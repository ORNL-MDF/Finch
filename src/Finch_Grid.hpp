/****************************************************************************
 * Copyright (c) 2024 by Oak Ridge National Laboratory                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of Finch. Finch is distributed under a                 *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

/*!
  \file Grid.hpp
  \brief Background grid on which to solve heat transport
*/

#ifndef Grid_H
#define Grid_H

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_Boundary.hpp>

namespace Finch
{

// Rank-zero stream used for decomposition reporting.
#define FINCH_INFO                                                             \
    if ( comm_rank == 0 )                                                      \
    std::cout

template <typename MemorySpace>
class Grid
{
  public:
    // Kokkos memory space
    using memory_space = MemorySpace;

    // Default Kokkos execution space for this memory space
    using exec_space = typename MemorySpace::execution_space;

    using entity_type = Cabana::Grid::Node;

    using mesh_type = Cabana::Grid::UniformMesh<double>;
    using local_mesh_type = Cabana::Grid::LocalMesh<memory_space, mesh_type>;

    using array_type =
        Cabana::Grid::Array<double, entity_type, mesh_type, memory_space>;
    using view_type = typename array_type::view_type;

    int comm_rank, comm_size;

    Grid( MPI_Comm comm, const double cell_size,
          std::array<double, 3> global_low_corner,
          std::array<double, 3> global_high_corner,
          std::array<int, 3> ranks_per_dim,
          const BoundaryConditions& boundary_conditions,
          const double initial_temperature,
          const exec_space& execution_space = exec_space{} )
        : exec_space_( execution_space )
        , boundary( Boundary( boundary_conditions, cell_size ) )
    {
        initialize( comm, cell_size, global_low_corner, global_high_corner,
                    ranks_per_dim, initial_temperature );

        boundary.create( local_grid, entity_type{} );

        // The material-aware physical boundary initialization and first halo
        // exchange occur when the solver is created.
    }

    void initialize( MPI_Comm comm, const double cell_size,
                     std::array<double, 3> global_low_corner,
                     std::array<double, 3> global_high_corner,
                     std::array<int, 3> ranks_per_dim,
                     const double initial_temperature )
    {
        // set up block decomposition
        MPI_Comm_size( comm, &comm_size );
        MPI_Comm_rank( comm, &comm_rank );

        MPI_Dims_create( comm_size, 3, ranks_per_dim.data() );
        Cabana::Grid::ManualBlockPartitioner<3> partitioner( ranks_per_dim );
        std::array<bool, 3> periodic = { false, false, false };

        std::array<int, 3> ranks_per_dim_manual =
            partitioner.ranksPerDimension( comm, ranks_per_dim );

        // print the created decomposition
        FINCH_INFO << "Ranks per dimension: ";
        for ( int d = 0; d < 3; ++d )
            FINCH_INFO << ranks_per_dim_manual[d] << " ";
        FINCH_INFO << std::endl;

        // create global mesh and global grid
        auto global_mesh = Cabana::Grid::createUniformGlobalMesh(
            global_low_corner, global_high_corner, cell_size );

        auto global_grid =
            createGlobalGrid( comm, global_mesh, periodic, partitioner );

        // create a local grid and local mesh with halo region
        local_grid = Cabana::Grid::createLocalGrid( global_grid, halo_width );

        // create layout for finite difference calculations
        auto layout =
            createArrayLayout( global_grid, halo_width, 1, entity_type() );

        T = Cabana::Grid::createArray<double, memory_space>( "temperature",
                                                             layout );
        Cabana::Grid::ArrayOp::assign( *T, initial_temperature,
                                       Cabana::Grid::Ghost() );

        // create an array to store previous temperature for explicit update
        // Note: this is an entirely separate array on purpose (no shallow copy)
        T0 = Cabana::Grid::createArray<double, memory_space>( "temperature",
                                                              layout );

        // create halo
        halo = createHalo( Cabana::Grid::FaceHaloPattern<3>(), halo_width, *T );
    }

    auto getLocalMesh()
    {
        return Cabana::Grid::createLocalMesh<memory_space>( *local_grid );
    }

    auto getLocalGrid() { return local_grid; }

    auto getIndexSpace()
    {
        return local_grid->indexSpace( Cabana::Grid::Own(), entity_type(),
                                       Cabana::Grid::Local() );
    }

    auto getGhostedIndexSpace()
    {
        return local_grid->indexSpace( Cabana::Grid::Ghost(), entity_type(),
                                       Cabana::Grid::Local() );
    }

    auto getTemperature() { return T->view(); }

    auto getPreviousTemperature() { return T0->view(); }

    const exec_space& executionSpace() const { return exec_space_; }

    // Start a new explicit step without copying the field. After the swap T0
    // is the completed previous state and T is the reusable output buffer.
    void swapTemperatureFields() { std::swap( T, T0 ); }

    void output( const int step, const double time )
    {
        Cabana::Grid::Experimental::BovWriter::writeTimeStep( step, time, *T );
        writeXdmfSidecar( step, time );
    }

    template <typename MaterialModel>
    void initializeBoundaries( const MaterialModel& material )
    {
        Kokkos::Profiling::ScopedRegion region( "Finch::physical_boundaries" );
        auto T_view = getTemperature();
        boundary.update( exec_space_, T_view, T_view, material );
    }

    template <typename MaterialModel>
    void updateBoundaries( const MaterialModel& material )
    {
        Kokkos::Profiling::ScopedRegion region( "Finch::physical_boundaries" );
        auto T_view = getTemperature();
        auto T0_view = getPreviousTemperature();
        boundary.update( exec_space_, T_view, T0_view, material );
    }

    void gather()
    {
        Kokkos::Profiling::ScopedRegion region( "Finch::halo_exchange" );
        halo->gather( exec_space_, *T );
    }

    // Point-centered visualization blocks overlap by one node so Fides can
    // construct cells across rank interfaces. Face halos already provide all
    // overlap values for a one-dimensional decomposition. With decomposition
    // in multiple dimensions, update edge and corner ghosts only when field
    // output is requested.
    void gatherForFieldOutput()
    {
        const auto& global_grid = local_grid->globalGrid();
        int decomposed_dimensions = 0;
        for ( int d = 0; d < 3; ++d )
            decomposed_dimensions += global_grid.dimNumBlock( d ) > 1;

        if ( decomposed_dimensions < 2 )
            return;

        Kokkos::Profiling::ScopedRegion region(
            "Finch::field_output_halo_exchange" );
        if ( !field_output_halo )
            field_output_halo = createHalo( Cabana::Grid::NodeHaloPattern<3>(),
                                            halo_width, *T );
        field_output_halo->gather( exec_space_, *T );
    }

    MPI_Comm getComm() { return local_grid->globalGrid().comm(); }

  protected:
    void writeXdmfSidecar( const int step, const double time ) const
    {
        if ( comm_rank != 0 )
            return;

        const auto& global_grid = T->layout()->localGrid()->globalGrid();
        const auto& global_mesh = global_grid.globalMesh();

        std::array<long, 3> dimensions;
        for ( int d = 0; d < 3; ++d )
            dimensions[d] =
                global_grid.globalNumEntity( Cabana::Grid::Cell(), d ) + 1;

        std::ostringstream prefix;
        prefix << "grid_" << T->label() << '_' << std::setfill( '0' )
               << std::setw( 6 ) << step;

        std::ofstream xdmf( prefix.str() + ".xmf" );
        xdmf << std::setprecision( std::numeric_limits<double>::max_digits10 );
        xdmf << "<?xml version=\"1.0\" ?>\n"
             << "<Xdmf Version=\"2.0\">\n"
             << "  <Domain>\n"
             << "    <Grid Name=\"" << T->label()
             << "\" GridType=\"Uniform\">\n"
             << "      <Time Value=\"" << time << "\"/>\n"
             << "      <Topology TopologyType=\"3DCoRectMesh\" "
                "Dimensions=\""
             << dimensions[2] << ' ' << dimensions[1] << ' ' << dimensions[0]
             << "\"/>\n"
             << "      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
             << "        <DataItem Dimensions=\"3\" NumberType=\"Float\" "
                "Precision=\"8\" Format=\"XML\">"
             << global_mesh.lowCorner( 0 ) << ' ' << global_mesh.lowCorner( 1 )
             << ' ' << global_mesh.lowCorner( 2 ) << "</DataItem>\n"
             << "        <DataItem Dimensions=\"3\" NumberType=\"Float\" "
                "Precision=\"8\" Format=\"XML\">"
             << global_mesh.cellSize( 0 ) << ' ' << global_mesh.cellSize( 1 )
             << ' ' << global_mesh.cellSize( 2 ) << "</DataItem>\n"
             << "      </Geometry>\n"
             << "      <Attribute Name=\"" << T->label()
             << "\" AttributeType=\"Scalar\" Center=\"Node\">\n"
             << "        <DataItem Dimensions=\"" << dimensions[2] << ' '
             << dimensions[1] << ' ' << dimensions[0]
             << "\" NumberType=\"Float\" Precision=\"8\" "
                "Endian=\"Little\" Format=\"Binary\">"
             << prefix.str() << ".dat</DataItem>\n"
             << "      </Attribute>\n"
             << "    </Grid>\n"
             << "  </Domain>\n"
             << "</Xdmf>\n";
    }

    // A single ordered execution-space instance is used for all operations on
    // this grid, including Cabana packing and unpacking.
    exec_space exec_space_;

    // Halo and stencil width;
    unsigned halo_width = 1;

    // Owned grid
    std::shared_ptr<Cabana::Grid::LocalGrid<Cabana::Grid::UniformMesh<double>>>
        local_grid;

    // Halo
    std::shared_ptr<Cabana::Grid::Halo<memory_space>> halo;

    // Full-connectivity halo created lazily only for point-field output on a
    // decomposition spanning multiple dimensions.
    std::shared_ptr<Cabana::Grid::Halo<memory_space>> field_output_halo;

    // Temperature field.
    std::shared_ptr<array_type> T;

    // Previous temperature field.
    std::shared_ptr<array_type> T0;

    //! Boundary conditions.
    Boundary boundary;
};

} // namespace Finch

#endif
