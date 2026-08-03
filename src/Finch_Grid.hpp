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
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_Boundary.hpp>

namespace Finch
{

// Info macro for writing on master
#define Info                                                                   \
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

    // Construct from components
    Grid( MPI_Comm comm, const double cell_size,
          std::array<double, 3> global_low_corner,
          std::array<double, 3> global_high_corner,
          std::array<int, 3> ranks_per_dim, std::array<std::string, 6> bc_types,
          Kokkos::Array<double, 6> bc_values, const double initial_temperature,
          const exec_space& execution_space = exec_space{} )
        : exec_space_( execution_space )
        , boundary( Boundary( bc_types, bc_values, cell_size ) )
    {
        initialize( comm, cell_size, global_low_corner, global_high_corner,
                    ranks_per_dim, initial_temperature );

        // Create boundaries
        boundary.create( local_grid, entity_type{} );

        // Initialize boundaries and halo
        updateBoundaries();
        gather();
    }

    // Constructor where no BC values are needed.
    Grid( MPI_Comm comm, const double cell_size,
          std::array<double, 3> global_low_corner,
          std::array<double, 3> global_high_corner,
          std::array<int, 3> ranks_per_dim, std::array<std::string, 6> bc_types,
          const double initial_temperature,
          const exec_space& execution_space = exec_space{} )
        : exec_space_( execution_space )
        , boundary( Boundary( bc_types, cell_size ) )
    {
        initialize( comm, cell_size, global_low_corner, global_high_corner,
                    ranks_per_dim, initial_temperature );

        // Create boundaries
        boundary.create( local_grid, entity_type{} );

        // Initialize boundaries and halo
        updateBoundaries();
        gather();
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
        Info << "Ranks per dimension: ";
        for ( int d = 0; d < 3; ++d )
            Info << ranks_per_dim_manual[d] << " ";
        Info << std::endl;

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

    auto getTemperature() { return T->view(); }

    auto getPreviousTemperature() { return T0->view(); }

    const exec_space& executionSpace() const { return exec_space_; }

    // Start a new explicit step without copying the field. After the swap T0
    // is the completed previous state and T is the reusable output buffer.
    void swapTemperatureFields() { std::swap( T, T0 ); }

    void output( const int step, const double time )
    {
        Cabana::Grid::Experimental::BovWriter::writeTimeStep( step, time, *T );
    }

    void updateBoundaries()
    {
        Kokkos::Profiling::ScopedRegion region( "Finch::physical_boundaries" );
        auto T_view = getTemperature();
        boundary.update( exec_space_, T_view );
    }

    void gather()
    {
        Kokkos::Profiling::ScopedRegion region( "Finch::halo_exchange" );
        halo->gather( exec_space_, *T );
    }

    MPI_Comm getComm() { return local_grid->globalGrid().comm(); }

  protected:
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

    // Temperature field.
    std::shared_ptr<array_type> T;

    // Previous temperature field.
    std::shared_ptr<array_type> T0;

    //! Boundary conditions.
    Boundary boundary;
};

} // namespace Finch

#endif
