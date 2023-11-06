#ifndef Grid_H
#define Grid_H

#include <Cajita.hpp>
#include <Kokkos_Core.hpp>

template <typename MemorySpace>
class Grid
{
  public:
    // Kokkos memory space.
    using memory_space = MemorySpace;
    // Default Kokkos execution space for this memory space.
    using exec_space = typename MemorySpace::execution_space;
    using entity_type = Cajita::Node;
    using mesh_type = Cajita::UniformMesh<double>;
    using array_type =
        Cajita::Array<double, entity_type, mesh_type, memory_space>;
    //! MPI rank and number of ranks
    int comm_rank, comm_size;

    Grid( MPI_Comm comm, const double cell_size,
          std::array<double, 3> global_low_corner,
          std::array<double, 3> global_high_corner,
          const double initial_temperature )
    {
        auto global_mesh = Cajita::createUniformGlobalMesh(
            global_low_corner, global_high_corner, cell_size );

        // set up block decomposition
        MPI_Comm_size( comm, &comm_size );
        MPI_Comm_rank( comm, &comm_rank );

        // Create the global grid
        std::array<int, 3> ranks_per_dim = { comm_size, 1, 1 };
        MPI_Dims_create( comm_size, 3, ranks_per_dim.data() );
        Cajita::ManualBlockPartitioner<3> partitioner( ranks_per_dim );
        std::array<bool, 3> periodic = { false, false, false };

        auto global_grid =
            createGlobalGrid( comm, global_mesh, periodic, partitioner );

        // Create a local grid and local mesh with halo region
        local_grid = Cajita::createLocalGrid( global_grid, halo_width );

        // Create cell array layout for finite difference calculations
        auto layout =
            createArrayLayout( global_grid, halo_width, 1, entity_type() );

        std::string name( "temperature" );
        T = Cajita::createArray<double, memory_space>( name, layout );
        Cajita::ArrayOp::assign( *T, initial_temperature, Cajita::Ghost() );

        // create an array to store previous temperature for explicit udpate
        // Note: this is an entirely separate array on purpose (no shallow copy)
        T0 = Cajita::createArray<double, memory_space>( name, layout );

        // Create halo
        halo = createHalo( Cajita::FaceHaloPattern<3>(), halo_width, *T );

        // Create boundaries
        createBoundaries();

        // Initialize boundaries and halo
        updateBoundaries();
        halo->gather( exec_space{}, *T );
    }

    auto getLocalMesh()
    {
        return Cajita::createLocalMesh<memory_space>( *local_grid );
    }

    auto getLocalGrid() { return local_grid; }

    auto getIndexSpace()
    {
        return local_grid->indexSpace( Cajita::Own(), entity_type(),
                                       Cajita::Local() );
    }

    auto getTemperature() { return T->view(); }

    auto getPreviousTemperature() { return T0->view(); }

    void output( const int step, const double time )
    {
        Cajita::Experimental::BovWriter::writeTimeStep( step, time, *T );
    }

    void createBoundaries()
    {
        // Generate the boundary condition index spaces.
        int count = 0;
        for ( int d = 0; d < 3; d++ )
        {
            for ( int dir = -1; dir < 2; dir += 2 )
            {
                boundary_planes[count] = { 0, 0, 0 };
                boundary_planes[count][d] = dir;

                // Get the boundary indices for this plane (each one is a
                // separate, contiguous index space).
                boundary_spaces[count] = local_grid->boundaryIndexSpace(
                    Cajita::Ghost(), entity_type(), boundary_planes[count][0],
                    boundary_planes[count][1], boundary_planes[count][2] );
                count++;
            }
        }
    }

    void updateBoundaries()
    {
        // Update the boundary on each face of the cube. Pass the vector of
        // index spaces to avoid launching 6 separate kernels.
        auto T_view = getTemperature();
        auto planes = boundary_planes;
        Cajita::grid_parallel_for(
            "boundary_update", exec_space{}, boundary_spaces,
            KOKKOS_LAMBDA( const int b, const int i, const int j,
                           const int k ) {
                T_view( i, j, k, 0 ) = T_view(
                    i - planes[b][0], j - planes[b][1], k - planes[b][2], 0 );
            } );
    }

    void gather() { halo->gather( exec_space{}, *T ); }

  protected:
    //! Halo and stencil width;
    unsigned halo_width = 1;
    //! Owned grid
    std::shared_ptr<Cajita::LocalGrid<Cajita::UniformMesh<double>>> local_grid;
    //! Halo
    std::shared_ptr<Cajita::Halo<memory_space>> halo;

    //! Temperature field.
    std::shared_ptr<array_type> T;
    //! Previous temperature field.
    std::shared_ptr<array_type> T0;

    //! Boundary indices for each plane.
    Kokkos::Array<Cajita::IndexSpace<3>, 6> boundary_spaces;
    // Boundary details.
    Kokkos::Array<Kokkos::Array<int, 3>, 6> boundary_planes;
};
#endif
