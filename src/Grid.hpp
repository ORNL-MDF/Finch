#ifndef Grid_H
#define Grid_H

#include <Cajita.hpp>
#include <Kokkos_Core.hpp>

#include <Boundary.hpp>

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

    using entity_type = Cajita::Node;

    using mesh_type = Cajita::UniformMesh<double>;
    using local_mesh_type = Cajita::LocalMesh<memory_space, mesh_type>;

    using array_type =
        Cajita::Array<double, entity_type, mesh_type, memory_space>;
    using view_type = typename array_type::view_type;

    int comm_rank, comm_size;

    // Construct from components
    Grid( MPI_Comm comm, const double cell_size,
          std::array<double, 3> global_low_corner,
          std::array<double, 3> global_high_corner,
          std::array<int, 3> ranks_per_dim, std::array<std::string, 6> bc_types,
          Kokkos::Array<double, 6> bc_values, const double initial_temperature )
        : boundary( Boundary( bc_types, bc_values ) )
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
          const double initial_temperature )
        : boundary( Boundary( bc_types ) )
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
        Cajita::ManualBlockPartitioner<3> partitioner( ranks_per_dim );
        std::array<bool, 3> periodic = { false, false, false };

        std::array<int, 3> ranks_per_dim_manual =
            partitioner.ranksPerDimension( MPI_COMM_WORLD, ranks_per_dim );

        // print the created decomposition
        Info << "Ranks per dimension: ";
        for ( int d = 0; d < 3; ++d )
            Info << ranks_per_dim_manual[d] << " ";
        Info << std::endl;

        // create global mesh and global grid
        auto global_mesh = Cajita::createUniformGlobalMesh(
            global_low_corner, global_high_corner, cell_size );

        auto global_grid =
            createGlobalGrid( comm, global_mesh, periodic, partitioner );

        // create a local grid and local mesh with halo region
        local_grid = Cajita::createLocalGrid( global_grid, halo_width );

        // create layout for finite difference calculations
        auto layout =
            createArrayLayout( global_grid, halo_width, 1, entity_type() );

        std::string name( "temperature" );
        T = Cajita::createArray<double, memory_space>( name, layout );
        Cajita::ArrayOp::assign( *T, initial_temperature, Cajita::Ghost() );

        // create an array to store previous temperature for explicit udpate
        // Note: this is an entirely separate array on purpose (no shallow copy)
        T0 = Cajita::createArray<double, memory_space>( name, layout );

        // create halo
        halo = createHalo( Cajita::FaceHaloPattern<3>(), halo_width, *T );
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

    void updateBoundaries()
    {
        auto T_view = getTemperature();
        boundary.update( exec_space{}, T_view );
    }

    void gather() { halo->gather( exec_space{}, *T ); }

  protected:
    // Halo and stencil width;
    unsigned halo_width = 1;

    // Owned grid
    std::shared_ptr<Cajita::LocalGrid<Cajita::UniformMesh<double>>> local_grid;

    // Halo
    std::shared_ptr<Cajita::Halo<memory_space>> halo;

    // Temperature field.
    std::shared_ptr<array_type> T;

    // Previous temperature field.
    std::shared_ptr<array_type> T0;

    //! Boundary conditions.
    Boundary boundary;
};
#endif
