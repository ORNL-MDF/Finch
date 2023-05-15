#include <Cajita.hpp>

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <array>

template <class array_t>
void setInternalField(array_t& field, double value)
{
    Cajita::ArrayOp::assign( field, value, Cajita::Ghost() );
    Cajita::ArrayOp::assign( field, value, Cajita::Own() );
}

void createGrid()
{
    int comm_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

    using exec_space = Kokkos::DefaultHostExecutionSpace;
    using device_type = exec_space::device_type;

    // set up block decomposition
    int comm_size;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
    std::array<int, 3> ranks_per_dim = { comm_size, 1, 1 };
    MPI_Dims_create( comm_size, 3, ranks_per_dim.data() );
    Cajita::ManualBlockPartitioner<3> partitioner( ranks_per_dim );
    std::array<bool, 3> periodic = { false, false, false };

    // Create the mesh and grid structures
    double cell_size = 0.1;
    std::array<int, 3> global_num_cell = { 40, 40, 40 };
    std::array<double, 3> global_low_corner = { 0, 0, 0 };
    std::array<double, 3> global_high_corner = {
        global_low_corner[0] + cell_size * global_num_cell[0],
        global_low_corner[1] + cell_size * global_num_cell[1],
        global_low_corner[2] + cell_size * global_num_cell[2] };

    auto global_mesh = Cajita::createUniformGlobalMesh(
        global_low_corner, global_high_corner, global_num_cell );

    auto global_grid = createGlobalGrid( MPI_COMM_WORLD, global_mesh,
                                         periodic, partitioner );

    // create local grid with halo region for mpi communication
    unsigned halo_width = 1;
    auto local_grid = Cajita::createLocalGrid(global_grid, halo_width);

    // create cell array layout
    auto layout =
        createArrayLayout( global_grid, halo_width, 1, Cajita::Cell() );

    std::string name( "temperature" );
    auto T = Cajita::createArray<double, device_type>( name, layout );

    setInternalField(*T, 0.0);

    auto T_view = T->view();

    // update boundary conditions
    // TODO : make a boudnary space that contains all boundaries
    auto boundary_space = local_grid->boundaryIndexSpace(
            Cajita::Own(), Cajita::Cell(), { 0, 0, 1 } );


    auto internal_space = local_grid->indexSpace(
            Cajita::Own(), Cajita::Cell(), Cajita::Local());

    auto halo =
        createHalo( Cajita::NodeHaloPattern<3>(), halo_width, *T );

    double alpha = 0.1;
    double dt = 0.1;
    int numSteps = 10000;
    
    for (int step = 0; step < numSteps; ++step)
    {
        // Update boundary conditions
        Cajita::grid_parallel_for
        (
            "boundary_grid_for",
            exec_space(),
            boundary_space,
            KOKKOS_LAMBDA( const int i, const int j, const int k )
            {
                T_view( i, j, k, 0 ) = 1.0;
            }
        );

        
        // Solve finite difference
        Cajita::grid_parallel_for
        (
            "local_grid_for",
            exec_space(),
            internal_space,
            KOKKOS_LAMBDA( const int i, const int j, const int k )
            {
                /*
                // if we want to avoid the boundaries
                if 
                (
                    (i >= internal_space.min(0)) 
                 && (j >= internal_space.min(1))
                 && (k >= internal_space.min(2))

                 && (i <= internal_space.max(0)) 
                 && (j <= internal_space.max(1)) 
                 && (k <= internal_space.max(2)) 
                )
                */
                {
                    double laplacian =
                      - 6.0*T_view(i, j, k, 0)
                      + T_view(i - 1, j, k, 0) + T_view(i + 1, j, k, 0)
                      + T_view(i, j - 1, k, 0) + T_view(i, j + 1, k, 0)
                      + T_view(i, j, k - 1, 0) + T_view(i, j, k + 1, 0);
                    
                    T_view( i, j, k, 0 ) += laplacian * alpha * dt;
                }
            }
        );

        // Exchange halo values
        halo->gather(  exec_space(), *T );
        //halo->scatter( exec_space(), Cajita::ScatterReduce::Replace(), *T );
    }

    Cajita::Experimental::BovWriter::writeTimeStep(0, 0, *T);
}

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    {
        Kokkos::ScopeGuard scope_guard( argc, argv );

        createGrid();
    }
    MPI_Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
