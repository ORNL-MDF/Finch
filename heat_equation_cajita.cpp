#include <Cajita.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <array>

// we are currently assuming each boundary is a uniform Dirichlet condition
template <class ExecutionSpace, class grid_t, class array_t>
void updateBoundaries(ExecutionSpace exec_space, grid_t local_grid, array_t& field)
{
    for ( int d = 0; d < 3; d++ )
    {
        for ( int dir = -1; dir < 2; dir += 2 )
        {
            std::array<int, 3> plane = { 0, 0, 0 };

            plane[d] = dir;

            auto boundary_space = local_grid->boundaryIndexSpace(
                Cajita::Own(), Cajita::Cell(), plane);

            Cajita::grid_parallel_for
            (
                "boundary_grid_for",
                exec_space,
                boundary_space,
                KOKKOS_LAMBDA( const int i, const int j, const int k )
                {
                    field( i, j, k, 0 ) = 1.0;
                }
            );
        }
    }
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
    double cell_size = 1;
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
    Cajita::ArrayOp::assign( *T, 0.0, Cajita::Ghost() );
    Cajita::ArrayOp::assign( *T, 0.0, Cajita::Own() );

    auto T_view = T->view();

    auto internal_space = local_grid->indexSpace(
            Cajita::Own(), Cajita::Cell(), Cajita::Local());

    auto halo =
        createHalo( Cajita::NodeHaloPattern<3>(), halo_width, *T );

    double alpha = 0.1;
    double dt = 0.1;
    int numSteps = 1000;

    updateBoundaries(exec_space(), local_grid, T_view);

    double alphaDtByDxSqr = alpha * dt / (cell_size * cell_size);

    for (int step = 0; step < numSteps; ++step)
    {
        // Solve finite difference
        Cajita::grid_parallel_for
        (
            "local_grid_for",
            exec_space(),
            internal_space,
            KOKKOS_LAMBDA( const int i, const int j, const int k )
            {
                {
                    double laplacian =
                      - 6.0*T_view(i, j, k, 0)
                      + T_view(i - 1, j, k, 0) + T_view(i + 1, j, k, 0)
                      + T_view(i, j - 1, k, 0) + T_view(i, j + 1, k, 0)
                      + T_view(i, j, k - 1, 0) + T_view(i, j, k + 1, 0);
                    
                    T_view( i, j, k, 0 ) += laplacian*alphaDtByDxSqr;
                }
            }
        );
        
        // update the boundary conditions
        updateBoundaries( exec_space(), local_grid, T_view);

        // Exchange halo values
        halo->gather(  exec_space(), *T );
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
