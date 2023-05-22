#include <Cajita.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <array>
#include <math.h>

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
                    // Dirichlet boundary condition example
                    //field( i, j, k, 0 ) = 1.0;

                    // Neumann boundary condition example
                    field( i, j, k, 0 ) = field( i - plane[0], j - plane[1], k - plane[2], 0 );
                }
            );
        }
    }
}

void createGrid()
{
    int comm_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

    using exec_space = Kokkos::DefaultExecutionSpace;
    using device_type = exec_space::device_type;
    using memory_space = device_type::memory_space;


    // set up block decomposition
    int comm_size;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );


    // Create the global mesh
    double cell_size = 1e-5;
    std::array<double, 3> global_low_corner  = { -5e-4, -5e-4, -5e-4 };
    std::array<double, 3> global_high_corner = { 5e-4, 5e-4, 0 };
    auto global_mesh = Cajita::createUniformGlobalMesh(
        global_low_corner, global_high_corner, cell_size );


    // Create the global grid
    std::array<int, 3> ranks_per_dim = { comm_size, 1, 1 };
    MPI_Dims_create( comm_size, 3, ranks_per_dim.data() );
    Cajita::ManualBlockPartitioner<3> partitioner( ranks_per_dim );
    std::array<bool, 3> periodic = { false, false, false };

    auto global_grid = createGlobalGrid( MPI_COMM_WORLD, global_mesh,
                                         periodic, partitioner );


    // Create a local grid and local mesh with halo region
    unsigned halo_width = 1;
    auto local_grid = Cajita::createLocalGrid( global_grid, halo_width );
    auto local_mesh = Cajita::createLocalMesh<device_type>( *local_grid );

    auto owned_space = local_grid->indexSpace(
            Cajita::Own(), Cajita::Cell(), Cajita::Local());

    // Create cell array layout for finite difference calculations
    auto layout =
        createArrayLayout( global_grid, halo_width, 1, Cajita::Cell() );

    std::string name( "temperature" );
    auto T = Cajita::createArray<double, device_type>( name, layout );
    Cajita::ArrayOp::assign( *T, 300.0, Cajita::Ghost() );
    Cajita::ArrayOp::assign( *T, 300.0, Cajita::Own() );

    auto T_view = T->view();

    auto T_halo =
        createHalo( Cajita::NodeHaloPattern<3>(), halo_width, *T );

    // update physical and processor boundaries
    updateBoundaries(exec_space(), local_grid, T_view);
    T_halo->gather(  exec_space(), *T );


    // create and array to store previous temperature for explicit udpate
    auto T0 = Cajita::createArray<double, device_type>( name, layout );
    auto T0_view = T0->view();


    // Solve heat conduction from point source
    double dt = 1e-6;
    double end_time = 1e-3;
    int numSteps = static_cast<int>(end_time / dt);

    // properties
    double rho = 1000.0;
    double specific_heat = 1000.0;
    double kappa = 10.0;

    double alpha = kappa / (rho * specific_heat);
    double alpha_dt_dx2 = alpha * dt / (cell_size * cell_size);
    double dt_rho_cp = dt / (rho * specific_heat);

    // gaussian heat source parameters (sigma is std dev of gaussian)
    double eta = 0.1;
    double power = 100.0;
    double sigma[3] = {50e-6, 50e-6, 25e-6};
    double r[3] = {sigma[0] / sqrt(2.0), sigma[1] / sqrt(2.0), sigma[2] / sqrt(2.0)};

    double I = 2.0 * eta*power / (M_PI * sqrt(M_PI) * r[0] * r[1] * r[2]);
    
    for (int step = 0; step < numSteps; ++step)
    {
        // store previous value for explicit update
        Kokkos::deep_copy(T0_view, T_view);

        // Solve finite difference
        Cajita::grid_parallel_for
        (
            "local_grid_for",
            exec_space(),
            owned_space,
            KOKKOS_LAMBDA( const int i, const int j, const int k )
            {
                double loc[3];
                int idx[3] = { i, j, k };
                local_mesh.coordinates( Cajita::Cell(), idx, loc);

                double f =
                    (loc[0] * loc[0] / r[0] / r[0])
                  + (loc[1] * loc[1] / r[1] / r[1])
                  + (loc[2] * loc[2] / r[2] / r[2]);

                double Q = I * exp(-f) * dt_rho_cp;

                double laplacian =
                (
                  - 6.0*T0_view(i, j, k, 0)
                  + T0_view(i - 1, j, k, 0) + T0_view(i + 1, j, k, 0)
                  + T0_view(i, j - 1, k, 0) + T0_view(i, j + 1, k, 0)
                  + T0_view(i, j, k - 1, 0) + T0_view(i, j, k + 1, 0)
                ) * alpha_dt_dx2;
                
                T_view( i, j, k, 0 ) = T0_view( i, j, k, 0) + laplacian + Q;
            }
        );
        
        // update the boundary conditions
        updateBoundaries( exec_space(), local_grid, T_view);

        // Exchange halo values
        T_halo->gather(  exec_space(), *T );
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
