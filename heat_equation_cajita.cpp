#include <array>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <string>

#include <Cajita.hpp>
#include <Kokkos_Core.hpp>

#include "movingBeam/moving_beam.hpp"

// we are currently assuming each boundary is a uniform Dirichlet condition
template <class ExecutionSpace, class grid_t, class array_t>
void updateBoundaries( ExecutionSpace exec_space, grid_t local_grid,
                       array_t& field )
{
    for ( int d = 0; d < 3; d++ )
    {
        for ( int dir = -1; dir < 2; dir += 2 )
        {
            std::array<int, 3> plane = { 0, 0, 0 };

            plane[d] = dir;

            auto boundary_space = local_grid->boundaryIndexSpace(
                Cajita::Own(), Cajita::Cell(), plane );

            Cajita::grid_parallel_for(
                "boundary_grid_for", exec_space, boundary_space,
                KOKKOS_LAMBDA( const int i, const int j, const int k ) {
                    // Apply boundary condition
                    field( i, j, k, 0 ) =
                        field( i - plane[0], j - plane[1], k - plane[2], 0 );
                } );
        }
    }
}

void createGrid()
{
    int comm_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

    using exec_space = Kokkos::DefaultExecutionSpace;
    using device_type = exec_space::device_type;

    // set up block decomposition
    int comm_size;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

    // Create the global mesh
    double cell_size = 25e-6;
    std::array<double, 3> global_low_corner = { -5e-4, -5e-4, -5e-4 };
    std::array<double, 3> global_high_corner = { 5.5e-3, 5.5e-3, 0 };
    auto global_mesh = Cajita::createUniformGlobalMesh(
        global_low_corner, global_high_corner, cell_size );

    // Create the global grid
    std::array<int, 3> ranks_per_dim = { comm_size, 1, 1 };
    MPI_Dims_create( comm_size, 3, ranks_per_dim.data() );
    Cajita::ManualBlockPartitioner<3> partitioner( ranks_per_dim );
    std::array<bool, 3> periodic = { false, false, false };

    auto global_grid =
        createGlobalGrid( MPI_COMM_WORLD, global_mesh, periodic, partitioner );

    // Create a local grid and local mesh with halo region
    unsigned halo_width = 1;
    auto local_grid = Cajita::createLocalGrid( global_grid, halo_width );
    auto local_mesh = Cajita::createLocalMesh<device_type>( *local_grid );

    auto owned_space = local_grid->indexSpace( Cajita::Own(), Cajita::Cell(),
                                               Cajita::Local() );

    // Create cell array layout for finite difference calculations
    auto layout =
        createArrayLayout( global_grid, halo_width, 1, Cajita::Cell() );

    std::string name( "temperature" );
    auto T = Cajita::createArray<double, device_type>( name, layout );
    Cajita::ArrayOp::assign( *T, 300.0, Cajita::Ghost() );
    Cajita::ArrayOp::assign( *T, 300.0, Cajita::Own() );

    auto T_view = T->view();

    auto T_halo = createHalo( Cajita::NodeHaloPattern<3>(), halo_width, *T );

    // update boundaries and halos
    updateBoundaries( exec_space(), local_grid, T_view );
    T_halo->gather( exec_space(), *T );

    // create and array to store previous temperature for explicit udpate
    auto T0 = Cajita::createArray<double, device_type>( name, layout );
    auto T0_view = T0->view();

    // Solve heat conduction from point source
    double dt = 1e-6;
    double end_time = 4e-3;
    int numSteps = static_cast<int>( end_time / dt );

    // properties
    double rho = 7600.0;
    double specific_heat = 750.0;
    double kappa = 30.0;

    double alpha = kappa / ( rho * specific_heat );
    double alpha_dt_dx2 = alpha * dt / ( cell_size * cell_size );
    double dt_rho_cp = dt / ( rho * specific_heat );

    // gaussian heat source parameters (sigma is std dev of gaussian)
    double eta = 0.35;
    double power = 195.0;
    double sigma[3] = { 50e-6, 50e-6, 60e-6 };
    double r[3] = { sigma[0] / sqrt( 2.0 ), sigma[1] / sqrt( 2.0 ),
                    sigma[2] / sqrt( 2.0 ) };

    // volumetric intensity of the heat source
    double I = 2.0 * eta * power / ( M_PI * sqrt( M_PI ) * r[0] * r[1] * r[2] );

    // parameters for a moving beam in x-direction
    movingBeam beam;

    // write 10 files to create a movie
    int num_output_steps = 10;
    int output_interval = static_cast<int>( numSteps / num_output_steps );
    double time = 0.0;

    for ( int step = 0; step < numSteps; ++step )
    {
        time += dt;
        // update beam position
        beam.move( time );
        std::vector<double> beam_pos = beam.position();
        // FIXME: not used
        double beam_power = beam.power();

        // store previous value for explicit update
        Kokkos::deep_copy( T0_view, T_view );

        // Solve finite difference
        Cajita::grid_parallel_for(
            "local_grid_for", exec_space(), owned_space,
            KOKKOS_LAMBDA( const int i, const int j, const int k ) {
                double loc[3];
                int idx[3] = { i, j, k };
                local_mesh.coordinates( Cajita::Cell(), idx, loc );

                loc[0] = fabs( loc[0] - beam_pos[0] );
                loc[1] = fabs( loc[1] - beam_pos[1] );
                loc[2] = fabs( loc[2] - beam_pos[2] );

                double f = ( loc[0] * loc[0] / r[0] / r[0] ) +
                           ( loc[1] * loc[1] / r[1] / r[1] ) +
                           ( loc[2] * loc[2] / r[2] / r[2] );

                double Q = I * exp( -f ) * dt_rho_cp;

                double laplacian =
                    ( -6.0 * T0_view( i, j, k, 0 ) + T0_view( i - 1, j, k, 0 ) +
                      T0_view( i + 1, j, k, 0 ) + T0_view( i, j - 1, k, 0 ) +
                      T0_view( i, j + 1, k, 0 ) + T0_view( i, j, k - 1, 0 ) +
                      T0_view( i, j, k + 1, 0 ) ) *
                    alpha_dt_dx2;

                T_view( i, j, k, 0 ) = T0_view( i, j, k, 0 ) + laplacian + Q;
            } );

        // update boundaries and halos
        updateBoundaries( exec_space(), local_grid, T_view );
        T_halo->gather( exec_space(), *T );

        // Write the current temperature field
        if ( step % output_interval == 0 )
        {
            Cajita::Experimental::BovWriter::writeTimeStep( step, step * dt,
                                                            *T );
        }
    }

    // Write the final temperature field
    Cajita::Experimental::BovWriter::writeTimeStep( numSteps, numSteps * dt,
                                                    *T );
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
