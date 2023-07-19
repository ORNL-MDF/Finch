#include <array>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <string>

#include <Cajita.hpp>
#include <Kokkos_Core.hpp>

#include "grid.hpp"
#include "movingBeam/moving_beam.hpp"

void run()
{
    using exec_space = Kokkos::DefaultExecutionSpace;
    using device_type = exec_space::device_type;

    // Create the global mesh
    double cell_size = 25e-6;
    std::array<double, 3> global_low_corner = { -5e-4, -5e-4, -5e-4 };
    std::array<double, 3> global_high_corner = { 5.5e-3, 5.5e-3, 0 };
    double init_temp = 300.0;

    Grid<device_type> grid( MPI_COMM_WORLD, cell_size, global_low_corner,
                            global_high_corner, init_temp );

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
    double sigma[3] = { 50e-6, 50e-6, 60e-6 };
    double r[3] = { sigma[0] / sqrt( 2.0 ), sigma[1] / sqrt( 2.0 ),
                    sigma[2] / sqrt( 2.0 ) };

    // parameters for a moving beam in x-direction
    MovingBeam beam;

    // write 10 files to create a movie
    int num_output_steps = 10;
    int output_interval = static_cast<int>( numSteps / num_output_steps );
    double time = 0.0;

    for ( int step = 0; step < numSteps; ++step )
    {
        time += dt;
        // update beam position
        beam.move( time );
        double beam_power = beam.power();
        double beam_pos_x = beam.position( 0 );
        double beam_pos_y = beam.position( 1 );
        double beam_pos_z = beam.position( 2 );

        // Get temperature views;
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();

        // store previous value for explicit update
        Kokkos::deep_copy( T0, T );

        // Solve finite difference
        auto local_mesh = grid.getLocalMesh();
        Cajita::grid_parallel_for(
            "local_grid_for", exec_space(), grid.getIndexSpace(),
            KOKKOS_LAMBDA( const int i, const int j, const int k ) {
                double loc[3];
                int idx[3] = { i, j, k };
                local_mesh.coordinates( Cajita::Cell(), idx, loc );

                loc[0] = fabs( loc[0] - beam_pos_x );
                loc[1] = fabs( loc[1] - beam_pos_y );
                loc[2] = fabs( loc[2] - beam_pos_z );

                double f = ( loc[0] * loc[0] / r[0] / r[0] ) +
                           ( loc[1] * loc[1] / r[1] / r[1] ) +
                           ( loc[2] * loc[2] / r[2] / r[2] );

                // volumetric intensity of the heat source
                double I = 2.0 * eta * beam_power /
                           ( M_PI * sqrt( M_PI ) * r[0] * r[1] * r[2] );
                double Q = I * exp( -f ) * dt_rho_cp;

                double laplacian =
                    ( -6.0 * T0( i, j, k, 0 ) + T0( i - 1, j, k, 0 ) +
                      T0( i + 1, j, k, 0 ) + T0( i, j - 1, k, 0 ) +
                      T0( i, j + 1, k, 0 ) + T0( i, j, k - 1, 0 ) +
                      T0( i, j, k + 1, 0 ) ) *
                    alpha_dt_dx2;

                T( i, j, k, 0 ) = T0( i, j, k, 0 ) + laplacian + Q;
            } );

        // update boundaries and halos
        grid.updateBoundaries();
        grid.gather();

        // Write the current temperature field
        if ( step % output_interval == 0 )
        {
            grid.output( step, step * dt );
        }
    }

    // Write the final temperature field
    grid.output( numSteps, numSteps * dt );
}

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    {
        Kokkos::ScopeGuard scope_guard( argc, argv );

        run();
    }
    MPI_Finalize();

    return 0;
}
