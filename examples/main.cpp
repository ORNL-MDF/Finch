#include <array>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <string>

#include <Cajita.hpp>
#include <Kokkos_Core.hpp>

#include "grid.hpp"
#include "movingBeam/moving_beam.hpp"
#include "simulation.hpp"
#include "solidificationData.hpp"

void run( int argc, char* argv[] )
{
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space  = exec_space::memory_space;

    // initialize the simulation
    Simulation db( MPI_COMM_WORLD, argc, argv );

    // initialize a moving beam
    MovingBeam beam(db.source.scan_path_file);

    // create the global mesh
    Grid<memory_space> grid( MPI_COMM_WORLD,
                             db.space.cell_size,
                             db.space.global_low_corner,
                             db.space.global_high_corner,
                             db.space.initial_temperature );

    auto local_mesh = grid.getLocalMesh();

    // gaussian heat source parameters
    double r[3] = { db.source.two_sigma[0] / sqrt( 2.0 ),
                    db.source.two_sigma[1] / sqrt( 2.0 ),
                    db.source.two_sigma[2] / sqrt( 2.0 ) };

    double I0 = ( 2.0 * db.source.absorption ) 
              / ( M_PI * sqrt( M_PI ) * r[0] * r[1] * r[2] );

    // time stepping
    double& time = db.time.time;

    int numSteps = static_cast<int>(
                        ( db.time.end_time - db.time.start_time )
                      / ( db.time.time_step ));

    int output_interval = static_cast<int>(
                        ( numSteps / db.time.num_output_steps ));
    
    // class for storing solidification data
    SolidificationData<memory_space> solidification_data( grid, db );

    // reference to thermophysical properties
    double dt = db.time.time_step;
    double dx = db.space.cell_size;
    double rho = db.properties.density;
    double cp = db.properties.specific_heat;
    double Lf = db.properties.latent_heat;
    double solidus = db.properties.solidus;
    double liquidus = db.properties.liquidus;

    double rho_cp = rho * cp;

    double rho_cp_Lf = rho * cp + rho * Lf / ( liquidus - solidus );

    double k_by_dx2 = 
        ( db.properties.thermal_conductivity )
      / ( dx * dx );


    // update the temperature field
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
        Cajita::grid_parallel_for(
            "local_grid_for", exec_space(), grid.getIndexSpace(),
            KOKKOS_LAMBDA( const int i, const int j, const int k )
            {
                double x = T0( i, j, k, 0 );

                // calculate linearized effective specific heat
                double rho_cp_eff = 
                    ( x >= solidus && x <= liquidus ) ? rho_cp_Lf : rho_cp;

                double dt_by_rho_cp = dt / rho_cp_eff;

                // calculate diffusion term
                double laplacian =
                    ( -6.0 * T0( i, j, k, 0 ) + 
                      T0( i - 1, j, k, 0 ) + T0( i + 1, j, k, 0 ) +
                      T0( i, j - 1, k, 0 ) + T0( i, j + 1, k, 0 ) +
                      T0( i, j, k - 1, 0 ) + T0( i, j, k + 1, 0 ) ) *
                    k_by_dx2 * dt_by_rho_cp;

                // calculate heating source term
                double q_dot = 0.0;

                if ( beam_power )
                {
                    double cell_loc[3], cell_dist_to_beam[3];
                    int idx[3] = { i, j, k };
                    local_mesh.coordinates( Cajita::Cell(), idx, cell_loc );

                    cell_dist_to_beam[0] = fabs( cell_loc[0] - beam_pos_x );
                    cell_dist_to_beam[1] = fabs( cell_loc[1] - beam_pos_y );
                    cell_dist_to_beam[2] = fabs( cell_loc[2] - beam_pos_z );

                    double f = ( cell_dist_to_beam[0] * cell_dist_to_beam[0] /
                                 r[0] / r[0] ) +
                               ( cell_dist_to_beam[1] * cell_dist_to_beam[1] /
                                 r[1] / r[1] ) +
                               ( cell_dist_to_beam[2] * cell_dist_to_beam[2] /
                                 r[2] / r[2] );

                    q_dot = I0 * beam_power * exp( -f ) * dt_by_rho_cp;
                }

                T( i, j, k, 0 ) = T0( i, j, k, 0 ) + laplacian + q_dot;

            } );

        // update boundaries
        grid.updateBoundaries();

        // communicate halos
        grid.gather();

        solidification_data.update();

        // Write the current temperature field
        if ( output_interval && (step % output_interval == 0 ) )
        {
            grid.output( step, step * db.time.time_step );
        }
    }

    // Write the final temperature field
    grid.output( numSteps, numSteps * db.time.time_step );

    // Write the temperature data used by ExaCA/other post-processing
    solidification_data.write();
}

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    {
        Kokkos::ScopeGuard scope_guard( argc, argv );

        run( argc, argv );
    }
    MPI_Finalize();

    return 0;
}
