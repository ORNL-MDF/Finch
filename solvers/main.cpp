#include <array>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <string>

#include <Cajita.hpp>
#include <Kokkos_Core.hpp>

#include "Grid.hpp"
#include "MovingBeam/MovingBeam.hpp"
#include "Simulation.hpp"
#include "SolidificationData.hpp"

void run( int argc, char* argv[] )
{
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;

    // initialize the simulation
    Simulation db( MPI_COMM_WORLD, argc, argv );

    // initialize a moving beam
    MovingBeam beam( db.source.scan_path_file );

    // Define boundary condition details.
    std::array<std::string, 6> bc_types = { "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic" };

    // create the global mesh
    Grid<memory_space> grid(
        MPI_COMM_WORLD, db.space.cell_size, db.space.global_low_corner,
        db.space.global_high_corner, db.space.ranks_per_dim, bc_types,
        db.space.initial_temperature );

    auto local_mesh = grid.getLocalMesh();
    using entity_type = typename Grid<memory_space>::entity_type;

    // gaussian heat source parameters
    double r[3] = { db.source.two_sigma[0] / sqrt( 2.0 ),
                    db.source.two_sigma[1] / sqrt( 2.0 ),
                    db.source.two_sigma[2] / sqrt( 2.0 ) };

    double A_inv[3] = { 1.0 / r[0] / r[0], 1.0 / r[1] / r[1],
                        1.0 / r[2] / r[2] };

    double I0 = ( 2.0 * db.source.absorption ) /
                ( M_PI * sqrt( M_PI ) * r[0] * r[1] * r[2] );

    // cut off for 3 standard deviations from heat source center
    double f_max = log( 3 ) + 2 * log( 10 );

    // time stepping
    double& time = db.time.time;

    int num_steps = db.time.num_steps;

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

    double k_by_dx2 = ( db.properties.thermal_conductivity ) / ( dx * dx );

    // update the temperature field
    for ( int step = 0; step < num_steps; ++step )
    {
        db.time_monitor.update();

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
            KOKKOS_LAMBDA( const int i, const int j, const int k ) {
                double x = T0( i, j, k, 0 );

                // calculate linearized effective specific heat
                double rho_cp_eff =
                    ( x >= solidus && x <= liquidus ) ? rho_cp_Lf : rho_cp;

                double dt_by_rho_cp = dt / rho_cp_eff;

                // calculate diffusion term
                double laplacian =
                    ( -6.0 * T0( i, j, k, 0 ) + T0( i - 1, j, k, 0 ) +
                      T0( i + 1, j, k, 0 ) + T0( i, j - 1, k, 0 ) +
                      T0( i, j + 1, k, 0 ) + T0( i, j, k - 1, 0 ) +
                      T0( i, j, k + 1, 0 ) ) *
                    k_by_dx2 * dt_by_rho_cp;

                // calculate heating source term
                double q_dot = 0.0;

                if ( beam_power )
                {
                    double grid_loc[3];
                    double dist_to_beam[3];
                    int idx[3] = { i, j, k };

                    local_mesh.coordinates( entity_type(), idx, grid_loc );

                    dist_to_beam[0] = grid_loc[0] - beam_pos_x;
                    dist_to_beam[1] = grid_loc[1] - beam_pos_y;
                    dist_to_beam[2] = grid_loc[2] - beam_pos_z;

                    double f =
                        ( dist_to_beam[0] * dist_to_beam[0] * A_inv[0] ) +
                        ( dist_to_beam[1] * dist_to_beam[1] * A_inv[1] ) +
                        ( dist_to_beam[2] * dist_to_beam[2] * A_inv[2] );

                    if ( f < f_max )
                    {
                        q_dot = I0 * beam_power * exp( -f ) * dt_by_rho_cp;
                    }
                }

                T( i, j, k, 0 ) = T0( i, j, k, 0 ) + laplacian + q_dot;
            } );

        // update boundaries
        grid.updateBoundaries();

        // communicate halos
        grid.gather();

        solidification_data.update();

        // Write the current temperature field (or, if step = num_steps - 1,
        // write the final temperature field)
        if ( ( ( step + 1 ) % db.time.output_interval == 0 ) ||
             ( step == num_steps - 1 ) )
        {
            grid.output( step + 1, ( step + 1 ) * db.time.time_step );
        }

        // Update time monitoring (or, if step = num_steps - 1, write the final
        // solution time)
        if ( ( ( step + 1 ) % db.time.monitor_interval == 0 ) ||
             ( step == num_steps - 1 ) )
        {
            db.time_monitor.write( step + 1 );
        }
    }

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
