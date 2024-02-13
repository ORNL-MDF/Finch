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

#include <array>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <string>

#include <Kokkos_Core.hpp>

#include "Finch_Grid.hpp"
#include "Finch_Inputs.hpp"
#include "Finch_SolidificationData.hpp"
#include "Finch_Solver.hpp"
#include "MovingBeam/Finch_MovingBeam.hpp"

void run( int argc, char* argv[] )
{
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;

    // initialize the simulation
    Inputs db( MPI_COMM_WORLD, argc, argv );

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

    auto owned_space = grid.getIndexSpace();

    // time stepping
    double& time = db.time.time;
    int num_steps = db.time.num_steps;
    double dt = db.time.time_step;

    // Create the solver
    auto fd = createSolver( db, grid );

    // class for storing solidification data
    SolidificationData<memory_space> solidification_data( grid, db );

    // update the temperature field
    for ( int step = 0; step < num_steps; ++step )
    {
        db.time_monitor.update();

        time += dt;

        // update beam position
        beam.move( time );
        double beam_power = beam.power();
        double beam_pos[3];
        for ( std::size_t d = 0; d < 3; ++d )
            beam_pos[d] = beam.position( d );

        // Get temperature views;
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();

        // store previous value for explicit update
        Kokkos::deep_copy( T0, T );

        // Solve finite difference
        fd.solve( exec_space(), owned_space, T, T0, beam_power, beam_pos );

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
