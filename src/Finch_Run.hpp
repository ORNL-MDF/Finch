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

#ifndef Layer_H
#define Layer_H

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include "Finch_Grid.hpp"
#include "Finch_Inputs.hpp"
#include "Finch_SolidificationData.hpp"
#include "Finch_Solver.hpp"
#include "MovingBeam/Finch_MovingBeam.hpp"

namespace Finch
{

template <typename MemorySpace>
class Layer
{
  public:
    using memory_space = MemorySpace;
    using sampling_type = Finch::SolidificationData<memory_space>;
    sampling_type solidification_data_;

    Layer( Inputs& inputs, Grid<MemorySpace>& grid )
    {
        // Only construct if turned on - will otherwise default and immediately
        // return from any member functions
        if ( inputs.sampling.enabled )
            solidification_data_ = sampling_type( inputs, grid );
    }

    // Run the full timestepped loop
    template <typename ExecutionSpace, typename SolverType>
    void run( ExecutionSpace exec_space, Inputs& inputs,
              Grid<MemorySpace>& grid, MovingBeam& beam, SolverType& fd )
    {
        // time stepping
        double& time = inputs.time.time;
        int num_steps = inputs.time.num_steps;
        double dt = inputs.time.time_step;
        int output_interval = inputs.time.output_interval;

        // update the temperature field
        for ( int n = 0; n < num_steps; ++n )
        {
            inputs.time_monitor.update();

            step( exec_space, n, time, num_steps, dt, output_interval, grid,
                  beam, fd );

            // Update time monitoring (or, if n = num_steps - 1, write the
            // final solution time)
            if ( ( ( n + 1 ) % inputs.time.monitor_interval == 0 ) ||
                 ( n == num_steps - 1 ) )
            {
                inputs.time_monitor.write( n + 1 );
            }
        }
    }

    // Run a single timestep
    template <typename ExecutionSpace, typename SolverType>
    void step( ExecutionSpace exec_space, const int n, double& time,
               const int num_steps, const double dt, const int output_interval,
               Grid<MemorySpace> grid, MovingBeam& beam, SolverType& fd )
    {
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
        auto owned_space = grid.getIndexSpace();
        fd.solve( exec_space, owned_space, T, T0, beam_power, beam_pos );

        // update boundaries
        grid.updateBoundaries();

        // communicate halos
        grid.gather();

        solidification_data_.update( grid, time );

        // Write the current temperature field (or, if n = num_steps - 1,
        // write the final temperature field)
        if ( ( ( n + 1 ) % output_interval == 0 ) || ( n == num_steps - 1 ) )
        {
            grid.output( n + 1, ( n + 1 ) * dt );
        }
    }

    auto getSolidificationData() { return solidification_data_.get(); }

    auto writeSolidificationData() { return solidification_data_.write(); }
};

} // namespace Finch

#endif
