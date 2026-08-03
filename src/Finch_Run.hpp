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

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

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
        const double nominal_dt = inputs.time.time_step;

        exec_space.fence( "Finch run start" );
        MPI_Barrier( grid.getComm() );
        const double run_start = MPI_Wtime();
        inputs.time_monitor.reset();

        // update the temperature field
        for ( int n = 0; n < num_steps; ++n )
        {
            const double dt =
                std::min( nominal_dt, inputs.time.end_time - time );

            step( exec_space, time, dt, grid, beam, fd );

            // Update time monitoring
            if ( inputs.time.monitor.isDue( n + 1, num_steps ) )
            {
                exec_space.fence( "Finch progress monitor" );
                inputs.time_monitor.write( n + 1 );
            }

            // Write the current temperature field
            if ( inputs.time.output.isDue( n + 1, num_steps ) )
            {
                Kokkos::Profiling::ScopedRegion output_region(
                    "Finch::field_output" );
                exec_space.fence( "Finch field output" );
                grid.output( n + 1, time );
            }
        }

        exec_space.fence( "Finch run complete" );
        MPI_Barrier( grid.getComm() );
        const double local_elapsed = MPI_Wtime() - run_start;
        double min_elapsed = 0.0;
        double max_elapsed = 0.0;
        double sum_elapsed = 0.0;
        MPI_Reduce( &local_elapsed, &min_elapsed, 1, MPI_DOUBLE, MPI_MIN, 0,
                    grid.getComm() );
        MPI_Reduce( &local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0,
                    grid.getComm() );
        MPI_Reduce( &local_elapsed, &sum_elapsed, 1, MPI_DOUBLE, MPI_SUM, 0,
                    grid.getComm() );

        unsigned long long local_nodes = grid.getIndexSpace().size();
        unsigned long long global_nodes = 0;
        MPI_Reduce( &local_nodes, &global_nodes, 1, MPI_UNSIGNED_LONG_LONG,
                    MPI_SUM, 0, grid.getComm() );
        if ( grid.comm_rank == 0 )
        {
            const double throughput =
                max_elapsed > 0.0 ? static_cast<double>( global_nodes ) *
                                        num_steps / max_elapsed
                                  : 0.0;
            std::cout << "Performance summary: " << num_steps << " steps, "
                      << std::fixed << std::setprecision( 6 ) << max_elapsed
                      << " s wall time (rank min/avg/max " << min_elapsed << "/"
                      << sum_elapsed / grid.comm_size << "/" << max_elapsed
                      << "), " << std::scientific << throughput
                      << " node updates/s" << std::endl;
        }
    }

    // Run a single timestep
    template <typename ExecutionSpace, typename SolverType>
    void step( ExecutionSpace exec_space, double& time, const double dt,
               Grid<MemorySpace>& grid, MovingBeam& beam, SolverType& fd )
    {
        Kokkos::Profiling::ScopedRegion step_region( "Finch::timestep" );
        time += dt;

        // update beam position
        beam.move( time );
        double beam_power = beam.power();
        const auto& beam_pos = beam.position();

        // Ping-pong invariant: T0 is the completed previous field and T is
        // the output buffer for this step. T0 is not reused until sampling is
        // complete.
        grid.swapTemperatureFields();

        // Get temperature views;
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();

        // Solve finite difference
        auto owned_space = grid.getIndexSpace();
        fd.solve( exec_space, owned_space, T, T0, dt, beam_power, beam_pos );

        // update boundaries
        grid.updateBoundaries();

        // communicate halos
        grid.gather();

        solidification_data_.update( grid, time, dt );
    }

    auto getSolidificationData() { return solidification_data_.get(); }
    // Append next layer's solidification data to input_solidification_data
    void appendSolidificationData(
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>&
            input_solidification_data,
        std::vector<int>& first_value_finch, std::vector<int>& last_value_finch,
        int finch_file_num, const int num_finch_simulations )
    {
        // Time-temperature history from the Finch simulation performed for this
        // layer
        auto new_layer_data = solidification_data_.get();
        // Number of events and components in new layer time-temperature history
        const int events_this_layer = new_layer_data.extent( 0 );
        const int n_cmpts = new_layer_data.extent( 1 );
        // Number of events in currently stored time-temperature history-
        // first_value_finch provides offset from data stored for previous
        // layers if performing more than 1 finch simulation at a time
        int events_prev_layers;
        if ( ( finch_file_num == 0 ) || ( num_finch_simulations == 1 ) )
        {
            first_value_finch[finch_file_num] = 0;
            events_prev_layers = 0;
        }
        else
        {
            first_value_finch[finch_file_num] =
                last_value_finch[finch_file_num - 1];
            events_prev_layers = input_solidification_data.extent( 0 );
        }
        // Resize input_solidification_data to accommodate both any events
        // stored from previous layers and the events calculated from simulation
        // of this layer
        Kokkos::resize( input_solidification_data,
                        events_prev_layers + events_this_layer, n_cmpts );
        // Copy this layer's data into the return view
        for ( int i = 0; i < events_this_layer; i++ )
            for ( int j = 0; j < n_cmpts; j++ )
                input_solidification_data(
                    first_value_finch[finch_file_num] + i, j ) =
                    new_layer_data( i, j );
        // Set last_value_finch to bound the indices with time-temperature
        // history data for this layer
        last_value_finch[finch_file_num] =
            events_prev_layers + events_this_layer;
    }

    auto writeSolidificationData( Sampling sampling_inputs, MPI_Comm comm )
    {
        return solidification_data_.write( sampling_inputs, comm );
    }

    [[deprecated( "Use of getLowerSolidificationDataBounds() without a "
                  "communicator is deprecated." )]] std::array<double, 3>
    getLowerSolidificationDataBounds()
    {
        return solidification_data_.getLowerBounds( MPI_COMM_WORLD );
    }
    [[deprecated( "Use of getUpperSolidificationDataBounds() without a "
                  "communicator is deprecated." )]] std::array<double, 3>
    getUpperSolidificationDataBounds()
    {
        return solidification_data_.getUpperBounds( MPI_COMM_WORLD );
    }

    std::array<double, 3> getLowerSolidificationDataBounds( MPI_Comm comm )
    {
        return solidification_data_.getLowerBounds( comm );
    }
    std::array<double, 3> getUpperSolidificationDataBounds( MPI_Comm comm )
    {
        return solidification_data_.getUpperBounds( comm );
    }
};

} // namespace Finch

#endif
