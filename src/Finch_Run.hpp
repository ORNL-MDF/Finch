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
    bool srdf_format;

    Layer( Inputs& inputs, Grid<MemorySpace>& grid )
    {
        // Only construct if turned on - will otherwise default and immediately
        // return from any member functions
        if ( inputs.sampling.enabled )
        {
            srdf_format = inputs.sampling.format == "stork";
            solidification_data_ = sampling_type( inputs, grid, srdf_format );
        }
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
        int output_interval = inputs.time.output.interval;

        // update the temperature field
        for ( int n = 0; n < num_steps; ++n )
        {
            inputs.time_monitor.update();

            step( exec_space, time, dt, grid, beam, fd );

            // Update time monitoring
            if ( ( n + 1 ) % inputs.time.monitor.interval == 0 )
            {
                inputs.time_monitor.write( n );
            }

            // Write the current temperature field
            if ( ( n + 1 ) % output_interval == 0 )
            {
                grid.output( n, n * dt );
            }
        }
    }

    // Run a single timestep
    template <typename ExecutionSpace, typename SolverType>
    void step( ExecutionSpace exec_space, double& time, const double dt,
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

        // communicate halos
        grid.gather();
        // update boundaries
        grid.updateBoundaries();

        solidification_data_.update( grid, time );
    }

    auto getSolidificationData( Grid<MemorySpace> grid, MPI_Comm comm,
                                Sampling sampling_inputs, bool write_data )
    {
        return solidification_data_.get( grid, comm, sampling_inputs,
                                         write_data );
    }

    // Append next layer's solidification data to input_solidification_data
    void appendSolidificationData(
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>&
            input_solidification_data,
        MPI_Comm comm, Grid<MemorySpace> grid, Sampling sampling_inputs,
        std::vector<int>& first_value_finch, std::vector<int>& last_value_finch,
        int finch_file_num, const int num_finch_simulations )
    {
        // Time-temperature history from the Finch simulation performed for this
        // layer
        auto new_layer_data =
            solidification_data_.get( grid, comm, sampling_inputs, false );
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

    auto writeSolidificationData( Grid<memory_space>& grid,
                                  Sampling sampling_inputs, MPI_Comm comm )
    {
        return solidification_data_.write( sampling_inputs, grid, comm );
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
