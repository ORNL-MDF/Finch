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
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <variant>
#include <vector>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include "Finch_FieldOutput.hpp"
#include "Finch_Grid.hpp"
#include "Finch_Inputs.hpp"
#include "Finch_MeltPoolDimensions.hpp"
#include "Finch_SolidificationData.hpp"
#include "Finch_Solver.hpp"
#include "MovingBeam/Finch_MovingBeam.hpp"

namespace Finch
{

namespace TimeIntegration
{

struct Event
{
    double time = 0.0;
    bool field_output = false;
    bool melt_pool_dimensions = false;
};

struct Interval
{
    double start = 0.0;
    double end = 0.0;
    double time_step = 0.0;
    int steps = 0;
    bool field_output = false;
    bool melt_pool_dimensions = false;
};

inline double tolerance( const double a, const double b )
{
    return 64.0 * std::numeric_limits<double>::epsilon() *
           std::max( { 1.0, std::abs( a ), std::abs( b ) } );
}

inline int fittedStepCount( const double duration,
                            const double maximum_time_step )
{
    const double raw_steps = duration / maximum_time_step;
    const double roundoff = 128.0 * std::numeric_limits<double>::epsilon() *
                            std::max( 1.0, std::abs( raw_steps ) );
    if ( !std::isfinite( raw_steps ) ||
         raw_steps > static_cast<double>( std::numeric_limits<int>::max() ) )
        throw std::runtime_error(
            "Automatic time interval exceeds the supported step count" );
    return std::max( 1, static_cast<int>( std::ceil( raw_steps - roundoff ) ) );
}

inline void appendCountedEvents( std::vector<Event>& events,
                                 const double start_time, const double end_time,
                                 const int count, const bool field_output )
{
    if ( count == 0 )
        return;

    const double duration = end_time - start_time;
    for ( int output = 1; output <= count; ++output )
    {
        Event event;
        event.time = output == count
                         ? end_time
                         : start_time + duration *
                                            static_cast<double>( output ) /
                                            static_cast<double>( count );
        event.field_output = field_output;
        event.melt_pool_dimensions = !field_output;
        events.push_back( event );
    }
}

inline std::vector<Interval> createSchedule( const Inputs& inputs,
                                             const MovingBeam& beam )
{
    std::vector<Event> events;
    std::vector<double> scan_path_events;
    beam.appendDiscontinuityTimes( scan_path_events );
    events.reserve( scan_path_events.size() +
                    inputs.functions.field_output.execute.count +
                    inputs.functions.melt_pool_dimensions.execute.count + 1 );
    for ( const double event_time : scan_path_events )
        events.push_back( { event_time, false, false } );

    appendCountedEvents( events, inputs.time.start_time, inputs.time.end_time,
                         inputs.functions.field_output.execute.count, true );
    appendCountedEvents( events, inputs.time.start_time, inputs.time.end_time,
                         inputs.functions.melt_pool_dimensions.execute.count,
                         false );
    events.push_back( { inputs.time.end_time, false, false } );

    std::sort( events.begin(), events.end(),
               []( const Event& lhs, const Event& rhs )
               { return lhs.time < rhs.time; } );

    std::vector<Event> merged_events;
    merged_events.reserve( events.size() );
    for ( auto event : events )
    {
        if ( event.time <=
                 inputs.time.start_time +
                     tolerance( event.time, inputs.time.start_time ) ||
             event.time > inputs.time.end_time +
                              tolerance( event.time, inputs.time.end_time ) )
            continue;
        if ( std::abs( event.time - inputs.time.end_time ) <=
             tolerance( event.time, inputs.time.end_time ) )
            event.time = inputs.time.end_time;

        if ( !merged_events.empty() &&
             std::abs( event.time - merged_events.back().time ) <=
                 tolerance( event.time, merged_events.back().time ) )
        {
            merged_events.back().field_output |= event.field_output;
            merged_events.back().melt_pool_dimensions |=
                event.melt_pool_dimensions;
        }
        else
            merged_events.push_back( event );
    }

    std::vector<Interval> schedule;
    schedule.reserve( merged_events.size() );
    double interval_start = inputs.time.start_time;
    long long total_steps = 0;
    for ( const auto& event : merged_events )
    {
        const double duration = event.time - interval_start;
        const int steps =
            fittedStepCount( duration, inputs.time.maximum_time_step );
        total_steps += steps;
        if ( total_steps > std::numeric_limits<int>::max() )
            throw std::runtime_error(
                "Automatic time schedule exceeds the supported step count" );

        schedule.push_back( { interval_start, event.time, duration / steps,
                              steps, event.field_output,
                              event.melt_pool_dimensions } );
        interval_start = event.time;
    }
    return schedule;
}

} // namespace TimeIntegration

template <typename MemorySpace>
class Layer
{
  public:
    using memory_space = MemorySpace;
    using solidification_data_type = Finch::SolidificationData<memory_space>;
    solidification_data_type solidification_data_;
    using melt_pool_dimensions_type = Finch::MeltPoolDimensions<memory_space>;
    std::unique_ptr<melt_pool_dimensions_type> melt_pool_dimensions_;
    using field_output_type = Finch::FieldOutput<memory_space>;
    std::unique_ptr<field_output_type> field_output_;

    Layer( Inputs& inputs, Grid<MemorySpace>& grid )
    {
        if ( inputs.functions.solidification.enabled )
            solidification_data_ = solidification_data_type(
                inputs.functions.solidification, inputs, grid );
        if ( inputs.functions.melt_pool_dimensions.enabled )
            melt_pool_dimensions_ = std::make_unique<melt_pool_dimensions_type>(
                inputs.functions.melt_pool_dimensions,
                inputs.properties.solidus, inputs.properties.liquidus,
                grid.getComm() );
        if ( inputs.functions.field_output.enabled )
            field_output_ = std::make_unique<field_output_type>(
                inputs.functions.field_output, grid );
    }

    template <typename ExecutionSpace, typename... SolverTypes>
    void run( ExecutionSpace exec_space, Inputs& inputs,
              Grid<MemorySpace>& grid, MovingBeam& beam,
              std::variant<SolverTypes...>& solvers )
    {
        std::visit( [&]( auto& solver )
                    { run( exec_space, inputs, grid, beam, solver ); },
                    solvers );
    }

    // Run the full timestepped loop
    template <typename ExecutionSpace, typename SolverType>
    void run( ExecutionSpace exec_space, Inputs& inputs,
              Grid<MemorySpace>& grid, MovingBeam& beam, SolverType& fd )
    {
        // Build a host-only schedule once. It fits stable steps exactly between
        // scan-path transitions, physical output times, and the requested end
        // time without any per-step device reduction or MPI collective.
        double& time = inputs.time.time;
        const auto schedule = TimeIntegration::createSchedule( inputs, beam );
        int num_steps = 0;
        double minimum_dt = std::numeric_limits<double>::max();
        double maximum_dt = 0.0;
        for ( const auto& interval : schedule )
        {
            num_steps += interval.steps;
            minimum_dt = std::min( minimum_dt, interval.time_step );
            maximum_dt = std::max( maximum_dt, interval.time_step );
        }
        inputs.time.num_steps = num_steps;
        inputs.time_monitor.setNumSteps( num_steps );
        if ( grid.comm_rank == 0 )
            std::cout << "Automatic time schedule: " << num_steps
                      << " steps, dt min/max " << minimum_dt << "/"
                      << maximum_dt << " s, " << schedule.size()
                      << " fitted intervals" << std::endl;

        exec_space.fence( "Finch run start" );
        MPI_Barrier( grid.getComm() );
        const double run_start = MPI_Wtime();
        inputs.time_monitor.reset();

        int completed_steps = 0;
        for ( const auto& interval : schedule )
        {
            for ( int interval_step = 0; interval_step < interval.steps;
                  ++interval_step )
            {
                const double next_time =
                    interval_step + 1 == interval.steps
                        ? interval.end
                        : interval.start +
                              static_cast<double>( interval_step + 1 ) *
                                  interval.time_step;
                const double dt = next_time - time;

                step( exec_space, time, next_time, dt, grid, beam, fd );
                ++completed_steps;

                if ( interval_step + 1 == interval.steps )
                {
                    if ( melt_pool_dimensions_ &&
                         interval.melt_pool_dimensions )
                        melt_pool_dimensions_->update( exec_space, grid,
                                                       beam.direction(), time );

                    if ( interval.field_output )
                    {
                        Kokkos::Profiling::ScopedRegion output_region(
                            "Finch::field_output" );
                        exec_space.fence( "Finch field output" );
                        field_output_->write( grid, fd, beam, completed_steps,
                                              time );
                    }
                }

                if ( inputs.time.monitor.isDue( completed_steps, num_steps ) )
                {
                    exec_space.fence( "Finch progress monitor" );
                    const double integrated_source_power =
                        fd.integratedSourcePower( exec_space, grid, beam );
                    const auto source_state = fd.heatSourceState();
                    inputs.time_monitor.write( completed_steps );
                    if ( grid.comm_rank == 0 )
                    {
                        const double commanded_power = beam.power();
                        const double target_absorbed_power =
                            source_state.absorptivity * commanded_power;
                        std::cout
                            << "  Source power: commanded=" << std::defaultfloat
                            << commanded_power
                            << " W, absorptivity=" << source_state.absorptivity
                            << ", target absorbed=" << target_absorbed_power
                            << " W, integrated=" << integrated_source_power
                            << " W";
                        if ( target_absorbed_power > 0.0 )
                            std::cout << ", error="
                                      << 100.0 * ( integrated_source_power /
                                                       target_absorbed_power -
                                                   1.0 )
                                      << "%";
                        std::cout << std::endl;
                        if ( commanded_power > 0.0 &&
                             source_state.depth_feedback )
                        {
                            constexpr double radians_to_degrees =
                                57.2957795130823208768;
                            std::cout << "  Source geometry: measured depth="
                                      << source_state.measured_melt_depth
                                      << " m, effective depth="
                                      << source_state.effective_source_depth
                                      << " m, D4sigma equivalent/major/minor="
                                      << source_state.lateral_d4_sigma << "/"
                                      << source_state.lateral_d4_sigma_major
                                      << "/"
                                      << source_state.lateral_d4_sigma_minor
                                      << " m, profile azimuth="
                                      << source_state.profile_azimuth *
                                             radians_to_degrees
                                      << " deg, aspect ratio="
                                      << source_state.aspect_ratio << std::endl;
                        }
                    }
                }
            }
        }

        if ( field_output_ )
            field_output_->close();

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

        solidification_data_.write( grid.getComm() );
    }

    // Run a single timestep
    template <typename ExecutionSpace, typename SolverType>
    void step( ExecutionSpace exec_space, double& time, const double next_time,
               const double dt, Grid<MemorySpace>& grid, MovingBeam& beam,
               SolverType& fd )
    {
        Kokkos::Profiling::ScopedRegion step_region( "Finch::timestep" );
        time = next_time;

        // update beam position
        beam.move( time );
        double beam_power = beam.power();
        const auto& beam_pos = beam.position();
        const auto& beam_direction = beam.direction();

        // A transient source uses the completed field and current halos to
        // update its explicitly lagged melt-pool-depth closure.
        fd.prepareSource( exec_space, grid, beam );

        // Ping-pong invariant: T0 is the completed previous field and T is
        // the output buffer for this step. T0 is not reused until
        // solidification-event collection is complete.
        grid.swapTemperatureFields();

        // Get temperature views;
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();

        // Solve finite difference
        auto owned_space = grid.getIndexSpace();
        fd.solve( exec_space, owned_space, T, T0, dt, beam_power, beam_pos,
                  beam_direction );

        // update boundaries
        grid.updateBoundaries( fd.materialModel() );

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
