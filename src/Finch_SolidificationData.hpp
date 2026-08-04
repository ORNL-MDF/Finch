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

/*!
  \file SolidificationData.hpp
  \brief Class to output solidification information (e.g. for later
  microstructure) simulation
*/

#ifndef SolidificationData_H
#define SolidificationData_H

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <mpi.h>
#include <utility>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_Grid.hpp>
#include <Finch_Inputs.hpp>

namespace Finch
{

template <typename MemorySpace>
class SolidificationData
{
    using memory_space = MemorySpace;
    using exec_space = typename memory_space::execution_space;
    using view_int = Kokkos::View<int*, memory_space>;
    using view_double2D = Kokkos::View<double**, memory_space>;
    using view_double4D = Kokkos::View<double****, memory_space>;
    using view_type_coupled =
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>;

  private:
    // Needed for file output
    int mpi_rank_ = 0;
    double liquidus_ = 0.0;
    double cell_size_ = 0.0;
    bool enabled_ = false;
    SolidificationDataInput input_;

    exec_space exec_space_;

    view_int count;
    Kokkos::View<int*, Kokkos::HostSpace> count_host;

    int count_cached = 0;

    int capacity = 0;

    int nCmpts = 9;

    view_double2D events;

    view_double4D tm_view;

  public:
    // Default constructor
    SolidificationData() = default;
    // constructor
    SolidificationData( const SolidificationDataInput& input,
                        const Inputs& inputs, Grid<memory_space>& grid )
        : mpi_rank_( grid.comm_rank )
        , liquidus_( inputs.properties.liquidus )
        , cell_size_( inputs.space.cell_size )
        , enabled_( input.enabled )
        , input_( input )
        , exec_space_( grid.executionSpace() )
    {
        count = view_int( "count", 1 );
        count_host = Kokkos::View<int*, Kokkos::HostSpace>( "count_host", 1 );
        Kokkos::deep_copy( exec_space_, count, 0 );

        const auto owned_size = grid.getIndexSpace().size();
        if ( owned_size > std::numeric_limits<int>::max() )
            throw std::runtime_error(
                "Local grid is too large for solidification event indexing" );
        // Most domains melt only in a compact region. Start with a bounded
        // fraction of the local grid and grow on demand instead of reserving
        // nine doubles for every local node.
        const auto initial_capacity = std::min<long>(
            owned_size, std::max<long>( 1024, owned_size / 64 ) );
        capacity = std::max<int>( 1, static_cast<int>( initial_capacity ) );

        // components: x, y, z, tm, ts, R, Gx, Gy, Gz
        events =
            view_double2D( Kokkos::ViewAllocateWithoutInitializing( "events" ),
                           capacity, nCmpts );

        auto local_grid = grid.getLocalGrid();
        using entity_type = typename Grid<memory_space>::entity_type;
        auto layout =
            Cabana::Grid::createArrayLayout( local_grid, 1, entity_type() );
        auto tm =
            Cabana::Grid::createArray<double, memory_space>( "tm", layout );
        tm_view = tm->view();
        Kokkos::deep_copy( exec_space_, tm_view, inputs.time.start_time );
    }

    void updateEvents( Grid<memory_space>& grid, const double time,
                       const double dt )
    {
        // get local copies from grid
        auto local_mesh = grid.getLocalMesh();
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();
        auto count_view = count;
        auto events_view = events;
        auto melting_time_view = tm_view;
        const double liquidus = liquidus_;
        const double cell_size = cell_size_;
        const int event_capacity = capacity;

        using entity_type = typename Grid<memory_space>::entity_type;

        Cabana::Grid::grid_parallel_for(
            "Finch::solidification_detection", exec_space_,
            grid.getIndexSpace(),
            KOKKOS_LAMBDA( const int i, const int j, const int k ) {
                double temp = T( i, j, k, 0 );
                double temp0 = T0( i, j, k, 0 );

                if ( ( temp <= liquidus ) && ( temp0 > liquidus ) )
                {
                    int current_count =
                        Kokkos::atomic_fetch_add( &count_view( 0 ), 1 );

                    if ( current_count < event_capacity )
                    {
                        // event coordinates
                        double pt[3];
                        int idx[3] = { i, j, k };
                        local_mesh.coordinates( entity_type(), idx, pt );
                        events_view( current_count, 0 ) = pt[0];
                        events_view( current_count, 1 ) = pt[1];
                        events_view( current_count, 2 ) = pt[2];

                        // event melting time
                        events_view( current_count, 3 ) =
                            melting_time_view( i, j, k, 0 );

                        // event solidification time
                        double m = ( temp - liquidus ) / ( temp - temp0 );
                        m = fmin( fmax( m, 0.0 ), 1.0 );
                        events_view( current_count, 4 ) = time - m * dt;

                        // cooling rate
                        events_view( current_count, 5 ) = ( temp0 - temp ) / dt;

                        // temperature gradient components
                        events_view( current_count, 6 ) =
                            ( T( i + 1, j, k, 0 ) - T( i - 1, j, k, 0 ) ) /
                            ( 2.0 * cell_size );

                        events_view( current_count, 7 ) =
                            ( T( i, j + 1, k, 0 ) - T( i, j - 1, k, 0 ) ) /
                            ( 2.0 * cell_size );

                        events_view( current_count, 8 ) =
                            ( T( i, j, k + 1, 0 ) - T( i, j, k - 1, 0 ) ) /
                            ( 2.0 * cell_size );
                    }
                }
                else if ( ( temp > liquidus ) && ( temp0 <= liquidus ) )
                {
                    double m = ( temp - liquidus ) / ( temp - temp0 );
                    m = fmin( fmax( m, 0.0 ), 1.0 );
                    melting_time_view( i, j, k, 0 ) = time - m * dt;
                }
            } );
    }

    // Update the solidification data
    void update( Grid<memory_space>& grid, const double time, const double dt )
    {
        if ( !enabled_ )
        {
            return;
        }

        Kokkos::Profiling::ScopedRegion region( "Finch::solidification_data" );
        updateEvents( grid, time, dt );

        // This is the sole synchronization needed by event collection in a
        // timestep. Keep the prior count cached on the host so an overflow can
        // be replayed without copying the counter before the kernel.
        Kokkos::deep_copy( exec_space_, count_host, count );
        exec_space_.fence( "Finch solidification count" );

        int new_count = count_host( 0 );

        // more events were added than the current view capacity.
        // resize view and update events starting from the previous counter.
        if ( new_count > capacity )
        {
            capacity = std::max( 2 * capacity, new_count );

            Kokkos::resize( Kokkos::WithoutInitializing, events, capacity,
                            nCmpts );

            Kokkos::deep_copy( exec_space_, count, count_cached );

            updateEvents( grid, time, dt );
        }

        // view size is within 90% of capacity. double current size.
        else if ( static_cast<double>( new_count ) /
                      static_cast<double>( capacity ) >
                  0.9 )
        {
            capacity = std::max( 2 * capacity, new_count + 1 );

            Kokkos::resize( Kokkos::WithoutInitializing, events, capacity,
                            nCmpts );
        }
        count_cached = new_count;
    }

    // Return all data for the events that have been recorded during the
    // simulation
    auto get()
    {
        if ( !enabled_ )
            return view_type_coupled(
                Kokkos::ViewAllocateWithoutInitializing( "copied_data" ), 0,
                nCmpts );

        auto valid_events = Kokkos::subview(
            events, std::make_pair( 0, count_cached ), Kokkos::ALL );
        auto events_host = Kokkos::create_mirror_view( valid_events );
        Kokkos::deep_copy( exec_space_, events_host, valid_events );
        exec_space_.fence( "Finch copy solidification events" );
        // Create a View on the host with fixed layout for coupling.
        view_type_coupled copied_data(
            Kokkos::ViewAllocateWithoutInitializing( "copied_data" ),
            count_cached, nCmpts );
        Kokkos::deep_copy( copied_data, events_host );
        return copied_data;
    }

    // Write the solidification data to separate files for each MPI rank
    void write( MPI_Comm comm )
    {
        if ( !enabled_ )
        {
            return;
        }

        std::chrono::high_resolution_clock::time_point
            start_solidification_print_time =
                std::chrono::high_resolution_clock::now();
        auto valid_events = Kokkos::subview(
            events, std::make_pair( 0, count_cached ), Kokkos::ALL );
        auto events_host = Kokkos::create_mirror_view( valid_events );
        Kokkos::deep_copy( exec_space_, events_host, valid_events );
        exec_space_.fence( "Finch copy solidification events for output" );

        // Create a shared output directory once, before any rank opens files.
        int directory_error = 0;
        if ( mpi_rank_ == 0 )
        {
            std::error_code ec;
            std::filesystem::create_directories( input_.directory, ec );
            if ( ec )
            {
                std::cerr << "Cannot create solidification directory "
                          << input_.directory << ": " << ec.message()
                          << std::endl;
                directory_error = 1;
            }
        }
        MPI_Bcast( &directory_error, 1, MPI_INT, 0, comm );
        if ( directory_error )
            throw std::runtime_error(
                "Unable to create solidification output directory" );
        MPI_Barrier( comm );

        std::ofstream fout;
        std::string filename( input_.directory + "/data_" +
                              std::to_string( mpi_rank_ ) + ".csv" );

        fout.open( filename );
        int local_file_error = !fout;
        int global_file_error = 0;
        MPI_Allreduce( &local_file_error, &global_file_error, 1, MPI_INT,
                       MPI_MAX, comm );
        if ( global_file_error )
            throw std::runtime_error( "Cannot open solidification output " +
                                      filename );
        fout << std::fixed << std::setprecision( 10 );

        for ( int n = 0; n < count_cached; n++ )
        {
            fout << events_host( n, 0 ) << "," << events_host( n, 1 ) << ","
                 << events_host( n, 2 ) << "," << events_host( n, 3 ) << ","
                 << events_host( n, 4 ) << "," << events_host( n, 5 );

            if ( input_.format == "default" )
            {
                fout << "," << events_host( n, 6 ) << "," << events_host( n, 7 )
                     << "," << events_host( n, 8 );
            }

            fout << '\n';
        }

        fout.close();

        MPI_Barrier( comm );
        std::chrono::high_resolution_clock::time_point
            end_solidification_print_time =
                std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds =
            end_solidification_print_time - start_solidification_print_time;
        if ( mpi_rank_ == 0 )
            std::cout << "Solidification data written in " << std::fixed
                      << std::setprecision( 6 ) << elapsed_seconds.count()
                      << " seconds" << std::endl;
    }

    std::array<double, 3> getLowerBounds( MPI_Comm comm )
    {
        if ( !enabled_ )
        {
            std::array<double, 3> bounds;
            bounds.fill( std::numeric_limits<double>::quiet_NaN() );
            return bounds;
        }

        auto data = get();
        std::array<double, 3> local = {
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity() };
        for ( int n = 0; n < count_cached; ++n )
            for ( int d = 0; d < 3; ++d )
                local[d] = std::min( local[d], data( n, d ) );

        int global_count = 0;
        MPI_Allreduce( &count_cached, &global_count, 1, MPI_INT, MPI_SUM,
                       comm );
        std::array<double, 3> data_bounds_low;
        MPI_Allreduce( local.data(), data_bounds_low.data(), 3, MPI_DOUBLE,
                       MPI_MIN, comm );

        if ( global_count == 0 )
        {
            data_bounds_low.fill( std::numeric_limits<double>::quiet_NaN() );
            if ( mpi_rank_ == 0 )
                std::cout << "No melted/resolidified region was recorded."
                          << std::endl;
            return data_bounds_low;
        }

        if ( mpi_rank_ == 0 )
        {
            std::cout << "Min X bound of the melted/resolidified region was "
                      << data_bounds_low[0] << std::endl;
            std::cout << "Min Y bound of the melted/resolidified region was "
                      << data_bounds_low[1] << std::endl;
            std::cout << "Min Z bound of the melted/resolidified region was "
                      << data_bounds_low[2] << std::endl;
        }
        return data_bounds_low;
    }

    std::array<double, 3> getUpperBounds( MPI_Comm comm )
    {
        if ( !enabled_ )
        {
            std::array<double, 3> bounds;
            bounds.fill( std::numeric_limits<double>::quiet_NaN() );
            return bounds;
        }

        auto data = get();
        std::array<double, 3> local = {
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity() };
        for ( int n = 0; n < count_cached; ++n )
            for ( int d = 0; d < 3; ++d )
                local[d] = std::max( local[d], data( n, d ) );

        int global_count = 0;
        MPI_Allreduce( &count_cached, &global_count, 1, MPI_INT, MPI_SUM,
                       comm );
        std::array<double, 3> data_bounds_high;
        MPI_Allreduce( local.data(), data_bounds_high.data(), 3, MPI_DOUBLE,
                       MPI_MAX, comm );

        if ( global_count == 0 )
        {
            data_bounds_high.fill( std::numeric_limits<double>::quiet_NaN() );
            return data_bounds_high;
        }

        if ( mpi_rank_ == 0 )
        {
            std::cout << "Max X bound of the melted/resolidified region was "
                      << data_bounds_high[0] << std::endl;
            std::cout << "Max Y bound of the melted/resolidified region was "
                      << data_bounds_high[1] << std::endl;
            std::cout << "Max Z bound of the melted/resolidified region was "
                      << data_bounds_high[2] << std::endl;
        }
        return data_bounds_high;
    }
};

} // namespace Finch

#endif
