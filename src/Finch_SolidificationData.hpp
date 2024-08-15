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

#include <array>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>

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
    int mpi_rank_;
    std::string folder_name_;
    double liquidus_;
    double dt_;
    double cell_size_;
    bool enabled_;
    std::string format_;

    view_int count;

    int capacity;

    int nCmpts;

    view_double2D events;

    view_double4D tm_view;

  public:
    // Default constructor
    SolidificationData() {}
    // constructor
    SolidificationData( const Inputs& inputs, Grid<memory_space>& grid )
        : mpi_rank_( grid.comm_rank )
        , folder_name_( inputs.sampling.directory_name )
        , liquidus_( inputs.properties.liquidus )
        , dt_( inputs.time.time_step )
        , cell_size_( inputs.space.cell_size )
        , enabled_( inputs.sampling.enabled )
        , format_( inputs.sampling.format )
    {
        count = view_int( "count", 1 );

        capacity = round( grid.getIndexSpace().size() );

        // components: x, y, z, tm, ts, R, Gx, Gy, Gz
        nCmpts = 9;

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
    }

    void updateEvents( Grid<memory_space>& grid, const double time )
    {
        // get local copies from grid
        auto local_mesh = grid.getLocalMesh();
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();

        using entity_type = typename Grid<memory_space>::entity_type;

        Cabana::Grid::grid_parallel_for(
            "local_grid_for", exec_space(), grid.getIndexSpace(),
            KOKKOS_CLASS_LAMBDA( const int i, const int j, const int k ) {
                double temp = T( i, j, k, 0 );
                double temp0 = T0( i, j, k, 0 );

                if ( ( temp <= liquidus_ ) && ( temp0 > liquidus_ ) )
                {
                    int current_count =
                        Kokkos::atomic_fetch_add( &count( 0 ), 1 );

                    if ( current_count < capacity )
                    {
                        // event coordinates
                        double pt[3];
                        int idx[3] = { i, j, k };
                        local_mesh.coordinates( entity_type(), idx, pt );
                        events( current_count, 0 ) = pt[0];
                        events( current_count, 1 ) = pt[1];
                        events( current_count, 2 ) = pt[2];

                        // event melting time
                        events( current_count, 3 ) = tm_view( i, j, k, 0 );

                        // event solidification time
                        double m = ( temp - liquidus_ ) / ( temp - temp0 );
                        m = fmin( fmax( m, 0.0 ), 1.0 );
                        events( current_count, 4 ) = time - m * dt_;

                        // cooling rate
                        events( current_count, 5 ) = ( temp0 - temp ) / dt_;

                        // temperature gradient components
                        events( current_count, 6 ) =
                            ( T( i + 1, j, k, 0 ) - T( i - 1, j, k, 0 ) ) /
                            ( 2.0 * cell_size_ );

                        events( current_count, 7 ) =
                            ( T( i, j + 1, k, 0 ) - T( i, j - 1, k, 0 ) ) /
                            ( 2.0 * cell_size_ );

                        events( current_count, 8 ) =
                            ( T( i, j, k + 1, 0 ) - T( i, j, k - 1, 0 ) ) /
                            ( 2.0 * cell_size_ );
                    }
                }
                else if ( ( temp > liquidus_ ) && ( temp0 <= liquidus_ ) )
                {
                    double m = ( temp - liquidus_ ) / ( temp - temp0 );
                    m = fmin( fmax( m, 0.0 ), 1.0 );
                    tm_view( i, j, k, 0 ) = time - m * dt_;
                }
            } );
    }

    // Update the solidification data
    void update( Grid<memory_space>& grid, const double time )
    {
        if ( !enabled_ )
        {
            return;
        }

        auto count_old_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        updateEvents( grid, time );

        auto count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        int new_count = count_host( 0 );

        // more events were added than the current view capacity.
        // resize view and update events starting from the previous counter.
        if ( new_count >= capacity )
        {
            capacity = 2.0 * new_count;

            Kokkos::resize( Kokkos::WithoutInitializing, events, capacity,
                            nCmpts );

            Kokkos::deep_copy( count, count_old_host( 0 ) );

            updateEvents( grid, time );
        }

        // view size is within 90% of capacity. double current size.
        else if ( new_count / capacity > 0.9 )
        {
            capacity = 2.0 * new_count;

            Kokkos::resize( Kokkos::WithoutInitializing, events, capacity,
                            nCmpts );
        }
    }

    // Return all data for the events that have been recorded during the
    // simulation
    auto get()
    {
        auto events_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), events );
        auto count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );
        // Resize the host copy so only valid events get copied.
        Kokkos::resize( events_host, count_host( 0 ), nCmpts );
        // Create a View on the host with fixed layout for coupling.
        view_type_coupled copied_data(
            Kokkos::ViewAllocateWithoutInitializing( "copied_data" ),
            count_host( 0 ), nCmpts );
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

        double start_solidification_print_time = MPI_Wtime();
        auto events_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), events );

        auto count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        // create directory is not present, otherwise overwrite existing files
        if ( mkdir( folder_name_.c_str(), 0777 ) != -1 )
        {
            std::cout << "Creating directory: " << folder_name_ << std::endl;
        }

        std::ofstream fout;
        std::string filename( folder_name_ + "/data_" +
                              std::to_string( mpi_rank_ ) + ".csv" );

        fout.open( filename );
        fout << std::fixed << std::setprecision( 10 );

        for ( int n = 0; n < count_host( 0 ); n++ )
        {
            fout << events_host( n, 0 ) << "," << events_host( n, 1 ) << ","
                 << events_host( n, 2 ) << "," << events_host( n, 3 ) << ","
                 << events_host( n, 4 ) << "," << events_host( n, 5 );

            if ( format_ == "default" )
            {
                fout << "," << events_host( n, 6 ) << "," << events_host( n, 7 )
                     << "," << events_host( n, 8 );
            }

            fout << std::endl;
        }

        fout.close();

        MPI_Barrier( comm );
        double end_solidification_print_time = MPI_Wtime();
        if ( mpi_rank_ == 0 )
            std::cout << "Solidification data written in "
                      << end_solidification_print_time -
                             start_solidification_print_time
                      << " seconds" << std::endl;
    }

    std::array<double, 3> getLowerBounds( MPI_Comm comm )
    {
        // Local copies for lambda capture
        auto events_ = events;
        auto count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        // Iterate over list of events, getting the min bounds in each direction
        double x_min, y_min, z_min;
        Kokkos::parallel_reduce(
            "solidification_event_bounds", count_host( 0 ),
            KOKKOS_LAMBDA( const int& n, double& x_min_th, double& y_min_th,
                           double& z_min_th ) {
                double x_event = events_( n, 0 );
                double y_event = events_( n, 1 );
                double z_event = events_( n, 2 );
                if ( x_event < x_min_th )
                    x_min_th = x_event;
                if ( y_event < y_min_th )
                    y_min_th = y_event;
                if ( z_event < z_min_th )
                    z_min_th = z_event;
            },
            Kokkos::Min<double>( x_min ), Kokkos::Min<double>( y_min ),
            Kokkos::Min<double>( z_min ) );

        // Get the min bounds on each direction across all ranks
        std::array<double, 3> data_bounds_low;
        MPI_Allreduce( &x_min, &data_bounds_low[0], 1, MPI_DOUBLE, MPI_MIN,
                       comm );
        MPI_Allreduce( &y_min, &data_bounds_low[1], 1, MPI_DOUBLE, MPI_MIN,
                       comm );
        MPI_Allreduce( &z_min, &data_bounds_low[2], 1, MPI_DOUBLE, MPI_MIN,
                       comm );

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
        // Local copies for lambda capture
        auto events_ = events;
        auto count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        // Iterate over list of events, getting the max bounds in each direction
        double x_max, y_max, z_max;
        Kokkos::parallel_reduce(
            "solidification_event_bounds", count_host( 0 ),
            KOKKOS_LAMBDA( const int& n, double& x_max_th, double& y_max_th,
                           double& z_max_th ) {
                double x_event = events_( n, 0 );
                double y_event = events_( n, 1 );
                double z_event = events_( n, 2 );
                if ( x_event > x_max_th )
                    x_max_th = x_event;
                if ( y_event > y_max_th )
                    y_max_th = y_event;
                if ( z_event > z_max_th )
                    z_max_th = z_event;
            },
            Kokkos::Max<double>( x_max ), Kokkos::Max<double>( y_max ),
            Kokkos::Max<double>( z_max ) );

        // Get the min/max bounds on each direction across all ranks
        std::array<double, 3> data_bounds_high;
        MPI_Allreduce( &x_max, &data_bounds_high[0], 1, MPI_DOUBLE, MPI_MAX,
                       comm );
        MPI_Allreduce( &y_max, &data_bounds_high[1], 1, MPI_DOUBLE, MPI_MAX,
                       comm );
        MPI_Allreduce( &z_max, &data_bounds_high[2], 1, MPI_DOUBLE, MPI_MAX,
                       comm );

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
