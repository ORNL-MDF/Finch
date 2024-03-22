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

  private:
    // Needed for file output
    int mpi_rank_;
    std::string folder_name_;
    double liquidus_;
    double dt_;
    double time_;
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
        , time_( inputs.time.time )
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

    void updateEvents( Grid<memory_space>& grid )
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
                        events( current_count, 4 ) = time_ - m * dt_;

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
                    tm_view( i, j, k, 0 ) = time_ - m * dt_;
                }
            } );
    }

    // Update the solidification data
    void update( Grid<memory_space>& grid )
    {
        if ( !enabled_ )
        {
            return;
        }

        auto count_old_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        updateEvents( grid );

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

            updateEvents( grid );
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
        auto event_list = std::make_pair( 0, count_host( 0 ) );
        auto temp_components = std::make_pair( 0, nCmpts );
        auto input_temperature_data_host =
            Kokkos::subview( events_host, event_list, temp_components );
        return input_temperature_data_host;
    }

    // Write the solidification data to separate files for each MPI rank
    void write()
    {
        if ( !enabled_ )
        {
            return;
        }

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
    }
};

} // namespace Finch

#endif
