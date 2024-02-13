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

template <typename MemorySpace>
class SolidificationData
{
    using memory_space = MemorySpace;
    using exec_space = typename memory_space::execution_space;
    using view_int = Kokkos::View<int*, memory_space>;
    using view_double2D = Kokkos::View<double**, memory_space>;
    using view_double4D = Kokkos::View<double****, memory_space>;

  private:
    Grid<memory_space>& grid;
    const Inputs& db;

    view_int count;

    int capacity;

    int nCmpts;

    view_double2D events;

    view_double4D tm_view;

    std::string folder_name;

  public:
    // constructor
    SolidificationData( Grid<memory_space>& grid, Inputs& db )
        : grid( grid )
        , db( db )
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

        folder_name = db.sampling.directory_name;
    }

    void updateEvents()
    {
        // get local copies from grid
        auto local_mesh = grid.getLocalMesh();
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();

        // get values from simulation data bases
        auto liquidus = db.properties.liquidus;
        auto dt = db.time.time_step;
        auto time = db.time.time;
        auto cell_size = db.space.cell_size;

        // get local copies of member variables on device for KOKKOS_LAMBDA
        auto local_count = count;
        auto local_capacity = capacity;
        auto local_events = events;
        auto local_tm_view = tm_view;

        using entity_type = typename Grid<memory_space>::entity_type;

        Cabana::Grid::grid_parallel_for(
            "local_grid_for", exec_space(), grid.getIndexSpace(),
            KOKKOS_LAMBDA( const int i, const int j, const int k ) {
                double temp = T( i, j, k, 0 );
                double temp0 = T0( i, j, k, 0 );

                if ( ( temp <= liquidus ) && ( temp0 > liquidus ) )
                {
                    int current_count =
                        Kokkos::atomic_fetch_add( &local_count( 0 ), 1 );

                    if ( current_count < local_capacity )
                    {
                        // event coordinates
                        double pt[3];
                        int idx[3] = { i, j, k };
                        local_mesh.coordinates( entity_type(), idx, pt );
                        local_events( current_count, 0 ) = pt[0];
                        local_events( current_count, 1 ) = pt[1];
                        local_events( current_count, 2 ) = pt[2];

                        // event melting time
                        local_events( current_count, 3 ) =
                            local_tm_view( i, j, k, 0 );

                        // event solidification time
                        double m = ( temp - liquidus ) / ( temp - temp0 );
                        m = fmin( fmax( m, 0.0 ), 1.0 );
                        local_events( current_count, 4 ) = time - m * dt;

                        // cooling rate
                        local_events( current_count, 5 ) =
                            ( temp0 - temp ) / dt;

                        // temperature gradient components
                        local_events( current_count, 6 ) =
                            ( T( i + 1, j, k, 0 ) - T( i - 1, j, k, 0 ) ) /
                            ( 2.0 * cell_size );

                        local_events( current_count, 7 ) =
                            ( T( i, j + 1, k, 0 ) - T( i, j - 1, k, 0 ) ) /
                            ( 2.0 * cell_size );

                        local_events( current_count, 8 ) =
                            ( T( i, j, k + 1, 0 ) - T( i, j, k - 1, 0 ) ) /
                            ( 2.0 * cell_size );
                    }
                }
                else if ( ( temp > liquidus ) && ( temp0 <= liquidus ) )
                {
                    double m = ( temp - liquidus ) / ( temp - temp0 );
                    m = fmin( fmax( m, 0.0 ), 1.0 );
                    local_tm_view( i, j, k, 0 ) = time - m * dt;
                }
            } );
    }

    // Update the solidification data
    void update()
    {
        if ( !db.sampling.enabled )
        {
            return;
        }

        auto count_old_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        updateEvents();

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

            updateEvents();
        }

        // view size is within 90% of capacity. double current size.
        else if ( new_count / capacity > 0.9 )
        {
            capacity = 2.0 * new_count;

            Kokkos::resize( Kokkos::WithoutInitializing, events, capacity,
                            nCmpts );
        }
    }

    // Write the solidification data to separate files for each MPI rank
    void write()
    {
        if ( !db.sampling.enabled )
        {
            return;
        }

        auto events_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), events );

        auto count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );

        // create directory is not present, otherwise overwrite existing files
        if ( mkdir( folder_name.c_str(), 0777 ) != -1 )
        {
            std::cout << "Creating directory: " << folder_name << std::endl;
        }

        std::ofstream fout;
        std::string filename( folder_name + "/data_" +
                              std::to_string( grid.comm_rank ) + ".csv" );

        fout.open( filename );
        fout << std::fixed << std::setprecision( 10 );

        for ( int n = 0; n < count_host( 0 ); n++ )
        {
            fout << events_host( n, 0 ) << "," << events_host( n, 1 ) << ","
                 << events_host( n, 2 ) << "," << events_host( n, 3 ) << ","
                 << events_host( n, 4 ) << "," << events_host( n, 5 );

            if ( db.sampling.format == "default" )
            {
                fout << "," << events_host( n, 6 ) << "," << events_host( n, 7 )
                     << "," << events_host( n, 8 );
            }

            fout << std::endl;
        }

        fout.close();
    }
};

#endif
