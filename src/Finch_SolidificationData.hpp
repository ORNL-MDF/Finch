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
#include <chrono>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#ifdef Finch_ENABLE_STORK
#include <Stork_Core.hpp>
#endif

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
    using view_double2D_host = typename view_double2D::host_mirror_type;
    using view_double4D = Kokkos::View<double****, memory_space>;
    using view_type_coupled =
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
#ifdef Finch_ENABLE_STORK
    using DualSRDF = Stork::Structs::SRDF_Dual<double>;
    using DualRDF = Stork::Structs::RDF_Dual<double>;
#endif
  private:
    // Needed for file output
    int mpi_rank_;
    double liquidus_;
    double dt_;
    double cell_size_;
    bool enabled_;
    std::string format_;
    int fine_factor_;
    bool srdf_format_;
    view_int count;

    int capacity;

    double x_min_sampling, y_min_sampling, z_min_sampling, x_max_sampling,
        y_max_sampling, z_max_sampling;

    // Used for SRDF-based data collection
    int nx_solidification, ny_solidification, nz_solidification;
    view_int cellnum;
    view_double2D timesview, thermalsview;
    // Used for RDF-based data collection
    view_double4D tm_view;
    view_double2D events;
    int nCmpts;

  public:
    // Default constructor
    SolidificationData() {}
    // constructor
    SolidificationData( const Inputs& inputs, Grid<memory_space>& grid,
                        const bool srdf_format )
        : mpi_rank_( grid.comm_rank )
        , liquidus_( inputs.properties.liquidus )
        , dt_( inputs.time.time_step )
        , cell_size_( inputs.space.cell_size )
        , enabled_( inputs.sampling.enabled )
        , format_( inputs.sampling.format )
        , fine_factor_( inputs.sampling.fine_factor )
        , srdf_format_( srdf_format )
    {
        count = view_int( "count", 1 );

        capacity = round( grid.getIndexSpace().size() );

        auto local_mesh = grid.getLocalMesh();
        auto local_grid = grid.getLocalGrid();
        using entity_type = typename Grid<memory_space>::entity_type;
        auto layout =
            Cabana::Grid::createArrayLayout( local_grid, 1, entity_type() );

        if ( srdf_format_ )
        {
            // Allocate views for SRDF data
            cellnum =
                view_int( Kokkos::ViewAllocateWithoutInitializing( "cellnum" ),
                          capacity );
            // 2 components: previous, current time step
            timesview = view_double2D(
                Kokkos::ViewAllocateWithoutInitializing( "times" ), capacity,
                2 );
            thermalsview = view_double2D(
                Kokkos::ViewAllocateWithoutInitializing( "thermals" ), capacity,
                16 );
            nCmpts = 6;
        }
        else
        {
            // Allocate views for RDF data
            auto tm =
                Cabana::Grid::createArray<double, memory_space>( "tm", layout );
            tm_view = tm->view();
            Kokkos::deep_copy( tm_view, 0 );
            events = view_double2D(
                Kokkos::ViewAllocateWithoutInitializing( "events" ), capacity,
                9 );
            nCmpts = 9;
        }

        // Lower bounds of region considered for sampling, plus round off buffer
        x_min_sampling = inputs.sampling.global_low_corner[0];
        y_min_sampling = inputs.sampling.global_low_corner[1];
        z_min_sampling = inputs.sampling.global_low_corner[2];

        // CA grid - store node data in halo regions in positive x,y,z, but if
        // at a global domain boundary, do not store boundary node data
        if ( std::abs( local_mesh.highCorner( Cabana::Grid::Own(), 0 ) -
                       inputs.space.global_high_corner[0] ) < 1e-10 )
            nx_solidification = grid.num_points_x;
        else
            nx_solidification = grid.num_points_x + 1;
        if ( std::abs( local_mesh.highCorner( Cabana::Grid::Own(), 1 ) -
                       inputs.space.global_high_corner[1] ) < 1e-10 )
            ny_solidification = grid.num_points_y;
        else
            ny_solidification = grid.num_points_y + 1;
        if ( std::abs( local_mesh.highCorner( Cabana::Grid::Own(), 2 ) -
                       inputs.space.global_high_corner[2] ) < 1e-10 )
            nz_solidification = grid.num_points_z;
        else
            nz_solidification = grid.num_points_z + 1;

        // Upper bounds of region considered for sampling (1e-10 to avoid
        // floating point error in comparisons)
        x_max_sampling = std::min( inputs.space.global_high_corner[0],
                                   inputs.sampling.global_high_corner[0] ) -
                         1e-10;
        y_max_sampling = std::min( inputs.space.global_high_corner[1],
                                   inputs.sampling.global_high_corner[1] ) -
                         1e-10;
        z_max_sampling = std::min( inputs.space.global_high_corner[2],
                                   inputs.sampling.global_high_corner[2] ) -
                         1e-10;
    }

    // Update the solidification data in the reduced data format
    void updateRDFEvents( Grid<memory_space>& grid, const double time )
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

    // Update the solidification data in the stork-compatible data format
    void updateSRDFEvents( Grid<memory_space>& grid, const double time )
    {
        // get local copies from grid
        auto T = grid.getTemperature();
        auto T0 = grid.getPreviousTemperature();
        auto local_mesh = grid.getLocalMesh();
        using entity_type = typename Grid<memory_space>::entity_type;
        double x_min_sampling_ = x_min_sampling;
        double y_min_sampling_ = y_min_sampling;
        double z_min_sampling_ = z_min_sampling;
        double x_max_sampling_ = x_max_sampling;
        double y_max_sampling_ = y_max_sampling;
        double z_max_sampling_ = z_max_sampling;
        double dt = dt_;
        int capacity_ = capacity;
        int ny_solidification_ = ny_solidification;
        int nz_solidification_ = nz_solidification;
        Cabana::Grid::grid_parallel_for(
            "local_grid_for", exec_space(), grid.getIndexSpace(),
            KOKKOS_CLASS_LAMBDA( const int i, const int j, const int k ) {
                double pt[3];
                int idx[3] = { i, j, k };
                local_mesh.coordinates( entity_type(), idx, pt );
                // Loop over owned points, except for global x,y,z bound
                if ( ( pt[0] >= x_min_sampling_ ) &&
                     ( pt[1] >= y_min_sampling_ ) &&
                     ( pt[2] >= z_min_sampling_ ) &&
                     ( pt[0] < x_max_sampling_ ) &&
                     ( pt[1] < y_max_sampling_ ) &&
                     ( pt[2] < z_max_sampling_ ) )
                {
                    // Count number of vertices above the liquidus on this time
                    // step
                    int vert_above_liquidus = 0;
                    int vert_above_liquidus_old = 0;
                    for ( int n_index = 0; n_index < 8; ++n_index )
                    {
                        const int neighbor_zn = k + ( n_index & 1 );
                        const int neighbor_yn = j + ( ( n_index >> 1 ) & 1 );
                        const int neighbor_xn = i + ( ( n_index >> 2 ) & 1 );

                        vert_above_liquidus_old +=
                            T0( neighbor_xn, neighbor_yn, neighbor_zn, 0 ) >=
                            liquidus_;
                        vert_above_liquidus += T( neighbor_xn, neighbor_yn,
                                                  neighbor_zn, 0 ) >= liquidus_;
                    }
                    // store previous, current temperature state if:
                    // - between 1 and 7 of the vertices were above the liquidus
                    // on either the previous or current time step
                    // - all vertices were above the liquidus and now all
                    // vertices are below the liquidus
                    // - all vertices were below the liquidus and now all
                    // vertices are above the liquidus
                    if ( ( vert_above_liquidus_old % 8 ) ||
                         ( vert_above_liquidus % 8 ) ||
                         ( vert_above_liquidus != vert_above_liquidus_old ) )
                    {
                        auto counter =
                            Kokkos::atomic_fetch_add( &count( 0 ), 1 );
                        // Store 1D index of cell in a way consistent with Stork
                        // - if there's space in the structs
                        if ( counter < capacity_ )
                        {
                            cellnum( counter ) =
                                ( i - 1 ) * ny_solidification_ *
                                    nz_solidification_ +
                                ( j - 1 ) * nz_solidification_ + ( k - 1 );

                            // Store previous, current times
                            timesview( counter, 0 ) = time - dt;
                            timesview( counter, 1 ) = time;

                            // Store vertex temperatures at previous, current
                            // times
                            for ( int n_index = 0; n_index < 8; ++n_index )
                            {
                                const int neighbor_zn = k + ( n_index & 1 );
                                const int neighbor_yn =
                                    j + ( ( n_index >> 1 ) & 1 );
                                const int neighbor_xn =
                                    i + ( ( n_index >> 2 ) & 1 );
                                thermalsview( counter, n_index ) = T0(
                                    neighbor_xn, neighbor_yn, neighbor_zn, 0 );
                                thermalsview( counter, n_index + 8 ) = T(
                                    neighbor_xn, neighbor_yn, neighbor_zn, 0 );
                            }
                        }
                    }
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

        // Get the number of events or snapshots
        auto count_old_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );
        int count_old = count_old_host( 0 );

        // Update the events or snapshots and get the new total number
        if ( srdf_format_ )
            updateSRDFEvents( grid, time );
        else
            updateRDFEvents( grid, time );
        auto count_new_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );
        int new_count = count_new_host( 0 );

        // more events were added than the current view capacity.
        // resize view and update events starting from the previous counters
        if ( new_count >= capacity )
        {
            capacity = 2.0 * new_count;
            Kokkos::deep_copy( count, count_old );
            if ( srdf_format_ )
            {
                Kokkos::resize( Kokkos::WithoutInitializing, cellnum,
                                capacity );
                Kokkos::resize( Kokkos::WithoutInitializing, timesview,
                                capacity, 2 );
                Kokkos::resize( Kokkos::WithoutInitializing, thermalsview,
                                capacity, 16 );
                updateSRDFEvents( grid, time );
            }
            else
            {
                Kokkos::resize( Kokkos::WithoutInitializing, events, capacity,
                                9 );
                updateRDFEvents( grid, time );
            }
        }
        // view size is within 90% of capacity. double current size.
        else if ( new_count / capacity > 0.9 )
        {
            capacity = 2.0 * new_count;
            if ( srdf_format_ )
            {
                Kokkos::resize( Kokkos::WithoutInitializing, cellnum,
                                capacity );
                Kokkos::resize( Kokkos::WithoutInitializing, timesview,
                                capacity, 2 );
                Kokkos::resize( Kokkos::WithoutInitializing, thermalsview,
                                capacity, 16 );
            }
            else
                Kokkos::resize( Kokkos::WithoutInitializing, events, capacity,
                                9 );
        }
    }

    // Without stork, copy the existing device view to the host
    view_double2D_host getEvents()
    {
        view_double2D_host events_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), events );
        auto count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );
        // Resize the host copy so only valid events get copied.
        Kokkos::resize( events_host, count_host( 0 ), nCmpts );
        return events_host;
    }

#ifdef Finch_ENABLE_STORK
    // Using Stork, interpolate the temperature vertex data to obtain the
    // solidification events
    view_double2D_host getEvents( Grid<memory_space>& grid, MPI_Comm comm )
    {
        view_double2D_host events_host(
            Kokkos::ViewAllocateWithoutInitializing( "events_host" ), 1, 6 );
        auto local_mesh = grid.getLocalMesh();
        DualSRDF SRDF;
        Stork::Structs::RegularGrid_Header<double, Stork::host_space> header =
            SRDF.host_header;

        // Origin point
        header.global_x0() = local_mesh.lowCorner( Cabana::Grid::Own(), 0 );
        header.global_y0() = local_mesh.lowCorner( Cabana::Grid::Own(), 1 );
        header.global_z0() = local_mesh.lowCorner( Cabana::Grid::Own(), 2 );
        std::cout << "Rank " << mpi_rank_ << " low corner is at "
                  << header.global_x0() << ", " << header.global_y0() << ", "
                  << header.global_z0() << std::endl;

        // Each MPI rank interpolates only on the local grid - no awareness of
        // global grid needed
        header.global_i0() = 0;
        header.global_j0() = 0;
        header.global_k0() = 0;

        // Number of points in each direction - includes halo in positive
        // directions unless at the global domain bound
        header.local_inum() = nx_solidification;
        header.local_jnum() = ny_solidification;
        header.local_knum() = nz_solidification;

        // Cell size
        header.gridResolution() = cell_size_;

        // Liquidus temperature
        SRDF.T_critical = liquidus_;

        // Resize views based on final count of temperature snapshots
        auto snapshot_count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );
        Kokkos::resize( cellnum, snapshot_count_host( 0 ) );
        Kokkos::resize( timesview, snapshot_count_host( 0 ), 2 );
        Kokkos::resize( thermalsview, snapshot_count_host( 0 ), 16 );
        // Number of snapshots
        SRDF.numSnaps = snapshot_count_host( 0 );
        std::cout << "Rank " << mpi_rank_
                  << " number of snapshots: " << snapshot_count_host( 0 )
                  << std::endl;
        // Copy views to host
        auto cellnum_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), cellnum );
        auto timesview_host = Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace(), timesview );
        auto thermalsview_host = Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace(), thermalsview );

        // Create tuples to order by cell number
        std::vector<int> old_list_position( snapshot_count_host( 0 ) );
        for ( int n = 0; n < snapshot_count_host( 0 ); n++ )
            old_list_position[n] = n;
        std::vector<std::tuple<int, int>> snapshots;
        snapshots.reserve( snapshot_count_host( 0 ) );

        for ( int n = 0; n < snapshot_count_host( 0 ); n++ )
        {
            snapshots.push_back(
                std::make_tuple( cellnum_host( n ), old_list_position[n] ) );
        }
        // Sorting from low to high
        std::sort( snapshots.begin(), snapshots.end() );

        // Create empty SRDF views based on count
        SRDF.template Make_Data_Views<Stork::device_space>(
            snapshot_count_host( 0 ) );
        // Make data mirrors on host
        SRDF.template Make_Data_Mirrors<Stork::device_space,
                                        Stork::host_space>();

        // Get reference to SRDF views on host
        Stork::Structs::SRDF_Data<double, Stork::host_space>& data =
            SRDF.host_data;

        // Fill SRDF views from sorted Finch data
        for ( int n = 0; n < snapshot_count_host( 0 ); n++ )
        {
            data.cellNum_view( n ) = std::get<0>( snapshots[n] );
            const int old_list_pos = std::get<1>( snapshots[n] );
            data.times_view( 2 * n ) = timesview_host( old_list_pos, 0 );
            data.times_view( 2 * n + 1 ) = timesview_host( old_list_pos, 1 );
            for ( int vert = 0; vert < 16; vert++ )
            {
                data.thermals_view( 16 * n + vert ) =
                    thermalsview_host( old_list_pos, vert );
            }
        }

        // Interpolate to fine cell size
        DualRDF RDF =
            Stork::Run::Interpolate_SRDF_to_RDF<double, Stork::host_space,
                                                double, Stork::device_space>(
                SRDF, fine_factor_ );

        // Copy data from device back to host
        RDF.template Make_Data_Mirrors<Stork::device_space,
                                       Stork::host_space>();
        RDF.template Copy_All<Stork::device_space, Stork::host_space>();
        std::cout << "Rank " << mpi_rank_ << " new domain "
                  << RDF.host_header.local_inum() << ","
                  << RDF.host_header.local_jnum() << ","
                  << RDF.host_header.local_knum() << std::endl;

        // Output Data to File (don't use Stork::IO::Output_RDF_csv since that
        // will put the x,y,z,tm,tl,cr in all files)

        // Get reduced data format header and data for output
        Stork::Structs::RegularGrid_Header<double, Stork::host_space>&
            output_header = RDF.host_header;
        Stork::Structs::RDF_Data<double, Stork::host_space>& output_data =
            RDF.host_data;

        // Get total number of events
        const size_t numEvents = RDF.numEvents;
        // With the number of events now known, allocate events view on host
        Kokkos::realloc( events_host, numEvents, 6 );
        std::cout << "Rank " << mpi_rank_ << " storing " << numEvents
                  << " interpolated solidification events" << std::endl;
        for ( uint32_t n = 0; n < numEvents; n++ )
        {
            // Get global xyz from local p
            const uint32_t& local_p = output_data.p( n );
            double global_xyz[3];
            output_header.LOCAL_p_to_GLOBAL_xyz( global_xyz, local_p );
            for ( int comp = 0; comp < 3; comp++ )
                events_host( n, comp ) = global_xyz[comp];
            events_host( n, 3 ) = output_data.tm( n );
            events_host( n, 4 ) = output_data.tl( n );
            events_host( n, 5 ) = output_data.cr( n );
        }
        // Copy event data, event count to device
        Kokkos::realloc( events, numEvents, 6 );
        events =
            Kokkos::create_mirror_view_and_copy( memory_space(), events_host );
        auto event_count_host =
            Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), count );
        event_count_host( 0 ) = numEvents;
        count = Kokkos::create_mirror_view_and_copy( memory_space(),
                                                     event_count_host );
        MPI_Barrier( comm );
        return events_host;
    }
#endif

    // Return all data for the events that have been recorded during the
    // simulation
    auto get( [[maybe_unused]] Grid<memory_space>& grid,
              [[maybe_unused]] MPI_Comm comm, Sampling sampling_inputs,
              bool write_data )
    {
        view_double2D_host events_host;
#ifdef Finch_ENABLE_STORK
        if ( srdf_format_ )
            events_host = getEvents( grid, comm );
        else
            events_host = getEvents();
#else
        events_host = getEvents();
#endif
        const int num_events = events_host.extent( 0 );
        // Create a View on the host with fixed layout for coupling.
        view_type_coupled copied_data(
            Kokkos::ViewAllocateWithoutInitializing( "copied_data" ),
            num_events, nCmpts );
        Kokkos::deep_copy( copied_data, events_host );
        if ( write_data )
            write( copied_data, sampling_inputs, comm );
        return copied_data;
    }

    // Write the solidification data to separate files for each MPI rank
    void write( view_type_coupled& solidification_data,
                Sampling sampling_inputs, MPI_Comm comm )
    {
        if ( !enabled_ )
        {
            return;
        }

        std::chrono::high_resolution_clock::time_point
            start_solidification_print_time =
                std::chrono::high_resolution_clock::now();

        // create directory is not present, otherwise overwrite existing files
        if ( mkdir( sampling_inputs.directory_name.c_str(), 0777 ) != -1 )
        {
            std::cout << "Creating directory: "
                      << sampling_inputs.directory_name << std::endl;
        }

        std::string filename( sampling_inputs.directory_name + "/data_" +
                              std::to_string( mpi_rank_ ) + ".csv" );
        std::ofstream fout;
        fout.open( filename );
        fout << std::fixed << std::setprecision( 10 );
        const int num_events = solidification_data.extent( 0 );
        for ( int n = 0; n < num_events; n++ )
        {
            fout << solidification_data( n, 0 ) << ","
                 << solidification_data( n, 1 ) << ","
                 << solidification_data( n, 2 ) << ","
                 << solidification_data( n, 3 ) << ","
                 << solidification_data( n, 4 ) << ","
                 << solidification_data( n, 5 );

            if ( sampling_inputs.format == "default" )
            {
                fout << "," << solidification_data( n, 6 ) << ","
                     << solidification_data( n, 7 ) << ","
                     << solidification_data( n, 8 );
            }

            fout << std::endl;
        }

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
