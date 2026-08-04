/****************************************************************************
 * Copyright (c) 2026 by Oak Ridge National Laboratory                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of Finch. Finch is distributed under a                 *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

/*!
  \file Finch_MeltPoolDimensions.hpp
  \brief Single-pass solidus and liquidus melt-pool dimension diagnostics.
*/

#ifndef FINCH_MELT_POOL_DIMENSIONS_HPP
#define FINCH_MELT_POOL_DIMENSIONS_HPP

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <Finch_Grid.hpp>
#include <Finch_Inputs.hpp>
#include <Finch_Isotherm.hpp>

namespace Finch
{

// Store minima and negated maxima so all twelve values use a single minimum
// reduction, including the subsequent MPI collective.
struct MeltPoolBoundsValue
{
    double value[12];
};

template <typename TemperatureView>
struct MeltPoolBoundsFunctor
{
    using value_type = MeltPoolBoundsValue;

    TemperatureView temperature;
    double solidus = 0.0;
    double liquidus = 0.0;
    double spacing = 0.0;
    double origin[3] = { 0.0, 0.0, 0.0 };
    double direction_x = 1.0;
    double direction_y = 0.0;
    double minimum_threshold = 0.0;
    int scan_path_frame = 0;
    int enabled[2] = { 1, 1 };
    int owned_high[3] = { 0, 0, 0 };
    int on_high_boundary[3] = { 0, 0, 0 };

    KOKKOS_INLINE_FUNCTION
    void init( value_type& bounds ) const
    {
        for ( int n = 0; n < 12; ++n )
            bounds.value[n] = Kokkos::reduction_identity<double>::min();
    }

    KOKKOS_INLINE_FUNCTION
    void join( value_type& destination, const value_type& source ) const
    {
        for ( int n = 0; n < 12; ++n )
            destination.value[n] =
                Kokkos::fmin( destination.value[n], source.value[n] );
    }

    KOKKOS_INLINE_FUNCTION
    void transform( const double point[3], double transformed[3] ) const
    {
        if ( scan_path_frame )
        {
            transformed[0] = direction_x * point[0] + direction_y * point[1];
            transformed[1] = -direction_y * point[0] + direction_x * point[1];
        }
        else
        {
            transformed[0] = point[0];
            transformed[1] = point[1];
        }
        transformed[2] = point[2];
    }

    KOKKOS_INLINE_FUNCTION
    void addPoint( const int isotherm, const double point[3],
                   value_type& bounds ) const
    {
        double transformed[3];
        transform( point, transformed );
        const int offset = 6 * isotherm;
        for ( int d = 0; d < 3; ++d )
        {
            bounds.value[offset + d] =
                Kokkos::fmin( bounds.value[offset + d], transformed[d] );
            bounds.value[offset + 3 + d] =
                Kokkos::fmin( bounds.value[offset + 3 + d], -transformed[d] );
        }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const int i, const int j, const int k,
                     value_type& bounds ) const
    {
        const int index[3] = { i, j, k };
        const double value = temperature( i, j, k, 0 );
        const double thresholds[2] = { solidus, liquidus };
        bool candidate = value >= minimum_threshold;
        bool has_edge[3] = { false, false, false };
        double neighbor_value[3] = { 0.0, 0.0, 0.0 };

        // Each positive grid edge is visited once. At an MPI subdomain high
        // side, the endpoint comes from the already-current temperature halo.
        for ( int d = 0; d < 3; ++d )
        {
            if ( index[d] + 1 >= owned_high[d] && on_high_boundary[d] )
                continue;

            has_edge[d] = true;
            int neighbor[3] = { i, j, k };
            ++neighbor[d];
            neighbor_value[d] =
                temperature( neighbor[0], neighbor[1], neighbor[2], 0 );
            candidate = candidate || crossesIsotherm( value, neighbor_value[d],
                                                      minimum_threshold );
        }
        if ( !candidate )
            return;

        double point[3];
        for ( int d = 0; d < 3; ++d )
            point[d] = origin[d] + static_cast<double>( index[d] ) * spacing;

        for ( int iso = 0; iso < 2; ++iso )
            if ( enabled[iso] && value >= thresholds[iso] )
                addPoint( iso, point, bounds );

        for ( int d = 0; d < 3; ++d )
        {
            if ( !has_edge[d] )
                continue;
            for ( int iso = 0; iso < 2; ++iso )
            {
                if ( !enabled[iso] )
                    continue;
                if ( !crossesIsotherm( value, neighbor_value[d],
                                       thresholds[iso] ) )
                    continue;
                double crossing[3] = { point[0], point[1], point[2] };
                crossing[d] +=
                    spacing * isothermFraction( value, neighbor_value[d],
                                                thresholds[iso] );
                addPoint( iso, crossing, bounds );
            }
        }
    }
};

template <typename MemorySpace>
class MeltPoolDimensions
{
  public:
    MeltPoolDimensions( const MeltPoolDimensionsInput& input,
                        const double solidus, const double liquidus,
                        MPI_Comm comm )
        : input_( input )
        , solidus_( solidus )
        , liquidus_( liquidus )
        , comm_( comm )
    {
        MPI_Comm_rank( comm_, &rank_ );

        int output_error = 0;
        if ( rank_ == 0 )
        {
            std::error_code error;
            std::filesystem::create_directories( input_.directory, error );
            if ( error )
            {
                std::cerr << "Cannot create melt-pool dimensions directory "
                          << input_.directory << ": " << error.message()
                          << std::endl;
                output_error = 1;
            }
            else
            {
                output_.open( input_.directory + "/dimensions.csv" );
                if ( !output_ )
                    output_error = 1;
            }
        }
        MPI_Bcast( &output_error, 1, MPI_INT, 0, comm_ );
        if ( output_error )
            throw std::runtime_error(
                "Unable to initialize melt-pool dimensions output" );

        if ( rank_ == 0 )
        {
            output_ << std::setprecision( 17 );
            output_ << "time";
            const char* names[2] = { "solidus", "liquidus" };
            for ( int iso = 0; iso < 2; ++iso )
            {
                if ( !input_.isotherms[iso] )
                    continue;
                output_ << ',' << names[iso] << "_active";
                if ( input_.coordinate_frame == "scan_path" )
                    output_ << ',' << names[iso] << "_length," << names[iso]
                            << "_width," << names[iso] << "_depth";
                else
                    output_ << ',' << names[iso] << "_x_extent," << names[iso]
                            << "_y_extent," << names[iso] << "_z_extent";
            }
            output_ << '\n';
        }
    }

    template <typename ExecutionSpace>
    void update( const ExecutionSpace& execution_space, Grid<MemorySpace>& grid,
                 const std::array<double, 3>& beam_direction,
                 const double time )
    {
        Kokkos::Profiling::ScopedRegion region( "Finch::melt_pool_dimensions" );
        auto local_mesh = grid.getLocalMesh();
        auto owned = grid.getIndexSpace();

        MeltPoolBoundsFunctor<decltype( grid.getTemperature() )> functor;
        functor.temperature = grid.getTemperature();
        functor.solidus = solidus_;
        functor.liquidus = liquidus_;
        functor.minimum_threshold = input_.isotherms[0] ? solidus_ : liquidus_;
        functor.spacing =
            grid.getLocalGrid()->globalGrid().globalMesh().cellSize( 0 );
        for ( int d = 0; d < 3; ++d )
            functor.origin[d] =
                local_mesh.lowCorner( Cabana::Grid::Ghost(), d );
        functor.scan_path_frame =
            input_.coordinate_frame == "scan_path" ? 1 : 0;
        for ( int iso = 0; iso < 2; ++iso )
            functor.enabled[iso] = input_.isotherms[iso] ? 1 : 0;

        functor.direction_x = beam_direction[0];
        functor.direction_y = beam_direction[1];
        for ( int d = 0; d < 3; ++d )
        {
            functor.owned_high[d] = static_cast<int>( owned.max( d ) );
            functor.on_high_boundary[d] =
                local_mesh.onHighBoundary( d ) ? 1 : 0;
        }

        MeltPoolBoundsValue bounds;
        Cabana::Grid::grid_parallel_reduce( "Finch::melt_pool_bounds",
                                            execution_space, owned, functor,
                                            bounds );
        execution_space.fence( "Finch melt-pool bounds" );

        if ( grid.comm_size > 1 )
            MPI_Allreduce( MPI_IN_PLACE, bounds.value, 12, MPI_DOUBLE, MPI_MIN,
                           comm_ );

        if ( rank_ == 0 )
        {
            output_ << time;
            for ( int iso = 0; iso < 2; ++iso )
            {
                if ( !input_.isotherms[iso] )
                    continue;
                const int offset = 6 * iso;
                const bool active = bounds.value[offset] !=
                                    Kokkos::reduction_identity<double>::min();
                output_ << ',' << ( active ? 1 : 0 );
                for ( int d = 0; d < 3; ++d )
                {
                    const double extent = active
                                              ? -bounds.value[offset + 3 + d] -
                                                    bounds.value[offset + d]
                                              : 0.0;
                    output_ << ',' << extent;
                }
            }
            output_ << '\n';
        }
    }

  private:
    MeltPoolDimensionsInput input_;
    double solidus_ = 0.0;
    double liquidus_ = 0.0;
    MPI_Comm comm_ = MPI_COMM_NULL;
    int rank_ = 0;
    std::ofstream output_;
};

} // namespace Finch

#endif
