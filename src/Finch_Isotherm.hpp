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
  \file Finch_Isotherm.hpp
  \brief Structured-grid isotherm crossing and local depth queries.
*/

#ifndef FINCH_ISOTHERM_HPP
#define FINCH_ISOTHERM_HPP

#include <algorithm>
#include <array>
#include <cmath>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <Finch_Grid.hpp>

namespace Finch
{

KOKKOS_INLINE_FUNCTION
bool crossesIsotherm( const double first, const double second,
                      const double threshold )
{
    return ( first < threshold && second >= threshold ) ||
           ( second < threshold && first >= threshold );
}

KOKKOS_INLINE_FUNCTION
double isothermFraction( const double first, const double second,
                         const double threshold )
{
    const double fraction = ( threshold - first ) / ( second - first );
    return Kokkos::fmin( Kokkos::fmax( fraction, 0.0 ), 1.0 );
}

template <typename TemperatureView, typename LocalMesh>
struct LocalDepthFunctor
{
    TemperatureView temperature;
    LocalMesh local_mesh;
    double beam_x = 0.0;
    double beam_y = 0.0;
    double beam_z = 0.0;
    double radius_squared = 0.0;
    double threshold = 0.0;
    double spacing = 0.0;
    int owned_high_z = 0;
    int on_high_z_boundary = 0;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int i, const int j, const int k,
                     double& deepest_depth ) const
    {
        const int index[3] = { i, j, k };
        double point[3];
        local_mesh.coordinates( Cabana::Grid::Node(), index, point );
        const double dx = point[0] - beam_x;
        const double dy = point[1] - beam_y;
        if ( dx * dx + dy * dy > radius_squared )
            return;

        const double value = temperature( i, j, k, 0 );
        if ( value >= threshold )
            deepest_depth = Kokkos::fmax( deepest_depth, beam_z - point[2] );

        const bool has_upper_edge = k + 1 < owned_high_z || !on_high_z_boundary;
        if ( !has_upper_edge )
            return;

        const double upper_z = point[2] + spacing;
        if ( upper_z > beam_z )
            return;
        const double upper_value = temperature( i, j, k + 1, 0 );
        if ( crossesIsotherm( value, upper_value, threshold ) )
        {
            const double fraction =
                isothermFraction( value, upper_value, threshold );
            const double crossing_z = point[2] + fraction * spacing;
            deepest_depth = Kokkos::fmax( deepest_depth, beam_z - crossing_z );
        }
    }
};

template <typename ExecutionSpace, typename MemorySpace>
double findLocalMeltPoolDepth( const ExecutionSpace& execution_space,
                               Grid<MemorySpace>& grid,
                               const std::array<double, 3>& beam_position,
                               const double search_radius,
                               const double threshold )
{
    Kokkos::Profiling::ScopedRegion region( "Finch::transient_depth_query" );
    auto local_mesh = grid.getLocalMesh();
    auto owned = grid.getIndexSpace();
    const double spacing =
        grid.getLocalGrid()->globalGrid().globalMesh().cellSize( 0 );
    const double ghost_x = local_mesh.lowCorner( Cabana::Grid::Ghost(), 0 );
    const double ghost_y = local_mesh.lowCorner( Cabana::Grid::Ghost(), 1 );
    const double ghost_z = local_mesh.lowCorner( Cabana::Grid::Ghost(), 2 );

    std::array<long, 3> lower = {
        std::max<long>(
            owned.min( 0 ),
            static_cast<long>( std::ceil(
                ( beam_position[0] - search_radius - ghost_x ) / spacing ) ) ),
        std::max<long>(
            owned.min( 1 ),
            static_cast<long>( std::ceil(
                ( beam_position[1] - search_radius - ghost_y ) / spacing ) ) ),
        owned.min( 2 ) };
    std::array<long, 3> upper = {
        std::min<long>(
            owned.max( 0 ),
            static_cast<long>( std::floor(
                ( beam_position[0] + search_radius - ghost_x ) / spacing ) ) +
                1 ),
        std::min<long>(
            owned.max( 1 ),
            static_cast<long>( std::floor(
                ( beam_position[1] + search_radius - ghost_y ) / spacing ) ) +
                1 ),
        std::min<long>( owned.max( 2 ),
                        static_cast<long>( std::floor(
                            ( beam_position[2] - ghost_z ) / spacing ) ) +
                            1 ) };
    for ( int d = 0; d < 3; ++d )
        upper[d] = std::max( upper[d], lower[d] );

    double local_depth = 0.0;
    const Cabana::Grid::IndexSpace<3> search_space( lower, upper );
    if ( search_space.size() > 0 )
    {
        LocalDepthFunctor<decltype( grid.getTemperature() ),
                          decltype( local_mesh )>
            functor;
        functor.temperature = grid.getTemperature();
        functor.local_mesh = local_mesh;
        functor.beam_x = beam_position[0];
        functor.beam_y = beam_position[1];
        functor.beam_z = beam_position[2];
        functor.radius_squared = search_radius * search_radius;
        functor.threshold = threshold;
        functor.spacing = spacing;
        functor.owned_high_z = static_cast<int>( owned.max( 2 ) );
        functor.on_high_z_boundary = local_mesh.onHighBoundary( 2 ) ? 1 : 0;

        Kokkos::Max<double> reducer( local_depth );
        Cabana::Grid::grid_parallel_reduce( "Finch::local_transient_depth",
                                            execution_space, search_space,
                                            functor, reducer );
        execution_space.fence( "Finch local transient depth" );
    }

    double global_depth = local_depth;
    if ( grid.comm_size > 1 )
        MPI_Allreduce( &local_depth, &global_depth, 1, MPI_DOUBLE, MPI_MAX,
                       grid.getComm() );
    return global_depth;
}

} // namespace Finch

#endif
