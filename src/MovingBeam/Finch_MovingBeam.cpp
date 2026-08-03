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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Finch_MovingBeam.hpp"

namespace Finch
{

MovingBeam::MovingBeam( const std::string& scan_path_file, MPI_Comm comm )
    : path( 1, Segment() )
    , index_( 0 )
    , position_( { 0.0, 0.0, 0.0 } )
    , power_( 0.0 )
    , endTime_( 0.0 )
    , current_time_( 0.0 )
    , comm_( comm )
{
    // read the scan path file
    pFile_ = scan_path_file;
    readPath();

    // find the beam end time, i.e. the last time power is on
    for ( std::size_t i = path.size() - 1; i > 0; i-- )
    {
        if ( path[i].power() > eps )
        {
            endTime_ = path[i].time();
            break;
        }
    }
}

void MovingBeam::readPath()
{
    std::string contents;
    std::string error;
    int rank = 0;
    if ( comm_ != MPI_COMM_NULL )
        MPI_Comm_rank( comm_, &rank );

    if ( rank == 0 )
    {
        std::ifstream input( pFile_ );
        if ( !input )
            error = "Cannot open scan-path file " + pFile_;
        else
        {
            std::ostringstream buffer;
            buffer << input.rdbuf();
            contents = buffer.str();
        }
    }

    if ( comm_ != MPI_COMM_NULL )
    {
        int error_size = static_cast<int>( error.size() );
        MPI_Bcast( &error_size, 1, MPI_INT, 0, comm_ );
        if ( error_size > 0 )
        {
            error.resize( error_size );
            MPI_Bcast( error.data(), error_size, MPI_CHAR, 0, comm_ );
            throw std::runtime_error( error );
        }

        int contents_size = static_cast<int>( contents.size() );
        MPI_Bcast( &contents_size, 1, MPI_INT, 0, comm_ );
        contents.resize( contents_size );
        if ( contents_size > 0 )
            MPI_Bcast( contents.data(), contents_size, MPI_CHAR, 0, comm_ );
    }
    else if ( !error.empty() )
        throw std::runtime_error( error );

    std::istringstream is( contents );

    std::string line;

    // skip the header line
    std::getline( is, line );

    while ( std::getline( is, line ) )
    {
        if ( line.empty() )
        {
            continue;
        }

        path.push_back( Segment( line ) );
    }

    if ( path.size() == 1 )
        throw std::runtime_error( "Scan-path file contains no path rows" );
    if ( path[1].mode() != 1 )
        throw std::runtime_error(
            "The first scan-path row must define a point source" );

    for ( std::size_t i = 1; i < path.size(); i++ )
    {
        if ( path[i].mode() == 1 )
        {
            path[i].setTime( path[i - 1].time() + path[i].parameter() );
        }
        else
        {
            const auto& p0 = path[i - 1].position();
            const auto& p1 = path[i].position();

            double d_ = sqrt( ( p0[0] - p1[0] ) * ( p0[0] - p1[0] ) +
                              ( p0[1] - p1[1] ) * ( p0[1] - p1[1] ) +
                              ( p0[2] - p1[2] ) * ( p0[2] - p1[2] ) );

            path[i].setTime( path[i - 1].time() + d_ / path[i].parameter() );
        }
    }
}

void MovingBeam::move( const double time )
{
    current_time_ = time;
    // turn off the laser power and stop position update at the end of the path
    if ( ( time - endTime_ ) > eps )
    {
        power_ = 0.0;
        return;
    }

    // update the current index of the path
    index_ = findIndex( time );

    const int i = index_;

    // update the beam center
    if ( path[i].mode() == 1 )
    {
        position_ = path[i].position();
    }
    else
    {
        double dt = path[i].time() - path[i - 1].time();

        if ( dt > 0 )
        {
            for ( int d = 0; d < 3; ++d )
            {
                const double dx =
                    path[i].position()[d] - path[i - 1].position()[d];
                position_[d] = path[i - 1].position()[d] +
                               dx * ( time - path[i - 1].time() ) / dt;
            }
        }
        else
            position_ = path[i].position();
    }

    // update the beam power
    if ( ( time - path[i - 1].time() ) > eps )
    {
        power_ = path[i].power();
    }
    else
    {
        power_ = path[i - 1].power();
    }
}

int MovingBeam::findIndex( const double time )
{
    const int n = path.size() - 1;

    // step back path index for safe updating
    int i = index_;
    for ( i = index_; i > 0 && path[i].time() > time; --i )
    {
    }

    // update the path index to the provided time
    for ( ; i < n && path[i].time() < time; ++i )
    {
    }

    // skip any point sources with zero time
    while ( i < n )
    {
        if ( path[i].mode() == 1 && path[i].parameter() == 0 )
        {
            ++i;
        }
        else
        {
            break;
        }
    }

    return std::min( std::max( i, 1 ), n );
}

bool MovingBeam::activePath() const { return current_time_ <= endTime_ + eps; }

} // namespace Finch
