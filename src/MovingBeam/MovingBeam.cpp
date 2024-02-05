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

#include "MovingBeam.hpp"

// a small floating point value used for consistent path update logic
const double MovingBeam::eps = 1e-10;

MovingBeam::MovingBeam( const std::string scan_path_file )
    : path( 1, Segment() )
    , index_( 0 )
    , power_( 0.0 )
    , endTime_( 0.0 )
{
    position_.resize( 3, 0.0 );

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
    std::ifstream is( pFile_ );

    if ( !is.good() )
    {
        std::string error = "Cannot find file " + pFile_;
        throw std::runtime_error( error );
    }

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

    for ( std::size_t i = 1; i < path.size(); i++ )
    {
        if ( path[i].mode() == 1 )
        {
            path[i].setTime( path[i - 1].time() + path[i].parameter() );
        }
        else
        {
            std::vector<double> p0 = path[i - 1].position();
            std::vector<double> p1 = path[i].position();

            double d_ = sqrt( ( p0[0] - p1[0] ) * ( p0[0] - p1[0] ) +
                              ( p0[1] - p1[1] ) * ( p0[1] - p1[1] ) +
                              ( p0[2] - p1[2] ) * ( p0[2] - p1[2] ) );

            path[i].setTime( path[i - 1].time() + d_ / path[i].parameter() );
        }
    }
}

void MovingBeam::move( const double time )
{
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
        std::vector<double> displacement( 3, 0 );

        double dt = path[i].time() - path[i - 1].time();

        if ( dt > 0 )
        {
            std::vector<double> dx( 3, 0 );
            dx[0] = path[i].position()[0] - path[i - 1].position()[0];
            dx[1] = path[i].position()[1] - path[i - 1].position()[1];
            dx[2] = path[i].position()[2] - path[i - 1].position()[2];
            displacement[0] = dx[0] * ( time - path[i - 1].time() ) / dt;
            displacement[1] = dx[1] * ( time - path[i - 1].time() ) / dt;
            displacement[2] = dx[2] * ( time - path[i - 1].time() ) / dt;
        }

        position_[0] = path[i - 1].position()[0] + displacement[0];
        position_[1] = path[i - 1].position()[1] + displacement[1];
        position_[2] = path[i - 1].position()[2] + displacement[2];
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

    return std::min( std::max( i, 0 ), n );
}
