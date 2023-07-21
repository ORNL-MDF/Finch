#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "moving_beam.hpp"

const double MovingBeam::eps = 1e-10;

MovingBeam::MovingBeam()
    : path( 1, BeamSegment() )
    , index_( 0 )
    , power_( 0.0 )
{
    position_.resize( 3, 0.0 );

    // read the scan path file
    readPath();

    std::cout << "initial path index: " << index_ << std::endl;
}

void MovingBeam::readPath()
{

    const std::string pFile_ = "scanPath.txt";

    std::ifstream is( pFile_ );

    if ( !is.good() )
    {
        std::string error = "Cannot find file " + pFile_;
        throw std::runtime_error( error );
    }
    else
    {
        std::cout << "Reading scan path from: " << pFile_ << std::endl;
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

        path.push_back( BeamSegment( line ) );
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
    int i1;
    for ( i1 = index_; i1 > 0 && path[i1].time() > time; --i1 )
    {
    }

    // update the path index to the provided time
    int i;
    for ( i = i1; i < n && path[i].time() < time; ++i )
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
