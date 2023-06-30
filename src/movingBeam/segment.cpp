#include <fstream>
#include <iostream>
#include <sstream>

#include "segment.hpp"

void segment::setTime( double time ) { time_ = time; }

void segment::setPosition( std::vector<double> position )
{
    position_ = position;
}

segment::segment()
    : mode_( 1 )
    , power_( 0.0 )
    , parameter_( 0.0 )
    , time_( 0.0 )
{
    position_.resize( 3, 0.0 );
}

segment::segment( std::string line )
{
    position_.resize( 3, 0.0 );
    std::stringstream lineStream( line );

    lineStream >> mode_ >> position_[0] >> position_[1] >> position_[2] >>
        power_ >> parameter_;
}

// ************************************************************************* //
