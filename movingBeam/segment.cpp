#include "segment.H"
#include <fstream>
#include <iostream>
#include <sstream>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Set the segment time to provided value
void segment::setTime( double time ) { time_ = time; }

void segment::setPosition( std::vector<double> position )
{
    position_ = position;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// construct default segment as a zero point source
segment::segment()
    : mode_( 1 )
    , power_( 0.0 )
    , parameter_( 0.0 )
    , time_( 0.0 )
{
    position_.resize( 3, 0.0 );
}

// set the segement properties given a space-delimited string
segment::segment( std::string line )
{
    position_.resize( 3, 0.0 );
    std::stringstream lineStream( line );

    lineStream >> mode_ >> position_[0] >> position_[1] >> position_[2] >>
        power_ >> parameter_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
