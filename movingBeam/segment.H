/*---------------------------------------------------------------------------*\
Class
    segment

Description
    Class for properties of a moving heat source

    mode:      0 or 1 (1 = point source, 0 = line source)
    position:  position of the heat source center
    power:     power of the heat source
    parameter: (mode = 1: time interval, mode = 0: scan velocity)

SourceFiles
    segment.C

\*---------------------------------------------------------------------------*/

#ifndef segment_H
#define segment_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include <string>
/*---------------------------------------------------------------------------*\
                      Class segment Declaration
\*---------------------------------------------------------------------------*/

class segment
{
  private:
    // Private data

    double mode_;

    std::vector<double> position_;

    double power_;

    double parameter_;

    double time_;

  public:
    // Constructors

    //- Default construction
    segment();

    //- Construct from space-delimited string
    segment( std::string );

    //- Destructor
    virtual ~segment() {}

    // Member Functions

    //- Set time to provided value
    void setTime( double time );

    //- Set position to provided value
    void setPosition( std::vector<double> position );

    double mode() { return mode_; }

    std::vector<double> position() { return position_; }

    double power() { return power_; }

    double parameter() { return parameter_; }

    double time() { return time_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "segment.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
