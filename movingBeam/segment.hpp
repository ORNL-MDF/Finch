/*!
  \file moving_beam.hpp
  \brief Class for properties of a moving heat source
*/

#ifndef segment_H
#define segment_H

#include <string>

class segment
{
  private:
    //! 0 or 1 (1 = point source, 0 = line source)
    double mode_;
    //! position of the heat source center
    std::vector<double> position_;
    //! power of the heat source
    double power_;
    //! (mode = 1: time interval, mode = 0: scan velocity)
    double parameter_;

    double time_;

  public:
    //! Default construction
    segment();

    //! Construct from space-delimited string
    segment( std::string );

    //! Destructor
    virtual ~segment() {}

    // Member Functions

    //! Set time to provided value
    void setTime( double time );

    //! Set position to provided value
    void setPosition( std::vector<double> position );

    double mode() { return mode_; }

    std::vector<double> position() { return position_; }

    double power() { return power_; }

    double parameter() { return parameter_; }

    double time() { return time_; }
};

#endif

// ************************************************************************* //
