/*!
  \file moving_beam.hpp
  \brief Class for properties of a moving heat source
*/

#ifndef Segment_H
#define Segment_H

#include <string>
#include <vector>

class Segment
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
    //! current time
    double time_;

  public:
    //! Default construction
    Segment();

    //! Construct from space-delimited string
    Segment( std::string );

    //! Destructor
    virtual ~Segment() {}

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
