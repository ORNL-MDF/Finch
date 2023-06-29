/*!
  \file moving_beam.hpp
  \brief Class for moving beam heat source used in additive manufacturing
*/

#ifndef movingBeam_H
#define movingBeam_H

#include "segment.H"

class movingBeam
{
  private:
    //! List of segments
    std::vector<segment> path;

    //! index of path
    int index_;

    //! end time of path
    double endTime_;

    //! current position of moving beam
    std::vector<double> position_;

    //! current power of moving beam
    double power_;

    //! tolerance for scan path intervals
    static const double eps;

  public:
    //! Default constructor
    movingBeam();

    //! Destructor
    virtual ~movingBeam() {}

    //! Read the path file
    void readPath();

    //! Move the beam to the provided time
    void move( const double time );

    //! Returns the path index at the provided time
    int findIndex( const double time );

    //! Returns true if the simulation time is less than path endTime
    bool activePath();

    //! Return the current path index
    int index() const { return index_; }

    //! Return end time of path
    int endTime() const { return endTime_; }

    //! Return current position of the moving beam
    std::vector<double> position() const { return position_; }

    //! Return current power of the moving beam
    double power() const { return power_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "movingBeam/movingBeam.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
