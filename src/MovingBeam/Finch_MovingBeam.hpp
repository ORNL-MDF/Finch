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

/*!
  \file MovingBeam.hpp
  \brief Class for moving beam heat source used in additive manufacturing
*/

#ifndef MovingBeam_H
#define MovingBeam_H

#include "MovingBeam/Finch_Segment.hpp"

#include <vector>

namespace Finch
{

class MovingBeam
{
  private:
    //! List of segments
    std::vector<Segment> path;

    //! Scan path file name
    std::string pFile_;

    //! index of path
    int index_;

    //! current position of moving beam
    std::vector<double> position_;

    //! current power of moving beam
    double power_;

    //! end time of path
    double endTime_;

    //! tolerance for scan path intervals
    static constexpr double eps = 1e-10;

  public:
    //! Default constructor
    MovingBeam( const std::string scan_path_file );

    //! Destructor
    virtual ~MovingBeam() {}

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

    //! Return current position component of the moving beam
    double position( const int dir ) const { return position_[dir]; }

    //! Return current power of the moving beam
    double power() const { return power_; }
};

} // namespace Finch

#endif
