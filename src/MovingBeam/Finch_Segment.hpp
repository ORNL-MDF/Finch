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
  \file Segment.hpp
  \brief Class for properties of a single segment of a moving heat source
*/

#ifndef Segment_H
#define Segment_H

#include <array>
#include <string>

namespace Finch
{

class Segment
{
  private:
    //! 0 or 1 (1 = point source, 0 = line source)
    int mode_;
    //! position of the heat source center
    std::array<double, 3> position_;
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
    explicit Segment( const std::string& line );

    //! Set time to provided value
    void setTime( double time );

    //! Set position to provided value
    void setPosition( const std::array<double, 3>& position );

    int mode() const { return mode_; }

    const std::array<double, 3>& position() const { return position_; }

    double power() const { return power_; }

    double parameter() const { return parameter_; }

    double time() const { return time_; }
};

} // namespace Finch

#endif
