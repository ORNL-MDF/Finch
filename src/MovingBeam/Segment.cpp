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

#include <fstream>
#include <iostream>
#include <sstream>

#include "Segment.hpp"

void Segment::setTime( double time ) { time_ = time; }

void Segment::setPosition( std::vector<double> position )
{
    position_ = position;
}

Segment::Segment()
    : mode_( 1 )
    , power_( 0.0 )
    , parameter_( 0.0 )
    , time_( 0.0 )
{
    position_.resize( 3, 0.0 );
}

Segment::Segment( std::string line )
{
    position_.resize( 3, 0.0 );
    std::stringstream lineStream( line );

    lineStream >> mode_ >> position_[0] >> position_[1] >> position_[2] >>
        power_ >> parameter_;
}
