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

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Finch_Segment.hpp"

namespace Finch
{

void Segment::setTime( double time ) { time_ = time; }

void Segment::setPosition( const std::array<double, 3>& position )
{
    position_ = position;
}

Segment::Segment()
    : mode_( 1 )
    , power_( 0.0 )
    , parameter_( 0.0 )
    , time_( 0.0 )
{
    position_.fill( 0.0 );
}

Segment::Segment( const std::string& line )
{
    position_.fill( 0.0 );
    std::stringstream lineStream( line );

    if ( !( lineStream >> mode_ >> position_[0] >> position_[1] >>
            position_[2] >> power_ >> parameter_ ) )
        throw std::runtime_error( "Invalid scan-path row: " + line );

    std::string trailing;
    if ( lineStream >> trailing )
        throw std::runtime_error( "Unexpected scan-path data: " + line );
    if ( mode_ != 0 && mode_ != 1 )
        throw std::runtime_error( "Scan-path mode must be 0 or 1" );
    if ( !std::isfinite( position_[0] ) || !std::isfinite( position_[1] ) ||
         !std::isfinite( position_[2] ) || !std::isfinite( power_ ) ||
         !std::isfinite( parameter_ ) || power_ < 0.0 || parameter_ < 0.0 ||
         ( mode_ == 0 && parameter_ <= 0.0 ) )
        throw std::runtime_error(
            "Scan-path coordinates, power, and parameter must be finite; "
            "power and dwell time must be nonnegative and line speed must be "
            "positive" );
}

} // namespace Finch
