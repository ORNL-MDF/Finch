/****************************************************************************
 * Copyright (c) 2026 by Oak Ridge National Laboratory                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of Finch. Finch is distributed under a                 *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef FINCH_BOUNDARY_CONDITIONS_HPP
#define FINCH_BOUNDARY_CONDITIONS_HPP

#include <array>
#include <string>

namespace Finch
{

inline constexpr std::array<const char*, 6> boundary_face_names = {
    "x_min", "x_max", "y_min", "y_max", "z_min", "z_max" };

struct BoundaryFaceCondition
{
    std::string type = "adiabatic";
    double value = 0.0;
    double convection_coefficient = 0.0;
    double emissivity = 0.0;
    double ambient_temperature = 0.0;
};

struct BoundaryConditions
{
    std::array<BoundaryFaceCondition, 6> faces;
};

} // namespace Finch

#endif
