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
  \file Version.hpp
  \brief Git version and hash
*/
#ifndef FINCH_VERSION_HPP
#define FINCH_VERSION_HPP

#include <Finch_Core_Config.hpp>

#include <string>

namespace Finch
{

//! Git version.
inline std::string version() { return Finch_VERSION_STRING; }

//! Git hash.
inline std::string commitHash() { return Finch_GIT_COMMIT_HASH; }

} // namespace Finch

#endif
