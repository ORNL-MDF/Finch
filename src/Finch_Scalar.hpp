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
  \file Finch_Scalar.hpp
  \brief Scalar-type dispatch so the solver can run on types other than double
*/

#ifndef Finch_Scalar_H
#define Finch_Scalar_H

#include <Kokkos_Core.hpp>

#include <cmath>
#include <type_traits>

namespace Finch
{
namespace Math
{

// Finch's kernels call these instead of Kokkos::<fn> directly so that the field
// scalar type is not required to be a built-in float/double.
//
// For arithmetic types the call forwards to Kokkos::<fn>, which is what the
// solver used before this indirection existed -- same device path, same
// architecture-specific implementation, no change in generated code.
//
// For any other type the call is made unqualified so that argument-dependent
// lookup finds an overload in the scalar type's own namespace. That is the
// extension point: a user-defined scalar type supplies its own exp/fmin/fmax
// and needs no edit here. Finch therefore carries no dependency on any
// particular AD library.

template <class T>
KOKKOS_INLINE_FUNCTION auto exp( const T& x )
{
    if constexpr ( std::is_arithmetic<T>::value )
    {
        return Kokkos::exp( x );
    }
    else
    {
        using std::exp;
        return exp( x );
    }
}

template <class T>
KOKKOS_INLINE_FUNCTION auto fmin( const T& a, const T& b )
{
    if constexpr ( std::is_arithmetic<T>::value )
    {
        return Kokkos::fmin( a, b );
    }
    else
    {
        using std::fmin;
        return fmin( a, b );
    }
}

template <class T>
KOKKOS_INLINE_FUNCTION auto fmax( const T& a, const T& b )
{
    if constexpr ( std::is_arithmetic<T>::value )
    {
        return Kokkos::fmax( a, b );
    }
    else
    {
        using std::fmax;
        return fmax( a, b );
    }
}

// Numeric value of a scalar. The trait exposes the underlying arithmetic type
// so an AD scalar backed by float remains float and one backed by double remains
// double. Call sites whose storage format requires double (file output and the
// solidification event records handed to downstream tools) convert at that
// boundary rather than forcing every scalar integration to double here.
//
// The primary template covers built-in types. A non-arithmetic scalar type
// specializes this to say which of its components is the "value" -- for a
// forward-mode AD type that is the real/zeroth coefficient.
template <class T, class Enable = void>
struct ScalarValue
{
    using value_type = T;

    KOKKOS_INLINE_FUNCTION static value_type value( const T& x ) { return x; }
};

template <class T>
using scalar_value_t = typename ScalarValue<T>::value_type;

template <class T>
KOKKOS_INLINE_FUNCTION scalar_value_t<T> value( const T& x )
{
    return ScalarValue<T>::value( x );
}

} // namespace Math
} // namespace Finch

#endif
