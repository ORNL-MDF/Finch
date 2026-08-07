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

/****************************************************************************
 * Sparrow OTI (order-truncated imaginary) scalar type glue for Finch.
 *
 * This header is the ONLY place where Finch and Sparrow meet. It lives in
 * integrations/sparrow rather than Finch's core or an application: the core
 * carries no dependency on any AD library, and the coupling is a single
 * ScalarValue specialization plus a parameter-seeding helper.
 ****************************************************************************/

#ifndef Finch_Sparrow_H
#define Finch_Sparrow_H

#include <Finch_Scalar.hpp>
#include <Finch_Solver.hpp>

#include "otinum/otinum.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace Finch
{
namespace Math
{

// Tell Finch which component of an OTI number is its numeric value. This is all
// the core library needs to know about the type; every other operation it
// performs (arithmetic, comparison, exp, fmin/fmax) is resolved by ADL into
// namespace oti.
template <int M, int N, class Coeff>
struct ScalarValue<oti::otinum<M, N, Coeff>>
{
    using value_type = Coeff;

    KOKKOS_INLINE_FUNCTION static value_type
    value( const oti::otinum<M, N, Coeff>& x )
    {
        return x.real();
    }
};

} // namespace Math

namespace Sensitivity
{

// The differentiated parameter set for this study.
//
// two_sigma is treated as a single parameter applied to all three source
// directions, which is what the input deck describes for an axisymmetric spot.
// Splitting it into three independent directions is a matter of adding two
// more slots below and seeding them separately.
enum Parameter
{
    Density = 0,
    SpecificHeat,
    ThermalConductivity,
    LatentHeat,
    Absorption,
    TwoSigma,
    NumParameters
};

inline const char* name( int p )
{
    switch ( p )
    {
    case Density:
        return "density";
    case SpecificHeat:
        return "specific_heat";
    case ThermalConductivity:
        return "thermal_conductivity";
    case LatentHeat:
        return "latent_heat";
    case Absorption:
        return "absorption";
    case TwoSigma:
        return "two_sigma";
    default:
        throw std::out_of_range(
            "Sensitivity parameter index is out of range" );
    }
}

inline const char* units( int p )
{
    switch ( p )
    {
    case Density:
        return "kg/m^3";
    case SpecificHeat:
        return "J/kg/K";
    case ThermalConductivity:
        return "W/m/K";
    case LatentHeat:
        return "J/kg";
    case Absorption:
        return "-";
    case TwoSigma:
        return "m";
    default:
        throw std::out_of_range(
            "Sensitivity parameter index is out of range" );
    }
}

// Nominal value of each parameter, read from the input deck.
inline std::array<double, NumParameters> nominal( const Inputs& db )
{
    // The single two_sigma slot presumes an axisymmetric source. Fail loudly
    // rather than silently differentiate only one direction.
    if ( db.source.two_sigma[0] != db.source.two_sigma[1] ||
         db.source.two_sigma[0] != db.source.two_sigma[2] )
        throw std::runtime_error(
            "Sensitivity: two_sigma is seeded as a single parameter, which "
            "requires the three components to be equal in the input deck." );

    std::array<double, NumParameters> p;
    p[Density] = db.properties.density;
    p[SpecificHeat] = db.properties.specific_heat;
    p[ThermalConductivity] = db.properties.thermal_conductivity;
    p[LatentHeat] = db.properties.latent_heat;
    p[Absorption] = db.source.absorption;
    p[TwoSigma] = db.source.two_sigma[0];
    return p;
}

// Assemble solver properties from an explicit parameter vector. Used for both
// the plain-double finite-difference runs (values perturbed) and the OTI run
// (values seeded as independent variables), so that both paths perturb exactly
// the same quantities.
template <class Scalar>
SolverParameters<Scalar> build( const std::array<Scalar, NumParameters>& p )
{
    SolverParameters<Scalar> params;
    params.density = p[Density];
    params.specific_heat = p[SpecificHeat];
    params.thermal_conductivity = p[ThermalConductivity];
    params.latent_heat = p[LatentHeat];
    params.absorption = p[Absorption];
    for ( int d = 0; d < 3; ++d )
        params.two_sigma[d] = p[TwoSigma];
    return params;
}

// Seed every parameter as an independent OTI variable at its nominal value.
template <class OTI>
std::array<OTI, NumParameters> seed( const Inputs& db )
{
    static_assert( OTI::nvars == NumParameters,
                   "OTI algebra must have one variable per differentiated "
                   "parameter" );
    auto nom = nominal( db );
    std::array<OTI, NumParameters> p;
    for ( int i = 0; i < NumParameters; ++i )
        p[i] = OTI::variable( i, nom[i] );
    return p;
}

} // namespace Sensitivity
} // namespace Finch

#endif
