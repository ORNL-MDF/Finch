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
  \file Solver.hpp
  \brief Main class for heat transport solve
*/

#ifndef Solver_H
#define Solver_H

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_Scalar.hpp>

namespace Finch
{

struct HostTag
{
};
struct DeviceTag
{
};

// The material and source inputs the solver treats as differentiable. Held in
// its own struct, templated on the scalar type, so a caller can hand the solver
// values of a type other than double (see makeProperties below). Quantities the
// solver only ever compares against -- solidus and liquidus -- stay double:
// they select a branch rather than entering the arithmetic.
template <typename Scalar>
struct MaterialProperties
{
    Scalar density;
    Scalar specific_heat;
    Scalar thermal_conductivity;
    Scalar latent_heat;
    Scalar absorption;
    Scalar two_sigma[3];
};

// Build the properties for a plain double solve directly from the input deck.
// Scalar defaults to double, so existing callers get exactly the previous
// behavior; callers wanting another scalar type request it explicitly and then
// overwrite the members they care about.
template <typename Scalar = double>
MaterialProperties<Scalar> makeProperties( const Inputs& db )
{
    MaterialProperties<Scalar> props;
    props.density = db.properties.density;
    props.specific_heat = db.properties.specific_heat;
    props.thermal_conductivity = db.properties.thermal_conductivity;
    props.latent_heat = db.properties.latent_heat;
    props.absorption = db.source.absorption;
    for ( std::size_t d = 0; d < 3; ++d )
        props.two_sigma[d] = db.source.two_sigma[d];
    return props;
}

template <typename ViewType, typename EntityType, typename LocalMeshType>
class Solver
{
  public:
    // The field scalar type is whatever the temperature view holds. No extra
    // template parameter is needed: making the Cabana array carry a different
    // value type is enough to change the arithmetic throughout the solver.
    using scalar_type = typename ViewType::non_const_value_type;

  protected:
    // temperature views are default constructed and updated every step.
    ViewType T_;
    ViewType T0_;

    LocalMeshType local_mesh_;

    // solution parameters
    double dt_;
    double solidus_;
    double liquidus_;
    scalar_type rho_cp_;
    scalar_type rho_Lf_by_dT_;
    scalar_type k_by_dx2_;

    // heat source parameters
    double power_;
    double position_[3];
    scalar_type r_[3];
    scalar_type A_inv_[3];
    scalar_type I0_;
    double w_max_;

  public:
    Solver( Inputs db, LocalMeshType local_mesh )
        : Solver( db, local_mesh, makeProperties<scalar_type>( db ) )
    {
    }

    Solver( Inputs db, LocalMeshType local_mesh,
            const MaterialProperties<scalar_type>& props )
        : local_mesh_( local_mesh )
        , power_( 0.0 )
    {
        // solution parameter constants
        double dx = db.space.cell_size;
        scalar_type rho = props.density;
        scalar_type cp = props.specific_heat;
        scalar_type Lf = props.latent_heat;

        dt_ = db.time.time_step;

        solidus_ = db.properties.solidus;

        liquidus_ = db.properties.liquidus;

        rho_cp_ = rho * cp;

        rho_Lf_by_dT_ = rho * Lf / ( liquidus_ - solidus_ );

        k_by_dx2_ = ( props.thermal_conductivity ) / ( dx * dx );

        // initialize beam position
        for ( std::size_t d = 0; d < 3; ++d )
        {
            position_[d] = 0.0;
        }

        // heat source parameter constants
        for ( std::size_t d = 0; d < 3; ++d )
        {
            r_[d] = props.two_sigma[d] / Kokkos::sqrt( 2.0 );
            A_inv_[d] = 1.0 / r_[d] / r_[d];
        }

        I0_ = ( 2.0 * props.absorption ) /
              ( M_PI * Kokkos::sqrt( M_PI ) * r_[0] * r_[1] * r_[2] );

        // cut off for 3 standard deviations from heat source center
        w_max_ = Kokkos::log( 3 ) + 2 * Kokkos::log( 10 );
    }

    // Function for temperature solve: forward time-centered space (FTCS) method
    template <class ExecSpace, class IndexSpaceType>
    void solve( ExecSpace exec_space, IndexSpaceType owned_space, ViewType& T,
                ViewType& T0, const double beam_power,
                const double beam_pos[3] )
    {
        // Update temperature views and beam parameters for current time step
        T_ = T;

        T0_ = T0;

        power_ = beam_power;

        for ( std::size_t d = 0; d < 3; ++d )
        {
            position_[d] = beam_pos[d];
        }

        // Tagged versions of temperature solver for architecture optimization
        using memory_space = typename ViewType::memory_space;

        if constexpr ( std::is_same<memory_space, Kokkos::HostSpace>::value )
        {
            Cabana::Grid::grid_parallel_for( "solve", exec_space, owned_space,
                                             HostTag{}, *this );
        }
        else
        {
            Cabana::Grid::grid_parallel_for( "solve", exec_space, owned_space,
                                             DeviceTag{}, *this );
        }
    }

    // Host tagged version of the temperature solver
    KOKKOS_INLINE_FUNCTION
    void operator()( HostTag tag, const int i, const int j, const int k ) const
    {
        scalar_type x = T0_( i, j, k, 0 );

        scalar_type dt_by_rho_cp = ( x >= solidus_ && x <= liquidus_ )
                                       ? dt_ / ( rho_cp_ + rho_Lf_by_dT_ )
                                       : dt_ / ( rho_cp_ );

        scalar_type rhs = laplacian( i, j, k ) + source( tag, i, j, k );

        T_( i, j, k, 0 ) = x + rhs * dt_by_rho_cp;
    }

    // Device tagged version of the temperature solver
    KOKKOS_INLINE_FUNCTION
    void operator()( DeviceTag tag, const int i, const int j,
                     const int k ) const
    {
        scalar_type x = T0_( i, j, k, 0 );

        scalar_type dt_by_rho_cp =
            dt_ / ( rho_cp_ +
                    ( x >= solidus_ ) * ( x <= liquidus_ ) * rho_Lf_by_dT_ );

        scalar_type rhs = laplacian( i, j, k ) + source( tag, i, j, k );

        T_( i, j, k, 0 ) = x + rhs * dt_by_rho_cp;
    }

    // First-order centered space laplacian stencil
    KOKKOS_INLINE_FUNCTION
    auto laplacian( const int i, const int j, const int k ) const
    {
        return ( T0_( i - 1, j, k, 0 ) + T0_( i + 1, j, k, 0 ) +
                 T0_( i, j - 1, k, 0 ) + T0_( i, j + 1, k, 0 ) +
                 T0_( i, j, k - 1, 0 ) + T0_( i, j, k + 1, 0 ) -
                 6.0 * T0_( i, j, k, 0 ) ) *
               k_by_dx2_;
    }

    // Normalized weight for the gaussian source term: x in exp(-x)
    KOKKOS_INLINE_FUNCTION
    scalar_type weight( const int i, const int j, const int k ) const
    {
        // Mesh coordinates and beam position are geometry, not field values,
        // and stay double whatever the field scalar type is. Only the
        // accumulation below picks up the scalar type, through A_inv_.
        double grid_loc[3];
        double dist_to_beam[3];
        int idx[3] = { i, j, k };

        local_mesh_.coordinates( EntityType(), idx, grid_loc );

        dist_to_beam[0] = grid_loc[0] - position_[0];
        dist_to_beam[1] = grid_loc[1] - position_[1];
        dist_to_beam[2] = grid_loc[2] - position_[2];

        return ( dist_to_beam[0] * dist_to_beam[0] * A_inv_[0] ) +
               ( dist_to_beam[1] * dist_to_beam[1] * A_inv_[1] ) +
               ( dist_to_beam[2] * dist_to_beam[2] * A_inv_[2] );
    }

    // Heating source term, device overload.
    KOKKOS_INLINE_FUNCTION
    scalar_type source( DeviceTag, const int i, const int j, const int k ) const
    {
        return I0_ * power_ * Math::exp( -weight( i, j, k ) );
    }

    // Heating source term, host overload.
    KOKKOS_INLINE_FUNCTION
    scalar_type source( HostTag, const int i, const int j, const int k ) const
    {
        // performance improvements on host: scoping the exponential
        if ( power_ )
        {
            scalar_type w = weight( i, j, k );

            if ( w < w_max_ )
            {
                return I0_ * power_ * Math::exp( -w );
            }
            else
            {
                return scalar_type( 0 );
            }
        }
        else
        {
            return scalar_type( 0 );
        }
    }
};

// Create a solver based on the grid details and simulation inputs.
template <typename MemorySpace, typename Scalar>
auto createSolver( Inputs db, Grid<MemorySpace, Scalar> grid )
{
    using grid_type = Grid<MemorySpace, Scalar>;
    using entity_type = typename grid_type::entity_type;
    using view_type = typename grid_type::view_type;
    using mesh_type = typename grid_type::local_mesh_type;

    auto local_mesh = grid.getLocalMesh();

    return Solver<view_type, entity_type, mesh_type>( db, local_mesh );
}

// Create a solver with explicitly supplied material properties. Used when the
// properties carry more than a value -- for example seeded AD variables.
template <typename MemorySpace, typename Scalar>
auto createSolver( Inputs db, Grid<MemorySpace, Scalar> grid,
                   const MaterialProperties<Scalar>& props )
{
    using grid_type = Grid<MemorySpace, Scalar>;
    using entity_type = typename grid_type::entity_type;
    using view_type = typename grid_type::view_type;
    using mesh_type = typename grid_type::local_mesh_type;

    auto local_mesh = grid.getLocalMesh();

    return Solver<view_type, entity_type, mesh_type>( db, local_mesh, props );
}

} // namespace Finch

#endif
