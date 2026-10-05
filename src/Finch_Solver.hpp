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

namespace Finch
{

struct HostTag
{
};
struct DeviceTag
{
};

template <typename ViewType, typename EntityType, typename LocalMeshType>
class Solver
{
  protected:
    // temperature views are default constructed and updated every step.
    ViewType T_;
    ViewType T0_;

    LocalMeshType local_mesh_;

    // solution parameters
    double dt_;
    double solidus_;
    double liquidus_;
    double rho_cp_;
    double rho_Lf_by_dT_;
    double k_by_dx2_;

    // heat source parameters
    double power_;
    double position_[3];

    // Original Gaussian heat source parameters
    double r_[3];
    double A_inv_[3];

    // Shared
    double I0_;
    double w_max_;

    // Dynamic super-Gaussian parameters (Coleman et al. 2024)
    // Only used when source.shape == "dynamic"
    bool dynamic_beam_;
    double sigma_[3]; // half-widths for super-Gaussian radial profile
    double k_;        // radial distribution parameter
    double m_;        // volumetric shape parameter
    double depth_;    // heat source depth (= two_sigma[2])

  public:
    Solver( Inputs db, LocalMeshType local_mesh )
        : local_mesh_( local_mesh )
        , power_( 0.0 )
    {
        // solution parameter constants
        double dx = db.space.cell_size;
        double rho = db.properties.density;
        double cp = db.properties.specific_heat;
        double Lf = db.properties.latent_heat;

        dt_ = db.time.time_step;

        solidus_ = db.properties.solidus;

        liquidus_ = db.properties.liquidus;

        rho_cp_ = rho * cp;

        rho_Lf_by_dT_ = rho * Lf / ( liquidus_ - solidus_ );

        k_by_dx2_ = ( db.properties.thermal_conductivity ) / ( dx * dx );

        // initialize beam position
        for ( std::size_t d = 0; d < 3; ++d )
        {
            position_[d] = 0.0;
        }

        dynamic_beam_ = ( db.source.shape == "dynamic" );

        if ( !dynamic_beam_ )
        {
            // Original isotropic 3D Gaussian
            for ( std::size_t d = 0; d < 3; ++d )
            {
                r_[d] = db.source.two_sigma[d] / Kokkos::sqrt( 2.0 );
                A_inv_[d] = 1.0 / r_[d] / r_[d];
            }

            I0_ = ( 2.0 * db.source.absorption ) /
                  ( M_PI * Kokkos::sqrt( M_PI ) * r_[0] * r_[1] * r_[2] );
        }
        else
        {
            // Dynamic super-Gaussian (Coleman et al. 2024, Eq. 2-6)
            k_ = db.source.k;
            m_ = db.source.m;

            // sigma_i = two_sigma_i / (2 * 2^(1/k))
            double two_to_1_over_k = Kokkos::pow( 2.0, 1.0 / k_ );
            for ( std::size_t d = 0; d < 3; ++d )
            {
                sigma_[d] = db.source.two_sigma[d] / ( 2.0 * two_to_1_over_k );
            }

            // Heat source depth = two_sigma[2]
            depth_ = db.source.two_sigma[2];

            // Areal normalization A0 (Eq. 3)
            double two_to_2_over_k = Kokkos::pow( 2.0, 2.0 / k_ );
            double rx0 = db.source.two_sigma[0] / two_to_2_over_k;
            double ry0 = db.source.two_sigma[1] / two_to_2_over_k;
            double A0 = M_PI * tgamma( 1.0 + 2.0 / k_ ) * rx0 * ry0;

            // Depth integral factor from Eq. 6
            double depth_factor =
                ( tgamma( 1.0 + 1.0 / m_ ) * tgamma( 1.0 + 2.0 / m_ ) ) /
                tgamma( 1.0 + 3.0 / m_ );

            // Volume normalization V0 = A0 * depth * depth_factor (Eq. 6)
            double V0 = A0 * depth_ * depth_factor;

            I0_ = db.source.absorption / V0;
        }

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
        double x = T0_( i, j, k, 0 );

        double dt_by_rho_cp = ( x >= solidus_ && x <= liquidus_ )
                                  ? dt_ / ( rho_cp_ + rho_Lf_by_dT_ )
                                  : dt_ / ( rho_cp_ );

        double rhs = laplacian( i, j, k ) + source( tag, i, j, k );

        T_( i, j, k, 0 ) = x + rhs * dt_by_rho_cp;
    }

    // Device tagged version of the temperature solver
    KOKKOS_INLINE_FUNCTION
    void operator()( DeviceTag tag, const int i, const int j,
                     const int k ) const
    {
        double x = T0_( i, j, k, 0 );

        double dt_by_rho_cp =
            dt_ / ( rho_cp_ +
                    ( x >= solidus_ ) * ( x <= liquidus_ ) * rho_Lf_by_dT_ );

        double rhs = laplacian( i, j, k ) + source( tag, i, j, k );

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

    // Normalized weight for the source term: returned value w used as exp(-w)
    // Branches on dynamic_beam_ set at construction time.
    KOKKOS_INLINE_FUNCTION
    auto weight( const int i, const int j, const int k ) const
    {
        double grid_loc[3];
        int idx[3] = { i, j, k };
        local_mesh_.coordinates( EntityType(), idx, grid_loc );

        double dx = grid_loc[0] - position_[0];
        double dy = grid_loc[1] - position_[1];
        double dz = grid_loc[2] - position_[2];

        if ( !dynamic_beam_ )
        {
            //  Original isotropic 3D Gaussian
            return ( dx * dx * A_inv_[0] ) + ( dy * dy * A_inv_[1] ) +
                   ( dz * dz * A_inv_[2] );
        }
        else
        {
            //  Dynamic super-Gaussian (Coleman et al. 2024, Eq. 4-5)
            // Normalized depth: only apply heat downward from beam center
            double dz_norm = Kokkos::fabs( dz ) / depth_;
            if ( dz_norm >= 1.0 )
                return w_max_ + 1.0; // below heat source depth: no contribution

            // Radius decays with depth: r(z) = sigma * (1 - |dz/d|^m)^(1/m)
            double shape_factor =
                Kokkos::pow( 1.0 - Kokkos::pow( dz_norm, m_ ), 1.0 / m_ );
            double rx = sigma_[0] * shape_factor;
            double ry = sigma_[1] * shape_factor;

            // Guard against division by zero near the tip
            if ( rx < 1e-15 || ry < 1e-15 )
                return w_max_ + 1.0;

            // Super-Gaussian radial exponent: (dx^2/rx^2 + dy^2/ry^2)^(k/2)
            double radial =
                ( dx * dx ) / ( rx * rx ) + ( dy * dy ) / ( ry * ry );
            return Kokkos::pow( radial, k_ / 2.0 );
        }
    }

    // Heating source term, device overload.
    KOKKOS_INLINE_FUNCTION
    auto source( DeviceTag, const int i, const int j, const int k ) const
    {
        return I0_ * power_ * Kokkos::exp( -weight( i, j, k ) );
    }

    // Heating source term, host overload.
    KOKKOS_INLINE_FUNCTION
    auto source( HostTag, const int i, const int j, const int k ) const
    {
        // performance improvements on host: scoping the exponential
        if ( power_ )
        {
            double w = weight( i, j, k );

            if ( w < w_max_ )
            {
                return I0_ * power_ * Kokkos::exp( -w );
            }
            else
            {
                return 0.0;
            }
        }
        else
        {
            return 0.0;
        }
    }
};

// Create a solver based on the grid details and simulation inputs.
template <typename MemorySpace>
auto createSolver( Inputs db, Grid<MemorySpace> grid )
{
    using entity_type = typename Grid<MemorySpace>::entity_type;
    using view_type = typename Grid<MemorySpace>::view_type;
    using mesh_type = typename Grid<MemorySpace>::local_mesh_type;

    auto local_mesh = grid.getLocalMesh();

    return Solver<view_type, entity_type, mesh_type>( db, local_mesh );
}

} // namespace Finch

#endif
