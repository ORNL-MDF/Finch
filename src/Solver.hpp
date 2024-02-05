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

#ifndef Solver_H
#define Solver_H

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

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
    double r_[3];
    double A_inv_[3];
    double I0_;
    double w_max_;

  public:
    Solver( Simulation db, LocalMeshType local_mesh )
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

        // heat source parameter constants
        for ( std::size_t d = 0; d < 3; ++d )
        {
            r_[d] = db.source.two_sigma[d] / Kokkos::sqrt( 2.0 );
            A_inv_[d] = 1.0 / r_[d] / r_[d];
        }

        I0_ = ( 2.0 * db.source.absorption ) /
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

    // Normalized weight for the gaussian source term: x in exp(-x)
    KOKKOS_INLINE_FUNCTION
    auto weight( const int i, const int j, const int k ) const
    {
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
auto createSolver( Simulation db, Grid<MemorySpace> grid )
{
    using entity_type = typename Grid<MemorySpace>::entity_type;
    using view_type = typename Grid<MemorySpace>::view_type;
    using mesh_type = typename Grid<MemorySpace>::local_mesh_type;

    auto local_mesh = grid.getLocalMesh();

    return Solver<view_type, entity_type, mesh_type>( db, local_mesh );
}

#endif
