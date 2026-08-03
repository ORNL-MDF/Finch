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

#include <algorithm>
#include <array>
#include <cmath>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_Grid.hpp>
#include <Finch_Inputs.hpp>

namespace Finch
{

struct DiffusionTag
{
};
struct SourceTag
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
    double cell_size_;

    // heat source parameters
    double power_;
    double position_[3];
    double r_[3];
    double cutoff_radius_[3];
    double A_inv_[3];
    double I0_;
    double w_max_;

  public:
    Solver( const Inputs& db, const LocalMeshType& local_mesh )
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
        cell_size_ = dx;

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

        constexpr double pi = 3.141592653589793238462643383279502884;
        I0_ = ( 2.0 * db.source.absorption ) /
              ( pi * Kokkos::sqrt( pi ) * r_[0] * r_[1] * r_[2] );

        // Truncate contributions below approximately 0.1% of peak intensity.
        w_max_ = Kokkos::log( 3 ) + 2 * Kokkos::log( 10 );
        for ( std::size_t d = 0; d < 3; ++d )
            cutoff_radius_[d] = r_[d] * std::sqrt( w_max_ );
    }

    // Function for temperature solve: forward time-centered space (FTCS) method
    template <class ExecSpace, class IndexSpaceType>
    void solve( ExecSpace exec_space, IndexSpaceType owned_space, ViewType& T,
                ViewType& T0, const double dt, const double beam_power,
                const std::array<double, 3>& beam_pos )
    {
        Kokkos::Profiling::ScopedRegion solve_region( "Finch::solve" );

        // Update temperature views and beam parameters for current time step
        T_ = T;

        T0_ = T0;

        dt_ = dt;
        power_ = beam_power;

        for ( std::size_t d = 0; d < 3; ++d )
        {
            position_[d] = beam_pos[d];
        }

        // Keep the numerical path identical on host and accelerator backends.
        // Splitting diffusion from the source avoids evaluating exp() over the
        // full domain and lets the source run only in its compact support.
        Cabana::Grid::grid_parallel_for( "Finch::diffusion", exec_space,
                                         owned_space, DiffusionTag{}, *this );

        if ( power_ > 0.0 && I0_ > 0.0 )
        {
            auto source_space = sourceIndexSpace( owned_space );
            if ( source_space.size() > 0 )
                Cabana::Grid::grid_parallel_for( "Finch::gaussian_source",
                                                 exec_space, source_space,
                                                 SourceTag{}, *this );
        }
    }

    // Explicit diffusion update.
    KOKKOS_INLINE_FUNCTION
    void operator()( DiffusionTag, const int i, const int j, const int k ) const
    {
        double x = T0_( i, j, k, 0 );

        double dt_by_rho_cp = ( x >= solidus_ && x <= liquidus_ )
                                  ? dt_ / ( rho_cp_ + rho_Lf_by_dT_ )
                                  : dt_ / ( rho_cp_ );

        T_( i, j, k, 0 ) = x + laplacian( i, j, k ) * dt_by_rho_cp;
    }

    // Add the compact Gaussian source to the completed diffusion update.
    KOKKOS_INLINE_FUNCTION
    void operator()( SourceTag, const int i, const int j, const int k ) const
    {
        double x = T0_( i, j, k, 0 );

        double dt_by_rho_cp = ( x >= solidus_ && x <= liquidus_ )
                                  ? dt_ / ( rho_cp_ + rho_Lf_by_dT_ )
                                  : dt_ / rho_cp_;

        const double w = weight( i, j, k );
        if ( w < w_max_ )
            T_( i, j, k, 0 ) += I0_ * power_ * Kokkos::exp( -w ) * dt_by_rho_cp;
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

    template <class IndexSpaceType>
    Cabana::Grid::IndexSpace<3>
    sourceIndexSpace( const IndexSpaceType& owned_space ) const
    {
        std::array<long, 3> source_min;
        std::array<long, 3> source_max;
        for ( int d = 0; d < 3; ++d )
        {
            const double ghost_low =
                local_mesh_.lowCorner( Cabana::Grid::Ghost(), d );

            source_min[d] = std::max<long>(
                owned_space.min( d ),
                static_cast<long>( std::ceil(
                    ( position_[d] - cutoff_radius_[d] - ghost_low ) /
                    cell_size_ ) ) );
            source_max[d] = std::min<long>(
                owned_space.max( d ),
                static_cast<long>( std::floor(
                    ( position_[d] + cutoff_radius_[d] - ghost_low ) /
                    cell_size_ ) ) +
                    1 );
            source_max[d] = std::max( source_max[d], source_min[d] );
        }
        return Cabana::Grid::IndexSpace<3>( source_min, source_max );
    }
};

// Create a solver based on the grid details and simulation inputs.
template <typename MemorySpace>
auto createSolver( const Inputs& db, Grid<MemorySpace>& grid )
{
    using entity_type = typename Grid<MemorySpace>::entity_type;
    using view_type = typename Grid<MemorySpace>::view_type;
    using mesh_type = typename Grid<MemorySpace>::local_mesh_type;

    auto local_mesh = grid.getLocalMesh();

    return Solver<view_type, entity_type, mesh_type>( db, local_mesh );
}

} // namespace Finch

#endif
