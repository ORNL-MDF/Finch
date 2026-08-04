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
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>
#include <variant>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_Grid.hpp>
#include <Finch_HeatSource.hpp>
#include <Finch_Inputs.hpp>
#include <Finch_Isotherm.hpp>
#include <Finch_Material.hpp>

namespace Finch
{

struct DiffusionTag
{
};
struct MaterialLookupTag
{
};
struct MaterialRefreshTag
{
};

template <typename ViewType, typename LocalMeshType, typename EntityType,
          typename SourceData, bool UsesLookupTable>
struct SourceApplyFunctor
{
    ViewType temperature;
    ViewType previous_temperature;
    ViewType inverse_heat_capacity;
    LocalMeshType local_mesh;
    SourceData source;
    double dt = 0.0;
    double solidus = 0.0;
    double liquidus = 0.0;
    double rho_cp = 0.0;
    double rho_latent_capacity = 0.0;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int i, const int j, const int k ) const
    {
        const int index[3] = { i, j, k };
        double point[3];
        local_mesh.coordinates( EntityType(), index, point );
        const double power_density = source.volumetricPower( point );
        if ( power_density == 0.0 )
            return;

        const double old_temperature = previous_temperature( i, j, k, 0 );
        double dt_by_heat_capacity;
        if constexpr ( UsesLookupTable )
            dt_by_heat_capacity = dt * inverse_heat_capacity( i, j, k, 0 );
        else
            dt_by_heat_capacity =
                ( old_temperature >= solidus && old_temperature <= liquidus )
                    ? dt / ( rho_cp + rho_latent_capacity )
                    : dt / rho_cp;
        temperature( i, j, k, 0 ) += power_density * dt_by_heat_capacity;
    }
};

template <typename LocalMeshType, typename EntityType, typename SourceData>
struct SourcePowerFunctor
{
    LocalMeshType local_mesh;
    SourceData source;
    Kokkos::Array<long, 3> owned_min;
    Kokkos::Array<long, 3> owned_max;
    Kokkos::Array<int, 3> on_low_boundary;
    Kokkos::Array<int, 3> on_high_boundary;
    double nodal_volume = 0.0;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int i, const int j, const int k,
                     double& integrated_power ) const
    {
        const int index[3] = { i, j, k };
        double point[3];
        local_mesh.coordinates( EntityType(), index, point );

        double quadrature_weight = 1.0;
        for ( int d = 0; d < 3; ++d )
        {
            if ( on_low_boundary[d] && index[d] == owned_min[d] )
                quadrature_weight *= 0.5;
            if ( on_high_boundary[d] && index[d] == owned_max[d] - 1 )
                quadrature_weight *= 0.5;
        }
        integrated_power +=
            quadrature_weight * nodal_volume * source.volumetricPower( point );
    }
};

template <typename LocalMeshType, typename EntityType, typename SourceData,
          typename OutputView>
struct SourceFieldFunctor
{
    LocalMeshType local_mesh;
    SourceData source;
    OutputView output;
    Kokkos::Array<long, 3> output_min;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int i, const int j, const int k ) const
    {
        const int index[3] = { i, j, k };
        double point[3];
        local_mesh.coordinates( EntityType(), index, point );
        output( i - output_min[0], j - output_min[1], k - output_min[2] ) =
            source.volumetricPower( point );
    }
};

template <typename ViewType, typename EntityType, typename LocalMeshType,
          typename MaterialModel>
class Solver
{
    using memory_space = typename ViewType::memory_space;
    using source_variant_type = HeatSourceVariant<memory_space>;

  protected:
    // temperature views are default constructed and updated every step.
    ViewType T_;
    ViewType T0_;

    LocalMeshType local_mesh_;
    MaterialModel material_;
    source_variant_type source_;
    ViewType thermal_potential_;
    ViewType inverse_heat_capacity_;
    Cabana::Grid::IndexSpace<3> ghosted_space_;
    bool material_state_initialized_ = false;

    // solution parameters
    double dt_;
    double solidus_;
    double liquidus_;
    double rho_cp_;
    double rho_Lf_by_dT_;
    double inverse_cell_size_squared_;
    double k_by_dx2_;
    double cell_size_;

  public:
    Solver( const Inputs& db, const LocalMeshType& local_mesh,
            MaterialModel material, source_variant_type source,
            ViewType thermal_potential = ViewType(),
            ViewType inverse_heat_capacity = ViewType(),
            Cabana::Grid::IndexSpace<3> ghosted_space =
                Cabana::Grid::IndexSpace<3>() )
        : local_mesh_( local_mesh )
        , material_( std::move( material ) )
        , source_( std::move( source ) )
        , thermal_potential_( thermal_potential )
        , inverse_heat_capacity_( inverse_heat_capacity )
        , ghosted_space_( ghosted_space )
    {
        // solution parameter constants
        double dx = db.space.cell_size;
        double rho = db.properties.density;
        double cp = db.properties.specific_heat.constantValue();
        double Lf = db.properties.latent_heat;

        dt_ = 0.0;
        solidus_ = db.properties.solidus;
        liquidus_ = db.properties.liquidus;
        rho_cp_ = rho * cp;
        rho_Lf_by_dT_ = rho * Lf / ( liquidus_ - solidus_ );
        inverse_cell_size_squared_ = 1.0 / ( dx * dx );
        k_by_dx2_ =
            db.properties.thermal_conductivity.constantValue() / ( dx * dx );
        cell_size_ = dx;
    }

    const MaterialModel& materialModel() const { return material_; }

    HeatSourceState heatSourceState() const
    {
        return std::visit( []( const auto& source ) { return source.state(); },
                           source_ );
    }

    // Integrate the exact discrete source used by the solver. This diagnostic
    // is called only at progress-monitor events, keeping reductions out of the
    // ordinary timestep path.
    template <typename ExecutionSpace, typename GridType, typename BeamType>
    double integratedSourcePower( const ExecutionSpace& execution_space,
                                  GridType& grid, const BeamType& beam ) const
    {
        if ( beam.power() <= 0.0 )
            return 0.0;

        Kokkos::Profiling::ScopedRegion region(
            "Finch::source_power_diagnostic" );
        const auto owned_space = grid.getIndexSpace();
        double local_power = 0.0;
        std::visit(
            [&]( const auto& source )
            {
                const auto source_space = sourceIndexSpace(
                    owned_space,
                    source.bounds( beam.position(), beam.direction() ) );
                if ( source_space.size() == 0 )
                    return;

                auto source_data = source.deviceData(
                    beam.power(), beam.position(), beam.direction() );
                SourcePowerFunctor<LocalMeshType, EntityType,
                                   decltype( source_data )>
                    functor;
                functor.local_mesh = local_mesh_;
                functor.source = source_data;
                functor.nodal_volume = cell_size_ * cell_size_ * cell_size_;
                for ( int d = 0; d < 3; ++d )
                {
                    functor.owned_min[d] = owned_space.min( d );
                    functor.owned_max[d] = owned_space.max( d );
                    functor.on_low_boundary[d] =
                        local_mesh_.onLowBoundary( d ) ? 1 : 0;
                    functor.on_high_boundary[d] =
                        local_mesh_.onHighBoundary( d ) ? 1 : 0;
                }

                Kokkos::Sum<double> reducer( local_power );
                Cabana::Grid::grid_parallel_reduce(
                    "Finch::integrated_source_power", execution_space,
                    source_space, functor, reducer );
                execution_space.fence( "Finch integrated source power" );
            },
            source_ );

        double global_power = 0.0;
        MPI_Reduce( &local_power, &global_power, 1, MPI_DOUBLE, MPI_SUM, 0,
                    grid.getComm() );
        return global_power;
    }

    // Materialize the instantaneous volumetric source only when requested by
    // field output. The dense output buffer is cleared once, then only the
    // compact source support is evaluated.
    template <typename ExecutionSpace, typename BeamType, typename OutputView>
    void
    fillVolumetricSourceField( const ExecutionSpace& execution_space,
                               const BeamType& beam,
                               const Cabana::Grid::IndexSpace<3>& output_space,
                               const OutputView& output ) const
    {
        Kokkos::deep_copy( execution_space, output, 0.0 );
        if ( beam.power() <= 0.0 )
            return;

        std::visit(
            [&]( const auto& source )
            {
                const auto source_space = sourceIndexSpace(
                    output_space,
                    source.bounds( beam.position(), beam.direction() ) );
                if ( source_space.size() == 0 )
                    return;

                auto source_data = source.deviceData(
                    beam.power(), beam.position(), beam.direction() );
                SourceFieldFunctor<LocalMeshType, EntityType,
                                   decltype( source_data ), OutputView>
                    functor;
                functor.local_mesh = local_mesh_;
                functor.source = source_data;
                functor.output = output;
                for ( int d = 0; d < 3; ++d )
                    functor.output_min[d] = output_space.min( d );
                Cabana::Grid::grid_parallel_for(
                    "Finch::volumetric_heat_source_output", execution_space,
                    source_space, functor );
            },
            source_ );
    }

    template <typename ExecutionSpace, typename GridType, typename BeamType>
    void prepareSource( const ExecutionSpace& execution_space, GridType& grid,
                        const BeamType& beam )
    {
        if ( beam.power() <= 0.0 )
            return;

        auto* source = std::get_if<TabulatedSource<memory_space>>( &source_ );
        if ( source == nullptr || !source->depthFeedbackEnabled() )
            return;

        const double detected_depth = findLocalMeltPoolDepth(
            execution_space, grid, beam.position(), source->searchRadius(),
            source->feedbackTemperature() );
        source->setDetectedDepth( detected_depth );
    }

    // Function for temperature solve: forward time-centered space (FTCS) method
    template <class ExecSpace, class IndexSpaceType>
    void solve( ExecSpace exec_space, IndexSpaceType owned_space, ViewType& T,
                ViewType& T0, const double dt, const double beam_power,
                const std::array<double, 3>& beam_pos,
                const std::array<double, 3>& beam_direction )
    {
        Kokkos::Profiling::ScopedRegion solve_region( "Finch::solve" );

        // Update temperature views and beam parameters for current time step
        T_ = T;

        T0_ = T0;

        dt_ = dt;

        if constexpr ( MaterialModel::uses_lookup_table )
        {
            if ( !material_state_initialized_ )
            {
                Cabana::Grid::grid_parallel_for(
                    "Finch::material_lookup_initial", exec_space,
                    ghosted_space_, MaterialLookupTag{}, *this );
                material_state_initialized_ = true;
            }
            else
                Cabana::Grid::grid_parallel_for(
                    "Finch::material_lookup_changed", exec_space,
                    ghosted_space_, MaterialRefreshTag{}, *this );
        }

        // Keep the numerical path identical on host and accelerator backends.
        // Splitting diffusion from the source avoids evaluating exp() over the
        // full domain and lets the source run only in its compact support.
        Cabana::Grid::grid_parallel_for( "Finch::diffusion", exec_space,
                                         owned_space, DiffusionTag{}, *this );

        if ( beam_power > 0.0 )
            std::visit(
                [&]( const auto& source )
                {
                    applySource( exec_space, owned_space, source, beam_power,
                                 beam_pos, beam_direction );
                },
                source_ );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( MaterialLookupTag, const int i, const int j,
                     const int k ) const
    {
        const auto state = material_.evaluate( T0_( i, j, k, 0 ) );
        thermal_potential_( i, j, k, 0 ) = state.thermal_potential;
        inverse_heat_capacity_( i, j, k, 0 ) =
            state.inverse_volumetric_heat_capacity;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( MaterialRefreshTag, const int i, const int j,
                     const int k ) const
    {
        // Before diffusion overwrites T, it still contains this node's value
        // from one step before T0. If the value did not change, the cached
        // material state is already exact and no lookup or write is needed.
        const double temperature = T0_( i, j, k, 0 );
        if ( temperature != T_( i, j, k, 0 ) )
        {
            const auto state = material_.evaluate( temperature );
            thermal_potential_( i, j, k, 0 ) = state.thermal_potential;
            inverse_heat_capacity_( i, j, k, 0 ) =
                state.inverse_volumetric_heat_capacity;
        }
    }

    // Explicit diffusion update.
    KOKKOS_INLINE_FUNCTION
    void operator()( DiffusionTag, const int i, const int j, const int k ) const
    {
        double x = T0_( i, j, k, 0 );

        if constexpr ( !MaterialModel::uses_lookup_table )
        {
            const double dt_by_rho_cp = ( x >= solidus_ && x <= liquidus_ )
                                            ? dt_ / ( rho_cp_ + rho_Lf_by_dT_ )
                                            : dt_ / rho_cp_;
            T_( i, j, k, 0 ) = x + laplacian( i, j, k ) * dt_by_rho_cp;
        }
        else
        {
            T_( i, j, k, 0 ) = x + conservativeDiffusion( i, j, k ) * dt_ *
                                       inverse_heat_capacity_( i, j, k, 0 );
        }
    }

    // A compile-time specialization keeps table lookup and face averaging out
    // of the constant-property kernel.
    KOKKOS_INLINE_FUNCTION
    auto laplacian( const int i, const int j, const int k ) const
    {
        return ( T0_( i - 1, j, k, 0 ) + T0_( i + 1, j, k, 0 ) +
                 T0_( i, j - 1, k, 0 ) + T0_( i, j + 1, k, 0 ) +
                 T0_( i, j, k - 1, 0 ) + T0_( i, j, k + 1, 0 ) -
                 6.0 * T0_( i, j, k, 0 ) ) *
               k_by_dx2_;
    }

    // The Kirchhoff thermal potential Phi(T) = integral(k(T) dT) turns
    // div(k grad(T)) into laplacian(Phi). Differences of Phi give an equal and
    // opposite flux on each shared face without six divisions per node.
    KOKKOS_INLINE_FUNCTION
    double conservativeDiffusion( const int i, const int j, const int k ) const
    {
        return inverse_cell_size_squared_ *
               ( thermal_potential_( i - 1, j, k, 0 ) +
                 thermal_potential_( i + 1, j, k, 0 ) +
                 thermal_potential_( i, j - 1, k, 0 ) +
                 thermal_potential_( i, j + 1, k, 0 ) +
                 thermal_potential_( i, j, k - 1, 0 ) +
                 thermal_potential_( i, j, k + 1, 0 ) -
                 6.0 * thermal_potential_( i, j, k, 0 ) );
    }

    template <class IndexSpaceType, class SourceType>
    void applySource( const typename memory_space::execution_space& exec_space,
                      const IndexSpaceType& owned_space,
                      const SourceType& source, const double beam_power,
                      const std::array<double, 3>& beam_position,
                      const std::array<double, 3>& beam_direction )
    {
        const auto source_space = sourceIndexSpace(
            owned_space, source.bounds( beam_position, beam_direction ) );
        if ( source_space.size() == 0 )
            return;

        auto source_data =
            source.deviceData( beam_power, beam_position, beam_direction );
        SourceApplyFunctor<ViewType, LocalMeshType, EntityType,
                           decltype( source_data ),
                           MaterialModel::uses_lookup_table>
            functor{ T_,           T0_,         inverse_heat_capacity_,
                     local_mesh_,  source_data, dt_,
                     solidus_,     liquidus_,   rho_cp_,
                     rho_Lf_by_dT_ };
        Cabana::Grid::grid_parallel_for( source.kernelLabel(), exec_space,
                                         source_space, functor );
    }

  private:
    template <class IndexSpaceType>
    Cabana::Grid::IndexSpace<3>
    sourceIndexSpace( const IndexSpaceType& owned_space,
                      const SourceBounds& bounds ) const
    {
        std::array<long, 3> source_min;
        std::array<long, 3> source_max;
        for ( int d = 0; d < 3; ++d )
        {
            const double ghost_low =
                local_mesh_.lowCorner( Cabana::Grid::Ghost(), d );

            const double low_index = ( bounds.low[d] - ghost_low ) / cell_size_;
            const double high_index =
                ( bounds.high[d] - ghost_low ) / cell_size_;
            constexpr double index_roundoff =
                32.0 * std::numeric_limits<double>::epsilon();
            const double low_tolerance =
                index_roundoff * std::max( 1.0, std::abs( low_index ) );
            const double high_tolerance =
                index_roundoff * std::max( 1.0, std::abs( high_index ) );

            source_min[d] = std::max<long>(
                owned_space.min( d ),
                static_cast<long>( std::ceil( low_index - low_tolerance ) ) );
            source_max[d] = std::min<long>(
                owned_space.max( d ),
                static_cast<long>( std::floor( high_index + high_tolerance ) ) +
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
    using constant_solver_type =
        Solver<view_type, entity_type, mesh_type, ConstantMaterial>;
    using table_material_type = UniformTableMaterial<MemorySpace>;
    using table_solver_type =
        Solver<view_type, entity_type, mesh_type, table_material_type>;
    using solver_variant =
        std::variant<constant_solver_type, table_solver_type>;

    auto local_mesh = grid.getLocalMesh();
    auto source = createHeatSource<MemorySpace>( db, grid.getComm(),
                                                 grid.executionSpace() );
    const auto initialize_boundaries = [&]( solver_variant solver )
    {
        std::visit(
            [&]( const auto& concrete_solver )
            { grid.initializeBoundaries( concrete_solver.materialModel() ); },
            solver );
        grid.gather();
        return solver;
    };
    if ( db.properties.isTemperatureDependent() )
    {
        table_material_type material( db, grid.executionSpace() );
        auto material_layout = Cabana::Grid::createArrayLayout(
            grid.getLocalGrid(), 1, entity_type() );
        auto thermal_potential = Cabana::Grid::createArray<double, MemorySpace>(
            "Finch thermal potential", material_layout );
        auto inverse_heat_capacity =
            Cabana::Grid::createArray<double, MemorySpace>(
                "Finch inverse heat capacity", material_layout );
        return initialize_boundaries( solver_variant(
            std::in_place_type<table_solver_type>, db, local_mesh,
            std::move( material ), std::move( source ),
            thermal_potential->view(), inverse_heat_capacity->view(),
            grid.getGhostedIndexSpace() ) );
    }

    ConstantMaterial material( db );
    return initialize_boundaries( solver_variant(
        std::in_place_type<constant_solver_type>, db, local_mesh,
        std::move( material ), std::move( source ) ) );
}

} // namespace Finch

#endif
