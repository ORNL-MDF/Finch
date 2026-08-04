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

/*!
  \file Finch_Material.hpp
  \brief Device-callable constant and tabulated thermal material models.
*/

#ifndef FINCH_MATERIAL_HPP
#define FINCH_MATERIAL_HPP

#include <algorithm>
#include <limits>

#include <Kokkos_Core.hpp>

#include <Finch_Inputs.hpp>

namespace Finch
{

struct MaterialState
{
    double thermal_potential;
    double inverse_volumetric_heat_capacity;
};

class ConstantMaterial
{
  public:
    static constexpr bool uses_lookup_table = false;

    explicit ConstantMaterial( const Inputs& inputs )
        : conductivity_(
              inputs.properties.thermal_conductivity.constantValue() )
    {
    }

    KOKKOS_INLINE_FUNCTION
    double conductivity( const double ) const { return conductivity_; }

  private:
    double conductivity_ = 0.0;
};

template <typename MemorySpace>
class UniformTableMaterial
{
  public:
    static constexpr bool uses_lookup_table = true;
    using memory_space = MemorySpace;
    // Use the backend's default layout: row-major on the current CPU build
    // and a coalescing-friendly layout on accelerator memory spaces.
    using table_type = Kokkos::View<double*[4], memory_space>;

    template <typename ExecutionSpace>
    UniformTableMaterial( const Inputs& inputs,
                          const ExecutionSpace& execution_space )
        : solidus_( inputs.properties.solidus )
        , liquidus_( inputs.properties.liquidus )
        , table_size_( inputs.properties.lookup_table_points )
    {
        const double density = inputs.properties.density;
        const double volumetric_latent_capacity =
            density * inputs.properties.latent_heat /
            ( inputs.properties.liquidus - inputs.properties.solidus );
        double minimum_temperature = std::numeric_limits<double>::max();
        double maximum_temperature = std::numeric_limits<double>::lowest();
        const auto include_range = [&]( const TemperatureProperty& property ) {
            if ( property.isTabulated() )
            {
                minimum_temperature =
                    std::min( minimum_temperature,
                              property.temperature.front() );
                maximum_temperature =
                    std::max( maximum_temperature,
                              property.temperature.back() );
            }
        };
        include_range( inputs.properties.specific_heat );
        include_range( inputs.properties.thermal_conductivity );

        minimum_temperature_ = minimum_temperature;
        maximum_temperature_ = maximum_temperature;
        temperature_spacing_ =
            ( maximum_temperature - minimum_temperature ) /
            static_cast<double>( table_size_ - 1 );
        inverse_temperature_spacing_ = 1.0 / temperature_spacing_;

        table_ = table_type( Kokkos::ViewAllocateWithoutInitializing(
                                 "Finch material lookup table" ),
                             table_size_ );
        auto host_table = Kokkos::create_mirror_view( table_ );
        for ( int i = 0; i < table_size_; ++i )
        {
            const double temperature =
                minimum_temperature +
                static_cast<double>( i ) * temperature_spacing_;
            const double specific_heat =
                inputs.properties.specific_heat.value( temperature );
            host_table( i, 0 ) = 1.0 / ( density * specific_heat );
            host_table( i, 1 ) =
                1.0 / ( density * specific_heat +
                        volumetric_latent_capacity );
            host_table( i, 2 ) =
                inputs.properties.thermal_conductivity.value( temperature );
            if ( i == 0 )
                host_table( i, 3 ) = 0.0;
            else
                host_table( i, 3 ) =
                    host_table( i - 1, 3 ) +
                    0.5 * temperature_spacing_ *
                        ( host_table( i - 1, 2 ) + host_table( i, 2 ) );
        }
        Kokkos::deep_copy( execution_space, table_, host_table );
        execution_space.fence( "Finch material table initialization" );
    }

    KOKKOS_INLINE_FUNCTION
    MaterialState evaluate( const double temperature ) const
    {
        int lower;
        double fraction;
        interval( temperature, lower, fraction );
        const double thermal_potential =
            potential( temperature, lower, fraction );
        const int capacity_component =
            ( temperature >= solidus_ && temperature <= liquidus_ ) ? 1 : 0;
        return { thermal_potential,
                 interpolate( lower, fraction, capacity_component ) };
    }

    KOKKOS_INLINE_FUNCTION
    double conductivity( const double temperature ) const
    {
        int lower;
        double fraction;
        interval( temperature, lower, fraction );
        return interpolate( lower, fraction, 2 );
    }

  private:
    KOKKOS_INLINE_FUNCTION
    void interval( const double temperature, int& lower,
                   double& fraction ) const
    {
        double location = ( temperature - minimum_temperature_ ) *
                          inverse_temperature_spacing_;
        if ( location <= 0.0 )
        {
            lower = 0;
            fraction = 0.0;
        }
        else if ( location >= static_cast<double>( table_size_ - 1 ) )
        {
            lower = table_size_ - 2;
            fraction = 1.0;
        }
        else
        {
            lower = static_cast<int>( location );
            fraction = location - static_cast<double>( lower );
        }
    }

    KOKKOS_INLINE_FUNCTION
    double interpolate( const int lower, const double fraction,
                        const int component ) const
    {
        const double low = table_( lower, component );
        return low + fraction * ( table_( lower + 1, component ) - low );
    }

    KOKKOS_INLINE_FUNCTION
    double potential( const double temperature, const int lower,
                      const double fraction ) const
    {
        if ( temperature < minimum_temperature_ )
            return table_( 0, 3 ) + table_( 0, 2 ) *
                                         ( temperature - minimum_temperature_ );

        if ( temperature > maximum_temperature_ )
            return table_( table_size_ - 1, 3 ) +
                   table_( table_size_ - 1, 2 ) *
                       ( temperature - maximum_temperature_ );

        // Conductivity is linear inside a lookup interval, so its integral is
        // quadratic. This evaluates Phi(T) exactly for the resampled table;
        // linearly interpolating the cumulative integral would instead use a
        // piecewise-constant effective conductivity.
        const double conductivity_low = table_( lower, 2 );
        const double conductivity_delta =
            table_( lower + 1, 2 ) - conductivity_low;
        return table_( lower, 3 ) +
               temperature_spacing_ *
                   ( conductivity_low * fraction +
                     0.5 * conductivity_delta * fraction * fraction );
    }

    table_type table_;
    double solidus_ = 0.0;
    double liquidus_ = 0.0;
    double minimum_temperature_ = 0.0;
    double maximum_temperature_ = 0.0;
    double temperature_spacing_ = 0.0;
    double inverse_temperature_spacing_ = 0.0;
    int table_size_ = 0;
};

} // namespace Finch

#endif
