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
  \file Finch_HeatSource.hpp
  \brief Host-selected, device-callable heat-source models.
*/

#ifndef FINCH_HEAT_SOURCE_HPP
#define FINCH_HEAT_SOURCE_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <Finch_Inputs.hpp>

namespace Finch
{

struct SourceBounds
{
    std::array<double, 3> low = { 0.0, 0.0, 0.0 };
    std::array<double, 3> high = { 0.0, 0.0, 0.0 };
};

struct HeatSourceState
{
    double absorptivity = 0.0;
    double measured_melt_depth = 0.0;
    double effective_source_depth = 0.0;
    double lateral_d4_sigma = 0.0;
    double lateral_d4_sigma_major = 0.0;
    double lateral_d4_sigma_minor = 0.0;
    double profile_azimuth = 0.0;
    double aspect_ratio = 0.0;
    bool depth_feedback = false;
    bool dynamic_absorption = false;
};

class ConstantAbsorption
{
  public:
    explicit ConstantAbsorption( const AbsorptionInput& input )
        : coefficient_( input.coefficient )
    {
    }

    double evaluate( const double ) const { return coefficient_; }

  private:
    double coefficient_ = 0.0;
};

class KellyAbsorption
{
  public:
    explicit KellyAbsorption( const AbsorptionInput& input )
        : fresnel_absorptivity_( input.fresnel_absorptivity )
        , conduction_absorptivity_( input.conduction_absorptivity )
        , transition_aspect_ratio_( input.transition_aspect_ratio )
        , cone_( input.geometry == "cone" )
    {
    }

    double evaluate( const double aspect_ratio ) const
    {
        if ( aspect_ratio <= transition_aspect_ratio_ )
            return conduction_absorptivity_;

        const double aspect_ratio_squared = aspect_ratio * aspect_ratio;
        double geometry_factor_f;
        double geometry_factor_g;
        if ( cone_ )
        {
            const double radial_factor =
                std::sqrt( 1.0 + aspect_ratio_squared );
            geometry_factor_f =
                1.0 / ( radial_factor * radial_factor * radial_factor );
            geometry_factor_g = 1.0 / ( 1.0 + radial_factor );
        }
        else
        {
            geometry_factor_f = 1.0 / ( 1.0 + aspect_ratio_squared );
            geometry_factor_g = 0.5 / ( 1.0 + aspect_ratio );
        }

        const double reflectivity = 1.0 - fresnel_absorptivity_;
        return fresnel_absorptivity_ *
               ( 1.0 + reflectivity *
                           ( geometry_factor_g - geometry_factor_f ) ) /
               ( 1.0 - reflectivity * ( 1.0 - geometry_factor_g ) );
    }

  private:
    double fresnel_absorptivity_ = 0.0;
    double conduction_absorptivity_ = 0.0;
    double transition_aspect_ratio_ = 1.0;
    bool cone_ = true;
};

using AbsorptionVariant = std::variant<ConstantAbsorption, KellyAbsorption>;

inline AbsorptionVariant createAbsorption( const AbsorptionInput& input )
{
    if ( input.type == "kelly" )
        return AbsorptionVariant( std::in_place_type<KellyAbsorption>, input );
    return AbsorptionVariant( std::in_place_type<ConstantAbsorption>, input );
}

struct GaussianSourceData
{
    Kokkos::Array<double, 3> position;
    Kokkos::Array<double, 3> inverse_radius_squared;
    double intensity = 0.0;
    double maximum_weight = 0.0;

    KOKKOS_INLINE_FUNCTION
    double volumetricPower( const double point[3] ) const
    {
        double weight = 0.0;
        for ( int d = 0; d < 3; ++d )
        {
            const double distance = point[d] - position[d];
            weight += distance * distance * inverse_radius_squared[d];
        }
        return weight < maximum_weight ? intensity * Kokkos::exp( -weight )
                                       : 0.0;
    }
};

class GaussianSource
{
  public:
    static constexpr bool supports_transient_depth = false;

    explicit GaussianSource( const Inputs& inputs )
    {
        constexpr double pi = 3.141592653589793238462643383279502884;
        maximum_weight_ = Kokkos::log( 3.0 ) + 2.0 * Kokkos::log( 10.0 );
        double radius_product = 1.0;
        for ( int d = 0; d < 3; ++d )
        {
            radius_[d] = inputs.source.two_sigma[d] / Kokkos::sqrt( 2.0 );
            inverse_radius_squared_[d] = 1.0 / ( radius_[d] * radius_[d] );
            cutoff_radius_[d] = radius_[d] * std::sqrt( maximum_weight_ );
            radius_product *= radius_[d];
        }
        absorptivity_ = inputs.source.absorption.coefficient;
        lateral_d4_sigma_major_ =
            2.0 * std::max( inputs.source.two_sigma[0],
                            inputs.source.two_sigma[1] );
        lateral_d4_sigma_minor_ =
            2.0 * std::min( inputs.source.two_sigma[0],
                            inputs.source.two_sigma[1] );
        lateral_d4_sigma_ = std::sqrt( lateral_d4_sigma_major_ *
                                       lateral_d4_sigma_minor_ );
        effective_source_depth_ = inputs.source.two_sigma[2];
        intensity_per_watt_ = ( 2.0 * absorptivity_ ) /
                              ( pi * std::sqrt( pi ) * radius_product );
    }

    SourceBounds bounds( const std::array<double, 3>& position,
                         const std::array<double, 3>& ) const
    {
        SourceBounds result;
        for ( int d = 0; d < 3; ++d )
        {
            result.low[d] = position[d] - cutoff_radius_[d];
            result.high[d] = position[d] + cutoff_radius_[d];
        }
        return result;
    }

    GaussianSourceData deviceData( const double power,
                                   const std::array<double, 3>& position,
                                   const std::array<double, 3>& ) const
    {
        GaussianSourceData data;
        for ( int d = 0; d < 3; ++d )
        {
            data.position[d] = position[d];
            data.inverse_radius_squared[d] = inverse_radius_squared_[d];
        }
        data.intensity = intensity_per_watt_ * power;
        data.maximum_weight = maximum_weight_;
        return data;
    }

    const char* kernelLabel() const { return "Finch::gaussian_source"; }

    HeatSourceState state() const
    {
        HeatSourceState result;
        result.absorptivity = absorptivity_;
        result.effective_source_depth = effective_source_depth_;
        result.lateral_d4_sigma = lateral_d4_sigma_;
        result.lateral_d4_sigma_major = lateral_d4_sigma_major_;
        result.lateral_d4_sigma_minor = lateral_d4_sigma_minor_;
        return result;
    }

  private:
    std::array<double, 3> radius_ = { 0.0, 0.0, 0.0 };
    std::array<double, 3> inverse_radius_squared_ = { 0.0, 0.0, 0.0 };
    std::array<double, 3> cutoff_radius_ = { 0.0, 0.0, 0.0 };
    double intensity_per_watt_ = 0.0;
    double maximum_weight_ = 0.0;
    double absorptivity_ = 0.0;
    double effective_source_depth_ = 0.0;
    double lateral_d4_sigma_ = 0.0;
    double lateral_d4_sigma_major_ = 0.0;
    double lateral_d4_sigma_minor_ = 0.0;
};

template <typename MemorySpace>
struct TabulatedSourceData
{
    using table_type = Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

    table_type table;
    Kokkos::Array<double, 3> position;
    Kokkos::Array<double, 2> scan_direction;
    double x0 = 0.0;
    double y0 = 0.0;
    double inverse_dx = 0.0;
    double inverse_dy = 0.0;
    double inverse_depth = 0.0;
    double cutoff_ratio = 0.0;
    double power_scale = 0.0;
    double axial_exponent = 0.0;
    int fast_exponent = 0;
    int nx = 0;
    int ny = 0;
    int scan_path_frame = 0;

    KOKKOS_INLINE_FUNCTION
    double volumetricPower( const double point[3] ) const
    {
        const double global_x = point[0] - position[0];
        const double global_y = point[1] - position[1];
        double profile_x = global_x;
        double profile_y = global_y;
        if ( scan_path_frame )
        {
            profile_x =
                global_x * scan_direction[0] + global_y * scan_direction[1];
            profile_y =
                -global_x * scan_direction[1] + global_y * scan_direction[0];
        }

        const double table_x = ( profile_x - x0 ) * inverse_dx;
        const double table_y = ( profile_y - y0 ) * inverse_dy;
        if ( table_x < 0.0 || table_x > static_cast<double>( nx - 1 ) ||
             table_y < 0.0 || table_y > static_cast<double>( ny - 1 ) )
            return 0.0;

        // The compact source bounds restrict evaluation to nodes at or below
        // the laser centroid. Clamp the signed centroid-to-node depth so
        // roundoff in a surface coordinate (for example, z = +5e-20 instead
        // of zero) does not incorrectly remove heating from the top node.
        const double normalized_depth = Kokkos::fmax(
            ( position[2] - point[2] ) * inverse_depth, 0.0 );
        if ( normalized_depth > cutoff_ratio )
            return 0.0;

        int i = static_cast<int>( table_x );
        int j = static_cast<int>( table_y );
        i = i < nx - 1 ? i : nx - 2;
        j = j < ny - 1 ? j : ny - 2;
        const double tx = table_x - static_cast<double>( i );
        const double ty = table_y - static_cast<double>( j );

        const double f00 = table( j, i );
        const double f10 = table( j, i + 1 );
        const double f01 = table( j + 1, i );
        const double f11 = table( j + 1, i + 1 );
        const double planar = ( 1.0 - tx ) * ( ( 1.0 - ty ) * f00 + ty * f01 ) +
                              tx * ( ( 1.0 - ty ) * f10 + ty * f11 );

        double axial_power;
        if ( fast_exponent == 1 )
            axial_power = normalized_depth;
        else if ( fast_exponent == 2 )
            axial_power = normalized_depth * normalized_depth;
        else if ( fast_exponent == 4 )
        {
            const double squared = normalized_depth * normalized_depth;
            axial_power = squared * squared;
        }
        else if ( fast_exponent == 8 )
        {
            const double squared = normalized_depth * normalized_depth;
            const double fourth = squared * squared;
            axial_power = fourth * fourth;
        }
        else
            axial_power = Kokkos::pow( normalized_depth, axial_exponent );

        return power_scale * planar * Kokkos::exp( -3.0 * axial_power );
    }
};

template <typename MemorySpace>
class TabulatedSource
{
  public:
    static constexpr bool supports_transient_depth = true;
    using memory_space = MemorySpace;
    using table_type = typename TabulatedSourceData<memory_space>::table_type;

    template <typename ExecutionSpace>
    TabulatedSource( const Inputs& inputs, MPI_Comm comm,
                     const ExecutionSpace& execution_space )
        : absorption_( createAbsorption( inputs.source.absorption ) )
        , minimum_depth_( inputs.source.minimum_depth )
        , current_depth_( inputs.source.minimum_depth )
        , exponent_slope_( inputs.source.exponent_slope )
        , exponent_intercept_( inputs.source.exponent_intercept )
        , feedback_temperature_( inputs.source.depth_temperature == "solidus"
                                     ? inputs.properties.solidus
                                     : inputs.properties.liquidus )
        , scan_path_frame_( inputs.source.coordinate_frame == "scan_path" )
        , depth_feedback_enabled_( inputs.source.transient_depth )
        , dynamic_absorption_( inputs.source.absorption.type == "kelly" )
    {
        readTable( inputs.source.profile_file, comm, execution_space );
        updateAxialState();
    }

    double feedbackTemperature() const { return feedback_temperature_; }
    double searchRadius() const { return search_radius_; }
    bool depthFeedbackEnabled() const { return depth_feedback_enabled_; }

    void setDetectedDepth( const double depth )
    {
        const double measured_depth = std::max( 0.0, depth );
        const double current_depth = std::max( minimum_depth_, measured_depth );
        const bool source_state_changed = current_depth != current_depth_;
        const bool absorption_state_changed =
            dynamic_absorption_ && measured_depth != measured_depth_;
        measured_depth_ = measured_depth;
        current_depth_ = current_depth;
        absorption_aspect_ratio_ =
            2.0 * measured_depth_ / lateral_d4_sigma_;
        if ( source_state_changed )
        {
            updateAxialState();
            return;
        }
        if ( absorption_state_changed )
            updateAbsorptionState();
    }

    SourceBounds bounds( const std::array<double, 3>& position,
                         const std::array<double, 3>& direction ) const
    {
        SourceBounds result;
        if ( !scan_path_frame_ )
        {
            result.low[0] = position[0] + x0_;
            result.high[0] = position[0] + x1_;
            result.low[1] = position[1] + y0_;
            result.high[1] = position[1] + y1_;
        }
        else
        {
            result.low[0] = result.low[1] = std::numeric_limits<double>::max();
            result.high[0] = result.high[1] =
                std::numeric_limits<double>::lowest();
            const double profile_x[2] = { x0_, x1_ };
            const double profile_y[2] = { y0_, y1_ };
            for ( int ix = 0; ix < 2; ++ix )
                for ( int iy = 0; iy < 2; ++iy )
                {
                    const double x = position[0] +
                                     profile_x[ix] * direction[0] -
                                     profile_y[iy] * direction[1];
                    const double y = position[1] +
                                     profile_x[ix] * direction[1] +
                                     profile_y[iy] * direction[0];
                    result.low[0] = std::min( result.low[0], x );
                    result.high[0] = std::max( result.high[0], x );
                    result.low[1] = std::min( result.low[1], y );
                    result.high[1] = std::max( result.high[1], y );
                }
        }
        result.low[2] = position[2] - cutoff_depth_;
        result.high[2] = position[2];
        return result;
    }

    TabulatedSourceData<memory_space>
    deviceData( const double power, const std::array<double, 3>& position,
                const std::array<double, 3>& direction ) const
    {
        TabulatedSourceData<memory_space> data;
        data.table = table_;
        for ( int d = 0; d < 3; ++d )
            data.position[d] = position[d];
        data.scan_direction[0] = direction[0];
        data.scan_direction[1] = direction[1];
        data.x0 = x0_;
        data.y0 = y0_;
        data.inverse_dx = inverse_dx_;
        data.inverse_dy = inverse_dy_;
        data.inverse_depth = 1.0 / current_depth_;
        data.cutoff_ratio = cutoff_depth_ / current_depth_;
        data.power_scale = power * normalization_per_watt_;
        data.axial_exponent = axial_exponent_;
        data.fast_exponent = fast_exponent_;
        data.nx = nx_;
        data.ny = ny_;
        data.scan_path_frame = scan_path_frame_ ? 1 : 0;
        return data;
    }

    const char* kernelLabel() const { return "Finch::tabulated_source"; }

    HeatSourceState state() const
    {
        HeatSourceState result;
        result.absorptivity = current_absorptivity_;
        result.measured_melt_depth = measured_depth_;
        result.effective_source_depth = current_depth_;
        result.lateral_d4_sigma = lateral_d4_sigma_;
        result.lateral_d4_sigma_major = d4_sigma_major_;
        result.lateral_d4_sigma_minor = d4_sigma_minor_;
        result.profile_azimuth = profile_azimuth_;
        result.aspect_ratio = absorption_aspect_ratio_;
        result.depth_feedback = depth_feedback_enabled_;
        result.dynamic_absorption = dynamic_absorption_;
        return result;
    }

  private:
    template <typename ExecutionSpace>
    void readTable( const std::string& filename, MPI_Comm comm,
                    const ExecutionSpace& execution_space )
    {
        int rank = 0;
        MPI_Comm_rank( comm, &rank );
        std::string error;
        std::vector<double> values;

        if ( rank == 0 )
        {
            try
            {
                std::ifstream input( filename );
                if ( !input )
                    throw std::runtime_error(
                        "Cannot open tabulated source profile " + filename );
                if ( !( input >> nx_ >> ny_ ) || nx_ < 2 || ny_ < 2 )
                    throw std::runtime_error(
                        "Tabulated source profile dimensions must be at "
                        "least 2 by 2" );
                if ( !( input >> x0_ >> y0_ >> dx_ >> dy_ ) ||
                     !std::isfinite( x0_ ) || !std::isfinite( y0_ ) ||
                     !std::isfinite( dx_ ) || !std::isfinite( dy_ ) ||
                     dx_ <= 0.0 || dy_ <= 0.0 )
                    throw std::runtime_error(
                        "Tabulated source origin and positive spacing are "
                        "required" );
                const long long table_size =
                    static_cast<long long>( nx_ ) * ny_;
                if ( table_size > INT_MAX )
                    throw std::runtime_error(
                        "Tabulated source profile is too large" );
                values.resize( static_cast<std::size_t>( table_size ) );
                for ( double& value : values )
                {
                    if ( !( input >> value ) || !std::isfinite( value ) ||
                         value < 0.0 )
                        throw std::runtime_error(
                            "Tabulated source values must be finite and "
                            "nonnegative" );
                }
                std::string trailing;
                if ( input >> trailing )
                    throw std::runtime_error(
                        "Unexpected trailing data in tabulated source "
                        "profile" );
                cropZeroBorder( values );
                integrateTable( values );
            }
            catch ( const std::exception& exception )
            {
                error = exception.what();
            }
        }

        int error_size = static_cast<int>( error.size() );
        MPI_Bcast( &error_size, 1, MPI_INT, 0, comm );
        if ( error_size > 0 )
        {
            error.resize( error_size );
            MPI_Bcast( error.data(), error_size, MPI_CHAR, 0, comm );
            throw std::runtime_error( error );
        }

        int dimensions[2] = { nx_, ny_ };
        double metadata[8] = { x0_,
                               y0_,
                               dx_,
                               dy_,
                               planar_integral_,
                               d4_sigma_major_,
                               d4_sigma_minor_,
                               profile_azimuth_ };
        MPI_Bcast( dimensions, 2, MPI_INT, 0, comm );
        MPI_Bcast( metadata, 8, MPI_DOUBLE, 0, comm );
        nx_ = dimensions[0];
        ny_ = dimensions[1];
        x0_ = metadata[0];
        y0_ = metadata[1];
        dx_ = metadata[2];
        dy_ = metadata[3];
        planar_integral_ = metadata[4];
        d4_sigma_major_ = metadata[5];
        d4_sigma_minor_ = metadata[6];
        profile_azimuth_ = metadata[7];
        lateral_d4_sigma_ =
            std::sqrt( d4_sigma_major_ * d4_sigma_minor_ );

        const int table_size = nx_ * ny_;
        if ( rank != 0 )
            values.resize( table_size );
        MPI_Bcast( values.data(), table_size, MPI_DOUBLE, 0, comm );

        inverse_dx_ = 1.0 / dx_;
        inverse_dy_ = 1.0 / dy_;
        x1_ = x0_ + static_cast<double>( nx_ - 1 ) * dx_;
        y1_ = y0_ + static_cast<double>( ny_ - 1 ) * dy_;
        search_radius_ =
            std::max( std::max( std::abs( x0_ ), std::abs( x1_ ) ),
                      std::max( std::abs( y0_ ), std::abs( y1_ ) ) );

        table_ = table_type( Kokkos::ViewAllocateWithoutInitializing(
                                 "Finch tabulated heat-source profile" ),
                             ny_, nx_ );
        auto host_table = Kokkos::create_mirror_view( table_ );
        for ( int j = 0; j < ny_; ++j )
            for ( int i = 0; i < nx_; ++i )
                host_table( j, i ) = values[i + nx_ * j];
        Kokkos::deep_copy( execution_space, table_, host_table );
        execution_space.fence( "Finch tabulated source initialization" );
    }

    void cropZeroBorder( std::vector<double>& values )
    {
        int low_i = nx_;
        int high_i = -1;
        int low_j = ny_;
        int high_j = -1;
        for ( int j = 0; j < ny_; ++j )
            for ( int i = 0; i < nx_; ++i )
                if ( values[i + nx_ * j] > 0.0 )
                {
                    low_i = std::min( low_i, i );
                    high_i = std::max( high_i, i );
                    low_j = std::min( low_j, j );
                    high_j = std::max( high_j, j );
                }
        if ( high_i < 0 )
            throw std::runtime_error(
                "Tabulated source profile must contain positive intensity" );

        low_i = std::max( 0, low_i - 1 );
        high_i = std::min( nx_ - 1, high_i + 1 );
        low_j = std::max( 0, low_j - 1 );
        high_j = std::min( ny_ - 1, high_j + 1 );
        if ( low_i == 0 && high_i == nx_ - 1 && low_j == 0 &&
             high_j == ny_ - 1 )
            return;

        const int old_nx = nx_;
        const int cropped_nx = high_i - low_i + 1;
        const int cropped_ny = high_j - low_j + 1;
        std::vector<double> cropped( cropped_nx * cropped_ny );
        for ( int j = 0; j < cropped_ny; ++j )
            for ( int i = 0; i < cropped_nx; ++i )
                cropped[i + cropped_nx * j] =
                    values[( i + low_i ) + old_nx * ( j + low_j )];
        x0_ += static_cast<double>( low_i ) * dx_;
        y0_ += static_cast<double>( low_j ) * dy_;
        nx_ = cropped_nx;
        ny_ = cropped_ny;
        values = std::move( cropped );
    }

    void integrateTable( const std::vector<double>& values )
    {
        planar_integral_ = 0.0;
        double first_moment_x = 0.0;
        double first_moment_y = 0.0;
        double second_moment_x = 0.0;
        double second_moment_y = 0.0;
        double cross_moment = 0.0;
        for ( int j = 0; j < ny_ - 1; ++j )
            for ( int i = 0; i < nx_ - 1; ++i )
            {
                const int index = i + nx_ * j;
                const double cell_x = x0_ + static_cast<double>( i ) * dx_;
                const double cell_y = y0_ + static_cast<double>( j ) * dy_;
                const double x_weights[3][2] = {
                    { 0.5 * dx_, 0.5 * dx_ },
                    { dx_ * ( 0.5 * cell_x + dx_ / 6.0 ),
                      dx_ * ( 0.5 * cell_x + dx_ / 3.0 ) },
                    { dx_ * ( 0.5 * cell_x * cell_x +
                              cell_x * dx_ / 3.0 + dx_ * dx_ / 12.0 ),
                      dx_ * ( 0.5 * cell_x * cell_x +
                              2.0 * cell_x * dx_ / 3.0 + dx_ * dx_ / 4.0 ) } };
                const double y_weights[3][2] = {
                    { 0.5 * dy_, 0.5 * dy_ },
                    { dy_ * ( 0.5 * cell_y + dy_ / 6.0 ),
                      dy_ * ( 0.5 * cell_y + dy_ / 3.0 ) },
                    { dy_ * ( 0.5 * cell_y * cell_y +
                              cell_y * dy_ / 3.0 + dy_ * dy_ / 12.0 ),
                      dy_ * ( 0.5 * cell_y * cell_y +
                              2.0 * cell_y * dy_ / 3.0 + dy_ * dy_ / 4.0 ) } };
                const double cell_values[2][2] = {
                    { values[index], values[index + nx_] },
                    { values[index + 1], values[index + nx_ + 1] } };
                for ( int x_node = 0; x_node < 2; ++x_node )
                    for ( int y_node = 0; y_node < 2; ++y_node )
                    {
                        const double value = cell_values[x_node][y_node];
                        planar_integral_ += value * x_weights[0][x_node] *
                                            y_weights[0][y_node];
                        first_moment_x += value * x_weights[1][x_node] *
                                          y_weights[0][y_node];
                        first_moment_y += value * x_weights[0][x_node] *
                                          y_weights[1][y_node];
                        second_moment_x += value * x_weights[2][x_node] *
                                           y_weights[0][y_node];
                        second_moment_y += value * x_weights[0][x_node] *
                                           y_weights[2][y_node];
                        cross_moment += value * x_weights[1][x_node] *
                                        y_weights[1][y_node];
                    }
            }
        if ( !std::isfinite( planar_integral_ ) || planar_integral_ <= 0.0 )
            throw std::runtime_error(
                "Tabulated source profile integral must be positive" );

        const double centroid_x = first_moment_x / planar_integral_;
        const double centroid_y = first_moment_y / planar_integral_;
        const double variance_x =
            second_moment_x / planar_integral_ - centroid_x * centroid_x;
        const double variance_y =
            second_moment_y / planar_integral_ - centroid_y * centroid_y;
        const double covariance =
            cross_moment / planar_integral_ - centroid_x * centroid_y;
        if ( !std::isfinite( variance_x ) || !std::isfinite( variance_y ) ||
             !std::isfinite( covariance ) || variance_x <= 0.0 ||
             variance_y <= 0.0 )
            throw std::runtime_error(
                "Tabulated source profile must have positive lateral width" );

        const double mean_variance = 0.5 * ( variance_x + variance_y );
        const double variance_radius = std::sqrt(
            0.25 * ( variance_x - variance_y ) *
                ( variance_x - variance_y ) +
            covariance * covariance );
        const double major_variance = mean_variance + variance_radius;
        const double minor_variance = mean_variance - variance_radius;
        if ( minor_variance <= 0.0 )
            throw std::runtime_error(
                "Tabulated source profile must have positive principal "
                "widths" );
        d4_sigma_major_ = 4.0 * std::sqrt( major_variance );
        d4_sigma_minor_ = 4.0 * std::sqrt( minor_variance );
        lateral_d4_sigma_ =
            std::sqrt( d4_sigma_major_ * d4_sigma_minor_ );
        profile_azimuth_ =
            variance_radius <=
                    64.0 * std::numeric_limits<double>::epsilon() *
                        mean_variance
                ? 0.0
                : 0.5 *
                      std::atan2( 2.0 * covariance,
                                  variance_x - variance_y );
    }

    void updateAxialState()
    {
        const double source_aspect_ratio =
            std::max( 2.0 * current_depth_ / lateral_d4_sigma_, 0.001 );
        const double exponent_power = std::max(
            0.0, std::min( exponent_slope_ *
                                   std::log2( source_aspect_ratio ) +
                               exponent_intercept_,
                           9.0 ) );
        axial_exponent_ = std::pow( 2.0, exponent_power );

        fast_exponent_ = 0;
        for ( const int candidate : { 1, 2, 4, 8 } )
            if ( std::abs( axial_exponent_ - candidate ) <=
                 16.0 * std::numeric_limits<double>::epsilon() * candidate )
                fast_exponent_ = candidate;

        const double inverse_exponent = 1.0 / axial_exponent_;
        const double axial_integral =
            current_depth_ * std::tgamma( inverse_exponent ) /
            ( axial_exponent_ * std::pow( 3.0, inverse_exponent ) );
        inverse_source_volume_ = 1.0 / ( planar_integral_ * axial_integral );
        updateAbsorptionState();

        constexpr double cutoff_tolerance = 1.0e-3;
        const double cutoff_ratio =
            std::pow( -std::log( cutoff_tolerance ) / 3.0, inverse_exponent );
        cutoff_depth_ = current_depth_ * cutoff_ratio;
    }

    void updateAbsorptionState()
    {
        current_absorptivity_ = std::visit(
            [&]( const auto& model )
            { return model.evaluate( absorption_aspect_ratio_ ); },
            absorption_ );
        normalization_per_watt_ =
            current_absorptivity_ * inverse_source_volume_;
    }

    table_type table_;
    AbsorptionVariant absorption_;
    double x0_ = 0.0;
    double x1_ = 0.0;
    double y0_ = 0.0;
    double y1_ = 0.0;
    double dx_ = 0.0;
    double dy_ = 0.0;
    double inverse_dx_ = 0.0;
    double inverse_dy_ = 0.0;
    double planar_integral_ = 0.0;
    double d4_sigma_major_ = 0.0;
    double d4_sigma_minor_ = 0.0;
    double lateral_d4_sigma_ = 0.0;
    double profile_azimuth_ = 0.0;
    double search_radius_ = 0.0;
    double minimum_depth_ = 0.0;
    double measured_depth_ = 0.0;
    double current_depth_ = 0.0;
    double exponent_slope_ = 0.0;
    double exponent_intercept_ = 0.0;
    double axial_exponent_ = 0.0;
    double inverse_source_volume_ = 0.0;
    double normalization_per_watt_ = 0.0;
    double cutoff_depth_ = 0.0;
    double feedback_temperature_ = 0.0;
    double absorption_aspect_ratio_ = 0.0;
    double current_absorptivity_ = 0.0;
    int nx_ = 0;
    int ny_ = 0;
    int fast_exponent_ = 0;
    bool scan_path_frame_ = false;
    bool depth_feedback_enabled_ = false;
    bool dynamic_absorption_ = false;
};

template <typename MemorySpace>
using HeatSourceVariant =
    std::variant<GaussianSource, TabulatedSource<MemorySpace>>;

template <typename MemorySpace, typename ExecutionSpace>
HeatSourceVariant<MemorySpace>
createHeatSource( const Inputs& inputs, MPI_Comm comm,
                  const ExecutionSpace& execution_space )
{
    if ( inputs.source.type == "tabulated" )
        return HeatSourceVariant<MemorySpace>(
            std::in_place_type<TabulatedSource<MemorySpace>>, inputs, comm,
            execution_space );
    return HeatSourceVariant<MemorySpace>( std::in_place_type<GaussianSource>,
                                           inputs );
}

} // namespace Finch

#endif
