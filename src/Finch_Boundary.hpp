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

#ifndef Boundary_H
#define Boundary_H

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Finch_BoundaryConditions.hpp>

namespace Finch
{

class Boundary
{
  public:
    Boundary( const BoundaryConditions& conditions, const double cell_size )
        : cell_size( cell_size )
    {
        for ( int face = 0; face < 6; ++face )
        {
            const auto& condition = conditions.faces[face];
            boundary_types[face] = condition.type;
            boundary_values[face] = condition.value;
            convection_coefficients[face] = condition.convection_coefficient;
            emissivities[face] = condition.emissivity;
            ambient_temperatures[face] = condition.ambient_temperature;
        }
        storeTypes();
    }

    // Use integer representation of BC types for device access.
    void storeTypes()
    {
        if ( !std::isfinite( cell_size ) || cell_size <= 0.0 )
            throw std::runtime_error(
                "Boundary cell size must be finite and positive." );

        has_mixed_boundary = false;
        for ( int d = 0; d < 6; d++ )
        {
            if ( !std::isfinite( boundary_values[d] ) )
                throw std::runtime_error( "Boundary values must be finite." );
            if ( boundary_types[d] == "dirichlet" )
                boundary_int[d] = 0;
            else if ( boundary_types[d] == "neumann" )
                boundary_int[d] = 1;
            else if ( boundary_types[d] == "adiabatic" )
                boundary_int[d] = 2;
            else if ( boundary_types[d] == "convection_radiation" )
            {
                if ( !std::isfinite( convection_coefficients[d] ) ||
                     convection_coefficients[d] < 0.0 ||
                     !std::isfinite( emissivities[d] ) ||
                     emissivities[d] < 0.0 || emissivities[d] > 1.0 ||
                     !std::isfinite( ambient_temperatures[d] ) ||
                     ambient_temperatures[d] <= 0.0 )
                    throw std::runtime_error(
                        "Invalid convection/radiation boundary parameters." );
                boundary_int[d] = 3;
                has_mixed_boundary = true;
            }
            else
                throw std::runtime_error( "Invalid boundary type: " +
                                          boundary_types[d] );
        }
    }

    // Create and store an array of each set of BC indices.
    template <typename LocalGridType, typename EntityType>
    void create( LocalGridType local_grid, EntityType )
    {
        // Generate the boundary condition index spaces.
        int count = 0;
        for ( int d = 0; d < 3; d++ )
        {
            for ( int dir = -1; dir < 2; dir += 2 )
            {
                boundary_planes[count] = { 0, 0, 0 };
                boundary_planes[count][d] = dir;

                // Get the boundary indices for this plane (each one is a
                // separate, contiguous index space).
                boundary_spaces[count] = local_grid->boundaryIndexSpace(
                    Cabana::Grid::Ghost(), EntityType(),
                    boundary_planes[count][0], boundary_planes[count][1],
                    boundary_planes[count][2] );
                count++;
            }
        }
    }

    template <typename ExecSpace, typename TempViewType,
              typename ConductivityModel>
    void update( ExecSpace exec_space, TempViewType T, TempViewType previous_T,
                 const ConductivityModel& conductivity_model )
    {
        if ( !has_mixed_boundary )
        {
            updateSimple( exec_space, T );
            return;
        }

        // Update the boundary on each face of the cube. Fuse iteration over all
        // index spaces to avoid launching 6 separate kernels.
        auto planes = boundary_planes;
        auto type = boundary_int;
        auto values = boundary_values;
        auto convection = convection_coefficients;
        auto emissivity = emissivities;
        auto ambient = ambient_temperatures;
        auto conductivity = conductivity_model;
        const double spacing = cell_size;
        constexpr double sigma = 5.67e-8;
        Cabana::Grid::grid_parallel_for(
            "Finch::boundary_update", exec_space, boundary_spaces,
            KOKKOS_LAMBDA( const int b, const int i, const int j,
                           const int k ) {
                const int inside_i = i - planes[b][0];
                const int inside_j = j - planes[b][1];
                const int inside_k = k - planes[b][2];
                const double inside_temperature =
                    T( inside_i, inside_j, inside_k, 0 );

                if ( type[b] == 0 )
                    T( i, j, k, 0 ) = values[b];
                else if ( type[b] == 1 )
                    T( i, j, k, 0 ) = inside_temperature + values[b] * spacing;
                else if ( type[b] == 2 )
                    T( i, j, k, 0 ) = inside_temperature;
                else
                {
                    // Linearize the radiative flux about the completed
                    // boundary temperature to avoid nonlinear iterations.
                    const double surface_temperature = previous_T( i, j, k, 0 );
                    const double ambient_temperature = ambient[b];
                    const double surface_squared =
                        surface_temperature * surface_temperature;
                    const double ambient_squared =
                        ambient_temperature * ambient_temperature;
                    const double radiative_coefficient =
                        sigma * emissivity[b] *
                        ( surface_squared + ambient_squared ) *
                        ( surface_temperature + ambient_temperature );
                    const double effective_coefficient = Kokkos::fmax(
                        convection[b] + radiative_coefficient, 1.0e-15 );
                    const double conductive_coefficient =
                        conductivity.conductivity( surface_temperature ) /
                        spacing;
                    T( i, j, k, 0 ) =
                        ( conductive_coefficient * inside_temperature +
                          effective_coefficient * ambient_temperature ) /
                        ( conductive_coefficient + effective_coefficient );
                }
            } );
    }

  protected:
    template <typename ExecSpace, typename TempViewType>
    void updateSimple( ExecSpace exec_space, TempViewType T )
    {
        auto planes = boundary_planes;
        auto type = boundary_int;
        auto values = boundary_values;
        const double spacing = cell_size;
        Cabana::Grid::grid_parallel_for(
            "Finch::boundary_update", exec_space, boundary_spaces,
            KOKKOS_LAMBDA( const int b, const int i, const int j,
                           const int k ) {
                if ( type[b] == 0 )
                    T( i, j, k, 0 ) = values[b];
                else if ( type[b] == 1 )
                    T( i, j, k, 0 ) = T( i - planes[b][0], j - planes[b][1],
                                         k - planes[b][2], 0 ) +
                                      values[b] * spacing;
                else
                    T( i, j, k, 0 ) = T( i - planes[b][0], j - planes[b][1],
                                         k - planes[b][2], 0 );
            } );
    }

    //! Boundary types for each plane.
    std::array<std::string, 6> boundary_types;
    //! Boundary values for each plane.
    Kokkos::Array<double, 6> boundary_values;
    Kokkos::Array<double, 6> convection_coefficients;
    Kokkos::Array<double, 6> emissivities;
    Kokkos::Array<double, 6> ambient_temperatures;
    //! Uniform grid spacing used to convert Neumann gradients to temperature.
    double cell_size;
    //! Boundary types for each plane, converted to int for device.
    Kokkos::Array<int, 6> boundary_int;
    bool has_mixed_boundary = false;
    //! Boundary indices for each plane.
    Kokkos::Array<Cabana::Grid::IndexSpace<3>, 6> boundary_spaces;
    // Boundary details.
    Kokkos::Array<Kokkos::Array<int, 3>, 6> boundary_planes;
};
} // namespace Finch

#endif
