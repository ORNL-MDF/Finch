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

namespace Finch
{

class Boundary
{
  public:
    // Constructor with BC types and values. Neumann values are outward
    // temperature gradients in K/m and are converted using cell_size.
    Boundary( std::array<std::string, 6> types, Kokkos::Array<double, 6> values,
              const double cell_size )
        : boundary_types( types )
        , boundary_values( values )
        , cell_size( cell_size )
    {
        storeTypes();
    }

    // Constructor where no values are needed.
    Boundary( std::array<std::string, 6> types, const double cell_size )
        : boundary_types( types )
        , boundary_values( { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } )
        , cell_size( cell_size )
    {
        for ( int d = 0; d < 6; d++ )
            if ( boundary_types[d] == "dirichlet" ||
                 boundary_types[d] == "neumann" )
                throw std::runtime_error(
                    "Boundary condition requires a separate value input." );

        storeTypes();
    }

    // Use integer representation of BC types for device access.
    void storeTypes()
    {
        if ( !std::isfinite( cell_size ) || cell_size <= 0.0 )
            throw std::runtime_error(
                "Boundary cell size must be finite and positive." );

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

    // Apply the boundary condition.
    template <typename ExecSpace, typename TempViewType>
    void update( ExecSpace exec_space, TempViewType T )
    {
        // Update the boundary on each face of the cube. Fuse iteration over all
        // index spaces to avoid launching 6 separate kernels.
        auto planes = boundary_planes;
        auto type = boundary_int;
        auto values = boundary_values;
        const double spacing = cell_size;
        Cabana::Grid::grid_parallel_for(
            "Finch::boundary_update", exec_space, boundary_spaces,
            KOKKOS_LAMBDA( const int b, const int i, const int j,
                           const int k ) {
                // Corresponds to the user-passed strings set up in the
                // constructor.
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

  protected:
    //! Boundary types for each plane.
    std::array<std::string, 6> boundary_types;
    //! Boundary values for each plane.
    Kokkos::Array<double, 6> boundary_values;
    //! Uniform grid spacing used to convert Neumann gradients to temperature.
    double cell_size;
    //! Boundary types for each plane, converted to int for device.
    Kokkos::Array<int, 6> boundary_int;
    //! Boundary indices for each plane.
    Kokkos::Array<Cabana::Grid::IndexSpace<3>, 6> boundary_spaces;
    // Boundary details.
    Kokkos::Array<Kokkos::Array<int, 3>, 6> boundary_planes;
};
} // namespace Finch

#endif
