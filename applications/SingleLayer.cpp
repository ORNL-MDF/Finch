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

#include <array>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <string>

#include <Kokkos_Core.hpp>

#include "Finch_Core.hpp"

void run( int argc, char* argv[] )
{
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;

    // initialize the simulation
    Finch::Inputs db( MPI_COMM_WORLD, argc, argv );

    // initialize a moving beam
    Finch::MovingBeam beam( db.source.scan_path_file );

    // Define boundary condition details.
    std::array<std::string, 6> bc_types = { "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic" };

    // create the global mesh
    Finch::Grid<memory_space> grid(
        MPI_COMM_WORLD, db.space.cell_size, db.space.global_low_corner,
        db.space.global_high_corner, db.space.ranks_per_dim, bc_types,
        db.space.initial_temperature );

    // Create the solver
    auto fd = Finch::createSolver( db, grid );

    // Run the full single layer problem
    Finch::Layer app( db, grid );
    app.run( exec_space(), db, grid, beam, fd );

    // Write the temperature data used by ExaCA/other post-processing
    app.writeSolidificationData( grid.getComm() );
    app.getLowerSolidificationDataBounds( grid.getComm() );
    app.getUpperSolidificationDataBounds( grid.getComm() );
}

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    Kokkos::initialize( argc, argv );

    run( argc, argv );

    Kokkos::finalize();
    MPI_Finalize();

    return 0;
}
