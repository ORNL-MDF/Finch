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
#include <cstdlib>
#include <math.h>
#include <mpi.h>
#include <string>

#include <Kokkos_Core.hpp>

#include "Finch_Core.hpp"

void run( int argc, char* argv[] )
{
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;
    exec_space execution_space;

    // initialize the simulation
    Finch::Inputs db( MPI_COMM_WORLD, argc, argv );

    if ( db.comm_rank == 0 )
        std::cout << "Kokkos execution space: " << exec_space::name()
                  << std::endl;

    if constexpr ( !Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                               memory_space>::accessible )
    {
        // Cabana 0.8 passes device-resident halo buffers directly to MPI.
        // There is no portable MPI capability query, so distributed device
        // runs require an explicit deployment assertion.
        const char* gpu_aware = std::getenv( "FINCH_GPU_AWARE_MPI" );
        int local_gpu_aware =
            gpu_aware != nullptr && std::string( gpu_aware ) == "1";
        int all_gpu_aware = 0;
        MPI_Allreduce( &local_gpu_aware, &all_gpu_aware, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD );
        if ( db.comm_size > 1 && !all_gpu_aware )
            throw std::runtime_error(
                "Distributed accelerator runs require GPU-aware MPI. Set "
                "FINCH_GPU_AWARE_MPI=1 after verifying the MPI deployment." );
    }

    // initialize a moving beam
    Finch::MovingBeam beam( db.source.scan_path_file, MPI_COMM_WORLD );

    // Define boundary condition details.
    std::array<std::string, 6> bc_types = { "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic",
                                            "adiabatic", "adiabatic" };

    // create the global mesh
    Finch::Grid<memory_space> grid(
        MPI_COMM_WORLD, db.space.cell_size, db.space.global_low_corner,
        db.space.global_high_corner, db.space.ranks_per_dim, bc_types,
        db.space.initial_temperature, execution_space );

    // Create the solver
    auto fd = Finch::createSolver( db, grid );

    // Run the full single layer problem
    Finch::Layer app( db, grid );
    app.run( execution_space, db, grid, beam, fd );

    // Write the temperature data used by ExaCA/other post-processing
    app.writeSolidificationData( db.sampling, grid.getComm() );
    app.getLowerSolidificationDataBounds( grid.getComm() );
    app.getUpperSolidificationDataBounds( grid.getComm() );
}

int main( int argc, char* argv[] )
{
    int provided = MPI_THREAD_SINGLE;
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided );
    if ( provided < MPI_THREAD_FUNNELED )
    {
        std::cerr << "Finch requires MPI_THREAD_FUNNELED or better."
                  << std::endl;
        MPI_Finalize();
        return 1;
    }

    Kokkos::initialize( argc, argv );

    int local_status = 0;
    try
    {
        run( argc, argv );
    }
    catch ( const std::exception& e )
    {
        int rank = 0;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        std::cerr << "Finch rank " << rank << ": " << e.what() << std::endl;
        local_status = 1;
    }

    int global_status = 0;
    MPI_Allreduce( &local_status, &global_status, 1, MPI_INT, MPI_MAX,
                   MPI_COMM_WORLD );

    Kokkos::finalize();
    MPI_Finalize();

    return global_status;
}
