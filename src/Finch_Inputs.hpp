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
  \file Simulation.hpp
  \brief Simulation inputs
*/

#ifndef Inputs_H
#define Inputs_H

#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unistd.h>

#include <nlohmann/json.hpp>

namespace Finch
{

// Info macro for writing on master
#define Info                                                                   \
    if ( comm_rank == 0 )                                                      \
    std::cout

struct Output
{
    int total_steps;
    int interval;

    void setInterval( const int num_steps )
    {
        // If total_output_steps = 0, set increment to greater than the number
        // of time steps to avoid printing output, otherwise bound
        // output_interval to be greater than 1 and no larger than the total
        // number of time steps
        if ( total_steps == 0 )
            interval = num_steps + 1;
        else
        {
            interval = static_cast<int>( ( num_steps / total_steps ) );
            interval = std::max( std::min( interval, num_steps ), 1 );
        }
    }
};

struct Time
{
    double Co;
    double start_time;
    double end_time;
    double time_step;
    double time;
    int num_steps;
    Output output;
    Output monitor;
};

struct Space
{
    double initial_temperature;
    double cell_size;
    std::array<double, 3> global_low_corner;
    std::array<double, 3> global_high_corner;
    std::array<int, 3> ranks_per_dim;
};

struct Source
{
    double absorption;
    std::array<double, 3> two_sigma;
    std::array<double, 3> r;
    std::string scan_path_file;
};

struct Properties
{
    double density;
    double specific_heat;
    double thermal_conductivity;
    double thermal_diffusivity;
    double latent_heat;
    double solidus;
    double liquidus;
};

struct Sampling
{
    std::string type;
    std::string format;
    std::string directory_name = "solidification";
    bool enabled;
};

struct TimeMonitor
{
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::duration<double> elapsed_seconds;
    double total_elapsed_time;
    int total_monitor_steps;
    int num_steps;
    int comm_rank;

    // Default constructor
    TimeMonitor(){};

    // Constructor with MPI_Comm
    TimeMonitor( MPI_Comm comm, Time& time )
        : total_elapsed_time( 0.0 )
    {
        start_time = std::chrono::high_resolution_clock::now();
        MPI_Comm_rank( comm, &comm_rank );
        total_monitor_steps = time.monitor.total_steps;
        num_steps = time.num_steps;
    }

    void update()
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end_time - start_time;
        total_elapsed_time += elapsed_seconds.count();

        start_time = std::chrono::high_resolution_clock::now();
    }

    void write( int step )
    {
        update();

        Info << "Time Step: " << step << "/" << num_steps << ", "
             << "Elapsed: " << std::fixed << std::setprecision( 6 )
             << elapsed_seconds.count() << " seconds, "
             << "Total: " << std::fixed << std::setprecision( 6 )
             << total_elapsed_time << " seconds" << std::endl;
    }
};

class Inputs
{
  public:
    Time time;
    Space space;
    Source source;
    Properties properties;
    Sampling sampling;
    TimeMonitor time_monitor;

    int comm_rank;
    int comm_size;

    // constructor for Finch run as standalone code
    Inputs( MPI_Comm comm, int argc, char* argv[] )
    {
        MPI_Comm_rank( comm, &comm_rank );
        MPI_Comm_size( comm, &comm_size );
        std::string filename = getFilename( argc, argv );
        parseInputFile( comm, filename );
    }
    // constructor for coupled run of Finch with ExaCA
    Inputs( MPI_Comm comm, std::string filename )
    {
        MPI_Comm_rank( comm, &comm_rank );
        MPI_Comm_size( comm, &comm_size );
        parseInputFile( comm, filename );
    }

    void write()
    {
        Info << "Simulation will be performed using parameters: " << std::endl;

        // Print time
        Info << "Time:" << std::endl;
        Info << "  Co: " << time.Co << std::endl;
        Info << "  Start Time: " << time.start_time << std::endl;
        Info << "  End Time: " << time.end_time << std::endl;
        Info << "  Num Output Steps: " << time.output.total_steps << std::endl;
        Info << "  Num Monitor Steps: " << time.monitor.total_steps
             << std::endl;

        // Print space
        Info << "Space:" << std::endl;
        Info << "  Initial temperature: " << space.initial_temperature
             << std::endl;
        Info << "  Cell Size: " << space.cell_size << std::endl;
        Info << "  Global Low Corner:" << std::endl;
        Info << "    X: " << space.global_low_corner[0] << std::endl;
        Info << "    Y: " << space.global_low_corner[1] << std::endl;
        Info << "    Z: " << space.global_low_corner[2] << std::endl;
        Info << "  Global High Corner:" << std::endl;
        Info << "    X: " << space.global_high_corner[0] << std::endl;
        Info << "    Y: " << space.global_high_corner[1] << std::endl;
        Info << "    Z: " << space.global_high_corner[2] << std::endl;

        // Print properties
        Info << "Properties:" << std::endl;
        Info << "  Density: " << properties.density << std::endl;
        Info << "  Specific Heat: " << properties.specific_heat << std::endl;
        Info << "  Thermal Conductivity: " << properties.thermal_conductivity
             << std::endl;
        Info << "  Latent Heat: " << properties.latent_heat << std::endl;
        Info << "  Solidus: " << properties.solidus << std::endl;
        Info << "  Liquidus: " << properties.liquidus << std::endl;

        // Print source
        Info << "Source:" << std::endl;
        Info << "  Absorption: " << source.absorption << std::endl;
        Info << "  two-sigma:" << std::endl;
        Info << "    X: " << source.two_sigma[0] << std::endl;
        Info << "    Y: " << source.two_sigma[1] << std::endl;
        Info << "    Z: " << source.two_sigma[2] << std::endl;
        Info << "  scan path file: " << source.scan_path_file << std::endl;

        // Print solidification output options
        Info << "Sampling:" << std::endl;
        if ( sampling.enabled )
        {
            Info << "  type: " << sampling.type << std::endl;
            Info << "  format:" << sampling.format << std::endl;
            Info << "  directory name:" << sampling.directory_name << std::endl;
        }
        else
        {
            Info << "Skipping optional sampling." << std::endl;
        }
    }

  private:
    std::string getFilename( int argc, char* argv[] )
    {
        const char* filename = nullptr;
        int option;

        while ( ( option = getopt( argc, argv, "i:" ) ) != -1 )
        {
            if ( option == 'i' )
            {
                filename = optarg;
            }
            else
            {
                std::string error_message =
                    "Error: the input file must be specified using -i "
                    "<input_json_file>";
                throw std::runtime_error( error_message );
            }
        }
        std::string filename_s( filename );
        return filename_s;
    }

    void parseInputFile( MPI_Comm comm, const std::string filename )
    {
        readInput( filename );

        write();

        // create auxilary properties
        properties.thermal_diffusivity =
            ( properties.thermal_conductivity ) /
            ( properties.density * properties.specific_heat );

        time.time_step = ( time.Co * space.cell_size * space.cell_size ) /
                         ( properties.thermal_diffusivity );

        Info << "Calculated time step: " << time.time_step << std::endl;

        time.time = time.start_time;

        time.num_steps = static_cast<int>( ( time.end_time - time.start_time ) /
                                           ( time.time_step ) );

        time.output.setInterval( time.num_steps );
        time.monitor.setInterval( time.num_steps );

        // initialize time monitoring
        time_monitor = TimeMonitor( comm, time );
    }

    void readInput( const std::string filename )
    {
        // parse input file
        std::ifstream db_stream( filename );
        nlohmann::json db = nlohmann::json::parse( db_stream );

        // Read time components
        time.Co = db["time"]["Co"];
        time.start_time = db["time"]["start_time"];
        time.end_time = db["time"]["end_time"];
        time.output.total_steps = db["time"]["total_output_steps"];
        time.monitor.total_steps = db["time"]["total_monitor_steps"];

        // Read space components
        space.initial_temperature = db["space"]["initial_temperature"];
        space.cell_size = db["space"]["cell_size"];
        space.global_low_corner = db["space"]["global_low_corner"];
        space.global_high_corner = db["space"]["global_high_corner"];

        /*
          Default block partitioner. This relies on MPI_Cart_create to
          balance the number of ranks in each direction. This partitioning
          is best only in the global mesh is a uniform cube.
        */
        std::array<int, 3> default_ranks_per_dim = { 0, 0, 0 };

        std::array<int, 3> ranks_per_dim = default_ranks_per_dim;
        if ( db["space"].contains( "ranks_per_dim" ) )
            ranks_per_dim = db["space"]["ranks_per_dim"];

        // Invalid partition strategy selected. Use Default block partioner.
        int product = ranks_per_dim[0] * ranks_per_dim[1] * ranks_per_dim[2];

        if ( product != comm_size )
        {
            ranks_per_dim = default_ranks_per_dim;
        }

        space.ranks_per_dim = ranks_per_dim;

        // Read properties components
        properties.density = db["properties"]["density"];
        properties.specific_heat = db["properties"]["specific_heat"];
        properties.thermal_conductivity =
            db["properties"]["thermal_conductivity"];
        properties.latent_heat = db["properties"]["latent_heat"];
        properties.solidus = db["properties"]["solidus"];
        properties.liquidus = db["properties"]["liquidus"];

        // Read heat source components
        source.absorption = db["source"]["absorption"];
        source.two_sigma = db["source"]["two_sigma"];

        source.two_sigma[0] = fabs( source.two_sigma[0] );
        source.two_sigma[1] = fabs( source.two_sigma[1] );
        source.two_sigma[2] = fabs( source.two_sigma[2] );

        source.scan_path_file = db["source"]["scan_path_file"];

        // Read sampling components
        sampling.enabled = false;
        if ( db.contains( "sampling" ) )
        {
            const std::string sampling_type = db["sampling"]["type"];

            if ( sampling_type == "solidification_data" )
            {
                sampling.type = sampling_type;
                sampling.enabled = true;
            }

            const std::string sampling_format = db["sampling"]["format"];

            if ( sampling_format == "exaca" )
            {
                sampling.format = sampling_format;
            }
            else
            {
                sampling.format = "default";
            }

            if ( db["sampling"].contains( "directory_name" ) )
            {
                sampling.directory_name = db["sampling"]["directory_name"];
            }
        }
    }
};

} // namespace Finch

#endif
