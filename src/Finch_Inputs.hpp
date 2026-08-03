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

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include <mpi.h>
#include <nlohmann/json.hpp>

#include <Finch_Version.hpp>

namespace Finch
{

// Info macro for writing on master
#define Info                                                                   \
    if ( comm_rank == 0 )                                                      \
    std::cout

struct Output
{
    int total_steps = 0;
    int interval = 1;

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

    bool isDue( const int completed_steps, const int num_steps ) const
    {
        if ( total_steps <= 0 || num_steps <= 0 || completed_steps <= 0 )
            return false;

        // Integer arithmetic gives evenly distributed output and always
        // includes the final step without accumulating floating-point error.
        const long long capped_outputs = std::min( total_steps, num_steps );
        return ( static_cast<long long>( completed_steps ) * capped_outputs ) /
                   num_steps >
               ( static_cast<long long>( completed_steps - 1 ) *
                 capped_outputs ) /
                   num_steps;
    }
};

struct Time
{
    double Co = 0.0;
    double start_time = 0.0;
    double end_time = 0.0;
    double time_step = 0.0;
    double time = 0.0;
    int num_steps = 0;
    Output output;
    Output monitor;
};

struct Space
{
    double initial_temperature = 0.0;
    double cell_size = 0.0;
    std::array<double, 3> global_low_corner = { 0.0, 0.0, 0.0 };
    std::array<double, 3> global_high_corner = { 0.0, 0.0, 0.0 };
    std::array<int, 3> ranks_per_dim = { 0, 0, 0 };
};

struct Source
{
    double absorption = 0.0;
    std::array<double, 3> two_sigma = { 0.0, 0.0, 0.0 };
    std::string scan_path_file;
};

struct Properties
{
    double density = 0.0;
    double specific_heat = 0.0;
    double thermal_conductivity = 0.0;
    double thermal_diffusivity = 0.0;
    double latent_heat = 0.0;
    double solidus = 0.0;
    double liquidus = 0.0;
};

struct Sampling
{
    std::string type;
    std::string format = "default";
    std::string directory_name = "solidification";
    bool enabled = false;
};

struct TimeMonitor
{
    std::chrono::steady_clock::time_point start_time;
    std::chrono::duration<double> elapsed_seconds;
    double total_elapsed_time = 0.0;
    int total_monitor_steps = 0;
    int num_steps = 0;
    int comm_rank = 0;

    // Default constructor
    TimeMonitor() = default;

    // Constructor with MPI_Comm
    TimeMonitor( MPI_Comm comm, Time& time )
    {
        start_time = std::chrono::steady_clock::now();
        MPI_Comm_rank( comm, &comm_rank );
        total_monitor_steps = time.monitor.total_steps;
        num_steps = time.num_steps;
    }

    void update()
    {
        auto end_time = std::chrono::steady_clock::now();
        elapsed_seconds = end_time - start_time;
        total_elapsed_time += elapsed_seconds.count();

        start_time = std::chrono::steady_clock::now();
    }

    void reset()
    {
        total_elapsed_time = 0.0;
        start_time = std::chrono::steady_clock::now();
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
        calcAuxiliaryProperties( comm );
    }
    // constructor for coupled run of Finch with ExaCA - inputs potentially
    // spread across both files
    Inputs( MPI_Comm comm, const std::string filename,
            const int input_file_number = 0 )
    {
        MPI_Comm_rank( comm, &comm_rank );
        MPI_Comm_size( comm, &comm_size );
        parseInputFile( comm, filename, input_file_number );
        calcAuxiliaryProperties( comm );
    }

    void write()
    {
        Info << "Finch version: " << version() << " (" << commitHash() << ")"
             << std::endl;
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
        if ( filename == nullptr )
            throw std::runtime_error(
                "Error: the input file must be specified using -i "
                "<input_json_file>" );

        return std::string( filename );
    }

    nlohmann::json readInputDocument( MPI_Comm comm,
                                      const std::string& filename )
    {
        std::string contents;
        std::string error;
        if ( comm_rank == 0 )
        {
            std::ifstream input_data_stream( filename );
            if ( !input_data_stream )
                error = "Error: cannot open Finch input file " + filename;
            else
            {
                std::ostringstream buffer;
                buffer << input_data_stream.rdbuf();
                contents = buffer.str();
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

        int contents_size = static_cast<int>( contents.size() );
        MPI_Bcast( &contents_size, 1, MPI_INT, 0, comm );
        contents.resize( contents_size );
        if ( contents_size > 0 )
            MPI_Bcast( contents.data(), contents_size, MPI_CHAR, 0, comm );

        try
        {
            return nlohmann::json::parse( contents );
        }
        catch ( const nlohmann::json::exception& e )
        {
            throw std::runtime_error( "Error parsing " + filename + ": " +
                                      e.what() );
        }
    }

    bool requiredSectionsFound( const std::vector<bool>& found ) const
    {
        // Sampling is optional; the other four sections are required.
        return std::all_of( found.begin(), found.begin() + 4,
                            []( const bool value ) { return value; } );
    }

    void parseInputFile( MPI_Comm comm, const std::string& filename,
                         const int input_file_number = 0 )
    {
        // Input file is either a Finch input file or an ExaCA input file with a
        // Finch object
        Info << "Parsing input file " << input_file_number << std::endl;
        nlohmann::json input_data_raw = readInputDocument( comm, filename );
        if ( !input_data_raw.contains( "Finch" ) )
        {
            // This is a Finch input file and should have all 5 sections
            std::vector<bool> found_sections = readSections( input_data_raw );
            if ( !requiredSectionsFound( found_sections ) )
                throw std::runtime_error(
                    "Error: Missing top-level sections of Finch input file, "
                    "see README for proper input file format" );
        }
        else
        {
            // This is an ExaCA input file - get data from the Finch object
            nlohmann::json top_level_input_data = input_data_raw["Finch"];
            // Sections to parse
            std::vector<std::string> input_file_sections = {
                "time", "space", "properties", "source", "sampling" };
            const int num_inp_file_sections = input_file_sections.size();
            if ( !top_level_input_data.contains( "layers" ) )
            {
                // All inputs should be present at the top level of the file
                std::vector<bool> found_sections =
                    readSections( top_level_input_data );
                if ( !requiredSectionsFound( found_sections ) )
                    throw std::runtime_error(
                        "Error: Missing top-level sections of Finch input "
                        "file, see README for proper input file format" );
            }
            else
            {
                // Check top level for inputs and parse
                std::vector<bool> found_sections_top_level =
                    readSections( top_level_input_data );
                // Check individual layer level for inputs and parse
                std::vector<bool> found_sections_layer_level = readSections(
                    top_level_input_data["layers"][input_file_number] );
                // Ensure each input was present at least once, warn about
                // redundant inputs (layer level takes priority)
                for ( int n = 0; n < num_inp_file_sections; n++ )
                {
                    if ( ( found_sections_top_level[n] ) &&
                         ( found_sections_layer_level[n] ) )
                    {
                        Info << "Warning: Finch input object "
                             << input_file_sections[n]
                             << " has multiple values given; values from "
                                "`layers` object will be used"
                             << std::endl;
                    }
                    else if ( n < 4 && ( !found_sections_top_level[n] ) &&
                              ( !found_sections_layer_level[n] ) )
                    {
                        std::string err_message = "Error: Finch input object " +
                                                  input_file_sections[n] +
                                                  " was not found";
                        throw std::runtime_error( err_message );
                    }
                }
            }
        }
        validateInputs();
        selectRankDecomposition();
        write();
    }

    void validateInputs() const
    {
        const auto finite = []( const double value )
        { return std::isfinite( value ); };

        if ( !finite( time.Co ) || time.Co <= 0.0 || time.Co > 1.0 / 6.0 )
            throw std::runtime_error(
                "Error: time.Co must be in (0, 1/6] for the 3D explicit "
                "diffusion stencil" );
        if ( !finite( time.start_time ) || !finite( time.end_time ) ||
             time.end_time <= time.start_time )
            throw std::runtime_error(
                "Error: end_time must be finite and greater than start_time" );
        if ( time.output.total_steps < 0 || time.monitor.total_steps < 0 )
            throw std::runtime_error(
                "Error: output and monitor step counts cannot be negative" );

        if ( !finite( space.initial_temperature ) ||
             !finite( space.cell_size ) || space.cell_size <= 0.0 )
            throw std::runtime_error(
                "Error: initial_temperature must be finite and cell_size "
                "must be positive" );
        for ( int d = 0; d < 3; ++d )
        {
            const double extent =
                space.global_high_corner[d] - space.global_low_corner[d];
            if ( !finite( space.global_low_corner[d] ) ||
                 !finite( space.global_high_corner[d] ) || extent <= 0.0 )
                throw std::runtime_error(
                    "Error: every global domain extent must be positive" );

            const double cells = extent / space.cell_size;
            if ( std::abs( cells - std::round( cells ) ) >
                 100.0 * std::numeric_limits<double>::epsilon() *
                     std::max( 1.0, std::abs( cells ) ) )
                throw std::runtime_error(
                    "Error: every global domain extent must be evenly "
                    "divisible by cell_size" );
            if ( space.ranks_per_dim[d] < 0 )
                throw std::runtime_error(
                    "Error: ranks_per_dim values cannot be negative" );
        }

        if ( !finite( properties.density ) || properties.density <= 0.0 ||
             !finite( properties.specific_heat ) ||
             properties.specific_heat <= 0.0 ||
             !finite( properties.thermal_conductivity ) ||
             properties.thermal_conductivity <= 0.0 ||
             !finite( properties.latent_heat ) ||
             properties.latent_heat < 0.0 || !finite( properties.solidus ) ||
             !finite( properties.liquidus ) ||
             properties.liquidus <= properties.solidus )
            throw std::runtime_error(
                "Error: material properties must be finite, density, heat "
                "capacity, and conductivity must be positive, latent heat "
                "must be nonnegative, and liquidus must exceed solidus" );

        if ( !finite( source.absorption ) || source.absorption < 0.0 ||
             source.absorption > 1.0 || source.scan_path_file.empty() )
            throw std::runtime_error(
                "Error: source absorption must be in [0,1] and "
                "scan_path_file cannot be empty" );
        for ( const double sigma : source.two_sigma )
            if ( !finite( sigma ) || sigma <= 0.0 )
                throw std::runtime_error(
                    "Error: every source two_sigma value must be positive" );

        if ( sampling.enabled && sampling.directory_name.empty() )
            throw std::runtime_error(
                "Error: sampling directory_name cannot be empty" );
    }

    void selectRankDecomposition()
    {
        const int requested_product = space.ranks_per_dim[0] *
                                      space.ranks_per_dim[1] *
                                      space.ranks_per_dim[2];
        std::array<int, 3> cells;
        for ( int d = 0; d < 3; ++d )
            cells[d] = static_cast<int>( std::llround(
                ( space.global_high_corner[d] - space.global_low_corner[d] ) /
                space.cell_size ) );

        if ( requested_product == comm_size )
        {
            for ( int d = 0; d < 3; ++d )
                if ( space.ranks_per_dim[d] <= 0 ||
                     space.ranks_per_dim[d] > cells[d] )
                    throw std::runtime_error(
                        "Error: ranks_per_dim would create an empty local "
                        "domain" );
            return;
        }

        if ( requested_product != 0 )
            Info << "Ignoring ranks_per_dim because its product does not "
                    "match the MPI communicator size; selecting a "
                    "geometry-aware decomposition."
                 << std::endl;

        double best_cost = std::numeric_limits<double>::max();
        std::array<int, 3> best = { 0, 0, 0 };
        for ( int px = 1; px <= comm_size; ++px )
        {
            if ( comm_size % px != 0 || px > cells[0] )
                continue;
            const int remainder = comm_size / px;
            for ( int py = 1; py <= remainder; ++py )
            {
                if ( remainder % py != 0 || py > cells[1] )
                    continue;
                const int pz = remainder / py;
                if ( pz > cells[2] )
                    continue;

                // Approximate the total internal face area communicated by a
                // Cartesian decomposition. Constant factors are omitted.
                const double cost = ( px - 1.0 ) * cells[1] * cells[2] +
                                    ( py - 1.0 ) * cells[0] * cells[2] +
                                    ( pz - 1.0 ) * cells[0] * cells[1];
                if ( cost < best_cost )
                {
                    best_cost = cost;
                    best = { px, py, pz };
                }
            }
        }

        if ( best[0] == 0 )
            throw std::runtime_error(
                "Error: MPI communicator is too large for the global grid" );
        space.ranks_per_dim = best;
    }

    void calcAuxiliaryProperties( MPI_Comm comm )
    {
        // create auxiliary properties
        properties.thermal_diffusivity =
            ( properties.thermal_conductivity ) /
            ( properties.density * properties.specific_heat );

        time.time_step = ( time.Co * space.cell_size * space.cell_size ) /
                         ( properties.thermal_diffusivity );

        Info << "Calculated time step: " << time.time_step << std::endl;

        time.time = time.start_time;

        const double duration = time.end_time - time.start_time;
        const double step_count = duration / time.time_step;
        if ( !std::isfinite( time.time_step ) || time.time_step <= 0.0 ||
             !std::isfinite( step_count ) ||
             step_count > std::numeric_limits<int>::max() - 1.0 )
            throw std::runtime_error( "Error: calculated timestep or step "
                                      "count is not representable" );
        time.num_steps = static_cast<int>(
            std::ceil( duration / time.time_step -
                       16.0 * std::numeric_limits<double>::epsilon() ) );

        time.output.setInterval( time.num_steps );
        time.monitor.setInterval( time.num_steps );

        // initialize time monitoring
        time_monitor = TimeMonitor( comm, time );
    }

    // Calls other read input functions to initialize variables, returning a
    // list of which sections were found
    std::vector<bool> readSections( nlohmann::json db )
    {
        std::vector<bool> found_sections( 5, false );
        if ( db.contains( "time" ) )
        {
            readInputTime( db );
            found_sections[0] = true;
        }
        if ( db.contains( "space" ) )
        {
            readInputSpace( db );
            found_sections[1] = true;
        }
        if ( db.contains( "properties" ) )
        {
            readInputProperties( db );
            found_sections[2] = true;
        }
        if ( db.contains( "source" ) )
        {
            readInputSource( db );
            found_sections[3] = true;
        }
        if ( db.contains( "sampling" ) )
        {
            readInputSampling( db );
            found_sections[4] = true;
        }
        return found_sections;
    }

    void readInputTime( nlohmann::json db )
    {
        // Read time components
        time.Co = db["time"]["Co"];
        time.start_time = db["time"]["start_time"];
        time.end_time = db["time"]["end_time"];
        time.output.total_steps = db["time"]["total_output_steps"];
        time.monitor.total_steps = db["time"]["total_monitor_steps"];
    }

    void readInputSpace( nlohmann::json db )
    {
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

        space.ranks_per_dim = ranks_per_dim;
    }

    void readInputProperties( nlohmann::json db )
    {
        // Read properties components
        properties.density = db["properties"]["density"];
        properties.specific_heat = db["properties"]["specific_heat"];
        properties.thermal_conductivity =
            db["properties"]["thermal_conductivity"];
        properties.latent_heat = db["properties"]["latent_heat"];
        properties.solidus = db["properties"]["solidus"];
        properties.liquidus = db["properties"]["liquidus"];
    }

    void readInputSource( nlohmann::json db )
    {
        // Read heat source components
        source.absorption = db["source"]["absorption"];
        source.two_sigma = db["source"]["two_sigma"];

        source.scan_path_file = db["source"]["scan_path_file"];
    }

    void readInputSampling( nlohmann::json db )
    {
        // Read sampling components
        sampling = Sampling{};
        if ( db.contains( "sampling" ) )
        {
            const std::string sampling_type = db["sampling"]["type"];

            if ( sampling_type == "solidification_data" )
            {
                sampling.type = sampling_type;
                sampling.enabled = true;
            }
            else
                throw std::runtime_error( "Error: unsupported sampling type " +
                                          sampling_type );

            const std::string sampling_format =
                db["sampling"].value( "format", "default" );

            if ( sampling_format == "exaca" )
            {
                sampling.format = sampling_format;
            }
            else if ( sampling_format == "default" )
            {
                sampling.format = "default";
            }
            else
                throw std::runtime_error(
                    "Error: sampling format must be default or exaca" );

            if ( db["sampling"].contains( "directory_name" ) )
            {
                sampling.directory_name = db["sampling"]["directory_name"];
            }
        }
    }
};

} // namespace Finch

#endif
