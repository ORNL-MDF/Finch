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

#include <Finch_BoundaryConditions.hpp>
#include <Finch_Version.hpp>

namespace Finch
{

// Rank-zero stream used for input and progress reporting.
#define FINCH_INFO                                                             \
    if ( comm_rank == 0 )                                                      \
    std::cout

struct Output
{
    int total_steps = 0;

    bool isDue( const int completed_steps, const int num_steps ) const
    {
        if ( total_steps == 0 )
            return false;

        // Integer arithmetic gives evenly distributed output and always
        // includes the final step without accumulating floating-point error.
        return ( static_cast<long long>( completed_steps ) * total_steps ) /
                   num_steps >
               ( static_cast<long long>( completed_steps - 1 ) * total_steps ) /
                   num_steps;
    }
};

struct Time
{
    // A three-dimensional forward-Euler diffusion update is stable at or below
    // 1/6. Keep some margin by default while allowing advanced users to choose
    // another stable value.
    double maximum_fourier_number = 0.125;
    double start_time = 0.0;
    double end_time = 0.0;
    double maximum_time_step = 0.0;
    double time = 0.0;
    int num_steps = 0;
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

struct AbsorptionInput
{
    std::string type;
    double coefficient = 0.0;
    std::string geometry;
    double fresnel_absorptivity = 0.0;
    double conduction_absorptivity = 0.0;
    double transition_aspect_ratio = 1.0;
};

struct Source
{
    std::string type = "gaussian";
    AbsorptionInput absorption;
    std::array<double, 3> two_sigma = { 0.0, 0.0, 0.0 };
    std::string scan_path_file;

    // Tabulated planar profile and projected axial distribution.
    std::string profile_file;
    std::string coordinate_frame = "global";
    double minimum_depth = 0.0;
    double exponent_slope = 0.0;
    double exponent_intercept = 1.0;

    // Optional explicit, one-step-lagged isotherm-depth feedback.
    bool transient_depth = false;
    std::string depth_temperature = "liquidus";
};

struct TemperatureProperty
{
    double constant = 0.0;
    std::vector<double> temperature;
    std::vector<double> values;

    bool isTabulated() const { return !temperature.empty(); }

    bool isConstant() const
    {
        if ( !isTabulated() )
            return true;
        return !values.empty() &&
               std::all_of( values.begin(), values.end(),
                            [&]( const double entry )
                            { return entry == values.front(); } );
    }

    double constantValue() const
    {
        return isTabulated() ? values.front() : constant;
    }

    double value( const double query_temperature ) const
    {
        if ( !isTabulated() )
            return constant;
        if ( query_temperature <= temperature.front() )
            return values.front();
        if ( query_temperature >= temperature.back() )
            return values.back();

        const auto upper = std::upper_bound(
            temperature.begin(), temperature.end(), query_temperature );
        const std::size_t high =
            static_cast<std::size_t>( upper - temperature.begin() );
        const std::size_t low = high - 1;
        const double fraction = ( query_temperature - temperature[low] ) /
                                ( temperature[high] - temperature[low] );
        return values[low] + fraction * ( values[high] - values[low] );
    }

    double minimumValue() const
    {
        if ( !isTabulated() )
            return constant;
        return *std::min_element( values.begin(), values.end() );
    }

    double maximumValue() const
    {
        if ( !isTabulated() )
            return constant;
        return *std::max_element( values.begin(), values.end() );
    }
};

struct Properties
{
    double density = 0.0;
    TemperatureProperty specific_heat;
    TemperatureProperty thermal_conductivity;
    // Conservative nonlinear-stencil bound. Conductivity and heat capacity
    // can occur at different temperatures on neighboring nodes, so this is
    // k_max / C_min rather than max(k/C).
    double thermal_diffusivity_upper_bound = 0.0;
    double latent_heat = 0.0;
    double solidus = 0.0;
    double liquidus = 0.0;
    static constexpr int lookup_table_points = 256;

    bool isTemperatureDependent() const
    {
        return !specific_heat.isConstant() ||
               !thermal_conductivity.isConstant();
    }
};

enum class FunctionControl
{
    every_step,
    output_count,
    execute,
    end
};

struct FunctionSchedule
{
    FunctionControl control = FunctionControl::end;
    int count = 0;
};

struct SolidificationDataInput
{
    std::string name;
    std::string type;
    std::string format = "default";
    std::string directory = "solidification";
    FunctionSchedule execute;
    FunctionSchedule write;
    bool enabled = false;
};

struct MeltPoolDimensionsInput
{
    std::string name;
    bool enabled = false;
    FunctionSchedule execute;
    FunctionSchedule write;
    std::array<bool, 2> isotherms = { true, true };
    std::string coordinate_frame = "global";
    std::string directory = "melt_pool_dimensions";
};

enum class FieldOutputField
{
    temperature,
    volumetric_heat_source
};

inline const char* fieldOutputFieldName( const FieldOutputField field )
{
    return field == FieldOutputField::temperature ? "temperature"
                                                  : "volumetric_heat_source";
}

struct FieldOutputInput
{
    std::string name;
    std::string format = "bov";
    bool enabled = false;
    FunctionSchedule execute;
    std::vector<FieldOutputField> fields = { FieldOutputField::temperature };
};

struct Functions
{
    SolidificationDataInput solidification;
    MeltPoolDimensionsInput melt_pool_dimensions;
    FieldOutputInput field_output;
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

    void setNumSteps( const int steps ) { num_steps = steps; }

    void write( int step )
    {
        update();

        FINCH_INFO << "Time Step: " << step << "/" << num_steps << ", "
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
    BoundaryConditions boundary;
    Functions functions;
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
        FINCH_INFO << "Finch version: " << version() << " (" << commitHash()
                   << ")" << std::endl;
        FINCH_INFO << "Simulation will be performed using parameters: "
                   << std::endl;

        // Print time
        FINCH_INFO << "Time:" << std::endl;
        FINCH_INFO << "  Maximum Fourier Number: "
                   << time.maximum_fourier_number << std::endl;
        FINCH_INFO << "  Start Time: " << time.start_time << std::endl;
        FINCH_INFO << "  End Time: " << time.end_time << std::endl;
        FINCH_INFO << "  Num Monitor Steps: " << time.monitor.total_steps
                   << std::endl;

        // Print space
        FINCH_INFO << "Space:" << std::endl;
        FINCH_INFO << "  Initial temperature: " << space.initial_temperature
                   << std::endl;
        FINCH_INFO << "  Cell Size: " << space.cell_size << std::endl;
        FINCH_INFO << "  Global Low Corner:" << std::endl;
        FINCH_INFO << "    X: " << space.global_low_corner[0] << std::endl;
        FINCH_INFO << "    Y: " << space.global_low_corner[1] << std::endl;
        FINCH_INFO << "    Z: " << space.global_low_corner[2] << std::endl;
        FINCH_INFO << "  Global High Corner:" << std::endl;
        FINCH_INFO << "    X: " << space.global_high_corner[0] << std::endl;
        FINCH_INFO << "    Y: " << space.global_high_corner[1] << std::endl;
        FINCH_INFO << "    Z: " << space.global_high_corner[2] << std::endl;

        // Print properties
        FINCH_INFO << "Properties:" << std::endl;
        FINCH_INFO << "  Density: " << properties.density << std::endl;
        if ( properties.specific_heat.isTabulated() )
        {
            FINCH_INFO << "  Specific Heat: "
                       << properties.specific_heat.values.size()
                       << " input points from "
                       << properties.specific_heat.temperature.front() << " to "
                       << properties.specific_heat.temperature.back() << " K"
                       << std::endl;
        }
        else
        {
            FINCH_INFO << "  Specific Heat: "
                       << properties.specific_heat.constant << std::endl;
        }
        if ( properties.thermal_conductivity.isTabulated() )
        {
            FINCH_INFO << "  Thermal Conductivity: "
                       << properties.thermal_conductivity.values.size()
                       << " input points from "
                       << properties.thermal_conductivity.temperature.front()
                       << " to "
                       << properties.thermal_conductivity.temperature.back()
                       << " K" << std::endl;
        }
        else
        {
            FINCH_INFO << "  Thermal Conductivity: "
                       << properties.thermal_conductivity.constant << std::endl;
        }
        FINCH_INFO << "  Latent Heat: " << properties.latent_heat << std::endl;
        FINCH_INFO << "  Solidus: " << properties.solidus << std::endl;
        FINCH_INFO << "  Liquidus: " << properties.liquidus << std::endl;

        FINCH_INFO << "Boundary conditions:" << std::endl;
        for ( int face = 0; face < 6; ++face )
        {
            const auto& condition = boundary.faces[face];
            FINCH_INFO << "  " << boundary_face_names[face] << ": "
                       << condition.type;
            if ( condition.type == "dirichlet" || condition.type == "neumann" )
            {
                FINCH_INFO << " (" << condition.value << ")";
            }
            else if ( condition.type == "convection_radiation" )
            {
                FINCH_INFO << " (h=" << condition.convection_coefficient
                           << ", emissivity=" << condition.emissivity
                           << ", ambient=" << condition.ambient_temperature
                           << " K)";
            }
            FINCH_INFO << std::endl;
        }

        // Print source
        FINCH_INFO << "Source:" << std::endl;
        FINCH_INFO << "  Type: " << source.type << std::endl;
        FINCH_INFO << "  Absorption: " << source.absorption.type;
        if ( source.absorption.type == "constant" )
        {
            FINCH_INFO << " (coefficient=" << source.absorption.coefficient
                       << ")";
        }
        else
        {
            FINCH_INFO << " (geometry=" << source.absorption.geometry
                       << ", fresnel=" << source.absorption.fresnel_absorptivity
                       << ", conduction="
                       << source.absorption.conduction_absorptivity
                       << ", transition aspect ratio="
                       << source.absorption.transition_aspect_ratio << ")";
        }
        FINCH_INFO << std::endl;
        if ( source.type == "gaussian" )
        {
            FINCH_INFO << "  two-sigma:" << std::endl;
            FINCH_INFO << "    X: " << source.two_sigma[0] << std::endl;
            FINCH_INFO << "    Y: " << source.two_sigma[1] << std::endl;
            FINCH_INFO << "    Z: " << source.two_sigma[2] << std::endl;
        }
        else
        {
            FINCH_INFO << "  Profile file: " << source.profile_file
                       << std::endl;
            FINCH_INFO << "  Coordinate frame: " << source.coordinate_frame
                       << std::endl;
            FINCH_INFO << "  Minimum depth: " << source.minimum_depth
                       << std::endl;
            FINCH_INFO << "  Axial exponent slope/intercept: "
                       << source.exponent_slope << "/"
                       << source.exponent_intercept << std::endl;
            if ( source.transient_depth )
                FINCH_INFO << "  Transient depth: " << source.depth_temperature
                           << " isotherm" << std::endl;
        }
        FINCH_INFO << "  scan path file: " << source.scan_path_file
                   << std::endl;

        FINCH_INFO << "Functions:" << std::endl;
        if ( functions.solidification.enabled )
        {
            FINCH_INFO
                << "  " << functions.solidification.name
                << ": solidification_data every step, write at end, format "
                << functions.solidification.format << ", directory "
                << functions.solidification.directory << std::endl;
        }
        if ( functions.melt_pool_dimensions.enabled )
        {
            FINCH_INFO << "  " << functions.melt_pool_dimensions.name
                       << ": melt_pool_dimensions, "
                       << functions.melt_pool_dimensions.execute.count
                       << " executions in the "
                       << functions.melt_pool_dimensions.coordinate_frame
                       << " frame, directory "
                       << functions.melt_pool_dimensions.directory << std::endl;
        }
        if ( functions.field_output.enabled )
        {
            FINCH_INFO << "  " << functions.field_output.name
                       << ": field_output, "
                       << functions.field_output.execute.count
                       << " executions, format "
                       << functions.field_output.format << ", fields";
            for ( const auto field : functions.field_output.fields )
                FINCH_INFO << ' ' << fieldOutputFieldName( field );
            FINCH_INFO << std::endl;
        }
        if ( !functions.solidification.enabled &&
             !functions.melt_pool_dimensions.enabled &&
             !functions.field_output.enabled )
            FINCH_INFO << "  None" << std::endl;
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
        // Time, space, properties, and source are required. Boundary and
        // functions are optional.
        return std::all_of( found.begin(), found.begin() + 4,
                            []( const bool value ) { return value; } );
    }

    void parseInputFile( MPI_Comm comm, const std::string& filename,
                         const int input_file_number = 0 )
    {
        // Input file is either a Finch input file or an ExaCA input file with a
        // Finch object
        FINCH_INFO << "Parsing input file " << input_file_number << std::endl;
        nlohmann::json input_data_raw = readInputDocument( comm, filename );
        if ( !input_data_raw.contains( "Finch" ) )
        {
            // This is a standalone Finch input file.
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
                "time",   "space",    "properties",
                "source", "boundary", "functions" };
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
                        FINCH_INFO << "Warning: Finch input object "
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

        const auto validate_temperature_property =
            [&]( const TemperatureProperty& property, const std::string& name )
        {
            if ( !property.isTabulated() )
            {
                if ( !finite( property.constant ) || property.constant <= 0.0 )
                    throw std::runtime_error( "Error: properties." + name +
                                              " must be finite and positive" );
                return;
            }

            if ( property.temperature.size() < 2 ||
                 property.temperature.size() != property.values.size() )
                throw std::runtime_error(
                    "Error: tabulated properties." + name +
                    " requires equally sized temperature and values "
                    "arrays with at least two entries" );
            for ( std::size_t i = 0; i < property.temperature.size(); ++i )
            {
                if ( !finite( property.temperature[i] ) ||
                     !finite( property.values[i] ) ||
                     property.values[i] <= 0.0 )
                    throw std::runtime_error(
                        "Error: tabulated properties." + name +
                        " temperatures must be finite and values must be "
                        "finite and positive" );
                if ( i > 0 &&
                     property.temperature[i] <= property.temperature[i - 1] )
                    throw std::runtime_error(
                        "Error: tabulated properties." + name +
                        " temperatures must be strictly increasing" );
            }
        };

        if ( !finite( time.maximum_fourier_number ) ||
             time.maximum_fourier_number <= 0.0 ||
             time.maximum_fourier_number > 1.0 / 6.0 )
            throw std::runtime_error(
                "Error: time.maximum_fourier_number must be in (0, 1/6] for "
                "the 3D explicit diffusion stencil" );
        if ( !finite( time.start_time ) || !finite( time.end_time ) ||
             time.end_time <= time.start_time )
            throw std::runtime_error(
                "Error: end_time must be finite and greater than start_time" );
        if ( time.monitor.total_steps < 0 )
            throw std::runtime_error(
                "Error: monitor step count cannot be negative" );

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
             !finite( properties.latent_heat ) ||
             properties.latent_heat < 0.0 || !finite( properties.solidus ) ||
             !finite( properties.liquidus ) ||
             properties.liquidus <= properties.solidus )
            throw std::runtime_error(
                "Error: density must be finite and positive, latent heat must "
                "be finite and nonnegative, and finite liquidus must exceed "
                "solidus" );
        validate_temperature_property( properties.specific_heat,
                                       "specific_heat" );
        validate_temperature_property( properties.thermal_conductivity,
                                       "thermal_conductivity" );
        if ( source.scan_path_file.empty() )
            throw std::runtime_error(
                "Error: source.scan_path_file cannot be empty" );

        const auto valid_absorptivity = [&]( const double value )
        { return finite( value ) && value >= 0.0 && value <= 1.0; };
        if ( source.absorption.type == "constant" )
        {
            if ( !valid_absorptivity( source.absorption.coefficient ) )
                throw std::runtime_error(
                    "Error: constant absorption coefficient must be in "
                    "[0,1]" );
        }
        else if ( source.absorption.type == "kelly" )
        {
            if ( source.absorption.geometry != "cone" &&
                 source.absorption.geometry != "cylinder" )
                throw std::runtime_error(
                    "Error: Kelly absorption geometry must be cone or "
                    "cylinder" );
            if ( !finite( source.absorption.fresnel_absorptivity ) ||
                 source.absorption.fresnel_absorptivity <= 0.0 ||
                 source.absorption.fresnel_absorptivity > 1.0 ||
                 !valid_absorptivity(
                     source.absorption.conduction_absorptivity ) ||
                 !finite( source.absorption.transition_aspect_ratio ) ||
                 source.absorption.transition_aspect_ratio < 0.0 )
                throw std::runtime_error(
                    "Error: Kelly absorptivities must be in [0,1] with "
                    "positive fresnel_absorptivity, and "
                    "transition_aspect_ratio must be nonnegative" );
            if ( source.type != "tabulated" || !source.transient_depth ||
                 source.depth_temperature != "liquidus" )
                throw std::runtime_error(
                    "Error: Kelly absorption requires a tabulated source "
                    "with liquidus transient_depth" );
        }
        else
            throw std::runtime_error(
                "Error: source.absorption.type must be constant or kelly" );

        if ( source.type == "gaussian" )
        {
            for ( const double sigma : source.two_sigma )
                if ( !finite( sigma ) || sigma <= 0.0 )
                    throw std::runtime_error(
                        "Error: every source two_sigma value must be "
                        "positive" );
        }
        else if ( source.type == "tabulated" )
        {
            if ( source.profile_file.empty() )
                throw std::runtime_error(
                    "Error: tabulated source profile_file cannot be empty" );
            if ( source.coordinate_frame != "global" &&
                 source.coordinate_frame != "scan_path" )
                throw std::runtime_error(
                    "Error: source coordinate_frame must be global or "
                    "scan_path" );
            if ( !finite( source.minimum_depth ) ||
                 source.minimum_depth <= 0.0 ||
                 !finite( source.exponent_slope ) ||
                 !finite( source.exponent_intercept ) )
                throw std::runtime_error(
                    "Error: tabulated source minimum_depth must be positive "
                    "and axial-profile coefficients must be finite" );
            if ( source.transient_depth )
            {
                if ( source.depth_temperature != "solidus" &&
                     source.depth_temperature != "liquidus" )
                    throw std::runtime_error(
                        "Error: transient_depth.temperature must be solidus "
                        "or liquidus" );
            }
        }
        else
            throw std::runtime_error(
                "Error: source type must be gaussian or tabulated" );

        if ( functions.solidification.enabled &&
             functions.solidification.directory.empty() )
            throw std::runtime_error(
                "Error: solidification_data directory cannot be empty" );

        if ( functions.melt_pool_dimensions.enabled )
        {
            const auto& melt_pool = functions.melt_pool_dimensions;
            if ( melt_pool.execute.count <= 0 )
                throw std::runtime_error(
                    "Error: melt_pool_dimensions output count must be "
                    "positive" );
            if ( melt_pool.coordinate_frame != "global" &&
                 melt_pool.coordinate_frame != "scan_path" )
                throw std::runtime_error(
                    "Error: melt_pool_dimensions.coordinate_frame must be "
                    "global or scan_path" );
            if ( melt_pool.directory.empty() )
                throw std::runtime_error(
                    "Error: melt_pool_dimensions directory cannot be "
                    "empty" );
        }

        if ( functions.field_output.enabled )
        {
            if ( functions.field_output.execute.count <= 0 )
                throw std::runtime_error(
                    "Error: field_output output count must be positive" );
            if ( functions.field_output.format != "bov" &&
                 functions.field_output.format != "adios2" )
                throw std::runtime_error(
                    "Error: field_output format must be bov or adios2" );
            if ( functions.field_output.format == "bov" &&
                 std::any_of( functions.field_output.fields.begin(),
                              functions.field_output.fields.end(),
                              []( const auto field ) {
                                  return field != FieldOutputField::temperature;
                              } ) )
                throw std::runtime_error(
                    "Error: BOV field_output supports only temperature; use "
                    "adios2 for derived fields" );
#if !Finch_ENABLE_ADIOS2
            if ( functions.field_output.format == "adios2" )
                throw std::runtime_error(
                    "Error: field_output format adios2 requires a Finch "
                    "build with Finch_ENABLE_ADIOS2=ON" );
#endif
        }

        for ( int face = 0; face < 6; ++face )
        {
            const auto& condition = boundary.faces[face];
            const std::string prefix =
                "Error: boundary." + std::string( boundary_face_names[face] ) +
                " ";
            if ( condition.type == "dirichlet" || condition.type == "neumann" )
            {
                if ( !finite( condition.value ) )
                    throw std::runtime_error( prefix + "value must be finite" );
            }
            else if ( condition.type == "convection_radiation" )
            {
                if ( !finite( condition.convection_coefficient ) ||
                     condition.convection_coefficient < 0.0 )
                    throw std::runtime_error(
                        prefix + "h must be finite and nonnegative" );
                if ( !finite( condition.emissivity ) ||
                     condition.emissivity < 0.0 || condition.emissivity > 1.0 )
                    throw std::runtime_error( prefix +
                                              "emissivity must be in [0,1]" );
                if ( !finite( condition.ambient_temperature ) ||
                     condition.ambient_temperature <= 0.0 )
                    throw std::runtime_error(
                        prefix +
                        "ambient_temperature must be positive and finite" );
            }
            else if ( condition.type != "adiabatic" )
                throw std::runtime_error(
                    prefix + "type must be adiabatic, dirichlet, neumann, or "
                             "convection_radiation" );
        }
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
            FINCH_INFO << "Ignoring ranks_per_dim because its product does not "
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
        // Use a conservative fixed stability bound across the property range.
        // This avoids a device reduction and MPI collective every timestep.
        const double minimum_heat_capacity =
            properties.specific_heat.minimumValue();
        const double maximum_conductivity =
            properties.thermal_conductivity.maximumValue();
        properties.thermal_diffusivity_upper_bound =
            maximum_conductivity /
            ( properties.density * minimum_heat_capacity );

        time.maximum_time_step = ( time.maximum_fourier_number *
                                   space.cell_size * space.cell_size ) /
                                 properties.thermal_diffusivity_upper_bound;

        FINCH_INFO << "Maximum diffusion-stable time step: "
                   << time.maximum_time_step << std::endl;

        time.time = time.start_time;

        const double duration = time.end_time - time.start_time;
        const double step_count = duration / time.maximum_time_step;
        if ( !std::isfinite( time.maximum_time_step ) ||
             time.maximum_time_step <= 0.0 || !std::isfinite( step_count ) ||
             step_count > std::numeric_limits<int>::max() - 1.0 )
            throw std::runtime_error( "Error: calculated timestep or step "
                                      "count is not representable" );
        const double step_tolerance = 128.0 *
                                      std::numeric_limits<double>::epsilon() *
                                      std::max( 1.0, std::abs( step_count ) );
        time.num_steps = std::max(
            1, static_cast<int>( std::ceil( step_count - step_tolerance ) ) );

        // initialize time monitoring
        time_monitor = TimeMonitor( comm, time );
    }

    // Calls other read input functions to initialize variables, returning a
    // list of which sections were found
    std::vector<bool> readSections( nlohmann::json db )
    {
        std::vector<bool> found_sections( 6, false );
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
        if ( db.contains( "boundary" ) )
        {
            readInputBoundary( db );
            found_sections[4] = true;
        }
        if ( db.contains( "functions" ) )
        {
            readInputFunctions( db );
            found_sections[5] = true;
        }
        return found_sections;
    }

    void readInputTime( nlohmann::json db )
    {
        // Read time components
        const auto& input = db.at( "time" );
        time.maximum_fourier_number =
            input.value( "maximum_fourier_number", 0.125 );
        time.start_time = input.at( "start_time" );
        time.end_time = input.at( "end_time" );
        time.monitor.total_steps = input.value( "total_monitor_steps", 0 );
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
        properties.specific_heat =
            readTemperatureProperty( db["properties"], "specific_heat" );
        properties.thermal_conductivity =
            readTemperatureProperty( db["properties"], "thermal_conductivity" );
        properties.latent_heat = db["properties"]["latent_heat"];
        properties.solidus = db["properties"]["solidus"];
        properties.liquidus = db["properties"]["liquidus"];
    }

    TemperatureProperty readTemperatureProperty( const nlohmann::json& input,
                                                 const std::string& name )
    {
        const auto& value = input.at( name );
        TemperatureProperty property;
        if ( value.is_number() )
        {
            property.constant = value.get<double>();
            return property;
        }
        if ( !value.is_object() )
            throw std::runtime_error( "Error: properties." + name +
                                      " must be a number or table object" );

        if ( !value.contains( "temperature" ) || !value.contains( "values" ) )
            throw std::runtime_error(
                "Error: tabulated properties." + name +
                " requires temperature and values arrays" );

        property.temperature =
            value.at( "temperature" ).get<std::vector<double>>();
        property.values = value.at( "values" ).get<std::vector<double>>();
        return property;
    }

    void readInputSource( nlohmann::json db )
    {
        // Read heat source components
        source = Source{};
        const auto& input = db.at( "source" );
        source.type = input.value( "type", "gaussian" );
        source.scan_path_file = input.at( "scan_path_file" );

        const auto& absorption = input.at( "absorption" );
        if ( !absorption.is_object() )
            throw std::runtime_error(
                "Error: source.absorption must be an object" );
        source.absorption.type = absorption.at( "type" );
        if ( source.absorption.type == "constant" )
            source.absorption.coefficient = absorption.at( "coefficient" );
        else if ( source.absorption.type == "kelly" )
        {
            source.absorption.geometry = absorption.at( "geometry" );
            source.absorption.fresnel_absorptivity =
                absorption.at( "fresnel_absorptivity" );
            source.absorption.conduction_absorptivity =
                absorption.at( "conduction_absorptivity" );
            source.absorption.transition_aspect_ratio =
                absorption.value( "transition_aspect_ratio", 1.0 );
        }

        if ( source.type == "gaussian" )
            source.two_sigma =
                input.at( "two_sigma" ).get<std::array<double, 3>>();
        else if ( source.type == "tabulated" )
        {
            source.profile_file = input.at( "profile_file" ).get<std::string>();
            source.coordinate_frame =
                input.value( "coordinate_frame", "global" );
            source.minimum_depth = input.at( "minimum_depth" );

            const auto& axial = input.at( "axial_profile" );
            source.exponent_slope = axial.at( "exponent_slope" );
            source.exponent_intercept = axial.at( "exponent_intercept" );

            if ( input.contains( "transient_depth" ) )
            {
                const auto& transient = input.at( "transient_depth" );
                if ( !transient.is_object() )
                    throw std::runtime_error(
                        "Error: source.transient_depth must be an object" );
                source.transient_depth = true;
                source.depth_temperature =
                    transient.value( "temperature", "liquidus" );
            }
        }
    }

    FunctionSchedule readFunctionSchedule( const nlohmann::json& function,
                                           const std::string& key,
                                           const FunctionControl fallback,
                                           const std::string& context )
    {
        FunctionSchedule result;
        result.control = fallback;
        if ( !function.contains( key ) )
            return result;

        const auto& schedule = function.at( key );
        if ( !schedule.is_object() )
            throw std::runtime_error( "Error: " + context + "." + key +
                                      " must be an object" );

        const std::string control = schedule.at( "control" );
        if ( control == "every_step" )
            result.control = FunctionControl::every_step;
        else if ( control == "output_count" )
        {
            result.control = FunctionControl::output_count;
            result.count = schedule.at( "count" );
        }
        else if ( control == "execute" )
            result.control = FunctionControl::execute;
        else if ( control == "end" )
            result.control = FunctionControl::end;
        else
            throw std::runtime_error( "Error: unsupported control " + control +
                                      " in " + context + "." + key );
        return result;
    }

    void readInputBoundary( const nlohmann::json& db )
    {
        boundary = BoundaryConditions{};
        const auto& input = db.at( "boundary" );
        if ( !input.is_object() )
            throw std::runtime_error( "Error: boundary must be a JSON object" );

        for ( int face = 0; face < 6; ++face )
        {
            const char* name = boundary_face_names[face];
            if ( !input.contains( name ) )
                continue;

            const auto& face_input = input.at( name );
            if ( !face_input.is_object() )
                throw std::runtime_error( "Error: boundary." +
                                          std::string( name ) +
                                          " must be a JSON object" );

            auto& condition = boundary.faces[face];
            condition.type = face_input.at( "type" ).get<std::string>();

            if ( condition.type == "dirichlet" )
                condition.value = face_input.at( "value" ).get<double>();
            else if ( condition.type == "neumann" )
                condition.value = face_input.at( "gradient" ).get<double>();
            else if ( condition.type == "convection_radiation" )
            {
                condition.convection_coefficient =
                    face_input.at( "h" ).get<double>();
                condition.emissivity =
                    face_input.at( "emissivity" ).get<double>();
                condition.ambient_temperature =
                    face_input.at( "ambient_temperature" ).get<double>();
            }
        }
    }

    void readInputFunctions( const nlohmann::json& db )
    {
        functions = Functions{};
        const auto& input = db.at( "functions" );
        if ( !input.is_object() )
            throw std::runtime_error(
                "Error: functions must be a JSON object" );

        for ( const auto& item : input.items() )
        {
            const std::string& name = item.key();
            const auto& value = item.value();
            if ( !value.is_object() )
                throw std::runtime_error( "Error: functions." + name +
                                          " must be an object" );

            const std::string context = "functions." + name;
            const std::string type = value.at( "type" );
            if ( type == "solidification_data" )
            {
                auto& function = functions.solidification;
                if ( function.enabled )
                    throw std::runtime_error(
                        "Error: functions may contain one "
                        "solidification_data entry" );
                function.enabled = true;
                function.name = name;
                function.type = type;
                function.execute = readFunctionSchedule(
                    value, "execute", FunctionControl::every_step, context );
                function.write = readFunctionSchedule(
                    value, "write", FunctionControl::end, context );
                if ( function.execute.control != FunctionControl::every_step ||
                     function.write.control != FunctionControl::end )
                    throw std::runtime_error(
                        "Error: solidification_data execute control must be "
                        "every_step and write control must be end" );
                function.format = value.value( "format", "default" );
                if ( function.format != "default" &&
                     function.format != "exaca" )
                    throw std::runtime_error(
                        "Error: solidification_data format must be default "
                        "or exaca" );
                function.directory =
                    value.value( "directory", "solidification" );
            }
            else if ( type == "melt_pool_dimensions" )
            {
                auto& function = functions.melt_pool_dimensions;
                if ( function.enabled )
                    throw std::runtime_error(
                        "Error: functions may contain one "
                        "melt_pool_dimensions entry" );
                function.enabled = true;
                function.name = name;
                function.execute = readFunctionSchedule(
                    value, "execute", FunctionControl::output_count, context );
                function.write = readFunctionSchedule(
                    value, "write", FunctionControl::execute, context );
                if ( function.execute.control !=
                         FunctionControl::output_count ||
                     function.write.control != FunctionControl::execute )
                    throw std::runtime_error(
                        "Error: melt_pool_dimensions execute control must be "
                        "output_count and write control must be execute" );
                function.isotherms = { false, false };
                const auto isotherms = value.value(
                    "isotherms",
                    std::vector<std::string>{ "solidus", "liquidus" } );
                for ( const auto& isotherm : isotherms )
                {
                    if ( isotherm == "solidus" )
                        function.isotherms[0] = true;
                    else if ( isotherm == "liquidus" )
                        function.isotherms[1] = true;
                    else
                        throw std::runtime_error(
                            "Error: melt_pool_dimensions isotherms must be "
                            "solidus or liquidus" );
                }
                if ( !function.isotherms[0] && !function.isotherms[1] )
                    throw std::runtime_error(
                        "Error: melt_pool_dimensions requires at least one "
                        "isotherm" );
                function.coordinate_frame =
                    value.value( "coordinate_frame", "global" );
                function.directory =
                    value.value( "directory", "melt_pool_dimensions" );
            }
            else if ( type == "field_output" )
            {
                auto& function = functions.field_output;
                if ( function.enabled )
                    throw std::runtime_error(
                        "Error: functions may contain one field_output entry" );
                function.enabled = true;
                function.name = name;
                function.execute = readFunctionSchedule(
                    value, "execute", FunctionControl::output_count, context );
                if ( function.execute.control != FunctionControl::output_count )
                    throw std::runtime_error(
                        "Error: field_output execute control must be "
                        "output_count" );
                if ( value.contains( "write" ) &&
                     readFunctionSchedule( value, "write",
                                           FunctionControl::execute, context )
                             .control != FunctionControl::execute )
                    throw std::runtime_error(
                        "Error: field_output writes when it executes" );
                const auto fields = value.value(
                    "fields", std::vector<std::string>{ "temperature" } );
                function.fields.clear();
                for ( const auto& field : fields )
                {
                    if ( field == "temperature" )
                        function.fields.push_back(
                            FieldOutputField::temperature );
                    else if ( field == "volumetric_heat_source" )
                        function.fields.push_back(
                            FieldOutputField::volumetric_heat_source );
                    else
                        throw std::runtime_error(
                            "Error: unsupported field_output field " + field );
                }
                if ( function.fields.empty() )
                    throw std::runtime_error(
                        "Error: field_output fields cannot be empty" );
#if Finch_ENABLE_ADIOS2
                const std::string default_format = "adios2";
#else
                const std::string default_format = "bov";
#endif
                function.format = value.value( "format", default_format );
            }
            else
                throw std::runtime_error( "Error: unsupported function type " +
                                          type );
        }
    }
};

} // namespace Finch

#endif
