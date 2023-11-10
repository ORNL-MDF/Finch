#ifndef Simulation_H
#define Simulation_H

#include "yaml-cpp/yaml.h"
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unistd.h>

// Info macro for writing on master
#define Info                                                                   \
    if ( comm_rank == 0 )                                                      \
    std::cout

struct Time
{
    double Co;
    double start_time;
    double end_time;
    double time_step;
    double time;
    int total_output_steps;
    int total_monitor_steps;
    int num_steps;
    int output_interval;
    int monitor_interval;
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
        total_monitor_steps = time.total_monitor_steps;
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

class Simulation
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

    // constructor
    Simulation( MPI_Comm comm, int argc, char* argv[] )
    {
        MPI_Comm_rank( comm, &comm_rank );
        MPI_Comm_size( comm, &comm_size );

        readInput( argc, argv );

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

        time.output_interval =
            static_cast<int>( ( time.num_steps / time.total_output_steps ) );

        time.output_interval =
            std::max( std::min( time.output_interval, time.num_steps ), 1 );

        time.monitor_interval =
            static_cast<int>( ( time.num_steps / time.total_monitor_steps ) );

        time.monitor_interval =
            std::max( std::min( time.monitor_interval, time.num_steps ), 1 );

        // initialize time monitoring
        time_monitor = TimeMonitor( comm, time );
    }

    void write()
    {
        Info << "  Co: " << time.Co << std::endl;
        Info << "  Start Time: " << time.start_time << std::endl;
        Info << "  End Time: " << time.end_time << std::endl;
        Info << "  Num Output Steps: " << time.total_output_steps << std::endl;
        Info << "  Num Monitor Steps: " << time.total_monitor_steps
             << std::endl;
        Info << "Simulation will be performed using parameters: " << std::endl;

        // Print time
        Info << "Time:" << std::endl;
        Info << "  Co: " << time.Co << std::endl;
        Info << "  Start Time: " << time.start_time << std::endl;
        Info << "  End Time: " << time.end_time << std::endl;
        Info << "  Num Output Steps: " << time.total_output_steps << std::endl;
        Info << "  Num Monitor Steps: " << time.total_monitor_steps
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
    void readInput( int argc, char* argv[] )
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
                std::cerr << "Usage: " << argv[0] << " -i <input_yaml_file>"
                          << std::endl;
                return;
            }
        }

        // parse input file
        {
            YAML::Node db = YAML::LoadFile( filename );

            // Read time components
            time.Co = db["time"]["Co"].as<double>();
            time.start_time = db["time"]["start_time"].as<double>();
            time.end_time = db["time"]["end_time"].as<double>();
            time.total_output_steps =
                db["time"]["total_output_steps"].as<int>();
            time.total_monitor_steps =
                db["time"]["total_monitor_steps"].as<int>();

            // Read space components
            space.initial_temperature =
                db["space"]["initial_temperature"].as<double>();
            space.cell_size = db["space"]["cell_size"].as<double>();
            space.global_low_corner =
                db["space"]["global_low_corner"].as<std::array<double, 3>>();
            space.global_high_corner =
                db["space"]["global_high_corner"].as<std::array<double, 3>>();

            /*
              Default block partitioner. This relies on MPI_Cart_create to
              balance the number of ranks in each direction. This partitioning
              is best only in the global mesh is a uniform cube.
            */
            std::array<int, 3> default_ranks_per_dim = { 0, 0, 0 };

            std::array<int, 3> ranks_per_dim =
                db["space"]["ranks_per_dim"].as<std::array<int, 3>>(
                    default_ranks_per_dim );

            // Invalid partition strategy selected. Use Default block partioner.
            int product =
                ranks_per_dim[0] * ranks_per_dim[1] * ranks_per_dim[2];

            if ( product != comm_size )
            {
                ranks_per_dim = default_ranks_per_dim;
            }

            space.ranks_per_dim = ranks_per_dim;

            // Read properties components
            properties.density = db["properties"]["density"].as<double>();
            properties.specific_heat =
                db["properties"]["specific_heat"].as<double>();
            properties.thermal_conductivity =
                db["properties"]["thermal_conductivity"].as<double>();
            properties.latent_heat =
                db["properties"]["latent_heat"].as<double>();
            properties.solidus = db["properties"]["solidus"].as<double>();
            properties.liquidus = db["properties"]["liquidus"].as<double>();

            // Read heat source components
            source.absorption = db["source"]["absorption"].as<double>();
            source.two_sigma =
                db["source"]["two_sigma"].as<std::array<double, 3>>();

            source.two_sigma[0] = fabs( source.two_sigma[0] );
            source.two_sigma[1] = fabs( source.two_sigma[1] );
            source.two_sigma[2] = fabs( source.two_sigma[2] );

            source.scan_path_file =
                db["source"]["scan_path_file"].as<std::string>();

            // Read sampling components
            sampling.enabled = false;
            if ( db["sampling"] )
            {
                const std::string sampling_type =
                    db["sampling"]["type"].as<std::string>();

                if ( sampling_type == "solidification_data" )
                {
                    sampling.type = sampling_type;
                    sampling.enabled = true;
                }

                const std::string sampling_format =
                    db["sampling"]["format"].as<std::string>();

                if ( sampling_format == "exaca" )
                {
                    sampling.format = sampling_format;
                }
                else
                {
                    sampling.format = "default";
                }

                if ( db["sampling"]["directory_name"] )
                {
                    sampling.directory_name =
                        db["sampling"]["directory_name"].as<std::string>();
                }
            }
        }
    }
};
#endif
