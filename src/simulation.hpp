#ifndef Simulation_H
#define Simulation_H

#include <iostream>
#include <fstream>
#include <array>
#include <unistd.h>
#include "yaml-cpp/yaml.h"

// Info macro for writing on master
#define Info if (rank == 0) std::cout

struct Time
{
    double Co;
    double start_time;
    double end_time;
    double time_step;
    double time;
    int num_output_steps;
    int output_interval;
};

struct Space
{
    double initial_temperature;
    double cell_size;
    std::array<double, 3> global_low_corner;
    std::array<double, 3> global_high_corner;
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
    double liquidus;
};

struct Sampling
{
    std::string type;
    std::string format;
    bool enabled;
};

class Simulation
{
public:
    //SolidificationOutput solidificationOutput;
    Time time;
    Space space;
    Source source;
    Properties properties;
    Sampling sampling;

    int rank;
    int size;

    // constructor
    Simulation(MPI_Comm comm, int argc, char* argv[])
    {
        MPI_Comm_rank( comm, &rank );
        MPI_Comm_size( comm, &size );

        readInput(argc, argv);
        
        properties.thermal_diffusivity = 
            ( properties.thermal_conductivity )
          / ( properties.density * properties.specific_heat );

        time.time_step = 
            ( time.Co * space.cell_size * space.cell_size )
          / ( properties.thermal_diffusivity );


        time.time = time.start_time;

        write();

        Info << "Calculated time step: " << time.time_step << std::endl;
    }

    void write()
    {
        Info << "  Co: " << time.Co << std::endl;
        Info << "  Start Time: " << time.start_time << std::endl;
        Info << "  End Time: " << time.end_time << std::endl;
        Info << "  Num Output Steps: " << time.num_output_steps << std::endl;
        Info << "Simulation will be performed using parameters: " << std::endl;

        // Print time
        Info << "Time:" << std::endl;
        Info << "  Co: " << time.Co << std::endl;
        Info << "  Start Time: " << time.start_time << std::endl;
        Info << "  End Time: " << time.end_time << std::endl;
        Info << "  Num Output Steps: " << time.num_output_steps << std::endl;

        // Print space
        Info << "Space:" << std::endl;
        Info << "  Initial temperature: " << space.initial_temperature << std::endl;
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
        Info << "  Thermal Conductivity: " << properties.thermal_conductivity << std::endl;
        Info << "  Latent Heat: " << properties.latent_heat << std::endl;
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
        if (sampling.enabled)
        {
            Info << "  type: " << sampling.type << std::endl;
            Info << "  format:" << sampling.format << std::endl;
        }
        else
        {
            Info<< "Skipping optional sampling." << std::endl;
        }
    }

private:
    void readInput(int argc, char* argv[])
    {
        const char* filename = nullptr;
        int option;

        while ((option = getopt(argc, argv, "i:")) != -1)
        {
            if (option == 'i') 
            {
                filename = optarg;
            }
            else
            {
                std::cerr << "Usage: " << argv[0]
                          << " -i <input_yaml_file>" << std::endl;
                return;
            }
        }

        // parse input file
        {
            YAML::Node db = YAML::LoadFile(filename);

            // Read time components
            time.Co =
                db["time"]["Co"].as<double>();
            time.start_time =
                db["time"]["start_time"].as<double>();
            time.end_time =
                db["time"]["end_time"].as<double>();
            time.num_output_steps =
                db["time"]["num_output_steps"].as<int>();

            // Read space components
            space.initial_temperature =
                db["space"]["initial_temperature"].as<double>();
            space.cell_size =
                db["space"]["cell_size"].as<double>();
            space.global_low_corner =
                db["space"]["global_low_corner"].as<std::array<double, 3>>();
            space.global_high_corner =
                db["space"]["global_high_corner"].as<std::array<double, 3>>();

            // Read properties components
            properties.density =
                db["properties"]["density"].as<double>();
            properties.specific_heat =
                db["properties"]["specific_heat"].as<double>();
            properties.thermal_conductivity =
                db["properties"]["thermal_conductivity"].as<double>();
            properties.latent_heat =
                db["properties"]["latent_heat"].as<double>();
            properties.liquidus =
                db["properties"]["liquidus"].as<double>();

            // Read heat source components
            source.absorption =
                db["source"]["absorption"].as<double>();
            source.two_sigma =
                db["source"]["two_sigma"].as<std::array<double, 3>>();

            source.two_sigma[0] = fabs( source.two_sigma[0] );
            source.two_sigma[1] = fabs( source.two_sigma[1] );
            source.two_sigma[2] = fabs( source.two_sigma[2] );

            source.scan_path_file =
                db["source"]["scan_path_file"].as<std::string>();

            // Read sampling components
            sampling.enabled = false;
            if (db["sampling"])
            {
                const std::string samplingType =
                    db["sampling"]["type"].as<std::string>();

                if (samplingType == "solidification_data")
                {
                    sampling.type = samplingType;
                    sampling.enabled = true;
                }

                const std::string samplingFormat =
                    db["sampling"]["format"].as<std::string>();

                if (samplingFormat == "exaca")
                {
                    sampling.format = samplingFormat;
                }
                else
                {
                    sampling.format = "default";
                }
            }
        }
    }
};
#endif
