#include <iostream>
#include <fstream>
#include <array>
#include <unistd.h>
#include "yaml-cpp/yaml.h"

struct Time
{
    double Co;
    double start_time;
    double end_time;
    double time_step;
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
};

struct Properties
{
    double density;
    double specific_heat;
    double thermal_conductivity;
    double thermal_diffusivity;
    double latent_heat;
};

class Simulation
{
public:
    Time time;
    Space space;
    Source source;
    Properties properties;

    // constructor
    Simulation(int argc, char* argv[])
    {
        readInput(argc, argv);
        
        properties.thermal_diffusivity = 
            ( properties.thermal_conductivity )
          / ( properties.density * properties.specific_heat );

        time.time_step = 
            ( time.Co * space.cell_size * space.cell_size )
          / ( properties.thermal_diffusivity );

        write();

        std::cout << "Calculated time step: " << time.time_step << std::endl;
    }

    void write()
    {
        std::cout << "Simulation will be performed using parameters: " << std::endl;

        // Print time
        std::cout << "Time:" << std::endl;
        std::cout << "  Co: " << time.Co << std::endl;
        std::cout << "  Start Time: " << time.start_time << std::endl;
        std::cout << "  End Time: " << time.end_time << std::endl;
        std::cout << "  Num Output Steps: " << time.num_output_steps << std::endl;

        // Print space
        std::cout << "Space:" << std::endl;
        std::cout << "  Initial temperature: " << space.initial_temperature << std::endl;
        std::cout << "  Cell Size: " << space.cell_size << std::endl;
        std::cout << "  Global Low Corner:" << std::endl;
        std::cout << "    X: " << space.global_low_corner[0] << std::endl;
        std::cout << "    Y: " << space.global_low_corner[1] << std::endl;
        std::cout << "    Z: " << space.global_low_corner[2] << std::endl;
        std::cout << "  Global High Corner:" << std::endl;
        std::cout << "    X: " << space.global_high_corner[0] << std::endl;
        std::cout << "    Y: " << space.global_high_corner[1] << std::endl;
        std::cout << "    Z: " << space.global_high_corner[2] << std::endl;

        // Print properties
        std::cout << "Properties:" << std::endl;
        std::cout << "  Density: " << properties.density << std::endl;
        std::cout << "  Specific Heat: " << properties.specific_heat << std::endl;
        std::cout << "  Thermal Conductivity: " << properties.thermal_conductivity << std::endl;
        std::cout << "  Latent Heat: " << properties.latent_heat << std::endl;

        // Print source
        std::cout << "Source:" << std::endl;
        std::cout << "  Absorption: " << source.absorption << std::endl;
        std::cout << "  two-sigma:" << std::endl;
        std::cout << "    X: " << source.two_sigma[0] << std::endl;
        std::cout << "    Y: " << source.two_sigma[1] << std::endl;
        std::cout << "    Z: " << source.two_sigma[2] << std::endl;
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

        if (!filename)
        {
            std::cerr << "Input file required." << std::endl
                      << "Usage: " << argv[0]
                      << " -i <input_yaml_file>" << std::endl;
            return;
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

            // Read heat source components
            source.absorption =
                db["source"]["absorption"].as<double>();
            source.two_sigma =
                db["source"]["two_sigma"].as<std::array<double, 3>>();

            source.two_sigma[0] = fabs( source.two_sigma[0] );
            source.two_sigma[1] = fabs( source.two_sigma[1] );
            source.two_sigma[2] = fabs( source.two_sigma[2] );
        }
    }
};
