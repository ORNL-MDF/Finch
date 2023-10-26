#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <unistd.h>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "CreateScanPaths.hpp"

int main( int argc, char* argv[] )
{
    // Read input file
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
            return 0;
        }
    }

    YAML::Node config = YAML::LoadFile( filename );

    Point minPoint;
    minPoint.x = config["min_point"][0].as<double>();
    minPoint.y = config["min_point"][1].as<double>();

    Point maxPoint;
    maxPoint.x = config["max_point"][0].as<double>();
    maxPoint.y = config["max_point"][1].as<double>();

    double angle = config["angle"].as<double>();
    double hatch = config["hatch"].as<double>();
    int num_rotations = config["num_rotations"].as<int>();

    double power = config["power"].as<double>();
    double speed = config["speed"].as<double>();
    double dwell_time = config["dwell_time"].as<double>();

    bool bi_direction = true;
    if ( config["bi_direction"] )
    {
        bi_direction = config["bi_direction"].as<bool>();
    }

    // Create bounding box for scan vectors
    boundBox boundingBox( minPoint, maxPoint );

    // Rotate and write scan vectors to file
    double rotation = 0.0;

    for ( int i = 0; i < num_rotations; ++i )
    {
        // create new path
        Path path( boundingBox, hatch, rotation );
        path.power = power;
        path.speed = speed;
        path.dwell_time = dwell_time;

        // write new path
        std::ostringstream oss;
        oss << std::fixed << std::setprecision( 0 ) << rotation;
        std::string rotationString = oss.str();
        std::string filename = "path_" + rotationString + ".txt";

        path.write( filename, bi_direction );

        rotation += angle;
    }

    return 0;
}
