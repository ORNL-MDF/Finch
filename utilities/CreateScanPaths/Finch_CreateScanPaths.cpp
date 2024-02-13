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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <unistd.h>
#include <vector>

#include <nlohmann/json.hpp>

#include "Finch_CreateScanPaths.hpp"

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
            std::cerr << "Usage: " << argv[0] << " -i <input_json_file>"
                      << std::endl;
            return 0;
        }
    }

    // parse input file
    std::ifstream config_stream( filename );
    nlohmann::json config = nlohmann::json::parse( config_stream );

    Finch::Point minPoint;
    minPoint.x = config["min_point"][0];
    minPoint.y = config["min_point"][1];

    Finch::Point maxPoint;
    maxPoint.x = config["max_point"][0];
    maxPoint.y = config["max_point"][1];

    double angle = config["angle"];
    double hatch = config["hatch"];
    int num_rotations = config["num_rotations"];

    double power = config["power"];
    double speed = config["speed"];
    double dwell_time = config["dwell_time"];

    bool bi_direction = true;
    if ( config["bi_direction"] )
    {
        bi_direction = config["bi_direction"];
    }

    // Create bounding box for scan vectors
    Finch::boundBox boundingBox( minPoint, maxPoint );

    // Rotate and write scan vectors to file
    double rotation = 0.0;

    for ( int i = 0; i < num_rotations; ++i )
    {
        // create new path
        Finch::Path path( boundingBox, hatch, rotation );
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
