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
#include <stdexcept>
#include <unistd.h>
#include <vector>

#include <nlohmann/json.hpp>

#include "Finch_CreateScanPaths.hpp"

int main( int argc, char* argv[] )
{
    try
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
                return 1;
            }
        }

        if ( filename == nullptr )
            throw std::runtime_error( "Usage: " + std::string( argv[0] ) +
                                      " -i <input_json_file>" );

        // parse input file
        std::ifstream config_stream( filename );
        if ( !config_stream )
            throw std::runtime_error( "Cannot open input file " +
                                      std::string( filename ) );
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

        bool bi_direction = config.value( "bi_direction", true );

        if ( !std::isfinite( minPoint.x ) || !std::isfinite( minPoint.y ) ||
             !std::isfinite( maxPoint.x ) || !std::isfinite( maxPoint.y ) ||
             maxPoint.x <= minPoint.x || maxPoint.y <= minPoint.y ||
             !std::isfinite( hatch ) || hatch <= 0.0 ||
             !std::isfinite( angle ) || num_rotations <= 0 ||
             !std::isfinite( power ) || power < 0.0 ||
             !std::isfinite( speed ) || speed <= 0.0 ||
             !std::isfinite( dwell_time ) || dwell_time < 0.0 )
            throw std::runtime_error(
                "Invalid scan-path inputs: bounds must have positive area, "
                "hatch, speed, and rotation count must be positive, and power "
                "and dwell time must be nonnegative" );

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
            std::string output_filename = "path_" + rotationString + ".txt";

            path.write( output_filename, bi_direction );

            rotation += angle;
        }

        return 0;
    }
    catch ( const std::exception& e )
    {
        std::cerr << "create_scan_paths: " << e.what() << std::endl;
        return 1;
    }
}
