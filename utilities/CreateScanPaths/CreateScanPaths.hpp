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

struct Point
{
    double x;
    double y;

    // Default constructor
    Point()
        : x( 0.0 )
        , y( 0.0 )
    {
    }

    // Construct from x and y
    Point( double xCoord, double yCoord )
        : x( xCoord )
        , y( yCoord )
    {
    }

    // Function to rotate a point around a specified origin by a given angle
    Point rotate( const Point& origin, double degrees ) const
    {
        double angle = degrees * ( M_PI / 180.0 );

        double s = sin( angle );
        double c = cos( angle );
        double translatedX = x - origin.x;
        double translatedY = y - origin.y;
        double newX = translatedX * c - translatedY * s + origin.x;
        double newY = translatedX * s + translatedY * c + origin.y;
        return Point( newX, newY );
    }

    // Define a custom operator<< to print a Point
    friend std::ostream& operator<<( std::ostream& os, const Point& point )
    {
        os << "(" << point.x << ", " << point.y << ")";
        return os;
    }
};

// Function to calculate the distance between two points
double distance( const Point& p1, const Point& p2 )
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt( dx * dx + dy * dy );
}

struct Line
{
    Point start;
    Point end;

    // Default constructor
    Line()
        : start( Point() )
        , end( Point() )
    {
    }

    // Construct from two points
    Line( Point startPoint, Point endPoint )
        : start( startPoint )
        , end( endPoint )
    {
    }

    // Function to rotate the line around a specified origin by a given angle
    void rotate( const Point& origin, double angle )
    {
        start = start.rotate( origin, angle );
        end = end.rotate( origin, angle );
    }

    bool isFinite()
    {
        return ( !std::isnan( start.x ) && !std::isnan( start.y ) &&
                 !std::isnan( end.x ) && !std::isnan( end.y ) );
    }

    // Define a custom operator<< to print a Line
    friend std::ostream& operator<<( std::ostream& os, const Line& line )
    {
        os << "(" << line.start.x << ", " << line.start.y << "), "
           << "(" << line.end.x << ", " << line.end.y << ")";
        return os;
    }
};

struct boundBox
{
    Point minPoint;
    Point maxPoint;
    Point midPoint;

    std::vector<Line> edges;

    // Construct from two points
    boundBox( Point minP, Point maxP )
        : minPoint( minP )
        , maxPoint( maxP )
        , midPoint( ( minP.x + maxP.x ) / 2.0, ( minP.y + maxP.y ) / 2.0 )
    {
        // left, right, tiop, bottom
        edges = { { Point( minPoint.x, minPoint.y ),
                    Point( minPoint.x, maxPoint.y ) },
                  { Point( maxPoint.x, minPoint.y ),
                    Point( maxPoint.x, maxPoint.y ) },
                  { Point( minPoint.x, maxPoint.y ),
                    Point( maxPoint.x, maxPoint.y ) },
                  { Point( minPoint.x, minPoint.y ),
                    Point( maxPoint.x, minPoint.y ) } };
    }

    // Function to check if a point is inside the bounding box
    bool isInside( Point p )
    {
        return p.x >= minPoint.x && p.x <= maxPoint.x && p.y >= minPoint.y &&
               p.y <= maxPoint.y;
    }

    Line cropLine( const Line& line )
    {
        std::vector<Point> intersections;

        for ( const Line& edge : edges )
        {
            Point intersection = intersect( edge, line );

            if ( !std::isnan( intersection.x ) &&
                 !std::isnan( intersection.y ) )
            {
                intersections.push_back( intersection );
            }
        }

        if ( intersections.empty() )
        {
            return Line( Point( NAN, NAN ), Point( NAN, NAN ) );
        }

        // Sort intersection points based on position along the original line
        std::sort( intersections.begin(), intersections.end(),
                   [&line]( const Point& p1, const Point& p2 )
                   {
                       double d1 = distance( line.start, p1 );
                       double d2 = distance( line.start, p2 );
                       return d1 < d2;
                   } );

        Line intersectedLine( intersections.front(), intersections.back() );

        return intersectedLine;
    }

    // Function to find the intersection point between two lines
    Point intersect( const Line& line1, const Line& line2 )
    {
        Point intersection( std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::quiet_NaN() );

        // components of first line
        double x1 = line1.start.x;
        double y1 = line1.start.y;
        double x2 = line1.end.x;
        double y2 = line1.end.y;

        // components of second line
        double x3 = line2.start.x;
        double y3 = line2.start.y;
        double x4 = line2.end.x;
        double y4 = line2.end.y;

        double denominator =
            ( x1 - x2 ) * ( y3 - y4 ) - ( y1 - y2 ) * ( x3 - x4 );

        // Lines are parallel or colinear, no intersection
        if ( denominator == 0 )
        {
            return intersection;
        }

        double t = ( ( x1 - x3 ) * ( y3 - y4 ) - ( y1 - y3 ) * ( x3 - x4 ) ) /
                   denominator;

        double u = -( ( x1 - x2 ) * ( y1 - y3 ) - ( y1 - y2 ) * ( x1 - x3 ) ) /
                   denominator;

        if ( t >= 0 && t <= 1 && u >= 0 && u <= 1 )
        {
            intersection.x = x1 + t * ( x2 - x1 );
            intersection.y = y1 + t * ( y2 - y1 );
            return intersection;
        }

        return intersection;
    }
};

struct Path
{
    std::vector<Line> lines;
    double power;
    double speed;
    double dwell_time;

    // Construct the path from a bounding box and hatch spacing
    Path( boundBox bbox, double step, double angle )
    {
        int numLines = numberOfLines( bbox, step );

        // create a pad of infinitely long, equally parallel lines
        std::vector<Line> pathLines;

        const float great = 1e10;

        // Create lines in the negative direction, excluding the midpoint line
        for ( int i = numLines - 1; i > 0; --i )
        {
            double height = bbox.midPoint.y - i * step;
            Line currentLine( Point( -great, height ), Point( great, height ) );
            pathLines.push_back( currentLine );
        }

        // Create lines in the positive direction, including the midpoint
        for ( int i = 0; i < numLines; ++i )
        {
            double height = bbox.midPoint.y + i * step;
            Line currentLine( Point( -great, height ), Point( great, height ) );
            pathLines.push_back( currentLine );
        }

        // apply rotation and cropping to the initial lines
        for ( const Line& pathLine : pathLines )
        {
            // Rotate the endpoints by the specified angle
            Line rotatedLine = pathLine;
            rotatedLine.rotate( bbox.midPoint, angle );

            Line croppedLine = bbox.cropLine( rotatedLine );

            if ( croppedLine.isFinite() )
            {
                lines.push_back( croppedLine );
            }
        }
    }

    // Function to find the number of scan vectors in the bounding box
    int numberOfLines( const boundBox& bbox, double step ) const
    {
        double rangeX = bbox.maxPoint.x - bbox.minPoint.x;
        double rangeY = bbox.maxPoint.y - bbox.minPoint.y;

        int nX = 0;
        int nY = 0;

        for ( double x = bbox.minPoint.x; x <= bbox.maxPoint.x; x += step )
        {
            nX++;
        }

        for ( double y = bbox.minPoint.y; y <= bbox.maxPoint.y; y += step )
        {
            nY++;
        }

        return ( nX > nY ) ? nX : nY;
    }

    void write( const std::string& filename,
                const bool bi_direction = true ) const
    {
        std::ofstream file( filename );

        file << "Mode\tX(m)\tY(m)\tZ(m)\tPower(W)\ttParam" << std::endl;

        for ( size_t i = 0; i < lines.size(); ++i )
        {
            const Line& line = lines[i];

            Point first = line.start;
            Point second = line.end;

            // reverse the odd lines for bi_directional
            if ( bi_direction && i % 2 == 1 )
            {
                first = line.end;
                second = line.start;
            }

            // hatch (with skywrite)
            if ( i == 0 )
            {
                // no initial dwell
                file << "1\t" << first.x << "\t" << first.y << "\t0\t0\t0\n";
            }
            else
            {
                file << "1\t" << first.x << "\t" << first.y << "\t0\t" << 0
                     << "\t" << dwell_time << "\n";
            }

            // raster
            file << "0\t" << second.x << "\t" << second.y << "\t0\t" << power
                 << "\t" << speed << "\n";
        }

        file.close();
    }
};
