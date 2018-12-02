#include <iostream>
#include <ctime>
#include <sstream>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

namespace bislip {

namespace variables {

std::string getCurrentDateTime( bool useLocalTime ) {
    std::stringstream currentDateTime;

    // current date/time based on current system
    time_t ttNow = time(0);
    tm * ptmNow;

    if (useLocalTime)
        ptmNow = localtime(&ttNow);
    else
        ptmNow = gmtime(&ttNow);

    currentDateTime << 1900 + ptmNow->tm_year;

    //month
    if (ptmNow->tm_mon < 9)
        //Fill in the leading 0 if less than 10
        currentDateTime << "0" << 1 + ptmNow->tm_mon;
    else
        currentDateTime << (1 + ptmNow->tm_mon);

    //day
    if (ptmNow->tm_mday < 10)
        currentDateTime << "0" << ptmNow->tm_mday << "_";
    else
        currentDateTime <<  ptmNow->tm_mday << "_";

    //hour
    if (ptmNow->tm_hour < 10)
        currentDateTime << "0" << ptmNow->tm_hour;
    else
        currentDateTime << ptmNow->tm_hour;

    //min
    if (ptmNow->tm_min < 10)
        currentDateTime << "0" << ptmNow->tm_min;
    else
        currentDateTime << ptmNow->tm_min;

    //sec
    if (ptmNow->tm_sec < 10)
        currentDateTime << "0" << ptmNow->tm_sec;
    else
        currentDateTime << ptmNow->tm_sec;

    //std::cout << "Got time: " << currentDateTime.str() << std::endl;

    return currentDateTime.str();
}

std::vector< std::string > getDataString ( const std::string &filename )
{
    std::ifstream inputdata;
    inputdata.open( filename.c_str( ) );
    if( inputdata.fail( ) )
    {
        std::cout << "File '" << filename << "' failed to open";
        std::cin.ignore( 100, '\n' );
        std::exit( 1 );
    }
    else{

        std::string  var;
        std::vector< std::string > stuff;

        while( std::getline( inputdata, var ) )
        {
            if( var.size( ) > 0 )
            {
                stuff.push_back( var );
            }
        }
        inputdata.close( );
        return stuff;
    }
}

std::vector< double > getDataNumeri ( const std::string &filename )
{
    std::ifstream inputdata;
    inputdata.open( filename.c_str( ) );
    if( inputdata.fail( ) )
    {
        std::cout << "File '" << filename << "' failed to open";
        std::cin.ignore( 100, '\n' );
        std::exit( 1 );
    }
    else{

        double var;
        std::vector< double > stuff;

        while (!inputdata.fail() && !inputdata.eof())
        {
            inputdata >> var;
            stuff.push_back(var);
        }

        inputdata.close( );
        return stuff;
    }
}

double computeAngularDistance (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    const double angular_distance = std::acos( std::sin( lat_c ) * std::sin( lat_f ) + std::cos( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

    return angular_distance;
}

double computeHeadingToTarget (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    const double heading_to_target = std::atan2( std::sin( lon_f - lon_c ) * std::cos( lat_f ) , std::cos( lat_c ) * std::sin( lat_f ) - std::sin( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

    return heading_to_target;
}

double computeHeadingError (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f,
        const double &heading)
{
    const double heading_to_target = computeHeadingToTarget( lat_c, lon_c, lat_f, lon_f );

    return heading_to_target -  heading;
}

double computeSpecificEnergy (
        const double &height,
        const double &airspeed)
{
    return 9.80665 * height + 0.5 * airspeed * airspeed;
}

double computeNormalizedSpecificEnergy (
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    return computeSpecificEnergy( height, airspeed ) / E_max;
}

double computeAngleofAttack (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_AoA,
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    return interpolator_AoA->interpolate( computeNormalizedSpecificEnergy( height, airspeed, E_max ) );
}

double computeThrottleSetting (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_throttle,
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    return interpolator_throttle->interpolate( computeNormalizedSpecificEnergy( height, airspeed, E_max ) );
}

double computeThrustElevationAngle (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_eps_T,
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    return interpolator_eps_T->interpolate( computeNormalizedSpecificEnergy( height, airspeed, E_max )  );
}

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass)
{
    bool currentEngineStatus = true;

    if ( currentMass <= landingMass )
    {
        currentEngineStatus = false;
    }

    return currentEngineStatus;
}
} // namespace variables
} // namespace bislip
