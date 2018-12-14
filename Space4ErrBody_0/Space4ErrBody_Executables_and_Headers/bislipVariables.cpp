#include <iostream>
#include <ctime>
#include <cmath>
#include <sstream>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/Aerodynamics/equilibriumWallTemperature.h>
#include <Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h>



#include "updateGuidance.h"

namespace bislip {

namespace variables {

std::string getCurrentDateTime( bool useLocalTime ) {
    std::stringstream currentDateTime;

    // current date/time based on current system
    time_t ttNow = time(nullptr);
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
    //const double angular_distance = std::acos( std::sin( lat_c ) * std::sin( lat_f ) + std::cos( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

    return std::acos( std::sin( lat_c ) * std::sin( lat_f ) + std::cos( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );
}

double computeHeadingToTarget (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    // const double heading_to_target = std::atan2( std::sin( lon_f - lon_c ) * std::cos( lat_f ) , std::cos( lat_c ) * std::sin( lat_f ) - std::sin( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

    return std::atan2( std::sin( lon_f - lon_c ) * std::cos( lat_f ) , std::cos( lat_c ) * std::sin( lat_f ) - std::sin( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

}

double computeHeadingError (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f,
        const double &heading)
{
    // const double heading_to_target = computeHeadingToTarget( lat_c, lon_c, lat_f, lon_f );

    return computeHeadingToTarget( lat_c, lon_c, lat_f, lon_f ) - heading;
}

double computeSpecificEnergy (
        const double &height,
        const double &airspeed)
{
    return tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * height + 0.5 * airspeed * airspeed;
}

double computeNormalizedSpecificEnergy (
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    return computeSpecificEnergy( height, airspeed ) / E_max;
}

//! https://doi.org/10.1016/j.cam.2017.09.049
std::vector< double > HermiteDerivatives( const Eigen::VectorXd &E_mapped, const Eigen::VectorXd &y, const long &nodes )
{
    Eigen::VectorXd h( nodes - 1 ), dely( nodes - 1 ), b( nodes ), x( nodes );
    Eigen::MatrixXd A( nodes, nodes );
    std::vector< double > x_vect;
    double mu, lambda;

    for ( long i = 0; i < nodes - 1; ++i ) { h( i ) = E_mapped( i + 1 ) - E_mapped( i ); dely( i ) = ( y( i + 1 ) - y( i ) ) / h( i ); }

    A = Eigen::MatrixXd::Zero( nodes, nodes );
    A( 0, 0 ) = 4.0;
    A( 0, 1 ) = -1.0;
    A( nodes - 1, nodes - 2 ) = -1.0;
    A( nodes - 1, nodes - 1 ) = 4.0;

    for ( long i = 1; i < nodes - 1; ++i )
    {
        lambda = h( i ) / ( h( i - 1 ) + h( i ) );
        mu = 1 - lambda;
        A( i , i - 1 ) = -mu;
        A( i , i ) = 4.0;
        A( i , i + 1 ) = -lambda * mu;
    }

    b( 0 ) = dely( 0 );
    b( nodes - 1 ) = dely( nodes - 2 );

    for ( long i = 1; i < nodes - 1; ++i ) { b( i ) = 3 * ( y( i + 1 ) - y( i - 1 ) ) / ( h( i - 1 ) + h( i ) ); }

    x = A.fullPivHouseholderQr().solve( b );

    for ( long i = 0; i < nodes; ++i ) { x_vect.push_back( x( i ) ); }

    return x_vect;
}

std::string passGuidanceParameter (
        const std::string &parameter)
{
    return parameter;
}

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator(
        const double &flight_path_angle,
        const std::string &parameter,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator;

    if ( flight_path_angle >  0 )
    {
        if (parameter == "Angle of Attack")
        {
            interpolator = vehicleSystems->getThrustElevationAngleInterpolator_Ascent();
        }
        if (parameter == "Bank Angle")
        {
            interpolator = vehicleSystems->getBankAngleInterpolator_Ascent();
        }
        if (parameter == "Thrust Elevation Angle")
        {
            interpolator = vehicleSystems->getThrustElevationAngleInterpolator_Ascent();
        }
        if (parameter == "Thrust Azimuth Angle")
        {
            interpolator = vehicleSystems->getThrustAzimuthAngleInterpolator_Ascent();
        }
        if (parameter == "Throttle Setting")
        {
            interpolator = vehicleSystems->getThrottleInterpolator_Ascent();
        }
    }
    else
    {
        if (parameter == "Angle of Attack")
        {
            interpolator = vehicleSystems->getThrustElevationAngleInterpolator_Descent();
        }
        if (parameter == "Bank Angle")
        {
            interpolator = vehicleSystems->getBankAngleInterpolator_Descent();
        }
        if (parameter == "Thrust Elevation Angle")
        {
            interpolator = vehicleSystems->getThrustElevationAngleInterpolator_Descent();
        }
        if (parameter == "Thrust Azimuth Angle")
        {
            interpolator = vehicleSystems->getThrustAzimuthAngleInterpolator_Descent();
        }
        if (parameter == "Throttle Setting")
        {
            interpolator = vehicleSystems->getThrottleInterpolator_Descent();
        }
    }
    return interpolator;
}

std::pair < double, double > chooseGuidanceBounds (
        const double &flight_path_angle,
        const std::string &parameter,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    std::pair < double, double > bounds;

    if ( flight_path_angle >  0 )
    {
        if (parameter == "Angle of Attack")
        {
            bounds = vehicleSystems->getParameterBounds( "ascent", parameter );
        }
        if (parameter == "Bank Angle")
        {
            bounds = vehicleSystems->getParameterBounds( "ascent", parameter );
        }
        if (parameter == "Thrust Elevation Angle")
        {
            bounds = vehicleSystems->getParameterBounds( "ascent", parameter );
        }
        if (parameter == "Thrust Azimuth Angle")
        {
            bounds = vehicleSystems->getParameterBounds( "ascent", parameter );
        }
        if (parameter == "Throttle Setting")
        {
            bounds = vehicleSystems->getParameterBounds( "ascent", parameter );
        }
    }
    else
    {
        if (parameter == "Angle of Attack")
        {
            bounds = vehicleSystems->getParameterBounds( "descent", parameter );
        }
        if (parameter == "Bank Angle")
        {
            bounds = vehicleSystems->getParameterBounds( "descent", parameter );
        }
        if (parameter == "Thrust Elevation Angle")
        {
            bounds = vehicleSystems->getParameterBounds( "descent", parameter );
        }
        if (parameter == "Thrust Azimuth Angle")
        {
            bounds = vehicleSystems->getParameterBounds( "descent", parameter );
        }
        if (parameter == "Throttle Setting")
        {
            bounds = vehicleSystems->getParameterBounds( "descent", parameter );
        }
    }
    return bounds;
}

double evaluateGuidanceInterpolator (
        const double &flightpathangle,
        const std::string &parameter,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    //! Select parameter interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator = bislip::variables::chooseGuidanceInterpolator( flightpathangle, parameter, vehicleSystems );

    //! Select parameter bounds.
    std::pair < double, double > bounds =  bislip::variables::chooseGuidanceBounds( flightpathangle, parameter, vehicleSystems );

    //! Evaluate interpolator.
    double evaluation = interpolator->interpolate( bislip::variables::computeNormalizedSpecificEnergy( height, airspeed, E_max ) );

    //! Impose bounds.
    if ( evaluation < bounds.first ){ evaluation = bounds.first; }
    if ( evaluation > bounds.second ){ evaluation = bounds.second; }

    return evaluation;
}

Eigen::Vector3d computeBodyFixedThrustDirection (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    // std::cout << "Computing Body Fixed Thrust Direction" << std::endl;
    const double eps_T =  evaluateGuidanceInterpolator (
                flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle ),
                "Thrust Elevation Angle",
                vehicleSystems,
                flightConditions->getCurrentAltitude(),
                flightConditions->getCurrentAirspeed(),
                vehicleSystems->getE_max() );

    const double phi_T = evaluateGuidanceInterpolator (
                flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle ),
                "Thrust Azimuth Angle",
                vehicleSystems,
                flightConditions->getCurrentAltitude(),
                flightConditions->getCurrentAirspeed(),
                vehicleSystems->getE_max() );

    Eigen::Vector3d bodyFixedThrustDirection;
    bodyFixedThrustDirection( 0 ) = std::cos( tudat::unit_conversions::convertDegreesToRadians( eps_T ) ) * std::cos( tudat::unit_conversions::convertDegreesToRadians( phi_T ) );
    bodyFixedThrustDirection( 1 ) = std::cos( tudat::unit_conversions::convertDegreesToRadians( eps_T ) ) * std::sin( tudat::unit_conversions::convertDegreesToRadians( phi_T ) );
    bodyFixedThrustDirection( 2 ) = std::sin( tudat::unit_conversions::convertDegreesToRadians( eps_T ) );

    return bodyFixedThrustDirection;
}

double computeThrustMagnitude (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    const double throttle = bislip::variables::evaluateGuidanceInterpolator (
                flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle ),
                "Throttle Setting",
                vehicleSystems,
                flightConditions->getCurrentAltitude(),
                flightConditions->getCurrentAirspeed(),
                vehicleSystems->getE_max() );

    return throttle * vehicleSystems->getMaxThrust();
}

Eigen::Vector3d computeBodyFixedThrustVector (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    const Eigen::Vector3d BodyFixedThrustDirection = bislip::variables::computeBodyFixedThrustDirection( flightConditions, vehicleSystems );
    const double thrustMagnitude = bislip::variables::computeThrustMagnitude( flightConditions, vehicleSystems );

    Eigen::Vector3d BodyFixedThrustVector;

    BodyFixedThrustVector( 0 ) = thrustMagnitude * BodyFixedThrustDirection( 0 );
    BodyFixedThrustVector( 1 ) = thrustMagnitude * BodyFixedThrustDirection( 1 );
    BodyFixedThrustVector( 2 ) = thrustMagnitude * BodyFixedThrustDirection( 2 );

    return BodyFixedThrustVector;
}

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass)
{
    bool currentEngineStatus = true;
    if ( currentMass <= landingMass ){ currentEngineStatus = false; }
    return currentEngineStatus;
}

double computeHeatingRate (
        const double &airdensity,
        const double &airspeed,
        const double &C,
        const double &N,
        const double &M)
{
    return C * std::pow( airdensity, N ) * std::pow( airspeed, M );
}

double computeStagnationHeat (const double &airdensity,
                              const double &airspeed,
                              const double &C_s,
                              const double &N,
                              const double &M,
                              const double &adiabaticWallTemperature,
                              const double &WallTemperature)
{
    double C = C_s * ( 1 - ( WallTemperature / adiabaticWallTemperature ) ) ;
    double q_dot_s = computeHeatingRate ( airdensity, airspeed, C, N, M) ;
    //std::cout << "q_dot_s: " << q_dot_s << std::endl;

    return q_dot_s;
}

double computeStagnationHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    const double airDensity = flightConditions->getCurrentDensity( );
    const double airSpeed = flightConditions->getCurrentAirspeed( );
    const double airTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double machNumber = flightConditions->getCurrentMachNumber( );
    const double noseRadius = vehicleSystems->getNoseRadius( );
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );

    const double nose_term = std::pow( noseRadius, -0.5 );
    const double C_s = ( 1.83 ) * nose_term * 1E-4;
    const double M = 3.0;
    const double N = 0.5;

    // Compute adiabatic wall temperature.
    double adiabaticWallTemperature
            = tudat::aerodynamics::computeAdiabaticWallTemperature( airTemperature , machNumber );

    std::function< double( const double ) > heatTransferFunction = std::bind(
                &computeStagnationHeat, airDensity, airSpeed, C_s, N, M, adiabaticWallTemperature, std::placeholders::_1 );

    return tudat::aerodynamics::computeEquilibriumHeatflux( heatTransferFunction, wallEmissivity, adiabaticWallTemperature );
}

double computeFlatPlateHeat (const double &airdensity,
                             const double &airspeed,
                             const double &C_FP_1,
                             const double &C_FP_2,
                             const double &adiabaticWallTemperature,
                             const double &WallTemperature)
{
    double C;
    double M;
    double N;

    if ( airspeed <= 3962.0 )
    {
        M = 3.37;
        N = 0.8;
        C = C_FP_1 * std::pow( 556 / WallTemperature, 1.0 / 4.0 ) * ( 1.0 - 1.11 * ( WallTemperature / adiabaticWallTemperature ) );
    }
    if ( airspeed > 3962.0 )
    {
        M = 3.7;
        N = 0.8;
        C = C_FP_2 * ( 1 - 1.11 * ( WallTemperature / adiabaticWallTemperature ) );
    }

    return computeHeatingRate ( airdensity, airspeed, C, N, M);;
}

double computeFlatPlateHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    const double airDensity = flightConditions->getCurrentDensity( );
    const double airSpeed = flightConditions->getCurrentAirspeed( );
    const double airTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double machNumber = flightConditions->getCurrentMachNumber( );
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );
    const double phi = vehicleSystems->getLocalBodyAngle( );
    const double x_T = vehicleSystems->getTransitionDistance( );

    const double sin_phi_term = std::pow( std::sin( phi ), 1.6 ) ;
    const double cos_phi = std::cos( phi );
    const double x_T_term = std::pow( x_T, -1.0 / 5.0 ) ;
    const double C_FP_1 = ( 3.35 ) * sin_phi_term * std::pow( cos_phi, 1.78 ) * x_T_term * 1E-4;
    const double C_FP_2 = ( 2.2 ) * sin_phi_term * std::pow( cos_phi, 2.08 ) * x_T_term * 1E-5;

    // Compute adiabatic wall temperature.
    double adiabaticWallTemperature
            = tudat::aerodynamics::computeAdiabaticWallTemperature( airTemperature , machNumber );

    std::function< double( const double ) > heatTransferFunction = std::bind(
                &computeFlatPlateHeat, airDensity, airSpeed, C_FP_1, C_FP_2, adiabaticWallTemperature, std::placeholders::_1 );

    return tudat::aerodynamics::computeEquilibriumHeatflux( heatTransferFunction, wallEmissivity, adiabaticWallTemperature );
}

double computeHeatingRateTauber (
        const double &q_dot_s,
        const double &q_dot_FP,
        const double &lambda)
{
    const double sin_lambda = std::sin( lambda );
    const double cos_lambda = std::cos( lambda );

    double q_dot_LE = std::pow ( 0.5 * q_dot_s * q_dot_s * cos_lambda * cos_lambda + q_dot_FP * q_dot_FP * sin_lambda * sin_lambda, 0.5 );
    //  std::cout << "q_dot_FP: " << q_dot_FP << std::endl;
    // std::cout << "q_dot_LE: " << q_dot_LE << std::endl;

    return q_dot_LE;
}
/*

Eigen::MatrixXd getDependentVariableMatrix( const tudat::propagators::SingleArcDynamicsSimulator< double > simulatedDynamics, const double simulationStartEpoch, const double fixedStepSize )
{

    const std::shared_ptr< tudat::propagators::SingleArcDynamicsSimulator< double, double > >& singleArcDynamicsSimulator,


    //! Extract map of dependent variables.
    const std::map< double, Eigen::VectorXd > dependentVariableMap = simulatedDynamics.SingleArcDynamicsSimulator::getDependentVariableHistory( );

    //! Declare number of rows (time) / columns (variables).
    unsigned long rows = dependentVariableMap.size();
    unsigned long columns = ( ( simulatedDynamics.getDependentVariableHistory( ).begin() )->second ).size();

    //! Declare and initialize dependent variable matrix.
    Eigen::MatrixXd dependentVariableMatrix( rows, columns );
    dependentVariableMatrix = Eigen::MatrixXd::Zero( rows, columns );

    //! Loop to populate the matrix with the extraction of the map.
    for ( unsigned long i = 0; i < rows; i++ )
    {
        dependentVariableMatrix.row( i ) = dependentVariableMap.at( simulationStartEpoch + i * fixedStepSize );

        std::cout << "dependentVariableMatrix.row( " << i << " ) : " << dependentVariableMatrix.row( i ) << std::endl;

    }

    return dependentVariableMatrix;
}

*/


} // namespace variables
} // namespace bislip
