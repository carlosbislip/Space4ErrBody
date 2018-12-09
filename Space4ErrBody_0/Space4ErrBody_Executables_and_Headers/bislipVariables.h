#ifndef BISLIP_VARIABLES_H
#define BISLIP_VARIABLES_H

#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <utility>

#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>

namespace bislip {

namespace variables {

enum GuidanceParameter
{
    AngleOfAttack,
    BankAngle,
    ThrustElevationAngle,
    ThrustAzimuthAngle,
    ThrottleSetting
};

std::string getCurrentDateTime ( const bool useLocalTime = false );

std::vector< std::string > getDataString ( const std::string &filename );

std::vector< double > getDataNumeri ( const std::string &filename );

double computeAngularDistance (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f);

double computeHeadingToTarget (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f);

double computeHeadingError (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f,
        const double &heading);

double computeSpecificEnergy (
        const double &height,
        const double &airspeed);

double computeNormalizedSpecificEnergy(
        const double &height,
        const double &airspeed,
        const double &E_max);

std::string passGuidanceParameter (
        const std::string &parameter);

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator(
        const double &flight_path_angle,
        const std::string &parameter,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems);

std::pair < double, double > chooseGuidanceBounds (
        const double &flight_path_angle,
        const std::string &parameter,
        const std::shared_ptr<tudat::system_models::VehicleSystems> &vehicleSystems);

double evaluateGuidanceInterpolator (
        const double &flightpathangle,
        const std::string &parameter,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const double &height,
        const double &airspeed,
        const double &E_max);

Eigen::Vector3d computeBodyFixedThrustDirection (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems);

double computeThrustMagnitude (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems);

Eigen::Vector3d computeBodyFixedThrustVector (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems);

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass);

double computeHeatingRate (
        const double &airdensity,
        const double &airspeed,
        const double &C,
        const double &N,
        const double &M);

double computeStagnationHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems);

double computeFlatPlateHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems);

double computeStagnationHeat (
        const double &airdensity,
        const double &airspeed,
        const double &C_s,
        const double &N,
        const double &M,
        const double &adiabaticWallTemperature,
        const double &WallTemperature);

double computeFlatPlateHeat (
        const double &airdensity,
        const double &airspeed,
        const double &C_FP_1,
        const double &C_FP_2,
        const double &adiabaticWallTemperature,
        const double &WallTemperature);

double computeHeatingRateTauber (
        const double &q_dot_s,
        const double &q_dot_FP,
        const double &lambda);

} // namespace variables
} // namespace bislip

#endif // BISLIP_VARIABLES_H
