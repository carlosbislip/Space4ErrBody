#ifndef BISLIP_VARIABLES_H
#define BISLIP_VARIABLES_H

#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <utility>

#include <Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h>
#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>
#include <Tudat/Astrodynamics/SystemModels/bislipSystems.h>
#include <Tudat/Astrodynamics/SystemModels/bislipOptimizationParameters.h>

namespace bislip {

namespace variables {

/*
enum OptimizationParameter
{
    AngleOfAttack,
    BankAngle,
    ThrustElevationAngle,
    ThrustAzimuthAngle,
    ThrottleSetting
};
*/


std::string getCurrentDateTime ( const bool useLocalTime = false );

unsigned int millis_since_midnight ( );

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

double computeNormalizedSpecificEnergy (
        const double &height,
        const double &airspeed,
        const double &E_max);

std::vector< double > HermiteDerivatives (
        const Eigen::VectorXd &mappedNormalizedSpecificEnergy,
        const Eigen::VectorXd &y);

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > createOneDimensionalHermiteInterpolator (
        const Eigen::VectorXd &parameterValues,
        const Eigen::VectorXd &normalizedSpecificEnergy,
        const std::map< double, double > &mapped_data,
        const std::shared_ptr< tudat::interpolators::InterpolatorSettings > &interpolatorSettings );

//std::string passOptimizationParameter (
//        const std::string &parameter);

bislip::variables::OptimizationParameter passOptimizationParameter (
        const bislip::variables::OptimizationParameter &parameter);

//std::string passDirection (
  //      const std::string &direction);

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator (
        const bislip::variables::OptimizationParameter &parameter,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems);

std::pair < double, double > chooseGuidanceBounds (
        const bislip::variables::OptimizationParameter &parameter,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems);

double evaluateGuidanceInterpolator (
        const bislip::variables::OptimizationParameter &parameter,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems,
        const double &height,
        const double &airspeed,
        const double &E_max);

Eigen::Vector6d computeCurrentCoefficients (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems);

Eigen::Vector3d computeBodyFixedThrustDirection (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems);

double computeThrustMagnitude (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems);

Eigen::Vector3d computeBodyFixedThrustVector (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems);

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass);

double getBodyFlapDeflection( const double &del_C_m_b );

Eigen::Vector2d getGravs (
        const double &r,
        const double &latitude);

double computeEquilibriumGlideLimit (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems,
        const double &currentMass );

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
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems);

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

double computeBendingMoment (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions);

double computePenalty (
        const Eigen::VectorXd &dependentVariable_TimeHistory,
        const long &startIterator,
        const long &endIterator,
        const double &constraint, const double &fixedStepSize, const double &tof, const bool &direct );





} // namespace variables
} // namespace bislip

#endif // BISLIP_VARIABLES_H
