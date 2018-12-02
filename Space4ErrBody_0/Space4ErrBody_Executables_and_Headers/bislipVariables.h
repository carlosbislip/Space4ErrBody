#ifndef BISLIP_VARIABLES_H
#define BISLIP_VARIABLES_H

#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <utility>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

namespace bislip {

namespace variables {

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

double computeAngleofAttack (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_AoA,
        const double &height,
        const double &airspeed,
        const double &E_max);

double computeThrottleSetting (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_throttle,
        const double &height,
        const double &airspeed,
        const double &E_max);

double computeThrustElevationAngle (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_eps_T,
        const double &height,
        const double &airspeed,
        const double &E_max);

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass);

} // namespace variables
} // namespace bislip

#endif // BISLIP_VARIABLES_H
