#include "updateGuidance.h"

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace tudat { namespace aerodynamics {


//! Trying to implement a simple aerodynamic guidance. Initially taken from the
//! TUDAT website.
//! http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/aerodynamicGuidance.html#FlightConditions
Eigen::Vector6d MyAerodynamicGuidance::getCoefficients( const double AoA, const double Mach )
{
    std::vector< double > CoefficientsInput_;
    Eigen::Vector6d newCoefficients;

    // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
    CoefficientsInput_.push_back( AoA );
    CoefficientsInput_.push_back( Mach );

    // Update and retrieve current aerodynamic coefficients
    coefficientInterface_->updateCurrentCoefficients( CoefficientsInput_ );
    newCoefficients = coefficientInterface_->getCurrentAerodynamicCoefficients( );


return newCoefficients;
}
} // namespace aerodynamics
} // namespace tudat
