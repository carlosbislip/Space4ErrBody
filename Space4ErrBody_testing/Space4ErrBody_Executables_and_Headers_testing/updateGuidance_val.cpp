#include "Space4ErrBody.h"
#include "updateGuidance_val.h"

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace tudat
{

namespace aerodynamics
{


//! Trying to implement a simple aerodynamic guidance. Initially taken from the
//! TUDAT website.
//! http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/aerodynamicGuidance.html#FlightConditions
void ValidationAerodynamicGuidance::updateGuidance( const double currentTime )
{
if( currentTime < 264 )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_ = 0.0;
}
else if ( ( currentTime >= 264 ) && ( currentTime <= 290 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_ = 0.0;
}
else if ( ( currentTime >= 290 ) && ( currentTime <= 554 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_ = 79.6;
}
else if ( ( currentTime >= 554 ) && ( currentTime <= 686 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_ = 56.0;
}
else if ( ( currentTime >= 686 ) && ( currentTime <= 924 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_ = 59.8;
}
else if ( ( currentTime >= 924 ) && ( currentTime <= 1319 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_ = 59.8;
}
else
{
    currentAngleOfAttack_ = 11.5;
    currentBankAngle_ = 54.0;
}

currentAngleOfAttack_ = unit_conversions::convertDegreesToRadians( currentAngleOfAttack_ ) ;
currentBankAngle_ = unit_conversions::convertDegreesToRadians( currentBankAngle_ ) ;



    //using namespace tudat;
    //using namespace tudat::aerodynamics;
    //using namespace tudat::basic_mathematics;
   // using namespace tudat::mathematical_constants;
   // using namespace tudat::simulation_setup;
   // using namespace tudat::unit_conversions;


    // Define aerodynamic coefficient interface/flight conditions (typically retrieved from body map; may also be a member variable)
    //boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_->getAerodynamicCoefficientInterface( );
  //  boost::shared_ptr< aerodynamics::FlightConditions > FlightConditions_;//->getFlightConditions( );
/*
    // Compute angles of attack and sideslip
    //currentAngleOfAttack_ = ...
    //currentAngleOfSideslip_ = ...



   //boost::shared_ptr< aerodynamics::FlightConditions > flightConditions_ = double currentFlightPathAngle = flightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );

   // double currentFlightPathAngle = flightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );




   // currentAngleOfAttack_ = unit_conversions::convertDegreesToRadians( 35.0 );
   // currentBankAngle_ = unit_conversions::convertDegreesToRadians( 35.0 );
    //double constantAngleOfAttack = unit_conversions::convertDegreesToRadians( 30 );
    //double constantBankAngle = unit_conversions::convertDegreesToRadians( 85 );
    //



    if( FlightConditions_->getCurrentAltitude( ) > 60.0E3 )
    {
        currentAngleOfAttack_ = unit_conversions::convertDegreesToRadians( 35.0 ) ;
    }
    else if( FlightConditions_->getCurrentAltitude( ) < 25.0E3 )
    {
        currentAngleOfAttack_ = unit_conversions::convertDegreesToRadians( 5.0 ) ;
    }
    else
    {
        currentAngleOfAttack_ = unit_conversions::convertDegreesToRadians(
                    ( 5.0 + 30.0 * ( FlightConditions_->getCurrentAltitude( ) - 25.0E3 ) / 35.0E3 ) );
    }

    currentAngleOfSideslip_ = 0.0;

    if( FlightConditions_->getCurrentMachNumber( ) > 8 )
    //   if( FlightConditions_->getCurrentAltitude( ) < 80.0E3 )
    {
        currentBankAngle_ = unit_conversions::convertDegreesToRadians( 10.0 ) ;
    }
    else if( ( FlightConditions_->getCurrentMachNumber( ) > 3) && ( FlightConditions_->getCurrentMachNumber( ) < 8) )
    {
        currentAngleOfAttack_ = unit_conversions::convertDegreesToRadians( 45.0 ) ;
    }
    else
    {
        currentBankAngle_ = unit_conversions::convertDegreesToRadians( 85.0 );
    }
*/
        /*
        // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
        std::vector< double > currentAerodynamicCoefficientsInput_;
        currentAerodynamicCoefficientsInput_.push_back( currentAngleOfAttack_ );
        currentAerodynamicCoefficientsInput_.push_back( currentAngleOfSideslip_ );
        currentAerodynamicCoefficientsInput_.push_back( currentBankAngle_ );
        currentAerodynamicCoefficientsInput_.push_back( FlightConditions_->getCurrentMachNumber( ) );

        // Update and retrieve current aerodynamic coefficients
        coefficientInterface_->updateCurrentCoefficients( currentAerodynamicCoefficientsInput_ );
        Eigen::Vector3d currentAerodynamicCoefficients = coefficientInterface_->getCurrentForceCoefficients( );

        // Compute bank angle
        currentBankAngle_ =  some function of currentAerodynamicCoefficients

        */
}


} // namespace aerodynamics

} // namespace tudat
