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
void ValidationAerodynamicGuidance::updateGuidance(const double currentTime)
{



//////////////// Calc Heading error
//! Extract current latitude
//! //const double lat_c_rad = FlightConditions_->getCurrentGeodeticLatitude( );
const double lat_c_rad = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::latitude_angle );

//! Extract current longitude
//! How equivalent are these functions?
//const double lon_c_rad = FlightConditions_->getCurrentLongitude( );
const double lon_c_rad = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::longitude_angle );

//! Extract current heading angle
const double chi_c_rad = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::heading_angle );
const double chi_c_deg = unit_conversions::convertRadiansToDegrees( chi_c_rad );

//! Calculate required heading angle: https://www.movable-type.co.uk/scripts/latlong.html
//! Maybe consider Vicenty's formulation: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
const double chi_req_rad = std::atan2( std::sin( lon_f_rad_ - lon_c_rad ) * std::cos( lat_f_rad_ ) , std::cos( lat_c_rad ) * std::sin( lat_f_rad_ ) - std::sin( lat_c_rad ) * std::cos( lat_f_rad_ ) * std::cos( lon_f_rad_ - lon_c_rad ) );
const double chi_req_deg = unit_conversions::convertRadiansToDegrees( chi_req_rad );

//! Calculate heading error
const double chi_err_deg = chi_c_deg - chi_req_deg ;
const double abs_chi_err_deg = abs( chi_err_deg );

//////////////// Calc Distance to target
//! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
//! Maybe consider Vicenty's formulation: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
const double a = std::sin( (lat_f_rad_ - lat_c_rad) / 2) * std::sin( (lat_f_rad_ - lat_c_rad) / 2) + std::cos( lat_c_rad ) * std::cos( lon_c_rad ) * std::sin( (lon_f_rad_ - lon_c_rad) / 2) * std::sin( (lon_f_rad_ - lon_c_rad) / 2);
const double d_rad = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
const double d_deg = unit_conversions::convertRadiansToDegrees( d_rad );

int sign = 1;
if ( currentBankAngle_ > 0.0 )
{
    sign = 1;
}
else if ( currentBankAngle_ < 0.0 )
{
    sign = -1;
}

double reversal_conditional = chi_err_deg * sign;
bool reversal = false;
//////////////// Determine SIGN of bank angle
//! I believe I may need some time history to properly determine if a sign
//! change is required. Or should this also be somehow included in the
//! determination of the bank angle?
if ( ( d_deg >= 30 ) && ( abs_chi_err_deg > 15 ) && ( reversal_conditional < 0 ) )
{
    reversal = true;
}
else if ( ( d_deg >= 30 ) && ( abs_chi_err_deg > 15 ) && ( reversal_conditional > 0 ) )
{
    reversal = false;
}
else if ( ( d_deg >= 30 ) && ( abs_chi_err_deg < 15 ) && ( reversal_conditional > 0 ) )
{
    reversal = false;
}
else if ( ( d_deg >= 30 ) && ( abs_chi_err_deg < 15 ) && ( reversal_conditional < 0 ) )
{
    reversal = false;
}
else if ( ( d_deg < 30 ) && ( d_deg > 10 ) && ( abs_chi_err_deg > abs( 23 + (d_deg - 10) * ( 15 - 23 ) / ( 30 - 10 )) ) && ( reversal_conditional < 0 ) )
{
    reversal = true;
}
else if ( ( d_deg < 30 ) && ( d_deg > 10 ) && ( abs_chi_err_deg > abs( 23 + (d_deg - 10) * ( 15 - 23 ) / ( 30 - 10 )) ) && ( reversal_conditional > 0 ) )
{
    reversal = false;
}
else if ( ( d_deg < 30 ) && ( d_deg > 10 ) && ( abs_chi_err_deg < abs( 23 + (d_deg - 10) * ( 15 - 23 ) / ( 30 - 10 )) ) && ( reversal_conditional > 0 ) )
{
    reversal = false;
}


//else if ( ( currentTime >= 924 ) && ( currentTime <= 1319 ) )


    /*
00.0		& 23.0 &  0.85776  & 14.35438\\
10.0		& 23.0 &  0.86476  & 15.20725\\
30.0		& 15.0 &  0.90739  & 19.62079\\
60.0		& 15.0 &  0.91452  & 20.36779\\
            &      &  0.94331  & 23.00000\\ \bottomrule


  */

//const double d_haversine = c * spice_interface::getAverageRadius( "Earth" );
double currentBankAngle_1;


//! Determine value of aerodynamic angles
if( currentTime < 264 )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_1 = 0.0;
}
else if ( ( currentTime >= 264 ) && ( currentTime <= 290 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_1 = 0.0;
}
else if ( ( currentTime >= 290 ) && ( currentTime <= 554 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_1 = 79.6;
}
else if ( ( currentTime >= 554 ) && ( currentTime <= 686 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_1 = 56.0;
}
else if ( ( currentTime >= 686 ) && ( currentTime <= 924 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_1 = 59.8;
}
else if ( ( currentTime >= 924 ) && ( currentTime <= 1319 ) )
{
    currentAngleOfAttack_ = 40.0;
    currentBankAngle_1 = 59.8;
}
else
{
    currentAngleOfAttack_ = 11.5;
    currentBankAngle_1 = 54.0;
}
currentBankAngle_ = currentBankAngle_1;

if ( reversal == true )
{
   currentBankAngle_ = -currentBankAngle_1;
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
