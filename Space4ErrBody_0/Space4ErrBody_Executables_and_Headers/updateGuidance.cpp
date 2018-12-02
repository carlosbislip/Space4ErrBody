//#include "Space4ErrBody.h"
#include "updateGuidance.h"
//#include "multi_dimensional_root_finding.hpp"
//#include "getStuff.h"
//#include "search_Coefficients.h"
#include "bislipVariables.h"

//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <Eigen/Dense>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"

namespace bislip {

void MyGuidance::updateGuidance( const double currentTime )
{
    std::cout << "Starting Aerodynamic guidance for this evaluation" << std::endl;

    if( currentTime == currentTime )
    {

        if( FlightConditions_ == nullptr )
        {
            FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                        bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
        }


        // E_hat = getE_hat();
        //std::cout << "E_hat: " << E_hat << std::endl;

        //       if( coefficientInterface_ == nullptr )
        //       {
        //           coefficientInterface_ = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
        //                       bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( ) );
        //       }


        //! Set of parameters that I am yet to figure out how to pass/extract them around.
        const double b_ref = 13;
        const double c_ref = 23;
        const double del_x = 3;
        const double del_x_T = 0;
        const double del_z_T = -5;
        const double n = 1.2;
        const double R_E = 6.378137e6;
        const double mu = 3.986004418e14;
        const double J2 = 1082.626523e-6;
        const double J3 = 2.532153e-7;
        const double J4 = 1.6109876e-7;
        //std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
        //                bodyMap_.at( vehicleName_ )->getFlightConditions( ) );

        //! Extract various parameters form current flight conditions
        double current_height = FlightConditions_->getCurrentAltitude( );
        //std::cout << "current_h:  " << current_h << std::endl;
        double current_rho = FlightConditions_->getCurrentDensity( );
        //std::cout << "current_rho:  " << current_rho << std::endl;
        double current_M = FlightConditions_->getCurrentMachNumber( );
        //std::cout << "current_M:  " << current_M << std::endl;
        double current_gamma = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        //std::cout << "current_gamma:  " << current_gamma << std::endl;
        double current_V = FlightConditions_->getCurrentAirspeed();
        //std::cout << "current_V:  " << current_V << std::endl;
        double current_AoA = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );
        //std::cout << "current_AoA:  " << current_AoA << std::endl;
        double S_ref = bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
        //std::cout << "S_ref:  " << S_ref << std::endl;
        double currentMass = bodyMap_.at( vehicleName_ )->getBodyMass( );//bodyMap_.at( vehicleName_ )->getCurrentMass( );
        //std::cout << "m:  " << m << std::endl;
        const double delta_rad = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        const double chi_rad = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        const double omega = 7.292115*1E-5;//bodyMap_.at( vehicleName_ )->getCentralBodyRotationRate( );
        const Eigen::Vector6d current_state = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getCurrentBodyCenteredBodyFixedState( );
        const Eigen::Vector3d r = current_state.segment( 0, 3 );
        const Eigen::Vector3d V = current_state.segment( 3, 3 );
        const double r_norm = r.norm( );
        const double V_norm = current_V;

        double AoA = bislip::variables::computeThrottleSetting(
                    bodyMap_.at( vehicleName_ )->getAoAInterpolator(),
                    FlightConditions_->getCurrentAltitude(),
                    FlightConditions_->getCurrentAirspeed(),
                    bodyMap_.at( vehicleName_ )->getE_max( ) );

        if ( AoA < parameterBounds_[ 2 ] )
        {
            AoA = parameterBounds_[ 2 ];
        }
        if ( AoA > parameterBounds_[ 3 ] )
        {
            AoA = parameterBounds_[ 3 ];
        }

        currentAngleOfAttack_ = tudat::unit_conversions::convertDegreesToRadians( AoA );
        currentBankAngle_ = 0;
        currentAngleOfSideslip_ = 0;

        /*        //currentSpecificImpulse_ = Isp_;

        //! Simplified expressions becuase thrust azimuth is known to be zero. I.e. phi_T = 0
       // currentbodyFixedThrustDirection_( 0 ) = std::cos( tudat::unit_conversions::convertDegreesToRadians( eps_T ) );
       // currentbodyFixedThrustDirection_( 1 ) = 0.0;
       // currentbodyFixedThrustDirection_( 2 ) = std::sin( tudat::unit_conversions::convertDegreesToRadians( eps_T ) );

        //currentThrustMagnitude_ = throttle * maxThrust_;

        //std::cout << "currentMass: " << currentMass << std::endl;
        //std::cout << "currentbodyFixedThrustDirection_: " << currentbodyFixedThrustDirection_ << std::endl;
        //std::cout << "currentThrustMagnitude_: " << currentThrustMagnitude_ << std::endl;
        //std::cout << "currentEngineStatus_: " << currentEngineStatus_ << std::endl;
        //std::cout << "currentSpecificImpulse_: " << currentSpecificImpulse_ << std::endl;
        //std::cout << "----------------------------"<< std::endl;

        // double alpha_deg = interpolator_alpha_deg_->interpolate( E );
        // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
        //Eigen::VectorXd coefficient_input ( 2 );
        // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
        std::vector< double > coefficient_input;
        coefficient_input.push_back( currentAngleOfAttack_ );
        coefficient_input.push_back( current_M );

        //       Eigen::Vector6d newCoefficients = MyGuidance::getCoefficients( coefficient_input );

        std::cout << std::fixed << std::setprecision(10) <<
                     std::setw(15) << "  E_hat:  " <<
                     std::setw(16) << E_hat <<
                     std::setw(15) << "  alpha_deg:  " <<
                     std::setw(16) << interpolator_alpha_deg_->interpolate( E_hat ) <<
                     //std::setw(15) << "  eps_T_deg: " <<
                     //std::setw(16) << eps_T_deg <<
                     //std::setw(15) << "  eps_T_deg: " <<
                     //std::setw(16) << eps_T_deg <<
                     //std::setw(8) << "new_V: " <<
                     //std::setw(16) << new_V(0)  << "  " <<new_V(1)<< "  " <<new_V(2) <<
                     //std::setw(15) << "  throttle:   "<<
                     //std::setw(16) << throttle_ <<
                     std::setw(10) << "  CD:    " <<
                     std::setw(16) << newCoefficients[0] <<
                     std::setw(10) << "  CL:    " <<
                     std::setw(16) << newCoefficients[2] <<
                     std::setw(10) << "  Cm:    " <<
                     std::setw(16) << newCoefficients[4] <<
                     std::setw(10) << "  Vx:    " <<
                     std::setw(16) << V[0] <<
                     std::setw(10) << "  Vy:    " <<
                     std::setw(16) << V[1] <<
                     std::setw(10) << "  Vz:    " <<
                     std::setw(16) << V[2] << std::endl;

*/
    }

}
/*
Eigen::Vector3d  MyGuidance::getCurrentBodyFixedThrustDirection( )
{

    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    double eps_T = bislip::variables::computeThrustElevationAngle(
                bodyMap_.at( vehicleName_ )->getThrustElevationAngleInterpolator(),
                FlightConditions_->getCurrentAltitude(),
                FlightConditions_->getCurrentAirspeed(),
                bodyMap_.at( vehicleName_ )->getE_max( ) );

    if ( eps_T < parameterBounds_[ 4 ] )
    {
        eps_T = parameterBounds_[ 4 ];
    }
    if ( eps_T > parameterBounds_[ 5 ] )
    {
        eps_T = parameterBounds_[ 5 ];
    }

    //! Simplified expressions becuase thrust azimuth is known to be zero. I.e. phi_T = 0
    currentbodyFixedThrustDirection_( 0 ) = std::cos( tudat::unit_conversions::convertDegreesToRadians( eps_T ) );
    currentbodyFixedThrustDirection_( 1 ) = 0.0;
    currentbodyFixedThrustDirection_( 2 ) = std::sin( tudat::unit_conversions::convertDegreesToRadians( eps_T ) );

    return currentbodyFixedThrustDirection_;
}

double  MyGuidance::getCurrentThrustMagnitude( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    double throttle = bislip::variables::computeThrottleSetting(
                bodyMap_.at( vehicleName_ )->getThrottleInterpolator(),
                FlightConditions_->getCurrentAltitude(),
                FlightConditions_->getCurrentAirspeed(),
                bodyMap_.at( vehicleName_ )->getE_max( ) );

    if ( throttle < 0.0 )
    {
        throttle = 0.0;
    }
    if ( throttle > 1.0 )
    {
        throttle = 1.0;
    }

  //  currentThrustMagnitude_ = throttle * maxThrust_;

    return throttle * maxThrust_;
}
double  MyGuidance::getCurrentSpecificImpulse( )
{
    //currentSpecificImpulse_ = Isp_;

    return Isp_;
}
bool  MyGuidance::getCurrentEngineStatus( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    //currentEngineStatus_ = bislip::variables::determineEngineStatus(  bodyMap_.at( vehicleName_ )->getBodyMass(), bodyMap_.at( vehicleName_ )->getLandingMass() );

    return bislip::variables::determineEngineStatus(  bodyMap_.at( vehicleName_ )->getBodyMass(), bodyMap_.at( vehicleName_ )->getLandingMass() );
}*/



} // namespace bislip

