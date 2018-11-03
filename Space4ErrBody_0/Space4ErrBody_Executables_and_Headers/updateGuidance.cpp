#include "Space4ErrBody.h"
#include "updateGuidance.h"
//#include "multi_dimensional_root_finding.hpp"
#include "getStuff.h"
#include "search_Coefficients.h"

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

namespace bislip { //namespace aerodynamics {

void MyGuidance::updateGuidance( const double currentTime )
{
    std::cout << "Starting Aerodynamic guidance for this evaluation" << std::endl;

    //! Set of parameters that I am yet to figure out how to pass/extract them around.
    const double c_ref = 13;
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
    double current_h = FlightConditions_->getCurrentAltitude( );
    std::cout << "current_h:  " << current_h << std::endl;
    double current_rho = FlightConditions_->getCurrentDensity( );
    double current_M = FlightConditions_->getCurrentMachNumber( );
    std::cout << "current_M:  " << current_M << std::endl;
    double current_gamma = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    double current_V = FlightConditions_->getCurrentAirspeed();
    std::cout << "current_V:  " << current_V << std::endl;
    double current_AoA = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );
    double S_ref = bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
    double m = bodyMap_.at( vehicleName_ )->getBodyMass( );//bodyMap_.at( vehicleName_ )->getCurrentMass( );
    const double delta_rad = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double chi_rad = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
    const double omega = bodyMap_.at( vehicleName_ )->getCentralBodyRotationRate( );
    const Eigen::Vector6d current_state = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getCurrentBodyCenteredBodyFixedState( );
    const Eigen::Vector3d r = current_state.segment( 0, 3 );
    const Eigen::Vector3d V = current_state.segment( 3, 3 );
    const double r_norm = r.norm( );
    const double V_norm = current_V;


    //interpolator_alpha_deg_ = interpolator_alpha_deg;
    //interpolator_eps_T_deg_ = interpolator_eps_T_deg;
    //interpolator_throttle_ = interpolator_throttle;
    //interpolator_alpha_rad_ = interpolator_alpha_rad;
    //interpolator_eps_T_rad_ = interpolator_eps_T_rad;

    //! Determine current gravity components. Recalculated because their particular form
    //! is used in the equations that will be solved further down. If these terms are analogous
    //! to what TUDAT generates, then this should be replaced.
    //Eigen::VectorXd gravs ( 2 );
    //
    const double g0 = bislip::getGravs ( mu, J2, J3, J4, R_E, r_norm, delta_rad ).norm();
    //const double g_n = gravs(0);
    //const double g_d = gravs(1);
    //const double g0 = g0_vector.norm();

    E_hat_ = ( g0 * current_h + 0.5 * current_V * current_V ) / E_max_;
    std::cout << "E_hat_" << E_hat_ << std::endl;

    //double throttle = interpolator_throttle_->interpolate( E_hat );
    currentAngleOfAttack_ = tudat::unit_conversions::convertDegreesToRadians( interpolator_alpha_deg_->interpolate( E_hat_ ) );
    currentBankAngle_ = 0;
    currentAngleOfSideslip_ = 0;

    // double alpha_deg = interpolator_alpha_deg_->interpolate( E );
    // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
    //Eigen::VectorXd coefficient_input ( 2 );
    // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
    std::vector< double > coefficient_input;
    coefficient_input.push_back( currentAngleOfAttack_ );
    coefficient_input.push_back( current_M );

    Eigen::Vector6d newCoefficients = MyGuidance::getCoefficients( coefficient_input );

    std::cout << std::fixed << std::setprecision(10) <<
                 std::setw(15) << "  E_hat:  " <<
                 std::setw(16) << E_hat_ <<
                 std::setw(15) << "  alpha_deg:  " <<
                 std::setw(16) << interpolator_alpha_deg_->interpolate( E_hat_ ) <<
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


}

void MyThrustGuidance::updateGuidance( const double currentTime )
{
    std::cout << "Starting Thrust guidance for this evaluation" << std::endl;

    currentEngineStatus_ = true;
    //std::cout << "currentEngineStatus_:  " << currentEngineStatus_ << std::endl;

    const double R_E = 6.378137e6;
    const double mu = 3.986004418e14;
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );

    double current_h = FlightConditions_->getCurrentAltitude( );
    std::cout << "current_h:  " << current_h << std::endl;

    double current_V = FlightConditions_->getCurrentAirspeed();
    std::cout << "current_V:  " << current_V << std::endl;
    const Eigen::Vector6d current_state = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getCurrentBodyCenteredBodyFixedState( );
    std::cout << "current_state:  " << current_state << std::endl;
    //const Eigen::Vector3d r = current_state.segment( 0, 3 );
    //std::cout << "r:  " << r << std::endl;
    const double r_norm = current_state.segment( 0, 3 ).norm( );
    std::cout << "r_norm:  " << r_norm << std::endl;
    const double delta_rad = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    std::cout << "delta_rad:  " << delta_rad << std::endl;

    const double g0 = bislip::getGravs ( mu, J2, J3, J4, R_E, r_norm, delta_rad ).norm();
    //const double g_n = gravs(0);
    //const double g_d = gravs(1);
    //const double g0 = g0_vector.norm();
    std::cout << "g0:  " << g0 << std::endl;

    double E_hat = ( g0 * current_h + 0.5 * current_V * current_V ) / E_max_;

    currentEpsilon_ = tudat::unit_conversions::convertDegreesToRadians( interpolator_eps_T_deg_->interpolate( E_hat ) );
    double throttle = interpolator_throttle_->interpolate( E_hat );
    std::cout << "throttle:  " << throttle << std::endl;

    double thrustMagnitude = 2500.0;

    currentThrustMagnitude_ = throttle * thrustMagnitude;

    currentSpecificImpulse_ = 472;

    //! Simplified expressions becuase thrust azimuth is known to be zero. I.e. phi_T = 0
    bodyFixedThrustDirection_( 0 ) = std::cos( currentEpsilon_ );
    bodyFixedThrustDirection_( 1 ) = 0.0;
    bodyFixedThrustDirection_( 2 ) = std::sin( currentEpsilon_ );

    std::cout << std::fixed << std::setprecision(10) <<
                 std::setw(15) << "  E_hat:  " <<
                 std::setw(16) << E_hat <<
                 std::setw(15) << "  alpha_deg:  " <<
                 std::setw(16) << interpolator_eps_T_deg_->interpolate( E_hat ) <<
                 std::setw(10) << "  ux:    " <<
                 std::setw(16) << bodyFixedThrustDirection_[0] <<
                 std::setw(10) << "  uy:    " <<
                 std::setw(16) << bodyFixedThrustDirection_[2] <<
                 std::setw(10) << "  uz:    " <<
                 std::setw(16) << bodyFixedThrustDirection_[4] <<
                 std::setw(10) << "  throttle: " <<
                 std::setw(16) << throttle <<
                 std::setw(10) << "  T:    " <<
                 std::setw(16) << currentThrustMagnitude_ << std::endl;

}

} // namespace aerodynamics
//} // namespace tudat
