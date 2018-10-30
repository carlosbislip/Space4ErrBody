#include "Space4ErrBody.h"
#include "updateGuidance.h"
#include "multi_dimensional_root_finding.hpp"
#include "getStuff.h"
#include "search_Coefficients.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Dense>

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"

namespace tudat { namespace aerodynamics {
void MyAerodynamicGuidance::updateGuidance( const double currentTime )
{
    //! Set of parameters that I am yet to figure out how to pass/extract them around.
    const double c_ref = 13;
    const double Isp = 472;
    const double del_x = 3;
    const double del_x_T = 0;
    const double del_z_T = -5;
    const double n = 1.2;
    const double R_E = 6.378137e6;
    const double mu = 3.986004418e14;
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;

    //! Extract various parameters form current flight conditions
    double current_h = FlightConditions_->getCurrentAltitude( );
    double current_rho = FlightConditions_->getCurrentDensity( );
    double current_M = FlightConditions_->getCurrentMachNumber( );
    double current_gamma = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );
    double current_V = FlightConditions_->getCurrentAirspeed();
    double current_AoA = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::angle_of_attack );
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

    /*
    alpha_rad_;
    eps_T_rad_ = eps_T_rad;
    m_dot_ = m_dot;
    gamma_rad_ = gamma_rad;
    sigma_rad_ = sigma_rad;
    height_ = height;
    Va_ = Va;
    */


    //! Determine current gravity components. Recalculated because their particular form
    //! is used in the equations that will be solved further down. If these terms are analogous
    //! to what TUDAT generates, then this should be replaced.
    Eigen::VectorXd gravs ( 2 );
    gravs << bislip::getGravs ( mu, J2, J3, J4, R_E, r_norm, delta_rad );
    const double g_n = gravs(0);
    const double g_d = gravs(1);

    //! Determine surface gravity.
    Eigen::VectorXd g0_vector ( 2 );
    g0_vector << bislip::getGravs ( mu, J2, J3, J4, R_E, R_E, delta_rad );
    const double g0 = g0_vector.norm();

    //! Calculate the flight-path angle required for a heading rate of zero. the function uses
    //! a root finder with derivatives. Various options are available:
    //!     Newton-Raphson - 1st order
    //!     Halley         - 2nd order
    //!     Schroder       - 2nd order
    //!
    //! Values are implicitly bounded to [0,pi/2]. The result may actually be < 0, but the
    //! returned value is always positive.
    const double s_delta   = std::sin( delta_rad );
    const double c_delta   = std::cos( delta_rad );
    const double s_chi     = std::sin( chi_rad );
    const double c_chi     = std::cos( chi_rad );
    const double A         = 2 * ( omega * r_norm * c_delta / current_V ) * sqrt( 1 - pow( c_delta * s_chi , 2 ) ) / ( (s_delta / c_delta ) * s_chi );
    const double B         = 1 + ( g_n * r_norm ) / ( current_V * current_V  * (s_delta / c_delta ) ) + pow( omega * r_norm * c_delta / current_V , 2 );
    const double phi       = std::atan2( 2 * omega * current_V * s_delta , -2 * omega * current_V * c_delta * c_chi );
    const double new_gamma = bislip::getGamma( A, B, phi );

    //! Rotate current velocity vector to what it has to be. Not sure what to do with this...
    //! but I do know that this is the orientation I want for the vehicle. Maybe its the BODY
    //! frame that must be rotated. Dont know.
    const double dif_gamma = new_gamma - current_gamma;
    Eigen::Matrix3d C2_gamma;
    const double s_dif_gamma = std::sin( dif_gamma );
    const double c_dif_gamma = std::cos( dif_gamma );
    C2_gamma<< c_dif_gamma,0,-s_dif_gamma,0,1,0,s_dif_gamma,0,c_dif_gamma;
    const Eigen::Vector3d new_V = C2_gamma * V;

    //! Find the angle of attack that generates a Cm = 0. Is there a better way to search the table?
    double eps_c = 1e-9;
    double a, b, c;
    double func_a, func_b, func_c;
    Eigen::Vector6d a_vect, b_vect, c_vect;
    a = 0;
    b = 50;

    //! Bisection Method!!!
    do {

        a_vect = MyAerodynamicGuidance::getCoefficients( a , current_M );
        b_vect = MyAerodynamicGuidance::getCoefficients( b , current_M );
        func_a = a_vect[ 4 ];
        func_b = b_vect[ 4 ];

        //! Find middle point
        c = ( a + b ) / 2;
        c_vect = MyAerodynamicGuidance::getCoefficients( c , current_M );
        func_c = c_vect[ 4 ];

        // Check if middle point is root
        if ( func_c == 0.0 )
            break;

        // Decide the side to repeat the steps
        else if ( func_c * func_a < 0 )
            b = c;
        else
            a = c;
    } while ( ( b - a ) >= eps_c );

    //! Re-assign values to readable variables.ß
    const double AoA = c;
    const Eigen::Vector6d newCoefficients = c_vect;

    //! Calculate various constant terms used in the EoM.
    const double q_d = current_rho * current_V * current_V / 2;
    const double s_new_gamma = std::sin( new_gamma );
    const double c_new_gamma = std::cos( new_gamma );
    const double extra1 = omega * omega * r_norm * c_delta * ( s_new_gamma * c_delta - c_new_gamma * s_delta * c_chi );
    const double extra2 = ( g_d * s_new_gamma + g_n * c_new_gamma * c_chi );
    const double extra = m * ( extra1 - extra2 );

    //! Copy parameters to new vector. This cleans up the clutter and avoids typing mistakes.
    Eigen::VectorXd p ( 12 );
    p << del_x, del_x_T, del_z_T, S_ref, c_ref, Isp, m, q_d, current_M, n, g0, extra;

    //! Get required thrust elevation angle. This function incorporates various concepts and is
    //! supposed satisfy various contraints:
    //!     Flight-path angle such that Heading angle remains constant
    //!     Angle of attack such that Cm = 0
    //!     Mass rate that ensures pitch trim
    //!
    //! I want to be able to pass the previous eps_T. I believe this may help guide the search
    //! for the current one.
    const double eps_T = bislip::getEps_T( AoA, p, newCoefficients );

    //std::cout << "gamma: " << new_gamma * 180 / tudat::mathematical_constants::PI<< std::endl;
    //std::cout << "AoA:   " << AoA * 180 / tudat::mathematical_constants::PI<< std::endl;
    //std::cout << "Cm:    " << newCoefficients[4] << std::endl;
    //std::cout << "eps_T: " << eps_T * 180 / tudat::mathematical_constants::PI << std::endl;

    //! Calculate new Mass rate
    const double m_dot = ( ( q_d * S_ref ) / ( g0 * Isp ) ) * ( del_x * ( newCoefficients[0] * std::sin( AoA ) - newCoefficients[2] * std::cos( AoA ) ) - newCoefficients[4] * c_ref ) / ( del_x_T * std::sin( eps_T ) + del_z_T * std::cos( eps_T ) );


    std::cout << std::fixed << std::setprecision(10) <<
                 std::setw(8) << "time:  " <<
                 std::setw(16) << currentTime <<
                 std::setw(8) << "gamma: " <<
                   std::setw(16) << new_gamma * 180 / tudat::mathematical_constants::PI <<
                   std::setw(8) << "new_V: " <<
                   std::setw(16) << new_V(0)  << "  " <<new_V(1)<< "  " <<new_V(2) <<
                   std::setw(8) << "AoA:   "<<
                   std::setw(16) << AoA * 180 / tudat::mathematical_constants::PI<<
                   std::setw(8) << "Cm:    " <<
                   std::setw(16) << newCoefficients[4] <<
                   std::setw(8) << "eps_T: " <<
                   std::setw(16) << eps_T * 180 / tudat::mathematical_constants::PI <<
                 std::setw(8) << "m_dot: " <<
                 std::setw(16) << m_dot  <<std::endl;



//Eigen::Matrix3d C1_VB;
//C1_VB<< 1,0,0,0,std::cos(  ),std::sin(  ), 0,-std::cos(  ),std::cos(  );
//Eigen::Matrix3d C2;
//C2<< std::cos(  ),0,-std::sin(  ),0,1,0,std::sin(  ),0,std::cos(  );
    //Eigen::Matrix3d C3;
    //C3<< std::cos(  ),0,-std::sin(  ),0,1,0,std::sin(  ),0,std::cos(  );




/*
    Eigen::Matrix3d C2_AB;
    C2_AB<< std::cos( -AoA ),0,-std::sin( -AoA ),0,1,0,std::sin( -AoA ),0,std::cos( -AoA );
    Eigen::Matrix3d C3_AB;
    C3_AB<< std::cos( 0 ),0,-std::sin( 0 ),0,1,0,std::sin( 0 ),0,std::cos( 0 );

    Eigen::Matrix3d C2_TV;
    C2_TV<< std::cos( new_gamma ),0,-std::sin( new_gamma ),0,1,0,std::sin( new_gamma ),0,std::cos( new_gamma );
    Eigen::Matrix3d C3_TV;
    C3_TV<< c_chi,0,-s_chi,0,1,0,s_chi,0,c_chi;

    Eigen::Matrix3d C1_TA;
    C1_TA<< 1,0,0,0,std::cos( 0 ),std::sin( 0 ), 0,-std::cos( 0 ),std::cos( 0 );

    Eigen::Matrix3d CAB, CTA, CTV, CVT, CVA;

    CAB = ( C3_AB.inverse() ) * ( C2_AB.inverse() );
    CTA = C1_TA;
    CTV = C2_TV * C3_TV;
    CVT = CTV.inverse();
    CVA = CVT * CTA;
*/


    currentAngleOfAttack_ = 10 * tudat::mathematical_constants::PI / 180;
    currentBankAngle_ = 0;
    currentAngleOfSideslip_ = 0;
    //currentRotationFromBodyToTrajectoryFrame_ = CTV * CVA * CAB;


/*

    do {

        for ( unsigned i = 0; i < 2 ; ++i ){
            for ( unsigned j = 0; j < 2 ; j++ ){
                func_plus  = getF ( i, j, p,  eps, x_old );
                func_minus = getF ( i, j, p, -eps, x_old );
                J( i, j ) = ( func_plus - func_minus ) / ( 2 * eps );
            }
        }

        f(0) = getF ( 0, 0, p, 0, x_old );
        f(1) = getF ( 1, 0, p, 0, x_old );

        invJ = J.inverse();
        invJf = invJ * f;
        x_new = x_old - invJf;
        dif = x_new - x_old;
        x_old = x_new;

        std::cout << "Jacobian: " << J(0,0) << "   " << J(0,1) << "   "<< J(1,0) << "   "<< J(1,1) << std::endl;
        std::cout << "Inverse:  " << invJ(0,0) << "   " << invJ(0,1) << "   "<< invJ(1,0) << "   "<< invJ(1,1) << std::endl;
        std::cout << "invJ * f: " << invJf(0) << "   " << invJf(1) << std::endl;
        std::cout << "f:        " << f[0] << "   " << f[1]<< std::endl;
        std::cout << "x_new:    " << x_new[0] << "   " << x_new[1] << std::endl;
        std::cout << "dif:      " << dif[0] << "   " << dif[1] << std::endl;

iter = iter + 1;
std::cout << "iteration: " << iter << std::endl;
std::cout << "dif_norm: " << dif.norm()  << std::endl;


sleep(2);

if ( x_old(0) < 0 )
{
    x_old(0) = 0;
}
if ( x_old(0) > 45.0 * tudat::mathematical_constants::PI / 180.0 )
{
    x_old(0) = 45.0 * tudat::mathematical_constants::PI / 180.0;
}

if ( std::abs( x_old(1) ) > 20.0 * tudat::mathematical_constants::PI / 180.0 )
{
    if ( x_old(1) < 0 )
    {
        sgn = -1;
    }
    x_old(1) = sgn * 20.0 * tudat::mathematical_constants::PI / 180.0;
}

    } while( dif.norm() > 1e-1 );
*/


    //vehicleSystems_->setCurrentControlSurfaceDeflection( "Elevon", elevonDeflection );
    //vehicleSystems_->setCurrentControlSurfaceDeflection( "Aileron1", aileron1Deflection );
    //vehicleSystems_->setCurrentControlSurfaceDeflection( "Aileron2", aileron2Deflection );


    // Define aerodynamic coefficient interface/flight conditions (typically retrieved from body map; may also be a member variable)
    //boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_->getAerodynamicCoefficientInterface( );
  //  boost::shared_ptr< aerodynamics::FlightConditions > FlightConditions_;//->getFlightConditions( );
/*

   //boost::shared_ptr< aerodynamics::FlightConditions > flightConditions_ = double currentFlightPathAngle = flightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );

   // double currentFlightPathAngle = flightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );


        // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
        std::vector< double > CoefficientsInput_;
        CoefficientsInput_.push_back( currentAngleOfAttack_ );
        CoefficientsInput_.push_back( FlightConditions_->getCurrentMachNumber( ) );

        // Update and retrieve current aerodynamic coefficients
        coefficientInterface_->updateCurrentCoefficients( CoefficientsInput_ );
        Eigen::Vector3d currentAerodynamicCoefficients = coefficientInterface_->getCurrentForceCoefficients( );

        // Compute bank angle
        currentBankAngle_ =  some function of currentAerodynamicCoefficients

        */
}


} // namespace aerodynamics

} // namespace tudat
