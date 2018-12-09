#include "updateGuidance.h"
#include "bislipVariables.h"

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
    //std::cout << "Starting Aerodynamic guidance for this evaluation" << std::endl;

    if( currentTime == currentTime )
    {

        if( FlightConditions_ == nullptr )
        {
            FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                        bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
        }

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

        std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap_.at( vehicleName_ )->getVehicleSystems();

        const double AoA = bislip::variables::evaluateGuidanceInterpolator (
                    FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle ),
                    "Angle of Attack",
                    vehicleSystems,
                    FlightConditions_->getCurrentAltitude(),
                    FlightConditions_->getCurrentAirspeed(),
                    vehicleSystems->getE_max() );

        const double sigma = bislip::variables::evaluateGuidanceInterpolator (
                    FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle ),
                    "Bank Angle",
                    vehicleSystems,
                    FlightConditions_->getCurrentAltitude(),
                    FlightConditions_->getCurrentAirspeed(),
                    vehicleSystems->getE_max() );

        currentAngleOfAttack_ = tudat::unit_conversions::convertDegreesToRadians( AoA );
        currentBankAngle_ = tudat::unit_conversions::convertDegreesToRadians( sigma );
        currentAngleOfSideslip_ = 0;

    }

} // MyGuidance
} // namespace bislip

