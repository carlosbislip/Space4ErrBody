#include "updateGuidance.h"

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

        if( coefficientInterface_ == nullptr )
        {
            coefficientInterface_ = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                        bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( ) );
        }

        //! Set of parameters that I am yet to figure out how to pass/extract them around.
        const double b_ref = 13.0;
        const double c_ref = 23.0;
        const double del_x = 3.0;
        const double del_x_T = 0.0;
        const double del_z_T = -5.0;
        //const double R_E = 6.378137e6;
        //const double mu = 3.986004418e14;
        //const double J2 = 1082.626523e-6;
        //const double J3 = 2.532153e-7;
        //const double J4 = 1.6109876e-7;
        const double n = 1.2;


        //! Extract various parameters form current flight conditions
        double current_height = FlightConditions_->getCurrentAltitude( );
        //std::cout << "current_h:  " << current_h << std::endl;
        double currentDensity = FlightConditions_->getCurrentDensity( );
        //std::cout << "currentDensity:  " << currentDensity << std::endl;
        double currentDynamicPressure = FlightConditions_->getCurrentDynamicPressure();


        double current_M = FlightConditions_->getCurrentMachNumber( );
        //std::cout << "current_M:  " << current_M << std::endl;
        double currentFlightPathAngle = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        //std::cout << "currentFlightPathAngle:  " << currentFlightPathAngle << std::endl;
        double currentAirspeed = FlightConditions_->getCurrentAirspeed();
        //std::cout << "currentAirspeed:  " << currentAirspeed << std::endl;
        double current_AoA = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );
        //std::cout << "current_AoA:  " << current_AoA << std::endl;
        double surfaceArea = bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
        //std::cout << "surfaceArea:  " << surfaceArea << std::endl;
        double currentMass = bodyMap_.at( vehicleName_ )->getBodyMass( );//bodyMap_.at( vehicleName_ )->getCurrentMass( );
        //std::cout << "m:  " << m << std::endl;
        //double q_dyn = ( 0.5 ) * currentDensity * currentAirspeed * currentAirspeed;

        const double currentLatitude = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        const double currentHeading = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        //const double omega = 7.292115*1E-5;//bodyMap_.at( vehicleName_ )->getCentralBodyRotationRate( );
        const Eigen::Vector6d current_state = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getCurrentBodyCenteredBodyFixedState( );
        const Eigen::Vector3d r = current_state.segment( 0, 3 );
        const Eigen::Vector3d V = current_state.segment( 3, 3 );
        const double r_norm = r.norm( );
        const double V_norm = currentAirspeed;

        std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap_.at( vehicleName_ )->getVehicleSystems();
        std::shared_ptr< bislip::VehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();

        double angleOfAttack = tudat::unit_conversions::convertDegreesToRadians(
                    bislip::Variables::evaluateGuidanceInterpolator (
                        bislip::Parameters::Optimization::AngleOfAttack,
                        bislipSystems,
                        FlightConditions_->getCurrentAltitude(),
                        FlightConditions_->getCurrentAirspeed(),
                        bislipSystems->getE_max() ) );

        double bankAngle = tudat::unit_conversions::convertDegreesToRadians(
                    bislip::Variables::evaluateGuidanceInterpolator (
                        bislip::Parameters::Optimization::BankAngle,
                        bislipSystems,
                        FlightConditions_->getCurrentAltitude(),
                        FlightConditions_->getCurrentAirspeed(),
                        bislipSystems->getE_max() ) );

        double thrustElevationAngle = tudat::unit_conversions::convertDegreesToRadians(
                    bislip::Variables::evaluateGuidanceInterpolator (
                        bislip::Parameters::Optimization::ThrustElevationAngle,
                        bislipSystems,
                        FlightConditions_->getCurrentAltitude(),
                        FlightConditions_->getCurrentAirspeed(),
                        bislipSystems->getE_max() ) );

        double thrustAzimuthAngle = tudat::unit_conversions::convertDegreesToRadians(
                    bislip::Variables::evaluateGuidanceInterpolator (
                        bislip::Parameters::Optimization::ThrustAzimuthAngle,
                        bislipSystems,
                        FlightConditions_->getCurrentAltitude(),
                        FlightConditions_->getCurrentAirspeed(),
                        bislipSystems->getE_max() ) );

        //! Update and retrieve current aerodynamic coefficients
        Eigen::Vector6d currentCoefficients = getCurrentCoefficients( );

        double maxBankAngle = bislip::Variables::computeEquilibriumGlideLimit( bodyMap_, vehicleName_, centralBodyName_ );
                  //  FlightConditions_,
                  //  coefficientInterface_,
                  //  bodyMap_ );
                    //bodyMap_.at( vehicleName_ )->getBislipSystems( ),
                    //bodyMap_.at( vehicleName_ )->getBodyMass( ) );

        if ( bankAngle > maxBankAngle ){ bankAngle = maxBankAngle; }

        currentAngleOfAttack_ = angleOfAttack ;
        currentBankAngle_ =  bankAngle;
        currentAngleOfSideslip_ = 0;


        double del_C_m_b = getBodyflapCmIncrement( );

        double  bodyFlap = bislip::Variables::getBodyFlapDeflection( del_C_m_b );
        if ( currentFlightPathAngle > 0 ){ bodyFlap = 0; }
        if ( currentFlightPathAngle < 0 ){ bodyFlap = 20; }

        double bodyFlapDeflection = tudat::unit_conversions::convertDegreesToRadians( bodyFlap );

        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bodyFlapDeflection );



    }

} // MyGuidance
} // namespace bislip
