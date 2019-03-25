#include "updateGuidance.h"

namespace bislip {

void MyGuidance::updateGuidance( const double currentTime )
{
    //std::cout << "Starting Aerodynamic guidance for this evaluation" << std::endl;
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap_.at( vehicleName_ )->getVehicleSystems();
    double startingEpoch = bislipSystems->getStartingEpoch();
    double fixedStepSize = bislipSystems->getFixedStepSize();
    double samplingRatio = bislipSystems->getSamplingRatio();
    std::string currentTrajectoryPhase = bislipSystems->getCurrentTrajectoryPhase();

    if ( currentTime == currentTime )
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
    }

    if ( ( currentTime >= startingEpoch ) && ( ( currentTime < startingEpoch + fixedStepSize * samplingRatio ) ) )
    {

        evaluateGuidanceFunctions( bislipSystems, vehicleSystems, currentTrajectoryPhase );

        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", 0.0 );

    }
    else if ( ( std::fmod( currentTime - startingEpoch, fixedStepSize * samplingRatio ) == 0.0 ) && ( currentTime > startingEpoch ) )
    {
        evaluateGuidanceFunctions( bislipSystems, vehicleSystems, currentTrajectoryPhase );

    }
    else if ( ( std::fmod( currentTime - startingEpoch, fixedStepSize * samplingRatio ) != 0.0 ) && ( ( currentTime > startingEpoch + fixedStepSize * samplingRatio ) ) )
    {

        currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack();
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getCurrentBodyFlapAngle() )  );

        currentBankAngle_ = bislipSystems->getCurrentBankAngle();

        currentAngleOfSideslip_ = 0.0;
        //std::cout << "Third one: " << currentTime - startingEpoch<< " | "<< currentTrajectoryPhase <<" | "<< tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack_ ) << " | "<< tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) <<std::endl;
        // bankLimit = tudat::unit_conversions::convertRadiansToDegrees( bislip::Variables::computeSkipSuppressionLimit( bodyMap_, vehicleName_, centralBodyName_ ) );
        //   std::cout << "currentBankAngle_: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
        //std::cout << "bankLimit: " <<  bankLimit << std::endl;

    }


} // updateGuidance


void MyGuidance::evaluateGuidanceFunctions(
        std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
        std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::string &currentTrajectoryPhase )
{
    bislipSystems->setCurrentAngleOfAttack( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::AngleOfAttack, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setEvaluatedBankAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::BankAngle, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setCurrentThrustElevationAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setCurrentThrustAzimuthAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustAzimuthAngle, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setCurrentThrottleSetting( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrottleSetting, bodyMap_, vehicleName_ ) );
    bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName_ )->getBodyMass(), bodyMap_.at( vehicleName_ )->getVehicleSystems()->getDryMass() ) );
    bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName_ ) );
    bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap_, vehicleName_) );

    currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack();

    bislipSystems->setCurrentBodyFlapAngle( bislip::Variables::computeBodyFlapDeflection( bodyMap_, vehicleName_ ) );

    vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getCurrentBodyFlapAngle() ) );

    double bankAngleForSkipSuppression = 0.0;
    bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );

    if ( ( currentTrajectoryPhase == "Ascent" ) ) { currentBankAngle_ = 0.0; }
    else
    {
        currentBankAngle_ = bislip::Variables::computeBankAngle( bodyMap_, vehicleName_ );

        bankAngleForSkipSuppression = double( bislip::Variables::determineBankAngleSign( currentBankAngle_ ) ) * std::abs( bislip::Variables::computeSkipSuppressionLimit( bodyMap_, vehicleName_, centralBodyName_ ) );
       //  std::cout << "Current (old) bank angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentBankAngle( ) ) << "  |  "<< "bankLimit: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankLimit ) << "  |  "<< "Calculated (new) Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
        bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );


        if ( std::isnan( bankAngleForSkipSuppression ) == false )
        {
            //  std::cout << "bankLimit: " <<  bankLimit << std::endl;
            //std::cout << "currentBankAngle_: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
            if ( std::abs( bankAngleForSkipSuppression ) > std::abs( currentBankAngle_ ) ){ currentBankAngle_ = bankAngleForSkipSuppression; }
         //   std::cout << "NEW currentBankAngle_: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
        }
    }
    bislipSystems->setCurrentBankAngle( currentBankAngle_ );

    // std::cout << "currentThrustMagnitude: " <<  bislipSystems->getCurrentThrustMagnitude() << std::endl;
    // std::cout << "currentBankAngle_: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
    // std::cout << "bankLimit: " <<  bankLimit << std::endl;

    currentAngleOfSideslip_ = 0.0;

}

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();

    return bislipSystems->getCurrentBodyFixedThrustDirection();
}

bool MyGuidance::getCurrentEngineStatus( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();

    return bislipSystems->getCurrentEngineStatus();
}

double MyGuidance::getCurrentThrustMagnitude( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();

    return bislipSystems->getCurrentThrustMagnitude();
}



double MyGuidance::getBodyFlapDeflection( )
{
    return bislip::Variables::computeBodyFlapDeflection( bodyMap_, vehicleName_ );
}

double MyGuidance::getBankAngle( )
{
    return bislip::Variables::computeBankAngle( bodyMap_, vehicleName_ );
}

Eigen::Vector6d MyGuidance::getPartialCurrentCoefficients( )
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

    if( bislipSystems_ == nullptr )
    {
        bislipSystems_ = std::dynamic_pointer_cast< bislip::BislipVehicleSystems >(
                    bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
    }

    return bislip::Variables::computePartialCurrentCoefficients( bodyMap_, vehicleName_ );
}




} // namespace bislip
