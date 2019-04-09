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

        if ( bislipSystems->getTempBankAngle() * currentBankAngle_ < 0.0 ) { bislipSystems->setBankAngleReversalTrigger( true ); }
        else { bislipSystems->setBankAngleReversalTrigger( false ); }


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

    }


} // updateGuidance


void MyGuidance::evaluateGuidanceFunctions(
        std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
        std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::string &currentTrajectoryPhase )
{
    //! Determine and set values by evaluating interpolators.
   // std::cout << "Determine and set values by evaluating interpolators." << std::endl;
    bislipSystems->setCurrentAngleOfAttack( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::AngleOfAttack, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setEvaluatedBankAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::BankAngle, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setCurrentThrustElevationAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setCurrentThrustAzimuthAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustAzimuthAngle, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setCurrentThrottleSetting( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrottleSetting, bodyMap_, vehicleName_ ) );

    //! Determine and set values by evaluating various functions.
    //std::cout << "Determine and set values by evaluating various functions." << std::endl;
    bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName_ )->getBodyMass(), bodyMap_.at( vehicleName_ )->getVehicleSystems()->getDryMass() ) );
    bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName_ ) );
    bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap_, vehicleName_) );

    //! Bound and assign the current angle of attack.
    //std::cout << "Bound and assign the current angle of attack." << std::endl;
    double angleOfAttackUpperBound = ( bislipSystems->getAlphaMachEnvelopeUBInterpolator() )->interpolate( FlightConditions_->getCurrentMachNumber() ) ;
    //std::cout << "Angle of attack Upper Bound: " << angleOfAttackUpperBound << std::endl;
    double angleOfAttackLowerBound = ( bislipSystems->getAlphaMachEnvelopeLBInterpolator() )->interpolate( FlightConditions_->getCurrentMachNumber() ) ;
    //std::cout << "Angle of attack Lower Bound: " << angleOfAttackLowerBound << std::endl;
    if ( bislipSystems->getCurrentAngleOfAttack() > angleOfAttackUpperBound )
    {
        bislipSystems->setCurrentAngleOfAttack( angleOfAttackUpperBound );
        currentAngleOfAttack_ = angleOfAttackUpperBound;
    }
    else if ( bislipSystems->getCurrentAngleOfAttack() < angleOfAttackLowerBound )
    {
        bislipSystems->setCurrentAngleOfAttack( angleOfAttackLowerBound );
        currentAngleOfAttack_ = angleOfAttackLowerBound;
    }
    else
    { currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack(); }


    //! Compute and assign the current bodyflap angle
    //std::cout << "Compute and assign the current bodyflap angle." << std::endl;
    bislipSystems->setCurrentBodyFlapAngle( bislip::Variables::computeBodyFlapDeflection( bodyMap_, vehicleName_ ) );
    vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getCurrentBodyFlapAngle() ) );


    double bankAngleForSkipSuppression = 0.0;
    bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );


    currentBankAngle_ = 0.0;
    if ( currentTrajectoryPhase != "Ascent" )
    {
        if ( FlightConditions_->getCurrentDynamicPressure() > bislipSystems->getMinimumDynamicPressureforControlSurface() )
        {
            currentBankAngle_ = bislip::Variables::computeBankAngle( bodyMap_, vehicleName_ );


            bankAngleForSkipSuppression = double( bislip::Variables::determineBankAngleSign( currentBankAngle_ ) ) * std::abs( bislip::Variables::computeSkipSuppressionLimit( bodyMap_, vehicleName_, centralBodyName_ ) );
            //  std::cout << "Current (old) bank angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentBankAngle( ) ) << "  |  "<< "bankLimit: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankLimit ) << "  |  "<< "Calculated (new) Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
            bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );

            if ( std::isnan( bankAngleForSkipSuppression ) == false )
            {
                //  std::cout << "bankLimit: " <<  bankLimit << std::endl;
                //std::cout << "currentBankAngle_: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
               // if ( std::abs( bankAngleForSkipSuppression ) > std::abs( currentBankAngle_ ) ){ currentBankAngle_ = bankAngleForSkipSuppression; }
                //   std::cout << "NEW currentBankAngle_: " <<  tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) << std::endl;
            }
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
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    return bislipSystems->getCurrentBodyFixedThrustDirection();
}

bool MyGuidance::getCurrentEngineStatus( )
{
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    return bislipSystems->getCurrentEngineStatus();
}

double MyGuidance::getCurrentThrustMagnitude( )
{
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    return bislipSystems->getCurrentThrustMagnitude();
}

/*
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
*/



} // namespace bislip
