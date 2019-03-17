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
        //std::cout << "First one" << currentTime << std::endl;
        currentAngleOfAttack_ = getAngleofAttack();
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", 0.0 );
        currentBankAngle_ = 0.0;
        currentAngleOfSideslip_ = 0.0;
       std::cout << "First one: " << currentTime- startingEpoch << " | "<< currentTrajectoryPhase <<" | "<< currentAngleOfAttack_ << " | "<< currentBankAngle_ <<std::endl;

    }
    else if ( ( std::fmod( currentTime - startingEpoch, fixedStepSize * samplingRatio ) == 0.0 ) && ( currentTime > startingEpoch ) )
{
        currentAngleOfAttack_ = getAngleofAttack();
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", tudat::unit_conversions::convertDegreesToRadians( getBodyFlapDeflection( ) ) );
        if ( ( currentTrajectoryPhase == "Ascent" ) ) { currentBankAngle_ = 0.0; }
        else { currentBankAngle_ = getBankAngle(); }

        std::cout << "Second one: " << currentTime- startingEpoch << " | "<< currentTrajectoryPhase <<" | "<< currentAngleOfAttack_ << " | "<< currentBankAngle_ <<std::endl;

        //std::cout << currentTime << " | "<< currentTrajectoryPhase <<" | "<< currentAngleOfAttack_ << currentBankAngle_ <<std::endl;

        //double maxBankAngle = tudat::unit_conversions::convertDegreesToRadians(
        //          bislip::Variables::computeEquilibriumGlideLimit( bodyMap_, vehicleName_, centralBodyName_ ) );

        //if ( std::isnan( maxBankAngle ) == false ){ if ( currentBankAngle_ > maxBankAngle ){ currentBankAngle_ = maxBankAngle; } }

        // }

        currentAngleOfSideslip_ = 0.0;

    }
    else if ( ( std::fmod( currentTime - startingEpoch, fixedStepSize * samplingRatio ) != 0.0 ) && ( ( currentTime > startingEpoch + fixedStepSize * samplingRatio ) ) )
    {
        // std::cout << "Third one:  " << currentTime << std::endl;

       // if( FlightConditions_ == nullptr )
       // {
            FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                        bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
       // }

        currentAngleOfAttack_ = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", tudat::unit_conversions::convertDegreesToRadians( getBodyFlapDeflection( ) ) );
        currentBankAngle_ = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::bank_angle );
        currentAngleOfSideslip_ = 0.0;
        std::cout << "Third one: " << currentTime - startingEpoch<< " | "<< currentTrajectoryPhase <<" | "<< currentAngleOfAttack_ << " | "<< currentBankAngle_ <<std::endl;


    }


} // updateGuidance

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName_ );
}

bool MyGuidance::getCurrentEngineStatus( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName_ )->getBodyMass(), bodyMap_.at( vehicleName_ )->getVehicleSystems()->getDryMass() );
}

double MyGuidance::getCurrentThrustMagnitude( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::computeThrustMagnitude( bodyMap_, vehicleName_ );
}

double MyGuidance::getAngleofAttack( )
{
    return bislip::Variables::computeAngleofAttack( bodyMap_, vehicleName_ );
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
