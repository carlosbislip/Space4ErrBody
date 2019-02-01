#include "updateGuidance.h"

namespace bislip {

Eigen::Vector6d MyGuidance::getCurrentCoefficients( )
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
        bislipSystems_ = std::dynamic_pointer_cast< bislip::VehicleSystems >(
                    bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
    }

    return bislip::Variables::computeCurrentCoefficients( FlightConditions_, coefficientInterface_, bislipSystems_ );
}

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::computeBodyFixedThrustDirection( FlightConditions_, bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
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

    return bislip::Variables::computeThrustMagnitude( FlightConditions_, bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
}

double MyGuidance::getBodyflapCmIncrement( )
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
        bislipSystems_ = std::dynamic_pointer_cast< bislip::VehicleSystems >(
                    bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
    }

    return bislip::Variables::computeBodyflapCmIncrement( FlightConditions_, coefficientInterface_, bislipSystems_  );
}




} // namespace bislip
