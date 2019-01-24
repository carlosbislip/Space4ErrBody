#include "updateGuidance.h"


namespace bislip {

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection(  )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::variables::computeBodyFixedThrustDirection( FlightConditions_, bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
}
} // namespace bislip
