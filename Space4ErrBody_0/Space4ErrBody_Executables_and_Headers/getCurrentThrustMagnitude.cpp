#include "updateGuidance.h"



namespace bislip {

double MyGuidance::getCurrentThrustMagnitude( )
{
    double throttle = interpolator_throttle_->interpolate( E_hat_ );

   if ( throttle < 0.0 )
    {
        throttle = 0.0;
    }
   if ( throttle > 1.0 )
    {
        throttle = 1.0;
    }

    currentThrustMagnitude_ = throttle * maxThrust_;

    return currentThrustMagnitude_;
}
} // namespace aerodynamics
//} // namespace tudat
