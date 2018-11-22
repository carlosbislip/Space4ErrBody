#include "updateGuidance.h"

#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

namespace bislip {

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection(  )
{
    double eps_T = interpolator_eps_T_deg_->interpolate( E_hat_ );

    if ( eps_T < parameterBounds_[ 4 ] )
    {
        eps_T = parameterBounds_[ 4 ];
    }
    if ( eps_T > parameterBounds_[ 5 ] )
    {
        eps_T = parameterBounds_[ 5 ];
    }

    currentEpsilon_ = tudat::unit_conversions::convertDegreesToRadians( eps_T );

    //! Simplified expressions becuase thrust azimuth is known to be zero. I.e. phi_T = 0
    bodyFixedThrustDirection_( 0 ) = std::cos( currentEpsilon_ );
    bodyFixedThrustDirection_( 1 ) = 0.0;
    bodyFixedThrustDirection_( 2 ) = std::sin( currentEpsilon_ );


    return bodyFixedThrustDirection_;
}
} // namespace aerodynamics
//} // namespace tudat
