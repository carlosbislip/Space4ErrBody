#include "updateGuidance.h"

#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

namespace bislip {

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection(  )
{
    return bodyFixedThrustDirection_;
}
} // namespace aerodynamics
//} // namespace tudat
