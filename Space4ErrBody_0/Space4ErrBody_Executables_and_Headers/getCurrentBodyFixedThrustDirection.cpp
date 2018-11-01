#include "updateGuidance.h"

#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

namespace bislip {

Eigen::Vector3d MyThrustGuidance::getCurrentBodyFixedThrustDirection(  )
{
    return bodyFixedThrustDirection_;
}
} // namespace aerodynamics
//} // namespace tudat
