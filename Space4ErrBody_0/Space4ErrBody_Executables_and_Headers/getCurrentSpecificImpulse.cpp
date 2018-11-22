#include "updateGuidance.h"



namespace bislip {

double MyGuidance::getCurrentSpecificImpulse( )
{
    currentSpecificImpulse_ = Isp_;

    return currentSpecificImpulse_;
}
} // namespace aerodynamics
//} // namespace tudat
