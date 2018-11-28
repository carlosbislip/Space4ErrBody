#include "updateGuidance.h"
#include "getStuff.h"

namespace bislip  {

double MyGuidance::setE_hat(
        const double &height,
        const double &airspeed,
        const double &E_max )
{
    return  ( 9.80665 * height + 0.5 * airspeed * airspeed ) / E_max;
}

} // namespace bislip
