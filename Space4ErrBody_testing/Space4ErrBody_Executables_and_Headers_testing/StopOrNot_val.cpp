#include "Space4ErrBody.h"
#include "StopOrNot_val.h"


#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
namespace tudat
{

bool StopOrNot( const simulation_setup::NamedBodyMap& bodyMap, const std::string vehicleName ) {
    //! Extract current latitude
    //! //const double lat_c_rad = FlightConditions_->getCurrentGeodeticLatitude( );
    const double lat_c_rad = bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::latitude_angle );

    //! Extract current longitude
    const double lon_c_rad = bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::longitude_angle );

    //! Extract target coordinates
    const double lat_f_rad = bodyMap.at( vehicleName )->getTargetLat();
    const double lon_f_rad = bodyMap.at( vehicleName )->getTargetLon();

    //////////////// Calc Distance to target
    //! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
    //! Maybe consider Vicenty's formulation: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    //const double a = std::sin( (lat_f_rad_ - lat_c_rad) / 2) * std::sin( (lat_f_rad_ - lat_c_rad) / 2) + std::cos( lat_c_rad ) * std::cos( lon_c_rad ) * std::sin( (lon_f_rad_ - lon_c_rad) / 2) * std::sin( (lon_f_rad_ - lon_c_rad) / 2);
    //const double d_rad = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
    const double d_rad =  std::acos( std::sin(lat_c_rad) * std::sin(lat_f_rad) + std::cos(lat_c_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_c_rad) );
    const double d_deg = unit_conversions::convertRadiansToDegrees( d_rad );

    bool done;

    if ( d_deg <= 0.75)
    {
        done = true;
    }
    else
    {
        done = false;
    }

    //! How can I incorporate when the distance increases after hitting the initial mark?


return done;
}

} // namespace tudat
