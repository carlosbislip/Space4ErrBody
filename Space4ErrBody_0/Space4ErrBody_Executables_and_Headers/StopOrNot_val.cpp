#include "Space4ErrBody.h"
#include "StopOrNot_val.h"
#include "getStuff.h"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>


bool StopOrNot( const tudat::simulation_setup::NamedBodyMap& bodyMap,
                const std::string &vehicleName,
                const std::vector< double > term_cond ) {

    //! Extract current latitude
    //! //const double lat_c_rad = FlightConditions_->getCurrentGeodeticLatitude( );
    const double lat_c_rad = bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );

    //! Extract current longitude
    const double lon_c_rad = bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );

    //! Extract initial coordinates
    const double lat_i_rad = bodyMap.at( vehicleName )->getInitialLat();
    const double lon_i_rad = bodyMap.at( vehicleName )->getInitialLon();

    //! Extract target coordinates
    const double lat_f_rad = bodyMap.at( vehicleName )->getTargetLat();
    const double lon_f_rad = bodyMap.at( vehicleName )->getTargetLon();

    //////////////// Calc Distance to target
    //! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
    //! Maybe consider Vicenty's formulation: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    //const double a = std::sin( (lat_f_rad_ - lat_c_rad) / 2) * std::sin( (lat_f_rad_ - lat_c_rad) / 2) + std::cos( lat_c_rad ) * std::cos( lon_c_rad ) * std::sin( (lon_f_rad_ - lon_c_rad) / 2) * std::sin( (lon_f_rad_ - lon_c_rad) / 2);
    //const double d_rad = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
    const double d_togo_rad = bislip::getAngularDistance( lat_c_rad , lon_c_rad , lat_f_rad , lon_f_rad );
    //std::acos( std::sin(lat_c_rad) * std::sin(lat_f_rad) + std::cos(lat_c_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_c_rad) );
    //const double d_togo_deg = unit_conversions::convertRadiansToDegrees( d_togo_rad );

    //! Calculate total distance traveled
    const double d_traveled_rad = bislip::getAngularDistance( lat_i_rad , lon_i_rad , lat_c_rad , lon_c_rad );
    //std::acos( std::sin(lat_c_rad) * std::sin(lat_f_rad) + std::cos(lat_c_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_c_rad) );
    //const double d_traveled_deg = unit_conversions::convertRadiansToDegrees( d_traveled_rad );

    //! Extract Initial distance to target
    const double initial_d_to_target = bodyMap.at( vehicleName )->getInitialDistanceToTarget();

    //! Extract current altitude
    const double current_altitude = bodyMap.at( vehicleName )->getFlightConditions( )->getCurrentAltitude();


    //! Declare boolean that is returns to terminate the simulation.
    bool done;

    //! Currently 3 conditions could temrinate a simulation
    //!     When distance to target is less than a specified number
    //!     When altitude is less than a specified number
    //!     When angular distance traveled is larger than angular distance among endpoints
    if ( ( d_togo_rad <= term_cond[0] ) || ( current_altitude <= term_cond[1] ) || ( d_traveled_rad >= initial_d_to_target ) )
    {
        done = true;
    }
    else
    {
        done = false;
    }
    
return done;
}
