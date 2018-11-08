#include "Space4ErrBody.h"
#include "StopOrNot.h"
#include "getStuff.h"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>

namespace bislip {

bool StopOrNot( const tudat::simulation_setup::NamedBodyMap& bodyMap,
                const std::string &vehicleName,
                const std::vector< double > term_cond,
                const std::vector< double > additional_data ) {

    //! Extract current latitude
    //! //const double lat_c_rad = FlightConditions_->getCurrentGeodeticLatitude( );
    const double lat_c_rad = bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );

    //! Extract current longitude
    const double lon_c_rad = bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );

    //! Extract initial coordinates
    const double lat_i_rad = additional_data[ 0 ];
    const double lon_i_rad = additional_data[ 1 ];

    //! Extract target coordinates
    const double lat_f_rad = additional_data[ 2 ];
    const double lon_f_rad = additional_data[ 3 ];

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
    const double initial_d_to_target = additional_data[ 4 ];

    //! Extract current altitude
    const double current_altitude = bodyMap.at( vehicleName )->getFlightConditions( )->getCurrentAltitude();


    //! Declare boolean that is returns to terminate the simulation.
    bool done;

    //! Currently 3 conditions could temrinate a simulation
    //!     When distance to target is less than a specified number
    //!     When altitude is more than a specified number
    //!     When altitude is less than a specified number
    //!     When angular distance traveled is larger than angular distance among endpoints
    if ( ( d_togo_rad <= tudat::unit_conversions::convertDegreesToRadians( term_cond[ 0 ] )  ) || ( current_altitude >= term_cond[ 1 ] ) || ( d_traveled_rad >= initial_d_to_target ) )
    {
        done = true;
        //std::cout << "d_togo_rad:  " << d_togo_rad << std::endl;
        //std::cout << "term_cond[ 0 ] :  " << term_cond[ 0 ]  << std::endl;
        //std::cout << "current_altitude:  " << current_altitude << std::endl;
        //std::cout << "term_cond[ 1 ] :  " << term_cond[ 1 ]  << std::endl;
        //std::cout << "d_traveled_rad:  " << d_traveled_rad << std::endl;
        //std::cout << "initial_d_to_target:  " << initial_d_to_target << std::endl;
        //std::cout << "Stopping propagation" << std::endl;

    }
    else if  ( current_altitude <= term_cond[ 2 ] )
    {
        done = true;
        //std::cout << "current_altitude:  " << current_altitude << std::endl;
        //std::cout << "term_cond[ 2 ] :  " << term_cond[ 2 ]  << std::endl;
        //std::cout << "Stopping propagation" << std::endl;

    }
    else
    {
        done = false;
    }
    
return done;
}
}
