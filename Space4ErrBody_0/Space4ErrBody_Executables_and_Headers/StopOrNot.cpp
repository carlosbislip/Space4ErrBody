#include "Space4ErrBody.h"
#include "StopOrNot.h"
#include "getStuff.h"
//##include "updateGuidance.h"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

namespace bislip {

bool StopOrNot( const tudat::simulation_setup::NamedBodyMap& bodyMap,
                const std::string &vehicleName,
                const std::vector< double > &vehicleParameterValues,
                const std::vector< double > &parameterBounds,
                const std::vector< double > &terminationConditionsValues,
                const std::vector< double > &additional_data,
                const double &E_max,
                const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_throttle) {

    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );
    //! Extract current latitude
    //! //const double lat_c_rad = FlightConditions_->getCurrentGeodeticLatitude( );
    const double lat_c_rad =  FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );

    //! Extract current longitude
    const double lon_c_rad =  FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );

    //! Extract initial coordinates
    const double lat_i_rad = additional_data[ 0 ];
    const double lon_i_rad = additional_data[ 1 ];


    //! Assign final coordinates and additional termination conditions
    const double lat_f_deg = terminationConditionsValues[ 0 ];
    const double lon_f_deg = terminationConditionsValues[ 1 ];
    const double d_f_deg = terminationConditionsValues[ 2 ];
    const double h_UP = terminationConditionsValues[ 3 ];
    const double h_DN = terminationConditionsValues[ 4 ];
    const double V_max = terminationConditionsValues[ 5 ];
    const double n_max = terminationConditionsValues[ 6 ];
    const double q_dot_max = terminationConditionsValues[ 7 ];
    const double q_dyn_max = terminationConditionsValues[ 8 ];




    //! Extract target coordinates
    const double lat_f_rad = tudat::unit_conversions::convertDegreesToRadians( lat_f_deg );;
    const double lon_f_rad = tudat::unit_conversions::convertDegreesToRadians( lon_f_deg );

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
    const double initial_d_to_target = bislip::getAngularDistance( lat_i_rad,lon_i_rad,lat_f_rad,lon_f_rad);

    //! Extract current altitude
    const double current_height = FlightConditions->getCurrentAltitude();

    //! Extract current airspeed
    const double current_V = FlightConditions->getCurrentAirspeed();

    //! Extract current angle of attack
    const double current_AoA = FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );

    //! Extract current flight-path angle
    const double current_gamma = FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );

    //! Exract current data
    const double currentMass = bodyMap.at( vehicleName )->getBodyMass( );

    //! Calculate current normalized mass-normalized energy.
    const double g0 = 9.80665;
    const double E_hat = ( g0 * current_height + 0.5 * current_V * current_V ) / E_max;

    //! Evaluate current throttle setting and thrust elevation angle.
    const double throttle = interpolator_throttle->interpolate( E_hat );
   // const double eps_T_deg = interpolator_eps_T_deg->interpolate( E_hat );

    //const double AoA_LB = tudat::unit_conversions::convertDegreesToRadians( parameterBounds[ 2 ] );
   // const double AoA_UB = tudat::unit_conversions::convertDegreesToRadians( parameterBounds[ 3 ] );
   // const double eps_T_deg_LB = parameterBounds[ 4 ];
   // const double eps_T_deg_UB = parameterBounds[ 5 ];

    const double finalMass = vehicleParameterValues[ 13 ];

    const double current_rho = FlightConditions->getCurrentDensity( );
    const double current_q_dot = FlightConditions->getCurrentAerodynamicHeatRate( );
    const double current_q_dot_calc = 0.5 * current_V * current_V * current_V * current_rho;
    const double current_q_dyn = FlightConditions->getCurrentDynamicPressure( );
    const double current_q_dyn_calc = 0.5 * current_V * current_V * current_rho;

    /*std::vector< double > term_cond;
    term_cond.push_back( d_f_deg );
    term_cond.push_back( h_UP );
    term_cond.push_back( h_DN );
    */

    //! Declare boolean that is returns to terminate the simulation.
    bool done = false;

    //! Currently 3 conditions could terminate a simulation
    //!     When distance to target is less than a specified number
    //!     When altitude is more than a specified number
    //!     When altitude is less than a specified number
    //!     When angular distance traveled is larger than angular distance among endpoints
    if ( ( d_togo_rad <= tudat::unit_conversions::convertDegreesToRadians( d_f_deg )  ) || ( current_height <= h_DN ) || ( current_height >= h_UP ) || ( d_traveled_rad >= initial_d_to_target ) )
    {
        done = true;
        std::cout << "d_togo_rad:  " << d_togo_rad << std::endl;
        std::cout << "d_f_deg :  " << d_f_deg  << std::endl;
        std::cout << "current_height:  " << current_height << std::endl;
        std::cout << "h_DN :  " << h_DN  << std::endl;
        std::cout << "h_UP :  " << h_UP  << std::endl;
        std::cout << "d_traveled_rad:  " << d_traveled_rad << std::endl;
        std::cout << "initial_d_to_target:  " << initial_d_to_target << std::endl;
        std::cout << "Stopping propagation" << std::endl;

    }
    /*else if  ( current_height <= terminationConditionsValues_[ 2 ] )
    {
        done = true;
        std::cout << "current_height:  " << current_height << std::endl;
        std::cout << "terminationConditionsValues_[ 2 ] :  " << terminationConditionsValues_[ 2 ]  << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }*/
    else if ( current_q_dot > q_dot_max*1 )
    {
        done = true;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - additional_data[ 2 ]<< std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "Expected thrust:  " << throttle * additional_data[ 3 ]<< std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "current_V:  " << current_V << std::endl;
        std::cout << "current_rho:  " << current_rho << std::endl;
        std::cout << "current_q_dot:  " << current_q_dot << std::endl;
        std::cout << "current_q_dot_calc:  " << current_q_dot_calc << std::endl;
        std::cout << "q_dot_max:  " << q_dot_max << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;

    }
    else if ( current_q_dyn > q_dyn_max*1 )
    {
        done = true;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - additional_data[ 2 ]<< std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "Expected thrust:  " << throttle * additional_data[ 3 ]<< std::endl;

        std::cout << "current_V:  " << current_V << std::endl;
        std::cout << "current_rho:  " << current_rho << std::endl;
        std::cout << "current_q_dyn:  " << current_q_dyn << std::endl;
        std::cout << "current_q_dyn_calc:  " << current_q_dyn_calc << std::endl;
        std::cout << "q_dyn_max:  " << q_dyn_max << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    else if ( currentMass < finalMass )
    {
        done = true;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - additional_data[ 2 ]<< std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "Expected thrust:  " << throttle * additional_data[ 3 ]<< std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    else if ( current_V * std::sin( current_gamma ) < 0 )
    {
        done = true;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - additional_data[ 2 ]<< std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "Expected thrust:  " << throttle * additional_data[ 3 ]<< std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "current_V * std::sin( current_gamma ):  " << current_V * std::sin( current_gamma ) << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    /*else if ( ( throttle < 0 ) || ( throttle < parameterBounds[ 6 ] ) || ( throttle > parameterBounds[ 7 ] ) )
    {
        //! This condition terminates the simulation if thrust goes negative or beyond engine capabilities.

        done = true;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    else if ( ( ( throttle >= 0 ) && ( throttle <= 1.0 ) )
              && ( ( ( current_gamma < 0 ) || ( current_gamma > tudat::unit_conversions::convertDegreesToRadians( 90 ) ) )
                   || ( ( current_AoA < AoA_LB ) || ( current_AoA > AoA_UB ) )
                   ||  ( ( eps_T_deg < eps_T_deg_LB ) || ( eps_T_deg > eps_T_deg_UB ) ) ) )
    */
    else if ( ( ( throttle >= 0 ) && ( throttle <= 1.0 ) )
              &&  ( ( current_gamma < 0 ) || ( current_gamma > tudat::mathematical_constants::PI / 2.0 ) ) )
    {
        //! This condition terminates the simulation if flight-path angle is negative  or beyond 90 deg while the engine is on.
        //if ( ( current_gamma < 0 ) || ( current_gamma > tudat::unit_conversions::convertDegreesToRadians( 90 ) ) )
        // {
        done = true;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "current_gamma:  " << tudat::unit_conversions::convertRadiansToDegrees( current_gamma ) << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    else if ( E_hat > 1.0 )
    {
        done = true;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }

    /*// else if ( ( current_AoA < tudat::unit_conversions::convertDegreesToRadians( parameterBounds[ 2 ] ) ) || ( current_AoA > tudat::unit_conversions::convertDegreesToRadians( parameterBounds[ 3 ] ) ) )
       // {
            //! This condition terminates the simulation if angle of attack is beyond the specified bounds while the engine is on.
       //     done = true;
       //     std::cout << "throttle:  " << throttle << std::endl;
            std::cout << "current_AoA:  " << tudat::unit_conversions::convertRadiansToDegrees( current_AoA ) << std::endl;
        //    std::cout << "Stopping propagation" << std::endl;
        //}
        //else if ( ( eps_T_deg < parameterBounds[ 4 ] ) || ( eps_T_deg > parameterBounds[ 5 ] ) )
       // {
            //! This condition terminates the simulation if thrust elevation angle is beyond the specified bounds while the engine is on.
        //    done = true;
        //    std::cout << "throttle:  " << throttle << std::endl;
            std::cout << "eps_T_deg:  " << eps_T_deg << std::endl;
            std::cout << "Stopping propagation" << std::endl;
        }

    //}*/
    else
    {
        done = false;
    }




    return done;
}
}
