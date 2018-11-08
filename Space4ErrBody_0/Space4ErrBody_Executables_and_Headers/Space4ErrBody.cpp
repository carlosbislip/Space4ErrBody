/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved.
 * qqwqwqw
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 *    Notes
 *
 */

#include <iostream>
//#include <iomanip>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <cstring>
//#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
//#include <sstream>
#include <utility>// std::pair, std::get
#include <boost/filesystem/operations.hpp>

//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include <boost/random/variate_generator.hpp>

#include <Eigen/Core>
//#include <boost/bind.hpp>
//#include <boost/make_shared.hpp>
//#include <boost/shared_ptr.hpp>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"



//! Do I need a costumized one for HORUS? Probably.... Maybe not. I believe it
//! would be sufficient with the aerodynamic coefficient tables.
//#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>

//! Mine
#include "Space4ErrBody.h"
#include "applicationOutput_tudat.h"
#include "getStuff.h"
#include "updateGuidance.h"
//#include "updateGuidance_val.h"
#include "StopOrNot.h"
//#include "getAngularDistance.h"

//#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
//#include <Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h>
//#include <Tudat/Astrodynamics/Ephemerides/ephemeris.h>
//#include <Tudat/JsonInterface/Propagation/variable.h>
//#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
//#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
//#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>
//#include <Tudat/Astrodynamics/Ephemerides/frameManager.h>

Space4ErrBody::Space4ErrBody(
        const std::vector< std::vector< double > > &bounds,
        const std::string &problem_name,
        const std::string &vehicle_name,
        const std::vector< std::string > &parameterList,
        const std::vector< double > &vehicleParameterValues,
        const std::vector< std::string > &aeroCoeffFileList,
        const std::vector< double > &simulation_settingsValues,
        const std::vector< double > &initialConditionsValues,
        const std::vector< double > &terminationConditionsValues,
        const std::vector< double > &output_settingsValues,
        const std::string &outputSubFolder ):
    problemBounds_( bounds ),
    problem_name_( problem_name ),
    vehicle_name_( vehicle_name ),
    parameterList_( parameterList ),
    vehicleParameterValues_( vehicleParameterValues ),
    aeroCoeffFileList_( aeroCoeffFileList ),
    simulation_settingsValues_( simulation_settingsValues ),
    initialConditionsValues_( initialConditionsValues ),
    terminationConditionsValues_( terminationConditionsValues ),
    output_settingsValues_( output_settingsValues ),
    outputSubFolder_( outputSubFolder ){ }

//        const bool useExtendedDynamics) :
//    useExtendedDynamics_( useExtendedDynamics ){ }

//! Descriptive name of the problem
std::string Space4ErrBody::get_name() const
{
    return problem_name_;
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > Space4ErrBody::get_bounds() const
{
    return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
}




//! Implementation of the fitness function
//! For this ballistic case, the fitness function returns delta, which is a
//! vector containing differences in GOAL and calculated latitudes and
//! longitudes, their "norm", and time of flight.
//!
//! This function seems to be able to get quite long. It is possible to create a
//! series of functions that could create a more compact display? How feasible
//! would that be?
std::vector<double> Space4ErrBody::fitness( const std::vector< double > &x )  const
{

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat::ephemerides;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::mathematical_constants;
    using namespace tudat::input_output;
    using namespace tudat::unit_conversions;
    using namespace tudat::reference_frames;
    using namespace tudat;
    using namespace tudat_applications;
    using namespace tudat::aerodynamics;
    using namespace bislip;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  std::cout << "Unpacking data" << std::endl;

    //! Set simulation start epoch.
    const double simulationStartEpoch = simulation_settingsValues_[ 0 ]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Set simulation end epoch.
    // const double simulationEndEpoch = simulationStartEpoch + 18000.0; // 5 hours
    const double simulationEndEpoch = simulationStartEpoch + simulation_settingsValues_[ 1 ];

    //! Set numerical integration fixed step size.
    const double fixedStepSize = simulation_settingsValues_[ 2 ];

    //! Set number of control nodes.
    const int nodes = int( simulation_settingsValues_[ 3 ] );

    //! Assign Reference area
    const double S_ref = vehicleParameterValues_[0]; //m^2

    //! Assign Reference length: wingspan
    const double b_ref = vehicleParameterValues_[0]; //m^2

    //! Assign Reference length: chord
    const double c_ref = vehicleParameterValues_[0]; //m^2

    //! Assign position vector of moment reference center
    const Eigen::Vector3d R_mrc( vehicleParameterValues_[4], vehicleParameterValues_[5], vehicleParameterValues_[6] ); //m^2

    //! Assign position vector of center of mass
    const Eigen::Vector3d R_com( vehicleParameterValues_[7], vehicleParameterValues_[8], vehicleParameterValues_[9] ); //m^2

    //! Assign position vector of center of thrust
    const Eigen::Vector3d R_cot( vehicleParameterValues_[10], vehicleParameterValues_[11], vehicleParameterValues_[12] ); //m^2

    //! Assign initial mass
    const double M_i = vehicleParameterValues_[13]; // kg

    //! Assign propellant mass
    const double M_p = vehicleParameterValues_[14]; // kg

    //! Assign specific impulse
    const double Isp = vehicleParameterValues_[15]; // kg

    //! Assign initial height
    const double h_i = initialConditionsValues_[2]; // kg

    //! Assign initial Mach
    const double Mach_i = initialConditionsValues_[3]; // kg

    //! Assign initial Mach
    const double gamma_i_deg = initialConditionsValues_[4]; // kg

    //! Assign initial coordinates
    const double lat_i_deg = initialConditionsValues_[ 0 ];
    const double lon_i_deg = initialConditionsValues_[ 1 ];

    //! Assign final coordinates and additional termination conditions
    const double lat_f_deg = terminationConditionsValues_[ 0 ];
    const double lon_f_deg = terminationConditionsValues_[ 1 ];
    const double d_f_deg = terminationConditionsValues_[ 2 ];
    const double h_UP = terminationConditionsValues_[ 3 ];
    const double h_DN = terminationConditionsValues_[ 4 ];
    const double n_max = terminationConditionsValues_[ 5 ];
    const double q_dot_max = terminationConditionsValues_[ 6 ];
    const double q_d_max = terminationConditionsValues_[ 7 ];

    //! Convert angles from degrees to radians
    const double lat_i_rad = unit_conversions::convertDegreesToRadians( lat_i_deg );
    const double lon_i_rad = unit_conversions::convertDegreesToRadians( lon_i_deg );
    const double lat_f_rad = unit_conversions::convertDegreesToRadians( lat_f_deg );
    const double lon_f_rad = unit_conversions::convertDegreesToRadians( lon_f_deg );
    const double gamma_i_rad = unit_conversions::convertDegreesToRadians( gamma_i_deg );

    //! Declare and allocate vectors of interest.
    Eigen::VectorXd xn_interval( nodes - 1 );
    Eigen::VectorXd xn( nodes );
    Eigen::VectorXd alpha_deg( nodes );
    Eigen::VectorXd eps_T_deg( nodes );
    Eigen::VectorXd throttle( nodes );
    Eigen::VectorXd alpha_rad( nodes );
    Eigen::VectorXd eps_T_rad( nodes );
    Eigen::VectorXd E( nodes );
    Eigen::VectorXd E_hat( nodes );
    Eigen::VectorXd E_mapped( nodes );

    //! Various parameters
    const double R_E = 6.378137e6;
    const double mu = 3.986004418e14;
    const double omega_E = 7.292115*1E-5;
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;

    Eigen::Vector2d gravs = getGravs ( mu, J2, J3, J4, R_E, R_E, 0.0000001 );
    double rho, a, g0;
    g0 = gravs.norm();

    //! Re-allocate desicion vector values into workable vectors.
    for( int i = 0; i < nodes; i++ )
    {
        if ( i < nodes )
        {
            xn_interval( i ) = x[ i ];
        }
        alpha_deg( i )   = x[ i + nodes-1 ];
        eps_T_deg( i )   = x[ i + 2*nodes-1 ];
        throttle( i )    = x[ i + 3*nodes-1 ];

        alpha_rad( i ) = unit_conversions::convertDegreesToRadians( alpha_deg( i ) );
        eps_T_rad( i ) = unit_conversions::convertDegreesToRadians( eps_T_deg( i ) );
    }

    //! Create vector of node locations
    xn( 0 ) = 0;
    for( int i = 1; i < nodes + 1; i++ )
    {
        xn( i )        = xn( i - 1 ) + xn_interval( i - 1 );
    }

    /*
    std::cout << "x =  " << std::endl;
    for( int i = 0; i < int(x.size()); i++ )
    {
        std::cout << i << "   " << x[i] << std::endl;
    }
    std::cout << "-------" << std::endl;
    std::cout << "xn_interval =  " << xn_interval << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "xn =  " << xn << std::endl;
    std::cout << "-------" << std::endl;

    std::cout << "alpha_deg =  " << alpha_deg << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "eps_T_deg =  " << eps_T_deg << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "throttle =  " << throttle << std::endl;
    std::cout << "-------" << std::endl;
*/

    //! Impose constraints on first and last energy nodes
    a = 301.7;//NRLMSISE00Atmosphere::getSpeedOfSound( R_E + height( 0 ), 0, 0, simulationStartEpoch );
    double V_i = a * Mach_i;
    double V_f = 0.99 * sqrt( mu / ( R_E + h_UP ) );
    double E_min = g0 * h_i + 0.5 * V_i * V_i;
    double E_max = g0 * h_UP + 0.5 * V_f * V_f;

    //! Normalize fist and last Energy nodes
    double E_hat_min = E_min/E_max;
    double E_hat_max = E_max/E_max;

    //! Map energy levels to control node locations
    E_mapped =  ( E_hat_max - E_hat_min ) * xn.array() + E_hat_min;

    //! Associate decision vector values to mapped energy levels within data maps
    std::map< double, double > map_alpha_deg, map_eps_T_deg, map_throttle, map_alpha_rad, map_eps_T_rad;
    for ( unsigned int i = 0; i < E_mapped.size( ); ++i )
    {
        map_alpha_deg[ E_mapped( i ) ] = alpha_deg( i );
        map_eps_T_deg[ E_mapped( i ) ] = eps_T_deg( i );
        map_throttle[ E_mapped( i ) ] = throttle( i );
        map_alpha_rad[ E_mapped( i ) ] = alpha_rad( i );
        map_eps_T_rad[ E_mapped( i ) ] = eps_T_rad( i );
    }

    //! Create interpolator
    std::shared_ptr< InterpolatorSettings > interpolatorSettings =
            std::make_shared< tudat::interpolators::InterpolatorSettings >( cubic_spline_interpolator );
    std::shared_ptr< OneDimensionalInterpolator< double, double > > interpolator_alpha_deg,interpolator_eps_T_deg,interpolator_throttle,interpolator_alpha_rad,interpolator_eps_T_rad;

    interpolator_alpha_deg = createOneDimensionalInterpolator< double, double >( map_alpha_deg, interpolatorSettings );
    interpolator_eps_T_deg = createOneDimensionalInterpolator< double, double >( map_eps_T_deg, interpolatorSettings );
    interpolator_throttle  = createOneDimensionalInterpolator< double, double >( map_throttle,  interpolatorSettings );
    interpolator_alpha_rad = createOneDimensionalInterpolator< double, double >( map_alpha_rad, interpolatorSettings );
    interpolator_eps_T_rad = createOneDimensionalInterpolator< double, double >( map_eps_T_rad, interpolatorSettings );

    double lat_f_deg_calc = lat_i_deg;
    double lon_f_deg_calc = lon_i_deg;
    double d_i_deg        = unit_conversions::convertRadiansToDegrees( getAngularDistance(lat_i_rad,lon_i_rad,lat_f_rad,lon_f_rad) );
    double d_f_deg_calc   = d_i_deg;
    //std::cout << d_f_deg_calc << std::endl;
    double h_UP_calc       = h_i;
    double tof = simulation_settingsValues_[ 1 ];
    std::string simulation_file_name_suffix;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //   std::cout << "Creating environment" << std::endl;

    //! Assign central body name.
    const std::string centralBodyName = "Earth";

    //! Declare/define simulation body settings data map.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { centralBodyName },
                                    simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );

    //! Define atmospheric model.
    bodySettings[ centralBodyName ]->atmosphereSettings = std::make_shared< AtmosphereSettings >(
                nrlmsise00 );

    //! Define ephemeris model settings.
    //! This is an acceptable 'cheat' were Earth is placed at the barycenter.
    //! Use only when there arent any third body perturbations (Moon, Sun, etc.)
    bodySettings[ centralBodyName ]->ephemerisSettings =
            std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    //! Reset ephemeris to J2000.
    bodySettings[ centralBodyName ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    //! Define Earth's radius. Using spice here. Is there a way to get a 'radius
    //! field'? Im interested in includig Earth's flattening, yet am unsure how
    //! to properly do it.
    const double radius_Earth = spice_interface::getAverageRadius( centralBodyName );
    double radius_Earth_i = radius_Earth;
    double radius_Earth_UP = radius_Earth;

    //!---------------------------- Still playing around with this section

    //! https://en.wikipedia.org/wiki/Flattening#Numerical_values_for_planets
    //! https://www.vcalc.com/wiki/vCalc/WGS-84+Earth+flattening+factor
    double flattening_Earth = 1 / 1;//nput_data_.back();
    if (flattening_Earth != 1)
    {
        bodySettings[ centralBodyName ]->shapeModelSettings = std::make_shared< OblateSphericalBodyShapeSettings >( radius_Earth, flattening_Earth );

        // Declare variable in which raw result is to be put by Spice function.
        double radii[ 3 ];

        // Call Spice function to retrieve gravitational parameter.
        SpiceInt numberOfReturnedParameters;
        bodvrd_c( "Earth", "RADII", 3, &numberOfReturnedParameters, radii );
        radii[0] = radii [0] * 1e3 ;
        radii[1] = radii [1] * 1e3 ;
        radii[1] = radii [2] * 1e3 ;
        // Compute average and convert from km to m.
        //return unit_conversions::convertKilometersToMeters< double >(
        //           radii[ 0 ] + radii[ 1 ] + radii[ 2 ] ) / 3.0;

        const double radius_Equator =  6378137.0;
        const double radius_pole = 6356752.3142;
        const double polar_angle_i_deg = 90 - lat_i_deg;
        const double polar_angle_f_deg = 90 - lat_f_deg;
        const double polar_angle_i_rad = polar_angle_i_deg * mathematical_constants::PI / 180;
        const double polar_angle_f_rad = polar_angle_f_deg * mathematical_constants::PI / 180;

        const double radius_Earth_i_x = radii[ 0 ] * std::cos(polar_angle_i_rad) * std::cos(lon_i_rad);
        const double radius_Earth_i_y = radii[ 1 ] * std::cos(polar_angle_i_rad) * std::sin(lon_i_rad);
        const double radius_Earth_i_z = radii[ 2 ] * std::sin(polar_angle_i_rad);
        radius_Earth_i = std::sqrt(pow(radius_Earth_i_x,2) + pow(radius_Earth_i_y,2) + pow(radius_Earth_i_z,2));
        const double radius_Earth_UP_x = radii[ 0 ] * std::cos(polar_angle_f_rad) * std::cos(lon_f_rad);
        const double radius_Earth_UP_y = radii[ 1 ] * std::cos(polar_angle_f_rad) * std::sin(lon_f_rad);
        const double radius_Earth_UP_z = radii[ 2 ] * std::sin(polar_angle_f_rad);
        radius_Earth_UP = std::sqrt(pow(radius_Earth_UP_x,2) + pow(radius_Earth_UP_y,2) + pow(radius_Earth_UP_z,2));

    }
    //!--------------------------------------------------------------

    //! Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    //! Create Earth's rotational ephemeris.
    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( centralBodyName )->getRotationalEphemeris( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Assign coefficient file names.
    //const std::string dragCoefficients = aeroCoeffFileList_[0];
    //const std::string liftCoefficients = aeroCoeffFileList_[1];
    // const std::string dragControlSurfaceForceCoefficients = "HORUS_CD_CS.txt";
    // const std::string liftControlSurfaceForceCoefficients = "HORUS_CL_CS.txt";
    //const std::string momentCoefficients = aeroCoeffFileList_[2];
    // const std::string controlSurfaceMomentCoefficients = "HORUS_CM_CS.txt";


    ///////// Start: Vehicle Aerodynamics Section

    //! Define physical meaning of independent variables
    //!     Angle of attack
    //!     Mach number
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );

    //! Define list of files for force coefficients. Entry denotes direction.
    //!     0 : x-direction (C ~D~/C ~X~)
    //!     1 : y-direction (C ~S~/C ~Y~)
    //!     2 : z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles;
    // std::map< int, std::string > controlSurfaceForceCoefficientFiles;
    std::map< int, std::string > momentCoefficientFiles;
    // std::map< int, std::string > controlSurfaceMomentCoefficientFiles;

    forceCoefficientFiles[ 0 ] = aeroCoeffFileList_[0]; // Set drag coefficient file
    forceCoefficientFiles[ 2 ] = aeroCoeffFileList_[1]; // Set lift coefficient file
    // controlSurfaceForceCoefficientFiles[ 0 ] = dragControlSurfaceForceCoefficients; // Set drag coefficient increment file
    // controlSurfaceForceCoefficientFiles[ 2 ] = liftControlSurfaceForceCoefficients; // Set lift coefficient increment file
    momentCoefficientFiles[ 1 ] = aeroCoeffFileList_[2]; // Set pitch moment coefficient file
    // controlSurfaceMomentCoefficientFiles[ 1 ] = controlSurfaceMomentCoefficients; // Set pitch moment coefficient increment file

    //       std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > controlSurfaceInterface =
    //              std::make_shared< CustomControlSurfaceIncrementAerodynamicInterface >(
    //                   &dummyControlIncrements,
    //                   boost::assign::list_of( angle_of_attack_dependent )( control_surface_deflection_dependent ) );
    //       std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
    //       controlSurfaceList[ "bodyflap" ] = controlSurfaceInterface;

    //! Define reference frame in which the loaded coefficients are defined.
    //! Have to get some more background info here to properly understand it.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;


    //! Load and parse coefficient files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles,
                momentCoefficientFiles,
                b_ref,
                S_ref,
                c_ref,
                R_mrc,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame,
                areCoefficientsInNegativeAxisDirection );

    //! Create shared pointer for aerodynamic coefficient increments.
    // std::shared_ptr< system_models::VehicleSystems > systemsModels = std::make_shared< system_models::VehicleSystems >( );

    //bodyMap.at( vehicle_name_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ) );

    ///////// End: Vehicle Aerodynamics Section

    //! Create vehicle objects.
    bodyMap[ vehicle_name_ ] = std::make_shared< simulation_setup::Body >( );

    //! Set body Mass.
    bodyMap[ vehicle_name_ ]->setConstantBodyMass( M_i );

    //! Set vehicle aerodynamic coefficients.
    bodyMap[ vehicle_name_ ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface(
                    aerodynamicCoefficientSettings,
                    vehicle_name_ ) );

    //! Set vehicle system models for aerodynamic coefficient increments.
    //bodyMap[ vehicle_name_ ]->setVehicleSystems( systemsModels );
    //bodyMap[ vehicle_name_ ]->getAerodynamicCoefficientInterface( )->setControlSurfaceIncrements( controlSurfaceList );

    //! Finalize body creation. Not entirely sure what this does. I believe it
    //! may place the body (Earth) at the barycenter. Probably another hack to
    //! facilitate bringing in everything into the inertial frame.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    //! Set spherical elements for vehicle's initial state.
    //! ARTBITRARILY CHOSEN to be in Earth-Fixed frame.
    Eigen::Vector6d EntryState_spherical;
    EntryState_spherical( SphericalOrbitalStateElementIndices::radiusIndex )       = radius_Earth_i + h_i;
    EntryState_spherical( SphericalOrbitalStateElementIndices::latitudeIndex )     = lat_i_rad;
    EntryState_spherical( SphericalOrbitalStateElementIndices::longitudeIndex )    = lon_i_rad;
    EntryState_spherical( SphericalOrbitalStateElementIndices::speedIndex )        = V_i;
    EntryState_spherical( SphericalOrbitalStateElementIndices::flightPathIndex )   = gamma_i_rad;
    EntryState_spherical( SphericalOrbitalStateElementIndices::headingAngleIndex ) = getHeadingToTarget( lat_i_rad,lon_i_rad,lat_f_rad,lon_f_rad );

    //std::cout << "radius_Earth_i + h_i:  " << radius_Earth_i + h_i << std::endl;
    //std::cout << "V_i:  " << V_i << std::endl;
    //std::cout << "gamma_i_rad:  " << gamma_i_rad << std::endl;

    //! Two things are being done here
    //!     Converting state vector from spherical to Cartesian elements
    //!     Transforming state vector from Earth-Fixed frame to Inertial frame.
    const Eigen::Vector6d systemInitialState = transformStateToGlobalFrame(
                convertSphericalOrbitalToCartesianState( EntryState_spherical ),
                simulationStartEpoch,
                earthRotationalEphemeris );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    //! Declare acceleration data map.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > vehicleAccelerations;

    //! Define gravitational model. The central body acts this force on itself.
    //! Arbitrary maximum degree/order. Equivalent functionality to Cartesian with corresponding maximum
    //! degree/order.
    vehicleAccelerations[ centralBodyName ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    //! Define aerodynamic accelerations. The atmosphere of the central body acts this force on the vehicle.
    vehicleAccelerations[ centralBodyName ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

    std::shared_ptr< MyGuidance > ThrustGuidance = std::make_shared< MyGuidance >(
                bodyMap,
                vehicle_name_,
                E_max,
                interpolator_alpha_deg,
                interpolator_eps_T_deg,
                interpolator_throttle );

    //double thrustMagnitude = 250000.0;
    //double specificImpulse = 500.0;

    std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< ThrustDirectionGuidanceSettings >(thrust_direction_from_existing_body_orientation, "Earth" );
    
    std::function< double( const double ) > thrustMagnitudeFunction =
            std::bind( &MyGuidance::getCurrentThrustMagnitude, ThrustGuidance );

    std::function< double( const double ) > specificImpulseFunction =
            std::bind( &MyGuidance::getCurrentSpecificImpulse, ThrustGuidance );

    std::function< bool( const double ) > isEngineOnFunction =
            std::bind( &MyGuidance::getCurrentEngineStatus, ThrustGuidance );

    std::function< Eigen::Vector3d( ) > BodyFixedThrustDirection =
            std::bind( &MyGuidance::getCurrentBodyFixedThrustDirection, ThrustGuidance );

    std::function< void( const double ) > customThrustResetFunction =
            std::bind( &MyGuidance::updateGuidance, ThrustGuidance, std::placeholders::_1 );

    std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction,
                specificImpulseFunction,
                isEngineOnFunction,
                BodyFixedThrustDirection,
                customThrustResetFunction );
    // std::make_shared< ConstantThrustMagnitudeSettings >( thrustMagnitude, specificImpulse );

    //! Define thrust acceleration settings. The vehicle acts this force on itself.
    vehicleAccelerations[ vehicle_name_ ].push_back( std::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );

    //! Declare acceleration settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    //! Assign acceleration map from the vehicleAccelerations data map.
    accelerationMap[ vehicle_name_ ] = vehicleAccelerations;

    //! Define bodies that will be propagated. Only 1.
    bodiesToPropagate.push_back( vehicle_name_ );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( centralBodyName );

    //! Set acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationMap,
                bodiesToPropagate,
                centralBodies );

    //! Aerodynamic guidance is set AFTER the accelerations and BEFORE propagating.
    //! Declare and assign aerodynamic guidance functions.
    std::shared_ptr< MyGuidance > AeroGuidance = std::make_shared< MyGuidance >(
                bodyMap,
                vehicle_name_,
                E_max,
                interpolator_alpha_deg,
                interpolator_eps_T_deg,
                interpolator_throttle);

    //! Set Guidance angle functions.
    setGuidanceAnglesFunctions( AeroGuidance, bodyMap.at( vehicle_name_ ) );

    //! Define constant orientation
    //double constantAngleOfAttack = unit_conversions::convertDegreesToRadians( 10 );
    //double constantBankAngle = unit_conversions::convertDegreesToRadians( 0 );
    //double constantSideSlipeAngle = unit_conversions::convertDegreesToRadians( 0 );
    //bodyMap.at( vehicle_name_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            [ = ]( ){ return constantAngleOfAttack; }, [ = ]( ){ return constantSideSlipeAngle; }, [ = ]( ){ return constantBankAngle; } );
    //std::cout << "Creating vehicle: Guidance is set" << std::endl;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE TERMINATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /* std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
             std::make_shared< SingleDependentVariableSaveSettings >(
                 altitude_dependent_variable,
                 vehicle_name_,
                 centralBodyName );

    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationDependentVariableTerminationSettings >(
                terminationDependentVariable,
                h_UP,
                true );
*/

    //In case your custom function requires more inputs (e.g., it may depend on the position of the spacecraft or other variables that are not the current time), you can use boost::bind to add more inputs.
    //As an example, the case where the state of the spacecraft is added as an input is shown below:
    // std::function< Eigen::Vector6d( ) > HORUS_StateFunction =
    //       std::bind( &Body::getState, bodyMap.at( vehicle_name_ ) );

    //! Define vector with termination conditions.
    std::vector< double > term_cond;
    term_cond.push_back( d_f_deg );
    term_cond.push_back( h_UP );
    term_cond.push_back( h_DN );

    std::vector< double > additional_data;
    //! Set initial coordinates. Earth-Fixed.
    additional_data.push_back( lat_i_rad );
    additional_data.push_back( lon_i_rad );

    //! Set target coordinates. Earth-Fixed.
    additional_data.push_back( lat_f_rad );
    additional_data.push_back( lon_f_rad );

    //! Set initial distance to target.
    additional_data.push_back( getAngularDistance( lat_i_rad,lon_i_rad,lat_f_rad,lon_f_rad) );

    //! Define termination settings.
    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationCustomTerminationSettings >(
                boost::bind( &bislip::StopOrNot, bodyMap, vehicle_name_, term_cond, additional_data ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          CREATE LIST OF DEPENDENT VARIABLES        ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  std::cout << "Creating propagation settings" << std::endl;

    //! Define list of dependent variables to save.
    //! The file that prints out the text saying what is saved has been modified
    //! to mute this message. During optimization this would clutter the
    //! terminal screen and is incredibly annoying.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dep_varList;
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_cartesian_position,
                    vehicle_name_,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_spherical_position,
                    vehicle_name_,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicle_name_,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic,
                    vehicle_name_,
                    centralBodyName,
                    1 ) );
    //dep_varList.push_back(
    //            std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
    //                vehicle_name_,
    //                latitude_angle,
    //                centralBodyName) );
    //dep_varList.push_back(
    //            std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
    //                vehicle_name_,
    //                longitude_angle,
    //                centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicle_name_,
                    heading_angle,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicle_name_,
                    flight_path_angle,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicle_name_,
                    angle_of_attack,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicle_name_,
                    angle_of_sideslip,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicle_name_,
                    bank_angle,
                    centralBodyName) );
    //dep_varList.push_back(
    //            std::make_shared< SingleDependentVariableSaveSettings >(
    //                total_mass_rate_dependent_variables,
    //                vehicle_name_ ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dep_varToSave =
            std::make_shared< DependentVariableSaveSettings >( dep_varList );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Create translational propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies,
              accelerationModelMap,
              bodiesToPropagate,
              systemInitialState,
              terminationSettings,
              cowell,
              dep_varToSave );


    ///////////////////////             CREATE MASS RATE SETTINGS              ////////////////////////////////////////////

    // Create mass rate models
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );
    //std::map< std::string, std::vector< std::shared_ptr< MassRateModelSettings > > > massRateModelSettings =

    std::map< std::string, std::shared_ptr< MassRateModel > > massRateModels;
    massRateModels[ vehicle_name_ ] = createMassRateModel(
                vehicle_name_, massRateModelSettings, bodyMap, accelerationModelMap );

    // Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( vehicle_name_ );

    //! Set initial mass of vehicle.
    Eigen::VectorXd initialBodyMasses( 1 );
    initialBodyMasses( 0 ) = M_i;

    //! Create mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses, terminationSettings, dep_varToSave );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // !Create list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    //! Create propagation settings for both mass and translational dynamics.
    //std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
    //        std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings, dep_varToSave );

    //! Create propagation settings for ONLY translational dynamics.
      std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            translationalPropagatorSettings;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Create integrator settings.
    std::shared_ptr< IntegratorSettings<  > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              fixedStepSize, 1, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::cout << "Starting propagation" << std::endl;

    //! Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< double > dynamicsSimulator(
                bodyMap,
                integratorSettings,
                propagatorSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        CLEAN UP AND CALCULATE FITNESS                //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Retrieve results
    //! Maybe this command would allow the use of certain results to guide
    //! further propagations within the same run. Sort of stitching them together.
    // std::map< double, Eigen::VectorXd > integrationResult =
    //        dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    //! Extract Final Epoch
    //! Trying to find the actual END time of the simulation via key extraction.
    //double simulationEndEpoch_calc = ( --dynamicsSimulator.getEquationsOfMotionNumericalSolution().end() )->first;

    //! Extracting the contents of the last timestep via extracted key
    //const Eigen::Vector6d systemFinalState = integrationResult[simulationEndEpoch_calc];

    //! Told by Dominic that this here gives the final epoch directly
    const double simulationEndEpoch_calc =  ( dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin() )->first;

    //! Told by Dominic that this here gives the final state directly
    //Eigen::Vector6d systemFinalState =  ( dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin() )->second;

    //! Calculate Time of Flight
    tof = simulationEndEpoch_calc - simulationStartEpoch;

    //! Transform final state from Inertial Frame to Earth-Fixed Frame
    //Eigen::Vector6d systemFinalState_EARTh_UPIXED = transformStateToTargetFrame( systemFinalState, simulationEndEpoch_calc, earthRotationalEphemeris );

    //! Transform expected final state to Inertial frame. This is where the
    //! destination is located at the termination of the simulation.
    //! Transformation is required because coordinates are Earth-fixed, similar
    //! to input values.
    //systemFinalStateGOAL = transformStateToGlobalFrame(systemFinalStateGOAL,simulationEndEpoch_calc,earthRotationalEphemeris );

    //! Calculate latitude and longitude of the GOAL state: Inertial Frame.
    //const double altitude_f_GOAL_calc = std::sqrt( pow(systemFinalStateGOAL[0],2) +
    //       pow(systemFinalStateGOAL[1],2) + pow(systemFinalStateGOAL[2],2) ) ;
    //double lon_f_rad_GOAL_calc = std::atan2(systemFinalStateGOAL[1] , systemFinalStateGOAL[0]);
    //double lat_f_rad_GOAL_calc = std::asin(systemFinalStateGOAL[2] / altitude_f_GOAL_calc) ;

    //! Convert lat/lon of GOAL state to degrees: Inertial Frame
    //const double lon_f_deg_GOAL_calc = lon_f_rad_GOAL_calc * 180 / mathematical_constants::PI;
    //const double lat_f_deg_GOAL_calc = lat_f_rad_GOAL_calc * 180 / mathematical_constants::PI;

    //! Calculate latitude and longitude of the final state: Earth-Fixed Frame.
    //! Dominic brought up the possibility of using:
    //!               body_fixed_relative_spherical_position
    //! as a dependent variable. If possible, convert this calulation to the
    //! extraction of the final entry. It generates radius, lat and lon.
    //!
    //!
    //!

    Eigen::VectorXd dep_var_FINAL_STATE = ( dynamicsSimulator.getDependentVariableHistory( ).rbegin() )->second;
    //std::cout << "dep_var_FINAL_STATE: " << dep_var_FINAL_STATE << std::endl;


    //const double altitude_f_calc = dep_var_FINAL_STATE[3];
    const double lat_f_rad_calc  = dep_var_FINAL_STATE[ 4 ];
    const double lon_f_rad_calc  = dep_var_FINAL_STATE[ 5 ];
    h_UP_calc                    = dep_var_FINAL_STATE[ 7 ];
    //const double altitude_f_calc = std::sqrt( pow(systemFinalState_EARTh_UPIXED[0],2) +
    //       pow(systemFinalState_EARTh_UPIXED[1],2) + pow(systemFinalState_EARTh_UPIXED[2],2) ) ;
    //const double lon_f_rad_calc = std::atan2(systemFinalState_EARTh_UPIXED[1] , systemFinalState_EARTh_UPIXED[0]);
    //const double lat_f_rad_calc = std::asin(systemFinalState_EARTh_UPIXED[2] / altitude_f_calc) ;

    //! Convert coordinates of final state to degrees: Earth-Fixed Frame
    lat_f_deg_calc = unit_conversions::convertRadiansToDegrees( lat_f_rad_calc );
    lon_f_deg_calc = unit_conversions::convertRadiansToDegrees( lon_f_rad_calc );

    //! Calculate angular distance of final state from target coordinates.
    const double d_rad = getAngularDistance( lat_f_rad_calc, lon_f_rad_calc, lat_f_rad, lon_f_rad );
    double d_deg = unit_conversions::convertRadiansToDegrees( d_rad );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Get time stamp for this specific simulation. This avoids overwriting the
    //! file if another individual with the same properties shows up in other
    //! evolutions.
    std::string simulation_save_time = getCurrentDateTime( false );

    //! Create unique filename that cannot be overwritten due to the timestamp.
    simulation_file_name_suffix =
            // std::to_string( Va_ordered( 0 ) ) + "_" +
            // std::to_string( gamma_rad_ordered( 0 ) ) + "_" +
            // std::to_string( chi_i_deg_ordered( 0 ) ) + "_" +
            std::to_string( tof ) + "_" +
            //std::to_string(x[3]) + "_" +
            std::to_string( lat_f_deg_calc ) + "_" +
            std::to_string( lon_f_deg_calc ) + "_" +
            simulation_save_time;

    std::string complete_file_name_Prop = "HORUSPropHistory_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_DepVar = "HORUSDepVar_" + simulation_file_name_suffix + ".dat";

    //! Will print out depending on some input values. Each entry corresponds to
    //! a different type of output. Entries are binary, 0 or 1.
    if ( output_settingsValues_[ 1 ] == 1 )//&& check == nodes )
    {

        //! Write HORUS propagation history to file.
        writeDataMapToTextFile( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                complete_file_name_Prop,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
    }
    if ( output_settingsValues_[ 2 ] == 1 )//&& check == nodes )
    {
        //! Write HORUS dependent variables' history to file.
        writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory( ),
                                complete_file_name_DepVar,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
    }


    //  }

    //! Calculate offset of final state from GOAL state: Earth-Fixed Frame
    //const double dif_lat_rad = lat_f_rad - lat_f_rad_calc;
    //const double dif_lon_rad = lon_f_rad - lon_f_rad_calc;

    //! Calculate offsets of final state in degrees: Earth-Fixed Frame
    const double dif_lat_deg = lat_f_deg - lat_f_deg_calc;
    const double dif_lon_deg = lon_f_deg - lon_f_deg_calc;

    //! Calculate "norm" of offsets. This is an arbitrary function I have
    //! implemented to pass on as an 'objective function'. It relates the
    //! differences such that when minimizing the offsets, there is an additional
    //! unsigned value that always goes to zero. Most definitely unsure about how
    //!  'proper' it is, yet is what works for the current BALLISTIC case.
    const double dif_norm = std::sqrt( ( dif_lat_deg * dif_lat_deg ) + ( dif_lon_deg * dif_lon_deg ) );

    //! Calculate offset of final angular distance to termination condition distance.
    const double dif_d_deg = d_deg - terminationConditionsValues_[ 2 ];
    //! Calculate offset from goal elevation.
    const double dif_h = h_UP - h_UP_calc;



    double dif_xn = xn( nodes - 1 ) - 1;
    double dif_E_mapped = E_mapped( nodes - 1 ) - 1;

    //! Assign values to Fitness vector! At the moment these are all 'objective
    //! functions'. No constraints have been implemented. To modify this I have
    //! change the header file and define how many are elements are objective
    //! function, equality contraints, and inequlity constraints. This vector
    //! here must contain them is that exact order: nOF, nEC, nIC.
    std::vector< double > delta;
    delta.push_back( dif_norm );
    delta.push_back( dif_lat_deg );
    delta.push_back( dif_lon_deg );
    delta.push_back( dif_d_deg );
    delta.push_back( dif_h );
    delta.push_back( tof );  // Not sure yet how this one affects the optimization. Included for completion.
    delta.push_back( dif_xn );
    delta.push_back( dif_E_mapped );

    //std::cout << "h_UP: " << h_UP << std::endl;
    //std::cout << "h_UP_calc: " << h_UP_calc << std::endl;

    //! Print results to terminal. Used to gauge progress.
    std::cout << std::fixed << std::setprecision(10) <<
                 std::setw(15) << "dif_xn = " <<
                 std::setw(14) << dif_xn <<
                 std::setw(15) << " dif_E_mapped = " <<
                 std::setw(16) <<  dif_E_mapped <<
                 std::setw(15) << " dif_norm = " <<
                 std::setw(16) << dif_norm <<
                 std::setw(15) << " dif_lat_deg = " <<
                 std::setw(16) << dif_lat_deg <<
                 std::setw(15) << " dif_lon_deg = " <<
                 std::setw(16) << dif_lon_deg <<
                 std::setw(13) << "dif_d_deg = " <<
                 std::setw(16) << dif_d_deg <<
                 std::setw(9) << "dif_h = " <<
                 std::setw(16) << dif_h <<
                 std::setw(7) << "tof = " <<
                 std::setw(16) << tof <<
                 std::setw(120) << simulation_file_name_suffix << std::endl;

    // std::cout <<dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) << std::endl;




    return delta;

} // Fitness function.
