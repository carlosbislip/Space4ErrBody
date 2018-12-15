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
#include <chrono>
#include <thread>

#include <Eigen/Core>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>
#include <Tudat/Mathematics/Interpolators/interpolator.h>
#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>

//! Mine
#include "Space4ErrBody.h"
#include "applicationOutput_tudat.h"
#include "getStuff.h"
#include "updateGuidance.h"
#include "bislipVariables.h"
//#include "updateGuidance_val.h"
#include "StopOrNot.h"

Space4ErrBody::Space4ErrBody(
        const std::vector< std::vector< double > > &bounds,
        const std::string &problem_name,
        const std::string &vehicle_name,
        const std::vector< std::string > &parameterList_Ascent,
        const std::vector< double > &parameterBounds_Ascent,
        const std::vector< std::string > &parameterList_Descent,
        const std::vector< double > &parameterBounds_Descent,
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
    parameterList_Ascent_( parameterList_Ascent ),
    parameterBounds_Ascent_( parameterBounds_Ascent ),
    parameterList_Descent_( parameterList_Descent ),
    parameterBounds_Descent_( parameterBounds_Descent ),
    vehicleParameterValues_( vehicleParameterValues ),
    aeroCoeffFileList_( aeroCoeffFileList ),
    simulation_settingsValues_( simulation_settingsValues ),
    initialConditionsValues_( initialConditionsValues ),
    terminationConditionsValues_( terminationConditionsValues ),
    output_settingsValues_( output_settingsValues ),
    outputSubFolder_( outputSubFolder ){ }

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
    //using namespace std::this_thread; // sleep_for, sleep_until
    //using namespace std::chrono; // nanoseconds, system_clock, seconds

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Unpacking data" << std::endl;

    //! Declare and initialize simulation start epoch.
    const double simulationStartEpoch = simulation_settingsValues_[ 0 ]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Declare and initialize simulation end epoch.
    const double simulationEndEpoch = simulationStartEpoch + simulation_settingsValues_[ 1 ];

    //! Declare and initialize numerical integration fixed step size.
    const double fixedStepSize = simulation_settingsValues_[ 2 ];

    //! Declare and initialize number of control nodes.
    const unsigned long nodes_Ascent = simulation_settingsValues_[ 3 ];
    const unsigned long nodes_Descent = simulation_settingsValues_[ 4 ];

    //! Declare and initialize Reference area
    const double S_ref = vehicleParameterValues_[ 0 ]; // m^2

    //! Declare and initialize Reference length: wingspan
    const double b_ref = vehicleParameterValues_[ 1 ]; // m^2

    //! Declare and initialize Reference length: chord
    const double c_ref = vehicleParameterValues_[ 2 ]; // m^2

    //! Declare and initialize position vector of moment reference center
    const Eigen::Vector3d R_mrc( vehicleParameterValues_[ 3 ], vehicleParameterValues_[ 4 ], vehicleParameterValues_[ 5 ] ); // m

    //! Declare and initialize position vector of center of mass
    const Eigen::Vector3d R_com( vehicleParameterValues_[ 6 ], vehicleParameterValues_[ 7 ], vehicleParameterValues_[ 8 ] ); // m

    //! Declare and initialize position vector of center of thrust
    const Eigen::Vector3d R_cot( vehicleParameterValues_[ 9 ], vehicleParameterValues_[ 10 ], vehicleParameterValues_[ 11 ] ); // m

    //! Declare and initialize initial mass
    const double initialMass = vehicleParameterValues_[ 12 ]; // kg

    //! Declare and initialize landing mass
    const double landingMass = vehicleParameterValues_[ 13 ]; // kg

    //! Declare and initialize specific impulse
    const double specificImpulse = vehicleParameterValues_[ 14 ]; // kg

    //! Declare and initialize maximum engine thrust
    const double maxThrust = vehicleParameterValues_[ 15 ]; // N

    //! Declare and initialize starting height
    const double h_i = initialConditionsValues_[ 2 ]; // m

    //! Declare and initialize initial Mach
    const double Mach_i = initialConditionsValues_[ 3 ]; // -

    //! Declare and initialize initial flight-path angle
    const double gamma_i_deg = initialConditionsValues_[ 4 ]; // deg

    //! Declare and initialize starting position coordinates.
    const double initialLat_deg = initialConditionsValues_[ 0 ];
    const double initialLon_deg = initialConditionsValues_[ 1 ];

    //! Declare and initialize final position coordinates and additional termination conditions
    const double targetLat_deg = terminationConditionsValues_[ 0 ];
    const double targetLon_deg = terminationConditionsValues_[ 1 ];

    //! Declare and initialize various termination conditions
    const double finalDistanceToTarget_deg = terminationConditionsValues_[ 2 ];
    //const double h_UP = terminationConditionsValues_[ 3 ];
    const double h_DN = terminationConditionsValues_[ 4 ];
    //const double V_UP = terminationConditionsValues_[ 5 ];
    const double n_max = terminationConditionsValues_[ 6 ];
    const double q_dot_max = terminationConditionsValues_[ 7 ];
    const double q_dyn_max = terminationConditionsValues_[ 8 ];

    //! Still working on these
    const double R_N = 3.0;
    const double lambda = unit_conversions::convertDegreesToRadians( 30.0 );
    const double epsilon = 0.7;
    const double x_T = 1E-9;
    const double phi = unit_conversions::convertDegreesToRadians( 3.0 );

    //! Convert angles from degrees to radians
    const double initialLat_rad = unit_conversions::convertDegreesToRadians( initialLat_deg );
    const double initialLon_rad = unit_conversions::convertDegreesToRadians( initialLon_deg );
    const double targetLat_rad = unit_conversions::convertDegreesToRadians( targetLat_deg );
    const double targetLon_rad = unit_conversions::convertDegreesToRadians( targetLon_deg );
    const double initialFlightPathAngle_rad= unit_conversions::convertDegreesToRadians( gamma_i_deg );

    //! Pre-define various variables used to determine fitness.
    double targetLat_deg_calc = initialLat_deg;
    double targetLon_deg_calc = initialLon_deg;
    double initialDistanceToTarget_rad = bislip::variables::computeAngularDistance( initialLat_rad, initialLon_rad, targetLat_rad, targetLon_rad );
    double initialDistanceToTarget_deg = unit_conversions::convertRadiansToDegrees( initialDistanceToTarget_rad );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //   std::cout << "Creating environment" << std::endl;

    //! Declare and initialize central body name.
    const std::string centralBodyName = "Earth";

    //! Declare and initialize simulation body settings data map.
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
    if ( flattening_Earth != 1 )
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
       /* const double polar_angle_i_deg = 90 - lat_i_deg;
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
*/
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

    //std::cout << "Creating vehicle" << std::endl;

    //! Assign coefficient file names.
    //const std::string dragCoefficients = aeroCoeffFileList_[0];
    //const std::string liftCoefficients = aeroCoeffFileList_[1];
    const std::string dragControlSurfaceForceCoefficients = "HORUS_CD_CS.txt";
    const std::string liftControlSurfaceForceCoefficients = "HORUS_CL_CS.txt";
    //const std::string momentCoefficients = aeroCoeffFileList_[2];
    const std::string controlSurfaceMomentCoefficients = "HORUS_CM_CS.txt";

    ///////// Start: Vehicle Aerodynamics Section

    //! Define physical meaning of independent variables
    //!     Angle of attack
    //!     Mach number
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );

    // Define physical meaning of independent variables for control surface increments, in this case Mach number, angle of attack and control surface deflection
    //std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > controlSurfaceIndependentVariableNames;
    //controlSurfaceIndependentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    //controlSurfaceIndependentVariableNames.push_back( aerodynamics::mach_number_dependent );
    //controlSurfaceIndependentVariableNames.push_back( aerodynamics::control_surface_deflection_dependent );


    //! Define list of files for force coefficients. Entry denotes direction.
    //!     0 : x-direction (C ~D~/C ~X~)
    //!     1 : y-direction (C ~S~/C ~Y~)
    //!     2 : z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles;
    //std::map< int, std::string > controlSurfaceForceCoefficientFiles;
    std::map< int, std::string > momentCoefficientFiles;
    //std::map< int, std::string > controlSurfaceMomentCoefficientFiles;

    std::string ELEVON_L = "ElevonLeft";
    std::string ELEVON_R = "ElevonRight";
    std::string BODYFLAP = "BodyFlap";


    forceCoefficientFiles[ 0 ] = aeroCoeffFileList_[0]; // Set drag coefficient file
    forceCoefficientFiles[ 2 ] = aeroCoeffFileList_[1]; // Set lift coefficient file

    momentCoefficientFiles[ 1 ] = aeroCoeffFileList_[2]; // Set pitch moment coefficient file
    // controlSurfaceForceCoefficientFiles[ 0 ] = dragControlSurfaceForceCoefficients; // Set drag coefficient increment file
    //controlSurfaceForceCoefficientFiles[ 2 ] = liftControlSurfaceForceCoefficients; // Set lift coefficient increment file
    //controlSurfaceMomentCoefficientFiles[ 1 ] = controlSurfaceMomentCoefficients; // Set pitch moment coefficient increment file

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

    // Add settings for control surface increments to main aerodynamic coefficients
    // aerodynamicCoefficientSettings->setControlSurfaceSettings(
    //             readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
    //                 controlSurfaceForceCoefficientFiles,
    //                 controlSurfaceMomentCoefficientFiles,
    //                 controlSurfaceIndependentVariableNames, ELEVON_L ) );



    //bodyMap.at( vehicle_name_ )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ) );

    ///////// End: Vehicle Aerodynamics Section

    //! Create vehicle objects.
    bodyMap[ vehicle_name_ ] = std::make_shared< simulation_setup::Body >( );

    //! Set body Mass.
    bodyMap[ vehicle_name_ ]->setConstantBodyMass( initialMass );

    //! Create vehicle systems and initialize by setting the landing (dry) Mass.
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared<system_models::VehicleSystems >( landingMass );

    //!
    bodyMap.at( vehicle_name_ )->setVehicleSystems( vehicleSystems );

    //! Set nose radius.
    vehicleSystems->setNoseRadius( R_N );

    //! Set wing sweep angle.
    vehicleSystems->setWingSweepAngle( lambda );

    //! Set wall emmisivity.
    vehicleSystems->setWallEmissivity( epsilon );

    //! Set local body angle.
    vehicleSystems->setLocalBodyAngle( phi );

    //! Set transition length.
    vehicleSystems->setTransitionDistance( x_T );


    //! Set vehicle aerodynamic coefficients.
    bodyMap[ vehicle_name_ ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface(
                    aerodynamicCoefficientSettings,
                    vehicle_name_ ) );

    //std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
    //controlSurfaceList[ ELEVON_L ] = controlSurfaceInterface;

    //! Set vehicle system models for aerodynamic coefficient increments.
    //bodyMap[ vehicle_name_ ]->setVehicleSystems( systemsModels );
    //bodyMap[ vehicle_name_ ]->getAerodynamicCoefficientInterface( )->setControlSurfaceIncrements( controlSurfaceList );

    //! Finalize body creation. Not entirely sure what this does. I believe it
    //! may place the body (Earth) at the barycenter. Probably another hack to
    //! facilitate bringing in everything into the inertial frame.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    bodyMap[ vehicle_name_ ]->setEphemeris(
                std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ), centralBodyName ) );
    /*
     *Thr problem with the mass propagation can be solved (-ish), by turning the throw.... of this error into a std::cerr (printing the error, but not terminating the program). Alternatively, it should also work if you use:
     *
     *   bodyMap[ "Spaceplane" ]->setEphemeris(
     *               std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
     *                   std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ),   "Earth" ) );
     *
     *What I think is going on is the following:
     * - When setting up the propagation, the code checks if all required environment models are available.
     * - The mass propagation requires the current state of the spacecraft
     * - The current state of the spacecraft requires that either the spacecraft has an ephemeris, or this state is propagated.
     * - When propagating the mass and state, the state is propagated. However, when making this check for the mass propagation, the code does not know this, resulting in this error.
     *
     *Either one of the fixes I give above should work. Let me know if it works out,
     */

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE NODAL STRUCTURE             /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // std::cout << "Creating nodal structure" << std::endl;

    //! Declare and allocate vectors of interest for Ascent phase.
    Eigen::VectorXd xn_interval_Ascent( nodes_Ascent - 1 );
    Eigen::VectorXd xn_Ascent( nodes_Ascent );
    Eigen::VectorXd alpha_deg_Ascent( nodes_Ascent );
    Eigen::VectorXd sigma_deg_Ascent( nodes_Ascent );
    Eigen::VectorXd eps_T_deg_Ascent( nodes_Ascent );
    Eigen::VectorXd phi_T_deg_Ascent( nodes_Ascent );
    Eigen::VectorXd throttle_Ascent( nodes_Ascent );

    //! Declare and allocate vectors of interest for Descent phase.
    Eigen::VectorXd xn_interval_Descent( nodes_Descent - 1 );
    Eigen::VectorXd xn_Descent( nodes_Descent );
    Eigen::VectorXd alpha_deg_Descent( nodes_Descent );
    Eigen::VectorXd sigma_deg_Descent( nodes_Descent );
    Eigen::VectorXd eps_T_deg_Descent( nodes_Descent );
    Eigen::VectorXd phi_T_deg_Descent( nodes_Descent );
    Eigen::VectorXd throttle_Descent( nodes_Descent );

    //! Various parameters
    const double R_E = radius_Earth;
    const double omega_E = 7.292115*1E-5;
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;

    //Eigen::Vector2d gravs = getGravs ( mu, J2, J3, J4, R_E, R_E, 0.0000001 );
    double g0 = tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION;

    //std::cout << "Re-allocate Ascent DVs" << std::endl;
    //! Re-allocate decision vector values into workable vectors for Ascent phase.
    for( unsigned int i = 0; i < nodes_Ascent; i++ )
    {
        if ( i < ( nodes_Ascent - 1 ) ){ xn_interval_Ascent( i ) = x[ i ]; }
        alpha_deg_Ascent( i ) = x[ i + 1*nodes_Ascent - 1 ];
        sigma_deg_Ascent( i ) = x[ i + 2*nodes_Ascent - 1 ];
        eps_T_deg_Ascent( i ) = x[ i + 3*nodes_Ascent - 1 ];
        phi_T_deg_Ascent( i ) = x[ i + 4*nodes_Ascent - 1 ];
        throttle_Ascent( i )  = x[ i + 5*nodes_Ascent - 1 ];
    }

    //! Declare and initialize number of parameters exclusive to Ascent phase.
    const unsigned long  N = ( parameterList_Ascent_.size() - 3 ) * nodes_Ascent - 1;

    //! Declare and initialize various parameters common to the entire trajectory.
    //!     Initial velocity.
    const double V_i = x[ N ];
    //!     Velocity at the maximum height.
    const double V_UP = x[ N + 1 ];
    //!     Maximum height.
    const double h_UP = x[ N + 2 ];


    // std::cout << "Re-allocate Descent DVs" << std::endl;
    //! Re-allocate decision vector values into workable vectors for Descent phase.
    for( unsigned int i = 0; i < nodes_Descent; i++ )
    {
        if ( i < ( nodes_Descent - 1) ){ xn_interval_Descent( i ) = x[ ( N + 3 ) + i ]; }
        alpha_deg_Descent( i ) = x[ ( N + 3 ) + i + 1*nodes_Descent - 1 ];
        sigma_deg_Descent( i ) = x[ ( N + 3 ) + i + 2*nodes_Descent - 1 ];
        eps_T_deg_Descent( i ) = x[ ( N + 3 ) + i + 3*nodes_Descent - 1 ];
        phi_T_deg_Descent( i ) = x[ ( N + 3 ) + i + 4*nodes_Descent - 1 ];
        throttle_Descent( i )  = x[ ( N + 3 ) + i + 5*nodes_Descent - 1 ];
    }

    //! Declare and initialize number of parameters exclusive to Descent phase.
    const unsigned long NN = ( parameterList_Descent_.size() - 1 ) * nodes_Descent - 1;

    //! Declare and initialize last parameter common to the entire trajectory.
    //!     Final velocity.
    const double V_DN = x[ N + NN + 3 ];

    //std::cout << "Create vector of node locations for ascent" << std::endl;
    //! Create vector of node locations for Ascent phase.
    //!     This vector starts from zero, where it's last element is the sum of
    //!     all intervals.
    xn_Ascent( 0 ) = 0;
    for( unsigned int i = 0; i < nodes_Ascent - 1; i++ )
    {
        xn_Ascent( i + 1 )        = xn_Ascent( i ) + xn_interval_Ascent( i );
    }

    // std::cout << "Create vector of node locations for descent" << std::endl;
    //! Create vector of node locations for Descent phase.
    //!     This vector starts from the sum of all intervals, where it's last element
    //!     is zero.
    xn_Descent( 0 ) = xn_interval_Descent.sum();
    for( unsigned int i = 0; i < nodes_Descent - 1; i++ )
    {
        xn_Descent( i + 1 )        = xn_Descent( i ) - xn_interval_Descent( i );
    }


    // xn_Descent = xn_Descent.reverse().eval();

    /*
 *
 *
    std::cout << "x =  " << std::endl;
    for( int i = 0; i < int( x.size( ) ); i++ )
    {
        std::cout << i << "   " << x[i] << std::endl;
    }
    std::cout << "-------" << std::endl;
    std::cout << "xn_interval_Ascent =  " << xn_interval_Ascent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "xn_Ascent=  " << xn_Ascent<< std::endl;
    std::cout << "-------" << std::endl;

    std::cout << "alpha_deg_Ascent =  " << alpha_deg_Ascent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "eps_T_deg =  " << eps_T_deg << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "throttle =  " << throttle << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "xn_interval_Descent =  " << xn_interval_Descent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "xn_Descent=  " << xn_Descent<< std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "alpha_deg_Descent =  " << alpha_deg_Descent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "sigma_deg_Descent =  " << sigma_deg_Descent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "V_i =  " << V_i << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "V_UP =  " << V_UP << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "V_DN =  " << V_DN << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "h_UP =  " << h_UP << std::endl;
    std::cout << "-------" << std::endl;
*/

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE KNOWN STATES              ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating known states" << std::endl;

    //! Impose constraints on first and last energy nodes_Ascent
    //const double mu = spice_interface::getBodyGravitationalParameter( centralBodyName );
    //a = 301.7;//NRLMSISE00Atmosphere::getSpeedOfSound( R_E + height( 0 ), 0, 0, simulationStartEpoch );
    //double V_i = a * Mach_i;
    //double V_DN = 0.99 * sqrt( mu / ( R_E + h_UP ) );

    //! Set spherical elements for vehicle's initial state.
    //! ARTBITRARILY CHOSEN to be in Earth-Fixed frame.
    Eigen::Vector6d EntryState_spherical;
    EntryState_spherical( SphericalOrbitalStateElementIndices::radiusIndex )       = radius_Earth_i + h_i;
    EntryState_spherical( SphericalOrbitalStateElementIndices::latitudeIndex )     = initialLat_rad;
    EntryState_spherical( SphericalOrbitalStateElementIndices::longitudeIndex )    = initialLon_rad;
    EntryState_spherical( SphericalOrbitalStateElementIndices::speedIndex )        = V_i;
    EntryState_spherical( SphericalOrbitalStateElementIndices::flightPathIndex )   = initialFlightPathAngle_rad;
    EntryState_spherical( SphericalOrbitalStateElementIndices::headingAngleIndex ) = bislip::variables::computeHeadingToTarget( initialLat_rad,initialLon_rad,targetLat_rad,targetLon_rad );

    //! Two things are being done here
    //!     Converting state vector from spherical to Cartesian elements
    //!     Transforming state vector from Earth-Fixed frame to Inertial frame.
    const Eigen::Vector6d systemInitialState = transformStateToGlobalFrame(
                convertSphericalOrbitalToCartesianState( EntryState_spherical ),
                simulationStartEpoch,
                earthRotationalEphemeris );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE INTERPOLATORS             ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating interpolators" << std::endl;

    //! Calculate initial, maximum, and final specific energy levels.
    //const double E_min = bislip::variables::computeSpecificEnergy( h_i, V_i );
    const double E_max = bislip::variables::computeSpecificEnergy( h_UP, V_UP );
    //const double E_f = bislip::variables::computeSpecificEnergy( h_DN, V_DN );

    //! Set maximum Energy level.
    vehicleSystems->setE_max( E_max );

    //! Normalize initial, maximum, and final specific energy levels.
    const double E_hat_i = bislip::variables::computeNormalizedSpecificEnergy( h_i, V_i, E_max );
    const double E_hat_max = bislip::variables::computeNormalizedSpecificEnergy( h_UP, V_UP, E_max );
    const double E_hat_f = bislip::variables::computeNormalizedSpecificEnergy( h_DN, V_DN, E_max );

    //! Map normalized specific energy levels to control node locations.
    Eigen::VectorXd E_mapped_Ascent( xn_Ascent.size() );
    Eigen::VectorXd E_mapped_Descent( xn_Descent.size() );
    E_mapped_Ascent = ( E_hat_max - E_hat_i ) * xn_Ascent.array() + E_hat_i;
    E_mapped_Descent = ( E_hat_max - E_hat_f ) * xn_Descent.array() + E_hat_f;
    //E_mapped_Descent = E_mapped_Descent.reverse().eval();
    //std::cout << "xn_Descent=  " << xn_Descent<< std::endl;
    //std::cout << "E_mapped_Ascent=  " << E_mapped_Ascent<< std::endl;
    std::cout << "E_mapped_Descent=  " << E_mapped_Descent<< std::endl;

    //! Declare data maps used for decision vector values related to Ascent phase.
    std::map< double, double > map_alpha_deg_Ascent, map_eps_T_deg_Ascent, map_phi_T_deg_Ascent, map_throttle_Ascent, map_sigma_deg_Ascent;
    std::map< double, Eigen::VectorXd > map_DV_mapped_Ascent;
    Eigen::VectorXd DV_mapped_Ascent ( 6 );

    //std::cout << "Mapping Ascent DVs" << std::endl;
    //! Associate decision vector values to mapped normalized specific energy levels within data maps.
    for ( unsigned int i = 0; i < E_mapped_Ascent.size( ); ++i )
    {
        map_alpha_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 0 ) = alpha_deg_Ascent( i );
        map_sigma_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 1 ) = sigma_deg_Ascent( i );
        map_eps_T_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 2 ) = eps_T_deg_Ascent( i );
        map_phi_T_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 3 ) = phi_T_deg_Ascent( i );
        map_throttle_Ascent[ E_mapped_Ascent( i ) ]  = DV_mapped_Ascent( 4 ) = throttle_Ascent( i );
        DV_mapped_Ascent( 5 ) = xn_Ascent( i );
        map_DV_mapped_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent;
    }

    //! Declare data maps used for decision vector values related to Descent phase.
    std::map< double, double > map_alpha_deg_Descent, map_eps_T_deg_Descent, map_phi_T_deg_Descent, map_throttle_Descent, map_sigma_deg_Descent;
    std::map< double, Eigen::VectorXd > map_DV_mapped_Descent;
    Eigen::VectorXd DV_mapped_Descent ( 6 );

    //std::cout << "Mapping Descent DVs" << std::endl;
    //! Associate decision vector values to mapped normalized specific energy levels within data maps.
    for ( unsigned int i = 0; i < E_mapped_Descent.size( ); ++i )
    {
        map_alpha_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 0 ) = alpha_deg_Descent( i );
        map_sigma_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 1 ) = sigma_deg_Descent( i );
        map_eps_T_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 2 ) = eps_T_deg_Descent( i );
        map_phi_T_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 3 ) = phi_T_deg_Descent( i );
        map_throttle_Descent[ E_mapped_Descent( i ) ]  = DV_mapped_Descent( 4 ) = throttle_Descent( i );
        DV_mapped_Descent( 5 ) = xn_Descent( i );
        map_DV_mapped_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent;
    }

    //! Calculate first derivatives for Hermite Spline Interpolation according to https://doi.org/10.1016/j.cam.2017.09.049.
    std::vector< double > alpha_deg_Ascent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, alpha_deg_Ascent, E_mapped_Ascent.size() );
    std::vector< double > sigma_deg_Ascent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, sigma_deg_Ascent, E_mapped_Ascent.size() );
    std::vector< double > eps_T_deg_Ascent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, eps_T_deg_Ascent, E_mapped_Ascent.size() );
    std::vector< double > phi_T_deg_Ascent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, phi_T_deg_Ascent, E_mapped_Ascent.size() );
    std::vector< double > throttle_Ascent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, throttle_Ascent, E_mapped_Ascent.size() );
    std::vector< double > alpha_deg_Descent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Descent, alpha_deg_Descent, E_mapped_Descent.size() );
    std::vector< double > sigma_deg_Descent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, sigma_deg_Descent, E_mapped_Descent.size() );
    std::vector< double > eps_T_deg_Descent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, eps_T_deg_Descent, E_mapped_Descent.size() );
    std::vector< double > phi_T_deg_Descent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, phi_T_deg_Descent, E_mapped_Descent.size() );
    std::vector< double > throttle_Descent_derivatives = bislip::variables::HermiteDerivatives( E_mapped_Ascent, throttle_Descent, E_mapped_Descent.size() );

    //! Declare and initialize interpolator settings.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = std::make_shared< interpolators::InterpolatorSettings >( hermite_spline_interpolator );

    //std::cout << "Creating Ascent Interpolators" << std::endl;
    //! Declare and initialize interpolators for Ascent phase.
    std::shared_ptr< OneDimensionalInterpolator< double, double > > interpolator_alpha_deg_Ascent, interpolator_sigma_deg_Ascent, interpolator_eps_T_deg_Ascent, interpolator_phi_T_deg_Ascent, interpolator_throttle_Ascent;
    interpolator_alpha_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_alpha_deg_Ascent, interpolatorSettings, std::make_pair( alpha_deg_Ascent( 0 ), alpha_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), alpha_deg_Ascent_derivatives );
    interpolator_sigma_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_sigma_deg_Ascent, interpolatorSettings, std::make_pair( sigma_deg_Ascent( 0 ), sigma_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), sigma_deg_Ascent_derivatives );
    interpolator_eps_T_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_eps_T_deg_Ascent, interpolatorSettings, std::make_pair( eps_T_deg_Ascent( 0 ), eps_T_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), eps_T_deg_Ascent_derivatives );
    interpolator_phi_T_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_phi_T_deg_Ascent, interpolatorSettings, std::make_pair( phi_T_deg_Ascent( 0 ), phi_T_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), phi_T_deg_Ascent_derivatives );
    interpolator_throttle_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_throttle_Ascent, interpolatorSettings, std::make_pair( throttle_Ascent( 0 ), throttle_Ascent( E_mapped_Ascent.size() - 1 ) ), throttle_Ascent_derivatives );

    //std::cout << "Creating Descent Interpolators" << std::endl;
    //! Declare and initialize interpolators for Descent phase.
    std::shared_ptr< OneDimensionalInterpolator< double, double > > interpolator_alpha_deg_Descent, interpolator_sigma_deg_Descent, interpolator_eps_T_deg_Descent, interpolator_phi_T_deg_Descent, interpolator_throttle_Descent;
    interpolator_alpha_deg_Descent = interpolators::createOneDimensionalInterpolator< double, double >( map_alpha_deg_Descent, interpolatorSettings, std::make_pair( alpha_deg_Descent( 0 ), alpha_deg_Descent( E_mapped_Descent.size() - 1 ) ), alpha_deg_Descent_derivatives );
    interpolator_sigma_deg_Descent = interpolators::createOneDimensionalInterpolator< double, double >( map_sigma_deg_Descent, interpolatorSettings, std::make_pair( sigma_deg_Descent( 0 ), sigma_deg_Descent( E_mapped_Descent.size() - 1 ) ), sigma_deg_Descent_derivatives );
    interpolator_eps_T_deg_Descent = interpolators::createOneDimensionalInterpolator< double, double >( map_eps_T_deg_Descent, interpolatorSettings, std::make_pair( eps_T_deg_Descent( 0 ), eps_T_deg_Descent( E_mapped_Descent.size() - 1 ) ), eps_T_deg_Descent_derivatives );
    interpolator_phi_T_deg_Descent = interpolators::createOneDimensionalInterpolator< double, double >( map_phi_T_deg_Descent, interpolatorSettings, std::make_pair( phi_T_deg_Descent( 0 ), phi_T_deg_Descent( E_mapped_Descent.size() - 1 ) ), phi_T_deg_Descent_derivatives );
    interpolator_throttle_Descent  = interpolators::createOneDimensionalInterpolator< double, double >( map_throttle_Descent,  interpolatorSettings, std::make_pair( throttle_Descent( 0 ), throttle_Descent( E_mapped_Descent.size() - 1 ) ), throttle_Descent_derivatives );

    // std::cout << "Pass Ascent Interpolators to vehicle systems" << std::endl;
    //! Pass Ascent phase interpolators to vehicle systems.
    vehicleSystems->setAoAInterpolator_Ascent( interpolator_alpha_deg_Ascent );
    vehicleSystems->setBankAngleInterpolator_Ascent( interpolator_sigma_deg_Ascent );
    vehicleSystems->setThrustElevationAngleInterpolator_Ascent( interpolator_eps_T_deg_Ascent );
    vehicleSystems->setThrustAzimuthAngleInterpolator_Ascent( interpolator_phi_T_deg_Ascent );
    vehicleSystems->setThrottleInterpolator_Ascent( interpolator_throttle_Ascent );

    //std::cout << "Pass Descent Interpolators to vehicle systems" << std::endl;
    //! Pass Descent phase interpolators to vehicle systems.
    vehicleSystems->setAoAInterpolator_Descent( interpolator_alpha_deg_Descent );
    vehicleSystems->setBankAngleInterpolator_Descent( interpolator_sigma_deg_Descent );
    vehicleSystems->setThrustElevationAngleInterpolator_Descent( interpolator_eps_T_deg_Descent );
    vehicleSystems->setThrustAzimuthAngleInterpolator_Descent( interpolator_phi_T_deg_Descent );
    vehicleSystems->setThrottleInterpolator_Descent( interpolator_throttle_Descent );

    //! Declare vectors containing interpolated values.
    Eigen::VectorXd interpolated_values_Ascent( 5 ), interpolated_values_Descent( 5 );

    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::VectorXd > interpolators_Ascent, interpolators_Descent;

    //! Declare evaluation variable.
    double eval;

    //std::cout << "Evaluate Interpolators and store to print out" << std::endl;
    //! Loop to populate vectors of interpolated values and then pass to data map.
    //!     Number of evaluations has been arbitrarily selected.
    double pp = 0;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        eval = pp * E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) / 1000;
        interpolated_values_Ascent( 0 ) = interpolator_alpha_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 1 ) = interpolator_sigma_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 2 ) = interpolator_eps_T_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 3 ) = interpolator_phi_T_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 4 ) = interpolator_throttle_Ascent->interpolate( eval );

        interpolators_Ascent[ eval ] = interpolated_values_Ascent;

        eval = pp * E_mapped_Descent( 0 ) / 1000;
        interpolated_values_Descent( 0 ) = interpolator_alpha_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 1 ) = interpolator_sigma_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 2 ) = interpolator_eps_T_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 3 ) = interpolator_phi_T_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 4 ) = interpolator_throttle_Descent->interpolate( eval );

        interpolators_Descent[ eval ] = interpolated_values_Descent;

        pp += 1;
    }

    //std::cout << "Print out Interpolators and DVs" << std::endl;
    //! Print data maps containing vectors of evaluation of interpolators.
    writeDataMapToTextFile( interpolators_Ascent,
                            "interpolators_pre_Ascent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
    writeDataMapToTextFile( map_DV_mapped_Ascent,
                            "DV_mapped_pre_Ascent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
    writeDataMapToTextFile( interpolators_Descent,
                            "interpolators_pre_Descent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
    writeDataMapToTextFile( map_DV_mapped_Descent,
                            "DV_mapped_pre_Descent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
    /*

    std::map< double, double > map_sin, map_sincos, map_x_squared, map_x_cubed, map_sin_x_squared, map_sin_of_x_squared, map_exp, map_exp_sin;
    std::map< double, Eigen::VectorXd > map_interpolatorValidation;
    Eigen::VectorXd interpolatorValidation ( 8 );

    Eigen::VectorXd eval_vector(100), sin(100), sincos(100), x_squared(100), x_cubed(100), sin_x_squared(100), sin_of_x_squared(100), exp(100), exp_sin(100);

    //! Associate decision vector values to mapped normalized specific energy levels within data maps.
    for ( unsigned int i = 0; i < 100; ++i )
    {
        eval = (double)i / 10;
        eval_vector(i) = eval;
        map_sin[ eval ]                    = interpolatorValidation( 0 ) = sin(i) = std::sin( eval );
        map_sincos[ eval ]                 = interpolatorValidation( 1 ) = sincos(i) = std::sin( eval ) * std::cos( eval );
        map_x_squared[ eval ]              = interpolatorValidation( 2 ) = x_squared(i) = eval * eval;
        map_x_cubed[ eval ]                = interpolatorValidation( 3 ) = x_cubed(i) = eval * eval * eval;
        map_sin_x_squared[ eval ]          = interpolatorValidation( 4 ) = sin_x_squared(i) = std::sin( eval ) * eval * eval;
        map_sin_of_x_squared[ eval ]       = interpolatorValidation( 5 ) = sin_of_x_squared(i) = std::sin( eval * eval );
        map_exp[ eval ]                    = interpolatorValidation( 6 ) = exp(i) = std::exp( eval );
        map_exp_sin[ eval ]                = interpolatorValidation( 7 ) = exp_sin(i) =  std::sin( eval ) * std::exp( eval );
        map_interpolatorValidation[ eval ] = interpolatorValidation;
    }

    std::vector< double > sin_derivatives              = bislip::variables::HermiteDerivatives( eval_vector, sin, 100 );
    std::vector< double > sincos_derivatives           = bislip::variables::HermiteDerivatives( eval_vector, sincos, 100 );
    std::vector< double > x_squared_derivatives        = bislip::variables::HermiteDerivatives( eval_vector, x_squared, 100 );
    std::vector< double > x_cubed_derivatives          = bislip::variables::HermiteDerivatives( eval_vector, x_cubed, 100 );
    std::vector< double > sin_x_squared_derivatives    = bislip::variables::HermiteDerivatives( eval_vector, sin_x_squared, 100 );
    std::vector< double > sin_of_x_squared_derivatives = bislip::variables::HermiteDerivatives( eval_vector, sin_of_x_squared, 100 );
    std::vector< double > exp_derivatives              = bislip::variables::HermiteDerivatives( eval_vector, exp, 100 );
    std::vector< double > exp_sin_derivatives          = bislip::variables::HermiteDerivatives( eval_vector, exp_sin, 100 );


    std::shared_ptr< OneDimensionalInterpolator< double, double > > interpolator_sin, interpolator_sincos, interpolator_x_squared, interpolator_x_cubed, interpolator_sin_x_squared, interpolator_sin_of_x_squared, interpolator_exp, interpolator_exp_sin;

    interpolator_sin              = interpolators::createOneDimensionalInterpolator< double, double >( map_sin, interpolatorSettings, std::make_pair( sin( 0 ), sin( 100 ) ), sin_derivatives);
    interpolator_sincos           = interpolators::createOneDimensionalInterpolator< double, double >( map_sincos, interpolatorSettings, std::make_pair( sincos( 0 ), sincos( 100 ) ), sincos_derivatives);
    interpolator_x_squared        = interpolators::createOneDimensionalInterpolator< double, double >( map_x_squared, interpolatorSettings, std::make_pair( x_squared( 0 ), x_squared( 100 ) ), x_squared_derivatives);
    interpolator_x_cubed          = interpolators::createOneDimensionalInterpolator< double, double >( map_x_cubed, interpolatorSettings, std::make_pair( x_cubed( 0 ), x_cubed( 100 ) ), x_cubed_derivatives);
    interpolator_sin_x_squared    = interpolators::createOneDimensionalInterpolator< double, double >( map_sin_x_squared, interpolatorSettings, std::make_pair( sin_x_squared( 0 ), sin_x_squared( 100 ) ), sin_x_squared_derivatives);
    interpolator_sin_of_x_squared = interpolators::createOneDimensionalInterpolator< double, double >( map_sin_of_x_squared, interpolatorSettings, std::make_pair( sin_of_x_squared( 0 ), sin_of_x_squared( 100) ), sin_of_x_squared_derivatives);
    interpolator_exp              = interpolators::createOneDimensionalInterpolator< double, double >( map_exp, interpolatorSettings, std::make_pair( exp( 0 ), exp( 100 ) ), exp_derivatives);
    interpolator_exp_sin          = interpolators::createOneDimensionalInterpolator< double, double >( map_exp_sin, interpolatorSettings, std::make_pair( exp_sin( 0 ), exp_sin( 100 ) ), exp_sin_derivatives);



    //! Declare vectors containing interpolated values.
    Eigen::VectorXd interpolater_Validation( 8 );

    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::VectorXd > interpolators_Validation;

    for ( unsigned int i = 0; i < 100001; ++i )
    {
        eval = (double)i / 1000;
        interpolater_Validation( 0 ) = interpolator_sin->interpolate( eval );
        interpolater_Validation( 1 ) = interpolator_sincos->interpolate( eval );
        interpolater_Validation( 2 ) = interpolator_x_squared->interpolate( eval );
        interpolater_Validation( 3 ) = interpolator_x_cubed->interpolate( eval );
        interpolater_Validation( 4 ) = interpolator_sin_x_squared->interpolate( eval );
        interpolater_Validation( 5 ) = interpolator_sin_of_x_squared->interpolate( eval );
        interpolater_Validation( 6 ) = interpolator_exp->interpolate( eval );
        interpolater_Validation( 7 ) = interpolator_exp_sin->interpolate( eval );

        interpolators_Validation[ eval ] = interpolater_Validation;


    }

    //! Print data maps containing vectors of evaluation of interpolators.
    writeDataMapToTextFile( interpolators_Validation,
                            "interpolators_Validation",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );



    writeDataMapToTextFile( map_interpolatorValidation,
                            "map_interpolatorValidation",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

*/
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE PARAMETER BOUNDS          ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Declare data maps of strings and pairs to store parameter bounds.
    //!     Strings are used for the parameter names.
    //!     Pairs are used for lower and upper bounds.
    std::map< std::string, std::pair < double, double > > Bounds_Ascent, Bounds_Descent;

    //! Declare and initialize counter.
    //!     Information is extracted directly from vectors created from input data, which requires the counter for the upper bound value.
    unsigned long p = 0;

    //! Loop to populate the data map with pairs associated to string values related to Ascent phase.
    for( unsigned long i = 0; i < parameterList_Ascent_.size(); i++)
    {
        Bounds_Ascent[ parameterList_Ascent_[ i ] ] = std::make_pair( parameterBounds_Ascent_[ p ], parameterBounds_Ascent_[ p + 1 ] );
        p += 2;
    }

    //! Re-initialize counter.
    p = 0;

    //! Loop to populate the data map with pairs associated to string values related to Descent phase.
    for( unsigned long i = 0; i < parameterList_Descent_.size(); i++)
    {
        Bounds_Descent[ parameterList_Descent_[ i ] ] = std::make_pair( parameterBounds_Descent_[ p ], parameterBounds_Descent_[ p + 1 ] );
        p += 2;
    }

    //! Pass parameter bounds to vehicle systems.
    vehicleSystems->setParameterBounds( Bounds_Ascent, "ascent" );
    vehicleSystems->setParameterBounds( Bounds_Descent, "descent" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Setting Accelerations" << std::endl;

    //! Declare acceleration data map.
    SelectedAccelerationMap accelerationSettingsMap;

    //! Define gravitational model.
    //!     The central body acts this force on the vehicle.
    //!     Arbitrary maximum degree/order. Equivalent functionality to Cartesian with corresponding maximum degree/order.
    accelerationSettingsMap[ vehicle_name_ ][ centralBodyName ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    //! Define aerodynamic accelerations.
    //!     The atmosphere of the central body acts this force on the vehicle.
    accelerationSettingsMap[ vehicle_name_ ][ centralBodyName ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

    //! Declare and initialize thrust guidance pointer.
    std::shared_ptr< MyGuidance > ThrustGuidance = std::make_shared< MyGuidance >(
                bodyMap,
                vehicle_name_);

    //! Declare and initialize thrust guidance direction guidance settings.
    std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< ThrustDirectionGuidanceSettings >( thrust_direction_from_existing_body_orientation, centralBodyName );

    //! Declare and initialize thrust magnitude function.
    std::function< double( const double ) > thrustMagnitudeFunction =
            std::bind( &bislip::MyGuidance::getCurrentThrustMagnitude, ThrustGuidance );

    //! Declare and initialize thrust specific impulse function.
    std::function< double( const double ) > specificImpulseFunction =
            std::bind( &bislip::MyGuidance::getCurrentSpecificImpulse, ThrustGuidance );

    //! Declare and initialize engine status function.
    std::function< bool( const double ) > isEngineOnFunction =
            std::bind( &bislip::MyGuidance::getCurrentEngineStatus, ThrustGuidance );

    //! Declare and initialize thrust direction vector function.
    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection =
            std::bind( &bislip::MyGuidance::getCurrentBodyFixedThrustDirection, ThrustGuidance );

    //! Declare and initialize thrust guidance update function.
    std::function< void( const double ) > customThrustResetFunction =
            std::bind( &bislip::MyGuidance::updateGuidance, ThrustGuidance, std::placeholders::_1 );

    //! Declare and initialize thrust magnitude settings.
    std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction,
                specificImpulseFunction,
                isEngineOnFunction,
                bodyFixedThrustDirection,
                customThrustResetFunction );

    //! Pass engine capabilities to vehicle systems.
    vehicleSystems->setMaxThrust( maxThrust );
    vehicleSystems->setSpecificImpulse( specificImpulse );

    //! Define thrust acceleration settings.
    //!     The vehicle acts this force on itself.
    accelerationSettingsMap[ vehicle_name_ ][ vehicle_name_ ].push_back( std::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );

    //std::cout << "Accelerations Set" << std::endl;

    //! Declare vectors that will contain the names of relevant bodies.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;

    //! Define bodies that will be propagated. Only 1.
    bodiesToIntegrate.push_back( vehicle_name_ );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( centralBodyName );

    //! Set acceleration models.
    basic_astrodynamics::AccelerationMap accelerationsMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationSettingsMap,
                bodiesToIntegrate,
                centralBodies );

    //! **************************************************************************************
    //! ********** Aerodynamic guidance is set AFTER the accelerations and BEFORE propagating.
    //! **************************************************************************************

    //! Declare and assign aerodynamic guidance functions.
    std::shared_ptr< MyGuidance > AeroGuidance = std::make_shared< MyGuidance >(
                bodyMap,
                vehicle_name_);

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
    ///////////////////////          CREATE LIST OF DEPENDENT VARIABLES        ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating list of dependent variables" << std::endl;

    //! Create vector that will contian the list of dependent variables to save/output.
    //!     The file that prints out the text saying what is saved has been modified
    //!     to mute this message. During optimization this would clutter the
    //!     terminal screen and is incredibly annoying.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dep_varList;

    //! Add dependent variables to list.
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
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicle_name_,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    airspeed_dependent_variable,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_aerodynamic_g_load_variable,
                    vehicle_name_,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic,
                    vehicle_name_,
                    centralBodyName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components.
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    spherical_harmonic_gravity,
                    vehicle_name_,
                    centralBodyName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration,
                    vehicle_name_,
                    vehicle_name_,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components.
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_dynamic_pressure_dependent_variable,
                    vehicle_name_,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_aerodynamic_heat_rate_dependent_variable,
                    vehicle_name_,
                    centralBodyName) );
    /*dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    stagnation_point_heat_flux_dependent_variable,
                    vehicle_name_,
                    centralBodyName) );*/
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_mass_rate_dependent_variables,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    current_mass,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    specific_energy,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    normalized_specific_energy,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    throttle_setting,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    thrust_elevation_angle,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    engine_status,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    angular_distance_traveled,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    angular_distance_to_go,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heading_to_target,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heading_error,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heat_rate_leading_edge,
                    vehicle_name_ ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_thrust_vector,
                    vehicle_name_ ) );

    //! Create and initialize object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dep_varToSave =
            std::make_shared< DependentVariableSaveSettings >( dep_varList );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE TERMINATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    vehicleSystems->setInitialLat( initialLat_rad );
    vehicleSystems->setInitialLon( initialLon_rad );
    vehicleSystems->setInitialCoordinates( std::make_pair( vehicleSystems->getInitialLat( ), vehicleSystems->getInitialLon( ) ) );
    vehicleSystems->setTargetLat( targetLat_rad );
    vehicleSystems->setTargetLon( targetLon_rad );
    vehicleSystems->setTargetCoordinates( std::make_pair( vehicleSystems->getTargetLat( ), vehicleSystems->getTargetLon( ) ) );
    vehicleSystems->setInitialDistanceToTarget( initialDistanceToTarget_rad );
    vehicleSystems->setStartingEpoch( simulationStartEpoch );

    //! Define CUSTOM termination settings.
    std::shared_ptr< PropagationTerminationSettings > customTermination =
            std::make_shared< PropagationCustomTerminationSettings >(
                boost::bind( &bislip::StopOrNot, bodyMap, vehicle_name_, vehicleParameterValues_, terminationConditionsValues_ ) );

    //! Define dependent variable termination settings.
    std::shared_ptr< PropagationTerminationSettings > thrustTerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration,
                    vehicle_name_,
                    vehicle_name_,
                    true ), maxThrust, false );

    std::shared_ptr< PropagationTerminationSettings > d_to_go_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    angular_distance_to_go,
                    vehicle_name_), vehicleSystems->getFinalDistanceToTarget() , true );

    std::shared_ptr< PropagationTerminationSettings > d_traveled_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    angular_distance_traveled,
                    vehicle_name_), vehicleSystems->getInitialDistanceToTarget(), false );

    std::shared_ptr< PropagationTerminationSettings > mass_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    current_mass,
                    vehicle_name_), vehicleSystems->getDryMass(), true );

    std::shared_ptr< PropagationTerminationSettings > E_hat_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    normalized_specific_energy,
                    vehicle_name_), 1.0 , false );

    std::shared_ptr< PropagationTerminationSettings > q_dyn_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_dynamic_pressure_dependent_variable,
                    vehicle_name_,
                    centralBodyName), q_dyn_max , false );

    std::shared_ptr< PropagationTerminationSettings > q_dot_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_aerodynamic_heat_rate_dependent_variable,
                    vehicle_name_,
                    centralBodyName), q_dot_max , false );
    std::shared_ptr< PropagationTerminationSettings > height_UP_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicle_name_,
                    centralBodyName), h_UP , false );
    std::shared_ptr< PropagationTerminationSettings > height_DN_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicle_name_,
                    centralBodyName), h_DN , true );
    std::shared_ptr< PropagationTerminationSettings > flight_path_angle_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicle_name_,
                    flight_path_angle,
                    centralBodyName) , tudat::mathematical_constants::PI / 2.0 , false );


    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList;
    // terminationSettingsList.push_back( customTermination );
    //terminationSettingsList.push_back( thrustTerminationSettings );
    terminationSettingsList.push_back( d_to_go_TerminationSettings );
    terminationSettingsList.push_back( d_traveled_TerminationSettings );
    terminationSettingsList.push_back( mass_TerminationSettings );
    terminationSettingsList.push_back( E_hat_TerminationSettings );
    terminationSettingsList.push_back( q_dyn_TerminationSettings );
    terminationSettingsList.push_back( q_dot_TerminationSettings );
    terminationSettingsList.push_back( height_UP_TerminationSettings );
    terminationSettingsList.push_back( height_DN_TerminationSettings );
    terminationSettingsList.push_back( flight_path_angle_TerminationSettings );

    // PropagationHybridTerminationSettings( terminationSettingsList,
    //                                       true,
    //                                       false );

    std::shared_ptr< PropagationTerminationSettings > terminationSettings = std::make_shared<
            propagators::PropagationHybridTerminationSettings >( terminationSettingsList, true );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Create translational propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies,
              accelerationsMap,
              bodiesToIntegrate,
              systemInitialState,
              terminationSettings,
              cowell,
              dep_varToSave );

    ///////////////////////             CREATE MASS RATE SETTINGS              ////////////////////////////////////////////
    //std::cout << "Create mass rate models" << std::endl;

    //! Declare and initialize mass rate model settings.
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );

    //! Declare and initialize mass rate model settings.
    std::map< std::string, std::shared_ptr< MassRateModel > > massRateModels;
    massRateModels[ vehicle_name_ ] = createMassRateModel(
                vehicle_name_, massRateModelSettings, bodyMap, accelerationsMap );

    //! Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( vehicle_name_ );

    //! Declare and initialize starting mass of vehicle.
    Eigen::VectorXd initialBodyMasses( 1 );
    initialBodyMasses( 0 ) = initialMass;

    //! Declare and initialize mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses, terminationSettings, dep_varToSave );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Declare and initialize list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    //! Declare and initialize propagation settings for both mass and translational dynamics.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings, dep_varToSave );


    //! **************************************************************************************
    //! ********** Create propagation settings for ONLY translational dynamics.

    //std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
    //        translationalPropagatorSettings;

    //! **************************************************************************************



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Declare and initialize integrator settings.
    std::shared_ptr< IntegratorSettings<  > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              fixedStepSize, 1, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          PROPAGATE TRAJECTORY          ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // std::cout << "Starting propagation" << std::endl;

    //! Propagate trajectory.
    SingleArcDynamicsSimulator< double > dynamicsSimulator(
                bodyMap,
                integratorSettings,
                propagatorSettings );

    //std::cout << "Propagation done" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        CLEAN UP AND CALCULATE FITNESS                //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Initialize variables used objectives functions.
    double finalDistanceToTarget_deg_calc = initialDistanceToTarget_deg;
    //std::cout << d_f_deg_calc << std::endl;
    double h_DN_calc = h_DN;
    double h_UP_calc = h_UP;
    double tof = simulation_settingsValues_[ 1 ];

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

    //! Extract Dependent variables of final state.
    const Eigen::VectorXd dep_var_FINAL_STATE = ( dynamicsSimulator.getDependentVariableHistory( ).rbegin() )->second;


    //! Extract Dependent variables of final state.
    const std::map< double, Eigen::VectorXd > dep_var_historyMap = dynamicsSimulator.getDependentVariableHistory( );
    // std::cout << "dep_var_FINAL_STATE: " << dep_var_FINAL_STATE << std::endl;

    //! Declare number of elements in time history.
    unsigned long rows = dep_var_historyMap.size();
    unsigned long columns = ( ( dynamicsSimulator.getDependentVariableHistory( ).begin() )->second ).size();

    Eigen::MatrixXd dep_var_history( rows, columns );
    dep_var_history = Eigen::MatrixXd::Zero( rows, columns );

    //    for ( std::map< double, Eigen::VectorXd >::iterator it = dep_var_history.begin(); it != dep_var_history.end(); ++it )
    for ( unsigned long i = 0; i < rows; i++ )
    {
        dep_var_history.row( i ) = dep_var_historyMap.at( simulationStartEpoch + i * fixedStepSize );
    }

    Eigen::VectorXd dep_var_height = dep_var_history.col( 11 );

    //const double altitude_f_calc = dep_var_FINAL_STATE[3];
    const double targetLat_rad_calc  = dep_var_FINAL_STATE[ 4 ];
    const double targetLon_rad_calc  = dep_var_FINAL_STATE[ 5 ];
    std::ptrdiff_t ind;
    h_UP_calc                    = dep_var_height.maxCoeff(&ind);
    h_DN_calc                    = dep_var_FINAL_STATE[ 11 ];
    //std::cout << "lat_f_rad_calc: " << lat_f_rad_calc << std::endl;
    //std::cout << "lon_f_rad_calc: " << lon_f_rad_calc << std::endl;
    //std::cout << "h_UP_calc     : " << h_UP_calc << std::endl;
    //std::cout << "h_DN_calc     : " << h_DN_calc << std::endl;

    //const double altitude_f_calc = std::sqrt( pow(systemFinalState_EARTh_UPIXED[0],2) +
    //       pow(systemFinalState_EARTh_UPIXED[1],2) + pow(systemFinalState_EARTh_UPIXED[2],2) ) ;
    //const double lon_f_rad_calc = std::atan2(systemFinalState_EARTh_UPIXED[1] , systemFinalState_EARTh_UPIXED[0]);
    //const double lat_f_rad_calc = std::asin(systemFinalState_EARTh_UPIXED[2] / altitude_f_calc) ;

    //! Convert coordinates of final state to degrees: Earth-Fixed Frame
    targetLat_deg_calc = unit_conversions::convertRadiansToDegrees( targetLat_rad_calc );
    targetLon_deg_calc = unit_conversions::convertRadiansToDegrees( targetLon_rad_calc );

    //! Calculate angular distance of final state from target coordinates.
    // const double d_rad = bislip::variables::computeAngularDistance( lat_f_rad_calc, lon_f_rad_calc, vehicleSystems->getTargetLat(), vehicleSystems->getTargetLon() );
    const double d_to_go_deg = unit_conversions::convertRadiansToDegrees( dep_var_FINAL_STATE[ 34 ] );
    //const double d_deg = unit_conversions::convertRadiansToDegrees( d_rad );

    //! Calculate offset of final state from GOAL state: Earth-Fixed Frame
    //const double dif_lat_rad = targetLat_rad - lat_f_rad_calc;
    //const double dif_lon_rad = targetLon_rad - lon_f_rad_calc;

    //! Calculate offsets of final state in degrees: Earth-Fixed Frame
    const double dif_lat_deg = targetLat_deg - targetLat_deg_calc;
    const double dif_lon_deg = targetLon_deg - targetLon_deg_calc;

    //! Calculate "norm" of offsets. This is an arbitrary function I have
    //! implemented to pass on as an 'objective function'. It relates the
    //! differences such that when minimizing the offsets, there is an additional
    //! unsigned value that always goes to zero. Most definitely unsure about how
    //!  'proper' it is, yet is what works for the current BALLISTIC case.
    const double dif_norm = std::sqrt( ( dif_lat_deg * dif_lat_deg ) + ( dif_lon_deg * dif_lon_deg ) );

    //! Calculate offset of final angular distance to termination condition distance.
    const double dif_d_deg = d_to_go_deg - terminationConditionsValues_[ 2 ];

    //! Calculate offset from maximum elevation.
    const double dif_h_UP = h_UP - h_UP_calc;

    //! Calculate offset from final elevation.
    const double dif_h_DN = h_DN - h_DN_calc;

    const double dif_xn_Ascent= xn_Ascent( E_mapped_Ascent.size() - 1 ) - 1;
    const double dif_E_mapped_Ascent = E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) - 1;
    const double dif_xn_Descent= xn_Descent( E_mapped_Descent.size() - 1 ) - 1;
    //const double dif_E_mapped_Descent_f = E_hat_f - E_mapped_Descent( nodes_Descent - 1 );
    const double dif_E_mapped_Descent_i = E_mapped_Descent( 0 ) - 1;

    //! Exract current mass
    double landingMass_calc = bodyMap.at( vehicle_name_ )->getBodyMass( );
    if ( isnan( landingMass_calc ) == 1 )
    {
        landingMass_calc = initialMass;
    }
    //! Calculate offset from goal mass.
    const double dif_mass = abs( landingMass - landingMass_calc );

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
    delta.push_back( dif_h_UP );
    delta.push_back( dif_h_DN );
    delta.push_back( tof );  // Not sure yet how this one affects the optimization. Included for completion.
    delta.push_back( dif_xn_Ascent);
    delta.push_back( dif_E_mapped_Ascent );
    delta.push_back( dif_xn_Descent);
   // delta.push_back( dif_E_mapped_Descent_i );
    delta.push_back( dif_E_mapped_Descent_i );
    delta.push_back( dif_mass );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PRINT SIMULATION OUTPUT TO FILE               //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Get time stamp for this specific simulation. This avoids overwriting the
    //! file if another individual with the same properties shows up in other
    //! evolutions.
    std::string simulation_save_time = bislip::variables::getCurrentDateTime( false );

    //! Create unique filename that cannot be overwritten due to the timestamp.
    std::string simulation_file_name_suffix =
            // std::to_string( Va_ordered( 0 ) ) + "_" +
            // std::to_string( gamma_rad_ordered( 0 ) ) + "_" +
            // std::to_string( chi_i_deg_ordered( 0 ) ) + "_" +
            std::to_string( tof ) + "_" +
            //std::to_string(x[3]) + "_" +
            std::to_string( targetLat_deg_calc ) + "_" +
            std::to_string( targetLon_deg_calc ) + "_" +
            simulation_save_time;

    //! Print results to terminal. Used to gauge progress.
    std::cout << std::fixed << std::setprecision(10) <<
                 std::setw(15) << "V_i = " <<
                 std::setw(14) << V_i <<
                 std::setw(15) << "dif_xn_Ascent= " <<
                 std::setw(14) << dif_xn_Ascent<<
                 std::setw(15) << " dif_E_mapped_Ascent = " <<
                 std::setw(16) <<  dif_E_mapped_Ascent <<
                 std::setw(15) << "dif_xn_Descent= " <<
                 std::setw(14) << dif_xn_Descent<<
                 std::setw(15) << " dif_E_mapped_Descent_i = " <<
                 std::setw(16) <<  dif_E_mapped_Descent_i <<
                 std::setw(15) << " dif_mass = " <<
                 std::setw(16) <<  dif_mass <<
                 std::setw(15) << " dif_norm = " <<
                 std::setw(16) << dif_norm <<
                 std::setw(15) << " dif_lat_deg = " <<
                 std::setw(16) << dif_lat_deg <<
                 std::setw(15) << " dif_lon_deg = " <<
                 std::setw(16) << dif_lon_deg <<
                 std::setw(13) << "dif_d_deg = " <<
                 std::setw(16) << dif_d_deg <<
                 std::setw(9) << "dif_h_UP = " <<
                 std::setw(16) << dif_h_UP <<
                 std::setw(7) << "tof = " <<
                 std::setw(16) << tof <<
                 std::setw(60) << simulation_file_name_suffix << std::endl;

    std::string complete_file_name_Prop = "HORUSPropHistory_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_DepVar = "HORUSDepVar_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_interpolators_Ascent = "interpolators_Ascent_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_map_DV_mapped_Ascent = "map_DV_mapped_Ascent_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_interpolators_Descent = "interpolators_Descent_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_map_DV_mapped_Descent = "map_DV_mapped_Descent_" + simulation_file_name_suffix + ".dat";

    //! Will print out depending on some input values. Each entry corresponds to
    //! a different type of output. Entries are binary, 0 or 1.
    if ( int( output_settingsValues_[ 1 ] ) == 1 )
    {
        //   std::cout << "Saving EQM " << std::endl;

        //! Write propagation history to file.
        writeDataMapToTextFile( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                complete_file_name_Prop,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //   std::cout << "EQM Saved" << std::endl;

    }
    if ( int( output_settingsValues_[ 2 ] ) == 1 )
    {
        //   std::cout << "Saving DepVar " << std::endl;

        //! Write dependent variables' history to file.
        writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory( ),
                                complete_file_name_DepVar,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //    std::cout << "DepVar Saved" << std::endl;


    }
    if ( int( output_settingsValues_[ 3 ] ) == 1 )
    {
        //   std::cout << "Saving interpolators_Ascent " << std::endl;

        writeDataMapToTextFile( interpolators_Ascent,
                                complete_file_name_interpolators_Ascent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //   std::cout << "interpolators_Ascent Saved" << std::endl;
        //   std::cout << "Saving map_DV_mapped_Ascent " << std::endl;

        writeDataMapToTextFile( map_DV_mapped_Ascent,
                                complete_file_name_map_DV_mapped_Ascent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //   std::cout << "map_DV_mapped_Ascent Saved" << std::endl;

        //   std::cout << "Saving interpolators_Descent " << std::endl;

        writeDataMapToTextFile( interpolators_Descent,
                                complete_file_name_interpolators_Descent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //  std::cout << "interpolators_Descent Saved" << std::endl;
        //   std::cout << "Saving map_DV_mapped_Descent " << std::endl;

        writeDataMapToTextFile( map_DV_mapped_Descent,
                                complete_file_name_map_DV_mapped_Descent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //std::cout << "map_DV_mapped_Descent Saved" << std::endl;


    }


    //  }



    //std::cout << "h_UP: " << h_UP << std::endl;
    //std::cout << "h_UP_calc: " << h_UP_calc << std::endl;



    // std::cout <<dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) << std::endl;
    //std::this_thread::sleep_for( std::chrono::nanoseconds( 10 ) );
    //using namespace std::this_thread; // sleep_for, sleep_until
    //using namespace std::chrono; // nanoseconds, system_clock, seconds




    return delta;

} // Fitness function.
