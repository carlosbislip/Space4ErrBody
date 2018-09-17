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

//#include <iostream>
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
#include <boost/filesystem/operations.hpp>

//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include <boost/random/variate_generator.hpp>

#include <Eigen/Core>
#include <boost/bind.hpp>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>



//! Do I need a costumized one for HORUS? Probably....
#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>

//! Mine
#include "Space4ErrBody.h"
#include "applicationOutput_tudat.h"
#include "getCurrentDateTime.h"
//#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
//#include <Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h>
//#include <Tudat/Astrodynamics/Ephemerides/ephemeris.h>
//#include <Tudat/JsonInterface/Propagation/variable.h>
//#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
//#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
//#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>
//#include <Tudat/Astrodynamics/Ephemerides/frameManager.h>

Space4ErrBody_Ballistic::Space4ErrBody_Ballistic(
        std::vector< std::vector< double > > &bounds,
        std::vector< double > &input_data,
        std::string &play):
    problemBounds_( bounds ),
    input_data_( input_data ),
    play_( play ){ }

//        const bool useExtendedDynamics) :
//    useExtendedDynamics_( useExtendedDynamics ){ }

//! Descriptive name of the problem
std::string Space4ErrBody_Ballistic::get_name() const
{
    return "AMS to IAD Ballistic Trajectory: 10/28/2018  11:00:00 AM";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > Space4ErrBody_Ballistic::get_bounds() const
{
    return { problemBounds_[0], problemBounds_[1] };
}

//! Implementation of the fitness function
//! For this ballistic case, the fitness function returns delta, which is a
//! vector containing differences in GOAL and calculated latitudes and
//! longitudes, their "norm", and time o flight.
//std::vector<double> Space4ErrBody_Ballistic::fitness(const std::vector<double> &x) const
std::vector<double> Space4ErrBody_Ballistic::fitness( const std::vector< double > &x )  const
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
    using namespace tudat::input_output;
    using namespace tudat::unit_conversions;
    using namespace tudat;
    using namespace pagmo;
    using namespace tudat_applications;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Set simulation start epoch.
    const double simulationStartEpoch = input_data_[0]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Set simulation end epoch.
    // const double simulationEndEpoch = simulationStartEpoch + 18000.0; // 5 hours
    const double simulationEndEpoch = simulationStartEpoch + input_data_[1];

    //! Set numerical integration fixed step size.
    const double fixedStepSize = input_data_[2];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" },
                                    simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );

    //! This is an acceptable 'cheat' were Earth is placed at the barycenter.
    //! Use only when there arent any third body perturbations (Moon, Sun, etc.)
    bodySettings[ "Earth" ]->ephemerisSettings =
            boost::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    double radius_Earth = spice_interface::getAverageRadius( "Earth" );
    double flattening_Earth = 1 / 298.257223563; //! https://www.vcalc.com/wiki/vCalc/WGS-84+Earth+flattening+factor
    bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< OblateSphericalBodyShapeSettings >( radius_Earth, flattening_Earth );

    //! Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    //! Create vehicle objects.
    bodyMap[ "HORUS" ] = boost::make_shared< simulation_setup::Body >( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Create vehicle aerodynamic coefficients. Gotta get my own for HORUS.
    bodyMap[ "HORUS" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );

    //! Create aerodynamic coefficient interface settings.
    double referenceArea = input_data_[3];//4.0*000000001;
    //double aerodynamicCoefficient = input_data_[4];//1.2 * 0.0000000000001;
   // boost::shared_ptr< simulation_setup::AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
     //       boost::make_shared< ConstantAerodynamicCoefficientSettings >(
       //         referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    //! Create and set aerodynamic coefficients object
   // bodyMap[ "HORUS" ]->setAerodynamicCoefficientInterface(
     //           createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "HORUS" ) );

    //! Set body Mass.
    bodyMap[ "HORUS" ]->setConstantBodyMass( input_data_[5] );

    //! Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    //! Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfHORUS;

    //! Gravitational model was stated to be Cartesian. I havent found this one
    //! in the available models. This means it must be developed, and this here
    //! modified.
    accelerationsOfHORUS[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    //! Aerodynamic accelerations are attached differently than gravitational. Why?
    accelerationsOfHORUS[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    //! Create Acceleration map from the accelerationsOfHORUS data map.
    accelerationMap[  "HORUS" ] = accelerationsOfHORUS;

    //! Define bodies that will be propagated. Only 1, HORUS.
    bodiesToPropagate.push_back( "HORUS" );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( "Earth" );

    //! Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap,
                accelerationMap,
                bodiesToPropagate,
                centralBodies );

    //! Define constant 30 degree angle of attack - ARBITRARY. Shouldn't this be
    //! a function of time or events? Would that require reading from a file?
    //! Would that require multiple propagations that are stitched together?
    double constantAngleOfAttack = 30.0 * mathematical_constants::PI / 180.0;
    bodyMap.at( "HORUS" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::lambda::constant( constantAngleOfAttack ) );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Unpack input data values.
    const double h_i         = input_data_[6]; //0 * 1E3; // m
    const double lat_i_deg   = input_data_[7]; //52.30805556; //52deg 18’29"N
    const double lon_i_deg   = input_data_[8]; //4.76416667 //4deg 45’51"E
    const double h_f         = input_data_[9]; //0 * 1E3; // m
    const double lat_f_deg   = input_data_[10]; //38.9444444444444; //38deg 56’40"N
    const double lon_f_deg   = input_data_[11]; //-77.4558333333 //77deg 27’21"W
    const double v_i         = x[0]; // OPTIMIZE IT!
    const double gamma_i_deg = x[1]; // OPTIMIZE IT!
    const double chi_i_deg   = x[2]; // OPTIMIZE IT!

    //! Convert angles from degrees to radians
    const double lat_i_rad   = unit_conversions::convertDegreesToRadians( lat_i_deg );
    const double lon_i_rad   = unit_conversions::convertDegreesToRadians( lon_i_deg );
    const double lat_f_rad   = unit_conversions::convertDegreesToRadians( lat_f_deg );
    const double lon_f_rad   = unit_conversions::convertDegreesToRadians( lon_f_deg );
    const double gamma_i_rad = unit_conversions::convertDegreesToRadians( gamma_i_deg );
    const double chi_i_rad   = chi_i_deg * mathematical_constants::PI / 180;

    //! Calculate initial heading angle: https://www.movable-type.co.uk/scripts/latlong.html
    //std::atan2( std::sin( lon_f_rad - lon_i_rad ) * std::cos( lat_f_rad ) , std::cos( lat_i_rad ) * std::sin( lat_f_rad ) - std::sin( lat_i_rad ) * std::cos( lat_f_rad ) * std::cos( lon_f_rad - lon_i_rad ) );

    //! Set spherical elements for HORUS initial state.
    //! ARTBITRARILY CHOSEN to be in Earth-Fixed frame.
    Eigen::Vector6d horusSphericalEntryState;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex )       = radius_Earth + h_i;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex )     = lat_i_rad;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex )    = lon_i_rad;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex )        = v_i;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex )   = gamma_i_rad;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = chi_i_rad;//0.6 * 0;

    //! Set spherical elements for HORUS final state.
    //! ARTBITRARILY CHOSEN to be in Earth-Fixed frame.
    Eigen::Vector6d horusSphericalFINALState = horusSphericalEntryState;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::radiusIndex )       = radius_Earth + h_f;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::latitudeIndex )     = lat_f_rad;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::longitudeIndex )    =  lon_f_rad;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::speedIndex )        = 0;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::flightPathIndex )   = 0;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0;//0.6 * 0;

    //! Convert HORUS state vector from spherical to Cartesian elements
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                horusSphericalEntryState );
    Eigen::Vector6d systemFinalStateGOAL = convertSphericalOrbitalToCartesianState(
                horusSphericalFINALState );

    //! Create Earth's rotational ephemeris.
    boost::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( "Earth" )->getRotationalEphemeris( );

    //! Transform Earth-fixed frame to Inertial Frame.
    systemInitialState = transformStateToGlobalFrame( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

    //! Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable,
                    "HORUS" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    "HORUS",
                    "Earth" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic,
                    "HORUS",
                    "Earth",
                    1 ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable,
                    "HORUS" ) );

    // Create object with list of dependent variables
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


    // Define termination conditions
    boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable,
                "HORUS",
                "Earth" );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationDependentVariableTerminationSettings >(
                terminationDependentVariable,
                h_f,
                true );

    //! Create propagation settings.
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies,
              accelerationModelMap,
              bodiesToPropagate,
              systemInitialState,
              terminationSettings,
              cowell);
       //       dependentVariablesToSave );
    boost::shared_ptr< IntegratorSettings<  > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              fixedStepSize );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< double > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings );


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
    Eigen::Vector6d systemFinalState =  ( dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin() )->second;

    //! Calculate Time of Flight
    const double tof = simulationEndEpoch_calc - simulationStartEpoch;

    //! Transform final state from Inertial Frame to Earth-Fixed Frame
    Eigen::Vector6d systemFinalState_EARTH_FIXED = transformStateToTargetFrame( systemFinalState, simulationEndEpoch_calc, earthRotationalEphemeris );

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
   const double altitude_f_calc = std::sqrt( pow(systemFinalState_EARTH_FIXED[0],2) +
           pow(systemFinalState_EARTH_FIXED[1],2) + pow(systemFinalState_EARTH_FIXED[2],2) ) ;
   double lon_f_rad_calc = std::atan2(systemFinalState_EARTH_FIXED[1] , systemFinalState_EARTH_FIXED[0]);
   double lat_f_rad_calc = std::asin(systemFinalState_EARTH_FIXED[2] / altitude_f_calc) ;

   //! Calculate offset of final state from GOAL state: Earth-Fixed Frame
   const double dif_lat_rad = lat_f_rad - lat_f_rad_calc;
   const double dif_lon_rad = lon_f_rad - lon_f_rad_calc;

   //! Convert offsets of final state to degrees: Earth-Fixed Frame
   const double dif_lat_deg = dif_lat_rad * 180 / mathematical_constants::PI;
   const double dif_lon_deg = dif_lon_rad * 180 / mathematical_constants::PI;
   //const double dif_lon_rad = unit_conversions::convertDegreesToRadians( dif_lon_deg );
   //const double dif_lat_rad = unit_conversions::convertDegreesToRadians( dif_lat_deg );

   //! Calculate "norm" of offsets. This is an arbitrary function I have
   //! implemented to pass on as an 'objective function'. It relates the
   //! differences such that when minimizing the offsets, there is an additional
   //! unsigned value that always goes to zero. Most definitely unsure about how
   //!  'proper' it is, yet is what works for the current BALLISTIC case.
   const double dif_norm = std::sqrt(dif_lat_deg * dif_lat_deg + dif_lon_deg * dif_lon_deg);

   //! Assign values to Fitness vector! At the moment these are all 'objective
   //! functions'. No constraints have been implemented.
   std::vector< double > delta;
   delta.push_back(dif_norm);
   delta.push_back(dif_lat_deg);
   delta.push_back(dif_lon_deg);
   delta.push_back(tof);  // Not sure yet how this one affects the optimization. Included for completion.

   //! Print results to terminal. Used to gauge progress.
   std::cout << std::fixed << std::setprecision(10) <<
                  std::setw(7) << "tof = " <<
                  std::setw(16) << tof <<
                  std::setw(7) << "v_i = " <<
                  std::setw(16) << v_i <<
                  std::setw(15) << "gamma_i_deg = " <<
                  std::setw(14) << gamma_i_deg <<
                  std::setw(13) << "chi_i_deg = " <<
                  std::setw(16) << chi_i_deg <<
                  std::setw(12) << "latitude = " <<
                  std::setw(16) << lat_f_rad_calc * 180 / mathematical_constants::PI  <<
                  std::setw(13) << "longitude = " <<
                  std::setw(16) << lon_f_rad_calc * 180 / mathematical_constants::PI <<
                  std::setw(19) << "latitude offset = " <<
                  std::setw(16) << dif_lat_deg <<
                  std::setw(20) << "longitude offset = " <<
                  std::setw(16) << dif_lon_deg <<std::endl;

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //! Define output folder for THIS particular set of runs
   std::string outputSubFolder = "HORUS_OUTPUT_" + std::to_string(int(input_data_[18])) + "_" + std::to_string(int(input_data_[19])) + "_" + std::to_string(int(input_data_[20])) +  "_" + std::to_string(fixedStepSize) + "_" + play_ +  "/";

   //! Get time stamp for this specific simulation. This avoids overwriting the
   //! file if another individual with the same properties shows up in other
   //! evolutions.
   std::string simulation_save_time = getCurrentDateTime( false );

   //! Write HORUS propagation history to file.
   writeDataMapToTextFile( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                           "HORUSPropHistory_" + std::to_string(v_i) + "_" + std::to_string(gamma_i_deg) + "_" + std::to_string(chi_i_deg) + "_" + std::to_string(tof) + "_" + simulation_save_time + ".dat",
                           tudat_applications::getOutputPath( ) + outputSubFolder,
                           "",
                           std::numeric_limits< double >::digits10,
                           std::numeric_limits< double >::digits10,
                           "," );

 /*  writeMatrixToFile( systemFinalStateGOAL,
                           "HORUSsystemFinalStateGOAL_" + std::to_string(v_i) + "_" + std::to_string(gamma_i_deg) + "_" + std::to_string(chi_i_deg) + ".dat",
                      std::numeric_limits< double >::digits10,
                           tudat_applications::getOutputPath( ) + outputSubFolder,
                           "");
   writeMatrixToFile( systemFinalState_EARTH_FIXED,
                           "HORUSsystemFinalState_EARTH_FIXED" + std::to_string(v_i) + "_" + std::to_string(gamma_i_deg) + "_" + std::to_string(chi_i_deg) + ".dat",
                      std::numeric_limits< double >::digits10,
                           tudat_applications::getOutputPath( ) + outputSubFolder,
                           "");
*/
//}

   return delta;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    // Use Boost library to make a random number generator.
    boost::mt19937 randomNumbergenerator( time( 0 ) );
    boost::random::uniform_real_distribution< > uniformDistribution( 0.0, 10.0 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumbers( randomNumbergenerator, uniformDistribution );

    // Generate random altitude value between 0 and 10 km.
    //Artbitraty altitude
    const double altitudeKilometers = 25; //generateRandomNumbers( );

    // Use the Tudat Core library to convert km to m.
    const double altitudeMeters = tudat::unit_conversions::convertKilometersToMeters(
            altitudeKilometers );

    // Use the Eigen library to create position vectors.
    Eigen::Vector3d positionOfBodySubjectToAcceleration;
    positionOfBodySubjectToAcceleration << 6378136.6 + altitudeMeters, 0.0, 0.0;

    // Position vector of Central body exerting acceleration.
    const Eigen::Vector3d positionOfBodyExertingAcceleration = Eigen::Vector3d::Zero( );

    // Use the Tudat library to compute the acceleration vector.
    const Eigen::Vector3d gravitationalAcceleration =
            tudat::gravitation::computeGravitationalAcceleration(
                    positionOfBodySubjectToAcceleration, 4.9e12,
                    positionOfBodyExertingAcceleration );
   */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

