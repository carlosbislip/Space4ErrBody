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
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>



//! Do I need a costumized one for HORUS? Probably.... Maybe not. I believe it
//! would be sufficient with the aerodynamic coefficient tables.
//#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>

//! Mine
#include "Space4ErrBody.h"
#include "applicationOutput_tudat.h"
#include "getStuff.h"
#include "updateGuidance.h"
#include "updateGuidance_val.h"
#include "StopOrNot_val.h"
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
        const std::vector< double > &input_data,
        const std::vector< double > &output_settings,
        const std::string &outputSubFolder):
    problemBounds_( bounds ),
    input_data_( input_data ),
    output_settings_( output_settings ),
    outputSubFolder_( outputSubFolder ){ }

//        const bool useExtendedDynamics) :
//    useExtendedDynamics_( useExtendedDynamics ){ }

//! Descriptive name of the problem
std::string Space4ErrBody::get_name() const
{
    if ( int(input_data_[0]) == 0 )
    {
        return "AMS to IAD Ballistic Trajectory: 10/28/2018  11:00:00 AM";
    }
    else
    {
        return "Validation: HORUS Entry towards Kourou";
    }
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > Space4ErrBody::get_bounds() const
{
    return { problemBounds_[0], problemBounds_[1] };
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




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  std::cout << "Unpacking data" << std::endl;

    //! Set simulation start epoch.
    const double simulationStartEpoch = input_data_[1]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Set simulation end epoch.
    // const double simulationEndEpoch = simulationStartEpoch + 18000.0; // 5 hours
    const double simulationEndEpoch = simulationStartEpoch + input_data_[2];

    //! Set numerical integration fixed step size.
    const double fixedStepSize = input_data_[3];

    //! Assign Reference area
    const double Ref_area = input_data_[4]; //m^2

    //! Assign initial mass
    const double M_i = input_data_[6]; // kg

    //! Declare various variables to unpack from input_data depending on a
    //! conditional.
    double v_i, gamma_i_deg, chi_i_deg, AoA_deg, AoA_rad;
    double h_i, lat_i_deg,lon_i_deg;
    double h_f, lat_f_deg,lon_f_deg;

    if ( int(input_data_[0]) == 0 )
    {
        //! Assign the initial and final position conditions.
        h_i         = input_data_[7]; //0 * 1E3; // m
        lat_i_deg   = input_data_[8]; //52.30805556; //52deg 18’29"N
        lon_i_deg   = input_data_[9]; //4.76416667 //4deg 45’51"E
        h_f         = input_data_[10]; //0 * 1E3; // m
        lat_f_deg   = input_data_[11]; //38.9444444444444; //38deg 56’40"N
        lon_f_deg   = input_data_[12]; //282.544166666667 //77deg 27’21"W

        //! Variables to optimize
        v_i         = x[0]; // FIND IT!
        gamma_i_deg = x[1]; // FIND IT!
        chi_i_deg   = x[2]; // FIND IT!
        AoA_deg     = x[3]; // FIND IT!
        AoA_rad     = unit_conversions::convertDegreesToRadians( AoA_deg );
    }
    else
    {
        //! Assign the initial and final position conditions.
        h_i         = input_data_[10]; //122 * 1E3 m // given
        lat_i_deg   = input_data_[11]; //-22.37 deg // given
        lon_i_deg   = input_data_[12]; //-106.7 deg // given
        h_f         = input_data_[13]; //25 * 1E3 m //
        lat_f_deg   = input_data_[14]; //5.0 deg //5 deg N //
        lon_f_deg   = input_data_[15]; //-53 deg //53 deg W//

        //! Variables to optimize
        v_i         = x[0]; // FIND IT!
        gamma_i_deg = x[1]; // FIND IT!
        chi_i_deg   = x[2]; // FIND IT!
    }

    //! Convert angles from degrees to radians
    const double lat_i_rad   = unit_conversions::convertDegreesToRadians( lat_i_deg );
    const double lon_i_rad   = unit_conversions::convertDegreesToRadians( lon_i_deg );
    const double lat_f_rad   = unit_conversions::convertDegreesToRadians( lat_f_deg );
    const double lon_f_rad   = unit_conversions::convertDegreesToRadians( lon_f_deg );
    const double gamma_i_rad = unit_conversions::convertDegreesToRadians( gamma_i_deg );
    const double chi_i_rad   = unit_conversions::convertDegreesToRadians( chi_i_deg );

    //! Calculate initial heading angle: https://www.movable-type.co.uk/scripts/latlong.html
    //std::atan2( std::sin( lon_f_rad - lon_i_rad ) * std::cos( lat_f_rad ) , std::cos( lat_i_rad ) * std::sin( lat_f_rad ) - std::sin( lat_i_rad ) * std::cos( lat_f_rad ) * std::cos( lon_f_rad - lon_i_rad ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //   std::cout << "Creating environment" << std::endl;

    //! Declare/define simulation body settings data map.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" },
                                    simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );

    //! Define atmospheric model.
    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >(
                nrlmsise00 );

    //! Define ephemeris model settings.
    //! This is an acceptable 'cheat' were Earth is placed at the barycenter.
    //! Use only when there arent any third body perturbations (Moon, Sun, etc.)
    bodySettings[ "Earth" ]->ephemerisSettings =
            boost::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    //! Reset ephemeris to J2000.
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    //! Define Earth's radius. Using spice here. Is there a way to get a 'radius
    //! field'? Im interested in includig Earth's flattening, yet am unsure how
    //! to properly do it.
    const double radius_Earth = spice_interface::getAverageRadius( "Earth" );
    double radius_Earth_i = radius_Earth;
    double radius_Earth_f = radius_Earth;

    //!---------------------------- Still playing around with this section

    //! https://en.wikipedia.org/wiki/Flattening#Numerical_values_for_planets
    //! https://www.vcalc.com/wiki/vCalc/WGS-84+Earth+flattening+factor
    double flattening_Earth = 1 / input_data_.back();
    if (flattening_Earth != 1)
    {
        bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< OblateSphericalBodyShapeSettings >( radius_Earth, flattening_Earth );

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
        const double radius_Earth_f_x = radii[ 0 ] * std::cos(polar_angle_f_rad) * std::cos(lon_f_rad);
        const double radius_Earth_f_y = radii[ 1 ] * std::cos(polar_angle_f_rad) * std::sin(lon_f_rad);
        const double radius_Earth_f_z = radii[ 2 ] * std::sin(polar_angle_f_rad);
        radius_Earth_f = std::sqrt(pow(radius_Earth_f_x,2) + pow(radius_Earth_f_y,2) + pow(radius_Earth_f_z,2));

    }
    //!--------------------------------------------------------------

    //! Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    //! Create Earth's rotational ephemeris.
    boost::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( "Earth" )->getRotationalEphemeris( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    std::cout << "Creating vehicle" << std::endl;

    //! Assign reference area
//    const double referenceArea = input_data_[3];//4.0*000000001;

    //! Assign constant mass
//    const double vehicleMass = input_data_[5];//4.0*000000001;
    //double aerodynamicCoefficient = input_data_[4];//1.2 * 0.0000000000001;
    //boost::shared_ptr< simulation_setup::AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
    //       boost::make_shared< ConstantAerodynamicCoefficientSettings >(
    //         referenceArea,
    //         aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

 //   std::cout << "Creating vehicle Aerodynamics" << std::endl;

      ///////// Start: Vehicle Aerodynamics Section

    //! Define physical meaning of independent variables
    //!     Mach number
    //!     Angle of attack
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );

    //! Define list of files for force coefficients. Entry denotes direction.
    //!     0 : x-direction (C ~D~/C ~X~)
    //!     1 : y-direction (C ~S~/C ~Y~)
    //!     2 : z-direction (C ~L~/C ~Z~)
    //! Gotta get my own for HORUS.
    std::map< int, std::string > forceCoefficientFiles;
    forceCoefficientFiles[ 0 ] = "HORUS_CD.txt"; // Set drag coefficient file
    forceCoefficientFiles[ 2 ] = "HORUS_CL.txt"; // Set lift coefficient file

    //! Define reference frame in which the loaded coefficients are defined.
    //! Have to get some more background info here to properly understand it.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;

    //! Load and parse coefficient files; create coefficient settings.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
        readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles,
                Ref_area,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame,
                areCoefficientsInNegativeAxisDirection );

    //! Define constant angle of attack - ARBITRARY. Shouldn't this be
    //! a function of time or events? Would that require reading from a file?
    //! Would that require multiple propagations that are 'stitched' together?
   // const double constantAngleOfAttack = AoA_deg * mathematical_constants::PI / 180.0;
    //bodyMap.at( "HORUS" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ) );


    ///////// End: Vehicle Aerodynamics Section

    //std::cout << "Creating vehicle objects" << std::endl;

    //! Create vehicle objects.
    bodyMap[ "HORUS" ] = boost::make_shared< simulation_setup::Body >( );

    //! Set body Mass.
    bodyMap[ "HORUS" ]->setConstantBodyMass( M_i );

    //! Set vehicle aerodynamic coefficients.
    bodyMap[ "HORUS" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "HORUS" ) );

    //! Finalize body creation. Not entirely sure what this does. I believe it
    //! may place the body (Earth) at the barycenter. Probably another hack to
    //! facilitate bringing in everything into the inertial frame.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //   std::cout << "Creating accelerations" << std::endl;

    //! Declare propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    //! Declare acceleration data map.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfHORUS;

    //! Define gravitational model. Arbitrary aximum degree/order. According to
    //! Dominic, equivalent functionality to Cartesian with corresponding maximum
    //! degree/order.
    accelerationsOfHORUS[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    //! Define aerodynamic accelerations.
    //! Aerodynamic accelerations are attached differently than gravitational. Why?
    accelerationsOfHORUS[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    //! Assign acceleration map from the accelerationsOfHORUS data map.
    accelerationMap[ "HORUS" ] = accelerationsOfHORUS;

    //! Define bodies that will be propagated. Only 1, HORUS.
    bodiesToPropagate.push_back( "HORUS" );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( "Earth" );

    //! Set acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationMap,
                bodiesToPropagate,
                centralBodies );


    //std::cout << "Creating vehicle: Guidance is set AFTER acceleration models" << std::endl;

    //! Guidance is set AFTER the accelerations and BEFORE propagating.

    //! Set initial coordinates. Earth-Fixed.
    bodyMap[ "HORUS" ]->setInitialLat( lat_i_rad );
    bodyMap[ "HORUS" ]->setInitialLon( lon_i_rad );

    //! Set target coordinates. Earth-Fixed.
    bodyMap[ "HORUS" ]->setTargetLat( lat_f_rad );
    bodyMap[ "HORUS" ]->setTargetLon( lon_f_rad );

    //! Calculate initial distance to target.
    //const double initial_d_to_target_rad = std::acos( std::sin(lat_i_rad) * std::sin(lat_f_rad) + std::cos(lat_i_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_i_rad) );
    //const double initial_d_to_target_rad = getAngularDistance(lat_i_rad,lon_i_rad,lat_f_rad,lon_f_rad);
    //const double initial_d_to_target_deg = unit_conversions::convertRadiansToDegrees( initial_d_to_target_rad );

    //! Set initial distance to target.
    bodyMap[ "HORUS" ]->setInitialDistanceToTarget( getAngularDistance(lat_i_rad,lon_i_rad,lat_f_rad,lon_f_rad) );

    //! Pass starting epoch to body.
    bodyMap[ "HORUS" ]->setStartingEpoch( simulationStartEpoch );

    //! Declare and assign aerodynamic guidance functions.
    boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance;
    if ( int(input_data_[0]) == 0 )
    {
        aerodynamicGuidance =
            boost::make_shared< MyAerodynamicGuidance >(bodyMap, "HORUS");
    }
    else
    {
        aerodynamicGuidance =
            boost::make_shared< ValidationAerodynamicGuidance >(
                    bodyMap,
                    "HORUS");
    }

    //! Set Guidance angle functions.
    setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap.at( "HORUS" ) );

    //! Define constant orientation
    //double constantAngleOfAttack = unit_conversions::convertDegreesToRadians( 30 );
    //double constantBankAngle = unit_conversions::convertDegreesToRadians( 85 );
    //bodyMap.at( "HORUS" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ),boost::lambda::constant( 0 ),boost::lambda::constant( constantBankAngle ));
    //std::cout << "Creating vehicle: Guidance is set" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    std::cout << "Creating propagation settings" << std::endl;

    //! Set spherical elements for HORUS initial state.
    //! ARTBITRARILY CHOSEN to be in Earth-Fixed frame.
    Eigen::Vector6d horusSphericalEntryState;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex )       = radius_Earth_i + h_i;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex )     = lat_i_rad;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex )    = lon_i_rad;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex )        = v_i;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex )   = gamma_i_rad;
    horusSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = chi_i_rad;

    //! Set spherical elements for HORUS final state.
    //! ARTBITRARILY CHOSEN to be in Earth-Fixed frame.
    Eigen::Vector6d horusSphericalFINALState = horusSphericalEntryState;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::radiusIndex )       = radius_Earth_f + h_f;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::latitudeIndex )     = lat_f_rad;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::longitudeIndex )    = lon_f_rad;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::speedIndex )        = 0;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::flightPathIndex )   = 0;
    horusSphericalFINALState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0;;

    //! Convert HORUS state vector from spherical to Cartesian elements
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                horusSphericalEntryState );
    Eigen::Vector6d systemFinalStateGOAL = convertSphericalOrbitalToCartesianState(
                horusSphericalFINALState );

    //! Transform INITIAL state from Earth-fixed frame to Inertial Frame.
    systemInitialState = transformStateToGlobalFrame(
                systemInitialState,
                simulationStartEpoch,
                earthRotationalEphemeris );

    //! Define list of dependent variables to save.
    //! The file that prints out the text saying what is saved has been modified
    //! to mute this message. During optimization this would clutter the
    //! terminal screen and is incredibly annoying.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_cartesian_position,
                    "HORUS",
                    "Earth" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_spherical_position,
                    "HORUS",
                    "Earth" ) );
    //dependentVariablesList.push_back(
    //            boost::make_shared< SingleDependentVariableSaveSettings >(
    //                mach_number_dependent_variable,
    //                "HORUS" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    "HORUS",
                    "Earth" ) );
    //dependentVariablesList.push_back(
    //            boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
    //                aerodynamic,
    //                "HORUS",
    //                "Earth",
    //                1 ) );
    dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "HORUS",
                    latitude_angle,
                    "Earth") );
    dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "HORUS",
                    longitude_angle,
                    "Earth") );
    dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "HORUS",
                    heading_angle,
                    "Earth") );
    dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "HORUS",
                    flight_path_angle,
                    "Earth") );
    dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "HORUS",
                    angle_of_attack,
                    "Earth") );
    dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "HORUS",
                    angle_of_sideslip,
                    "Earth") );
    dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "HORUS",
                    bank_angle,
                    "Earth") );

    // Create object with list of dependent variables
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    //! Define termination conditions
    //! Dominic mentioned something about terminating at exactly the termination
    //! condition. Somewhere in the documentation. This works for now, though
    //! the altitude temrination condition may yield a "nearest neighbor" or
    //! dnet error with nrlmsise00 as the atmospheric model.
    //!
    //! Per Dominic:
    //! If I remember correctly, these errors (dnet....) are thrown by NRLMSISE
    //! when the altitude becomes < 0. The good thing is that the propagation
    //! results up until the error should be saved (numerical solution and
    //! dependent variables), so you can verify whether the last step gets
    //! close to zero. Note that the very last step, where it crashes, won't be
    //! saved, so you may find that the results don't actually show it getting
    //! below zero.

    boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable,
                "HORUS",
                "Earth" );


    //In case your custom function requires more inputs (e.g., it may depend on the position of the spacecraft or other variables that are not the current time), you can use boost::bind to add more inputs.

    //As an example, the case where the state of the spacecraft is added as an input is shown below:

   // boost::function< Eigen::Vector6d( ) > HORUS_StateFunction =
     //       boost::bind( &Body::getState, bodyMap.at( "HORUS" ) );

    std::vector< double > term_cond;
    term_cond.push_back( unit_conversions::convertDegreesToRadians( 0.75 ) );
    term_cond.push_back( 25000 );

    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationCustomTerminationSettings >(
                boost::bind( &StopOrNot, bodyMap, "HORUS", term_cond ) );

    //boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
    //        boost::make_shared< PropagationDependentVariableTerminationSettings >(
    //            terminationDependentVariable,
    //            h_f,
    //            true );

    //! Create propagation settings.
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies,
              accelerationModelMap,
              bodiesToPropagate,
              systemInitialState,
              terminationSettings,
              cowell,
              dependentVariablesToSave );
    boost::shared_ptr< IntegratorSettings<  > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              fixedStepSize );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //   std::cout << "Propagating" << std::endl;

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
   const double lon_f_rad_calc = std::atan2(systemFinalState_EARTH_FIXED[1] , systemFinalState_EARTH_FIXED[0]);
   const double lat_f_rad_calc = std::asin(systemFinalState_EARTH_FIXED[2] / altitude_f_calc) ;

   //! Convert coordinates of final state to degrees: Earth-Fixed Frame
   const double lon_f_deg_calc = unit_conversions::convertRadiansToDegrees( lon_f_rad_calc );
   const double lat_f_deg_calc = unit_conversions::convertRadiansToDegrees( lat_f_rad_calc );

   //! Calculate angular distance of final state from target coordinates.
   const double d_rad = getAngularDistance( lat_f_rad_calc, lon_f_rad_calc, lat_f_rad, lon_f_rad );
   const double d_deg = unit_conversions::convertRadiansToDegrees( d_rad );

   //! Calculate offset of final angular distance to termination condition distance.
   const double dif_d_deg = d_deg - unit_conversions::convertDegreesToRadians( term_cond[0] );

   //! Calculate offset of final state from GOAL state: Earth-Fixed Frame
   const double dif_lat_rad = lat_f_rad - lat_f_rad_calc;
   const double dif_lon_rad = lon_f_rad - lon_f_rad_calc;

   //! Calculate offsets of final state in degrees: Earth-Fixed Frame
   const double dif_lat_deg = lat_f_deg - lat_f_deg_calc;
   const double dif_lon_deg = lon_f_deg - lon_f_deg_calc;

   //! Calculate "norm" of offsets. This is an arbitrary function I have
   //! implemented to pass on as an 'objective function'. It relates the
   //! differences such that when minimizing the offsets, there is an additional
   //! unsigned value that always goes to zero. Most definitely unsure about how
   //!  'proper' it is, yet is what works for the current BALLISTIC case.
   const double dif_norm = std::sqrt( pow(dif_lat_deg,2) + pow(dif_lon_deg,2) );

   //! Calculate difference of final angular distance to
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
   delta.push_back( tof );  // Not sure yet how this one affects the optimization. Included for completion.

   //! Print results to terminal. Used to gauge progress.
   std::cout << std::fixed << std::setprecision(10) <<
                  std::setw(7) << "tof = " <<
                  std::setw(16) << tof <<
                  std::setw(7) << "v_i = " <<
                  std::setw(16) << v_i <<
                  std::setw(15) << "flight-path = " <<
                  std::setw(14) << gamma_i_deg <<
                  std::setw(13) << "heading = " <<
                  std::setw(16) << chi_i_deg <<
                  std::setw(7) << "lat = " <<
                  std::setw(16) << lat_f_deg_calc <<
                  std::setw(7) << "lon = " <<
                  std::setw(16) << lon_f_deg_calc <<
                  std::setw(14) << "lat offset = " <<
                  std::setw(16) << dif_lat_deg <<
                  std::setw(14) << "lon offset = " <<
                  std::setw(16) << dif_lon_deg <<
                  std::setw(17) << "angular dist. = " <<
                  std::setw(16) << d_deg <<std::endl;
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //! Get time stamp for this specific simulation. This avoids overwriting the
   //! file if another individual with the same properties shows up in other
   //! evolutions.
   std::string simulation_save_time = getCurrentDateTime( false );

   //! Create unique filename that cannot be overwritten due to the timestamp.
   std::string simulation_file_name_suffix = std::to_string(v_i) + "_" +
           std::to_string(gamma_i_deg) + "_" +
           std::to_string(chi_i_deg) + "_" +
           std::to_string(tof) + "_" +
           std::to_string(x[3]) + "_" +
           std::to_string(lat_f_deg_calc) + "_" +
           std::to_string(lon_f_deg_calc) + "_" +
           simulation_save_time;

   //! Will print out depending on some input values. Each entry corresponds to
   //! a different type of output. Entries are binary, 0 or 1.
   if (output_settings_[1] == 1)
   {

   //! Write HORUS propagation history to file.
   writeDataMapToTextFile( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                           "HORUSPropHistory_" + simulation_file_name_suffix + ".dat",
                           tudat_applications::getOutputPath( ) + outputSubFolder_,
                           "",
                           std::numeric_limits< double >::digits10,
                           std::numeric_limits< double >::digits10,
                           "," );
   }
   if (output_settings_[2] == 1)
   {
   //! Write HORUS dependent variables' history to file.
   writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory( ),
                           "HORUSDepVar_" + simulation_file_name_suffix + ".dat",
                           tudat_applications::getOutputPath( ) + outputSubFolder_,
                           "",
                           std::numeric_limits< double >::digits10,
                           std::numeric_limits< double >::digits10,
                           "," );
   }
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

   return delta;

} // Fitness function.
