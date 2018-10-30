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
        const std::vector< double > &problem_name,
        const std::vector< double > &vehicleParameterValues,
        const std::vector< std::string > &aeroCoeffFileList,
        const std::vector< double > &simulation_settingsValues,
        const std::vector< double > &initialConditionsValues,
        const std::vector< double > &terminationConditionsValues,
        const std::vector< double > &output_settingsValues,
        const std::string &outputSubFolder ):
    problemBounds_( bounds ),
    problem_name_( problem_name ),
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
    const double S_ref = input_data_[4]; //m^2

    //! Assign initial mass
    const double M_i = input_data_[6]; // kg

    //! Declare various variables to unpack from input_data depending on a
    //! conditional.
    double v_i, gamma_i_deg, chi_i_deg;//, AoA_deg, AoA_rad;
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
        //AoA_deg     = x[3]; // FIND IT!
        //AoA_rad     = unit_conversions::convertDegreesToRadians( AoA_deg );
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

    //! Assign central body name.
    const std::string centralBodyName = "Earth";

    //! Declare/define simulation body settings data map.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { centralBodyName },
                                    simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );

    //! Define atmospheric model.
    bodySettings[ centralBodyName ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >(
                nrlmsise00 );

    //! Define ephemeris model settings.
    //! This is an acceptable 'cheat' were Earth is placed at the barycenter.
    //! Use only when there arent any third body perturbations (Moon, Sun, etc.)
    bodySettings[ centralBodyName ]->ephemerisSettings =
            boost::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    //! Reset ephemeris to J2000.
    bodySettings[ centralBodyName ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    //! Define Earth's radius. Using spice here. Is there a way to get a 'radius
    //! field'? Im interested in includig Earth's flattening, yet am unsure how
    //! to properly do it.
    const double radius_Earth = spice_interface::getAverageRadius( centralBodyName );
    double radius_Earth_i = radius_Earth;
    double radius_Earth_f = radius_Earth;

    //!---------------------------- Still playing around with this section

    //! https://en.wikipedia.org/wiki/Flattening#Numerical_values_for_planets
    //! https://www.vcalc.com/wiki/vCalc/WGS-84+Earth+flattening+factor
    double flattening_Earth = 1 / input_data_.back();
    if (flattening_Earth != 1)
    {
        bodySettings[ centralBodyName ]->shapeModelSettings = boost::make_shared< OblateSphericalBodyShapeSettings >( radius_Earth, flattening_Earth );

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
            bodyMap.at( centralBodyName )->getRotationalEphemeris( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    std::cout << "Creating vehicle" << std::endl;
    //! Assign vehicle name.
    const std::string vehicleName = "HORUS";

    //! Assign coefficient file names.
    const std::string dragCoefficients = "HORUS_CD.txt";
    const std::string liftCoefficients = "HORUS_CL.txt";
   // const std::string dragControlSurfaceForceCoefficients = "HORUS_CD_CS.txt";
   // const std::string liftControlSurfaceForceCoefficients = "HORUS_CL_CS.txt";
    const std::string momentCoefficients = "HORUS_Cm.txt";
   // const std::string controlSurfaceMomentCoefficients = "HORUS_CM_CS.txt";
    double x_mrc = 13;
    double b_ref = 13;
    double c_ref = 13;

    Eigen::Vector3d R_mrc( x_mrc, 0, 0 );


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
    //! Gotta get my own for HORUS.
    std::map< int, std::string > forceCoefficientFiles;
   // std::map< int, std::string > controlSurfaceForceCoefficientFiles;
    std::map< int, std::string > momentCoefficientFiles;
   // std::map< int, std::string > controlSurfaceMomentCoefficientFiles;

    forceCoefficientFiles[ 0 ] = dragCoefficients; // Set drag coefficient file
    forceCoefficientFiles[ 2 ] = liftCoefficients; // Set lift coefficient file
   // controlSurfaceForceCoefficientFiles[ 0 ] = dragControlSurfaceForceCoefficients; // Set drag coefficient increment file
   // controlSurfaceForceCoefficientFiles[ 2 ] = liftControlSurfaceForceCoefficients; // Set lift coefficient increment file
    momentCoefficientFiles[ 1 ] = momentCoefficients; // Set pitch moment coefficient file
   // controlSurfaceMomentCoefficientFiles[ 1 ] = controlSurfaceMomentCoefficients; // Set pitch moment coefficient increment file

//       boost::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > controlSurfaceInterface =
//              boost::make_shared< CustomControlSurfaceIncrementAerodynamicInterface >(
//                   &dummyControlIncrements,
//                   boost::assign::list_of( angle_of_attack_dependent )( control_surface_deflection_dependent ) );
//       std::map< std::string, boost::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
//       controlSurfaceList[ "bodyflap" ] = controlSurfaceInterface;

    //! Define reference frame in which the loaded coefficients are defined.
    //! Have to get some more background info here to properly understand it.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;


    //! Load and parse coefficient files; create coefficient settings.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
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
   // boost::shared_ptr< system_models::VehicleSystems > systemsModels = boost::make_shared< system_models::VehicleSystems >( );

    //bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ) );

    ///////// End: Vehicle Aerodynamics Section

    //! Create vehicle objects.
    bodyMap[ vehicleName ] = boost::make_shared< simulation_setup::Body >( );

    //! Set body Mass.
    bodyMap[ vehicleName ]->setConstantBodyMass( M_i );

    //! Pass Earth's rotation rate
    bodyMap[ vehicleName ]->setCentralBodyRotationRate( 7.292115*1E-5 );

    //! Set initial coordinates. Earth-Fixed.
    bodyMap[ vehicleName ]->setInitialLat( lat_i_rad );
    bodyMap[ vehicleName ]->setInitialLon( lon_i_rad );

    //! Set target coordinates. Earth-Fixed.
    bodyMap[ vehicleName ]->setTargetLat( lat_f_rad );
    bodyMap[ vehicleName ]->setTargetLon( lon_f_rad );

    //! Set initial distance to target.
    bodyMap[ vehicleName ]->setInitialDistanceToTarget(
                getAngularDistance(
                    lat_i_rad,
                    lon_i_rad,
                    lat_f_rad,
                    lon_f_rad) );

    //! Pass starting epoch to body.
    bodyMap[ vehicleName ]->setStartingEpoch( simulationStartEpoch );

    //! Set vehicle aerodynamic coefficients.
    bodyMap[ vehicleName ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface(
                    aerodynamicCoefficientSettings,
                    vehicleName ) );

    //! Set vehicle system models for aerodynamic coefficient increments.
    //bodyMap[ vehicleName ]->setVehicleSystems( systemsModels );
    //bodyMap[ vehicleName ]->getAerodynamicCoefficientInterface( )->setControlSurfaceIncrements( controlSurfaceList );

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
    EntryState_spherical( SphericalOrbitalStateElementIndices::speedIndex )        = v_i;
    EntryState_spherical( SphericalOrbitalStateElementIndices::flightPathIndex )   = gamma_i_rad;
    EntryState_spherical( SphericalOrbitalStateElementIndices::headingAngleIndex ) = chi_i_rad;

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

    //! Declare propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    //! Declare acceleration data map.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > vehicleAccelerations;

    //! Define gravitational model. Arbitrary aximum degree/order. According to
    //! Dominic, equivalent functionality to Cartesian with corresponding maximum
    //! degree/order.
    vehicleAccelerations[ centralBodyName ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    //! Define aerodynamic accelerations.
    //! Aerodynamic accelerations are attached differently than gravitational. Why?
    vehicleAccelerations[ centralBodyName ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    //! Assign acceleration map from the vehicleAccelerations data map.
    accelerationMap[ vehicleName ] = vehicleAccelerations;

    //! Define bodies that will be propagated. Only 1.
    bodiesToPropagate.push_back( vehicleName );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( centralBodyName );

    //! Set acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationMap,
                bodiesToPropagate,
                centralBodies );

    //! Guidance is set AFTER the accelerations and BEFORE propagating.
    //! Declare and assign aerodynamic guidance functions.
    boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance;
    if ( int( input_data_[ 0 ] ) == 0 )
    {
        aerodynamicGuidance =
            boost::make_shared< MyAerodynamicGuidance >(
                    bodyMap,
                    vehicleName);
    }
    else
    {
        aerodynamicGuidance =
            boost::make_shared< ValidationAerodynamicGuidance >(
                    bodyMap,
                    vehicleName);
    }

    //! Set Guidance angle functions.
    setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap.at( vehicleName ) );

    //! Define constant orientation
    //double constantAngleOfAttack = unit_conversions::convertDegreesToRadians( 30 );
    //double constantBankAngle = unit_conversions::convertDegreesToRadians( 85 );
    //bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ),boost::lambda::constant( 0 ),boost::lambda::constant( constantBankAngle ));
    //std::cout << "Creating vehicle: Guidance is set" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    std::cout << "Creating propagation settings" << std::endl;

    //! Define list of dependent variables to save.
    //! The file that prints out the text saying what is saved has been modified
    //! to mute this message. During optimization this would clutter the
    //! terminal screen and is incredibly annoying.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dep_varList;
    dep_varList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_cartesian_position,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_spherical_position,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable,
                    vehicleName ) );
    dep_varList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic,
                    vehicleName,
                    centralBodyName,
                    1 ) );
    //dep_varList.push_back(
    //            boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
    //                vehicleName,
    //                latitude_angle,
    //                centralBodyName) );
    //dep_varList.push_back(
    //            boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
    //                vehicleName,
    //                longitude_angle,
    //                centralBodyName) );
    dep_varList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    heading_angle,
                    centralBodyName) );
    dep_varList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    flight_path_angle,
                    centralBodyName) );
    dep_varList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    angle_of_attack,
                    centralBodyName) );
    dep_varList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    angle_of_sideslip,
                    centralBodyName) );
    dep_varList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    bank_angle,
                    centralBodyName) );

    // Create object with list of dependent variables
    boost::shared_ptr< DependentVariableSaveSettings > dep_varToSave =
            boost::make_shared< DependentVariableSaveSettings >( dep_varList );

    //! Define termination conditions
    //! Dominic mentioned something about terminating at exactly the termination
    //! condition. Somewhere in the documentation. This works for now, though
    //! the altitude temrination condition may yield a "nearest neighbor" or
    //! dnet error with nrlmsise00 as the atmospheric model.
    //!
    //! Per Dominic:
    //! If I remember correctly, these errors (dnet....) are thrown by NRLMSISE
    //! when the altitude becomes out of range. The good thing is that the propagation
    //! results up until the error should be saved (numerical solution and
    //! dependent variables), so you can verify whether the last step gets
    //! close to zero. Note that the very last step, where it crashes, won't be
    //! saved, so you may find that the results don't actually show it getting
    //! below zero.

    /* boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
             boost::make_shared< SingleDependentVariableSaveSettings >(
                 altitude_dependent_variable,
                 vehicleName,
                 centralBodyName );

    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationDependentVariableTerminationSettings >(
                terminationDependentVariable,
                h_f,
                true );
*/

    //In case your custom function requires more inputs (e.g., it may depend on the position of the spacecraft or other variables that are not the current time), you can use boost::bind to add more inputs.
    //As an example, the case where the state of the spacecraft is added as an input is shown below:
   // boost::function< Eigen::Vector6d( ) > HORUS_StateFunction =
     //       boost::bind( &Body::getState, bodyMap.at( vehicleName ) );

    std::vector< double > term_cond;
    term_cond.push_back( unit_conversions::convertDegreesToRadians( 0.75 ) );
    term_cond.push_back( h_f );

    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationCustomTerminationSettings >(
                boost::bind( &StopOrNot, bodyMap, vehicleName, term_cond ) );


    //! Create propagation settings.
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies,
              accelerationModelMap,
              bodiesToPropagate,
              systemInitialState,
              terminationSettings,
              cowell,
              dep_varToSave );

    //! Create integrator settings.
    boost::shared_ptr< IntegratorSettings<  > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              fixedStepSize, 1, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    const double tof = simulationEndEpoch_calc - simulationStartEpoch;

    //! Transform final state from Inertial Frame to Earth-Fixed Frame
    //Eigen::Vector6d systemFinalState_EARTH_FIXED = transformStateToTargetFrame( systemFinalState, simulationEndEpoch_calc, earthRotationalEphemeris );

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
   //const double altitude_f_calc = dep_var_FINAL_STATE[3];
   const double lat_f_rad_calc  = dep_var_FINAL_STATE[4];
   const double lon_f_rad_calc  = dep_var_FINAL_STATE[5];
   const double h_f_calc        = dep_var_FINAL_STATE[6];
   //const double altitude_f_calc = std::sqrt( pow(systemFinalState_EARTH_FIXED[0],2) +
   //       pow(systemFinalState_EARTH_FIXED[1],2) + pow(systemFinalState_EARTH_FIXED[2],2) ) ;
   //const double lon_f_rad_calc = std::atan2(systemFinalState_EARTH_FIXED[1] , systemFinalState_EARTH_FIXED[0]);
   //const double lat_f_rad_calc = std::asin(systemFinalState_EARTH_FIXED[2] / altitude_f_calc) ;

   //! Convert coordinates of final state to degrees: Earth-Fixed Frame
   const double lat_f_deg_calc = unit_conversions::convertRadiansToDegrees( lat_f_rad_calc );
   const double lon_f_deg_calc = unit_conversions::convertRadiansToDegrees( lon_f_rad_calc );

   //! Calculate angular distance of final state from target coordinates.
   const double d_rad = getAngularDistance( lat_f_rad_calc, lon_f_rad_calc, lat_f_rad, lon_f_rad );
   const double d_deg = unit_conversions::convertRadiansToDegrees( d_rad );

   //! Calculate offset of final angular distance to termination condition distance.
   const double dif_d_deg = d_deg - unit_conversions::convertDegreesToRadians( term_cond[0] );

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

   //! Calculate offset from goal elevation.
   const double dif_h = h_f_calc - h_f;

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


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //! Get time stamp for this specific simulation. This avoids overwriting the
   //! file if another individual with the same properties shows up in other
   //! evolutions.
   std::string simulation_save_time = getCurrentDateTime( false );

   //! Create unique filename that cannot be overwritten due to the timestamp.
   std::string simulation_file_name_suffix =
           std::to_string(v_i) + "_" +
           std::to_string(gamma_i_deg) + "_" +
           std::to_string(chi_i_deg) + "_" +
           std::to_string(tof) + "_" +
           //std::to_string(x[3]) + "_" +
           std::to_string(lat_f_deg_calc) + "_" +
           std::to_string(lon_f_deg_calc) + "_" +
           simulation_save_time;

   std::string complete_file_name_Prop = "HORUSPropHistory_" + simulation_file_name_suffix + ".dat";
   std::string complete_file_name_DepVar = "HORUSDepVar_" + simulation_file_name_suffix + ".dat";

   //! Will print out depending on some input values. Each entry corresponds to
   //! a different type of output. Entries are binary, 0 or 1.
   if (output_settings_[1] == 1)
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
   if (output_settings_[2] == 1)
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



   //! Print results to terminal. Used to gauge progress.
   std::cout << std::fixed << std::setprecision(10) <<
                  std::setw(7) << "v_i = " <<
                  std::setw(16) << v_i <<
                  std::setw(15) << "flight-path = " <<
                  std::setw(14) << gamma_i_deg <<
                  std::setw(13) << "heading = " <<
                  std::setw(16) << chi_i_deg <<
                  std::setw(12) << "dif_norm = " <<
                  std::setw(16) << dif_norm <<
                  std::setw(14) << "dif_lat_deg = " <<
                  std::setw(16) << dif_lat_deg <<
                  std::setw(14) << "dif_lon_deg = " <<
                  std::setw(16) << dif_lon_deg <<
                  std::setw(13) << "dif_d_deg = " <<
                  std::setw(16) << dif_d_deg <<
                  std::setw(9) << "dif_h = " <<
                  std::setw(16) << dif_h <<
                  std::setw(7) << "tof = " <<
                  std::setw(16) << tof <<
                  std::setw(120) << simulation_file_name_suffix << std::endl;

   std::cout <<dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) << std::endl;




   return delta;

} // Fitness function.
