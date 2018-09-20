#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/filesystem/operations.hpp>

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
//#include "getFitness.h"
#include "updateGuidance.h"


//namespace tudat
//{
std::vector< double > getFitness(
        const std::vector< double > &input_data_,
        const std::vector< double > &output_settings_,
        const std::string &outputSubFolder_,
        const std::vector< double > &x)
{

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Reading 'using' statements" << std::endl;
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
    using namespace tudat;
  //  using namespace pagmo;
    using namespace tudat_applications;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK SOME DATA              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Unpacking data" << std::endl;

    //! Unpack input data values.
    //! Declare variables to unpack from input_data. This list may include
    //! variables that may not be used further.
    double Ref_area;
    double M_i;
    double v_i;
    double gamma_i_deg;
    double chi_i_deg;
    double AoA_deg;
    double h_i;
    double lat_i_deg;
    double lon_i_deg;
    double h_f;
    double lat_f_deg;
    double lon_f_deg;

    if ( int(input_data_[0]) == 0 )
    {
        //! Assign Reference area
        Ref_area    = input_data_[4]; //m^2

        //! Assign initial mass
        M_i         = input_data_[6]; // kg

        //! Assign the initial and final position conditions.
        h_i         = input_data_[7]; //0 * 1E3; // m
        lat_i_deg   = input_data_[8]; //52.30805556; //52deg 18’29"N
        lon_i_deg   = input_data_[9]; //4.76416667 //4deg 45’51"E
        h_f         = input_data_[10]; //0 * 1E3; // m
        lat_f_deg   = input_data_[11]; //38.9444444444444; //38deg 56’40"N
        lon_f_deg   = input_data_[12]; //282.544166666667 //77deg 27’21"W

        //! Variables to optimize
        v_i         = x[0]; // OPTIMIZE IT!
        gamma_i_deg = x[1]; // OPTIMIZE IT!
        chi_i_deg   = x[2]; // OPTIMIZE IT!
        AoA_deg     = x[3]; // OPTIMIZE IT!
        const double AoA_rad     = unit_conversions::convertDegreesToRadians( AoA_deg );

    }
    else
    {
        //! Assign Reference area
        Ref_area    = input_data_[4]; // ??? m^2 // dont have it yet

        //! Assign initial mass
        M_i         = input_data_[6]; //??? * 1E3 kg // dont have it yet

        //! Assign the initial and final position conditions.
        h_i         = input_data_[10]; //122 * 1E3 m // given
        lat_i_deg   = input_data_[11]; //-22.37 deg // given
        lon_i_deg   = input_data_[12]; //-106.7 deg // given
        h_f         = input_data_[13]; //25 * 1E3 m // arbitrary?
        lat_f_deg   = input_data_[14]; //5.237222 deg //5 deg 14’14"N // arbitrary?
        lon_f_deg   = input_data_[15]; //-52.760556 deg //52 deg 45’38"W// arbitrary?

        //! Variables to optimize
        v_i         = x[0]; // OPTIMIZE IT!
        gamma_i_deg = x[1]; // OPTIMIZE IT!
        chi_i_deg   = x[2]; // OPTIMIZE IT!

    }
    /*
    const double h_i         = input_data_[6]; //0 * 1E3; // m
    const double lat_i_deg   = input_data_[7]; //52.30805556; //52deg 18’29"N
    const double lon_i_deg   = input_data_[8]; //4.76416667 //4deg 45’51"E
    const double h_f         = input_data_[9]; //0 * 1E3; // m
    const double lat_f_deg   = input_data_[10]; //38.9444444444444; //38deg 56’40"N
    const double lon_f_deg   = input_data_[11]; //-77.4558333333 //77deg 27’21"W
    const double v_i         = x[0]; // OPTIMIZE IT!
    const double gamma_i_deg = x[1]; // OPTIMIZE IT!
    const double chi_i_deg   = x[2]; // OPTIMIZE IT!
*/
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
    ///////////////////////            DEFINE TIMEFRAME & TIMESTEP            /////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   std::cout << "Defining timeframe and timestep" << std::endl;

    //! Set simulation start epoch.
    const double simulationStartEpoch = input_data_[1]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Set simulation end epoch.
    // const double simulationEndEpoch = simulationStartEpoch + 18000.0; // 5 hours
    const double simulationEndEpoch = simulationStartEpoch + input_data_[2];

    //! Set numerical integration fixed step size.
    const double fixedStepSize = input_data_[3];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Creating environment" << std::endl;

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
    //!
    //! Per Dominic:
    //! The best way to do this, is to define geodetic coordinates (altitude,
    //!  geodetic latitude, longitude), which you can then transform to
    //! Cartesian elements (or any other set) to define your initial condition.
    //!  This will work for both a spherical and flattened Earth model.
    const double radius_Earth = spice_interface::getAverageRadius( "Earth" );


    //!---------------------------- Still playing around with this section

    double radius_Earth_i = radius_Earth;
    double radius_Earth_f = radius_Earth;
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
   // std::cout << "Creating vehicle" << std::endl;

    //! Assign reference area
    const double referenceArea = input_data_[3];//4.0*000000001;

    //! Assign constant mass
    const double vehicleMass = input_data_[5];//4.0*000000001;
    //double aerodynamicCoefficient = input_data_[4];//1.2 * 0.0000000000001;
    //boost::shared_ptr< simulation_setup::AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
    //       boost::make_shared< ConstantAerodynamicCoefficientSettings >(
    //         referenceArea,
    //         aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

      ///////// Start: Vehicle Aerodynamics Section

    //! Define physical meaning of independent variables
    //!     Mach number
    //!     Angle of attack
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );
    independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );

    //! Define list of files for force coefficients. Entry denotes direction.
    //!     0 : x-direction (C ~D~/C ~X~)
    //!     1 : y-direction (C ~S~/C ~Y~)
    //!     2 : z-direction (C ~L~/C ~Z~)
    //! Gotta get my own for HORUS.
    std::map< int, std::string > forceCoefficientFiles;
    forceCoefficientFiles[ 0 ] = "aurora_CD.txt"; // Set drag coefficient file
    forceCoefficientFiles[ 2 ] = "aurora_CL.txt"; // Set lift coefficient file

    //! Define reference frame in which the loaded coefficients are defined.
    //! Have to get some more background info here to properly understand it.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;

    //! Load and parse coefficient files; create coefficient settings.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
        readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles,
                referenceArea,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame,
                areCoefficientsInNegativeAxisDirection );
  //  std::cout << "Creating vehicle: Aerodynamic coefficients are done" << std::endl;

    //! Define constant angle of attack - ARBITRARY. Shouldn't this be
    //! a function of time or events? Would that require reading from a file?
    //! Would that require multiple propagations that are 'stitched' together?
   // const double constantAngleOfAttack = AoA_deg * mathematical_constants::PI / 180.0;
    //bodyMap.at( "HORUS" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ) );

    ///////// End: Vehicle Aerodynamics Section

    //! Create vehicle objects.
    bodyMap[ "HORUS" ] = boost::make_shared< simulation_setup::Body >( );

    //! Set body Mass.
    bodyMap[ "HORUS" ]->setConstantBodyMass( vehicleMass );

    //! Set vehicle aerodynamic coefficients.
    bodyMap[ "HORUS" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "HORUS" ) );


    //! Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Creating accelerations" << std::endl;

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

    //! Aerodynamic accelerations are attached differently than gravitational. Why?
    accelerationsOfHORUS[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    //! Create Acceleration map from the accelerationsOfHORUS data map.
    accelerationMap[ "HORUS" ] = accelerationsOfHORUS;

    //! Define bodies that will be propagated. Only 1, HORUS.
    bodiesToPropagate.push_back( "HORUS" );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( "Earth" );

    //! Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationMap,
                bodiesToPropagate,
                centralBodies );

    //! Guidance is set AFTER the accelerations and BEFORE propagating. Placed
    //! both statements here becuase I am unsure which one matters. Haven't
    //! investigated further, but it works.

    //! Create aerodynamic guidance. Does this have to be done AFTER acceleration models?
    boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance =
            boost::make_shared< FlightConditionsBasedAerodynamicGuidance >(bodyMap, "HORUS");
    //std::cout << "Creating vehicle: Guidance is done AFTER acceleration models" << std::endl;

    //boost::shared_ptr< aerodynamics::FlightConditions > flightConditions_ = double currentFlightPathAngle = flightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );

    //! Set Guidance angle functions.
    setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap.at( "HORUS" ) );

   // std::cout << "Creating vehicle: Guidance is set" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Creating propagation settings" << std::endl;


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
    //dependentVariablesList.push_back(
    //            boost::make_shared< SingleDependentVariableSaveSettings >(
    //                aerodynamic_force_coefficients_dependent_variable,
    //                "HORUS" ) );

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
              cowell,
              dependentVariablesToSave );
    boost::shared_ptr< IntegratorSettings<  > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              fixedStepSize );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE                  ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Propagating" << std::endl;

    //! Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< double > dynamicsSimulator(
                bodyMap,
                integratorSettings,
                propagatorSettings );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        CLEAN UP AND CALCULATE FITNESS                //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Clean up and calculate fitness" << std::endl;

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

   //! Assign values to Fitness vector! At the moment these are all 'objective
   //! functions'. No constraints have been implemented. To modify this I have
   //! change the header file and define how many are elements are objective
   //! function, equality contraints, and inequlity constraints. This vector
   //! here must contain them is that exact order: nOF, nEC, nIC.
   std::vector< double > delta;
   delta.push_back(dif_norm);
   delta.push_back(dif_lat_deg);
   delta.push_back(dif_lon_deg);
   delta.push_back(tof);  // Not sure yet how this one affects the optimization. Included for completion.

   //! Print results to terminal. Used to gauge progress.
   if ( int(input_data_[0]) == 0 )
{
       std::cout << std::fixed << std::setprecision(10) <<
                  std::setw(7) << "tof = " <<
                  std::setw(16) << tof <<
                  std::setw(7) << "v_i = " <<
                  std::setw(16) << v_i <<
                  std::setw(15) << "flight-path = " <<
                  std::setw(14) << gamma_i_deg <<
                  std::setw(13) << "heading = " <<
                  std::setw(16) << chi_i_deg <<
                  std::setw(7) << "AoA = " <<
                  std::setw(16) << AoA_deg <<
                  std::setw(7) << "lat = " <<
                  std::setw(16) << lat_f_deg_calc <<
                  std::setw(7) << "lon = " <<
                  std::setw(16) << lon_f_deg_calc <<
                  std::setw(14) << "lat offset = " <<
                  std::setw(16) << dif_lat_deg <<
                  std::setw(14) << "lon offset = " <<
                  std::setw(16) << dif_lon_deg <<std::endl;
   }
   else
   {
       std::cout << std::fixed << std::setprecision(10) <<
                  std::setw(7) << "tof = " <<
                  std::setw(16) << tof <<
                  std::setw(7) << "v_i = " <<
                  std::setw(16) << v_i <<
                  std::setw(15) << "flight-path = " <<
                  std::setw(14) << gamma_i_deg <<
                  std::setw(13) << "heading = " <<
                  std::setw(16) << chi_i_deg <<
                  std::setw(7) << "AoA = " <<
                  std::setw(16) << AoA_deg <<
                  std::setw(7) << "lat = " <<
                  std::setw(16) << lat_f_deg_calc <<
                  std::setw(7) << "lon = " <<
                  std::setw(16) << lon_f_deg_calc <<
                  std::setw(14) << "lat offset = " <<
                  std::setw(16) << dif_lat_deg <<
                  std::setw(14) << "lon offset = " <<
                  std::setw(16) << dif_lon_deg <<std::endl;
   }

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
}
//}// namespace tudat
