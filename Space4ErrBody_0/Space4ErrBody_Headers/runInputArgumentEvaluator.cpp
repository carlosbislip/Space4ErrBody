#include "runInputArgumentEvaluator.h"
#include "Space4ErrBodyProblem.h"
#include "getPagmoAlgorithm.h"
#include <Tudat/Bislip/bislipUtilities.h>
#include <Tudat/Bislip/bislipVariables.h>
#include <Tudat/Bislip/bislipConstraints.h>
#include <Tudat/Bislip/bislipVehicleParameters.h>


#include <Tudat/Bislip/bislipDecisionVectorEvaluation.h>
#include <Tudat/Bislip/createTUDATSettings.h>
#include <Tudat/Bislip/updateGuidance.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <pagmo/pagmo.hpp>

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
//using namespace tudat_applications;
using namespace pagmo;
using namespace bislip;
using namespace bislip::Variables;
using namespace bislip::parameters;
using namespace bislip::utilities;
using namespace boost;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Main Program                  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace bislip
{

void runInputArgumentEvaluator( const std::vector< std::string > &runInputArguments )
{
    // std::cout << "Declaring Problem Input pointer" << std::endl;
    std::shared_ptr< bislip::problem_input > problemInput = std::make_shared< bislip::problem_input >( );
    std::shared_ptr< bislip::problem_input > problemInputCopy;
    //! Get this run's time stamp. Will be used to create a unique output
    //! subfolder name, where all files created by the optimizer will be stored.
    std::chrono::time_point< std::chrono::system_clock > playTime = bislip::utilities::getDateTime( );
    std::string playTimeString = bislip::utilities::convertDateTimeToString( false, playTime );
    std::pair< std::chrono::time_point< std::chrono::system_clock >, std::string > playTimePair = std::make_pair( playTime, playTimeString );

    //std::cout << "Initializing Problem Input pointer" << std::endl;
    problemInput->setPlayTimePair( playTimePair );

    std::string inputCaseMatrixFileName = runInputArguments[ 0 ];
    std::string inputCaseMatrixFieldNamesFileName = runInputArguments[ 1 ];
    std::string caseToRun = runInputArguments[ 2 ];
    //std::string runDetails = runInputArguments[ 3 ];


    //std::string primaryFileName = runInputArguments[ 0 ];
    // std::string simulationType = runInputArguments[ 1 ];
    // std::string continuityType = runInputArguments[ 2 ];
    //std::string nodalParameterInputMatrixFileName = "NOTREQUIRED";

    // std::cout << "Read inputs from files present in: " << primaryFileName << std::endl;
    unsigned int p = 0;
    std::vector< std::string > primary;
    Eigen::MatrixXd densityParameterMatrixForJerk;
    std::vector< std::string > vehicleParameterData;
    std::vector< std::string > vehicleParameterList;
    std::vector< double > vehicleParameterValues;
    std::vector< std::string > aeroCoeffFileList;
    std::vector< std::string > headingErrorDeadBandFileList;
    Eigen::MatrixXd headingErrorDeadBandCoarse;
    Eigen::MatrixXd headingErrorDeadBandLowDistance;
    std::vector< std::string > alphaMachEnvelopeFileList;
    Eigen::MatrixXd alphaMachEnvelopeLB;
    Eigen::MatrixXd alphaMachEnvelopeUB;
    std::vector< std::string > simulationSettingsList;
    std::vector< double > simulationSettingsValues;
    std::vector< std::string > initialConditionsList;
    std::vector< double > initialConditionsValues;
    std::vector< std::string > constraintsList;
    std::vector< double > constraintsValues;
    std::vector< std::string > outputSettingsList;
    std::vector< double > outputSettingsValues;
    //Eigen::VectorXd  parameterBoundsValues;
    //std::vector< std::string > parameterList;
    std::vector< std::string > nonNodalParameterBoundsList;
    std::vector< std::string > optimizationSettingsList;
    Eigen::VectorXd optimizationSettings;
    Eigen::MatrixXd nonNodalParameterBounds;
    std::vector< std::string > nonNodalParameterBoundsData;
    std::vector< std::string >optimizerParameterList;
    Eigen::MatrixXd optimizerParameterCaseMatrix;
    std::vector< std::string > caseMatrix;
    std::vector< std::string > inputCaseMatrixFieldNames;


    unsigned int caseToRun_int = 100000000;
    std::string caseDetailsToRun_String;
    std::vector< std::string > caseDetailsToRun_Vector;
    std::vector< std::vector< std::string > > vehicleParameterData_Split;
    std::string mainSubFolder, simulationType, continuityType, problemName, vehicleName, trajectoryType, nodalParameterInputMatrixFileName, kourouTimeHistoriesFileName;
    Eigen::MatrixXd nodalParameterInputMatrix;
    Eigen::MatrixXd kourouTimeHistories;
    std::string resumeOptimizationFolderPath;
    bool importInitialPopulation = false;

    std::string initialPopulationFile;
    std::string initialFitnessFile;
    std::string initialExtremesAndConstraintsFile;

    bool validation = false;
    bool testCase = false;
    bool resumeOptimization = false;

    try
    {
        std::cout << "Loading Input Case Matrix: " << inputCaseMatrixFileName << std::endl;
        caseMatrix    = bislip::utilities::getDataString ( inputCaseMatrixFileName );

        try
        {
            std::cout << "  Case to Run: " << caseToRun << std::endl;
            caseToRun_int = static_cast< unsigned int >( std::stoi( caseToRun ) );
        }
        catch( const char *a )
        {
            throw std::runtime_error("Case string not convertible to integer.");
        }

        std::cout << "Extracting case details from Case Matrix" << std::endl;
        caseDetailsToRun_String = caseMatrix[ caseToRun_int ];

        std::cout << "  Splitting details and populating case details vector" << std::endl;
        boost::algorithm::split( caseDetailsToRun_Vector, caseDetailsToRun_String, boost::is_any_of( "\t" ),
                                 boost::algorithm::token_compress_on );

        std::cout << "  Assigning to variables and reading in data based on initial vector entries" << std::endl;
        mainSubFolder                     = caseDetailsToRun_Vector[ p ];
        simulationType                    = caseDetailsToRun_Vector[ p += 1 ];
        continuityType                    = caseDetailsToRun_Vector[ p += 1 ];
        resumeOptimizationFolderPath      = caseDetailsToRun_Vector[ p += 1 ];
        problemName                       = caseDetailsToRun_Vector[ p += 1 ];
        vehicleName                       = caseDetailsToRun_Vector[ p += 1 ];
        densityParameterMatrixForJerk     = tudat::input_output::readMatrixFromFile( caseDetailsToRun_Vector[ p += 1 ], " \t", "#" );
        vehicleParameterData              = bislip::utilities::getDataString( caseDetailsToRun_Vector[ p += 1 ] );
        aeroCoeffFileList                 = bislip::utilities::getDataString( caseDetailsToRun_Vector[ p += 1 ] );
        headingErrorDeadBandFileList      = bislip::utilities::getDataString( caseDetailsToRun_Vector[ p += 1 ] );
        alphaMachEnvelopeFileList         = bislip::utilities::getDataString( caseDetailsToRun_Vector[ p += 1 ] );
        nonNodalParameterBoundsData       = bislip::utilities::getDataString( caseDetailsToRun_Vector[ p += 1 ] );
        trajectoryType                    = caseDetailsToRun_Vector[ p += 1 ];
        nodalParameterInputMatrixFileName = caseDetailsToRun_Vector[ p += 1 ];
        kourouTimeHistoriesFileName       = caseDetailsToRun_Vector[ p += 1 ];
        importInitialPopulation           = bool( std::stoi( caseDetailsToRun_Vector[ p += 1 ] ) );
        initialPopulationFile             = caseDetailsToRun_Vector[ p += 1 ];
        initialFitnessFile                = caseDetailsToRun_Vector[ p += 1 ];
        initialExtremesAndConstraintsFile = caseDetailsToRun_Vector[ p += 1 ];

        std::cout << "  Splitting vehicle parameter data" << std::endl;
        std::vector< std::string > lineSplit_;
        std::vector< std::string > termsToRemove;
        termsToRemove.push_back( " " );
        std::string cleanLine;

        for ( unsigned int i = 0; i < vehicleParameterData.size(); i++ )
        {
            lineSplit_.clear( );
            cleanLine.clear();

            // Read current line and split into separate entries.
            boost::algorithm::split( lineSplit_, vehicleParameterData[ i ],
                                     boost::is_any_of( "," ),
                                     boost::algorithm::token_compress_on );

            cleanLine = lineSplit_[ 1 ];
            bislip::utilities::eraseSubStrings( cleanLine, termsToRemove );
            vehicleParameterList.push_back( lineSplit_[ 0 ] );
            vehicleParameterValues.push_back( std::stod( cleanLine ) );
        }
        nonNodalParameterBounds.resize( nonNodalParameterBoundsData.size(), 2 );

        std::cout << "  Populating static parameter bounds matrix" << std::endl;
        for ( unsigned int i = 0; i < nonNodalParameterBoundsData.size(); i++ )
        {
            lineSplit_.clear( );
            cleanLine.clear();

            // Read current line and split into separate entries.
            boost::algorithm::split( lineSplit_, nonNodalParameterBoundsData[ i ],
                                     boost::is_any_of( "," ),
                                     boost::algorithm::token_compress_on );

            nonNodalParameterBoundsList.push_back( lineSplit_[ 0 ] );

            cleanLine = lineSplit_[ 1 ];
            bislip::utilities::eraseSubStrings( cleanLine, termsToRemove );
            nonNodalParameterBounds( i, 0 ) = std::stod( cleanLine );

            cleanLine.clear();
            cleanLine = lineSplit_[ 2 ];
            bislip::utilities::eraseSubStrings( cleanLine, termsToRemove );
            nonNodalParameterBounds( i, 1 ) = std::stod( cleanLine );

        }

        headingErrorDeadBandCoarse        = tudat::input_output::readMatrixFromFile( headingErrorDeadBandFileList[ 0 ] , " \t", "#" );
        headingErrorDeadBandLowDistance   = tudat::input_output::readMatrixFromFile( headingErrorDeadBandFileList[ 1 ] , " \t", "#" );
        alphaMachEnvelopeLB               = tudat::input_output::readMatrixFromFile( alphaMachEnvelopeFileList[ 0 ] , " \t", "#" );
        alphaMachEnvelopeUB               = tudat::input_output::readMatrixFromFile( alphaMachEnvelopeFileList[ 1 ] , " \t", "#" );

        std::cout << "Simulation Type: " << simulationType << std::endl;
        if( simulationType == "thesis"  || simulationType == "Thesis" || simulationType == "THESIS" )
        {
            std::cout << "  Nodal Parameter Input Matrix File Name: " << nodalParameterInputMatrixFileName << std::endl;
            try
            {
                std::cout << "  Loading Nodal Parameter Input Matrix" << std::endl;
                nodalParameterInputMatrix = tudat::input_output::readMatrixFromFile( nodalParameterInputMatrixFileName , " \t", "#" );
            }
            catch( const char *something )
            {
                std::string strFromChar;

                std::string text = "Could not open " + nodalParameterInputMatrixFileName;
                std::string exception = "Exception: " + strFromChar.append( something );
                throw std::runtime_error( text + exception );
            }
        }
        else if( simulationType == "validation"  || simulationType == "Validation" || simulationType == "VALIDATION" )
        {
            validation = true;
            std::cout << "  Kourou Time Histories File Name: " << kourouTimeHistoriesFileName << std::endl;
            //std::cout << "  Kourou Bank Angle History File Name: " << primary.rbegin()[ 2 ] << std::endl;
            kourouTimeHistories       = tudat::input_output::readMatrixFromFile( kourouTimeHistoriesFileName, " \t", "#" );
        }
        else if( simulationType == "test"  || simulationType == "Test" || simulationType == "TEST" )
        {
            testCase = true;
            std::cout << "  Nodal Parameter Input Matrix File Name: " << nodalParameterInputMatrixFileName << std::endl;
            try
            {
                std::cout << "  Loading Nodal Parameter Input Matrix" << std::endl;
                nodalParameterInputMatrix = tudat::input_output::readMatrixFromFile( nodalParameterInputMatrixFileName , " \t", "#" );
            }
            catch( const char *something )
            {
                std::string strFromChar;

                std::string text = "Could not open " + nodalParameterInputMatrixFileName;
                std::string exception = "Exception: " + strFromChar.append( something );
                throw std::runtime_error( text + exception );
            }
        }

        //! If 'resume' type, where is the previous data?
        if( continuityType == "resume"  || continuityType == "Resume" || continuityType == "RESUME" )
        {
            resumeOptimization = true;
            std::cout << "Folder containing case to resume: " << resumeOptimizationFolderPath << std::endl;
        }

        std::cout << "  Erasing initial entries from case details vector" << std::endl;
        caseDetailsToRun_Vector.erase( caseDetailsToRun_Vector.begin(), caseDetailsToRun_Vector.begin() + p + 1 );

        inputCaseMatrixFieldNames    = bislip::utilities::getDataString ( inputCaseMatrixFieldNamesFileName );
        inputCaseMatrixFieldNames.erase( inputCaseMatrixFieldNames.begin(), inputCaseMatrixFieldNames.begin() + p  );

        std::cout << "  Extracting remaining (numeric) subvectors from case details vector and populating relevant vectors" << std::endl;
        simulationSettingsList   = bislip::utilities::extractContiguousSubVector( inputCaseMatrixFieldNames, 0, 13 );
        simulationSettingsValues = bislip::utilities::convertVectorOfStringsToVectorOfDoubles( bislip::utilities::extractContiguousSubVector( caseDetailsToRun_Vector, 0, 13 ) );
        outputSettingsList       = bislip::utilities::extractContiguousSubVector( inputCaseMatrixFieldNames, 13, 3 );
        outputSettingsValues     = bislip::utilities::convertVectorOfStringsToVectorOfDoubles( bislip::utilities::extractContiguousSubVector( caseDetailsToRun_Vector, 13, 3 ) );
        initialConditionsList    = bislip::utilities::extractContiguousSubVector( inputCaseMatrixFieldNames, 16, 7 );
        initialConditionsValues  = bislip::utilities::convertVectorOfStringsToVectorOfDoubles( bislip::utilities::extractContiguousSubVector( caseDetailsToRun_Vector, 16, 7 ) );
        constraintsList          = bislip::utilities::extractContiguousSubVector( inputCaseMatrixFieldNames, 23, 16 );
        constraintsValues        = bislip::utilities::convertVectorOfStringsToVectorOfDoubles( bislip::utilities::extractContiguousSubVector( caseDetailsToRun_Vector, 23, 16 ) );
        //parameterBoundsList      = bislip::utilities::extractContiguousSubVector( inputCaseMatrixFieldNames, 39, 18 );
        //parameterBoundsValues    = bislip::utilities::convertSTLVectorToEigenVectorXd( bislip::utilities::convertVectorOfStringsToVectorOfDoubles( bislip::utilities::extractContiguousSubVector( caseDetailsToRun_Vector, 39, 18 ) ) );
        optimizationSettingsList = bislip::utilities::extractContiguousSubVector( inputCaseMatrixFieldNames, 39, 18 );
        optimizationSettings     = bislip::utilities::convertSTLVectorToEigenVectorXd( bislip::utilities::convertVectorOfStringsToVectorOfDoubles( bislip::utilities::extractContiguousSubVector( caseDetailsToRun_Vector, 39, 18 ) ) );
    }
    catch( const char * something )
    {
        std::string strFromChar;

        std::string text = "Could not read file in position: " + std::to_string( p ) + " . ";
        std::string exception = "Exception: " + strFromChar.append( something );
        throw std::runtime_error( text + exception );
    }

    std::cout << "Set various Problem Inputs"<< std::endl;
    problemInput->setProblemName( problemName );
    problemInput->setOutputSettings( outputSettingsValues );
    problemInput->setSimulationSettings( simulationSettingsValues );
    problemInput->setVehicleName( vehicleName );
    problemInput->setVehicleParameterList( vehicleParameterList );
    problemInput->setVehicleParameters( vehicleParameterValues );
    problemInput->setAerodynamicDatabaseFileList( aeroCoeffFileList );
    problemInput->setNodalParameterInputMatrixStructure( nodalParameterInputMatrix );
    problemInput->setInitialConditions( initialConditionsValues );

    std::map< bislip::BislipVehicleSystems::body_names, std::string > bodyNamesMap;
    bodyNamesMap[ bislip::BislipVehicleSystems::body_names::vehicle_name ] = vehicleName;

    //std::cout << "Creating containers for the pinted population/fitness" << std::endl;
    std::map < std::string, Eigen::VectorXd > printedPopulation, printedFitness;

    Eigen::VectorXd dummy( 2 );
    dummy << 0, 0;
    /*printedPopulation[ std::to_string( 0 ) ] = dummy;
    printedPopulation.erase( std::to_string( 0 )  );
    printedFitness[ std::to_string( 0 ) ] = dummy;
    printedFitness.erase( std::to_string( 0 )  );

    //std::cout << "Passing containers to 'PROBLEM INPUT'" << std::endl;
    problemInput->setPrintedPopulation( printedPopulation );
    problemInput->setPrintedFitness( printedFitness );*/

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Declare and initialize simulation start epoch.
    const int simulationStartEpoch_CALENDAR_YEAR   = int( simulationSettingsValues[ 0 ] ); // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000
    const int simulationStartEpoch_CALENDAR_MONTH  = int( simulationSettingsValues[ 1 ] ); // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000
    const int simulationStartEpoch_CALENDAR_DAY    = int( simulationSettingsValues[ 2 ] ); // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000
    const int simulationStartEpoch_CALENDAR_HOUR   = int( simulationSettingsValues[ 3 ] ); // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000
    const int simulationStartEpoch_CALENDAR_MINUTE = int( simulationSettingsValues[ 4 ] ); // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000
    const int simulationStartEpoch_CALENDAR_SECOND = int( simulationSettingsValues[ 5 ] ); // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    const double simulationStartEpoch = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                simulationStartEpoch_CALENDAR_YEAR,
                simulationStartEpoch_CALENDAR_MONTH,
                simulationStartEpoch_CALENDAR_DAY,
                simulationStartEpoch_CALENDAR_HOUR,
                simulationStartEpoch_CALENDAR_MINUTE,
                simulationStartEpoch_CALENDAR_SECOND );

    //! Declare and initialize maximum simulation time of flight.
    const double max_tof = simulationSettingsValues[ 6 ];

    //! Declare and initialize simulation end epoch.
    const double simulationEndEpoch = simulationStartEpoch + max_tof;

    //! Declare and initialize numerical integration fixed step size.
    double propagationStepSize = simulationSettingsValues[ 7 ];

    //! Declare and initialize propagation to guidance sampling ratio.
    double guidanceStepSize = simulationSettingsValues[ 8 ];

    const unsigned int sphericalHarmonicDegree = static_cast< unsigned int >( simulationSettingsValues[ 9 ] );

    const unsigned int sphericalHarmonicOrder = static_cast< unsigned int >( simulationSettingsValues[ 10 ] );

    problemInput->setSimulationSettings( simulationSettingsValues );

    const bool maximizeThrust = int( simulationSettingsValues.rbegin()[ 1 ] );
     unsigned int debugInfo = static_cast< unsigned int >( simulationSettingsValues.rbegin()[ 0 ] );

    //! Declare and initialize Reference area
    const double S_ref = vehicleParameterValues[ 0 ]; // m^2

    //! Declare and initialize Reference length: wingspan
    const double b_ref = vehicleParameterValues[ 1 ]; // m

    //! Declare and initialize Reference length: chord
    const double c_ref = vehicleParameterValues[ 2 ]; // m

    Eigen::VectorXd referenceValues( 3 ); // m^2, m, m
    referenceValues( bislip::BislipVehicleSystems::reference_values::surface_area ) = S_ref;
    referenceValues( bislip::BislipVehicleSystems::reference_values::wingspan )     = b_ref;
    referenceValues( bislip::BislipVehicleSystems::reference_values::chord_length ) = c_ref;

    //! Declare and initialize position vector of moment reference center
    const Eigen::Vector3d R_mrc( vehicleParameterValues[ 3 ], vehicleParameterValues[ 4 ], vehicleParameterValues[ 5 ] ); // m

    //! Declare and initialize position vector of center of mass
    const Eigen::Vector3d R_com( vehicleParameterValues[ 6 ], vehicleParameterValues[ 7 ], vehicleParameterValues[ 8 ] ); // m

    //! Declare and initialize position vector of center of thrust
    const Eigen::Vector3d R_cot( vehicleParameterValues[ 9 ], vehicleParameterValues[ 10 ], vehicleParameterValues[ 11 ] ); // m

    //! Declare and initialize initial mass
    //const double initialMass_Ascent = vehicleParameterValues[ 12 ]; // kg

    //! Declare and initialize landing mass
    const double landingMass = vehicleParameterValues[ 13 ]; // kg

    //! Declare and initialize specific impulse
    //const double specificImpulse = vehicleParameterValues[ 14 ]; // kg

    //! Declare and initialize maximum engine thrust
    //const double maxThrust = vehicleParameterValues[ 15 ]; // N

    //! Declare and initialize relative chest forward angle
    //const double relativeChestForwardAngle_deg = vehicleParameterValues[ 16 ]; // deg

    //! Declare and initialize vertebral column inclination angle
    //const double vertebralColumnInclinationAngle_deg = vehicleParameterValues[ 17 ]; // deg

    //! Declare and initialize the minimum dynamic pressure required for bank angle assignment.
    const double minimumDynamicPressureforControlSurface = vehicleParameterValues[ 18 ];

    //! Declare and initialize nose radius
    const double noseRadius = vehicleParameterValues[ 19 ]; // m

    //! Declare and initialize nose radius
    //const double leadingEdgeRadius = vehicleParameterValues[ 20 ]; // m

    //! Declare and initialize nose radius
    const double wallEmissivity = vehicleParameterValues[ 21 ];

    //! Declare and initialize Bodyflap deflection bounds
    const double minimumBodyflapAngle = vehicleParameterValues[ 22 ];
    const double maximumBodyflapAngle = vehicleParameterValues[ 23 ];

    //! Declare and initialize Elevon deflection bounds
    const double minimumElevonAngle = vehicleParameterValues[ 24 ];
    const double maximumElevonAngle = vehicleParameterValues[ 25 ];

    //! Declare and initialize Thrust Elevation bounds
    const double minimumThrustElevationAngle = vehicleParameterValues[ 26 ];
    const double maximumThrustElevationAngle = vehicleParameterValues[ 27 ];
    //const double minimumThrustAzimuthAngle = vehicleParameterValues[ 28 ];
    //const double maximumThrustAzimuthAngle = vehicleParameterValues[ 29 ];

    //! Declare and initialize Throttle Setting bounds
    const double minimumThrottleSetting = vehicleParameterValues[ 30 ];
    const double maximumThrottleSetting = vehicleParameterValues[ 31 ];


    //! Declare and initialize starting position coordinates.
    const double initialLat_deg = initialConditionsValues[ 0 ];
    const double initialLon_deg = initialConditionsValues[ 1 ];

    //! Declare and initialize starting height
    const double initialHeight = initialConditionsValues[ 2 ]; // m

    //! Declare and initialize initial Mach
    const double initialMachNumber = initialConditionsValues[ 3 ]; // -

    //! Declare and initialize initial flight-path angle
    const double initialFlightPathAngle_deg = initialConditionsValues[ 4 ]; // deg

    std::vector< bislip::constraints::constraint_list > constraintList = bislip::constraints::getConstraintList();
    for( unsigned int i = 0; i < constraintList.size(); i++ )
    { problemInput->setConstraintValue( constraintList[ i ], constraintsValues[ i ] );}

    //! Convert angles from degrees to radians
    const double initialLat_rad                      = unit_conversions::convertDegreesToRadians( initialLat_deg );
    const double initialLon_rad                      = unit_conversions::convertDegreesToRadians( initialLon_deg );
    const double targetLat_rad                       = unit_conversions::convertDegreesToRadians( problemInput->getConstraintValue( bislip::constraints::constraint_list::final_latitude ) );
    const double targetLon_rad                       = unit_conversions::convertDegreesToRadians( problemInput->getConstraintValue( bislip::constraints::constraint_list::final_longitude ) );
    const double finalDistanceToTarget_rad           = unit_conversions::convertDegreesToRadians( problemInput->getConstraintValue( bislip::constraints::constraint_list::termination_angular_distance_to_go ) );
    const double initialFlightPathAngle_rad          = unit_conversions::convertDegreesToRadians( initialFlightPathAngle_deg );

    Eigen::VectorXd centralTargetCoordinates( 2 );
    centralTargetCoordinates( 0 ) = targetLat_rad;
    centralTargetCoordinates( 1 ) = targetLon_rad;

    //! Pre-define various variables used to determine fitness.
    //double targetLat_deg_calc          = initialLat_deg;
    //double targetLon_deg_calc          = initialLon_deg;
    double initialDistanceToTarget_rad = bislip::Variables::computeAngularDistance( initialLat_rad, initialLon_rad, targetLat_rad, targetLon_rad );
    //double initialDistanceToTarget_deg = unit_conversions::convertRadiansToDegrees( initialDistanceToTarget_rad );

    //! Still working on these
    //const double noseRadius = 3.0;
    //const double leadingEdgeRadius = 1.5;
    const double lambda = unit_conversions::convertDegreesToRadians( 30.0 );
    const double x_T = 1E-9;
    const double phi = unit_conversions::convertDegreesToRadians( 3.0 );

    std::cout << "Creating partial output subfolder filename suffix" << std::endl;
    std::cout << "    Extracting run details from input arguments: " << runInputArguments[ runInputArguments.size() - 1 ] << std::endl;

    /*if( debugInfo == 1 ){ std::cout << "    Extracting optimizer case from input arguments: " << runInputArguments[ runInputArguments.size() - 1 ] << std::endl; }
    const std::string optimizerParameterCase_string = runInputArguments[ runInputArguments.size() - 1 ];
    int optimizerParameterCase_int = 10000000;
    try
    {
        optimizerParameterCase_int = std::stoi( optimizerParameterCase_string );
    }
    catch(...)
    {
        throw("Optimizer case string not convertible to integer.");
    }*/

    //if( debugInfo == 1 ){ std::cout << "    Get optimizer settings from Case Matrix" << std::endl; }
    // Eigen::VectorXd optimizationSettings = optimizerParameterCaseMatrix.row( optimizerParameterCase_int );

    std::cout << "    Store optimizer parameters" << std::endl;
    std::vector< bislip::optimization_settings::optimization_settings_list > optimizationSettingsEnumList = bislip::optimization_settings::getOptimizationSettingsList();
    for( unsigned int i = 0; i < optimizationSettingsEnumList.size(); i++ )
    { problemInput->setOptimizationSettingsValue( optimizationSettingsEnumList[ i ], static_cast< unsigned int >( optimizationSettings[ i ] ) ); }

    std::cout << "    Convert Objective Function Case optimization setting to character and set" << std::endl;
    const char objectiveFunctionCase = static_cast< char >( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::objective_function_case ) );
    const std::string objectiveFunctionCaseString( 1, objectiveFunctionCase );

    const int optimizerIndex = int( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::optimizer ) );
    const bool optimizerCustomSettingsFlag = bool( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::custom_optimizer_parameters_flag ) );
    std::string algo_method;
    std::string algo_method_printout;
    switch( optimizerIndex )
    {
    case 0:
    { algo_method = "NSGA2"; algo_method_printout = algo_method; break; }
    case 1:
    { algo_method = "MOEA/D (DE)"; algo_method_printout = "MOEADDE"; break; }
    case 2:
    { algo_method = "IHS"; algo_method_printout = algo_method; break; }
    }

    /*std::vector< std::string > termsToRemoveFrom_primaryFileName;
    termsToRemoveFrom_primaryFileName.push_back( "primary" );
    termsToRemoveFrom_primaryFileName.push_back( ".txt" );

    if( debugInfo == 1 ){ std::cout << "Removing terms from primary file name" << std::endl; }
    bislip::utilities::eraseSubStrings( primaryFileName, termsToRemoveFrom_primaryFileName );*/

    std::vector< std::string > termsToRemove;
    termsToRemove.push_back( "primary" );
    termsToRemove.push_back( "THESIS_" );
    termsToRemove.push_back( "Thesis_" );
    termsToRemove.push_back( "thesis_" );
    termsToRemove.push_back( "VALIDATION_" );
    termsToRemove.push_back( "Validation_" );
    termsToRemove.push_back( "test_" );
    termsToRemove.push_back( "Test_" );
    termsToRemove.push_back( "TEST_" );
    termsToRemove.push_back( "THESISnodalParameters_" );
    termsToRemove.push_back( "nodalParameters_" );
    termsToRemove.push_back( "nodalParameters_" );
    termsToRemove.push_back( "NodalParameters_" );
    termsToRemove.push_back( "nodalParameterFiles" );
    termsToRemove.push_back( "HORUS_" );
    termsToRemove.push_back( "/" );
    termsToRemove.push_back( ".txt" );

    std::string this_run_settings =
            caseToRun + "_" +
            nodalParameterInputMatrixFileName + "_" +
            algo_method_printout + "_" +
            std::to_string( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::seed_generator_initializer ) ) + "_" +
            objectiveFunctionCase;
    bislip::utilities::eraseSubStrings( this_run_settings, termsToRemove );

    if( debugInfo == 1 ){ std::cout << "Create output subfolder filename. Based on arbitrary prefix, the" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    previously created partial suffix, and the run's time stamp." << std::endl; }
    std::string outputSubFolder = "OUTPUT_" + this_run_settings + "_" + playTimeString;

    if( resumeOptimization ) { outputSubFolder = resumeOptimizationFolderPath; }
    std::string outputPath = bislip::Variables::getOutputPath( mainSubFolder + "/" );
    problemInput->setOutputPath( outputPath );
    problemInput->setOutputSubFolder( outputSubFolder );

    tudat::input_output::writeDataMapToTextFile( bislip::Variables::createTargetCircle( targetLat_rad, targetLon_rad, finalDistanceToTarget_rad ),
                                                 "targetCircleCoordinates.dat",
                                                 outputPath + outputSubFolder,
                                                 "",
                                                 std::numeric_limits< double >::digits10,
                                                 std::numeric_limits< double >::digits10,
                                                 "," );

    std::cout << "Constructing Decision Vector Bounds Matrix" << std::endl;
    unsigned int n_totalParameters = nonNodalParameterBounds.rows();
    Eigen::MatrixXd parameterMatrix( n_totalParameters, 2 );
    parameterMatrix = nonNodalParameterBounds;

    std::map< bislip::trajectory_phases::trajectory_phase_list, std::map< bislip::parameters::nodal_parameter, std::map< bislip::parameters::node, Eigen::Vector2d > > > nodalParameterBoundsPerTrajectoryPhaseMap;
    bool ascentIsPresent, descentIsPresent, AWD, AWOD, DWOA;

    if( validation == true )
    {
        ascentIsPresent = false;
        descentIsPresent = true;
        AWD = false;
        AWOD = false;
        DWOA = true;

        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_is_present, ascentIsPresent );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_is_present, descentIsPresent );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_with_descent, AWD );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_without_descent, AWOD );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_without_ascent, DWOA );


    }
    else
    {
        std::cout << "Processing nodal parameters" << std::endl;

        //bislip::parameters::processNodalParameters( problemInput, nodalParameterInputMatrix, nonNodalParameterBounds, nodalParameterBoundsPerTrajectoryPhaseMap, parameterMatrix, n_totalParameters );

        /*
        void processNodalParameters(
                    const std::shared_ptr< bislip::problem_input > &problemInput,
                    const Eigen::MatrixXd &nodalParameterInputMatrix )
        {
*/    nodalParameterBoundsPerTrajectoryPhaseMap = bislip::parameters::createMapOfNodalParameterBounds( nodalParameterInputMatrix );

        problemInput->setNodalParameterBoundsPerTrajectoryPhaseMap( nodalParameterBoundsPerTrajectoryPhaseMap );

        unsigned int n_nodesAscent = 0, n_nodesDescent = 0, n_intervalsAscent = 0, n_intervalsDescent = 0, n_nodalVariablesAscent = 0, n_nodalVariablesDescent = 0;
        std::map< bislip::parameters::n_map, unsigned int > nodesPerAscentMap, nodesPerDescentMap;
        std::map< bislip::trajectory_phases::trajectory_phase_list, std::map< bislip::parameters::n_map, unsigned int > > nodesPerTrajectoryMap;


        if( nodalParameterBoundsPerTrajectoryPhaseMap.count( bislip::trajectory_phases::trajectory_phase_list::ascent ) )
        {
            n_nodesAscent = nodalParameterBoundsPerTrajectoryPhaseMap.at( bislip::trajectory_phases::trajectory_phase_list::ascent ).at( bislip::parameters::nodal_parameter::angle_of_attack ).size();
            n_intervalsAscent = n_nodesAscent - 1;
            n_nodalVariablesAscent = nodalParameterBoundsPerTrajectoryPhaseMap.at( bislip::trajectory_phases::trajectory_phase_list::ascent ).size();
            nodesPerAscentMap[ bislip::parameters::n_map::number_of_nodes ] = n_nodesAscent;
            nodesPerAscentMap[ bislip::parameters::n_map::number_of_intervals ] = n_intervalsAscent;
            nodesPerAscentMap[ bislip::parameters::n_map::nodal_variables ] = n_nodalVariablesAscent;
            nodesPerTrajectoryMap[ bislip::trajectory_phases::trajectory_phase_list::ascent ] = nodesPerAscentMap;
        }
        if( nodalParameterBoundsPerTrajectoryPhaseMap.count( bislip::trajectory_phases::trajectory_phase_list::descent ) )
        {
            n_nodesDescent = nodalParameterBoundsPerTrajectoryPhaseMap.at( bislip::trajectory_phases::trajectory_phase_list::descent ).at( bislip::parameters::nodal_parameter::angle_of_attack ).size();
            n_intervalsDescent = n_nodesDescent - 1;
            n_nodalVariablesDescent = nodalParameterBoundsPerTrajectoryPhaseMap.at( bislip::trajectory_phases::trajectory_phase_list::descent ).size();
            nodesPerDescentMap[ bislip::parameters::n_map::number_of_nodes ] = n_nodesDescent;
            nodesPerDescentMap[ bislip::parameters::n_map::number_of_intervals ] = n_intervalsDescent;
            nodesPerDescentMap[ bislip::parameters::n_map::nodal_variables ] = n_nodalVariablesDescent;
            nodesPerTrajectoryMap[ bislip::trajectory_phases::trajectory_phase_list::descent ] = nodesPerDescentMap;
        }

        problemInput->setNodesPerTrajectoryPhaseMap( nodesPerTrajectoryMap );

        const unsigned int n_nodalVariables = bislip::parameters::getNodalParameters().size();
        const unsigned int n_nodalParameters = static_cast< unsigned int >( nodalParameterInputMatrix.rows() ) + n_intervalsAscent + n_intervalsDescent;
        const unsigned int n_staticParameters = static_cast< unsigned int >( nonNodalParameterBounds.rows() );
        n_totalParameters += n_nodalParameters;

        //    problemInput->setNodesPerTrajectoryPhase( bislip::trajectory_phases::trajectory_phase_list::ascent, n_nodesAscent );
        //  problemInput->setNodesPerTrajectoryPhase( bislip::trajectory_phases::trajectory_phase_list::descent, n_nodesDescent );

        std::map< bislip::parameters::n_map, unsigned int > n_map;
        n_map[ bislip::parameters::n_map::nodes_ascent ]      = n_nodesAscent;
        n_map[ bislip::parameters::n_map::nodes_descent ]     = n_nodesDescent;
        n_map[ bislip::parameters::n_map::intervals_ascent ]  = n_intervalsAscent;
        n_map[ bislip::parameters::n_map::intervals_descent ] = n_intervalsDescent;
        n_map[ bislip::parameters::n_map::nodal_variables ]   = n_nodalVariables;
        n_map[ bislip::parameters::n_map::nodal_parameter_input_matrix_rows ] = static_cast< unsigned int >( nodalParameterInputMatrix.rows() );
        n_map[ bislip::parameters::n_map::nodal_parameters ]  = n_nodalParameters;
        n_map[ bislip::parameters::n_map::static_parameters ] = n_staticParameters;
        n_map[ bislip::parameters::n_map::total_parameters ]  = n_totalParameters;

        problemInput->setParametern_map( n_map );

        Eigen::MatrixXd nodalParameterBounds( nodalParameterInputMatrix.rows(), 2 );
        nodalParameterBounds.col( 0 ) = nodalParameterInputMatrix.col( nodalParameterInputMatrix.cols() - 2 );
        nodalParameterBounds.col( 1 ) = nodalParameterInputMatrix.col( nodalParameterInputMatrix.cols() - 1 );

        Eigen::MatrixXd nodeIntervalBounds( n_intervalsAscent + n_intervalsDescent, 2 );
        nodeIntervalBounds = Eigen::MatrixXd::Zero( n_intervalsAscent + n_intervalsDescent, 2 );
        nodeIntervalBounds.col( 1 ) = Eigen::MatrixXd::Ones( n_intervalsAscent + n_intervalsDescent, 1 );

        //! Resize and re-populate parameter matrix
        parameterMatrix.resize( nonNodalParameterBounds.rows() + nodalParameterBounds.rows() + nodeIntervalBounds.rows(), 2 );
        parameterMatrix << nonNodalParameterBounds, nodeIntervalBounds, nodalParameterBounds;
        Eigen::MatrixXd parameterMatrixTransposed( 2, parameterMatrix.rows() );
        parameterMatrixTransposed = parameterMatrix.transpose();

        ascentIsPresent = bool( nodalParameterBoundsPerTrajectoryPhaseMap.count( bislip::trajectory_phases::trajectory_phase_list::ascent ) );
        descentIsPresent = bool( nodalParameterBoundsPerTrajectoryPhaseMap.count( bislip::trajectory_phases::trajectory_phase_list::descent ) );
        AWD = ascentIsPresent && descentIsPresent;
        AWOD = ascentIsPresent && !descentIsPresent;
        DWOA = !ascentIsPresent && descentIsPresent;

        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_is_present, ascentIsPresent );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_is_present, descentIsPresent );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_with_descent, AWD );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_without_descent, AWOD );
        problemInput->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_without_ascent, DWOA );

        //}
        std::cout << "Nodal parameters processed" << std::endl;

    }


    //! Create vector containing decision vector bounds.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( n_totalParameters, 0.0 ) );
    std::vector< double > validationDecisionVector, testCaseDecisionVector;
    //!---------------------------------------   ^ for lb/up  (rows)       ^ for # of parameters

    for( unsigned int i = 0; i < n_totalParameters ; i++ )
    {
        bounds[ 0 ][ i ] = parameterMatrix( i, 0 );
        bounds[ 1 ][ i ] = parameterMatrix( i, 1 );
    }

    if( validation )
    {
        for( unsigned int i = 0; i < n_totalParameters ; i++ )
        {
            validationDecisionVector.push_back( ( bounds[ 0 ][ i ] + bounds[ 1 ][ i ] ) / 2.0 );
        }
    }
    if( testCase )
    {
        for( unsigned int i = 0; i < n_totalParameters ; i++ )
        {
            testCaseDecisionVector.push_back( ( bounds[ 0 ][ i ] + bounds[ 1 ][ i ] ) / 2.0 );
        }
    }
    /*
    for( unsigned int i = 0; i < n_totalParameters ; i++ )
    {
        std::cout << "bounds[ 0 ][ " << i << " ] = " << bounds[ 0 ][ i ] << std::endl;
        std::cout << "bounds[ 1 ][ " << i << " ] = " << bounds[ 1 ][ i ] << std::endl;
    }
    */

    if( debugInfo == 1 ){ std::cout << "Bounds Matrix Constructed" << std::endl; }

    problemInput->setDecisionVectorBounds( bounds );

    //! Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating environment" << std::endl;

    //! Declare and initialize central body name.
    const std::string centralBodyName = "Earth";
    problemInput->setCentralBodyName( centralBodyName );

    bodyNamesMap[ bislip::BislipVehicleSystems::body_names::central_body_name ] = centralBodyName;

    //! Declare and initialize simulation body settings data map.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { centralBodyName },
                                    simulationStartEpoch - 10000.0 * propagationStepSize,
                                    simulationEndEpoch + 10000.0 * propagationStepSize );

    // std::cout << "Define ephemeris model settings." << std::endl;
    //! Define ephemeris model settings.
    //! This is an acceptable 'cheat' were Earth is placed at the barycenter.
    //! Use only when there arent any third body perturbations (Moon, Sun, etc.)
    bodySettings[ centralBodyName ]->ephemerisSettings =
            std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    // std::cout << "Reset ephemeris to J2000." << std::endl;
    //! Reset ephemeris to J2000.
    bodySettings[ centralBodyName ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    //! Define Earth's radius. Using spice here. Is there a way to get a 'radius
    //! field'? Im interested in includig Earth's flattening, yet am unsure how
    //! to properly do it.
    const double radiusEarth = tudat::spice_interface::getAverageRadius( centralBodyName );
    //double radiusEarth_i = radiusEarth;
    //double radiusEartmaximumHeightAllowable = radiusEarth;

    //!---------------------------- Still playing around with this section

    //! https://en.wikipedia.org/wiki/Flattening#Numerical_values_for_planets
    //! https://www.vcalc.com/wiki/vCalc/WGS-84+Earth+flattening+factor
    double flattening_Earth = 1 / 1;//nput_data_.back();
    if ( flattening_Earth != 1 )
    {
        bodySettings[ centralBodyName ]->shapeModelSettings = std::make_shared< OblateSphericalBodyShapeSettings >( radiusEarth, flattening_Earth );

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

        //const double radius_Equator =  6378137.0;
        //const double radius_pole = 6356752.3142;
        /* const double polar_angle_i_deg = 90 - lat_i_deg;
        const double polar_angle_f_deg = 90 - lat_f_deg;
        const double polar_angle_i_rad = polar_angle_i_deg * mathematical_constants::PI / 180;
        const double polar_angle_f_rad = polar_angle_f_deg * mathematical_constants::PI / 180;

        const double radiusEarth_i_x = radii[ 0 ] * std::cos(polar_angle_i_rad) * std::cos(initialLon_rad);
        const double radiusEarth_i_y = radii[ 1 ] * std::cos(polar_angle_i_rad) * std::sin(initialLon_rad);
        const double radiusEarth_i_z = radii[ 2 ] * std::sin(polar_angle_i_rad);
        radiusEarth_i = std::sqrt(pow(radiusEarth_i_x,2) + pow(radiusEarth_i_y,2) + pow(radiusEarth_i_z,2));
        const double radiusEartmaximumHeightAllowable_x = radii[ 0 ] * std::cos(polar_angle_f_rad) * std::cos(targetLon_rad);
        const double radiusEartmaximumHeightAllowable_y = radii[ 1 ] * std::cos(polar_angle_f_rad) * std::sin(targetLon_rad);
        const double radiusEartmaximumHeightAllowable_z = radii[ 2 ] * std::sin(polar_angle_f_rad);
        radiusEartmaximumHeightAllowable = std::sqrt(pow(radiusEartmaximumHeightAllowable_x,2) + pow(radiusEartmaximumHeightAllowable_y,2) + pow(radiusEartmaximumHeightAllowable_z,2));
*/
    }
    //!--------------------------------------------------------------

    //const double initialAltitude = initialHeight + radiusEarth;

    double initialSpeedOfSound = 0.0;
    double initialDensity = 0.0;
    if( validation )
    {
        guidanceStepSize = propagationStepSize;

        // Define inputs
        std::string atmosphereFile = tudat::input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";

        std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables;
        dependentVariables.push_back( aerodynamics::density_dependent_atmosphere );
        dependentVariables.push_back( aerodynamics::pressure_dependent_atmosphere );
        dependentVariables.push_back( aerodynamics::temperature_dependent_atmosphere );

        double specificGasConstant = 287.0;

        double ratioOfSpecificHeats = 1.4;

        tudat::interpolators::BoundaryInterpolationType boundaryHandling = tudat::interpolators::use_default_value_with_warning;

        double defaultExtrapolationValue = TUDAT_NAN;

        // Set tabulated atmosphere
        bodySettings[ centralBodyName ]->atmosphereSettings = std::make_shared< tudat::simulation_setup::TabulatedAtmosphereSettings >( atmosphereFile, dependentVariables, specificGasConstant, ratioOfSpecificHeats, boundaryHandling, defaultExtrapolationValue );

        if( debugInfo == 1 ){ std::cout << "    Creating Atmospheric Model" << std::endl; }
        //! US Standard Atmosphere Model 1976
        //std::string atmosphereTableFile = tudat::input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
        tudat::aerodynamics::TabulatedAtmosphere atmosphericModel_US76( atmosphereFile );

        if( debugInfo == 1 ){std::cout << "    Evaluating US76 Atmosphere Model" << std::endl; }
        initialSpeedOfSound = atmosphericModel_US76.getSpeedOfSound( initialHeight );
        initialDensity = atmosphericModel_US76.getDensity( initialHeight );

        if( debugInfo == 1 ){std::cout << "         Initial Speed of Sound = " << initialSpeedOfSound << std::endl; }
        if( debugInfo == 1 ){std::cout << "         Initial Density        = " << initialDensity << std::endl; }

    }
    else
    {
        bodySettings[ centralBodyName ]->atmosphereSettings = std::make_shared< tudat::simulation_setup::AtmosphereSettings >( nrlmsise00 );

        if( debugInfo == 1 ){ std::cout << "    Creating Atmospheric Model" << std::endl; }
        //! NRLMSISE-00 Atmosphere Model
        std::string spaceWeatherFilePath = "spaceweather.txt";

        tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
                tudat::input_output::solar_activity::readSolarActivityData( spaceWeatherFilePath ) ;

        std::function< tudat::aerodynamics::NRLMSISE00Input ( double, double, double, double ) > atmosphericFunction_NRLMSISE00 =
                std::bind( &tudat::aerodynamics::nrlmsiseInputFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, solarActivityData , true , 0 );

        tudat::aerodynamics::NRLMSISE00Atmosphere atmosphericModel_NRLMSISE00( atmosphericFunction_NRLMSISE00 );

        if( debugInfo == 1 ){std::cout << "    Evaluating NRLMSISE-00 Atmosphere Model" << std::endl; }
        initialSpeedOfSound = atmosphericModel_NRLMSISE00.getSpeedOfSound( initialHeight, initialLon_rad, initialLat_rad, simulationStartEpoch );
        initialDensity = atmosphericModel_NRLMSISE00.getDensity( initialHeight, initialLon_rad, initialLat_rad, simulationStartEpoch );

        if( debugInfo == 1 ){std::cout << "         Initial Speed of Sound = " << initialSpeedOfSound << std::endl; }
        if( debugInfo == 1 ){std::cout << "         Initial Density        = " << initialDensity << std::endl; }
    }

    if( debugInfo == 1 ){ std::cout << "Creating map of desity parameters for Jerk computation" << std::endl; }

    std::map< int, Eigen::VectorXd > densityParameterMapForJerk;
    densityParameterMapForJerk[ 0 ] = dummy;
    densityParameterMapForJerk.erase( 0 );

    for( int i = 0; i < densityParameterMatrixForJerk.rows(); i++ )
    { densityParameterMapForJerk[ i ] = ( densityParameterMatrixForJerk.row( i ) ).tail( densityParameterMatrixForJerk.cols() - 1 ); }

    if( debugInfo == 1 ){ std::cout << "Creating Earth object" << std::endl; }

    //! Create Earth object
    tudat::simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );
    tudat::simulation_setup::NamedBodyMap bodyMapCopy;

    //std::cout << "Create Earth's rotational ephemeris" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating vehicle" << std::endl;

    ///////// Start: Vehicle Aerodynamics Section
    if( debugInfo == 1 ){ std::cout << "           Start: Vehicle Aerodynamics Section" << std::endl; }

    //! Define physical meaning of independent variables
    //!     Angle of attack
    //!     Mach number
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );

    // Define physical meaning of independent variables for control surface increments, in this case Mach number, angle of attack and control surface deflection
    //!     Angle of attack
    //!     Mach number
    //!     Control surface deflection
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > controlSurfaceIndependentVariableNames;
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::mach_number_dependent );
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::control_surface_deflection_dependent );

    //! Define physical meaning of independent variables for drag increment due to altitude
    //!     Angle of attack
    //!     Mach number
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > altitudeIndependentVariableNames;
    altitudeIndependentVariableNames.push_back( aerodynamics::altitude_dependent );
    altitudeIndependentVariableNames.push_back( aerodynamics::mach_number_dependent );

    //! Define list of files for force coefficients. Entry denotes direction.
    //!     0 : x-direction (C ~D~/C ~X~)
    //!     1 : y-direction (C ~S~/C ~Y~)
    //!     2 : z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles;
    std::map< int, std::string > forceCoefficientFiles_CS_B;
    std::map< int, std::string > forceCoefficientFiles_CS_EL;
    std::map< int, std::string > forceCoefficientFiles_CS_ER;
    std::map< int, std::string > forceCoefficientFiles_Altitude;
    std::map< int, std::string > momentCoefficientFiles;
    std::map< int, std::string > momentCoefficientFiles_CS_B;
    std::map< int, std::string > momentCoefficientFiles_CS_EL;
    std::map< int, std::string > momentCoefficientFiles_CS_ER;

    std::string BODYFLAP = "BodyFlap";
    std::string ELEVON_L = "ElevonLeft";
    std::string ELEVON_R = "ElevonRight";
    std::string ALTITUDE = "Altitude";

    std::map< bislip::BislipVehicleSystems::control_surface_list, std::string > controlSurfaceNameMap;
    controlSurfaceNameMap[ bislip::BislipVehicleSystems::control_surface_list::bodyflap ] = BODYFLAP;
    controlSurfaceNameMap[ bislip::BislipVehicleSystems::control_surface_list::elevon_left ] = ELEVON_L;
    controlSurfaceNameMap[ bislip::BislipVehicleSystems::control_surface_list::elevon_right ] = ELEVON_R;

    if( debugInfo == 1 ){ std::cout << "           Assigning files names for Aerodynamic Coefficient Database" << std::endl; }

    // Define list of files for force coefficients.
    forceCoefficientFiles[ 0 ] = aeroCoeffFileList[ 0 ]; // Set clean coefficient drag coefficient file
    forceCoefficientFiles[ 2 ] = aeroCoeffFileList[ 1 ]; // Set clean coefficient lift coefficient file

    // Define list of files for moment coefficients.
    momentCoefficientFiles[ 1 ] = aeroCoeffFileList[ 2 ]; // Set clean coefficient pitch moment coefficient file

    // Define list of files for force coefficients for control surfaces ( Bodyflap )
    forceCoefficientFiles_CS_B[ 0 ] = aeroCoeffFileList[ 3 ];
    forceCoefficientFiles_CS_B[ 2 ] = aeroCoeffFileList[ 4 ];

    // Define list of files for moment coefficients for control surfaces ( Bodyflap )
    momentCoefficientFiles_CS_B[ 1 ] = aeroCoeffFileList[ 5 ];

    // Define list of files for force coefficients for control surfaces ( Elevon - Left )
    forceCoefficientFiles_CS_EL[ 0 ] = aeroCoeffFileList[ 6 ];
    forceCoefficientFiles_CS_EL[ 2 ] = aeroCoeffFileList[ 7 ];

    // Define list of files for moment coefficients for control surfaces ( Elevon - Left )
    momentCoefficientFiles_CS_EL[ 1 ] = aeroCoeffFileList[ 8 ];

    // Define list of files for force coefficients for control surfaces ( Elevon - Right )
    forceCoefficientFiles_CS_ER[ 0 ] = aeroCoeffFileList[ 9 ];
    forceCoefficientFiles_CS_ER[ 2 ] = aeroCoeffFileList[ 10 ];

    // Define list of files for moment coefficients for control surfaces ( Elevon - Right )
    momentCoefficientFiles_CS_ER[ 1 ] = aeroCoeffFileList[ 11 ];

    // Define list of files for force coefficients for control surfaces ( Elevon - Right )
    forceCoefficientFiles_Altitude[ 0 ] = aeroCoeffFileList[ 12 ];

    //! Define reference frame in which the loaded coefficients are defined.
    //! Have to get some more background info here to properly understand it.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;


    std::vector< tudat::interpolators::BoundaryInterpolationType > aerodynamicCoefficientCleanConfigurationBoundaryHandling( 2 );

    aerodynamicCoefficientCleanConfigurationBoundaryHandling[ 0 ] = tudat::interpolators::BoundaryInterpolationType::extrapolate_at_boundary ;
    aerodynamicCoefficientCleanConfigurationBoundaryHandling[ 1 ] = tudat::interpolators::BoundaryInterpolationType::extrapolate_at_boundary ;
    // aerodynamicCoefficientCleanConfigurationBoundaryHandling[ 2 ] = tudat::interpolators::BoundaryInterpolationType::extrapolate_at_boundary ;


    std::shared_ptr< tudat::interpolators::InterpolatorSettings > multi3DLinearInterpolatorSettings =
            std::make_shared< tudat::interpolators::InterpolatorSettings >(
                tudat::interpolators::InterpolatorTypes::multi_linear_interpolator,
                tudat::interpolators::AvailableLookupScheme::huntingAlgorithm,
                false,
                aerodynamicCoefficientCleanConfigurationBoundaryHandling );

    if( debugInfo == 1 ){ std::cout << "           Load and parse coefficient files; create coefficient settings" << std::endl; }

    //! Load and parse coefficient files; create coefficient settings.
    std::shared_ptr< tudat::simulation_setup::AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            tudat::simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles,
                momentCoefficientFiles,
                b_ref,
                S_ref,
                c_ref,
                R_mrc,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame,
                areCoefficientsInNegativeAxisDirection,
                multi3DLinearInterpolatorSettings );

    if( debugInfo == 1 ){ std::cout << "           Add settings for control surface increments to main aerodynamic coefficients" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Bodyflap" << std::endl; }

    //! Add settings for control surface increments to main aerodynamic coefficients
    aerodynamicCoefficientSettings->setControlSurfaceSettings(
                tudat::simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles_CS_B,
                    momentCoefficientFiles_CS_B,
                    controlSurfaceIndependentVariableNames), BODYFLAP );

    if( debugInfo == 1 ){ std::cout << "                Left Elevon" << std::endl; }

    aerodynamicCoefficientSettings->setControlSurfaceSettings(
                tudat::simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles_CS_EL,
                    momentCoefficientFiles_CS_EL,
                    controlSurfaceIndependentVariableNames), ELEVON_L );

    if( debugInfo == 1 ){ std::cout << "                Right Elevon" << std::endl; }

    aerodynamicCoefficientSettings->setControlSurfaceSettings(
                tudat::simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles_CS_ER,
                    momentCoefficientFiles_CS_ER,
                    controlSurfaceIndependentVariableNames), ELEVON_R );

    // if( debugInfo == 1 ){ std::cout << "           --------> Altitude" << std::endl; }

    // aerodynamicCoefficientSettings->setControlSurfaceSettings(
    //             tudat::simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
    //                 forceCoefficientFiles_Altitude,
    //                 altitudeIndependentVariableNames), ALTITUDE );

    ///////// End: Vehicle Aerodynamics Section
    if( debugInfo == 1 ){ std::cout << "           End: Vehicle Aerodynamics Section" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Create vehicle object" << std::endl; }

    //! Create vehicle objects.
    bodyMap[ vehicleName ] = std::make_shared< tudat::simulation_setup::Body >( );

    if( debugInfo == 1 ){ std::cout << "    Create shared pointer for Vehicle Systems" << std::endl; }

    //! Create vehicle systems and initialize by setting the landing (dry) Mass.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = std::make_shared< tudat::system_models::VehicleSystems >( landingMass );

    if( debugInfo == 1 ){ std::cout << "    Create shared pointer for Bislip Systems" << std::endl; }
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = std::make_shared< bislip::BislipVehicleSystems >( landingMass );
    //! Pass shared pointer of vehicle and bislip systems.
    bodyMap.at( vehicleName )->setVehicleSystems( vehicleSystems );
    bodyMap.at( vehicleName )->setBislipSystems( bislipSystems );

    vehicleSystems->setNoseRadius( noseRadius );
    vehicleSystems->setWallEmissivity( wallEmissivity );


    bislipSystems->setOutputPath( outputPath );
    bislipSystems->setOutputSubFolder( outputSubFolder );
    bislipSystems->setBodyNames( bodyNamesMap );
    bislipSystems->setDebugInfo( debugInfo );

    bislipSystems->setStartingEpoch( simulationStartEpoch );
    bislipSystems->setWingSweepAngle( lambda );
    bislipSystems->setLocalBodyAngle( phi );
    bislipSystems->setTransitionDistance( x_T );
    bislipSystems->setCentralTargetLatitude( targetLat_rad );
    bislipSystems->setCentralTargetLongitude( targetLon_rad );
    bislipSystems->setAngularDistanceForTermination( finalDistanceToTarget_rad );
    bislipSystems->setCentralTargetCoordinates( centralTargetCoordinates );
    bislipSystems->setFinalDistanceToTarget( finalDistanceToTarget_rad );
    bislipSystems->setReferenceArea( S_ref );
    bislipSystems->setMinimumDynamicPressureforControlSurface( minimumDynamicPressureforControlSurface );
    bislipSystems->setAverageEarthRadius( radiusEarth );

    bislipSystems->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_is_present, ascentIsPresent );
    bislipSystems->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_is_present, descentIsPresent );
    bislipSystems->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_with_descent, AWD );
    bislipSystems->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_without_descent, AWOD );
    bislipSystems->setTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_without_ascent, DWOA );

    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::latitude_angle, initialLat_rad );
    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::longitude_angle, initialLon_rad );
    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::distance_to_central_target, initialDistanceToTarget_rad );
    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::flight_path_angle, initialFlightPathAngle_rad );
    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::height, initialHeight );
    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::speed_of_sound, initialSpeedOfSound );
    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::density, initialDensity );
    bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::mach_number, initialMachNumber );

    bislipSystems->setInitialCoordinates( std::make_pair( bislipSystems->getInitialValue( bislip::initial_values::initial_value_type::latitude_angle ), bislipSystems->getInitialValue( bislip::initial_values::initial_value_type::longitude_angle ) ) );
    bislipSystems->setDensityParameterMapForJerk( densityParameterMapForJerk );

    bislipSystems->setBodyFlapDeflectionLimits( std::make_pair( minimumBodyflapAngle, maximumBodyflapAngle ) );
    bislipSystems->setElevonDeflectionLimits( std::make_pair( minimumElevonAngle, maximumElevonAngle ) );
    bislipSystems->setThrustElevationLimits( std::make_pair( minimumThrustElevationAngle, maximumThrustElevationAngle ) );
    bislipSystems->setThrottleSettingLimits( std::make_pair( minimumThrottleSetting, maximumThrottleSetting ) );
    bislipSystems->setBankAngleReversalTimepoint( simulationStartEpoch );
    bislipSystems->setControlSurfaceNameMap( controlSurfaceNameMap );

    std::vector< bislip::vehicle_parameters::vehicle_parameter_list > vehicleParameters = bislip::vehicle_parameters::getVehicleParameterList();
    for( unsigned int i = 0; i < vehicleParameterList.size(); i++ )
    { bislipSystems->setVehicleParameter( vehicleParameters[ i ], vehicleParameterValues[ i ] );}

    //! Initialize all flags
    std::cout << "Initialize flags" << std::endl;
    bislipSystems->setFlag( bislip::flags::flag_list::validation, validation );
    bislipSystems->setFlag( bislip::flags::flag_list::test_case, testCase );
    bislipSystems->setFlag( bislip::flags::flag_list::evolution_evaluation, true );
    bislipSystems->setFlag( bislip::flags::flag_list::flight_corridor_evaluation, false );
    bislipSystems->setFlag( bislip::flags::flag_list::low_distance, false );
    bislipSystems->setFlag( bislip::flags::flag_list::bank_angle_reversal, false );
    bislipSystems->setFlag( bislip::flags::flag_list::thrust_maximizer, maximizeThrust );
    bislipSystems->setFlag( bislip::flags::flag_list::impose_parameter_bounds_on_interpolators, true );

    if( bislipSystems->getFlag( bislip::flags::flag_list::validation ) )
    {
        std::cout << "  Validation case --> Set thrust maximizer flag to false" << std::endl;
        bislipSystems->setFlag( bislip::flags::flag_list::thrust_maximizer, false );
        //bislipSystems->setFlag( bislip::flags::flag_list::evolution_evaluation, false );
    }
    else
    {
        std::cout << "  Set thrust maximizer flag" << std::endl;
        if( bislipSystems->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_is_present ) )
        {
            std::cout << "      Ascent is present" << std::endl;
            if( nodalParameterBoundsPerTrajectoryPhaseMap.at( bislip::trajectory_phases::trajectory_phase_list::ascent ).count( bislip::parameters::nodal_parameter::throttle_setting ) )
            {
                std::cout << "          Throttle setting is a nodal parameter --> Set thrust maximizer flag to false" << std::endl;
                bislipSystems->setFlag( bislip::flags::flag_list::thrust_maximizer, false );
            }

        }
        else if( bislipSystems->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_without_ascent ) )
        {
            std::cout << "      Desscent is present" << std::endl;
            if( nodalParameterBoundsPerTrajectoryPhaseMap.at( bislip::trajectory_phases::trajectory_phase_list::descent ).count( bislip::parameters::nodal_parameter::throttle_setting ) )
            {
                std::cout << "          Throttle setting is a nodal parameter --> Set thrust maximizer flag to false" << std::endl;
                bislipSystems->setFlag( bislip::flags::flag_list::thrust_maximizer, false );
            }
        }
    }

    if( bislipSystems->getFlag( bislip::flags::flag_list::validation ) )
    {
        double initialAirspeed_Ascent = initialMachNumber;
        bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::airspeed, initialAirspeed_Ascent );
        bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::mach_number, initialAirspeed_Ascent / bislipSystems->getInitialValue( bislip::initial_values::initial_value_type::speed_of_sound ) );
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "  Set initial Mach number" << std::endl; }
        bislipSystems->setInitialValue( bislip::initial_values::initial_value_type::mach_number, initialMachNumber );
    }

    if( debugInfo == 1 ){ std::cout << "Set constraint map" << std::endl; }
    bislipSystems->setConstraintMap( problemInput->getConstraintMap() );

    if( debugInfo == 1 ){ std::cout << "Create and set aerodynamic coefficient interface" << std::endl; }
    //! Set vehicle aerodynamic coefficients.
    bodyMap.at( vehicleName )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, vehicleName ) );

    //std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
    //controlSurfaceList[ ELEVON_L ] = controlSurfaceInterface;

    //! Set vehicle system models for aerodynamic coefficient increments.
    //bodyMap[ vehicleName ]->setVehicleSystems( systemsModels );
    //bodyMap[ vehicleName ]->getAerodynamicCoefficientInterface( )->setControlSurfaceIncrements( controlSurfaceList );

    //! Finalize body creation. Not entirely sure what this does. I believe it
    //! may place the body (Earth) at the barycenter. Probably another hack to
    //! facilitate bringing in everything into the inertial frame.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    bodyMap.at( vehicleName )->setEphemeris(
                std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ), centralBodyName ) );
    /*
     *The problem with the mass propagation can be solved (-ish), by turning the throw.... of this error into a std::cerr (printing the error, but not terminating the program). Alternatively, it should also work if you use:
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
    ///////////////////////            CREATE KNOWN STATES              ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating known states" << std::endl;
    std::cout << "Creating (partial) known states" << std::endl;

    //! Impose constraints on first and last energy nodesAscent
    //const double mu = spice_interface::getBodyGravitationalParameter( centralBodyName );
    //a = 301.7;//NRLMSISE00Atmosphere::getSpeedOfSound( R_E + height( 0 ), 0, 0, simulationStartEpoch );
    //double V_i = a * initialMachNumber;
    //double V_DN = 0.99 * sqrt( mu / ( R_E + maximumHeightAllowable ) );

    //! Set spherical elements for vehicle's initial state. Arbitrarily chosen to be in Earth-Fixed frame.
    //!     Transformation to Inertial Frame is done within fitness function, as the initial velocity is a
    //!     decision variable.
    Eigen::Vector6d initialState_spherical;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::radiusIndex )       = radiusEarth + initialHeight;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::latitudeIndex )     = initialLat_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::longitudeIndex )    = initialLon_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::speedIndex )        = 0.0;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::flightPathIndex )   = 0.0;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.0;

    problemInput->setInitialState_Spherical( initialState_spherical );

    //! Create Earth's rotational ephemeris.
    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( centralBodyName )->getRotationalEphemeris( );

    problemInput->setEarthRotationalEphemeris( earthRotationalEphemeris );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   std::cout << "Setting Accelerations" << std::endl;

    //! Declare acceleration data map.
    tudat::simulation_setup::SelectedAccelerationMap accelerationSettingsMap;

    //! Define gravitational model.
    //!     The central body acts this force on the vehicle.
    if( bislipSystems->getFlag( bislip::flags::flag_list::validation ) )
    {
        //accelerationSettingsMap[ vehicleName ][ centralBodyName ].push_back( std::make_shared< tudat::simulation_setup::AccelerationSettings >( tudat::basic_astrodynamics::AvailableAcceleration::point_mass_gravity ) );
        accelerationSettingsMap[ vehicleName ][ centralBodyName ].push_back( std::make_shared< tudat::simulation_setup::SphericalHarmonicAccelerationSettings >( sphericalHarmonicDegree, sphericalHarmonicOrder ) );
    }
    else
    {
        //!     Arbitrary maximum degree/order. Equivalent functionality to Cartesian with corresponding maximum degree/order.
        accelerationSettingsMap[ vehicleName ][ centralBodyName ].push_back( std::make_shared< tudat::simulation_setup::SphericalHarmonicAccelerationSettings >( sphericalHarmonicDegree, sphericalHarmonicOrder ) );
    }

    //! Define aerodynamic accelerations.
    //!     The atmosphere of the central body acts this force on the vehicle.
    accelerationSettingsMap[ vehicleName ][ centralBodyName ].push_back( std::make_shared< tudat::simulation_setup::AccelerationSettings >( tudat::basic_astrodynamics::AvailableAcceleration::aerodynamic ) );

    //! Declare and initialize thrust guidance pointer.
    std::shared_ptr< bislip::MyGuidance > ThrustGuidance = std::make_shared< bislip::MyGuidance >(
                bodyMap,
                vehicleName,
                centralBodyName,
                bislipSystems->getStartingEpoch() );

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
    //bislipSystems->setMaxThrust( maxThrust );
    //bislipSystems->setSpecificImpulse( specificImpulse );

    //! Define thrust acceleration settings.
    //!     The vehicle acts this force on itself.
    accelerationSettingsMap[ vehicleName ][ vehicleName ].push_back( std::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );

    //std::cout << "Accelerations Set" << std::endl;
    if( debugInfo == 1 ){ std::cout << "Accelerations Set" << std::endl; }

    //! Declare vectors that will contain the names of relevant bodies.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;

    //! Define bodies that will be propagated. Only 1.
    bodiesToIntegrate.push_back( vehicleName );
    problemInput->setBodiesToIntegrate( bodiesToIntegrate );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( centralBodyName );
    problemInput->setCentralBodies( centralBodies );

    //! Set acceleration models.
    tudat::basic_astrodynamics::AccelerationMap accelerationMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationSettingsMap,
                bodiesToIntegrate,
                centralBodies );

    problemInput->setAccelerationMap( accelerationMap );
    bislipSystems->setSphericalHarmonicsParameters( sphericalHarmonicDegree, sphericalHarmonicOrder );


    //! **************************************************************************************
    //! ********** Aerodynamic guidance is set AFTER the accelerations and BEFORE propagating.
    //! **************************************************************************************

    if( debugInfo == 1 ){ std::cout << "Declare and assign aerodynamic guidance functions" << std::endl; }

    //! Declare and assign aerodynamic guidance functions.
    std::shared_ptr< bislip::MyGuidance > AeroGuidance = std::make_shared< bislip::MyGuidance >(
                bodyMap,
                vehicleName,
                centralBodyName,
                bislipSystems->getStartingEpoch() );

    //! Set Guidance angle functions.
    setGuidanceAnglesFunctions( AeroGuidance, bodyMap.at( vehicleName ) );

    //! **************************************************************************************
    //! ********** Create alpha-Mach interpolators delineating aerodyanmic data range.
    //! **************************************************************************************

    if( debugInfo == 1 ){ std::cout << "     Create alpha-Mach interpolators delineating aerodynamic data range" << std::endl; }

    bislip::Variables::createAlphaMachBoundingInterpolators( bodyMap, vehicleName, alphaMachEnvelopeUB, alphaMachEnvelopeLB, outputPath, outputSubFolder );

    //! **************************************************************************************
    //! ********** Create heading error deadband interpolators used for bank angle reversals.
    //! **************************************************************************************

    if( debugInfo == 1 ){ std::cout << "     Create heading error deadband interpolators used for bank angle reversals" << std::endl; }

    bislip::Variables::createHeadingErrorDeadBandInterpolator( bodyMap, vehicleName, headingErrorDeadBandCoarse, headingErrorDeadBandLowDistance, outputPath, outputSubFolder);

    //! **************************************************************************************
    //! ********** Create aerodynamic guidance interpolators used for HORUS Validation.
    //! **************************************************************************************

    if( bislipSystems->getFlag( bislip::flags::flag_list::validation ) )
    {
        //bislipSystems->setMaxThrust( 0.0 );
        bislipSystems->setVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::maximum_thrust, 0.0 );

        if( debugInfo == 1 ){ std::cout << "     Create HORUS Re-entry to Kourou guidance interpolators" << std::endl; }
        bislip::Variables::createKourouGuidanceInterpolators( bodyMap, vehicleName, kourouTimeHistories, outputPath, outputSubFolder );
    }

    if( debugInfo == 1 ){ std::cout << "Guidance is set" << std::endl; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          CREATE LIST OF DEPENDENT VARIABLES        ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating list of dependent variables" << std::endl;

    bislip::createTUDATSettings::createDependentVariableList( problemInput, bodyMap );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE TERMINATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating Termination Settings" << std::endl;

    std::map< bislip::problem_input::termination_constraints, double > terminationConstraintsMap;
    terminationConstraintsMap[ bislip::problem_input::termination_constraints::minimum_mach_number ] = problemInput->getConstraintValue( bislip::constraints::constraint_list::minimum_mach_number );
    terminationConstraintsMap[ bislip::problem_input::termination_constraints::ascent_time_of_flight ] = 0.0;
    terminationConstraintsMap[ bislip::problem_input::termination_constraints::maximum_time_of_flight ] = max_tof;
    terminationConstraintsMap[ bislip::problem_input::termination_constraints::minimum_allowable_height_ascent ] = initialHeight * 9.0 / 10.0;
    terminationConstraintsMap[ bislip::problem_input::termination_constraints::minimum_allowable_height_descent ] = 10000;
    terminationConstraintsMap[ bislip::problem_input::termination_constraints::angular_distance_to_go ] = problemInput->getConstraintValue( bislip::constraints::constraint_list::termination_angular_distance_to_go );
    terminationConstraintsMap[ bislip::problem_input::termination_constraints::angular_distance_traveled ] = tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getInitialValue( bislip::initial_values::initial_value_type::distance_to_central_target ) );

    if( bislipSystems->getFlag( bislip::flags::flag_list::validation ) )
    { terminationConstraintsMap[ bislip::problem_input::termination_constraints::minimum_allowable_height_descent ] = 10000; }

    problemInput->setTerminationConstraintsMap( terminationConstraintsMap );

    bislip::createTUDATSettings::createPropagationTerminationSettings( problemInput, bodyMap );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE MASS RATE SETTINGS: ASCENT            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Calculate initial heading angle to target: https://www.movable-type.co.uk/scripts/latlong.html
    double chi_i_deg_calc = bislip::math_tools::convertNegativeAnglesInDegreesToPositive( unit_conversions::convertRadiansToDegrees( bislip::Variables::computeHeadingToTarget( initialLat_rad , initialLon_rad , targetLat_rad , targetLon_rad ) ) );

    //! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
    const double a = std::sin( (targetLat_rad - initialLat_rad) / 2) * std::sin( (targetLat_rad - initialLat_rad) / 2) + std::cos( initialLat_rad ) * std::cos( initialLon_rad ) * std::sin( (targetLon_rad - initialLon_rad) / 2) * std::sin( (targetLon_rad - initialLon_rad) / 2);
    const double d_angular = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
    const double d_spherical_law_cosines = unit_conversions::convertRadiansToDegrees( bislip::Variables::computeAngularDistance( initialLat_rad , initialLon_rad , targetLat_rad , targetLon_rad ) );
    //std::acos( std::sin(initialLat_rad) * std::sin(targetLat_rad) + std::cos(initialLat_rad) * std::cos(targetLat_rad) * std::cos(targetLon_rad-initialLon_rad) ) * spice_interface::getAverageRadius( "Earth" );

    //! Print to Terminal screen.
    std::cout << " " << std::endl;
    std::cout << problemName << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Non-nodal Parameters" << std::endl;
    for( unsigned int i = 0; i < static_cast< unsigned int >( nonNodalParameterBoundsList.size() ) ; i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.1f %-2s %-8.1f %-2s\n") % nonNodalParameterBoundsList[ i ] % "[" % nonNodalParameterBounds( i, 0 ) % "," % nonNodalParameterBounds( i, 1 ) % "]";
    }

    if( bislipSystems->getFlag( bislip::flags::flag_list::validation ) == false )
    {
        std::cout << "--------------------------------------------------------" << std::endl;
        std::cout << "Nodal Parameters" << std::endl;

        std::cout << nodalParameterInputMatrix << std::endl;

    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Vehicle Parameters" << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[0] << "      " <<  vehicleParameterValues[0] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[1] << "      " <<  vehicleParameterValues[1] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[2] << "      " <<  vehicleParameterValues[2] << std::endl;
    std::cout << std::setw(30) << "R_mrc" << "     [ " << vehicleParameterValues[3] << " , " << vehicleParameterValues[4] << " , " << vehicleParameterValues[5] << " ]" <<  std::endl;
    std::cout << std::setw(30) << "R_com" << "     [ " << vehicleParameterValues[6] << " , " << vehicleParameterValues[7] << " , " << vehicleParameterValues[8] << " ]" <<  std::endl;
    std::cout << std::setw(30) << "R_cot" << "     [ " << vehicleParameterValues[9] << " , " << vehicleParameterValues[10] << " , " << vehicleParameterValues[11] << " ]" <<  std::endl;
    std::cout << std::setw(30) << vehicleParameterList[12] << "      " <<  vehicleParameterValues[12] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[13] << "      " <<  vehicleParameterValues[13] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[14] << "      " <<  vehicleParameterValues[14] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[15] << "      " <<  vehicleParameterValues[15] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[16] << "      " <<  vehicleParameterValues[16] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[17] << "      " <<  vehicleParameterValues[17] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[18] << "      " <<  vehicleParameterValues[18] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[19] << "      " <<  vehicleParameterValues[19] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[20] << "      " <<  vehicleParameterValues[20] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[21] << "      " <<  vehicleParameterValues[21] << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Aerodynamic Coefficients File Names" << std::endl;
    for( unsigned int i = 0; i < static_cast< unsigned int >( aeroCoeffFileList.size() ) ; i++)
    {
        std::cout << "File " << i + 1 << ":      " <<  aeroCoeffFileList[ i ] << std::endl;
    }
    //  std::cout << "Coefficient File: " << "      " <<  aeroCoeffFileList[1] << std::endl;
    //  std::cout << "Coefficient File: " << "      " <<  aeroCoeffFileList[2] << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Optimization Settings" << std::endl;

    std::cout << boost::format("%-25s %-2s %-8.0f\n") % "Optimizer Selected" % "" % algo_method;

    for( unsigned int i = 0; i < static_cast< unsigned int >( optimizationSettings.size() ); i++)
    {
        if ( i == 7 )
            continue;
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % optimizationSettingsList[ i ] % "" % optimizationSettings[ i ];
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Simulation Settings" << std::endl;
    for( unsigned int i = 0; i < static_cast< unsigned int >( simulationSettingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % simulationSettingsList[i] % "" % simulationSettingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output Settings" << std::endl;
    for( unsigned int i = 0; i < static_cast< unsigned int >( outputSettingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % outputSettingsList[i] % "" % outputSettingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Initial Conditions" << std::endl;
    for( unsigned int i = 0; i < static_cast< unsigned int >( initialConditionsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % initialConditionsList[i] % "" % initialConditionsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Constraints / Termination Conditions" << std::endl;
    for( unsigned int i = 0; i < static_cast< unsigned int >( constraintsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % constraintsList[i] % "" % constraintsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Ground distance to cover " << std::endl;
    std::cout << "  Haversine Formula:               " << unit_conversions::convertRadiansToDegrees( d_angular ) << " degrees." << std::endl;
    std::cout << "  Spherical Law of Cosines:        " << d_spherical_law_cosines << " degrees." << std::endl;
    std::cout << "                                   " << d_spherical_law_cosines * radiusEarth / 1E3 << " km." << std::endl;
    std::cout << "  Initial heading angle to Target: " << chi_i_deg_calc << " degrees. Calculated."  << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output Location: " << std::endl;
    std::cout << "     '" << outputPath + outputSubFolder << "'" << std::endl;
    std::cout << " DV   " << "  Lower Boundary " << "  Upper Boundary " << std::endl;
    for( unsigned int i = 0; i < n_totalParameters ; i++ )
    {
        std::cout << std::fixed << std::setprecision( 10 ) <<
                     std::setw( 6 ) << i <<
                     std::setw( 17 ) << bounds[ 0 ][ i ] <<
                     std::setw( 17 ) << bounds[ 1 ][ i ] << std::endl;
    }


    //*********************************************************************************************************************


    bislipSystems->setPropagationStepSize( propagationStepSize );
    bislipSystems->setGuidanceStepSize( guidanceStepSize );

    problemInput->setObjectiveFunctionCase( objectiveFunctionCase );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PAGMO ISLAND(S) & EVOLVE IF APPLICABLE        ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::cout << "Starting Optimization Section" << std::endl;

    double fitnessMagnitudeAverage = TUDAT_NAN;
    double fitnessMagnitudeTopIndividualsVectorAverage = TUDAT_NAN;
    double fitnessMagnitudeTopIndividualsMagnitudeAverage = TUDAT_NAN;
    std::map< double, Eigen::VectorXd > fitnessMagnitudeAveragesOutputHistory;
    std::vector< double > fitnessMagnitudeAverageHistory;
    std::vector< double > fitnessMagnitudeAverageTopIndividualsVectorHistory, fitnessMagnitudeAverageTopIndividualsMagnitudeHistory;
    std::vector< double > fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory, fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory;

    Eigen::VectorXd temp( 9 );


    if( bislipSystems->getFlag( bislip::flags::flag_list::validation ) )
    {
        std::cout << "Running a Single Validation case." << std::endl;
        std::cout << "  Number of Evolutions    = " << 0 << std::endl;
        std::cout << "  Propagation Step Size   = " << propagationStepSize << std::endl;
        std::cout << "  Guidance Step Size      = " << guidanceStepSize << std::endl;
        std::cout << "  " << std::endl;

        problemInput->setIndividualNumber( 0 );
        problemInput->setGenerationNumber( 0 );
        bislipSystems->setFlag( bislip::flags::flag_list::evolution_evaluation, true );
        bislip::Space4ErrBody validationProblem( problemInput, bodyMap );

        std::vector< double > validationRun = validationProblem.fitness( validationDecisionVector );

        bislip::Variables::reruns( bodyMap, problemInput, bislip::Variables::selectTopIndividuals( bodyMap, problemInput, 1u ) );

        fitnessMagnitudeAverageHistory.push_back( fitnessMagnitudeAverage );
        //fitnessMagnitudeTopIndividuals = bislip::Variables::getFitnessMagnitudeOfTopIndividuals( bodyMap, problemInput, topIndividuals );
        fitnessMagnitudeAverageTopIndividualsVectorHistory.push_back( fitnessMagnitudeTopIndividualsVectorAverage );
        fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsVectorHistory ) );
        fitnessMagnitudeAverageTopIndividualsMagnitudeHistory.push_back( fitnessMagnitudeTopIndividualsMagnitudeAverage );
        fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory ) );

        temp( 0 ) = fitnessMagnitudeAverageHistory[ 0 ];
        temp( 1 ) = fitnessMagnitudeAverageTopIndividualsVectorHistory[ 0 ];
        temp( 2 ) = TUDAT_NAN;
        temp( 3 ) = fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ 0 ];
        temp( 4 ) = TUDAT_NAN;
        temp( 5 ) = fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ 0 ];
        temp( 6 ) = TUDAT_NAN;
        temp( 7 ) = fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ 0 ];
        temp( 8 ) = TUDAT_NAN;

        fitnessMagnitudeAveragesOutputHistory[ 0 ] = temp;

        //! Write population and fitness to file
        if( bool( int( outputSettingsValues[ 0 ] ) ) )
        {
            std::cout<< "Printing population and fitness" << std::endl;
            tudat::input_output::writeDataMapToTextFile( problemInput->getPopulationHistory(), "populationHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( problemInput->getFitnessHistory(), "fitnessHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( problemInput->getExtremesAndConstraintsHistory(), "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            std::cout<< " " << std::endl;
        }
    }
    else if( bislipSystems->getFlag( bislip::flags::flag_list::test_case ) )
    {
        std::cout << "Running a Test case." << std::endl;
        std::cout << "  Population Size         = " << problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::population_per_island ) << std::endl;
        std::cout << "  Number of Evolutions    = " << 0 << std::endl;
        std::cout << "  Propagation Step Size   = " << propagationStepSize << std::endl;
        std::cout << "  Guidance Step Size      = " << guidanceStepSize << std::endl;
        std::cout << "  " << std::endl;

        //! Assign population size.
        pagmo::population::size_type populationSizePerIsland = problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::population_per_island );
        if( debugInfo == 1 ){ std::cout << "  Initialize seed generator for reproducible results." << std::endl; }
        pagmo::random_device::set_seed( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::seed_generator_initializer ) );

        //! Assign number of evolutions
        unsigned int evolutions = problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::number_of_evolutions );

        if( debugInfo == 1 ){ std::cout << "  Evaluate PRNG and store seeds used throughout the evolutions" << std::endl; }
        std::vector< pagmo::detail::random_engine_type::result_type > seedSequenceVector;
        std::generate_n( std::back_inserter( seedSequenceVector ), evolutions + 1, pagmo::random_device::next );
        problemInput->setSeedSequenceVector( seedSequenceVector );

        //! Create object to compute the problem fitness.
        if( debugInfo == 1 ){ std::cout << "  Create Pagmo Problem" << std::endl; }
        pagmo::problem prob{ bislip::Space4ErrBody( problemInput, bodyMap ) };

        if( debugInfo == 1 ){ std::cout << "    Generate initial population (single individual)" << std::endl; }
        problemInput->setIndividualNumber( 0 );
        problemInput->setGenerationNumber( 0 );
        bislipSystems->setFlag( bislip::flags::flag_list::evolution_evaluation, true );
        pagmo::population population{ prob, populationSizePerIsland, static_cast< unsigned int >( seedSequenceVector[ 0 ] ) };

        unsigned int topIndividuals = static_cast< unsigned int >( populationSizePerIsland );
        if( debugInfo == 1 ){ std::cout << "    'Ranking' individuals" << std::endl; }

        bislip::Variables::rankFitnessVectors( bodyMap, problemInput, topIndividuals, fitnessMagnitudeAverage, fitnessMagnitudeTopIndividualsVectorAverage, fitnessMagnitudeTopIndividualsMagnitudeAverage );

        fitnessMagnitudeAverageHistory.push_back( fitnessMagnitudeAverage );
        fitnessMagnitudeAverageTopIndividualsVectorHistory.push_back( fitnessMagnitudeTopIndividualsVectorAverage );
        fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsVectorHistory ) );
        fitnessMagnitudeAverageTopIndividualsMagnitudeHistory.push_back( fitnessMagnitudeTopIndividualsMagnitudeAverage );
        fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory ) );

        temp( 0 ) = fitnessMagnitudeAverageHistory[ 0 ];
        temp( 1 ) = fitnessMagnitudeAverageTopIndividualsVectorHistory[ 0 ];
        temp( 2 ) = TUDAT_NAN;
        temp( 3 ) = fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ 0 ];
        temp( 4 ) = TUDAT_NAN;
        temp( 5 ) = fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ 0 ];
        temp( 6 ) = TUDAT_NAN;
        temp( 7 ) = fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ 0 ];
        temp( 8 ) = TUDAT_NAN;

        fitnessMagnitudeAveragesOutputHistory[ 0 ] = temp;

        bislip::Variables::reruns( bodyMap, problemInput, bislip::Variables::selectTopIndividuals( bodyMap, problemInput, topIndividuals ) );

        //! Write population and fitness to file
        if( bool( int( outputSettingsValues[ 0 ] ) ) )
        {
            std::cout<< "Printing population and fitness" << std::endl;
            tudat::input_output::writeDataMapToTextFile( problemInput->getPopulationHistory(), "populationHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( problemInput->getFitnessHistory(), "fitnessHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( problemInput->getExtremesAndConstraintsHistory(), "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            std::cout<< " " << std::endl;
        }

    }
    else
    {
        std::cout << "  Main Optimization Section" << std::endl;

        unsigned int k = 0;
        problemInput->setIndividualNumber( 0 );
        problemInput->setGenerationNumber( k );
        bislipSystems->setFlag( bislip::flags::flag_list::evolution_evaluation, true );

        bislip::createTUDATSettings::createDependentVariableList( problemInput, bodyMap );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PAGMO PROBLEM TO COMPUTE FITNESS            ////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //bislipSystems->setDebugInfo( 1 );
debugInfo = 1;
        if( debugInfo == 1 ){ std::cout << "Starting creation of Pagmo Problem" << std::endl; }

        //! Assign population size.
        pagmo::population::size_type populationSizePerIsland = problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::population_per_island );

        if( debugInfo == 1 ){ std::cout << "  Initialize seed generator for reproducible results." << std::endl; }
        pagmo::random_device::set_seed( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::seed_generator_initializer ) );

        //! Assign number of evolutions
        unsigned int evolutions = problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::number_of_evolutions );

        if( debugInfo == 1 ){ std::cout << "  Evaluate PRNG and store seeds used throughout the evolutions" << std::endl; }
        std::vector< pagmo::detail::random_engine_type::result_type > seedSequenceVector;
        std::generate_n( std::back_inserter( seedSequenceVector ), evolutions + 1, pagmo::random_device::next );
        problemInput->setSeedSequenceVector( seedSequenceVector );

        //for( unsigned int i = 0; i < seedSequenceVector.size(); i++ ) { std::cout << seedSequenceVector[ i ] << std::endl; }

        if( debugInfo == 1 ){ std::cout << "  Reset Seed Generator" << std::endl; }
        pagmo::random_device::set_seed( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::seed_generator_initializer ) );

        //! Create object to compute the problem fitness.
        if( debugInfo == 1 ){ std::cout << "  Create Pagmo Problem" << std::endl; }
        pagmo::problem prob{ bislip::Space4ErrBody( problemInput, bodyMap ) };

        std::cout << "--------------------------------------------------------" << std::endl;
        std::cout << "Creating initial population based off problem, propulation size and seed number" << std::endl;
        std::cout << " " << std::endl;
        //! Create initial population based off problem, propulation size and seed number
        if( int( outputSettingsValues[ 2 ] ) == 1 ) { populationSizePerIsland = 1; }
        std::cout << "  Population per island: " << populationSizePerIsland << std::endl;
        std::cout << "  Seed for individuals:  " << seedSequenceVector[ 0 ] << std::endl;
        std::cout<< "  " << std::endl;

        pagmo::population population{ prob, populationSizePerIsland, static_cast< unsigned int >( seedSequenceVector[ 0 ] ) };
        pagmo::population populationImport( population );
        std::cout<< "  " << std::endl;
        std::cout<< "Ranking Population" << std::endl;
        std::cout<< "  " << std::endl;

        unsigned int topIndividuals = static_cast< unsigned int >( outputSettingsValues[ 1 ] );
        if( topIndividuals == 0u ) { topIndividuals = 1u; }

        bislip::Variables::rankFitnessVectors( bodyMap, problemInput, topIndividuals, fitnessMagnitudeAverage, fitnessMagnitudeTopIndividualsVectorAverage, fitnessMagnitudeTopIndividualsMagnitudeAverage );
        std::cout<< "  " << std::endl;

        //! Retrieve algorithm. Three options available in the following function:
        //!        getMultiObjectiveAlgorithm
        //!               case 0 --> nsga2
        //!               case 1 --> moead
        //!               case 2 --> ihs
        //! moead.hpp should be modified such that the points are generated with a Halton sequence.
        //pagmo::algorithm algo{ getPagmoAlgorithm( optimizerSpecificParameterValues ) };
        std::cout<< "Getting algorithm" << std::endl;
        pagmo::algorithm algo{ bislip::getPagmoAlgorithm( optimizerCustomSettingsFlag, optimizationSettings.tail( optimizationSettings.size() - 7 ) ) };
        std::cout<< "  " << std::endl;

        //! Create island
        //pagmo::island island{ algo, population };

        //pagmo::island island{ algo, prob, populationSizePerIsland };

        //std::vector< std::vector < double > > fitnessMagnitudeAveragesOutputHistory( 1, std::vector< double >( 6 ) );


        // Store the string representation of p.
        std::stringstream ss_algo, ss_prob, ss_population;
        std::stringstream ss_populationHistory, ss_fitnessHistory;
        std::map< std::string, Eigen::VectorXd > populationHistory, fitnessHistory, extremesAndconstraintsHistory;
        std::map< std::string, Eigen::VectorXd > populationMap = problemInput->getPopulation();
        std::map< std::string, Eigen::VectorXd > clearPopulationMap = populationMap;
        clearPopulationMap.clear();
        std::map< std::string, Eigen::VectorXd > fitnessMap = problemInput->getFitness();
        std::map< std::string, Eigen::VectorXd > clearFitnessMap = fitnessMap;
        clearFitnessMap.clear();
        std::map< std::string, Eigen::VectorXd > extremesMap = problemInput->getExtremesAndConstraints();
        std::map< std::string, Eigen::VectorXd > clearExtremesMap = extremesMap;
        clearExtremesMap.clear();
        std::map< std::string, std::vector< double > > topIndividualsMap;


        if( resumeOptimization || importInitialPopulation )
        {

            std::vector< std::string > individualNames;
            unsigned int purePopulationColumns = static_cast< unsigned int >( ( populationMap.begin()->second ).size() ) - 5u;
            //unsigned int pureFitnessColumns = static_cast< unsigned int >( ( fitnessHistory.begin()->second ).size() ) - 5u;
            //std::cout << "Original Individual 1: " <<  populationMap.begin()->second << std::endl;


            problemInput->setPopulation( clearPopulationMap );
            problemInput->setFitness( clearFitnessMap );
            problemInput->setPopulationHistory( clearPopulationMap );
            problemInput->setFitnessHistory( clearFitnessMap );
            problemInput->setExtremesAndConstraints( clearExtremesMap );
            problemInput->setExtremesAndConstraintsHistory( clearExtremesMap );

            if( ( int( outputSettingsValues[ 2 ] ) == 0 ) && resumeOptimization )
            {
                std::cout << "Restarting an optimization run. Must read in various files and overwrite the population." << std::endl;
                std::cout << "  Deserializing (importing) PAGMO population." << std::endl;
                std::ifstream is( outputPath + outputSubFolder + "/archiveJSON_Population.dat", std::ios::binary );
                { // load
                    cereal::PortableBinaryInputArchive archiveLoad( is );
                    archiveLoad( populationImport );
                }

                std::cout << "      Replacing initialized population with imported population." << std::endl;
                for ( pagmo::population::size_type i = 0; i < populationSizePerIsland; i++ )
                {
                    population.set_xf( i, populationImport.get_x()[ i ], populationImport.get_f()[ i ] );
                }

                std::cout << "  Loading and setting Population/Fitness history" << std::endl;
                populationHistory = bislip::utilities::readFileToDataMapWithStringKey( "populationHistory.dat", outputPath + outputSubFolder );
                fitnessHistory = bislip::utilities::readFileToDataMapWithStringKey( "fitnessHistory.dat", outputPath + outputSubFolder );
                extremesAndconstraintsHistory = bislip::utilities::readFileToDataMapWithStringKey( "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder );

                if( debugInfo == 1 ){ std::cout << "Extract names" << std::endl; }
                individualNames.clear();
                for ( auto const& element : populationHistory ) { individualNames.push_back( element.first ); }

                purePopulationColumns = static_cast< unsigned int >( ( populationHistory.begin()->second ).size() ) - 5u;
                //pureFitnessColumns = static_cast< unsigned int >( ( fitnessHistory.begin()->second ).size() ) - 5u;

                problemInput->setFitnessHistory( fitnessHistory );
                problemInput->setFitness( fitnessHistory );
                problemInput->setPopulationHistory( populationHistory );
                problemInput->setPopulation( populationHistory );
                problemInput->setExtremesAndConstraints( extremesAndconstraintsHistory );
                problemInput->setExtremesAndConstraintsHistory( extremesAndconstraintsHistory );

                if( bool( int( outputSettingsValues[ 2 ] ) ) )
                {
                    std::cout << "  The interest is to re-evaluate with the purpose of printing them all. So we must make sure that the " << std::endl;
                    std::cout <<             "number of top individals is the same as the population size before 'ranking' them." << std::endl;
                    std::cout << "  " << std::endl;
                    std::cout << "  Rank imported population." << std::endl;
                    bislip::Variables::rankFitnessVectors( bodyMap, problemInput, populationHistory.size(), fitnessMagnitudeAverage, fitnessMagnitudeTopIndividualsVectorAverage, fitnessMagnitudeTopIndividualsMagnitudeAverage );
                }

                std::cout << "  Determine the generation number of loaded population" << std::endl;
                std::vector< double > distinctGenerationNumberList;

                for( unsigned int i = 0; i < static_cast< unsigned int >( fitnessHistory.size() ); i++ )
                { distinctGenerationNumberList.push_back( fitnessHistory.at( individualNames[ i ] )( 3 ) ); }
                std::sort( distinctGenerationNumberList.begin(), distinctGenerationNumberList.end() );
                auto last = std::unique( distinctGenerationNumberList.begin(), distinctGenerationNumberList.end() );
                distinctGenerationNumberList.erase( last, distinctGenerationNumberList.end() );
                k = static_cast< unsigned int >( *std::max_element( distinctGenerationNumberList.begin(), distinctGenerationNumberList.end() ) );

                if( bool( int( outputSettingsValues[ 2 ] ) ) )
                { pagmo::random_device::set_seed( seedSequenceVector[ k ] ); }

                std::cout << "  Loading fitness magnitude average history and populating associated vectors" << std::endl;
                fitnessMagnitudeAveragesOutputHistory = tudat::input_output::readEigenVectorMapFromFile( outputPath + outputSubFolder + "/fitnessMagnitudeAveragesOutputHistory.dat" );

                Eigen::MatrixXd fitnessMagnitudeAveragesOutputHistoryMatrix( fitnessMagnitudeAveragesOutputHistory.size(), ( fitnessMagnitudeAveragesOutputHistory.begin()->second ).size() );

                for ( unsigned int i = 0; i < fitnessMap.size(); i++ )
                { fitnessMagnitudeAveragesOutputHistoryMatrix.row( i ) = fitnessMagnitudeAveragesOutputHistory.at( i ); }

                fitnessMagnitudeAverageHistory.clear();
                fitnessMagnitudeAverageTopIndividualsVectorHistory.clear();
                fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory.clear();
                fitnessMagnitudeAverageTopIndividualsMagnitudeHistory.clear();
                fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory.clear();
                fitnessMagnitudeAverageHistory = bislip::utilities::convertEigenVectorXdToSTLVector( fitnessMagnitudeAveragesOutputHistoryMatrix.col( 0 ) );
                fitnessMagnitudeAverageTopIndividualsVectorHistory = bislip::utilities::convertEigenVectorXdToSTLVector( fitnessMagnitudeAveragesOutputHistoryMatrix.col( 1 ) );
                fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory = bislip::utilities::convertEigenVectorXdToSTLVector( fitnessMagnitudeAveragesOutputHistoryMatrix.col( 3 ) );
                fitnessMagnitudeAverageTopIndividualsMagnitudeHistory = bislip::utilities::convertEigenVectorXdToSTLVector( fitnessMagnitudeAveragesOutputHistoryMatrix.col( 5 ) );
                fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory = bislip::utilities::convertEigenVectorXdToSTLVector( fitnessMagnitudeAveragesOutputHistoryMatrix.col( 7 ) );

                if( bool( int( outputSettingsValues[ 2 ] ) ) )
                {
                    topIndividualsMap = bislip::Variables::selectTopIndividuals( bodyMap, problemInput, populationHistory.size() );

                    std::cout<< "Printing population, fitness, extremes, and fitness magnitude history with printed flag" << std::endl;
                    tudat::input_output::writeDataMapToTextFile( problemInput->getPopulationHistory(), "populationHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    tudat::input_output::writeDataMapToTextFile( problemInput->getFitnessHistory(), "fitnessHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    tudat::input_output::writeDataMapToTextFile( problemInput->getExtremesAndConstraintsHistory(), "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    std::cout<< " " << std::endl;

                    bislip::Variables::reruns( bodyMap, problemInput, topIndividualsMap);
                }
            }

            if( importInitialPopulation )
            {

                std::map< std::string, Eigen::VectorXd > initialPopulationImportMap;
                std::cout << " Loading and setting initial Population" << std::endl;
                std::cout << "    Initial Population path and file name: " << initialPopulationFile << std::endl;
                initialPopulationImportMap = bislip::utilities::readFileToDataMapWithStringKey( initialPopulationFile, "" );

                unsigned int pureImportedPopulationColumns = static_cast< unsigned int >( ( initialPopulationImportMap.begin()->second ).size() ) - 5u;


                Eigen::MatrixXd initialPopulationImportMatrix( initialPopulationImportMap.size(), pureImportedPopulationColumns );

                if( debugInfo == 1 ){ std::cout << "    Extract names" << std::endl; }
                individualNames.clear();
                for ( auto const& element : initialPopulationImportMap ) { individualNames.push_back( element.first ); }

                for ( unsigned int i = 0; i < initialPopulationImportMap.size(); i++ )
                { initialPopulationImportMatrix.row( i ) = initialPopulationImportMap.at( individualNames[ i ] ).tail( pureImportedPopulationColumns ); }

                //std::cout << "Imported Individual 1: " <<  initialPopulationImportMatrix.row( 0 ) << std::endl;

                if( pureImportedPopulationColumns != purePopulationColumns ) { bislip::utilities::removeColumnFromEigenMatrix( initialPopulationImportMatrix, 9 ); }

                //std::cout << "Imported Individual 1: " <<  initialPopulationImportMatrix.row( 0 ) << std::endl;

                //pureFitnessColumns = static_cast< unsigned int >( ( fitnessHistory.begin()->second ).size() ) - 5u;
                problemInput->setIndividualNumber( 0 );
                problemInput->setGenerationNumber( 0 );
                bislipSystems->setFlag( bislip::flags::flag_list::impose_parameter_bounds_on_interpolators, false );

                std::vector< double > individual;
                std::cout << "    Replacing initialized population with pre-determined initial population." << std::endl;
                //std::cout << "    ";

                for ( pagmo::population::size_type i = 0; i < populationSizePerIsland; i++ )
                {
                    individual.clear();
                    individual = bislip::utilities::convertEigenVectorXdToSTLVector( initialPopulationImportMatrix.row( i ) );
                    population.set_x( i, individual );
                }
                //std::cout << "100%"<< std::endl;
                bislipSystems->setFlag( bislip::flags::flag_list::impose_parameter_bounds_on_interpolators, true );

                bislip::Variables::rankFitnessVectors( bodyMap, problemInput, topIndividuals, fitnessMagnitudeAverage, fitnessMagnitudeTopIndividualsVectorAverage, fitnessMagnitudeTopIndividualsMagnitudeAverage );

                if( debugInfo == 1 ){ std::cout << "  Reset Seed Generator" << std::endl; }
                pagmo::random_device::set_seed( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::seed_generator_initializer ) );

                fitnessMagnitudeAverageHistory.push_back( fitnessMagnitudeAverage );
                //fitnessMagnitudeTopIndividuals = bislip::Variables::getFitnessMagnitudeOfTopIndividuals( bodyMap, problemInput, topIndividuals );
                fitnessMagnitudeAverageTopIndividualsVectorHistory.push_back( fitnessMagnitudeTopIndividualsVectorAverage );
                fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsVectorHistory ) );
                fitnessMagnitudeAverageTopIndividualsMagnitudeHistory.push_back( fitnessMagnitudeTopIndividualsMagnitudeAverage );
                fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory ) );

                temp( 0 ) = fitnessMagnitudeAverageHistory[ k ];
                temp( 1 ) = fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ];
                temp( 2 ) = TUDAT_NAN;
                temp( 3 ) = fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ];
                temp( 4 ) = TUDAT_NAN;
                temp( 5 ) = fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ];
                temp( 6 ) = TUDAT_NAN;
                temp( 7 ) = fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ];
                temp( 8 ) = TUDAT_NAN;

                fitnessMagnitudeAveragesOutputHistory[ k ] = temp;
                std::cout<< "  " << std::endl;

                std::cout<< "Identifying and Storing Top Individuals from this generation" << std::endl;
                std::cout<< "  " << std::endl;
                topIndividualsMap = bislip::Variables::selectTopIndividuals( bodyMap, problemInput, topIndividuals );

                //! Write original (unevolved) population and fitness to file
                if( bool( int( outputSettingsValues[ 0 ] ) ) )
                {
                    std::cout<< "Printing of initial population, fitness, extremes, and fitness magnitude history" << std::endl;
                    tudat::input_output::writeDataMapToTextFile( problemInput->getPopulationHistory(), "populationHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    tudat::input_output::writeDataMapToTextFile( problemInput->getFitnessHistory(), "fitnessHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    tudat::input_output::writeDataMapToTextFile( problemInput->getExtremesAndConstraintsHistory(), "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                    std::cout<< "Archiving initial population" << std::endl;
                    std::ofstream os( outputPath + outputSubFolder + "/archiveJSON_Population.dat", std::ios::binary );
                    { // save
                        cereal::PortableBinaryOutputArchive archiveSave( os );
                        archiveSave( population );
                    }
                    std::cout<< " " << std::endl;
                }

                bislip::Variables::reruns( bodyMap, problemInput, topIndividualsMap );
            }
        }
        else
        {
            /*unsigned int p = 0;
                    for( auto const& fitnessVector : fitnessMap)
                    {
                        fitnessMagnitude( p ) = ( fitnessVector.second ).tail( fitnessVector.second.size() - 5 ).norm();
                        p += 1;
                    }
        */
            fitnessMagnitudeAverageHistory.push_back( fitnessMagnitudeAverage );
            //fitnessMagnitudeTopIndividuals = bislip::Variables::getFitnessMagnitudeOfTopIndividuals( bodyMap, problemInput, topIndividuals );
            fitnessMagnitudeAverageTopIndividualsVectorHistory.push_back( fitnessMagnitudeTopIndividualsVectorAverage );
            fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsVectorHistory ) );
            fitnessMagnitudeAverageTopIndividualsMagnitudeHistory.push_back( fitnessMagnitudeTopIndividualsMagnitudeAverage );
            fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory ) );

            temp( 0 ) = fitnessMagnitudeAverageHistory[ k ];
            temp( 1 ) = fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ];
            temp( 2 ) = TUDAT_NAN;
            temp( 3 ) = fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ];
            temp( 4 ) = TUDAT_NAN;
            temp( 5 ) = fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ];
            temp( 6 ) = TUDAT_NAN;
            temp( 7 ) = fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ];
            temp( 8 ) = TUDAT_NAN;

            fitnessMagnitudeAveragesOutputHistory[ k ] = temp;
            std::cout<< "  " << std::endl;

            std::cout<< "Identifying and Storing Top Individuals from this generation" << std::endl;
            std::cout<< "  " << std::endl;
            topIndividualsMap = bislip::Variables::selectTopIndividuals( bodyMap, problemInput, topIndividuals );

            //! Write original (unevolved) population and fitness to file
            if( bool( int( outputSettingsValues[ 0 ] ) ) )
            {
                std::cout<< "Printing of initial population, fitness, extremes, and fitness magnitude history" << std::endl;
                tudat::input_output::writeDataMapToTextFile( problemInput->getPopulationHistory(), "populationHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                tudat::input_output::writeDataMapToTextFile( problemInput->getFitnessHistory(), "fitnessHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                tudat::input_output::writeDataMapToTextFile( problemInput->getExtremesAndConstraintsHistory(), "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                std::cout<< "Archiving initial population" << std::endl;
                std::ofstream os( outputPath + outputSubFolder + "/archiveJSON_Population.dat", std::ios::binary );
                { // save
                    cereal::PortableBinaryOutputArchive archiveSave( os );
                    archiveSave( population );
                }
                std::cout<< " " << std::endl;
            }


            if( evolutions == 0 )
            {
                problemInput->setPopulation( clearPopulationMap );
                problemInput->setFitness( clearFitnessMap );
                problemInput->setPopulationHistory( clearPopulationMap );
                problemInput->setFitnessHistory( clearFitnessMap );
                problemInput->setExtremesAndConstraints( clearExtremesMap );
                problemInput->setExtremesAndConstraintsHistory( clearExtremesMap );

            }

            bislip::Variables::reruns( bodyMap, problemInput, topIndividualsMap );


            //prob = pagmo::problem{ bislip::Space4ErrBody( problemInput, bodyMap ) };

            //algo = pagmo::algorithm{ getPagmoAlgorithm( optimizerCustomSettingsFlag, optimizationSettings.tail( optimizationSettings.size() - 7 ) ) };

            // population = pagmo::population{ prob, populationSizePerIsland };

            //island.set_algorithm( algo );
            //island.set_population( population );

            std::cout << " " << std::endl;
            std::cout << "Initial population created" << std::endl;
            std::cout << " " << std::endl;

        }

        if( evolutions == 0 )
        {
            std::cout << " " << std::endl;
            std::cout << "Number of evolutions to run is equal to 0. Nothing else to do." << std::endl;
            std::cout << " " << std::endl;

        }
        else
        {
            std::cout << " " << std::endl;
            std::cout << "Number of evolutions to run is greater than 0." << std::endl;
            std::cout << " " << std::endl;


            //! Assign archipelago number.
            pagmo::archipelago::size_type archipelagoSize = problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::archipelago_size );

            //! Check if PAGMO is to use a Single Island or an Archipelago
            bool singleIsland = bool( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::single_island_flag ) );

            //! Check if it is a MonteCarlo analysis
            bool monteCarloSimulation = bool( problemInput->getOptimizationSettingsValue( bislip::optimization_settings::optimization_settings_list::monte_carlo_flag ) );

            evolutions += k;
            k += 1;
            bool exitCondition1 = true;
            bool exitCondition2 = true;
            bool exitCondition3 = true;
            bool exitCondition4 = true;
            bool exitCondition5 = true;
            bool exitCondition6 = true;
            bool exitCondition7 = true;
            bool exitCondition8 = true;
            bool rerunCondition = false;

            std::cout << "  " << std::endl;
            std::cout << "Running an Optimization problem." << std::endl;
            std::cout << "  Individuals per Island  = " << populationSizePerIsland << std::endl;
            if( singleIsland == true ) { std::cout << "  Number of Islands       = 1"  << std::endl; }
            else { std::cout << "  Number of Islands       = " << archipelagoSize << std::endl; }
            std::cout << "  Starting Generation     = " << k - 1 << std::endl;
            std::cout << "  Ending Generation       = " << evolutions << std::endl;
            std::cout << "  Propagation Step Size   = " << propagationStepSize << std::endl;
            std::cout << "  Guidance Step Size      = " << guidanceStepSize << std::endl;
            std::cout << "  Objective Function Case = " << objectiveFunctionCase << std::endl;
            std::cout << "  " << std::endl;

            tudat::input_output::writeDataMapToTextFile( problemInput->getPopulationHistory(), "populationHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( problemInput->getFitnessHistory(), "fitnessHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( problemInput->getExtremesAndConstraintsHistory(), "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );

            if( singleIsland && optimizerIndex != 2 )
            {
                do
                {
                    problemInput->setIndividualNumber( 0 );
                    problemInput->setGenerationNumber( k );
                    problemInput->setPopulation( clearPopulationMap );
                    problemInput->setFitness( clearFitnessMap );
                    //problemInput->setEvolutionEvaluationFlag( true );
                    bislipSystems->setFlag( bislip::flags::flag_list::evolution_evaluation, true );
                    topIndividualsMap.clear();

                    bislip::createTUDATSettings::createDependentVariableList( problemInput, bodyMap );

                    population = algo.evolve( population );

                    if( ( problemInput->getListOfIndividualsToReplace() ).empty() == 0 )
                    {
                        for ( auto const& element : problemInput->getListOfIndividualsToReplace() ) { population.set_xf( element.first, element.second.first, element.second.second ); }
                        problemInput->clearListOfIndividualsToReplace();
                    }

                    std::cout<< "  " << std::endl;
                    std::cout<< "Ranking Population" << std::endl;

                    bislip::Variables::rankFitnessVectors( bodyMap, problemInput, topIndividuals, fitnessMagnitudeAverage, fitnessMagnitudeTopIndividualsVectorAverage, fitnessMagnitudeTopIndividualsMagnitudeAverage );

                    /*
                    island.evolve( );
                    while( island.status( ) != pagmo::evolve_status::idle &&
                           island.status( ) != pagmo::evolve_status::idle_error )
                    {
                        island.wait( );
                    }
                    island.wait_check( ); // Raises errors
                     */

                    std::cout<< "  " << std::endl;
                    std::cout<< "Generation " << k << "/" << evolutions <<" created" << std::endl;
                    std::cout<< "  " << std::endl;

                    fitnessMagnitudeAverageHistory.push_back( fitnessMagnitudeAverage );
                    fitnessMagnitudeAverageTopIndividualsVectorHistory.push_back( fitnessMagnitudeTopIndividualsVectorAverage );
                    fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsVectorHistory ) );
                    fitnessMagnitudeAverageTopIndividualsMagnitudeHistory.push_back( fitnessMagnitudeTopIndividualsMagnitudeAverage );
                    fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory ) );

                    temp( 0 ) = fitnessMagnitudeAverageHistory[ k ];
                    temp( 1 ) = fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ];
                    temp( 2 ) = 1000 * std::abs( ( fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ] - fitnessMagnitudeAverageTopIndividualsVectorHistory[ k - 1 ] ) / fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ] );
                    temp( 3 ) = fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ];
                    temp( 4 ) = 1000 * std::abs( ( fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ] - fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k - 1 ] ) / fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ] );
                    temp( 5 ) = fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ];
                    temp( 6 ) = 1000 * std::abs( ( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ] - fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k - 1 ] ) / fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ] );
                    temp( 7 ) = fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ];
                    temp( 8 ) = 1000 * std::abs( ( fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ] - fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k - 1 ] ) / fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ] );

                    fitnessMagnitudeAveragesOutputHistory[ k ] = temp;

                    std::cout<< "  " << std::endl;
                    std::cout<< "Identifying and Storing Top Individuals from this generation" << std::endl;
                    topIndividualsMap = bislip::Variables::selectTopIndividuals( bodyMap, problemInput, topIndividuals );

                    if  ( int( outputSettingsValues[ 0 ] ) == 1 )
                    {
                        std::cout<< "Printing of population, fitness, extremes, and fitness magnitude history" << std::endl;
                        bislip::utilities::writeDataMapToExistingTextFile( problemInput->getPopulation(),"populationHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                        bislip::utilities::writeDataMapToExistingTextFile( problemInput->getFitness(),"fitnessHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                        bislip::utilities::writeDataMapToExistingTextFile( problemInput->getExtremesAndConstraints(),"extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                        tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );

                        std::cout<< "Archiving population" << std::endl;
                        std::ofstream os( outputPath + outputSubFolder + "/archiveJSON_Population.dat", std::ios::binary );
                        { // save
                            cereal::PortableBinaryOutputArchive archiveSave( os );
                            archiveSave( population );
                        }
                    }

                    rerunCondition = ( k == static_cast< unsigned int >( evolutions / 4  + 0.5 ) ) || ( k == static_cast< unsigned int >( 2 * evolutions / 4 + 0.5 ) ) || ( k == static_cast< unsigned int >( 3 * evolutions / 4 + 0.5 ) );

                    if( rerunCondition )
                    {
                        std::cout<< "Trajectory data printing conditional met: " << "k = " << k << " == ";

                        if( k == static_cast< unsigned int >( evolutions / 4  + 0.5 ) )
                        { std::cout << static_cast< unsigned int >( evolutions / 4  + 0.5 ) << std::endl; }
                        else if( k == static_cast< unsigned int >( 2 * evolutions / 4  + 0.5 ) )
                        { std::cout << static_cast< unsigned int >( 2 * evolutions / 4  + 0.5 ) << std::endl; }
                        else if( k == static_cast< unsigned int >( 3 * evolutions / 4  + 0.5 ) )
                        { std::cout << static_cast< unsigned int >( 3 * evolutions / 4  + 0.5 ) << std::endl; }

                        bislip::Variables::reruns( bodyMap, problemInput, topIndividualsMap );

                    }

                    exitCondition1 = fitnessMagnitudeAveragesOutputHistory.at( k )( 2 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 2 ) ) == true )
                    { exitCondition2 = true; }
                    else
                    { exitCondition2 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 2 ) >= 1; }
                    exitCondition3 = fitnessMagnitudeAveragesOutputHistory.at( k )( 4 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 4 ) ) == true )
                    { exitCondition4 = true; }
                    else
                    { exitCondition4 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 4 ) >= 1; }
                    exitCondition5 = fitnessMagnitudeAveragesOutputHistory.at( k )( 6 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 6 ) ) == true )
                    { exitCondition6 = true; }
                    else
                    { exitCondition6 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 6 ) >= 1; }
                    exitCondition7 = fitnessMagnitudeAveragesOutputHistory.at( k )( 8 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 8 ) ) == true )
                    { exitCondition8 = true; }
                    else
                    { exitCondition8 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 8 ) >= 1; }


                    k += 1;

                    std::cout<< "  " << std::endl;
                    std::cout<< "  " << std::endl;
                } while( ( k < static_cast< unsigned int >( 3 * evolutions / 4 + 0.5 ) ) || ( ( k < evolutions + 1 ) && ( ( exitCondition1 && exitCondition2 ) || ( exitCondition3 && exitCondition4 ) || ( exitCondition5 && exitCondition6 ) || ( exitCondition7 && exitCondition8 ) ) ) );
            }
            else
            {
                evolutions *= populationSizePerIsland;
                std::vector< std::string > suffixes;

                do
                {
                    problemInput->setGenerationNumber( k );
                    bislipSystems->setFlag( bislip::flags::flag_list::evolution_evaluation, true );
                    topIndividualsMap.clear();

                    bislip::createTUDATSettings::createDependentVariableList( problemInput, bodyMap );

                    population = algo.evolve( population );

                    if( ( problemInput->getListOfIndividualsToReplace() ).empty() == 0 )
                    {
                        for ( auto const& element : problemInput->getListOfIndividualsToReplace() ) { population.set_xf( element.first, element.second.first, element.second.second ); }
                        problemInput->clearListOfIndividualsToReplace();
                    }

                    /*for (unsigned int i = 0; i != population.get_x().size(); ++i)
                    {
                        std::cout << "Individual " << i << " : ";
                        for (unsigned int j = 0; j != population.get_x()[i].size(); ++j )
                        { std::cout << population.get_x()[i][j] << " "; }

                        std::cout << std::endl;
                    }*/

                    for( auto const& element : problemInput->getPopulation() ) { suffixes.push_back( element.first ); }

                    for( unsigned int i = 0; i < suffixes.size(); i++ )
                    {
                        problemInput->updateIndividualGeneration( suffixes[ i ], k );
                    }

                    //! Rank individuals. They will be ranked by vector and by magnitude.
                    bislip::Variables::rankFitnessVectors( bodyMap, problemInput, topIndividuals, fitnessMagnitudeAverage, fitnessMagnitudeTopIndividualsVectorAverage, fitnessMagnitudeTopIndividualsMagnitudeAverage );

                    //std::cout << " " << std::endl;
                    //std::cout << "Remove lowest ranked individual" << std::endl;
                    for( unsigned int i = 0; i < suffixes.size(); i++ )
                    {
                        if( static_cast< unsigned int >( problemInput->getFitness().at( suffixes[ i ] )( 1 ) ) == populationSizePerIsland )
                        {
                            // std::cout << "  Lowest ranked individual found" << std::endl;
                            problemInput->eraseIndividual( suffixes[ i ] );
                            // std::cout << "  Lowest ranked individual removed" << std::endl;
                            //suffixes.erase( suffixes.begin() + i - 1 );
                            suffixes.erase( std::remove( suffixes.begin(), suffixes.end(), suffixes[ i ] ), suffixes.end() );
                            break;
                        }
                    }


                    // std::cout << " " << std::endl;
                    //  std::cout << "Re-rank individuals to re-order magnitude ranking" << std::endl;
                    bislip::Variables::rankFitnessVectors( bodyMap, problemInput, topIndividuals, fitnessMagnitudeAverage, fitnessMagnitudeTopIndividualsVectorAverage, fitnessMagnitudeTopIndividualsMagnitudeAverage );

                    /*
                    island.evolve( );
                    while( island.status( ) != pagmo::evolve_status::idle &&
                           island.status( ) != pagmo::evolve_status::idle_error )
                    {
                        island.wait( );
                    }
                    island.wait_check( ); // Raises errors
                     */

                    // std::cout<< "  " << std::endl;
                    // std::cout<< "Generation " << k << "/" << evolutions <<" created" << std::endl;
                    // std::cout<< "  " << std::endl;

                    fitnessMagnitudeAverageHistory.push_back( fitnessMagnitudeAverage );
                    fitnessMagnitudeAverageTopIndividualsVectorHistory.push_back( fitnessMagnitudeTopIndividualsVectorAverage );
                    fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsVectorHistory ) );
                    fitnessMagnitudeAverageTopIndividualsMagnitudeHistory.push_back( fitnessMagnitudeTopIndividualsMagnitudeAverage );
                    fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory.push_back( bislip::utilities::computeAverageOfSTLVector( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory ) );

                    temp( 0 ) = fitnessMagnitudeAverageHistory[ k ];
                    temp( 1 ) = fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ];
                    temp( 2 ) = 1000 * std::abs( ( fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ] - fitnessMagnitudeAverageTopIndividualsVectorHistory[ k - 1 ] ) / fitnessMagnitudeAverageTopIndividualsVectorHistory[ k ] );
                    temp( 3 ) = fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ];
                    temp( 4 ) = 1000 * std::abs( ( fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ] - fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k - 1 ] ) / fitnessMagnitudeAverageOfAverageTopIndividualsVectorHistory[ k ] );
                    temp( 5 ) = fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ];
                    temp( 6 ) = 1000 * std::abs( ( fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ] - fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k - 1 ] ) / fitnessMagnitudeAverageTopIndividualsMagnitudeHistory[ k ] );
                    temp( 7 ) = fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ];
                    temp( 8 ) = 1000 * std::abs( ( fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ] - fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k - 1 ] ) / fitnessMagnitudeAverageOfAverageTopIndividualsMagnitudeHistory[ k ] );

                    fitnessMagnitudeAveragesOutputHistory[ k ] = temp;

                    //std::cout<< "  " << std::endl;
                    // std::cout<< "Identifying and Storing Top Individuals from this generation" << std::endl;
                    //topIndividualsMap = bislip::Variables::selectTopIndividuals( bodyMap, problemInput, topIndividuals );

                    if( k % populationSizePerIsland == 0 )
                    {
                        // std::cout << " " << std::endl;
                        // std::cout << "Clear print status" << std::endl;
                        for( unsigned int i = 0; i < suffixes.size(); i++ )
                        { problemInput->updateIndividualPrintStatus( suffixes[ i ], 0 ); }


                        //std::cout<< "  " << std::endl;
                        //std::cout<< "Identifying and Storing Top Individuals from this generation" << std::endl;
                        topIndividualsMap = bislip::Variables::selectTopIndividuals( bodyMap, problemInput, topIndividuals );

                        if( int( outputSettingsValues[ 0 ] ) == 1 )
                        {
                            std::cout<< "Printing of population, fitness, extremes, and fitness magnitude history" << std::endl;

                            bislip::utilities::writeDataMapToExistingTextFile( problemInput->getPopulation(),"populationHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                            bislip::utilities::writeDataMapToExistingTextFile( problemInput->getFitness(),"fitnessHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                            bislip::utilities::writeDataMapToExistingTextFile( problemInput->getExtremesAndConstraints(),"extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                            tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );

                            std::cout<< "Archiving population" << std::endl;
                            std::ofstream os( outputPath + outputSubFolder + "/archiveJSON_Population.dat", std::ios::binary );
                            { // save
                                cereal::PortableBinaryOutputArchive archiveSave( os );
                                archiveSave( population );
                            }

                        }
                    }

                    rerunCondition = ( k == static_cast< unsigned int >( evolutions / 4  + 0.5 ) ) || ( k == static_cast< unsigned int >( 2 * evolutions / 4 + 0.5 ) ) || ( k == static_cast< unsigned int >( 3 * evolutions / 4 + 0.5 ) );

                    if( rerunCondition )
                    {
                        std::cout<< "Trajectory data printing conditional met: " << "k = " << k << " == ";

                        if( k == static_cast< unsigned int >( evolutions / 4  + 0.5 ) )
                        { std::cout << static_cast< unsigned int >( evolutions / 4  + 0.5 ) << std::endl; }
                        else if( k == static_cast< unsigned int >( 2 * evolutions / 4  + 0.5 ) )
                        { std::cout << static_cast< unsigned int >( 2 * evolutions / 4  + 0.5 ) << std::endl; }
                        else if( k == static_cast< unsigned int >( 3 * evolutions / 4  + 0.5 ) )
                        { std::cout << static_cast< unsigned int >( 3 * evolutions / 4  + 0.5 ) << std::endl; }

                        bislip::Variables::reruns( bodyMap, problemInput, topIndividualsMap );

                    }


                    exitCondition1 = fitnessMagnitudeAveragesOutputHistory.at( k )( 2 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 2 ) ) == true )
                    { exitCondition2 = true; }
                    else
                    { exitCondition2 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 2 ) >= 1; }
                    exitCondition3 = fitnessMagnitudeAveragesOutputHistory.at( k )( 4 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 4 ) ) == true )
                    { exitCondition4 = true; }
                    else
                    { exitCondition4 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 4 ) >= 1; }
                    exitCondition5 = fitnessMagnitudeAveragesOutputHistory.at( k )( 6 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 6 ) ) == true )
                    { exitCondition6 = true; }
                    else
                    { exitCondition6 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 6 ) >= 1; }
                    exitCondition7 = fitnessMagnitudeAveragesOutputHistory.at( k )( 8 ) >= 1;
                    if( std::isnan( fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 8 ) ) == true )
                    { exitCondition8 = true; }
                    else
                    { exitCondition8 = fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 8 ) >= 1; }


                    k += 1;
                    suffixes.clear();
                    // std::cout<< "  " << std::endl;
                    // std::cout<< "  " << std::endl;
                } while( ( k < static_cast< unsigned int >( 3 * evolutions / 4 + 0.5 ) ) || ( ( k < evolutions + 1 ) && ( ( exitCondition1 && exitCondition2 ) || ( exitCondition3 && exitCondition4 ) || ( exitCondition5 && exitCondition6 ) || ( exitCondition7 && exitCondition8 ) ) ) );
            }

            std::cout<< "exitCondition1 = " << exitCondition1 << std::endl;
            std::cout<< "exitCondition2 = " << exitCondition2 << std::endl;
            std::cout<< "exitCondition3 = " << exitCondition3 << std::endl;
            std::cout<< "exitCondition4 = " << exitCondition4 << std::endl;
            std::cout<< "exitCondition5 = " << exitCondition5 << std::endl;
            std::cout<< "exitCondition6 = " << exitCondition6 << std::endl;
            std::cout<< "exitCondition7 = " << exitCondition7 << std::endl;
            std::cout<< "exitCondition8 = " << exitCondition8 << std::endl;
            std::cout<< "k - 2 -> " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 0 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 1 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 2 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 3 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 4 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 5 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 6 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 7 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 3 )( 8 ) << std::endl;
            std::cout<< "k - 1 -> " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 0 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 1 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 2 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 3 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 4 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 5 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 6 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 7 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 2 )( 8 ) << std::endl;
            std::cout<< "k     -> " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 0 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 1 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 2 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 3 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 4 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 5 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 6 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 7 ) << "  " << fitnessMagnitudeAveragesOutputHistory.at( k - 1 )( 8 ) << std::endl;

            if( singleIsland && optimizerIndex == 2 )
            {

                std::vector< std::string > suffixes;
                for( auto const& element : problemInput->getPopulation() ) { suffixes.push_back( element.first ); }

                // std::cout << " " << std::endl;
                // std::cout << "Clear print status" << std::endl;
                for( unsigned int i = 0; i < suffixes.size(); i++ )
                { problemInput->updateIndividualPrintStatus( suffixes[ i ], 0 ); }

                //std::cout<< "  " << std::endl;
                //std::cout<< "Identifying and Storing Top Individuals from this generation" << std::endl;
                topIndividualsMap = bislip::Variables::selectTopIndividuals( bodyMap, problemInput, topIndividuals );

            }

            bislip::Variables::reruns( bodyMap, problemInput, topIndividualsMap );
            bislip::utilities::writeDataMapToExistingTextFile( problemInput->getPopulation(),"populationHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            bislip::utilities::writeDataMapToExistingTextFile( problemInput->getFitness(),"fitnessHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            bislip::utilities::writeDataMapToExistingTextFile( problemInput->getExtremesAndConstraints(),"extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
            tudat::input_output::writeDataMapToTextFile( fitnessMagnitudeAveragesOutputHistory, "fitnessMagnitudeAveragesOutputHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );


            /*if( int( outputSettingsValues[ 0 ] ) == 1 )
            {
                std::cout<< "Printing of population, fitness, extremes, and fitness magnitude history with printed flag" << std::endl;
                tudat::input_output::writeDataMapToTextFile( problemInput->getPopulationHistory(), "populationHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                tudat::input_output::writeDataMapToTextFile( problemInput->getFitnessHistory(), "fitnessHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                tudat::input_output::writeDataMapToTextFile( problemInput->getExtremesAndConstraintsHistory(), "extremesAndConstraintsHistory.dat", outputPath + outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
                std::cout<< " " << std::endl;
            }*/

        }
    }
}


}

// namespace bislip

