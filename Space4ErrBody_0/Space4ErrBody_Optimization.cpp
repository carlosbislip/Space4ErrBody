
#include <pagmo/problem.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>

#include "Space4ErrBody_Executables_and_Headers/Space4ErrBody.h"
#include "Space4ErrBody_Executables_and_Headers/getAlgorithm.h"
#include "Space4ErrBody_Executables_and_Headers/applicationOutput_tudat.h"
#include "Space4ErrBody_Executables_and_Headers/applicationOutput_pagmo.h"
#include "Space4ErrBody_Executables_and_Headers/saveOptimizationResults.h"
#include "Space4ErrBody_Executables_and_Headers/getStuff.h"

//#include "Space4ErrBody_Executables_and_Headers_testing/getConditions.h"
//#include "Space4ErrBody_Executables_and_Headers_testing/getOptimizationSettings.h"
#include <boost/format.hpp>

#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

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
using namespace tudat_applications;
using namespace pagmo;
using namespace bislip;
using namespace boost;
using namespace tudat_pagmo_applications;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Main Program                  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    //! Get this run's time stamp. Will be used to create a unique output
    //! subfolder name, where all files created by the optimizer will be stored.
    std::string play_time = getCurrentDateTime( false );

    //! Get all inputs.
    const std::vector< std::string > primary                   = getDataString ( "primary.txt" );
    const std::vector< std::string > parameterList             = getDataString ( primary[ 0 ] );
    const std::vector< double > parameterBounds                = getDataNumeri ( primary[ 1 ] );
    const std::vector< std::string > vehicleParameterList      = getDataString ( primary[ 2 ] );
    const std::vector< double > vehicleParameterValues         = getDataNumeri ( primary[ 3 ] );
    const std::vector< std::string > aeroCoeffFileList         = getDataString ( primary[ 4 ] );
    const std::vector< std::string > optimization_settingsList = getDataString ( primary[ 5 ] );
    const std::vector< double > optimization_settingsValues    = getDataNumeri ( primary[ 6 ] );
    const std::vector< std::string > simulation_settingsList   = getDataString ( primary[ 7 ] );
    const std::vector< double > simulation_settingsValues      = getDataNumeri ( primary[ 8 ] );
    const std::vector< std::string > output_settingsList       = getDataString ( primary[ 9 ] );
    const std::vector< double > output_settingsValues          = getDataNumeri ( primary[ 10 ] );
    const std::vector< std::string > initialConditionsList     = getDataString ( primary[ 11 ] );
    const std::vector< double > initialConditionsValues        = getDataNumeri ( primary[ 12 ] );
    const std::vector< std::string > terminationConditionsList = getDataString ( primary[ 13 ] );
    const std::vector< double > terminationConditionsValues    = getDataNumeri ( primary[ 14 ] );
    const std::string problem_name                             = primary.rbegin()[ 1 ];
    const std::string vehicle_name                             = primary.rbegin()[ 0 ];

    //! Set seed for reproducible results.
    pagmo::random_device::set_seed( optimization_settingsValues[ 4 ] );

    //! Create partial output subfolder filename suffix based on optimization
    //! settings and fixed time step.
    const std::string this_run_settings = std::to_string( int( optimization_settingsValues[ 0 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 1 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 2 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 3 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 4 ] ) ) + "_" +
            std::to_string( simulation_settingsValues[ 3 ] );

    //! Create output subfolder filename. Based on arbitrary prefix, the
    //! previously created partial suffix, and the run's time stamp.
    const std::string outputSubFolder = "OUTPUT_" + this_run_settings + "_" + play_time + "/";

    //! Determine number of parameters to vary based of size of parameterBounds vector.
    const int nodes = simulation_settingsValues[ 3 ];
    const int N = parameterList.size() * nodes - 1;

    //! Create vector containing decision vector bounds.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( N, 0.0 ) );
    //!---------------------------------------   ^ for lb/up  (rows)       ^ for # of parameters

    //! Loop to build the bounds matrix.
    int p;

    for( int i = 0; i < N; i++ )
    {
        if( i <  ( nodes -1 ) )
        {
            p = 0;
        }
        else if( ( i - ( nodes - 1 ) ) % nodes == 0 )
        {
            p += 2;
        }

        bounds[ 0 ][ i ] = parameterBounds[ p ];
        bounds[ 1 ][ i ] = parameterBounds[ p + 1 ];
    }

    //! Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    //! Create object to compute the problem fitness; no perturbations
    //! An example uses 'extended dynamics'. Some conditional in the original
    //! example to include gravitational pertubations from additional celestial
    //! bodies, such as the Moon, Sun, etc. Will have to change some things in
    //! "space4Errbody.h" If I want to use it, maybe the "extended dynamics'
    //! cases would then include aerodynamics and eventually thrust. Probably
    //! not. might be a bit too messy.
    problem prob{ Space4ErrBody( bounds,
                                 problem_name,
                                 vehicle_name,
                                 parameterList,
                                 vehicleParameterValues,
                                 aeroCoeffFileList,
                                 simulation_settingsValues,
                                 initialConditionsValues,
                                 terminationConditionsValues,
                                 output_settingsValues,
                                 outputSubFolder ) };

    //! Retrieve algorithm. Three options available in the following function:
    //!        getMultiObjectiveAlgorithm
    //!               case 0 --> nsga2
    //!               case 1 --> moead
    //!               case 2 --> ihs
    //! Selection is currently arbitrary. moead.hpp has been modified such that
    //! the points are generated with a Halton sequence.
    pagmo::algorithm algo{ getMultiObjectiveAlgorithm( optimization_settingsValues[ 0 ] ) };

    //! Assign population size.
    pagmo::population::size_type populationSize = optimization_settingsValues[ 1 ];

    //! Assign archipelago size.
    pagmo::archipelago::size_type archipelagoSize = optimization_settingsValues[ 2 ];

    //! Assign population per archipelago
    const int pop_per_archi = populationSize/archipelagoSize;

    //! Assign number of evolutions
    const int evolutions = optimization_settingsValues[ 3 ];

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    const int index = int( optimization_settingsValues[ 0 ] );
    std::string algo_method;
    switch( index )
    {
    case 0:
    {
        algo_method = "NSGA2";
        break;
    }
    case 1:
    {
        algo_method = "MOEAD";
        break;
    }
    case 2:
    {
        algo_method = "IHS";
        break;
    }
    }


    //! Convert angles from degrees to radians
    const double lat_i_rad = unit_conversions::convertDegreesToRadians( initialConditionsValues[ 0 ] );
    const double lon_i_rad = unit_conversions::convertDegreesToRadians( initialConditionsValues[ 1 ] );
    const double lat_f_rad = unit_conversions::convertDegreesToRadians( terminationConditionsValues[ 0 ] );
    const double lon_f_rad = unit_conversions::convertDegreesToRadians( terminationConditionsValues[ 1 ] );

    //! Calculate initial heading angle: https://www.movable-type.co.uk/scripts/latlong.html
    //const double chi_i_rad_calc = getHeadingToTarget( lat_i_rad , lon_i_rad , lat_f_rad , lon_f_rad );
    //std::atan2( std::sin( lon_f_rad - lon_i_rad ) * std::cos( lat_f_rad ) , std::cos( lat_i_rad ) * std::sin( lat_f_rad ) - std::sin( lat_i_rad ) * std::cos( lat_f_rad ) * std::cos( lon_f_rad - lon_i_rad ) );
    double chi_i_deg_calc = unit_conversions::convertRadiansToDegrees( getHeadingToTarget( lat_i_rad , lon_i_rad , lat_f_rad , lon_f_rad ) );
    //! If heading angle is negative, this may help visualize it.
    if (chi_i_deg_calc < 0)
    {
        chi_i_deg_calc = 360 + chi_i_deg_calc;
    }

    //! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
    const double a = std::sin( (lat_f_rad - lat_i_rad) / 2) * std::sin( (lat_f_rad - lat_i_rad) / 2) + std::cos( lat_i_rad ) * std::cos( lon_i_rad ) * std::sin( (lon_f_rad - lon_i_rad) / 2) * std::sin( (lon_f_rad - lon_i_rad) / 2);
    const double d_angular = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
    const double d_spherical_law_cosines = unit_conversions::convertRadiansToDegrees( getAngularDistance( lat_i_rad , lon_i_rad , lat_f_rad , lon_f_rad ) );
    //std::acos( std::sin(lat_i_rad) * std::sin(lat_f_rad) + std::cos(lat_i_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_i_rad) ) * spice_interface::getAverageRadius( "Earth" );

    //! Print to Terminal screen.
    std::cout << " " << std::endl;
    std::cout << problem_name << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Parameters to vary" << std::endl;
    for( int i = 0; i < int( parameterList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.1f %-2s %-8.1f %-2s\n") % parameterList[i] % "[" % bounds[0][i] % "," % bounds[1][i] % "]";
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
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Aerodynamic Coefficients" << std::endl;
    std::cout << "Coefficient File: " << "      " <<  aeroCoeffFileList[0] << std::endl;
    std::cout << "Coefficient File: " << "      " <<  aeroCoeffFileList[1] << std::endl;
    std::cout << "Coefficient File: " << "      " <<  aeroCoeffFileList[2] << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Optimizations Settings" << std::endl;
    for( int i = 0; i < int( optimization_settingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % optimization_settingsList[i] % "" % optimization_settingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Simulation Settings" << std::endl;
    for( int i = 0; i < int( simulation_settingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % simulation_settingsList[i] % "" % simulation_settingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output Settings" << std::endl;
    for( int i = 0; i < int( output_settingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % output_settingsList[i] % "" % output_settingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Initial Conditions" << std::endl;
    for( int i = 0; i < int( initialConditionsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % initialConditionsList[i] % "" % initialConditionsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Termination Conditions" << std::endl;
    for( int i = 0; i < int( terminationConditionsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % terminationConditionsList[i] % "" % terminationConditionsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Ground distance to cover " << std::endl;
    std::cout << "  Haversine Formula:        " << unit_conversions::convertRadiansToDegrees( d_angular ) << " degrees." << std::endl;
    std::cout << "  Spherical Law of Cosines: " << d_spherical_law_cosines << " degrees." << std::endl;
    std::cout << "Heading angle to Target:    " << chi_i_deg_calc << " degrees. Calculated."  << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output subfolder: " << std::endl;
    std::cout << "     '" << outputSubFolder <<"'" << std::endl;
    std::cout << " DV  " << " Lower Boundary " << " Upper Boundary " << std::endl;
    for( int i = 0; i < N; i++ )
    {
        std::cout << std::fixed << std::setprecision(10) <<
                     std::setw(6) << i <<
                     std::setw(17) << bounds[ 0 ][ i ] <<
                     std::setw(17) << bounds[ 1 ][ i ] << std::endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    //! Create an island with populationSize individuals
    pagmo::island isl{algo, prob, populationSize};

    //! Instantiate an archipelago
    /*  archipelago archi{archipelagoSize, algo, prob, pop_per_archi};

    //! Evolve
    archi.evolve(evolutions);

    //! Wait for the evolutions to be finished
    archi.wait_check();// Raises errors,

    for (const auto &isl : archi) {
        std::cout << isl.get_population().champion_f()[0] << '\n';
    }
    for (const auto &isl : archi) {
        std::cout << isl.get_population().champion_x()[0] << '\n';
    }
    for (const auto &isl : archi) {
     printPopulationToFile( isl.get_population().champion_f()[0],
             "mo_AMSIAD_BALLISTIC_" ,
             false );
     printPopulationToFile( isl.get_population().champion_x()[0],
             "mo_AMSIAD_BALLISTIC_" ,
             true );
    }
*/
    //! Write original (unevolved) population to file
    if  ( output_settingsValues[ 0 ] == 1 )
    {
     //   printPopulationToFile( isl_pert.get_population( ).get_x( ), this_run_settings_1,
     //                          outputSubFolder, false );
     //   printPopulationToFile( isl_pert.get_population( ).get_f( ), this_run_settings_1,
     //                          outputSubFolder, true );
    }

    for( int i = 0; i < evolutions; i++ )
    {
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        if  ( output_settingsValues[ 0 ] == 1 )
        {

            std::string this_run_settings_1 = std::to_string( i ) + "_" +
                    this_run_settings;

            //! Write current iteration results to file
            printPopulationToFile( isl.get_population( ).get_x( ),
                                   this_run_settings_1,
                                   outputSubFolder,
                                   false );
            printPopulationToFile( isl.get_population( ).get_f( ),
                                   this_run_settings_1,
                                   outputSubFolder,
                                   true ); //! true boolean is for the printPopulationToFile to print with a different prefix: fitness_
        }

        std::cout<<i<<std::endl;

    }

}

