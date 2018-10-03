
#include <pagmo/problem.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>

#include "Space4ErrBody_Executables_and_Headers_testing/Space4ErrBody.h"
#include "Space4ErrBody_Executables_and_Headers_testing/getAlgorithm.h"
#include "Space4ErrBody_Executables_and_Headers_testing/applicationOutput_tudat.h"
#include "Space4ErrBody_Executables_and_Headers_testing/applicationOutput_pagmo.h"
#include "Space4ErrBody_Executables_and_Headers_testing/saveOptimizationResults.h"
#include "Space4ErrBody_Executables_and_Headers_testing/getStuff.h"

//#include "Space4ErrBody_Executables_and_Headers_testing/getConditions.h"
//#include "Space4ErrBody_Executables_and_Headers_testing/getOptimizationSettings.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"
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
using namespace tudat_pagmo_applications;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Main Program                  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    //! Get this run's time stamp. Will be used to create a unique output
    //! subfolder name, where all files created by the optimizer will be stored.
    std::string play_time = getCurrentDateTime( false );

    //! Get list of files containing all relevant information for this run.
    //! Size of vector is NOT pre-defined.
    std::vector< std::string > file_list = getFiles ( "file_list.txt" );


    //! Get conditions that are passed to the propagator. Size of vector is NOT
    //! pre-defined. The structure defined here is that for actual run, not for
    //! validation cases.
    //!     Validation flag
    //!     initial Epoch
    //!     max additional time from initial epoch
    //!     fixed time step
    //!     reference area
    //!     aerodynamic coefficient (single value)
    //!     mass of vehicle (initial?)
    //!     initial height (relative to Earth's surface)
    //!     initial latitude (Earth-Fixed)
    //!     initial longitude (Earth-Fixed)
    //!     final height (relative to Earth's surface)
    //!     final latitude (Earth-Fixed)
    //!     final longitude (Earth-Fixed)
    std::vector< double > input_data = getData ( file_list[0] );

    //! Get optimization settings. Size of vector is NOT pre-defined.
    //!     population size
    //!     archipelago size
    //!     number of evolutions
    std::vector< double > opt_set = getData ( file_list[1] );

    //! Get decision variable bounds. Size of vector is NOT pre-defined.
    //!     lb - initial velocity (Earth-Fixed)
    //!     ub - initial velocity (Earth-Fixed)
    //!     lb - initial flight-path angle
    //!     ub - initial flight-path angle
    //!     lb - initial heading angle
    //!     ub - initial heading angle
    //!     lb - angle of attack
    //!     ub - angle of attack
    std::vector< double > dv_bounds = getData ( file_list[2] );

    //! Create partial output subfolder filename suffix based on optimization
    //! settings and fixed time step.
    std::string this_run_settings = std::to_string(int(opt_set[0])) + "_" +
            std::to_string(int(opt_set[1])) + "_" +
            std::to_string(int(opt_set[2])) + "_" +
            std::to_string(int(opt_set[3])) + "_" +
            std::to_string(input_data[3]);

    //! Get output settings. Size of vector is NOT pre-defined.
    //!     print fitness & population
    //!     print propagation history
    //!     print dependent varialbes history
    std::vector< double > output_settings = getData ( file_list[3] );


    //! Create output subfolder filename. Based on arbitrary prefix, the
    //! previously created partial suffix, and the run's time stamp.
    std::string outputSubFolder = "HORUS_OUTPUT_" + this_run_settings + "_" + play_time + "/";
    if ( int(input_data[0]) == 0 )
    {
        outputSubFolder = "HORUS_OUTPUT_" + this_run_settings + "_" + play_time + "/";
    }
    else
    {
        outputSubFolder = "HORUS_VALIDATION_" + this_run_settings + "_" + play_time + "/";
    }

    //! Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    //! Set seed for reproducible results. Should I make this an input value?
    //! It may be better to hardcode this value to avoid issues down the line.
    pagmo::random_device::set_seed(234);

    //! Determine number of parameters to vary based of size of dv_bounds input.
    const int N = dv_bounds.size()/2;

    //! Create vector containing decision vector bounds.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( N, 0.0 ) );
    //!---------------------------------------   ^ for lb/up  (rows)       ^ for # of parameters

    if ( int(input_data[0]) == 0 )
    {
        //! Assign bounds of decision variables. Loop unpacks data decision vector
        //! bounds into the bounds matrix that will be passed on for the optimizer.
        //!       Initial velocity (Earth-Fixed) - v_i
        //!       Initial flight-path angle      - gamma_i
        //!       Initial heading angle          - chi_i
        //!       Angle of attack                - AoA
        int p = 0;
        for( int i = 0; i < N; i++ )
        {
            bounds[ 0 ][ i ] = dv_bounds[ p ];
            bounds[ 1 ][ i ] = dv_bounds[ p + 1 ];
            p = p + 2;
        }
    }
    else
    {
        //! Same loop as before. However, the purpose of this distinction it
        //! that a validation exercise is to be done. Various initial conditions
        //! are known. The corresponding bounds of interest are then tightened
        //! to reflect this. The bounds are still defined in an input file.
        //!       Initial velocity (Earth-Fixed) - v_i = 7435.5 m/s
        //!                                              Not sure if the stated
        //!                                              velocity is Earth-Fixed
        //!                                              or Inertial.
        //!       Initial flight-path angle      - gamma_i = -1.437 deg
        //!       Initial heading angle          - chi_i = 70.757 deg
        int p = 0;
        for( int i = 0; i < N; i++ )
        {
            bounds[ 0 ][ i ] = dv_bounds[ p ];
            bounds[ 1 ][ i ] = dv_bounds[ p + 1 ];
            p = p + 2;
        }

    }

    //! Create object to compute the problem fitness; no perturbations
    //! An example uses 'extended dynamics'. Some conditional in the original
    //! example to include gravitational pertubations from additional celestial
    //! bodies, such as the Moon, Sun, etc. Will have to change some things in
    //! "space4Errbody.h" If I want to use it, maybe the "extended dynamics'
    //! cases would then include aerodynamics and eventually thrust. Probably
    //! not. might be a bit too messy.
    problem prob{Space4ErrBody( bounds, input_data, output_settings, outputSubFolder) };

    //! Retrieve algorithm. Three options available in the following function:
    //!        getMultiObjectiveAlgorithm
    //!               case 0 --> nsga2
    //!               case 1 --> moead
    //!               case 2 --> ihs
    //! Selection is currently arbitrary. moead.hpp has been modified such that
    //! the points are generated with a Halton sequence.
    algorithm algo{getMultiObjectiveAlgorithm( opt_set[0] )};

    //! Assign population size.
    pagmo::population::size_type populationSize = opt_set[1];

    //! Assign archipelago size.
    pagmo::archipelago::size_type archipelagoSize = opt_set[2];

    //! Assign population per archipelago
    const int pop_per_archi = populationSize/archipelagoSize;

    //! Assign number of evolutions
    const int evolutions = opt_set[3];

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    //! If we are here, we have read, assigned, created, and passed mostly all
    //! information. The next step is to actually start doing stuff. However,
    //! before that, let's print out some stuff to the screen to make sure we
    //! know what's supposed to happen.

    //! Assign Reference area
    const double Ref_area = input_data[4]; //m^2

    //! Assign initial mass
    const double M_i = input_data[6]; // kg

    //! Declare various variables to unpack from input_data depending on a
    //! conditional.
    double v_i, gamma_i_deg, chi_i_deg;
    double h_i, lat_i_deg,lon_i_deg;
    double h_f, lat_f_deg,lon_f_deg;

    if ( int(input_data[0]) == 0 )
     {
         //! Assign the initial and final position conditions.
         h_i         = input_data[7]; //0 * 1E3; // m
         lat_i_deg   = input_data[8]; //52.30805556; //52deg 18’29"N
         lon_i_deg   = input_data[9]; //4.76416667 //4deg 45’51"E
         h_f         = input_data[10]; //0 * 1E3; // m
         lat_f_deg   = input_data[11]; //38.9444444444444; //38deg 56’40"N
         lon_f_deg   = input_data[12]; //282.544166666667 //77deg 27’21"W
     }
     else
     {
         //! Assign expected initial values. Only used to display and remind
         //! what the goal values are supposed to be.
         v_i         = input_data[7]; //7435.5 m/s // given
         gamma_i_deg = input_data[8]; //-1.437 * 1E3 kg // given
         chi_i_deg   = input_data[9]; //70.757 deg // given

         //! Assign the initial and final position conditions.
         h_i         = input_data[10]; //122 * 1E3 m // given
         lat_i_deg   = input_data[11]; //-22.37 deg // given
         lon_i_deg   = input_data[12]; //-106.7 deg // given
         h_f         = input_data[13]; //25 * 1E3 m // arbitrary?
         lat_f_deg   = input_data[14]; //5.237222 deg //5 deg 14’14"N // 5.0 deg according to Mooij Dissertation
         lon_f_deg   = input_data[15]; //-52.760556 deg //52 deg 45’38"W// -53.0 deg according to Mooij Dissertation
     }

    //! Convert angles from degrees to radians
    const double lat_i_rad = unit_conversions::convertDegreesToRadians( lat_i_deg );
    const double lon_i_rad = unit_conversions::convertDegreesToRadians( lon_i_deg );
    const double lat_f_rad = unit_conversions::convertDegreesToRadians( lat_f_deg );
    const double lon_f_rad = unit_conversions::convertDegreesToRadians( lon_f_deg );

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
    const double d_haversine = d_angular * spice_interface::getAverageRadius( "Earth" );
    const double d_spherical_law_cosines = unit_conversions::convertRadiansToDegrees( getAngularDistance( lat_i_rad , lon_i_rad , lat_f_rad , lon_f_rad ) );
    //std::acos( std::sin(lat_i_rad) * std::sin(lat_f_rad) + std::cos(lat_i_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_i_rad) ) * spice_interface::getAverageRadius( "Earth" );

    const int index = int(opt_set[0]);
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

    //! Print to Terminal screen.
       std::cout << " " << std::endl;
       if ( int(input_data[0]) == 0 )
       {
           std::cout << "AMS to IAD Ballistic Trajectory: 10/28/2018  11:00:00 AM" << std::endl;
       }
       else
       {
           std::cout << "Validation: HORUS Entry towards Kourou" << std::endl;
       }
       std::cout << "--------------------------------------------------------" << std::endl;
       std::cout << "Starting Epoch:             " << input_data[1] << " seconds." << std::endl;
       std::cout << "Max. Epoch dif.:            " << input_data[2] << " seconds." << std::endl;
       if ( int(input_data[0]) == 0 )
       {
           std::cout << "Initial mass:               " << M_i << " kg." << std::endl;
           std::cout << "Reference area:             " << Ref_area << " m^2." << std::endl;
       }
       else
       {
           std::cout << "Initial mass:               " << M_i << " kg. Dont have it yet" << std::endl;
           std::cout << "Reference area:             " << Ref_area << " m^2. Dont have it yet" << std::endl;
       }
       std::cout << "Initial height:             " << h_i / 1E3 << " km above Earth's surface." << std::endl;
       std::cout << "Final height:               " << h_f / 1E3 << " km above Earth's surface." << std::endl;
       std::cout << "Initial latitude:           " << lat_i_deg << " degrees. Earth-Fixed" << std::endl;
       std::cout << "Initial longitude:          " << lon_i_deg << " degrees. Earth-Fixed" << std::endl;

       if ( int(input_data[0]) == 0 )
       {
           std::cout << "GOAL latitude:              " << lat_f_deg << " degrees. Earth-Fixed" << std::endl;
           std::cout << "GOAL longitude:             " << lon_f_deg << " degrees. Earth-Fixed" << std::endl;
       }
       else
       {
           std::cout << "GOAL latitude:              " << lat_f_deg << " degrees. Not sure if its this one." << std::endl;
           std::cout << "GOAL longitude:             " << lon_f_deg << " degrees. Not sure if its this one." << std::endl;
           std::cout << "Initial velocity:           " << v_i << " m/s. Unclear if Earth-Fixed." << std::endl;
           std::cout << "Initial flight-path angle:  " << gamma_i_deg << " degrees." << std::endl;
           std::cout << "Initial heading:            " << chi_i_deg << " degrees." << std::endl;
       }
       std::cout << "Ground distance to cover " << std::endl;
       std::cout << "  Haversine Formula:        " << d_haversine / 1E3 << " km." << std::endl;
       std::cout << "  Spherical Law of Cosines: " << d_spherical_law_cosines << " degrees." << std::endl;
       std::cout << "Initial Heading angle:      " << chi_i_deg_calc << " degrees. Calculated."  << std::endl;
       std::cout << "Optimization method:        " << algo_method << std::endl;
       std::cout << "Population size:            " << populationSize << std::endl;
       std::cout << "Num. of archipelagos:       " << archipelagoSize << std::endl;
       std::cout << "Num. of evolutions:         " << evolutions << std::endl;
       std::cout << "Output subfolder: " << std::endl;
       std::cout << "     '" << outputSubFolder <<"'" << std::endl;
       std::cout << " DV " << " Lower Boundary " << " Upper Boundary " << std::endl;
    for( int i = 0; i < N; i++ )
    {
        std::cout << std::fixed << std::setprecision(10) <<
                     std::setw(4) << i + 1 <<
                     std::setw(16) << bounds[ 0 ][ i ] <<
                     std::setw(16) << bounds[ 1 ][ i ] << std::endl;
    }


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    //! Create an island with populationSize individuals
    island isl{algo, prob, populationSize};

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




    for( int i = 0; i < evolutions; i++ )
    {
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        if (output_settings[0] == 1)
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
/*
    // Create object to compute the problem fitness; with perturbations
    problem prob_pert{Space4ErrBody_Ballistic( var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, true ) };


    // Instantiate a pagmo algorithm for the new problem.
    algorithm algo_pert{de1220( )};

    // Create an empty population for perturbed problem
    population population_pert = population( prob_pert, 0 );

    // Retrieve population of unperturbed problem, and instantiate population of perturbed problem
    std::vector<vector_double> original_population = isl.get_population( ).get_x( );
    for( unsigned int k = 0; k < populationSize; k++ )
    {
        population_pert.push_back( original_population.at( k ) );
    }

    // Create island for perturbed problem
    island isl_pert{algo_pert, population_pert};

    // Perform Grid Search for perturbed priblem and write results to file
    if( performGridSearch )
    {
        createGridSearch( isl_pert.get_population( ).get_problem( ), {{0.0, 0.0}, {360.0, 180.0}}, { 100, 50 }, "propagationTargetingGridSearch_pert" );
    }

    // Write original (unevolved) population to file
    printPopulationToFile( isl_pert.get_population( ).get_x( ), "targetingPropagation_pert_orig" , false );
    printPopulationToFile( isl_pert.get_population( ).get_f( ), "targetingPropagation_pert_orig" , true );


    // Evolve for 4 generations
    for( int i = 0; i < 1; i++ )
    {
        isl_pert.evolve( );
        while( isl_pert.status( ) != pagmo::evolve_status::idle &&
               isl_pert.status( ) != pagmo::evolve_status::idle_error )
        {
            isl_pert.wait( );
        }
        isl_pert.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( isl_pert.get_population( ).get_x( ), "targetingPropagation_pert_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
        printPopulationToFile( isl_pert.get_population( ).get_f( ), "targetingPropagation_pert_" + std::to_string( i ) + "_" + std::to_string( i ) , true );

        std::cout<<i<<std::endl;
    }

    */
}

