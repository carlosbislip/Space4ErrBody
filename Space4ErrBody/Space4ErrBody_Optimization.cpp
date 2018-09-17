
#include <pagmo/problem.hpp>
#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include <pagmo/algorithms/de.hpp>
#include <pagmo/algorithms/simulated_annealing.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>

#include "Space4ErrBody_Executables_and_Headers/Space4ErrBody.h"
#include "Space4ErrBody_Executables_and_Headers/getAlgorithm.h"
#include "Space4ErrBody_Executables_and_Headers/applicationOutput_tudat.h"
#include "Space4ErrBody_Executables_and_Headers/applicationOutput_pagmo.h"
#include "Space4ErrBody_Executables_and_Headers/saveOptimizationResults.h"
#include "Space4ErrBody_Executables_and_Headers/getCurrentDateTime.h"

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
//using namespace pagmo::problem;
using namespace tudat_pagmo_applications;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Read Input File               //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void readinputfile ( double *var1, double *var2, double *var3, double *var4,
                     double *var5, double *var6, double *var7, double *var8,
                     double *var9, double *var10, double *var11, double *var12,
                     double *var13, double *var14, double *var15, double *var16,
                     double *var17, double *var18, double *var19, double *var20,
                     double *var21)
//void readinputfile ( std::vector< double > *input_data )
{
    using namespace std;

 const std::string& filename = "input.txt";
    std::ifstream inputdata;
    inputdata.open(filename.c_str());
    if(inputdata.fail())
    {
        std::cout << "File '" << filename << "' failed to open";
        std::cin.ignore(100, '\n');
        std::exit(1);
//        return 1;
    }
    else{
        std::cout << "The following text had been read." << std::endl;
        std::cout << " " << std::endl;
        std::cout << inputdata.rdbuf();  // Prints content to screen.
        inputdata.seekg(std::ios_base::beg) ;
        std::cout << " " << std::endl;
        std::cout << "The program will now proceed to import the data." << std::endl;


        inputdata>>*var1;
        inputdata>>*var2;
        inputdata>>*var3;
        inputdata>>*var4;
        inputdata>>*var5;
        inputdata>>*var6;
        inputdata>>*var7;
        inputdata>>*var8;
        inputdata>>*var9;
        inputdata>>*var10;
        inputdata>>*var11;
        inputdata>>*var12;
        inputdata>>*var13;
        inputdata>>*var14;
        inputdata>>*var15;
        inputdata>>*var16;
        inputdata>>*var17;
        inputdata>>*var18;
        inputdata>>*var19;
        inputdata>>*var20;
        inputdata>>*var21;
        inputdata.close();
        /*std::cout << "Imported to local variables." << std::endl;
        std::cout << "Variable 1 is: " << var1 << std::endl;
        *input_data.push_back(var1);
        *input_data.push_back(var2);
        *input_data.push_back(var3);
        *input_data.push_back(var4);
        *input_data.push_back(var5);
        input_data.push_back(var6);
        input_data.push_back(var7);
        input_data.push_back(var8);
        input_data.push_back(var9);
        input_data.push_back(var10);
        input_data.push_back(var11);
        input_data.push_back(var12);
        input_data.push_back(var13);
        input_data.push_back(var14);
        std::cout << "Placed all variables" << std::endl;*/

        }
}






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Main Program                  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{

     std::string play = getCurrentDateTime( false );

    double var1;
    double var2;
    double var3;
    double var4;
    double var5;
    double var6;
    double var7;
    double var8;
    double var9;
    double var10;
    double var11;
    double var12;
    double var13;
    double var14;
    double var15;
    double var16;
    double var17;
    double var18;
    double var19;
    double var20;
    double var21;
    std::vector< double > input_data;

    //! Read input file. Provides some basic data that need not be hardcoded.
    readinputfile ( &var1, &var2, &var3, &var4, &var5, &var6, &var7,
                    &var8, &var9, &var10, &var11, &var12, &var13, &var14,
                    &var15, &var16, &var17, &var18, &var19, &var20, &var21);
    //readinputfile ( &input_data );

   //! Plop it all into a vector. Maybe push this to the readinputfile subroutine
    input_data.push_back(var1);
    input_data.push_back(var2);
    input_data.push_back(var3);
    input_data.push_back(var4);
    input_data.push_back(var5);
    input_data.push_back(var6);
    input_data.push_back(var7);
    input_data.push_back(var8);
    input_data.push_back(var9);
    input_data.push_back(var10);
    input_data.push_back(var11);
    input_data.push_back(var12);
    input_data.push_back(var13);
    input_data.push_back(var14);
    input_data.push_back(var15);
    input_data.push_back(var16);
    input_data.push_back(var17);
    input_data.push_back(var18);
    input_data.push_back(var19);
    input_data.push_back(var20);
    input_data.push_back(var21);

    //! Lazy hack to pass this info into the fitness function

  //  std::chrono::system_clock::now()

    //! Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    const double h_i       = input_data[6]; //0 * 1E3; // m
    const double lat_i_deg = input_data[7]; //52.30805556; //52deg 18’29"N
    const double lon_i_deg = input_data[8]; //4.76416667 //4deg 45’51"E    
    const double h_f       = input_data[9]; //0 * 1E3; // m
    const double lat_f_deg = input_data[10]; //38.9444444444444; //38deg 56’40"N
    const double lon_f_deg = input_data[11]; //282.544166666667 //77deg 27’21"W

    //! Convert angles from degrees to radians
    const double lat_i_rad = unit_conversions::convertDegreesToRadians( lat_i_deg );
    const double lon_i_rad = unit_conversions::convertDegreesToRadians( lon_i_deg );
    const double lat_f_rad = unit_conversions::convertDegreesToRadians( lat_f_deg );
    const double lon_f_rad = unit_conversions::convertDegreesToRadians( lon_f_deg );

    //! Calculate initial heading angle: https://www.movable-type.co.uk/scripts/latlong.html
    const double chi_i_rad = std::atan2( std::sin( lon_f_rad - lon_i_rad ) * std::cos( lat_f_rad ) , std::cos( lat_i_rad ) * std::sin( lat_f_rad ) - std::sin( lat_i_rad ) * std::cos( lat_f_rad ) * std::cos( lon_f_rad - lon_i_rad ) );
    double chi_i_deg = chi_i_rad * 180 / mathematical_constants::PI;

    if (chi_i_deg < 0)
    {
        chi_i_deg = 360 + chi_i_deg;
    }

    //! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
    const double a = std::sin( (lat_f_rad - lat_i_rad) / 2) * std::sin( (lat_f_rad - lat_i_rad) / 2) + std::cos( lat_i_rad ) * std::cos( lon_i_rad ) * std::sin( (lon_f_rad - lon_i_rad) / 2) * std::sin( (lon_f_rad - lon_i_rad) / 2);
    const double c = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
    const double d_haversine = c * spice_interface::getAverageRadius( "Earth" );
    const double d_spherical_law_cosines =  std::acos( std::sin(lat_i_rad) * std::sin(lat_f_rad) + std::cos(lat_i_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_i_rad) ) * spice_interface::getAverageRadius( "Earth" );

    //! Print out what is initiall being passed into the Spherical State vector
    std::cout << " " << std::endl;
    std::cout << "AMS to IAD Ballistic Trajectory: 10/28/2018  11:00:00 AM" << std::endl;
    std::cout << "Imported values that have been read and passed---------------" << std::endl;
    std::cout << "Starting Epoch: " << input_data[0] << " seconds." << std::endl;
    std::cout << "Initial height: " << h_i / 1E3 << " km above the Earth's surface." << std::endl;
    std::cout << "Final height: " << h_f / 1E3 << " km above the Earth's surface." << std::endl;
    std::cout << "Origin latitude: " << lat_i_deg << " degrees." << std::endl;
    std::cout << "Origin longitude: " << lon_i_deg << " degrees." << std::endl;
    std::cout << "Destination latitude: " << lat_f_deg << " degrees." << std::endl;
    std::cout << "Destination longitude: " << lon_f_deg << " degrees." << std::endl;
    std::cout << "Groud distance to cover " << std::endl;
    std::cout << "     Haversine Formula:        " << d_haversine / 1E3<< " km." << std::endl;
    std::cout << "     Spherical Law of Cosines: " << d_spherical_law_cosines / 1E3<< " km." << std::endl;
    std::cout << "Initial Heading angle " << chi_i_deg << " degrees."  << std::endl;

    //bool performGridSearch = false;

    //! Set seed for reproducible results. Should I make this an input value?
    pagmo::random_device::set_seed(234);

    //! For the ballistic case we have two decision variables each with a lower
    //! and upper bound. These are initial velocity and initial flight path
    //! angle. We now create a vector of vectors that will contain these.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 3, 0.0 ) );

    //! Define bounds (should I make these input values?)
    //!       Initial velocity will be varied from 1,000 to 9,000 m/s (inertial)
    //!       Initial flight-path angle with be varied form 20 to 80 deg.
    bounds[ 0 ][ 0 ] = input_data[12];
    bounds[ 1 ][ 0 ] = input_data[13];
    bounds[ 0 ][ 1 ] = input_data[14];
    bounds[ 1 ][ 1 ] = input_data[15];
    bounds[ 0 ][ 2 ] = input_data[16];
    bounds[ 1 ][ 2 ] = input_data[17];

    //! Create object to compute the problem fitness; no perturbations
    //! The original uses 'extended dynamics'. Some conditional in the original
    //! example to include gravitational pertubations from additional celestial
    //! bodies, such as the Moon, Sun, etc. Will have to change some things in
    //! "space4Errbody.h" If I want to use it, maybe the "extended dynamics'
    //! cases would then include aerodynamics and eventually thrust.
    //problem prob{Space4ErrBody_Ballistic( var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, false ) };
    problem prob{Space4ErrBody_Ballistic( bounds, input_data, play) };

    //! Perform Grid Search and write results to file
    //if( performGridSearch )
    //{
    //   createGridSearch( prob, {{0.0, 0.0}, {360.0, 180.0}}, { 100, 50 },
    //                      "propagationTargetingGridSearch_" );
    //}

    //! Retrieve algorithm. Three options available in the following function:
    //!        getMultiObjectiveAlgorithm
    //!               case 0 --> nsga2
    //!               case 1 --> moead
    //!               case 2 --> ihs
    //! Selection is currently arbitrary.
    algorithm algo{getMultiObjectiveAlgorithm( 0 )};

    //! Define population size.
    pagmo::population::size_type populationSize = input_data[18];

    //! Define archipelago size.
  //  pagmo::archipelago::size_type archipelagoSize = input_data[17];

    //! Define population per archipelago
  //  const int pop_per_archi = populationSize/archipelagoSize;

    //! Define number of evolutions
    const int evolutions = input_data[20];

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
    std::string outputSubFolder = "HORUS_OUTPUT_" + std::to_string(int(input_data[18])) + "_" + std::to_string(int(input_data[19])) + "_" + std::to_string(int(input_data[20])) +  "_" + std::to_string(input_data[2]) + "_" + play + "/";


    for( int i = 0; i < evolutions; i++ )
    {
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        //! Write current iteration results to file
        printPopulationToFile( isl.get_population( ).get_x( ),
                               std::to_string( i ) + "_" + std::to_string(int(input_data[18])) + "_" + std::to_string(int(input_data[19])) + "_" + std::to_string(int(input_data[20])) +  "_" + std::to_string(input_data[3]),
                               outputSubFolder,
                               false );
        printPopulationToFile( isl.get_population( ).get_f( ),
                               std::to_string( i ) + "_" + std::to_string(int(input_data[18])) + "_" + std::to_string(int(input_data[19])) + "_" + std::to_string(int(input_data[20])) +  "_" + std::to_string(input_data[3]),
                               outputSubFolder,
                               true ); //! true boolean is for the printPopulationToFile to print with a different prefix: fitness_

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

