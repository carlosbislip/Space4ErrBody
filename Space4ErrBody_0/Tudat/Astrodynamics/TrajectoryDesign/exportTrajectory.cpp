#include <cmath>
#include <fstream>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"

#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Return a vector of positions and times corresponding to a trajectory at a certain epoch.
void returnTrajectory( Eigen::VectorXd initialCartesianState,
                       double centralBodyGravitationalParameter,
                       double duration,
                       double maximumTimeStep,
                       std::vector < Eigen::Vector3d >& positionVector,
                       std::vector < double >& timeVector,
                       double startingTime )
{
    // Calculate the number of intermediate steps required and the corresponding time step.
    int numberOfSteps = std::ceil( duration / maximumTimeStep );
    double timeStep = duration / numberOfSteps;

    // Resize the position and time vectors. Note that the "+ 1" comes from the end points.
    positionVector.resize( numberOfSteps + 1 );
    timeVector.resize( numberOfSteps + 1 );

    // Initialize Cartesian and Keplerian coordinates vectors, for intermediate calculations.
    Eigen::Vector6d cartesianElements ( 6 ), keplerianElements ( 6 );
    cartesianElements = initialCartesianState;


    // Initiate the time variable.
    double time = 0.0 + startingTime;

    // Set the first position and time to that of the departure.
    positionVector[ 0 ] = cartesianElements.segment( 0, 3 );
    timeVector[ 0 ] = time;

    // Calculate and set the intermediate points
    for ( int counter = 1; counter <= numberOfSteps; counter++ )
    {
        // Convert the cartesian elements into keplerian elements.
        keplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                    cartesianElements, centralBodyGravitationalParameter );


        // Propagate the keplerian elements with one timestep.
        keplerianElements = orbital_element_conversions::propagateKeplerOrbit( keplerianElements,
                    timeStep, centralBodyGravitationalParameter );


        // Convert the keplerian elements back into Cartesian elements.
        cartesianElements = orbital_element_conversions::convertKeplerianToCartesianElements(
                    keplerianElements, centralBodyGravitationalParameter );


        // Update the time.
        time += timeStep;

        // Store the current time and position in their corresponding entry in their vector.
        positionVector[ counter ] = cartesianElements.segment( 0 , 3 );
        timeVector[ counter ] = time;
    }
}

//! Return a vector of positions and times corresponding to a 2D circular trajectory.
void returnCircularTrajectory( double maximumTimeStep,
                               std::vector < Eigen::Vector3d >& positionVector,
                               std::vector < double >& timeVector,
                               double semiMajorAxis,
                               double centralGravitationalParamater,
                               double startingTime )
{
    // Calculate the mean motion, orbital period and the resulting number of steps required.
    double meanMotion = std::sqrt( centralGravitationalParamater / std::pow( semiMajorAxis, 3 ) );
    double orbitalPeriod = 2 * mathematical_constants::PI / meanMotion;
    double numberOfSteps = std::ceil( orbitalPeriod / maximumTimeStep );
    double timeStep = orbitalPeriod / numberOfSteps;

    // Create time variable, for use in the coming loop.
    double time = 0;

    timeVector.resize( numberOfSteps + 1 );
    positionVector.resize( numberOfSteps + 1 );

    // Create a for loop in which the positions are set to the vectors.
    for ( int counter = 0; counter < numberOfSteps + 1; counter++ )
    {
        timeVector[ counter ] = time + startingTime;
        positionVector[ counter ] << semiMajorAxis * std::cos( time * meanMotion ),
                                     semiMajorAxis * std::sin( time * meanMotion ), 0.0;
        time += timeStep;
    }
}

//! Write a trajectory to a data file.
void writeTrajectoryToFile( std::vector < Eigen::Vector3d > positionVector,
                            std::vector < double > timeVector,
                            const std::string fileName)
{
    std::ofstream exportFile( fileName );
    for ( unsigned int counter = 0; counter < timeVector.size( ); counter++ )
    {
        exportFile << timeVector[ counter ] << ", " << positionVector[ counter ]( 0 ) << ", "
                   << positionVector[ counter ]( 1 ) << ", " << positionVector[ counter ]( 2 )
                   << std::endl;
    }
    exportFile.close( );
}

//! Write a trajectory to a data file.
void writeManeuversToFile( std::vector < Eigen::Vector3d > positionVector,
                           std::vector < double > timeVector,
                           std::vector < double > deltaVVector,
                           const char * fileName )
{
    std::ofstream exportFile( fileName );
    for ( unsigned int counter = 0; counter < timeVector.size( ); counter++ )
    {
        exportFile << timeVector[ counter ] << ", " << positionVector[ counter ]( 0 ) << ", "
                   << positionVector[ counter ]( 1 ) << ", " << positionVector[ counter ]( 2 )
                   << ", " << deltaVVector[ counter ] << std::endl;
    }
    exportFile.close( );
}

} // namespace transfer_trajectories
} // namespace tudat

