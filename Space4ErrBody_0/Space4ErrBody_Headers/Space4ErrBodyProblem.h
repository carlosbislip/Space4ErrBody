/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef SPACE4ERRBODY_HPP
#define SPACE4ERRBODY_HPP
/*
#include <utility>
#include <vector>
#include <limits>
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

#include <pagmo/problem.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/island.hpp>
#include <pagmo/serialization.hpp>
#include <pagmo/types.hpp>

//#include "getStuff.h"
//#include <Tudat/Bislip/Parameters/bislipOptimizationParameters.h>
//#include <Tudat/Bislip/Variables/bislipVariables.h>
//#include <Tudat/Bislip/SystemModels/bislipSystems.h>

#include <boost/format.hpp>


//#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
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
*/
//#include "updateGuidance.h"

#include <pagmo/pagmo.hpp>
#include <pagmo/utils/multi_objective.hpp>
#include <pagmo/serialization.hpp>
#include <pagmo/detail/constants.hpp>
#include <pagmo/exceptions.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "getOutputPath.h"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/InputOutput/solarActivityData.h>
#include <Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h>
#include <Tudat/Astrodynamics/Aerodynamics/nrlmsise00InputFunctions.h>
#include <Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h>

//#include "getPagmoAlgorithm.h"
//#include "bislipUtilities.h"
//#include "applicationOutput_tudat.h"
//#include "applicationOutput_pagmo.h"

#include <Tudat/Bislip/bislipProblemInput.h>
#include <Tudat/Bislip/bislipFitnessVector.h>

namespace bislip
{//using namespace pagmo::problem;

//! Define the problem PaGMO-style
//struct Space4ErrBodyProblem
class Space4ErrBody : public pagmo::problem
{
/*
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar)
    {
        ar & problemInput_;
        ar & bodyMap_;
    }
*/

public:

    //! Empty constructor
    Space4ErrBody( ) { }

    //! Default constructor, required for Pagmo compatibility
    Space4ErrBody(
            std::shared_ptr< bislip::problem_input > &problemInput,
            tudat::simulation_setup::NamedBodyMap& bodyMap ):
        problemInput_( problemInput ), bodyMap_( bodyMap )
    { }

    virtual ~Space4ErrBody( ) { }

    //! Calculate the fitness as a function of the parameter vector x
    std::vector< double > fitness( const std::vector< double >  &x ) const;

    //! Retrieve the allowable limits of the parameters to vary: pair
    //! containing minima and maxima of parameter values
    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    //! Retrieve the name of the problem
    std::string get_name( ) const;

    std::vector< double >::size_type get_nobj() const
    {
        char objectiveFunctionCase = problemInput_->getObjectiveFunctionCase();
        unsigned int nOF = 1u;

        std::vector< bislip::fitness_vector::cost > costList = bislip::fitness_vector::getCostList( objectiveFunctionCase );
        std::vector< bislip::fitness_vector::penalty > penaltyList = bislip::fitness_vector::getPenaltyList( objectiveFunctionCase );

        nOF = static_cast< unsigned int >( costList.size() );

        if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_with_descent ) )
        { nOF += 2 * static_cast< unsigned int >( penaltyList.size() ); }
        else
        { nOF += static_cast< unsigned int >( penaltyList.size() ); }

        problemInput_->setFitnessVectorLength( nOF );

        std::cout << "    Objective Function Case: " << objectiveFunctionCase << std::endl;
        std::cout << "    Size of Fitness Vector : " << std::to_string( nOF ) << std::endl;
        std::cout << "    " << std::endl;

        /*
        if( objectiveFunctionCase == 'A' ) { nOF = 3u; }
        else if( objectiveFunctionCase == 'B' ) { nOF = 3u; }
        else if( objectiveFunctionCase == 'C' ) { nOF = 5u; }
        else if( objectiveFunctionCase == 'D' ) { nOF = 7u; }
        else if( objectiveFunctionCase == 'E' ) { nOF = 9u; }
        else if( objectiveFunctionCase == 'F' ) { nOF = 5u; }
        else if( objectiveFunctionCase == 'G' ) { nOF = 3u; }
        else if( objectiveFunctionCase == 'H' ) { nOF = 7u; }
        else if( objectiveFunctionCase == 'I' ) { nOF = 9u; }
        else if( objectiveFunctionCase == 'J' ) { nOF = 7u; }
        else if( objectiveFunctionCase == 'K' ) { nOF = 5u; }
        else if( objectiveFunctionCase == 'L' ) { nOF = 11u; }
        else if( objectiveFunctionCase == 'M' ) { nOF = 31u; }

        if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_with_descent ) == false )
        { nOF = 1u + ( nOF - 1u ) / 2u; } */

        return nOF;
    }


    std::vector< double >::size_type get_nec() const
    {
        return 0u;
    }

    //const std::vector< std::vector< double > > problemBounds_ = problemInput_->getDecisionVectorBounds();
    /*
    //! Serialization function for Pagmo compatibility
    template < typename Archive >
    void serialize( Archive &ar )
    {
        ar( problemBounds_ );
    }
*/
private:

    const std::shared_ptr< bislip::problem_input > problemInput_;
    const tudat::simulation_setup::NamedBodyMap bodyMap_;



};
}
PAGMO_REGISTER_PROBLEM(bislip::Space4ErrBody)

#endif // SPACE4ERRBODY_OPTIMIZATION_HPP
