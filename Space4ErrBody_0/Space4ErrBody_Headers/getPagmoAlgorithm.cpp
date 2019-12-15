
#include <iostream>
#include "getPagmoAlgorithm.h"
/*
#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/pso.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/sga.hpp"
#include "pagmo/algorithms/sea.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/nlopt.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "pagmo/algorithms/xnes.hpp"
*/
//#include "pagmo/algorithms/nsga2.hpp"
//#include "pagmo/algorithms/moead.hpp"
//#include "pagmo/algorithms/ihs.hpp"

#include <pagmo/pagmo.hpp>

namespace bislip {

pagmo::algorithm getPagmoAlgorithm( const bool customParameters, const Eigen::VectorXd optimizerSpecificParameterValues )
{
    //! Default multi-objective algorithm. (Arbitrary)
    //pagmo::algorithm algo{ pagmo::nsga2( ) };

    switch( int( customParameters ) )
    {
    case 0 :
    {
        switch( int( optimizerSpecificParameterValues( 0 ) ) )
        {
        case 0:
        {
            pagmo::algorithm algo{ pagmo::nsga2( ) };
            return algo;
            // break;
        }
        case 1:
        {
            pagmo::algorithm algo{ pagmo::moead( ) };
            return algo;
            //break;
        }
        case 2:
        {
            pagmo::algorithm algo{ pagmo::ihs( ) };
            return algo;
            //break;
        }
        default:
        {
            throw std::runtime_error( "Error, multi-objective pagmo algorithm " + std::to_string( int( optimizerSpecificParameterValues( 0 ) ) ) + " was not found." );
        }
        }
    }
    case 1 :
    {
        switch( int( optimizerSpecificParameterValues( 0 ) ) )
        {
        case 0:
        {
            pagmo::algorithm algo{ pagmo::nsga2( 1u, optimizerSpecificParameterValues( 1 ), optimizerSpecificParameterValues( 2 ), optimizerSpecificParameterValues( 3 ), optimizerSpecificParameterValues( 4 ) ) };
            return algo;
            //break;
        }
        case 1:
        {
            const std::vector< std::string > weightGenerationVector = { "grid", "low discrepancy", "random" };
            const std::vector< std::string > decompositionVector = { "weighted", "tchebycheff", "bi" };
            const std::string weightGeneration = weightGenerationVector[ int( optimizerSpecificParameterValues( 1 ) ) ];
            const std::string decomposition = decompositionVector[ int( optimizerSpecificParameterValues( 2 ) ) ];
            const std::vector< std::vector< double > >::size_type neighbours = int( optimizerSpecificParameterValues( 3 ) );
            const double CR = optimizerSpecificParameterValues( 4 );
            const double F = optimizerSpecificParameterValues( 5 );
            const double eta_m = optimizerSpecificParameterValues( 6 );
            const bool preserve_diversity = bool( optimizerSpecificParameterValues( 7 ) );
            const double realb = optimizerSpecificParameterValues( 8 );
            const unsigned int limit = int( optimizerSpecificParameterValues( 9 ) );

            pagmo::algorithm algo{ pagmo::moead( 1u,
                                                 weightGeneration,
                                                 decomposition,
                                                 neighbours,
                                                 CR, F, eta_m, realb, limit, preserve_diversity ) };
            return algo;
            //break;
        }
        case 2:
        {
            const double phmcr = optimizerSpecificParameterValues( 1 );
            const double ppar_min = optimizerSpecificParameterValues( 2 );
            const double ppar_max = optimizerSpecificParameterValues( 3 );
            const double bw_min = optimizerSpecificParameterValues( 4 );
            const double bw_max = optimizerSpecificParameterValues( 5 );
            pagmo::algorithm algo{ pagmo::ihs( 1u,
                                               phmcr,
                                               ppar_min, ppar_max,
                                               bw_min, bw_max ) };
            return algo;
            //break;
        }
        default:
        {
            throw std::runtime_error( "Error, multi-objective pagmo algorithm " + std::to_string( int( optimizerSpecificParameterValues[ 0 ] ) ) + " was not found." );
        }
        }
        //      break;
    }
    default:
    {
        throw std::runtime_error( "Error, Conditional must be 0 (generic parameters) or 1 (custom parameters). Value provided was " + std::to_string( int( customParameters ) ) + "." );
    }
    }
}
} // namespace bislip
