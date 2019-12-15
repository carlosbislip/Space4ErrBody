#ifndef GETPAGMOALGORITHM_H
#define GETPAGMOALGORITHM_H

#include <iostream>

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

pagmo::algorithm getPagmoAlgorithm( const bool customParameters, const Eigen::VectorXd optimizerSpecificParameterValues );

} // namespace bislip
#endif // GETPAGMOALGORITHM_H
