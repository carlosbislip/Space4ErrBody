#ifndef CUSTOMTERMINATIONFUNCTION_VAL_H
#define CUSTOMTERMINATIONFUNCTION_VAL_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

namespace bislip {

bool StopOrNot ( const tudat::simulation_setup::NamedBodyMap& bodyMap,
                 const std::string &vehicleName,
                 const std::vector< double > &parameterBounds,
                 const std::vector< double > &terminationConditionsValues_,
                 const std::vector< double > &additional_data,
                 const double &E_max,
                 const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_throttle,
                 const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &interpolator_eps_T_deg);

}

#endif // CUSTOMTERMINATIONFUNCTION_VAL_H

