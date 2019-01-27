#ifndef CUSTOMTERMINATIONFUNCTION_VAL_H
#define CUSTOMTERMINATIONFUNCTION_VAL_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

//#include "Space4ErrBodyProblem.h"
#include <Tudat/Bislip/Variables/bislipVariables.h>


namespace bislip {

bool StopOrNot (const tudat::simulation_setup::NamedBodyMap& bodyMap,
                const std::string &vehicleName,
                const std::vector< double > &vehicleParameterValues,
                const std::vector< double > &terminationConditionsValues);

}

#endif // CUSTOMTERMINATIONFUNCTION_VAL_H

