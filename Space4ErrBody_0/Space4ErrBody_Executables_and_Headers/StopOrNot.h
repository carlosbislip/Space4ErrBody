#ifndef CUSTOMTERMINATIONFUNCTION_VAL_H
#define CUSTOMTERMINATIONFUNCTION_VAL_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

bool StopOrNot ( const tudat::simulation_setup::NamedBodyMap& bodyMap,
                 const std::string &vehicleName,
                 const std::vector< double > term_cond,
                 const std::vector< double > additional_data );



#endif // CUSTOMTERMINATIONFUNCTION_VAL_H

