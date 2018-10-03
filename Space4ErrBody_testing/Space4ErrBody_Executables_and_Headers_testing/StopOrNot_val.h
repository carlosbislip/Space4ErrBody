#ifndef CUSTOMTERMINATIONFUNCTION_VAL_H
#define CUSTOMTERMINATIONFUNCTION_VAL_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
namespace tudat
{

bool StopOrNot ( const simulation_setup::NamedBodyMap& bodyMap, const std::string vehicleName);

} // namespace tudat


#endif // CUSTOMTERMINATIONFUNCTION_VAL_H

