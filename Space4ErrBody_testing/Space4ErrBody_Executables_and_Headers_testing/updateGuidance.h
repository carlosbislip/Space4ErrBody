#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>



using namespace tudat::aerodynamics;
using namespace tudat::mathematical_constants;
using namespace tudat::simulation_setup; //! required for NamedBodyMap

class FlightConditionsBasedAerodynamicGuidance: public AerodynamicGuidance
{

    FlightConditionsBasedAerodynamicGuidance(
            const NamedBodyMap& bodyMap,
            const std::string vehicleName )
    {
        vehicleFlightConditions_ = bodyMap.at( vehicleName )->getFlightConditions( );
    }

    void updateGuidance( const double currentTime );

private:

    boost::shared_ptr< FlightConditions > vehicleFlightConditions_;
};
