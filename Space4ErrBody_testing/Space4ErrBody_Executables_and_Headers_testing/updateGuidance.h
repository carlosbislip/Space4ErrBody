#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

//! Trying to implement a simple aerodynamic guidance. Initially taken from the
//! TUDAT website.
//! http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/aerodynamicGuidance.html#FlightConditions
namespace tudat
{
namespace aerodynamics
{
class FlightConditionsBasedAerodynamicGuidance: public AerodynamicGuidance
{
public:
    FlightConditionsBasedAerodynamicGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName )
    {
        FlightConditions_ = boost::dynamic_pointer_cast< AtmosphericFlightConditions >(
                    bodyMap.at( vehicleName )->getFlightConditions( ) );
    }

    void updateGuidance( const double currentTime );

private:

    boost::shared_ptr< AtmosphericFlightConditions > FlightConditions_;
};

} // namespace aerodynamics

} // namespace tudat
