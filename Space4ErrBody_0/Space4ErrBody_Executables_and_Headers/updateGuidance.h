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
class MyAerodynamicGuidance: public AerodynamicGuidance
{
public:
    MyAerodynamicGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName,
            const Eigen::VectorXd xn,
            const Eigen::VectorXd alpha_rad,
            const Eigen::VectorXd eps_T_rad,
            const Eigen::VectorXd throttle,
            const Eigen::VectorXd E_mapped )
    {
        bodyMap_ = bodyMap;
        FlightConditions_ = boost::dynamic_pointer_cast< AtmosphericFlightConditions >(
                    bodyMap.at( vehicleName )->getFlightConditions( ) );
        vehicleName_ = vehicleName;
        xn_ = xn;
        alpha_rad_ = alpha_rad;
        eps_T_rad_ = eps_T_rad;
        throttle_ = throttle;
        E_mapped_ = E_mapped;
        vehicleSystems_ = bodyMap.at( vehicleName )->getVehicleSystems( );
        coefficientInterface_ = boost::dynamic_pointer_cast< AerodynamicCoefficientInterface >(
                    bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );
    }

    void updateGuidance( const double currentTime );
    //double getAoAforCmEqualToZero( const double AoA, const double Mach );
    Eigen::Vector6d getCoefficients( const double AoA, const double Mach );
    /* double getF (
            unsigned i,
            unsigned j,
            const Eigen::VectorXd &p,
            const double &eps,
            Eigen::Vector2d &x );*/

private:

    boost::shared_ptr< AtmosphericFlightConditions > FlightConditions_;
    boost::shared_ptr< system_models::VehicleSystems > vehicleSystems_;
    tudat::simulation_setup::NamedBodyMap bodyMap_;
    std::string vehicleName_;
    boost::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface_;
    Eigen::VectorXd xn_;
    Eigen::VectorXd alpha_rad_;
    Eigen::VectorXd eps_T_rad_;
    Eigen::VectorXd throttle_;
    Eigen::VectorXd E_mapped_;

};

} // namespace aerodynamics

} // namespace tudat
