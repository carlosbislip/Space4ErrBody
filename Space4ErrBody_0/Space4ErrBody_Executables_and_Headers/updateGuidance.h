#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

//! http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/aerodynamicGuidance.html#FlightConditions
namespace bislip { //namespace aerodynamics {

class MyAerodynamicGuidance: public tudat::aerodynamics::AerodynamicGuidance
{
public:
    MyAerodynamicGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName,
            const double E_max,
            const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_alpha_deg ):
        bodyMap_( bodyMap ),
        vehicleName_( vehicleName ),
        E_max_( E_max ),
        interpolator_alpha_deg_( interpolator_alpha_deg )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap.at( vehicleName )->getFlightConditions( ) );
        vehicleSystems_ = bodyMap.at( vehicleName )->getVehicleSystems( );
        coefficientInterface_ = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                    bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );
    }

    void updateGuidance( const double currentTime );

    //double getAoAforCmEqualToZero( const double AoA, const double Mach );

    Eigen::Vector6d getCoefficients( const std::vector< double > &coefficient_input );
    double getE_hat( );

private:
    tudat::simulation_setup::NamedBodyMap bodyMap_;
    std::string vehicleName_;
    double E_max_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_alpha_deg_;
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions_;
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems_;
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_;
    double E_hat_;

};

class MyThrustGuidance: public tudat::aerodynamics::AerodynamicGuidance
{
public:
    MyThrustGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName,
            const double E_max,
            const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_deg,
            const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle ):
        bodyMap_( bodyMap ),
        vehicleName_( vehicleName ),
        E_max_( E_max ),
        interpolator_eps_T_deg_( interpolator_eps_T_deg ),
        interpolator_throttle_( interpolator_throttle ) { }

    //! Destructor.
    ~MyThrustGuidance( ){ }

    void updateGuidance( const double currentTime );

    Eigen::Vector3d getCurrentBodyFixedThrustDirection( );

    double getCurrentThrustMagnitude( );

    double getCurrentSpecificImpulse( );

    bool getCurrentEngineStatus( );

private:

    tudat::simulation_setup::NamedBodyMap bodyMap_;
    std::string vehicleName_;
    double E_max_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_deg_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle_;
    double currentEpsilon_;
    Eigen::Vector3d bodyFixedThrustDirection_;
    double currentThrustMagnitude_;
    double currentSpecificImpulse_;
    bool currentEngineStatus_;
    double E_hat_;

};



} // namespace aerodynamics
//} // namespace tudat
