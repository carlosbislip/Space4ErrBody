#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

namespace bislip { //namespace aerodynamics {

class MyGuidance: public tudat::aerodynamics::AerodynamicGuidance
{
public:
    //! Constructor.
    MyGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName,
            const std::vector< double > parameterBounds,
            const double finalMass,
            const double E_max,
            const double Isp,
            const double maxThrust,
            const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_alpha_deg,
            const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_deg,
            const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle ):
        bodyMap_( bodyMap ),
        vehicleName_( vehicleName ),
        parameterBounds_( parameterBounds ),
        finalMass_( finalMass ),
        E_max_( E_max ),
        Isp_( Isp ),
        maxThrust_( maxThrust ),
        interpolator_alpha_deg_( interpolator_alpha_deg ),
        interpolator_eps_T_deg_( interpolator_eps_T_deg ),
        interpolator_throttle_( interpolator_throttle ){ }

    //! Destructor.
    ~MyGuidance( ){ }

    void updateGuidance( const double currentTime );

    //double getE_hat( );

    //double getAoAforCmEqualToZero( const double AoA, const double Mach );

    Eigen::Vector6d getCoefficients( const std::vector< double > &coefficient_input );

    Eigen::Vector3d getCurrentBodyFixedThrustDirection( )
    {
        std::cout << "bodyFixedThrustDirection_: " << bodyFixedThrustDirection_ << std::endl;

        return bodyFixedThrustDirection_;
    }
    double getCurrentThrustMagnitude( )
    {
        std::cout << "currentThrustMagnitude_: " << currentThrustMagnitude_ << std::endl;

        return currentThrustMagnitude_;
    }
    double getCurrentSpecificImpulse( )
    {
        std::cout << "currentSpecificImpulse_: " << currentSpecificImpulse_ << std::endl;

        return currentSpecificImpulse_;
    }
    bool getCurrentEngineStatus( )
    {
        std::cout << "currentEngineStatus_: " << currentEngineStatus_ << std::endl;

        return currentEngineStatus_;
    }



//protected:
    tudat::simulation_setup::NamedBodyMap bodyMap_;
    std::string vehicleName_;
    std::vector< double > parameterBounds_;
    double finalMass_;
    double E_max_;
    double Isp_;
    double maxThrust_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_alpha_deg_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_deg_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle_;
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions_;// = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems_;
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_;




};
/*
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
        interpolator_throttle_( interpolator_throttle ){}
//{       FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
//                bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
//}

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
    //std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions_;
    double currentEpsilon_;
    Eigen::Vector3d bodyFixedThrustDirection_;
    double currentThrustMagnitude_;
    double currentSpecificImpulse_;
    bool currentEngineStatus_;
    //double E_hat_;

};
*/


} // namespace aerodynamics
//} // namespace tudat
