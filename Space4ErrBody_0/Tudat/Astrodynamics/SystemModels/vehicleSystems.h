/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VEHICLESYSTEMS_H
#define TUDAT_VEHICLESYSTEMS_H

#include <map>
#include <iostream>

#include <memory>

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"

namespace tudat
{

namespace system_models
{

//! Wrapper class that contains the relevant hardware systems of a vehicle.
/*!
 *  Wrapper class that contains the relevant hardware systems of a vehicle. Not all member objects need to be set; nullptr
 *  member objects denote that the vehicle does not contain the associated hardware.
 */
class VehicleSystems
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dryMass Total dry mass of the vehicle (not defined; NaN by default).
     */
    VehicleSystems( const double dryMass = TUDAT_NAN ):
        dryMass_( dryMass ){ }

    //! Destructor
    ~VehicleSystems( ){ }

    //! Function to retrieve the engine models
    /*!
     * Function to retrieve the engine models
     * \return Named list of engine models in the vehicle
     */
    std::map< std::string, std::shared_ptr< EngineModel > > getEngineModels( )
    {
        return engineModels_;
    }

    //! Function to set a single engine in the vehicle
    /*!
     * Function to set a single engine in the vehicle. Each engine can be identified by a string. If only a single
     * engine is set, the default (empty string) can be used.
     * \param engineModel Model of engine that is to be set
     * \param engineName Reference id of the engine that is to be set.
     */
    void setEngineModel(
            const std::shared_ptr< EngineModel > engineModel, const std::string engineName = "" )
    {
        // Check if engine with this name already exists.
        if( engineModels_.count( engineName ) )
        {
            std::cerr << "Warning, engine model of name " << engineModel << " already exists, overriding old model" << std::endl;
        }

        engineModels_[ engineName ] = engineModel;
    }

    //! Function to retrieve the total dry mass of the vehicle
    /*!
     * Function to retrieve the  total dry mass of the vehicle
     * \return Total dry mass of the vehicle
     */
    double getDryMass( )
    {
        return dryMass_;
    }

    //! Function to set the current deflection of a single control surface
    /*!
     * Function to set the current deflection of a single control surface
     * \param controlSurfaceId Name of control surface for which deflection is to be set
     * \param deflectionAngle Current deflection of control surface that is to be set.
     */
    void setCurrentControlSurfaceDeflection(
            const std::string& controlSurfaceId, const double deflectionAngle )
    {
        currentControlSurfaceDeflections_[ controlSurfaceId ] =  deflectionAngle;
    }

    //! Function to retrieve the current deflection of a single control surface
    /*!
     * Function to retrieve the current deflection of a single control surface
     * \param controlSurfaceId Name of control surface for which deflection is to be set
     * \return Current deflection of control surface that is requested.
     */
    double getCurrentControlSurfaceDeflection(
            const std::string& controlSurfaceId )
    {
        return currentControlSurfaceDeflections_.at( controlSurfaceId );
    }

    //! Function to (re)set the vehicle nose radius
    /*!
     * Function to (re)set the vehicle nose radius
     * \param noseRadius The  vehicle nose radius that is to be set
     */
    void setNoseRadius( const double noseRadius )
    {
        noseRadius_ = noseRadius;
    }

    //! Function to retrieve the vehicle nose radius
    /*!
     * Function to retrieve the vehicle nose radius
     * \return The vehicle nose radius
     */
    double getNoseRadius( )
    {
        return noseRadius_;
    }

    //! Function to (re)set the vehicle wall emissivity
    /*!
     * Function to (re)set the vehicle wall emissivity
     * \param wallEmissivity The vehicle wall emissivity that is to be set
     */
    void setWallEmissivity( const double wallEmissivity )
    {
        wallEmissivity_ = wallEmissivity;
    }

    //! Function to retrieve the vehicle wall emissivity
    /*!
     * Function to retrieve the vehicle wall emissivity
     * \return The vehicle wall emissivity
     */
    double getWallEmissivity( )
    {
        return wallEmissivity_;
    }

    void setMaxThrust( const double maxThrust )
    {
        maxThrust_ = maxThrust;
    }

    double getMaxThrust( )
    {
        return maxThrust_;
    }

    void setSpecificImpulse( const double specificImpulse )
    {
        specificImpulse_ = specificImpulse;
    }

    double getSpecificImpulse( )
    {
        return specificImpulse_;
    }

    void setWingSweepAngle( const double lambda )
    {
        lambda_ = lambda;
    }

    double getWingSweepAngle( )
    {
        return lambda_;
    }

    void setLocalBodyAngle( const double phi )
    {
        phi_ = phi;
    }

    double getLocalBodyAngle( )
    {
        return phi_;
    }

    void setTransitionDistance( const double x_T )
    {
        x_T_ = x_T;
    }

    double getTransitionDistance( )
    {
        return x_T_;
    }

    void setE_max( const double E_max )
    {
        E_max_ = E_max;
    }

    double getE_max( )
    {
        return E_max_;
    }

    void setAoAInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_AoA_Ascent )
    {
        interpolator_AoA_Ascent_ = interpolator_AoA_Ascent;
    }

    void setAoAInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_AoA_Descent )
    {
        interpolator_AoA_Descent_ = interpolator_AoA_Descent;
    }

    void setBankAngleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_sigma_Ascent )
    {
        interpolator_sigma_Ascent_ = interpolator_sigma_Ascent;
    }

    void setBankAngleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_sigma_Descent )
    {
        interpolator_sigma_Descent_ = interpolator_sigma_Descent;
    }


    void setThrottleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle_Ascent )
    {
        interpolator_throttle_Ascent_ = interpolator_throttle_Ascent;
    }

    void setThrottleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle_Descent )
    {
        interpolator_throttle_Descent_ = interpolator_throttle_Descent;
    }
    void setThrustElevationAngleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_Ascent )
    {
        interpolator_eps_T_Ascent_ = interpolator_eps_T_Ascent;
    }
    void setThrustElevationAngleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_Descent )
    {
        interpolator_eps_T_Descent_ = interpolator_eps_T_Descent;
    }

    void setThrustAzimuthAngleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_phi_T_Ascent )
    {
        interpolator_phi_T_Ascent_ = interpolator_phi_T_Ascent;
    }

    void setThrustAzimuthAngleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_phi_T_Descent )
    {
        interpolator_phi_T_Descent_ = interpolator_phi_T_Descent;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getAoAInterpolator_Ascent( )
    {
        return interpolator_AoA_Ascent_;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getAoAInterpolator_Descent( )
    {
        return interpolator_AoA_Descent_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getBankAngleInterpolator_Ascent( )
    {
        return interpolator_sigma_Ascent_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getBankAngleInterpolator_Descent( )
    {
        return interpolator_sigma_Descent_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrottleInterpolator_Ascent( )
    {
        return interpolator_throttle_Ascent_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrottleInterpolator_Descent( )
    {
        return interpolator_throttle_Descent_;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustElevationAngleInterpolator_Ascent( )
    {
        return interpolator_eps_T_Ascent_;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustElevationAngleInterpolator_Descent( )
    {
        return interpolator_eps_T_Descent_;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustAzimuthAngleInterpolator_Ascent( )
    {
        return interpolator_phi_T_Ascent_;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustAzimuthAngleInterpolator_Descent( )
    {
        return interpolator_phi_T_Descent_;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getGuidanceInterpolator(
            const double &flight_path_angle,
            const std::string &parameter )
    {

        std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator;

        if ( flight_path_angle >  0 )
        {
            if (parameter == "Angle of Attack")
            {
                interpolator = getThrustElevationAngleInterpolator_Ascent();
            }
            if (parameter == "Bank Angle")
            {
                interpolator = getBankAngleInterpolator_Ascent();
            }
            if (parameter == "Thrust Elevation Angle")
            {
                interpolator = getThrustElevationAngleInterpolator_Ascent();
            }
            if (parameter == "Thrust Azimuth Angle")
            {
                interpolator = getThrustAzimuthAngleInterpolator_Ascent();
            }
            if (parameter == "Throttle Setting")
            {
                interpolator = getThrottleInterpolator_Ascent();
            }
        }
        else
        {
            if (parameter == "Angle of Attack")
            {
                interpolator = getThrustElevationAngleInterpolator_Descent();
            }
            if (parameter == "Bank Angle")
            {
                interpolator = getBankAngleInterpolator_Descent();
            }
            if (parameter == "Thrust Elevation Angle")
            {
                interpolator = getThrustElevationAngleInterpolator_Descent();
            }
            if (parameter == "Thrust Azimuth Angle")
            {
                interpolator = getThrustAzimuthAngleInterpolator_Descent();
            }
            if (parameter == "Throttle Setting")
            {
                interpolator = getThrottleInterpolator_Descent();
            }
        }

        return interpolator;
    }

    void setTargetLat( const double target_lat )
    {
        target_lat_ = target_lat;
    }

    void setTargetLon( const double target_lon )
    {
        target_lon_ = target_lon;
    }

    void setInitialLat( const double initial_lat )
    {
        initial_lat_ = initial_lat;
    }

    void setInitialLon( const double initial_lon )
    {
        initial_lon_ = initial_lon;
    }
    void setInitialDistanceToTarget( const double initial_d_to_target )
    {
        initial_d_to_target_ = initial_d_to_target;
    }
    void setFinalDistanceToTarget( const double final_d_to_target )
    {
        final_d_to_target_ = final_d_to_target;
    }
    double getTargetLat( )
    {
        return target_lat_;
    }
    double getTargetLon( )
    {
        return target_lon_;
    }
    double getInitialLat( )
    {
        return initial_lat_;
    }
    double getInitialLon( )
    {
        return initial_lon_;
    }
    double getInitialDistanceToTarget( )
    {
        return initial_d_to_target_;
    }
    double getFinalDistanceToTarget( )
    {
        return final_d_to_target_;
    }
    void setParameterBounds( const std::map< std::string, std::pair < double, double > > &parameterBounds, const std::string &direction )
    {
        if( direction == "ascent" ){ Bounds_Ascent_ = parameterBounds; }
        if( direction == "descent" ){ Bounds_Descent_ = parameterBounds; }
    }

    std::pair < double, double > getParameterBounds( const std::string &direction, const std::string &parameter )
    {
        std::map< std::string, std::pair < double, double > > Bounds;

        if( direction == "ascent" ){ Bounds = Bounds_Ascent_; }
        if( direction == "descent" ){ Bounds = Bounds_Descent_; }

        return Bounds[ parameter ];
    }

    void setInitialCoordinates( const std::pair < double, double > &initialCoordinates )
    {
        initialCoordinates_ = initialCoordinates;
    }

    std::pair < double, double > getInitialCoordinates( )
    {
        return initialCoordinates_;
    }

    void setTargetCoordinates( const std::pair < double, double > &targetCoordinates )
    {
        targetCoordinates_ = targetCoordinates;
    }

    std::pair < double, double > getTargetCoordinates( )
    {
        return targetCoordinates_;
    }
    void setStartingEpoch( const double startingEpoch )
    {
        startingEpoch_ = startingEpoch;
    }
    double getStartingEpoch( )
    {
        return startingEpoch_;
    }


private:

    std::pair < double, double > initialCoordinates_;
    std::pair < double, double > targetCoordinates_;

    double lambda_;
    double phi_;
    double x_T_;
    double E_max_;
    double maxThrust_;
    double specificImpulse_;
    double startingEpoch_;

    std::map< std::string, std::pair < double, double > > Bounds_Ascent_;
    std::map< std::string, std::pair < double, double > > Bounds_Descent_;

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_AoA_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_phi_T_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_sigma_Ascent_;

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_AoA_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_phi_T_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_sigma_Descent_;

    //! Passing target coordinates
    double target_lat_;
    double target_lon_;

    //! Passing initial coordinates
    double initial_lat_;
    double initial_lon_;

    //! Passing initial distance to target
    double initial_d_to_target_;

    //! Passing final distance to target
    double final_d_to_target_;

    //! Named list of engine models in the vehicle
    std::map< std::string, std::shared_ptr< EngineModel > > engineModels_;

    //! Total dry mass of the vehicle
    double dryMass_;

    //! List if current control surface deflections (with key the control surface id).
    std::map< std::string, double > currentControlSurfaceDeflections_;

    //! Nose radius of the vehicle (used for heating computations)
    double noseRadius_;

    //! Wall emissivity of the vehicle (used for heating computations)
    double wallEmissivity_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_VEHICLESYSTEMS_H
