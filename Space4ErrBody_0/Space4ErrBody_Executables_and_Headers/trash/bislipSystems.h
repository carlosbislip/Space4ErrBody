/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_BISLIPSYSTEMS_H
#define TUDAT_BISLIPSYSTEMS_H

#include <map>
#include <iostream>

#include <memory>

//#include <Tudat/Astrodynamics/SystemModels/engineModel.h>
#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>

#include "bislipVariables.h"

namespace tudat
{

namespace system_models
{

//! Wrapper class that contains the relevant hardware systems of a vehicle.
/*!
 *  Wrapper class that contains the relevant hardware systems of a vehicle. Not all member objects need to be set; nullptr
 *  member objects denote that the vehicle does not contain the associated hardware.
 */
class BislipSystems
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dryMass Total dry mass of the vehicle (not defined; NaN by default).
     */
    BislipSystems( const double landingMass = TUDAT_NAN ):
        landingMass_( landingMass ){ }

    //! Destructor
    ~BislipSystems( ){ }


    void setInitialMass( const double initialMass )
    {
        initialMass_ = initialMass;
    }


    //! Function to retrieve the initial mass of the vehicle
    /*!
     * Function to retrieve the  total dry mass of the vehicle
     * \return initial mass of the vehicle
     */
    double getInitialMass( )
    {
        return initialMass_;
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

    /*
           void setAoAInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Ascent )
           {
               interpolatorAngleOfAttack_Ascent_ = interpolatorAngleOfAttack_Ascent;
           }

           void setAoAInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Descent )
           {
               interpolatorAngleOfAttack_Descent_ = interpolatorAngleOfAttack_Descent;
           }

           void setBankAngleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleAscent )
           {
               interpolatorBankAngleAscent_ = interpolatorBankAngleAscent;
           }

           void setBankAngleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleDescent )
           {
               interpolatorBankAngleDescent_ = interpolatorBankAngleDescent;
           }

           void setThrustElevationAngleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleAscent )
           {
               interpolatorThrustElevationAngleAscent_ = interpolatorThrustElevationAngleAscent;
           }
           void setThrustElevationAngleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleDescent )
           {
               interpolatorThrustElevationAngleDescent_ = interpolatorThrustElevationAngleDescent;
           }

           void setThrustAzimuthAngleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Ascent )
           {
               interpolatorThrustAzimuthAngle_Ascent_ = interpolatorThrustAzimuthAngle_Ascent;
           }

           void setThrustAzimuthAngleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Descent )
           {
               interpolatorThrustAzimuthAngle_Descent_ = interpolatorThrustAzimuthAngle_Descent;
           }

           void setThrottleInterpolator_Ascent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Ascent )
           {
               interpolatorThrottleSetting_Ascent_ = interpolatorThrottleSetting_Ascent;
           }

           void setThrottleInterpolator_Descent( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Descent )
           {
               interpolatorThrottleSetting_Descent_ = interpolatorThrottleSetting_Descent;
           }

       */

    void setAngleOfAttackInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack )
    {
        interpolatorAngleOfAttack_ = interpolatorAngleOfAttack;
    }
    void setBankAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngle )
    {
        interpolatorBankAngle_ = interpolatorBankAngle;
    }
    void setThrustElevationAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngle )
    {
        interpolatorThrustElevationAngle_ = interpolatorThrustElevationAngle;
    }
    void setThrustAzimuthAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle )
    {
        interpolatorThrustAzimuthAngle_ = interpolatorThrustAzimuthAngle;
    }
    void setThrottleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting )
    {
        interpolatorThrottleSetting_ = interpolatorThrottleSetting;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getAngleOfAttackInterpolator( )
    {
        return interpolatorAngleOfAttack_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getBankAngleInterpolator( )
    {
        return interpolatorBankAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustElevationAngleInterpolator( )
    {
        return interpolatorThrustElevationAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustAzimuthAngleInterpolator( )
    {
        return interpolatorThrustAzimuthAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrottleInterpolator( )
    {
        return interpolatorThrottleSetting_;
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
    void setParameterBounds( const std::map< bislip::variables::OptimizationParameter, std::pair < double, double > > &Bounds )
    {

        Bounds_ = Bounds;
        //   if( direction == "ascent" ){ Bounds_Ascent_ = parameterBounds; }
        // if( direction == "descent" ){ Bounds_Descent_ = parameterBounds; }
    }

    //    std::pair < double, double > getParameterBounds( const std::string &parameter )
    std::pair < double, double > getParameterBounds( const bislip::variables::OptimizationParameter &parameter )
    {
        // std::map< std::string, std::pair < double, double > > Bounds;

        //        if( direction == "ascent" ){ Bounds = Bounds_Ascent_; }
        //       if( direction == "descent" ){ Bounds = Bounds_Descent_; }

        return Bounds_[ parameter ];
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
    // double initialMass_;

    //    std::map< std::string, std::pair < double, double > > Bounds_;
    std::map< std::string, std::pair < double, double > > Bounds_Ascent_;
    std::map< std::string, std::pair < double, double > > Bounds_Descent_;
    std::map< bislip::variables::OptimizationParameter, std::pair < double, double > > Bounds_;

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_;


    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleAscent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleAscent_;

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleDescent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleDescent_;
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

    //! Initial mass of the vehicle
    double initialMass_;

    //! Landing mass of the vehicle
    double landingMass_;



};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_BISLIPSYSTEMS_H
