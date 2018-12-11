/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ATMOSPHERE_MODEL_H
#define TUDAT_ATMOSPHERE_MODEL_H

#include <memory>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/Aerodynamics/windModel.h"

namespace tudat
{

namespace aerodynamics
{

//! Enum of all the possible independent variables on which the atmosphere can depend.
/*!
 * Enum of all the possible independent variables on which the atmosphere can depend.
 */
enum AtmosphereIndependentVariables
{
    altitude_dependent_atmosphere = 0,
    longitude_dependent_atmosphere = 1,
    latitude_dependent_atmosphere = 2,
    time_dependent_atmosphere = 3
};

//! Enum of all the possible dependent variables that an atmosphere can describe.
/*!
 * Enum of all the possible dependent variables that an atmosphere can describe.
 */
enum AtmosphereDependentVariables
{
    density_dependent_atmosphere = 0,
    pressure_dependent_atmosphere = 1,
    temperature_dependent_atmosphere = 2,
    gas_constant_dependent_atmosphere = 3,
    specific_heat_ratio_dependent_atmosphere = 4,
    molar_mass_dependent_atmosphere = 5
};

//! Atmosphere model class.
/*!
 * Base class for all atmosphere models.
 * To keep the function generic for both reference atmospheres and standard
 * atmospheres, altitude, longitude, latitude, and time are inputs for all functions.
 */
class AtmosphereModel
{
public:

    //! Default destructor.
    /*!
    * Default destructor.
    */
    virtual ~AtmosphereModel( ) { }

    //! Get local density.
    /*!
    * Returns the local density parameter of the atmosphere in kg per meter^3.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric density.
    */
    virtual double getDensity( const double altitude, const double longitude,
                               const double latitude, const double time ) = 0;

    //! Get local pressure.
    /*!
    * Returns the local pressure of the atmosphere parameter in Newton per meter^2.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric pressure.
    */
    virtual double getPressure( const double altitude, const double longitude,
                                const double latitude, const double time ) = 0;

    //! Get local temperature.
    /*!
    * Returns the local temperature of the atmosphere parameter in Kelvin.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric temperature.
    */
    virtual double getTemperature( const double altitude, const double longitude,
                                   const double latitude, const double time ) = 0;

    //! Get local speed of sound.
    /*!
    * Returns the local speed of sound of the atmosphere in m/s.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric speed of sound.
    */
    virtual double getSpeedOfSound( const double altitude, const double longitude,
                                    const double latitude, const double time ) = 0;

    //! Function to retrieve the model describing the wind velocity vector of the atmosphere
    /*!
     * Function to retrieve the model describing the wind velocity vector of the atmosphere
     * \return Model describing the wind velocity vector of the atmosphere
     */
    std::shared_ptr< WindModel > getWindModel( )
    {
        return windModel_;
    }

    //! Function to set the model describing the wind velocity vector of the atmosphere
    /*!
     * Function to set the model describing the wind velocity vector of the atmosphere
     * \param windModel New model describing the wind velocity vector of the atmosphere
     */
    void setWindModel( const std::shared_ptr< WindModel > windModel )
    {
        windModel_ = windModel;
    }

protected:

    //! Model describing the wind velocity vector of the atmosphere
    std::shared_ptr< WindModel > windModel_;

private:

};

//! Typedef for shared-pointer to AtmosphereModel object.
typedef std::shared_ptr< AtmosphereModel > AtmosphereModelPointer;

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_ATMOSPHERE_MODEL_H
