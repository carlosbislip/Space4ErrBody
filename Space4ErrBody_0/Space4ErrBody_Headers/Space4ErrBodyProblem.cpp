/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved.
 * qqwqwqw
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 *    Notes
 *
 */

#include "Space4ErrBodyProblem.h"
#include <Tudat/Bislip/bislipHeaders.h>
#include <Tudat/Bislip/bislipUtilities.h>
#include <Tudat/Bislip/bislipPenaltyFunctions.h>
#include <Tudat/Bislip/bislipConstraints.h>
#include <Tudat/Bislip/updateGuidance.h>
#include <Tudat/Bislip/bislipProblemInput.h>
#include <Tudat/Bislip/bislipDecisionVectorEvaluation.h>
#include <Tudat/Bislip/createTUDATSettings.h>

namespace bislip{


/*Space4ErrBodyProblem::Space4ErrBodyProblem(
        const std::shared_ptr< bislip::problem_input > &problemInput,
        const tudat::simulation_setup::NamedBodyMap& bodyMap ):
    problemInput_( problemInput ),
    bodyMap_( bodyMap )
{ }
*/
//! Descriptive name of the problem
std::string Space4ErrBody::get_name() const
{
    return problemInput_->getProblemName();
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > Space4ErrBody::get_bounds() const
{
    return { problemInput_->getDecisionVectorBounds()[ 0 ], problemInput_->getDecisionVectorBounds()[ 1 ] };
}

std::vector< double > Space4ErrBody::fitness( const std::vector< double > &x ) const
{

    const std::string vehicleName =  problemInput_->getVehicleName();

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName )->getBislipSystems( ) ;

    unsigned int debugInfo = bislipSystems->getDebugInfo();

    std::vector< double > x_mod = x;
    std::vector< std::vector< double > > decisionVectorBounds = problemInput_->getDecisionVectorBounds( );

    if( bislipSystems->getFlag( bislip::flags::flag_list::impose_parameter_bounds_on_interpolators ) )
    {
        if( debugInfo == 1 ){ std::cout << "Verifying that the decision vector is both within bounds and doesnt contain nan values" << std::endl; }
        for( unsigned int i = 0; i < x.size( ); i++ )
        {
            if( std::isnan( x_mod[ i ] ) )
            { x_mod[ i ] = ( decisionVectorBounds[ 0 ][ i ] + decisionVectorBounds[ 1 ][ i ] ) / 2; }
            else
            {
                if( x_mod[ i ] < decisionVectorBounds[ 0 ][ i ])
                { x_mod[ i ] = decisionVectorBounds[ 0 ][ i ]; }
                else if( x_mod[ i ] > decisionVectorBounds[ 1 ][ i ] )
                { x_mod[ i ] = decisionVectorBounds[ 1 ][ i ]; }
            }
        }
    }

    if( debugInfo == 1 ){ std::cout << "Evaluating Decision Vector" << std::endl; }
    Eigen::MatrixXd depVarTimeHistoryMatrix;
    unsigned int rowsAscent;
    bislip::Variables::decisionVectorEvaluation( x_mod , problemInput_, bodyMap_, depVarTimeHistoryMatrix, rowsAscent );

    const unsigned int rowsTotal = static_cast< unsigned int >( depVarTimeHistoryMatrix.rows() );

    const double propagationStepSize = bislipSystems->getPropagationStepSize();

    if( debugInfo == 1 ){ std::cout << "Extract various variable time histories." << std::endl; }


    Eigen::VectorXd depVar_TimeOfFlight                        = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::time_of_flight );
    Eigen::VectorXd depVar_Height                              = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::height_dependent_variable );
    Eigen::VectorXd depVar_Airspeed                            = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::airspeed_dependent_variable );
    Eigen::VectorXd depVar_NormalizedSpecificEnergy            = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::normalized_specific_energy );
    Eigen::VectorXd depVar_NormalizedSpecificEnergy_Ascent     = depVar_NormalizedSpecificEnergy.segment( 0 , rowsAscent );
    Eigen::VectorXd depVar_NormalizedSpecificEnergy_Descent    = depVar_NormalizedSpecificEnergy.segment( rowsAscent - 1 , rowsTotal - rowsAscent );
    Eigen::VectorXd depVar_DynamicPressure                     = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::local_dynamic_pressure_dependent_variable );
    Eigen::VectorXd depVar_BendingMoment                       = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::bending_moment );
    Eigen::VectorXd depVar_CurrentMass                         = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::current_mass );
    Eigen::VectorXd depVar_FlightPathAngle                     = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::flight_path_angle );
    Eigen::VectorXd depVar_CentralTargetAngularDistanceToGo    = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::central_target_angular_distance_to_go );
    Eigen::VectorXd depVar_BodyFrameTotalGLoad_x               = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::body_frame_total_g_load_x );
    Eigen::VectorXd depVar_BodyFrameTotalGLoad_z               = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::body_frame_total_g_load_z );
    Eigen::VectorXd depVar_BodyFrameMechanicalGLoadMag         = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::body_frame_total_mechanical_g_load_magnitude );
    Eigen::VectorXd depVar_PassengerFrameTotalGLoad_x          = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::passenger_frame_total_g_load_x );
    Eigen::VectorXd depVar_PassengerFrameTotalGLoad_z          = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::passenger_frame_total_g_load_z );
    Eigen::VectorXd depVar_HeatFluxChapmanNose                 = depVarTimeHistoryMatrix.col( bislip::dependent_variables::dependent_variable_list::heat_flux_chapman );
    Eigen::VectorXd depVar_PassengerFrame_Jerk_x, depVar_PassengerFrame_Jerk_y, depVar_PassengerFrame_Jerk_z, depVar_IntegratedHeatLoad;


    depVar_PassengerFrame_Jerk_x = depVarTimeHistoryMatrix.col( depVarTimeHistoryMatrix.cols() - 4u );
    depVar_PassengerFrame_Jerk_y = depVarTimeHistoryMatrix.col( depVarTimeHistoryMatrix.cols() - 3u );
    depVar_PassengerFrame_Jerk_z = depVarTimeHistoryMatrix.col( depVarTimeHistoryMatrix.cols() - 2u );
    depVar_IntegratedHeatLoad = depVarTimeHistoryMatrix.col( depVarTimeHistoryMatrix.cols() - 1u );

    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }
    if( debugInfo == 1 ){ std::cout << "angularDistanceToGo_TerminationSettings         = " << depVar_CentralTargetAngularDistanceToGo( depVar_CentralTargetAngularDistanceToGo.size() - 1 ) << std::endl; }
    //if( debugInfo == 1 ){ std::cout << "netAngularDistanceTravelled_TerminationSettings = " << depVar_AngularDistanceTravelled( depVar_AngularDistanceTravelled.size() - 1 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "minimumHeightAllowable_TerminationSettings      = " << depVar_Height( depVar_Height.size() - 1 ) << std::endl; }
    //if( debugInfo == 1 ){ std::cout << "maximumTimeOfFlight_TerminationSettings         = " << depVar_TimeOfFlight( depVar_TimeOfFlight.size() - 1 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }

    //! Extract Dependent variables of final state.
    const Eigen::VectorXd dependentVariableFinalState_Ascent  = depVarTimeHistoryMatrix.row( rowsAscent - 1 );
    const Eigen::VectorXd dependentVariableFinalState_Descent = depVarTimeHistoryMatrix.row( rowsTotal - 1 );

    //std::cout << "dependentVariableFinalState_Ascent  = " << dependentVariableFinalState_Ascent.transpose()  << std::endl;
    //std::cout << "dependentVariableFinalState_Descent = " << dependentVariableFinalState_Descent.transpose()  << std::endl;
    //const double altitude_f_calc = depVar_FINALSTATE[3];
    const double targetLat_deg_calc  = depVarTimeHistoryMatrix( rowsTotal - 1, bislip::dependent_variables::dependent_variable_list::latitude_angle );
    const double targetLon_deg_calc  = depVarTimeHistoryMatrix( rowsTotal - 1, bislip::dependent_variables::dependent_variable_list::longitude_angle );
    //const double E_hat_f_calc        = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 29 );
    //std::ptrdiff_t ind_maximum_height;
    //objectiveHeight_Ascent_calc      = depVar_Height.maxCoeff( &ind_maximum_height );
    //objectiveAirspeed_Ascent_calc    = depVar_Airspeed( ind_maximum_height );
    //objectiveHeight_Ascent_calc      = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 11 );
    //objectiveAirspeed_Descent_calc   = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 13 );

    std::ptrdiff_t index_MinimumCentralTargetAngularDistanceToGo;
    std::ptrdiff_t index_MaximumNormalizedSpecificEnergy_Ascent;
    std::ptrdiff_t index_MinimumNormalizedSpecificEnergy_Ascent;
    std::ptrdiff_t index_MaximumNormalizedSpecificEnergy_Descent;
    std::ptrdiff_t index_MinimumNormalizedSpecificEnergy_Descent;
    std::ptrdiff_t index_MaximumAirspeed;
    std::ptrdiff_t index_maximumBodyFrameMechanicalLoad;
    std::ptrdiff_t index_MaximumDynamicPressure;
    std::ptrdiff_t index_MaximumBendingMoment;
    std::ptrdiff_t index_MaximumHeatFluxChapmanNose;
    std::ptrdiff_t index_MaximumIntegratedHeatLoad;
    std::ptrdiff_t index_MaximumBodyFrameTotalGLoad_x;
    std::ptrdiff_t index_MinimumBodyFrameTotalGLoad_x;
    std::ptrdiff_t index_MaximumPassengerFrameTotalGLoad_x;
    std::ptrdiff_t index_MinimumPassengerFrameTotalGLoad_x;
    std::ptrdiff_t index_MaximumPassengerFrameJerk_x;
    std::ptrdiff_t index_MinimumPassengerFrameJerk_x;
    std::ptrdiff_t index_MaximumBodyFrameTotalGLoad_z;
    std::ptrdiff_t index_MinimumBodyFrameTotalGLoad_z;
    std::ptrdiff_t index_MaximumPassengerFrameTotalGLoad_z;
    std::ptrdiff_t index_MinimumPassengerFrameTotalGLoad_z;
    std::ptrdiff_t index_MaximumPassengerFrameJerk_z;
    std::ptrdiff_t index_MinimumPassengerFrameJerk_z;

    const double minimum_CentralTargetAngularDistanceToGo = depVar_CentralTargetAngularDistanceToGo.minCoeff( &index_MinimumCentralTargetAngularDistanceToGo );

    double maximum_NormalizedSpecificEnergy_Ascent = TUDAT_NAN;// = depVar_NormalizedSpecificEnergy_Ascent.maxCoeff( &index_MaximumNormalizedSpecificEnergy_Ascent );
    double maximum_NormalizedSpecificEnergy_Descent = TUDAT_NAN;// = depVar_NormalizedSpecificEnergy_Descent.maxCoeff( &index_MaximumNormalizedSpecificEnergy_Descent );
    double minimum_NormalizedSpecificEnergy_Ascent = TUDAT_NAN;// = depVar_NormalizedSpecificEnergy_Ascent.maxCoeff( &index_MaximumNormalizedSpecificEnergy_Ascent );
    double minimum_NormalizedSpecificEnergy_Descent = TUDAT_NAN;// = depVar_NormalizedSpecificEnergy_Descent.maxCoeff( &index_MaximumNormalizedSpecificEnergy_Descent );

    if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_is_present ) )
    {
        maximum_NormalizedSpecificEnergy_Ascent = depVar_NormalizedSpecificEnergy_Ascent.maxCoeff( &index_MaximumNormalizedSpecificEnergy_Ascent );
        minimum_NormalizedSpecificEnergy_Ascent = depVar_NormalizedSpecificEnergy_Ascent.minCoeff( &index_MinimumNormalizedSpecificEnergy_Ascent );
    }

    if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_is_present ) )
    {
        maximum_NormalizedSpecificEnergy_Descent = depVar_NormalizedSpecificEnergy_Descent.maxCoeff( &index_MaximumNormalizedSpecificEnergy_Descent );
        minimum_NormalizedSpecificEnergy_Descent = depVar_NormalizedSpecificEnergy_Descent.minCoeff( &index_MinimumNormalizedSpecificEnergy_Descent );
    }

    const double maximumAirspeed_Ascent            = depVar_Airspeed.head( rowsAscent ).maxCoeff( &index_MaximumAirspeed );
    const double maximumAirspeed                   = depVar_Airspeed.maxCoeff( &index_MaximumAirspeed );
    const double maximumBodyFrameMechanicalLoad    = depVar_BodyFrameMechanicalGLoadMag.maxCoeff( &index_maximumBodyFrameMechanicalLoad );
    //const double maximumHeadingError              = depVar_HeadingToDynamicTargetError.maxCoeff( &index_MaximumHeadingError );
    const double maximumDynamicPressure            = depVar_DynamicPressure.maxCoeff( &index_MaximumDynamicPressure );
    const double maximumBendingMoment              = depVar_BendingMoment.maxCoeff( &index_MaximumBendingMoment );
    const double maximumHeatFluxChapmanNose        = depVar_HeatFluxChapmanNose.maxCoeff( &index_MaximumHeatFluxChapmanNose );
    const double maximumIntegratedHeatLoad         = depVar_IntegratedHeatLoad.maxCoeff( &index_MaximumIntegratedHeatLoad );
    const double maximumBodyFrameTotalGLoad_x      = depVar_BodyFrameTotalGLoad_x.maxCoeff( &index_MaximumBodyFrameTotalGLoad_x );
    const double minimumBodyFrameTotalGLoad_x      = depVar_BodyFrameTotalGLoad_x.minCoeff( &index_MinimumBodyFrameTotalGLoad_x );
    const double maximumBodyFrameTotalGLoad_z      = depVar_BodyFrameTotalGLoad_z.maxCoeff( &index_MaximumBodyFrameTotalGLoad_z );
    const double minimumBodyFrameTotalGLoad_z      = depVar_BodyFrameTotalGLoad_z.minCoeff( &index_MinimumBodyFrameTotalGLoad_z );
    const double maximumPassengerFrameTotalGLoad_x = depVar_PassengerFrameTotalGLoad_x.maxCoeff( &index_MaximumPassengerFrameTotalGLoad_x );
    const double minimumPassengerFrameTotalGLoad_x = depVar_PassengerFrameTotalGLoad_x.minCoeff( &index_MinimumPassengerFrameTotalGLoad_x );
    const double maximumPassengerFrameTotalGLoad_z = depVar_PassengerFrameTotalGLoad_z.maxCoeff( &index_MaximumPassengerFrameTotalGLoad_z );
    const double minimumPassengerFrameTotalGLoad_z = depVar_PassengerFrameTotalGLoad_z.minCoeff( &index_MinimumPassengerFrameTotalGLoad_z );
    const double maximumPassengerFrameJerk_x       = depVar_PassengerFrame_Jerk_x.segment( 1 , rowsTotal - 2 ).maxCoeff( &index_MaximumPassengerFrameJerk_x );
    const double minimumPassengerFrameJerk_x       = depVar_PassengerFrame_Jerk_x.segment( 1 , rowsTotal - 2 ).minCoeff( &index_MinimumPassengerFrameJerk_x );
    const double maximumPassengerFrameJerk_z       = depVar_PassengerFrame_Jerk_z.segment( 1 , rowsTotal - 2 ).maxCoeff( &index_MaximumPassengerFrameJerk_z );
    const double minimumPassengerFrameJerk_z       = depVar_PassengerFrame_Jerk_z.segment( 1 , rowsTotal - 2 ).minCoeff( &index_MinimumPassengerFrameJerk_z );

    Eigen::VectorXd extremeValues( 24 );
    extremeValues << minimum_CentralTargetAngularDistanceToGo,
            maximum_NormalizedSpecificEnergy_Ascent,
            minimum_NormalizedSpecificEnergy_Ascent,
            maximum_NormalizedSpecificEnergy_Descent,
            minimum_NormalizedSpecificEnergy_Descent,
            maximumAirspeed_Ascent,
            maximumAirspeed,
            maximumBodyFrameMechanicalLoad,
            maximumDynamicPressure,
            maximumBendingMoment,
            maximumHeatFluxChapmanNose,
            maximumIntegratedHeatLoad,
            maximumBodyFrameTotalGLoad_x,
            minimumBodyFrameTotalGLoad_x,
            maximumBodyFrameTotalGLoad_z,
            minimumBodyFrameTotalGLoad_z,
            maximumPassengerFrameTotalGLoad_x,
            minimumPassengerFrameTotalGLoad_x,
            maximumPassengerFrameTotalGLoad_z,
            minimumPassengerFrameTotalGLoad_z,
            maximumPassengerFrameJerk_x,
            minimumPassengerFrameJerk_x,
            maximumPassengerFrameJerk_z,
            minimumPassengerFrameJerk_z;


    std::vector< bislip::trajectory_phases::trajectory_phase_list > trajectoryPhaseList = bislip::trajectory_phases::getTrajectoryPhaseList();
    std::vector< bislip::parameters::nodal_parameter > nodalParametersList = bislip::parameters::getNodalParameters();
    std::vector< bislip::interpolators::interpolator_list > nodalInterpolatorList = bislip::interpolators::getNodalInterpolatorList();
    std::vector< bislip::parameters::node > nodeList = bislip::parameters::getNodeList();
    std::map< bislip::trajectory_phases::trajectory_phase_list,
            std::map< bislip::parameters::nodal_parameter,
            std::map< bislip::parameters::node, Eigen::Vector2d > > > nodalParameterBoundsPerTrajectoryPhaseMap =
            problemInput_->getNodalParameterBoundsPerTrajectoryPhaseMap( );
    
    
    const double originalInitialMass                   = bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::initial_mass ); // kg
    const double originalLandingMass                   = bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::landing_mass ); // kg
    const double actualLandingMass                     = bodyMap_.at( vehicleName )->getVehicleSystems()->getDryMass(); // kg
    const double initialMass_Ascent                    = depVar_CurrentMass[ 0 ];
    const double finalMass_Ascent                      = depVar_CurrentMass[ rowsAscent - 1 ];
    const double initialMass_Descent                   = finalMass_Ascent;
    const double finalMass_Descent                     = depVar_CurrentMass[ rowsTotal - 1 ];
    const double initialAirspeed_Ascent                = depVar_Airspeed[ 0 ];
    const double finalAirspeed_Ascent                  = depVar_Airspeed[ rowsAscent - 1 ];
    const double initialAirspeed_Descent               = finalAirspeed_Ascent;
    const double finalAirspeed_Descent                 = depVar_Airspeed[ rowsTotal - 1 ];

    const double timeOfFlight_Ascent                   = depVar_TimeOfFlight[ rowsAscent - 1 ];
    const double finalNormalizedSpecificEnergy_Descent = depVar_NormalizedSpecificEnergy[ rowsTotal - 1 ];
    const double finalHeight                           = depVar_Height[ rowsTotal - 1 ];
    const double finalBodyFrameMechanicalGLoadMag      = depVar_BodyFrameMechanicalGLoadMag[ rowsTotal - 1 ];
    const double timeOfFlightTotal                     = depVar_TimeOfFlight[ rowsTotal - 1 ];
    const double timeOfFlight_Descent                  = timeOfFlightTotal - timeOfFlight_Ascent;

    const double constraint_FinalHeight            = problemInput_->getConstraintValue( bislip::constraints::constraint_list::minimum_height_allowable );
    const double constraint_AscentMechanicalLoad   = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_mechanical_load_ascent );
    const double constraint_DescentMechanicalLoad  = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_mechanical_load_descent );
    const double constraint_ChapmanHeatFlux        = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_heat_flux );
    const double constraint_DynamicPressure        = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_dynamic_pressure );
    const double constraint_PitchMomentCoefficient = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_absolute_pitch_moment_coefficient );
    const double constraint_BendingMoment          = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_bending_moment );
    const double constraint_PassengerPosZLoad      = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_passenger_positive_z_load );
    const double constraint_PassengerNegZLoad      = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_passenger_negative_z_load );
    const double constraint_PassengerJerk          = problemInput_->getConstraintValue( bislip::constraints::constraint_list::maximum_passenger_jerk_load );


    Eigen::VectorXd constraints( 10 );
    constraints << constraint_FinalHeight,
            constraint_AscentMechanicalLoad,
            constraint_DescentMechanicalLoad,
            constraint_ChapmanHeatFlux,
            constraint_DynamicPressure,
            constraint_PitchMomentCoefficient,
            constraint_BendingMoment,
            constraint_PassengerPosZLoad,
            constraint_PassengerNegZLoad,
            constraint_PassengerJerk;

    //debugInfo = 1;

    if( debugInfo == 1 ){ std::cout << "Penalty: Angular Distance To Go" << std::endl; }
    const double angularDistanceToGoRatio = depVar_CentralTargetAngularDistanceToGo( rowsTotal - 1 ) / depVar_CentralTargetAngularDistanceToGo( 0 );
    const double costAngularDistanceToGo = 100 * ( angularDistanceToGoRatio );

    if( debugInfo == 1  ){ std::cout << "Extract Objective Function Case" << std::endl; }
    char objectiveFunctionCase = problemInput_->getObjectiveFunctionCase();

    std::map< bislip::fitness_vector::cost, double > costMap;
    costMap[ bislip::fitness_vector::cost::angular_distance_to_go ] = costAngularDistanceToGo;
    std::map< bislip::fitness_vector::penalty, double > penaltyMapAscent;
    std::map< bislip::fitness_vector::penalty, double > penaltyMapDescent;


    if( debugInfo == 1  ){ std::cout << "Select and Populate Fitness Vector" << std::endl; }

    if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_with_descent ) )
    {
        switch ( objectiveFunctionCase )
        {
        case 'A' :
        {
            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = 0.0;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = 0.0;
            break;
        }
        case 'B' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            break;
        }
        case 'C' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }


            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;

            break;
        }
        case 'D' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;

            break;
        }
        case 'E' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            break;
        }
        case 'F' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyHeatFluxChapmanAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyHeatFluxChapmanDescent;
            break;
        }
        case 'G' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateDescent;
            break;
        }
        case 'H' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;


            break;
        }
        case 'I' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            double logTerm = finalMass_Ascent / initialMass_Ascent;
            if( logTerm < 1 ) { logTerm = 1; }
            double velocityDifference = std::abs( initialAirspeed_Ascent - maximumAirspeed_Ascent );
            if( velocityDifference == 0.0 ) { velocityDifference = 1e-10; }
            const double costFuelMass = ( initialMass_Ascent - originalInitialMass ) / originalInitialMass + bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::specific_impulse ) * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * std::log( logTerm ) / velocityDifference;
            if( debugInfo == 1  ){ std::cout << "     costFuelMassAscent: " << costFuelMass << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;
            costMap[ bislip::fitness_vector::cost::fuel_mass ] = costFuelMass;

            break;
        }
        case 'J' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;

            break;
        }
        case 'K' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;

            break;
        }
        case 'L' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleAscent + penaltyMonotonicEnergyStateAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleDescent + penaltyMonotonicEnergyStateDescent;

            break;
        }
        case 'M' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleAscent + penaltyMonotonicEnergyStateAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleDescent + penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;

            break;
        }
        case 'N' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleAscent + penaltyMonotonicEnergyStateAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleDescent + penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            double logTerm = finalMass_Ascent / initialMass_Ascent;
            if( logTerm < 1 ) { logTerm = 1; }
            double velocityDifference = std::abs( initialAirspeed_Ascent - maximumAirspeed_Ascent );
            if( velocityDifference == 0.0 ) { velocityDifference = 1e-10; }
            const double costFuelMass = ( initialMass_Ascent - originalInitialMass ) / originalInitialMass + bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::specific_impulse ) * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * std::log( logTerm ) / velocityDifference;
            if( debugInfo == 1  ){ std::cout << "     costFuelMassAscent: " << costFuelMass << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;
            costMap[ bislip::fitness_vector::cost::fuel_mass ] = costFuelMass;

            break;
        }
        }
    }
    else if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_without_descent ) )
    {
        switch ( objectiveFunctionCase )
        {
        case 'A' :
        {
            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = 0.0;
            break;
        }
        case 'B' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            break;
        }
        case 'C' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;

            break;
        }
        case 'D' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;

            break;
        }
        case 'E' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            break;
        }
        case 'F' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyHeatFluxChapmanAscent;
            break;
        }
        case 'G' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }


            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateAscent;
            break;
        }
        case 'H' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateAscent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;


            break;
        }
        case 'I' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateAscent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            double logTerm = finalMass_Ascent / initialMass_Ascent;
            if( logTerm < 1 ) { logTerm = 1; }
            double velocityDifference = std::abs( initialAirspeed_Ascent - maximumAirspeed_Ascent );
            if( velocityDifference == 0.0 ) { velocityDifference = 1e-10; }
            const double costFuelMass = ( initialMass_Ascent - originalInitialMass ) / originalInitialMass + bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::specific_impulse ) * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * std::log( logTerm ) / velocityDifference;
            if( debugInfo == 1  ){ std::cout << "     costFuelMassAscent: " << costFuelMass << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;
            costMap[ bislip::fitness_vector::cost::fuel_mass ] = costFuelMass;

            break;
        }
        case 'J' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;

            break;
        }
        case 'K' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleAscent;

            break;
        }
        case 'L' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleAscent + penaltyMonotonicEnergyStateAscent;

            break;
        }
        case 'M' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleAscent + penaltyMonotonicEnergyStateAscent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;

            break;
        }
        case 'N' :
        {
            const double penaltyMechanicalLoadAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.head( rowsAscent ), constraint_AscentMechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_AscentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadAscent:  " <<  penaltyMechanicalLoadAscent << std::endl; }

            const double penaltyDynamicPressureAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl; }

            const double penaltyBendingMomentAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl; }

            const double penaltyHeatFluxChapmanAscent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl; }

            const double penaltyFlightPathAngleAscent = bislip::penalty_functions::computeConstraintViolationPenalty( -depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleAscent: " << penaltyFlightPathAngleAscent << std::endl; }

            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadAscent + penaltyDynamicPressureAscent + penaltyBendingMomentAscent + penaltyHeatFluxChapmanAscent;
            penaltyMapAscent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleAscent + penaltyMonotonicEnergyStateAscent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            double logTerm = finalMass_Ascent / initialMass_Ascent;
            if( logTerm < 1 ) { logTerm = 1; }
            double velocityDifference = std::abs( initialAirspeed_Ascent - maximumAirspeed_Ascent );
            if( velocityDifference == 0.0 ) { velocityDifference = 1e-10; }
            const double costFuelMass = ( initialMass_Ascent - originalInitialMass ) / originalInitialMass + bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::specific_impulse ) * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * std::log( logTerm ) / velocityDifference;
            if( debugInfo == 1  ){ std::cout << "     costFuelMassAscent: " << costFuelMass << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;
            costMap[ bislip::fitness_vector::cost::fuel_mass ] = costFuelMass;

            break;
        }
        }
    }
    else if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::descent_without_ascent ) )
    {
        switch ( objectiveFunctionCase )
        {
        case 'A' :
        {
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = 0.0;
            break;
        }
        case 'B' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            break;
        }
        case 'C' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;

            break;
        }
        case 'D' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;

            break;
        }
        case 'E' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            break;
        }
        case 'F' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyHeatFluxChapmanDescent;
            break;
        }
        case 'G' :
        {
            const double penaltyMonotonicEnergyStateAscent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateAscent:  " <<  penaltyMonotonicEnergyStateAscent << std::endl; }

            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapAscent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateAscent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateDescent;
            break;
        }
        case 'H' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;


            break;
        }
        case 'I' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::mechanical_load ] = penaltyMechanicalLoadDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::dynamic_pressure ] = penaltyDynamicPressureDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::bending_moment ] = penaltyBendingMomentDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::heat_flux ] = penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::monotonic_energy ] = penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            double logTerm = finalMass_Ascent / initialMass_Ascent;
            if( logTerm < 1 ) { logTerm = 1; }
            double velocityDifference = std::abs( initialAirspeed_Ascent - maximumAirspeed_Ascent );
            if( velocityDifference == 0.0 ) { velocityDifference = 1e-10; }
            const double costFuelMass = ( initialMass_Ascent - originalInitialMass ) / originalInitialMass + bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::specific_impulse ) * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * std::log( logTerm ) / velocityDifference;
            if( debugInfo == 1  ){ std::cout << "     costFuelMassAscent: " << costFuelMass << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;
            costMap[ bislip::fitness_vector::cost::fuel_mass ] = costFuelMass;

            break;
        }
        case 'J' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;

            break;
        }
        case 'K' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::flight_path_angle ] = penaltyFlightPathAngleDescent;

            break;
        }
        case 'L' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleDescent + penaltyMonotonicEnergyStateDescent;

            break;
        }
        case 'M' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleDescent + penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;

            break;
        }
        case 'N' :
        {
            const double penaltyMechanicalLoadDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BodyFrameMechanicalGLoadMag.tail( rowsTotal - rowsAscent ), constraint_DescentMechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_DescentMechanicalLoad );
            if( debugInfo == 1  ){ std::cout << "penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl; }

            const double penaltyDynamicPressureDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
            if( debugInfo == 1  ){ std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl; }

            const double penaltyBendingMomentDescent =  bislip::penalty_functions::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
            if( debugInfo == 1  ){ std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl; }

            const double penaltyHeatFluxChapmanDescent = bislip::penalty_functions::computeCompoundViolationPenalty( depVar_HeatFluxChapmanNose.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl; }

            const double penaltyFlightPathAngleDescent = bislip::penalty_functions::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );
            if( debugInfo == 1  ){ std::cout << "penaltyFlightPathAngleDescent: " << penaltyFlightPathAngleDescent << std::endl; }

            const double penaltyMonotonicEnergyStateDescent = bislip::penalty_functions::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );
            if( debugInfo == 1  ){ std::cout << "penaltyMonotonicEnergyStateDescent: " << penaltyMonotonicEnergyStateDescent << std::endl; }

            penaltyMapDescent[ bislip::fitness_vector::penalty::kinetics ] = penaltyMechanicalLoadDescent + penaltyDynamicPressureDescent + penaltyBendingMomentDescent + penaltyHeatFluxChapmanDescent;
            penaltyMapDescent[ bislip::fitness_vector::penalty::kinematics ] = penaltyFlightPathAngleDescent + penaltyMonotonicEnergyStateDescent;

            const double costHeatLoadChapman = bislip::math_tools::computeSumOfEigenVectorXd( depVar_HeatFluxChapmanNose.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
            if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapman: " << costHeatLoadChapman << std::endl; }

            double logTerm = finalMass_Ascent / initialMass_Ascent;
            if( logTerm < 1 ) { logTerm = 1; }
            double velocityDifference = std::abs( initialAirspeed_Ascent - maximumAirspeed_Ascent );
            if( velocityDifference == 0.0 ) { velocityDifference = 1e-10; }
            const double costFuelMass = ( initialMass_Ascent - originalInitialMass ) / originalInitialMass + bislipSystems->getVehicleParameter( bislip::vehicle_parameters::vehicle_parameter_list::specific_impulse ) * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * std::log( logTerm ) / velocityDifference;
            if( debugInfo == 1  ){ std::cout << "     costFuelMassAscent: " << costFuelMass << std::endl; }

            costMap[ bislip::fitness_vector::cost::heat_load ] = costHeatLoadChapman;
            costMap[ bislip::fitness_vector::cost::fuel_mass ] = costFuelMass;

            break;
        }
        }
    }


    //! Assign values to Fitness vector! At the moment these are all 'objective
    //! functions'. To modify this I have change the header file and define how
    //! many elements are objective functions, equality contraints, and
    //! inequality constraints. This vector here must contain them is that exact order: nOF, nEC, nIC.
    std::vector< double > fitness;

    bislip::utilities::convertMapToVector( costMap, fitness );
    bislip::utilities::convertMapToVector( penaltyMapAscent, fitness );
    bislip::utilities::convertMapToVector( penaltyMapDescent, fitness );

    if( objectiveFunctionCase != 'A' )
    {
        if( problemInput_->getTrajectoryTypeToSimulate( bislip::trajectory_phases::trajectory_type_to_simulate::ascent_without_descent ) == false )
        {
            if( debugInfo == 1 ){ std::cout << "Add '1.0'  to Fitness vector elements and multiply by multiple of 10^n if (n+1)*10% > Angular Distance To Go >= n*10%" << std::endl; }
            if( debugInfo == 1 ){ std::cout << "    This is an attempt to avoid having these individuals end up at the top simply becuase most elements";
                std::cout << "     in the fitness vector will inevitably be close  or equal to zero (0)" << std::endl; }
            if( angularDistanceToGoRatio >= 0.1 )
            {
                double multiplier = 1.0;
                if( angularDistanceToGoRatio >= 0.9 ) { multiplier = 1000000000; }
                else if( angularDistanceToGoRatio >= 0.8 ) { multiplier = 100000000; }
                else if( angularDistanceToGoRatio >= 0.7 ) { multiplier = 10000000; }
                else if( angularDistanceToGoRatio >= 0.6 ) { multiplier = 1000000; }
                else if( angularDistanceToGoRatio >= 0.5 ) { multiplier = 100000; }
                else if( angularDistanceToGoRatio >= 0.4 ) { multiplier = 10000; }
                else if( angularDistanceToGoRatio >= 0.3 ) { multiplier = 1000; }
                else if( angularDistanceToGoRatio >= 0.2 ) { multiplier = 100; }
                else if( angularDistanceToGoRatio >= 0.1 ) { multiplier = 10; }

                std::transform( std::next( fitness.begin() ), fitness.end(), std::next( fitness.begin() ), std::bind( std::plus< double >(), std::placeholders::_1, 1.0 ) );
                std::transform( std::next( fitness.begin() ), fitness.end(), std::next( fitness.begin() ), std::bind( std::multiplies< double >(), std::placeholders::_1, multiplier ) );
            }

            if( debugInfo == 1 ){ std::cout << "Add '1.0'  to Fitness vector elements and multiply by '10.0' if Final Height >= 120% the value of the constraint" << std::endl; }
            if( debugInfo == 1 ){ std::cout << "    This is an attempt to avoid giving preference to individuals that may reach the destination but are";
                std::cout << "    too high" << std::endl; }
            if( finalHeight / constraint_FinalHeight >= 1.2 )
            {
                std::transform( std::next( fitness.begin() ), fitness.end(), std::next( fitness.begin() ), std::bind( std::plus< double >(), std::placeholders::_1, 1.0 ) );
                std::transform( std::next( fitness.begin() ), fitness.end(), std::next( fitness.begin() ), std::bind( std::multiplies< double >(), std::placeholders::_1, 10.0 ) );
            }

        }
        else
        {
            if( debugInfo == 1 ){ std::cout << "Add '1.0'  to Fitness vector elements and multiply by multiple of 10^n if (n+1)*10% > Angular Distance To Go >= n*10%" << std::endl; }
            if( debugInfo == 1 ){ std::cout << "    This is an attempt to avoid having these individuals end up at the top simply becuase most elements";
                std::cout << "     in the fitness vector will inevitably be close  or equal to zero (0)" << std::endl; }
            if( angularDistanceToGoRatio >= 0.7 )
            {
                double multiplier = 1.0;
                if( angularDistanceToGoRatio >= 0.9 ) { multiplier = 1000; }
                else if( angularDistanceToGoRatio >= 0.8 ) { multiplier = 100; }
                else if( angularDistanceToGoRatio >= 0.7 ) { multiplier = 10; }

                std::transform( std::next( fitness.begin() ), fitness.end(), std::next( fitness.begin() ), std::bind( std::plus< double >(), std::placeholders::_1, 1.0 ) );
                std::transform( std::next( fitness.begin() ), fitness.end(), std::next( fitness.begin() ), std::bind( std::multiplies< double >(), std::placeholders::_1, multiplier ) );
            }

        }
    }

    if( bislipSystems->getFlag( bislip::flags::flag_list::evolution_evaluation ) )
    {

        if( debugInfo == 1 ){ std::cout << "Printing to terminal" << std::endl; }
        unsigned int generationNumber = problemInput_->getGenerationNumber( );
        unsigned int individualNumber = problemInput_->getIndividualNumber( );

        std::cout  << generationNumber << " - " << individualNumber << " | " ;
        for( unsigned int i = 0; i < fitness.size() - 1; i++)
        {
            std::cout  << fitness[ i ] << " | " ;
        }
        //std::cout  << fitness[ fitness.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof << "  ||  " << " Theoretical delV = " << theoreticaldelV_Ascent << "  ||  " << " Goal delV = " << goaldelV_Ascent << "  ||  " << " Actual delV = " << actualdelV_Ascent << std::endl;
        //    std::cout  << fitness[ fitness.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof <<  std::endl;
        //std::cout  << fitness[ fitness.size() - 1 ] << " ||  " << timeOfFlight << " ||  " << saveTrajectoryOutput << " ||  " << rowsAscent << std::endl;
        std::cout << fitness[ fitness.size() - 1 ] << " ||  " << timeOfFlightTotal  << std::endl;


        //! Get time stamp for this specific simulation. This avoids overwriting the
        //! file if another individual with the same properties shows up in other
        //! evolutions.
        std::chrono::time_point< std::chrono::system_clock > printTime = bislip::utilities::getDateTime( );

        std::string printTimeString = bislip::utilities::convertDateTimeToString( false, printTime );

        unsigned int millisSincePlayTime = static_cast< unsigned int >( std::chrono::duration_cast< std::chrono::milliseconds >( printTime - ( problemInput_->getPlayTimePair() ).first ).count() );

        //! Create unique filename that cannot be overwritten due to the timestamp.
        std::string simulation_file_name_suffix =
                ( problemInput_->getPlayTimePair() ).second + "_" +
                printTimeString + "_" +
                std::to_string( millisSincePlayTime ) + "_" +
                std::to_string( timeOfFlightTotal ) + "_" +
                std::to_string( finalHeight ) + "_" +
                std::to_string( targetLat_deg_calc ) + "_" +
                std::to_string( targetLon_deg_calc ) + "_" +
                std::to_string( costAngularDistanceToGo ) + "_" +
                std::to_string( maximumBodyFrameMechanicalLoad ) + "_" +
                std::to_string( finalBodyFrameMechanicalGLoadMag );

        //! Convert decision vector from std::vector< double > to Eigen::VectorXd
        const Eigen::VectorXd x_Eigen = Eigen::Map< const Eigen::VectorXd, Eigen::Unaligned >( x.data(), static_cast< unsigned int >( x.size() ) );

        //! Create new Eigen::VectorXd to include flag that the individual was printed or not
        Eigen::VectorXd x_EigenFlagged( x_Eigen.size() + 5u );
        x_EigenFlagged << 0, 0, 0, generationNumber, individualNumber, x_Eigen;

        problemInput_->setPopulationIndividual( simulation_file_name_suffix, x_EigenFlagged );
        problemInput_->setPopulationHistoryIndividual( simulation_file_name_suffix, x_EigenFlagged );

        //! Convert fitness vector from std::vector< double > to Eigen::VectorXd
        const Eigen::VectorXd fitness_Eigen = Eigen::Map< const Eigen::VectorXd, Eigen::Unaligned >( fitness.data(), static_cast< unsigned int >( fitness.size() ) );

        //! Create new Eigen::VectorXd to include flag that the individual was printed or not
        Eigen::VectorXd fitness_EigenFlagged( fitness_Eigen.size() + 5u );
        fitness_EigenFlagged << 0, 0, 0, generationNumber, individualNumber, fitness_Eigen;

        problemInput_->setFitnessIndividual( simulation_file_name_suffix, fitness_EigenFlagged );
        problemInput_->setFitnessHistoryIndividual( simulation_file_name_suffix, fitness_EigenFlagged );


        //! Create new Eigen::VectorXd to include flag that the individual was printed or not
        Eigen::VectorXd extremesAndConstraints( extremeValues.size() + constraints.size() + 5u );
        extremesAndConstraints << 0, 0, 0, generationNumber, individualNumber, extremeValues, constraints;
        problemInput_->setExtremesAndConstraintsIndividual( simulation_file_name_suffix, extremesAndConstraints );
        problemInput_->setExtremesAndConstraintsHistoryIndividual( simulation_file_name_suffix, extremesAndConstraints );

        if( x != x_mod )
        {
             if( bislipSystems->getFlag( bislip::flags::flag_list::evolution_evaluation ) )
             {
                 problemInput_->setIndividualsToReplace( problemInput_->getIndividualNumber(), x_mod, fitness );
             }
        }
        problemInput_->setIndividualNumber( individualNumber + 1 );

    }
    return fitness;
}// Fitness function.

} //namespace bislip
