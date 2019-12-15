function [ evolutions ] = initializeExtremesAndConstraintsFields( evolutions, k, i )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


    extremesAndConstraintsList = readcell( 'extremesAndConstraintsList.dat','Delimiter',{','} );

  for j = 1:size(extremesAndConstraintsList,1)
      evolutions(k).population(i).extremesAndConstraints.(extremesAndConstraintsList{j,1}).figureSaveNameContent = extremesAndConstraintsList{j,1};
      evolutions(k).population(i).extremesAndConstraints.(extremesAndConstraintsList{j,1}).variableLabel = extremesAndConstraintsList{j,2};
      evolutions(k).population(i).extremesAndConstraints.(extremesAndConstraintsList{j,1}).units = extremesAndConstraintsList{j,3};
      evolutions(k).population(i).extremesAndConstraints.(extremesAndConstraintsList{j,1}).scalingFactor = extremesAndConstraintsList{j,4};
      evolutions(k).population(i).extremesAndConstraints.(extremesAndConstraintsList{j,1}).limits = [extremesAndConstraintsList{j,5} extremesAndConstraintsList{j,6}];
      evolutions(k).population(i).extremesAndConstraints.(extremesAndConstraintsList{j,1}).tick = extremesAndConstraintsList{j,7};
      evolutions(k).population(i).extremesAndConstraints.(extremesAndConstraintsList{j,1}).value = nan;
  end

%{
evolutions(k).population(i).extremesAndConstraints.minimum_CentralTargetAngularDistanceToGo.variableLabel = 'Minimum Central Target Angular Distance To Go';
evolutions(k).population(i).extremesAndConstraints.minimum_CentralTargetAngularDistanceToGo.figureSaveNameContent = 'minimumCentralTargetAngularDistanceToGo';
evolutions(k).population(i).extremesAndConstraints.minimum_CentralTargetAngularDistanceToGo.units = '(deg)';
evolutions(k).population(i).extremesAndConstraints.minimum_CentralTargetAngularDistanceToGo.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimum_CentralTargetAngularDistanceToGo.limits = [0 60];
evolutions(k).population(i).extremesAndConstraints.minimum_CentralTargetAngularDistanceToGo.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Ascent.variableLabel = 'Maximum Normalized Specific Energy - Ascent';
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Ascent.figureSaveNameContent = 'maximumNormalizedSpecificEnergyAscent';
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Ascent.units = '(-)';
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Ascent.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Ascent.limits = [0 2];
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Ascent.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Ascent.variableLabel = 'Minimum Normalized Specific Energy - Ascent';
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Ascent.figureSaveNameContent = 'minimumNormalizedSpecificEnergyAscent';
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Ascent.units = '(-)';
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Ascent.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Ascent.limits = [0 2];
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Ascent.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Descent.variableLabel =  'Maximum Normalized Specific Energy - Descent';
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Descent.figureSaveNameContent = 'maximumNormalizedSpecificEnergyDescent';
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Descent.units = '(-)';
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Descent.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Descent.limits = [0 2];
evolutions(k).population(i).extremesAndConstraints.maximum_NormalizedSpecificEnergy_Descent.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Descent.variableLabel =  'Minimum Normalized Specific Energy - Descent';
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Descent.figureSaveNameContent = 'minimumNormalizedSpecificEnergyDescent';
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Descent.units = '(-)';
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Descent.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Descent.limits = [0 2];
evolutions(k).population(i).extremesAndConstraints.minimum_NormalizedSpecificEnergy_Descent.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumAirspeed_Ascent.variableLabel = 'Maximum Airspeed - Ascent';
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed_Ascent.figureSaveNameContent = 'maximumAirspeedAscent';
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed_Ascent.units = '(m/s)';
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed_Ascent.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed_Ascent.limits = [0 8000];
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed_Ascent.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumAirspeed.variableLabel = 'Maximum Airspeed';
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed.figureSaveNameContent = 'maximumAirspeed';
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed.units = '(m/s)';
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed.limits = [0 8000];
evolutions(k).population(i).extremesAndConstraints.maximumAirspeed.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameMechanicalLoad.variableLabel = 'Maximum Body Frame Mechanical Load';
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameMechanicalLoad.figureSaveNameContent = 'maximumBodyFrameMechanicalLoad';
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameMechanicalLoad.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameMechanicalLoad.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameMechanicalLoad.limits = [0 5];
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameMechanicalLoad.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumDynamicPressure.variableLabel = 'Maximum Dynamic Pressure';
evolutions(k).population(i).extremesAndConstraints.maximumDynamicPressure.figureSaveNameContent = 'maximumDynamicPressure';
evolutions(k).population(i).extremesAndConstraints.maximumDynamicPressure.units = '(kPa)';
evolutions(k).population(i).extremesAndConstraints.maximumDynamicPressure.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumDynamicPressure.limits = [0 25000];
evolutions(k).population(i).extremesAndConstraints.maximumDynamicPressure.scalingFactor = 1e3;

evolutions(k).population(i).extremesAndConstraints.maximumBendingMoment.variableLabel = 'Maximum Bending Moment';
evolutions(k).population(i).extremesAndConstraints.maximumBendingMoment.figureSaveNameContent = 'maximumBendingMoment';
evolutions(k).population(i).extremesAndConstraints.maximumBendingMoment.units = '(Pa-rad)';
evolutions(k).population(i).extremesAndConstraints.maximumBendingMoment.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumBendingMoment.limits = [0 5000];
evolutions(k).population(i).extremesAndConstraints.maximumBendingMoment.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumHeatFluxChapmanNose.variableLabel = 'Maximum Chapman Heat Flux - Nose';
evolutions(k).population(i).extremesAndConstraints.maximumHeatFluxChapmanNose.figureSaveNameContent = 'maximumChapmanHeatFlux';
evolutions(k).population(i).extremesAndConstraints.maximumHeatFluxChapmanNose.units = '(kW/m^2)';
evolutions(k).population(i).extremesAndConstraints.maximumHeatFluxChapmanNose.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumHeatFluxChapmanNose.limits = [0 700];
evolutions(k).population(i).extremesAndConstraints.maximumHeatFluxChapmanNose.scalingFactor = 1e3;

evolutions(k).population(i).extremesAndConstraints.maximumIntegratedHeatLoad.variableLabel = 'Maximum Integrated Heat Flux';
evolutions(k).population(i).extremesAndConstraints.maximumIntegratedHeatLoad.figureSaveNameContent = 'maximumIntegratedHeatFlux';
evolutions(k).population(i).extremesAndConstraints.maximumIntegratedHeatLoad.units = '(MJ/m^2)';
evolutions(k).population(i).extremesAndConstraints.maximumIntegratedHeatLoad.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumIntegratedHeatLoad.limits = [0 500];
evolutions(k).population(i).extremesAndConstraints.maximumIntegratedHeatLoad.scalingFactor = 1e6;

evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_x.variableLabel = strcat('Maximum Body Frame Total g-load - x');
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_x.figureSaveNameContent = 'maximumBodyFrameTotalGLoad_x';
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_x.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_x.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_x.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_x.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_x.variableLabel = strcat('Minimum Body Frame Total g-load - x');
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_x.figureSaveNameContent = 'minimumBodyFrameTotalGLoad_x';
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_x.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_x.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_x.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_x.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_z.variableLabel = strcat('Maximum Body Frame Total g-load - z');
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_z.figureSaveNameContent = 'maximumBodyFrameTotalGLoad_z';
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_z.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_z.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_z.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.maximumBodyFrameTotalGLoad_z.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_z.variableLabel = strcat('Minimum Body Frame Total g-load - z');
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_z.figureSaveNameContent = 'minimumBodyFrameTotalGLoad_z';
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_z.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_z.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_z.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.minimumBodyFrameTotalGLoad_z.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_x.variableLabel = strcat('Maximum Passenger Frame Total g-load - x');
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_x.figureSaveNameContent = 'maximumPassengerFrameTotalGLoad_x';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_x.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_x.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_x.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_x.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_x.variableLabel = strcat('Minimum Passenger Frame Total g-load - x');
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_x.figureSaveNameContent = 'minimumPassengerFrameTotalGLoad_x';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_x.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_x.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_x.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_x.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_z.variableLabel = strcat('Maximum Passenger Frame Total g-load - z');
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_z.figureSaveNameContent = 'maximumPassengerFrameTotalGLoad_z';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_z.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_z.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_z.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameTotalGLoad_z.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_z.variableLabel = strcat('Minimum Passenger Frame Total g-load - z');
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_z.figureSaveNameContent = 'minimumPassengerFrameTotalGLoad_z';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_z.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_z.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_z.limits = [-3 3];
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameTotalGLoad_z.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_x.variableLabel = strcat('Maximum Passenger Frame Jerk - x');
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_x.figureSaveNameContent = 'maximumPassengerJerk_x';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_x.units = '(m/s^3)';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_x.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_x.limits = [-5 5];
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_x.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_x.variableLabel = strcat('Minimum Passenger Frame Jerk - x');
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_x.figureSaveNameContent = 'maximumPassengerJerk_x';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_x.units = '(m/s^3)';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_x.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_x.limits = [-5 5];
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_x.scalingFactor = 1;


evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_z.variableLabel = strcat('Maximum Passenger Frame Jerk - z');
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_z.figureSaveNameContent = 'maximumPassengerJerk_z';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_z.units = '(m/s^3)';
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_z.value = nan;
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_z.limits = [-5 5];
evolutions(k).population(i).extremesAndConstraints.maximumPassengerFrameJerk_z.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_z.variableLabel = strcat('Minimum Passenger Frame Jerk - z');
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_z.figureSaveNameContent = 'minimumPassengerJerk_z';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_z.units = '(m/s^3)';
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_z.value = nan;
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_z.limits = [-5 5];
evolutions(k).population(i).extremesAndConstraints.minimumPassengerFrameJerk_z.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_FinalHeight.variableLabel = 'Final Height';
evolutions(k).population(i).extremesAndConstraints.constraint_FinalHeight.figureSaveNameContent = 'Final Height';
evolutions(k).population(i).extremesAndConstraints.constraint_FinalHeight.units = 'km';
evolutions(k).population(i).extremesAndConstraints.constraint_FinalHeight.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_FinalHeight.limits = [0 200];
evolutions(k).population(i).extremesAndConstraints.constraint_FinalHeight.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_AscentMechanicalLoad.variableLabel = 'Constraint: Mechanical Load - Ascent';
evolutions(k).population(i).extremesAndConstraints.constraint_AscentMechanicalLoad.figureSaveNameContent = 'constraintMechanicalLoadAscent';
evolutions(k).population(i).extremesAndConstraints.constraint_AscentMechanicalLoad.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.constraint_AscentMechanicalLoad.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_AscentMechanicalLoad.limits = [0 5];
evolutions(k).population(i).extremesAndConstraints.constraint_AscentMechanicalLoad.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_DescentMechanicalLoad.variableLabel = 'Constraint: Mechanical Load - Descent';
evolutions(k).population(i).extremesAndConstraints.constraint_DescentMechanicalLoad.figureSaveNameContent = 'constraintMechanicalLoadDescent';
evolutions(k).population(i).extremesAndConstraints.constraint_DescentMechanicalLoad.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.constraint_DescentMechanicalLoad.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_DescentMechanicalLoad.limits = [0 5];
evolutions(k).population(i).extremesAndConstraints.constraint_DescentMechanicalLoad.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_ChapmanHeatFlux.variableLabel = 'Constraint: Chapman Heat Flux';
evolutions(k).population(i).extremesAndConstraints.constraint_ChapmanHeatFlux.figureSaveNameContent = 'constraintChapmanHeatFlux';
evolutions(k).population(i).extremesAndConstraints.constraint_ChapmanHeatFlux.units = '(kW/m^2)';
evolutions(k).population(i).extremesAndConstraints.constraint_ChapmanHeatFlux.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_ChapmanHeatFlux.limits = [0 700];
evolutions(k).population(i).extremesAndConstraints.constraint_ChapmanHeatFlux.scalingFactor = 1e3;

evolutions(k).population(i).extremesAndConstraints.constraint_DynamicPressure.variableLabel = 'Constraint: Dynamic Pressure';
evolutions(k).population(i).extremesAndConstraints.constraint_DynamicPressure.figureSaveNameContent = 'constraintDynamicPressure';
evolutions(k).population(i).extremesAndConstraints.constraint_DynamicPressure.units = '(kPa)';
evolutions(k).population(i).extremesAndConstraints.constraint_DynamicPressure.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_DynamicPressure.limits = [0 25000];
evolutions(k).population(i).extremesAndConstraints.constraint_DynamicPressure.scalingFactor = 1e3;

evolutions(k).population(i).extremesAndConstraints.constraint_PitchMomentCoefficient.variableLabel = 'Constraint: Pitch Moment Coefficient';
evolutions(k).population(i).extremesAndConstraints.constraint_PitchMomentCoefficient.figureSaveNameContent = 'constraintPitchMomentCoefficient';
evolutions(k).population(i).extremesAndConstraints.constraint_PitchMomentCoefficient.units = '(-)';
evolutions(k).population(i).extremesAndConstraints.constraint_PitchMomentCoefficient.vale = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_PitchMomentCoefficient.limits = [-1 1];
evolutions(k).population(i).extremesAndConstraints.constraint_PitchMomentCoefficient.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_BendingMoment.variableLabel = 'Constraint: Bending Moment';
evolutions(k).population(i).extremesAndConstraints.constraint_BendingMoment.figureSaveNameContent = 'constraintBendingMoment';
evolutions(k).population(i).extremesAndConstraints.constraint_BendingMoment.units = '(Pa-rad)';
evolutions(k).population(i).extremesAndConstraints.constraint_BendingMoment.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_BendingMoment.limits = [0 5000];
evolutions(k).population(i).extremesAndConstraints.constraint_BendingMoment.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_PassengerPosZLoad.variableLabel = strcat('Constraint: Passenger Frame Positive Load - z');
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerPosZLoad.figureSaveNameContent = 'constraintPassengerFramePositiveLoad_z';
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerPosZLoad.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerPosZLoad.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerPosZLoad.limits = [0 5];
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerPosZLoad.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_PassengerNegZLoad.variableLabel = strcat('Constraint: Passenger Frame Negative Load - z');
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerNegZLoad.figureSaveNameContent = 'constraintPassengerFrameNegativeLoad_z';
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerNegZLoad.units = '(g_0)';
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerNegZLoad.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerNegZLoad.limits = [-5 0];
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerNegZLoad.scalingFactor = 1;

evolutions(k).population(i).extremesAndConstraints.constraint_PassengerJerk.variableLabel = 'Constraint: Passenger Frame Jerk';
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerJerk.figureSaveNameContent = 'constraintPassengerFrameJerk';
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerJerk.units = '(m/s^3)';
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerJerk.value = nan;
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerJerk.limits = [0 5];
evolutions(k).population(i).extremesAndConstraints.constraint_PassengerJerk.scalingFactor = 1;
%}


end

