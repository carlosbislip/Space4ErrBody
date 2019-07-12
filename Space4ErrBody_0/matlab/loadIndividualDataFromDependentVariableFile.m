function [ dependentVariableTimeHistory ] = loadIndividualDataFromDependentVariableFile( dependentVariableNames, depvar, interp_Ascent, interp_Descent, DV_mapped_Ascent, DV_mapped_Descent, headingErrorDeadbandBounds, alphaMachBounds )

%% Dependent Variables


        
        % Loop through each relevant type of parameter
        for i = 1:length(dependentVariableNames)
            dependentVariableTimeHistory.(dependentVariableNames{i}) = depvar(:,i);
        end
dependentVariableTimeHistory.acc_x                                = depvar(:,17) + depvar(:,20) + depvar(:,23);
dependentVariableTimeHistory.acc_y                                = depvar(:,18) + depvar(:,21) + depvar(:,24);
dependentVariableTimeHistory.acc_z                                = depvar(:,19) + depvar(:,22) + depvar(:,25);
dependentVariableTimeHistory.acc_aero_M                           = sqrt(depvar(:,18).^2 + depvar(:,19).^2 + depvar(:,20).^2);
dependentVariableTimeHistory.acc_grav_M                           = sqrt(depvar(:,21).^2 + depvar(:,22).^2 + depvar(:,23).^2);
dependentVariableTimeHistory.acc_thru_M                           = sqrt(depvar(:,24).^2 + depvar(:,25).^2 + depvar(:,26).^2);



% 
% individual.dependentVariableTimeHistory.t                                    = depvar(:,1);
% individual.dependentVariableTimeHistory.x_R                                  = depvar(:,2);
% individual.dependentVariableTimeHistory.y_R                                  = depvar(:,3);
% individual.dependentVariableTimeHistory.z_R                                  = depvar(:,4);
% individual.dependentVariableTimeHistory.altitude                             = depvar(:,5);
% individual.dependentVariableTimeHistory.latitude_angle                       = depvar(:,6);
% individual.dependentVariableTimeHistory.longitude_angle                      = depvar(:,7);
% individual.dependentVariableTimeHistory.heading_angle                        = depvar(:,8);
% individual.dependentVariableTimeHistory.flight_path_angle                    = depvar(:,9);
% individual.dependentVariableTimeHistory.angle_of_attack                      = depvar(:,10);
% individual.dependentVariableTimeHistory.angle_of_sideslip                    = depvar(:,11);
% individual.dependentVariableTimeHistory.bank_angle                           = depvar(:,12);
% individual.dependentVariableTimeHistory.height                               = depvar(:,13);
% individual.dependentVariableTimeHistory.mach                                 = depvar(:,14);
% individual.dependentVariableTimeHistory.airspeed                             = depvar(:,15);
% individual.dependentVariableTimeHistory.total_aero_g_load                    = depvar(:,16);
% individual.dependentVariableTimeHistory.acc_aero_x                           = depvar(:,17);
% individual.dependentVariableTimeHistory.acc_aero_y                           = depvar(:,18);
% individual.dependentVariableTimeHistory.acc_aero_z                           = depvar(:,19);
% individual.dependentVariableTimeHistory.acc_grav_x                           = depvar(:,21);
% individual.dependentVariableTimeHistory.acc_grav_y                           = depvar(:,21);
% individual.dependentVariableTimeHistory.acc_grav_z                           = depvar(:,22);
% individual.dependentVariableTimeHistory.acc_thru_x                           = depvar(:,23);
% individual.dependentVariableTimeHistory.acc_thru_y                           = depvar(:,24);
% individual.dependentVariableTimeHistory.acc_thru_z                           = depvar(:,25);
% individual.dependentVariableTimeHistory.dynamicPressure                      = depvar(:,26);
% individual.dependentVariableTimeHistory.heat_rate_TUDAT_nose                 = depvar(:,27);
% individual.dependentVariableTimeHistory.mass_rate                            = depvar(:,28);
% individual.dependentVariableTimeHistory.mass                                 = depvar(:,29);
% individual.dependentVariableTimeHistory.specificEnergy                       = depvar(:,30);
% individual.dependentVariableTimeHistory.normalizedSpecificEnergy             = depvar(:,31);
% individual.dependentVariableTimeHistory.evaluated_throttle_setting           = depvar(:,32);
% individual.dependentVariableTimeHistory.evaluated_thrust_elevation_angle     = depvar(:,33);
% individual.dependentVariableTimeHistory.evaluated_thrust_azimuth_angle       = depvar(:,34);
% individual.dependentVariableTimeHistory.evaluated_angle_of_attack            = depvar(:,35);
% individual.dependentVariableTimeHistory.evaluated_bank_angle                 = depvar(:,36);
% individual.dependentVariableTimeHistory.engine_status                        = depvar(:,37);
% individual.dependentVariableTimeHistory.distance_traveled                    = depvar(:,38);
% individual.dependentVariableTimeHistory.distance_to_go                       = depvar(:,39);
% individual.dependentVariableTimeHistory.heading_to_target                    = depvar(:,40);
% individual.dependentVariableTimeHistory.heading_error                        = depvar(:,41);
% individual.dependentVariableTimeHistory.heat_flux_tauber_leadingedge         = depvar(:,42);
% individual.dependentVariableTimeHistory.body_fixed_thrust_load_x             = depvar(:,43);
% individual.dependentVariableTimeHistory.body_fixed_thrust_load_y             = depvar(:,44);
% individual.dependentVariableTimeHistory.body_fixed_thrust_load_z             = depvar(:,45);
% individual.dependentVariableTimeHistory.bending_moment                       = depvar(:,46);
% individual.dependentVariableTimeHistory.local_gravity_1                      = depvar(:,47);
% individual.dependentVariableTimeHistory.local_gravity_2                      = depvar(:,48);
% individual.dependentVariableTimeHistory.local_gravity_3                      = depvar(:,49);
% individual.dependentVariableTimeHistory.skip_suppression_limit               = depvar(:,50);
% individual.dependentVariableTimeHistory.bodyflap_deflection                  = depvar(:,51);
% individual.dependentVariableTimeHistory.increment_Cm_bodyflap                = depvar(:,52);
% individual.dependentVariableTimeHistory.increment_Cm_bodyflap_dif            = depvar(:,53);
% individual.dependentVariableTimeHistory.body_fixed_total_load_x              = depvar(:,54);
% individual.dependentVariableTimeHistory.body_fixed_total_load_y              = depvar(:,55);
% individual.dependentVariableTimeHistory.body_fixed_total_load_z              = depvar(:,56);
% individual.dependentVariableTimeHistory.body_fixed_total_g_load_x            = depvar(:,57);
% individual.dependentVariableTimeHistory.body_fixed_total_g_load_y            = depvar(:,58);
% individual.dependentVariableTimeHistory.body_fixed_total_g_load_z            = depvar(:,59);
% individual.dependentVariableTimeHistory.body_fixed_total_g_load_mag          = depvar(:,60);
% individual.dependentVariableTimeHistory.body_fixed_aero_load_x               = depvar(:,61);
% individual.dependentVariableTimeHistory.body_fixed_aero_load_y               = depvar(:,62);
% individual.dependentVariableTimeHistory.body_fixed_aero_load_z               = depvar(:,63);
% individual.dependentVariableTimeHistory.aero_force_coefficient_C_D           = depvar(:,64);
% individual.dependentVariableTimeHistory.aero_force_coefficient_C_S           = depvar(:,65);
% individual.dependentVariableTimeHistory.aero_force_coefficient_C_L           = depvar(:,66);
% individual.dependentVariableTimeHistory.aero_moment_coefficient_C_l          = depvar(:,67);
% individual.dependentVariableTimeHistory.aero_moment_coefficient_C_m          = depvar(:,68);
% individual.dependentVariableTimeHistory.aero_moment_coefficient_C_n          = depvar(:,69);
% individual.dependentVariableTimeHistory.heat_flux_chapman_nose               = depvar(:,70);
% individual.dependentVariableTimeHistory.localDensity                         = depvar(:,71);
% individual.dependentVariableTimeHistory.passenger_fixed_total_g_load_x       = depvar(:,72);
% individual.dependentVariableTimeHistory.passenger_fixed_total_g_load_y       = depvar(:,73);
% individual.dependentVariableTimeHistory.passenger_fixed_total_g_load_z       = depvar(:,74);
% individual.dependentVariableTimeHistory.commanded_throttle_setting           = depvar(:,75);
% individual.dependentVariableTimeHistory.commanded_thrust_elevation_angle     = depvar(:,76);
% individual.dependentVariableTimeHistory.commanded_thrust_azimuth_angle       = depvar(:,77);
% individual.dependentVariableTimeHistory.commanded_angle_of_attack            = depvar(:,78);
% individual.dependentVariableTimeHistory.commanded_bank_angle                 = depvar(:,79);
% individual.dependentVariableTimeHistory.currentLiftForce                     = depvar(:,80);
% individual.dependentVariableTimeHistory.headingErrorDeadband                 = depvar(:,81);
% individual.dependentVariableTimeHistory.tempBankAngle                        = depvar(:,82);
% individual.dependentVariableTimeHistory.reversalConditional                  = depvar(:,83);
% individual.dependentVariableTimeHistory.bank_angle_reversal_trigger          = depvar(:,84);
% individual.dependentVariableTimeHistory.wall_temperature_chapman             = depvar(:,85);
% individual.dependentVariableTimeHistory.wall_temperature_tauber_stagnation   = depvar(:,86);
% individual.dependentVariableTimeHistory.wall_temperature_tauber_flatplate    = depvar(:,87);
% individual.dependentVariableTimeHistory.heat_flux_tauber_stagnation          = depvar(:,88);
% individual.dependentVariableTimeHistory.heat_flux_tauber_flatplate           = depvar(:,89);
% individual.dependentVariableTimeHistory.groundtrack_covered                  = depvar(:,90);
% individual.dependentVariableTimeHistory.groundtrack_difference               = depvar(:,91);
% individual.dependentVariableTimeHistory.time_vector                          = depvar(:,92);
% individual.dependentVariableTimeHistory.flight_path_angle_rate               = depvar(:,93);
% individual.dependentVariableTimeHistory.cumulative_distance_travelled        = depvar(:,94);
% individual.dependentVariableTimeHistory.elevon_deflection                    = depvar(:,95);
% individual.dependentVariableTimeHistory.thrustMagnitude                      = depvar(:,96);
% individual.dependentVariableTimeHistory.speedOfSound                         = depvar(:,97);
% individual.dependentVariableTimeHistory.adiabaticWallTemperature             = depvar(:,98);
% individual.dependentVariableTimeHistory.freestreamTemperature                = depvar(:,99);
% individual.dependentVariableTimeHistory.currentDragForce                     = depvar(:,100);
% individual.dependentVariableTimeHistory.estimatedFightPathAngle              = depvar(:,101);
% individual.dependentVariableTimeHistory.aerodyamicFrameAerodynamicLoad_x     = depvar(:,102);
% individual.dependentVariableTimeHistory.aerodyamicFrameAerodynamicLoad_y     = depvar(:,103);
% individual.dependentVariableTimeHistory.aerodyamicFrameAerodynamicLoad_z     = depvar(:,104);
% individual.dependentVariableTimeHistory.aerodyamicFrameTotalLoad_x           = depvar(:,105);
% individual.dependentVariableTimeHistory.aerodyamicFrameTotalLoad_y           = depvar(:,106);
% individual.dependentVariableTimeHistory.aerodyamicFrameTotalLoad_z           = depvar(:,107);
% individual.dependentVariableTimeHistory.aerodyamicFrameTotalAcceleration_x   = depvar(:,108);
% individual.dependentVariableTimeHistory.aerodyamicFrameTotalAcceleration_y   = depvar(:,109);
% individual.dependentVariableTimeHistory.aerodyamicFrameTotalAcceleration_z   = depvar(:,110);
% 
% individual.dependentVariableTimeHistory.passengerFrameTotalLoad_x            = depvar(:,111);
% individual.dependentVariableTimeHistory.passengerFrameTotalLoad_y            = depvar(:,112);
% individual.dependentVariableTimeHistory.passengerFrameTotalLoad_z            = depvar(:,113);
% individual.dependentVariableTimeHistory.passengerFrameTotalAcceleration_x    = depvar(:,114);
% individual.dependentVariableTimeHistory.passengerFrameTotalAcceleration_y    = depvar(:,115);
% individual.dependentVariableTimeHistory.passengerFrameTotalAcceleration_z    = depvar(:,116);
% 
% individual.dependentVariableTimeHistory.passengerFrameJerk_x                 = depvar(:,117);
% individual.dependentVariableTimeHistory.passengerFrameJerk_y                 = depvar(:,118);
% individual.dependentVariableTimeHistory.passengerFrameJerk_z                 = depvar(:,119);
% 
% individual.dependentVariableTimeHistory.trajectoryphase                      = depvar(:,120);
% individual.dependentVariableTimeHistory.passengerFrameJerk_x_calc            = depvar(:,121);
% individual.dependentVariableTimeHistory.passengerFrameJerk_y_calc            = depvar(:,122);
% individual.dependentVariableTimeHistory.passengerFrameJerk_z_calc            = depvar(:,123);


individual.dependentVariableTimeHistory.acc_x                                = depvar(:,17) + depvar(:,20) + depvar(:,23);
individual.dependentVariableTimeHistory.acc_y                                = depvar(:,18) + depvar(:,21) + depvar(:,24);
individual.dependentVariableTimeHistory.acc_z                                = depvar(:,19) + depvar(:,22) + depvar(:,25);
individual.dependentVariableTimeHistory.acc_aero_M                           = sqrt(depvar(:,18).^2 + depvar(:,19).^2 + depvar(:,20).^2);
individual.dependentVariableTimeHistory.acc_grav_M                           = sqrt(depvar(:,21).^2 + depvar(:,22).^2 + depvar(:,23).^2);
individual.dependentVariableTimeHistory.acc_thru_M                           = sqrt(depvar(:,24).^2 + depvar(:,25).^2 + depvar(:,26).^2);


individual.DV_E_mapped_Ascent               = DV_mapped_Ascent(:,1);
individual.DV_angle_of_attack_Ascent        = DV_mapped_Ascent(:,2);
individual.DV_bank_angle_Ascent             = DV_mapped_Ascent(:,3);
individual.DV_thrust_elevation_angle_Ascent = DV_mapped_Ascent(:,4);
individual.DV_thrust_azimuth_angle_Ascent   = DV_mapped_Ascent(:,5);
individual.DV_throttle_setting_Ascent       = DV_mapped_Ascent(:,6);
individual.DV_node_location_Ascent          = DV_mapped_Ascent(:,7);

individual.DV_E_mapped_Descent               = DV_mapped_Descent(:,1);
individual.DV_angle_of_attack_Descent        = DV_mapped_Descent(:,2);
individual.DV_bank_angle_Descent             = DV_mapped_Descent(:,3);
individual.DV_thrust_elevation_angle_Descent = DV_mapped_Descent(:,4);
individual.DV_thrust_azimuth_angle_Descent   = DV_mapped_Descent(:,5);
individual.DV_throttle_setting_Descent       = DV_mapped_Descent(:,6);
individual.DV_node_location_Descent          = DV_mapped_Descent(:,7);


individual.interpolators.Ascent.InputData.normalizedSpecificEnergy = DV_mapped_Ascent(:,1);
individual.interpolators.Ascent.InputData.angleOfAttack        = DV_mapped_Ascent(:,2);
individual.interpolators.Ascent.InputData.bankAngle            = DV_mapped_Ascent(:,3);
individual.interpolators.Ascent.InputData.thrustElevationAngle = DV_mapped_Ascent(:,4);
individual.interpolators.Ascent.InputData.thrustAzimuthAngle   = DV_mapped_Ascent(:,5);
individual.interpolators.Ascent.InputData.throttleSetting      = DV_mapped_Ascent(:,6);
individual.interpolators.Ascent.InputData.nodeLocation         = DV_mapped_Ascent(:,7);

individual.interpolators.Descent.InputData.normalizedSpecificEnergy = DV_mapped_Descent(:,1);
individual.interpolators.Descent.InputData.angleOfAttack        = DV_mapped_Descent(:,2);
individual.interpolators.Descent.InputData.bankAngle            = DV_mapped_Descent(:,3);
individual.interpolators.Descent.InputData.thrustElevationAngle = DV_mapped_Descent(:,4);
individual.interpolators.Descent.InputData.thrustAzimuthAngle   = DV_mapped_Descent(:,5);
individual.interpolators.Descent.InputData.throttleSetting      = DV_mapped_Descent(:,6);
individual.interpolators.Descent.InputData.nodeLocation         = DV_mapped_Descent(:,7);

individual.interpolators.Ascent.Evaluation.normalizedSpecificEnergy = interp_Ascent(:,1);
individual.interpolators.Ascent.Evaluation.angleOfAttack        = interp_Ascent(:,2);
individual.interpolators.Ascent.Evaluation.bankAngle            = interp_Ascent(:,3);
individual.interpolators.Ascent.Evaluation.thrustElevationAngle = interp_Ascent(:,4);
individual.interpolators.Ascent.Evaluation.thrustAzimuthAngle   = interp_Ascent(:,5);
individual.interpolators.Ascent.Evaluation.throttleSetting      = interp_Ascent(:,6);

individual.interpolators.Descent.Evaluation.normalizedSpecificEnergy = interp_Descent(:,1);
individual.interpolators.Descent.Evaluation.angleOfAttack        = interp_Descent(:,2);
individual.interpolators.Descent.Evaluation.bankAngle            = interp_Descent(:,3);
individual.interpolators.Descent.Evaluation.thrustElevationAngle = interp_Descent(:,4);
individual.interpolators.Descent.Evaluation.thrustAzimuthAngle   = interp_Descent(:,5);
individual.interpolators.Descent.Evaluation.throttleSetting      = interp_Descent(:,6);
%% Interpolators
individual.interp_E_mapped_Ascent               = interp_Ascent(:,1);
individual.interp_angle_of_attack_Ascent        = interp_Ascent(:,2);
individual.interp_bank_angle_Ascent             = interp_Ascent(:,3);
individual.interp_thrust_elevation_angle_Ascent = interp_Ascent(:,4);
individual.interp_thrust_azimuth_angle_Ascent   = interp_Ascent(:,5);
individual.interp_throttle_setting_Ascent       = interp_Ascent(:,6);

individual.interpolators.Ascent.Evaluation.normalizedSpecificEnergy = interp_Ascent(:,1);
individual.interpolators.Ascent.Evaluation.angleOfAttack        = interp_Ascent(:,2);
individual.interpolators.Ascent.Evaluation.bankAngle            = interp_Ascent(:,3);
individual.interpolators.Ascent.Evaluation.thrustElevationAngle = interp_Ascent(:,4);
individual.interpolators.Ascent.Evaluation.thrustAzimuthAngle   = interp_Ascent(:,5);
individual.interpolators.Ascent.Evaluation.throttleSetting      = interp_Ascent(:,6);


individual.interp_E_mapped_Descent               = interp_Descent(:,1);
individual.interp_angle_of_attack_Descent        = interp_Descent(:,2);
individual.interp_bank_angle_Descent             = interp_Descent(:,3);
individual.interp_thrust_elevation_angle_Descent = interp_Descent(:,4);
individual.interp_thrust_azimuth_angle_Descent   = interp_Descent(:,5);
individual.interp_throttle_setting_Descent       = interp_Descent(:,6);


individual.interpolators.Ascent.Evaluation.normalizedSpecificEnergy = interp_Ascent(:,1);
individual.interpolators.Ascent.Evaluation.angleOfAttack        = interp_Ascent(:,2);
individual.interpolators.Ascent.Evaluation.bankAngle            = interp_Ascent(:,3);
individual.interpolators.Ascent.Evaluation.thrustElevationAngle = interp_Ascent(:,4);
individual.interpolators.Ascent.Evaluation.thrustAzimuthAngle   = interp_Ascent(:,5);
individual.interpolators.Ascent.Evaluation.throttleSetting      = interp_Ascent(:,6);

individual.interpolators.Descent.Evaluation.normalizedSpecificEnergy = interp_Descent(:,1);
individual.interpolators.Descent.Evaluation.angleOfAttack        = interp_Descent(:,2);
individual.interpolators.Descent.Evaluation.bankAngle            = interp_Descent(:,3);
individual.interpolators.Descent.Evaluation.thrustElevationAngle = interp_Descent(:,4);
individual.interpolators.Descent.Evaluation.thrustAzimuthAngle   = interp_Descent(:,5);
individual.interpolators.Descent.Evaluation.throttleSetting      = interp_Descent(:,6);



individual.alphaMachBounds_Mach        = alphaMachBounds(:,1);
individual.alphaMachBounds_LB          = alphaMachBounds(:,2);
individual.alphaMachBounds_UB          = alphaMachBounds(:,3);

individual.headingErrorDeadBand_distance    = headingErrorDeadbandBounds(:,1);
individual.headingErrorDeadBand_LB          = headingErrorDeadbandBounds(:,3);
individual.headingErrorDeadBand_UP          = headingErrorDeadbandBounds(:,2);


end