function [evolutions,output] = Load_Data(option,p,objectiveFunctionCase,Folder_Path_List,npop_files)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
disp('LOAD_DATA')


%% Load lists of analyzed simulations

switch option
    
    case 1
        
        % Figure out output name based on what is in folder path
        if contains(convertCharsToStrings(Folder_Path_List{p}),'OUTPUT') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'T_'),'.mat'));
        elseif contains(convertCharsToStrings(Folder_Path_List{p}),'VALIDATION') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'N_'),'.mat'));
        end
        
        % Open file if it exists. Create variable otherwise.
        if exist(output,'file') == 2
            load(output)
        else
            disp('Run again with option 2 or 3. No .mat file found.');
            
            %             k = npop_files + 1;
            %             evolutions(k).evolution(npop_files)    = nan;
            %             evolutions(k).v_i(npop_files)          = nan;
            %             evolutions(k).gamma_i(npop_files)      = nan;
            %             evolutions(k).chi_i(npop_files)        = nan;
            %             evolutions(k).population.trajectory(npop_files) = nan;
            %             evolutions(k).fitness(npop_files)      = nan;
            %             evolutions(k).lat(npop_files)          = nan;
            %             evolutions(k).lon(npop_files)          = nan;
            %             evolutions(k).tof(npop_files)          = nan;
        end
        
        idx = [];
        
    case {2, 3}
        
        % Figure out output name based on what is in folder path
        if contains(convertCharsToStrings(Folder_Path_List),'OUTPUT') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List,'T_'),'.mat'));
        elseif contains(convertCharsToStrings(Folder_Path_List),'VALIDATION') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List,'N_'),'.mat'));
        end
        
        % Open file if it exists. Create variable otherwise.
        %if exist(output,'file') == 2
        %    load(output)
        % else
        %if contains(pop_path, 'monteCarlo') == true
        %k = npop_files;
        %else
        %   k = npop_files + 1;
        %end
        for k = 1: npop_files
            evolutions(k).evolution                                          = nan;
            
            evolutions(k).population.size.collective                         = nan;
            evolutions(k).population.size.printed                            = nan;
            evolutions(k).population.size.nonPrinted                         = nan;
            evolutions(k).population.size.fitness                            = nan;
            evolutions(k).population.indices.printed                         = nan;
            evolutions(k).population.indices.nonPrinted                      = nan;
            evolutions(k).population.indices.nonDominatedFront               = nan;
            evolutions(k).population.indices.trajectoryPhaseChange           = nan;
            
            evolutions(k).population.decisionVector.collective               = nan;
            evolutions(k).population.decisionVector.printed                  = nan;
            evolutions(k).population.decisionVector.nonPrinted               = nan;
            
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.figureTitleContent             =        'Node Interval';
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.figureSaveNameContent          =        'ascentNodeInterval';
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.units                          =        '(-)';
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.one                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.two                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.three                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.four                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.five                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.six                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.seven                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.eight                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.nodeInterval.data.nine                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.figureTitleContent            =        'Angle of Attack';
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.figureSaveNameContent         =        'ascentAngleOfAttack';
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.units                         =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.one                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.two                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.three                    =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.four                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.five                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.six                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.seven                    =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.eight                    =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.nine                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.angleOfAttack.data.ten                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.figureTitleContent                =        'Bank Angle';
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.figureSaveNameContent             =        'ascentBankAngle';
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.units                             =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.one                          =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.two                          =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.three                        =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.four                         =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.five                         =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.six                          =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.seven                        =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.eight                        =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.nine                         =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.bankAngle.data.ten                          =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.figureTitleContent     =        'Thrust Elevation Angle';
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.figureSaveNameContent  =        'ascentThrustElevationAngle';
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.units                          =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.one                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.two                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.three                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.four                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.five                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.six                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.seven                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.eight                     =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.nine                      =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustElevationAngle.data.ten                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.figureTitleContent       =        'Thrust Azimuth Angle';
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.figureSaveNameContent    =        'ascentThrustAzimuthAngle';
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.units                    =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.one                         =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.two                         =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.three                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.four                        =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.five                        =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.six                         =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.seven                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.eight                       =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.nine                        =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.ten                         =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.figureTitleContent          =        'Throttle Setting';
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.figureSaveNameContent       =        'ascentThrottleSetting';
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.units                       =        '(-)';
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.one                            =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.two                            =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.three                          =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.four                           =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.five                           =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.six                            =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.seven                          =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.eight                          =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.nine                           =       nan;
            evolutions(k).population.decisionVector.parameters.Ascent.throttleSetting.data.ten                            =       nan;
            evolutions(k).population.decisionVector.parameters.Common.initialFlightPathAngle.figureTitleContent           =        'Initial Flight-Path Angle';
            evolutions(k).population.decisionVector.parameters.Common.initialFlightPathAngle.figureSaveNameContent        =        'initialFlightPathAngle';
            evolutions(k).population.decisionVector.parameters.Common.initialFlightPathAngle.units                        =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Common.initialFlightPathAngle.data.one                         =       nan;
            evolutions(k).population.decisionVector.parameters.Common.initialLaunchHeadingAngle.figureTitleContent        =        'Initial Launch Heading Angle';
            evolutions(k).population.decisionVector.parameters.Common.initialLaunchHeadingAngle.figureSaveNameContent     =        'initialLaunchHeadingAngle';
            evolutions(k).population.decisionVector.parameters.Common.initialLaunchHeadingAngle.units                     =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Common.initialLaunchHeadingAngle.data.one                      =       nan;
            evolutions(k).population.decisionVector.parameters.Common.initialVelocity.figureTitleContent                  =        'Initial Velocity';
            evolutions(k).population.decisionVector.parameters.Common.initialVelocity.figureSaveNameContent               =        'initialVelocity';
            evolutions(k).population.decisionVector.parameters.Common.initialVelocity.units                               =        '(m/s)';
            evolutions(k).population.decisionVector.parameters.Common.initialVelocity.data.one                                =       nan;
            evolutions(k).population.decisionVector.parameters.Common.maximumVelocity.figureTitleContent                  =        'Maximum Velocity';
            evolutions(k).population.decisionVector.parameters.Common.maximumVelocity.figureSaveNameContent               =        'maximumVelocity';
            evolutions(k).population.decisionVector.parameters.Common.maximumVelocity.units                               =        '(m/s)';
            evolutions(k).population.decisionVector.parameters.Common.maximumVelocity.data.one                                =       nan;
            evolutions(k).population.decisionVector.parameters.Common.maximumHeight.figureTitleContent                    =        'Maximum Height';
            evolutions(k).population.decisionVector.parameters.Common.maximumHeight.figureSaveNameContent                 =        'maximumHeight';
            evolutions(k).population.decisionVector.parameters.Common.maximumHeight.units                                 =        '(m)';
            evolutions(k).population.decisionVector.parameters.Common.maximumHeight.data.one                                  =       nan;
            evolutions(k).population.decisionVector.parameters.Common.additionalMass.figureTitleContent                   =        'Additional Mass';
            evolutions(k).population.decisionVector.parameters.Common.additionalMass.figureSaveNameContent                =        'additionalMass';
            evolutions(k).population.decisionVector.parameters.Common.additionalMass.units                                =        '(kg)';
            evolutions(k).population.decisionVector.parameters.Common.additionalMass.data.one                                 =       nan;
            evolutions(k).population.decisionVector.parameters.Common.terminationDistanceRatio.figureTitleContent         =        'Termination Distance Ratio';
            evolutions(k).population.decisionVector.parameters.Common.terminationDistanceRatio.figureSaveNameContent      =        'terminationDistanceRatio';
            evolutions(k).population.decisionVector.parameters.Common.terminationDistanceRatio.units                      =        '(-)';
            evolutions(k).population.decisionVector.parameters.Common.terminationDistanceRatio.data.one                       =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.figureTitleContent            =        'Node Interval';
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.figureSaveNameContent         =        'descentNodeInterval';
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.units                         =        '(-)';
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.one                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.two                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.three                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.four                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.five                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.six                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.seven                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.eight                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.nodeInterval.data.nine                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.figureTitleContent           =        'Angle of Attack';
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.figureSaveNameContent        =        'descentAngleOfAttack';
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.units                        =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.one                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.two                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.three                   =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.four                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.five                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.six                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.seven                   =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.eight                   =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.nine                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.angleOfAttack.data.ten                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.figureTitleContent       =        'Bank Angle';
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.figureSaveNameContent    =        'descentBankAngle';
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.units                    =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.one                         =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.two                         =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.three                       =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.four                        =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.five                        =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.six                         =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.seven                       =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.eight                       =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.nine                        =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.bankAngle.data.ten                         =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.figureTitleContent            =        'Thrust Elevation Angle';
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.figureSaveNameContent         =        'descentThrustElevationAngle';
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.units                         =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.one                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.two                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.three                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.four                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.five                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.six                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.seven                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.eight                    =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.nine                     =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustElevationAngle.data.ten                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.figureTitleContent      =        'Thrust Azimuth Angle';
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.figureSaveNameContent   =        'descentThrustAzimuthAngle';
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.units                   =        '(deg)';
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.one                        =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.two                        =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.three                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.four                       =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.five                       =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.six                        =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.seven                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.eight                      =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.nine                       =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.ten                        =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.figureTitleContent         =        'Throttle Setting';
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.figureSaveNameContent      =        'descentThrottleSetting';
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.units                      =        '(-)';
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.one                           =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.two                           =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.three                         =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.four                          =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.five                          =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.six                           =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.seven                         =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.eight                         =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.nine                          =       nan;
            evolutions(k).population.decisionVector.parameters.Descent.throttleSetting.data.ten                           =       nan;
            evolutions(k).population.decisionVector.parameters.Common.finalVelocity.figureTitleContent                    =        'Final Velocity';
            evolutions(k).population.decisionVector.parameters.Common.finalVelocity.figureSaveNameContent                 =        'finalVelocity';
            evolutions(k).population.decisionVector.parameters.Common.finalVelocity.units                                 =        '(m/s)';
            evolutions(k).population.decisionVector.parameters.Common.finalVelocity.data.one                                  =       nan;
            evolutions(k).population.decisionVector.parameters.Common.skipSuppressionTriggerTime.figureTitleContent       =        'Skip Suppression Trigger Time';
            evolutions(k).population.decisionVector.parameters.Common.skipSuppressionTriggerTime.figureSaveNameContent    =        'skipSuppressionTriggerTime';
            evolutions(k).population.decisionVector.parameters.Common.skipSuppressionTriggerTime.units                    =        '(s)';
            evolutions(k).population.decisionVector.parameters.Common.skipSuppressionTriggerTime.data.one                     =       nan;
            evolutions(k).population.decisionVector.parameters.Common.maximumMechanicalLoad.figureTitleContent            =        'Maximum Mechanical Load';
            evolutions(k).population.decisionVector.parameters.Common.maximumMechanicalLoad.figureSaveNameContent         =        'maximumMechanicalLoad';
            evolutions(k).population.decisionVector.parameters.Common.maximumMechanicalLoad.units                         =        '(g_0)';
            evolutions(k).population.decisionVector.parameters.Common.maximumMechanicalLoad.data.one                          =       nan;
            
            evolutions(k).population.dependentVariableTimeHistory.t                                    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.x_R                                  	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.y_R                                  	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.z_R                                  	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.altitude                             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.latitudeAngle                       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.longitudeAngle                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.headingAngle                        	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.flightPathAngle                    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.angleOfAttack                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.angleOfSideslip                    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bankAngle                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.height                               	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.machNumber                            =	nan;
            evolutions(k).population.dependentVariableTimeHistory.airspeed                             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodynamicFrameAerodynamicGLoad      =	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_aero_x                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_aero_y                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_aero_z                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_grav_x                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_grav_y                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_grav_z                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_thru_x                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_thru_y                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_thru_z                           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.dynamicPressure                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.heatFluxTUDATNose                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.massRate                            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.mass                                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.specificEnergy                       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.normalizedSpecificEnergy             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.evaluatedThrottleSetting           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.evaluatedThrustElevationAngle     	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.evaluatedThrustAzimuthAngle       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.evaluatedAngleOfAttack            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.evaluatedBankAngle                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.engine_status                        	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.angularDistanceTraveled               =	nan;
            evolutions(k).population.dependentVariableTimeHistory.angularDistanceToGo                   =	nan;
            evolutions(k).population.dependentVariableTimeHistory.headingToTarget                    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.headingError                        	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.heatFluxTauberLeadingEdge         	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameThrustLoad_x             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameThrustLoad_y             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameThrustLoad_z             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bendingMoment                       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.localGravity_1                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.localGravity_2                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.localGravity_3                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.skipSuppressionBankAngleLimit         =	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyflapDeflectionAngle               =	nan;
            evolutions(k).population.dependentVariableTimeHistory.increment_Cm_bodyflap                	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.increment_Cm_bodyflap_dif            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameTotalLoad_x              	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameTotalLoad_y              	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameTotalLoad_z              	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameTotalGLoad_x            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameTotalGLoad_y            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameTotalGLoad_z            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude          	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameAerodynamicLoad_x               	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameAerodynamicLoad_y               	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bodyFrameAerodynamicLoad_z               	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodynamicCoefficient_CD           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodynamicCoefficient_CS           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodynamicCoefficient_CL           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodynamicCoefficient_Cl          	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodynamicCoefficient_Cm          	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodynamicCoefficient_Cn          	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.heatFluxChapmanNose               	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.localDensity                         	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalGLoad_x       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalGLoad_y       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalGLoad_z       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.commandedThrottleSetting           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.commandedThrustElevationAngle     	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.commandedThrustAzimuthAngle       	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.commandedAngleOfAttack            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.commandedBankAngle                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.currentLiftForce                     	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.headingErrorDeadband                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.tempBankAngle                        	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.reversalConditional                  	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.bankAngleReversalTrigger          	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.equilibriumWallTemperatureChapmanNose             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.equilibriumWallTemperatureTauberStagnationLeadingEdge   	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.equilibriumWallTemperatureTauberFlatPlateLeadingEdge    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.heatFluxTauberStagnationLeadingEdge          	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.heatFluxTauberFlatPlateLeadingEdge           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.cumulativeAngularDistanceTravelled                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.groundtrackDifference               	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.timeOfFlight                          	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.flightPathAngleRate               	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.cumulativeCartesianDistanceTravelled        	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.elevonDeflectionAngle                    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.thrustMagnitude                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.speedOfSound                         	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.adiabaticWallTemperature             	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.freestreamTemperature                	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.currentDragForce                     	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.estimatedFightPathAngle              	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameAerodynamicLoad_x     	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameAerodynamicLoad_y     	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameAerodynamicLoad_z     	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameTotalLoad_x           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameTotalLoad_y           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameTotalLoad_z           	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameTotalAcceleration_x   	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameTotalAcceleration_y   	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.aerodyamicFrameTotalAcceleration_z   	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalLoad_x            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalLoad_y            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalLoad_z            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalAcceleration_x    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalAcceleration_y    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameTotalAcceleration_z    	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameJerk_x                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameJerk_y                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameJerk_z                 	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.trajectoryPhase                      	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.angularDistanceCoveredRatio        	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameJerk_x_calc            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameJerk_y_calc            	=	nan;
            evolutions(k).population.dependentVariableTimeHistory.passengerFrameJerk_z_calc            	=	nan;
            
            
            evolutions(k).population.dependentVariableTimeHistory.acc_x                            	    = nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_y                           	    = nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_z                           	    = nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_aero_M                      	    = nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_grav_M                          	= nan;
            evolutions(k).population.dependentVariableTimeHistory.acc_thru_M                       	    = nan;
            
            evolutions(k).population.interpolators.Ascent.InputData.normalizedSpecificEnergy 	=	nan;
            evolutions(k).population.interpolators.Ascent.InputData.angleOfAttack        	=	nan;
            evolutions(k).population.interpolators.Ascent.InputData.bankAngle            	=	nan;
            evolutions(k).population.interpolators.Ascent.InputData.thrustElevationAngle 	=	nan;
            evolutions(k).population.interpolators.Ascent.InputData.thrustAzimuthAngle   	=	nan;
            evolutions(k).population.interpolators.Ascent.InputData.throttleSetting      	=	nan;
            evolutions(k).population.interpolators.Ascent.InputData.nodeLocation         	=	nan;
            evolutions(k).population.interpolators.Ascent.Evaluation.normalizedSpecificEnergy 	=	nan;
            evolutions(k).population.interpolators.Ascent.Evaluation.angleOfAttack        	=	nan;
            evolutions(k).population.interpolators.Ascent.Evaluation.bankAngle            	=	nan;
            evolutions(k).population.interpolators.Ascent.Evaluation.thrustElevationAngle 	=	nan;
            evolutions(k).population.interpolators.Ascent.Evaluation.thrustAzimuthAngle   	=	nan;
            evolutions(k).population.interpolators.Ascent.Evaluation.throttleSetting      	=	nan;
            
            evolutions(k).population.interpolators.Descent.InputData.normalizedSpecificEnergy 	=	nan;
            evolutions(k).population.interpolators.Descent.InputData.angleOfAttack        	=	nan;
            evolutions(k).population.interpolators.Descent.InputData.bankAngle            	=	nan;
            evolutions(k).population.interpolators.Descent.InputData.thrustElevationAngle 	=	nan;
            evolutions(k).population.interpolators.Descent.InputData.thrustAzimuthAngle   	=	nan;
            evolutions(k).population.interpolators.Descent.InputData.throttleSetting      	=	nan;
            evolutions(k).population.interpolators.Descent.InputData.nodeLocation         	=	nan;
            evolutions(k).population.interpolators.Descent.Evaluation.normalizedSpecificEnergy 	=	nan;
            evolutions(k).population.interpolators.Descent.Evaluation.angleOfAttack        	=	nan;
            evolutions(k).population.interpolators.Descent.Evaluation.bankAngle            	=	nan;
            evolutions(k).population.interpolators.Descent.Evaluation.thrustElevationAngle 	=	nan;
            evolutions(k).population.interpolators.Descent.Evaluation.thrustAzimuthAngle   	=	nan;
            evolutions(k).population.interpolators.Descent.Evaluation.throttleSetting      	=	nan;
            
            if objectiveFunctionCase == 'A'
                evolutions(k).population.size.fitness = 3;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
            elseif objectiveFunctionCase == 'B'
                evolutions(k).population.size.fitness = 5;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energy   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.energy   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
            elseif objectiveFunctionCase == 'C'
                evolutions(k).population.size.fitness = 5;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.monotonicEnergy   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.monotonicEnergy   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
            elseif objectiveFunctionCase == 'D'
                evolutions(k).population.size.fitness = 5;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.monotonicEnergyAndFlightPath   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.monotonicEnergyAndFlightPath   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
            elseif objectiveFunctionCase == 'E'
                evolutions(k).population.size.fitness = 5;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.energyAndFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
            elseif objectiveFunctionCase == 'F'
                evolutions(k).population.size.fitness = 7;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.fuelMassCost   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.energyAndFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.fuelMassCost   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
            elseif objectiveFunctionCase == 'G'
                evolutions(k).population.size.fitness = 5;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
            elseif objectiveFunctionCase == 'H'
                evolutions(k).population.size.fitness = 7;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Ascent.passengerAccelerationsZ   = nan;
                evolutions(k).population.fitnessVector.Descent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.passengerAccelerationsZ   = nan;
            elseif objectiveFunctionCase == 'I'
                evolutions(k).population.size.fitness = 9;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Ascent.passengerAccelerationsZ   = nan;
                evolutions(k).population.fitnessVector.Ascent.passengerJerkXYZ   = nan;
                evolutions(k).population.fitnessVector.Descent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.passengerAccelerationsZ   = nan;
                evolutions(k).population.fitnessVector.Descent.passengerJerkXYZ   = nan;
            elseif objectiveFunctionCase == 'J'
                evolutions(k).population.size.fitness = 7;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Ascent.passengerAccelZJerkXYZ   = nan;
                evolutions(k).population.fitnessVector.Descent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.passengerAccelZJerkXYZ   = nan;
            elseif objectiveFunctionCase == 'K'
                evolutions(k).population.size.fitness = 5;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.allDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.energyMassFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.allDynamics   = nan;
            elseif objectiveFunctionCase == 'L'
                evolutions(k).population.size.fitness = 11;
                evolutions(k).population.fitnessVector.Common.angularDistanceToGo   = nan;
                evolutions(k).population.fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Ascent.fuelMassCost   = nan;
                evolutions(k).population.fitnessVector.Ascent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Ascent.passengerAccelerationsZ   = nan;
                evolutions(k).population.fitnessVector.Ascent.passengerJerkXYZ   = nan;
                evolutions(k).population.fitnessVector.Descent.energyAndFlightPathAngle   = nan;
                evolutions(k).population.fitnessVector.Descent.fuelMassCost   = nan;
                evolutions(k).population.fitnessVector.Descent.basicDynamics   = nan;
                evolutions(k).population.fitnessVector.Descent.passengerAccelerationsZ   = nan;
                evolutions(k).population.fitnessVector.Descent.passengerJerkXYZ   = nan;
            end
            
            evolutions(k).population.fitnessVector.Cost                             = nan;
            evolutions(k).population.fitnessVector.DominationSet                    = nan;
            evolutions(k).population.fitnessVector.DominatedCount                   = nan;
            evolutions(k).population.fitnessVector.Rank                             = nan;
            evolutions(k).population.fitnessVector.Common.magnitude                 = nan;
            evolutions(k).population.fitnessVector.Common.ranking                   = nan;
            evolutions(k).population.fitnessVector.Common.scalingFactor             = nan;
            
            evolutions(k).Common.Bounds.headingError.angularDistanceToGo  = nan;
            evolutions(k).Common.Bounds.headingError.upperBound  = nan;
            evolutions(k).Common.Bounds.headingError.lowerBound  = nan;
            evolutions(k).Common.Bounds.AngleOfAttack.machNumber  = nan;
            evolutions(k).Common.Bounds.AngleOfAttack.upperBound  = nan;
            evolutions(k).Common.Bounds.AngleOfAttack.lowerBound  = nan;
            
            evolutions(k).max_tof                     = nan;
            %         evolutions(k).max_interp_E_mapped_Ascent  = nan;
            %         evolutions(k).max_interp_E_mapped_Descent = nan;
            %         evolutions(k).population.individual       = nan;
            %         evolutions(k).individuals.v_i             = nan;
            %         evolutions(k).individuals.gamma_i         = nan;
            %         evolutions(k).individuals.chi_i           = nan;
            %         evolutions(k).fitness.dif_norm            = nan;
            %         evolutions(k).fitness.dif_lat             = nan;
            %         evolutions(k).fitness.dif_lon             = nan;
            %         evolutions(k).fitness.dif_d_deg           = nan;
            %         evolutions(k).fitness.dif_h               = nan;
            %         evolutions(k).fitness.tof                 = nan;
            %         evolutions(k).best(3).criteria            = nan;
            %         evolutions(k).best(3).v_i                 = nan;
            %         evolutions(k).best(3).gamma_i             = nan;
            %         evolutions(k).best(3).chi_i               = nan;
            %         evolutions(k).best(3).dif_d_deg           = nan;
            %         evolutions(k).best(3).dif_h               = nan;
            %         evolutions(k).best(3).tof                 = nan;
            %         evolutions(k).best(3).max_tof             = nan;
        end
end


end

