function [ evolutions ] = initializeDecisionVectorNodalParameterFields( simulationDetails, evolutions, k, i )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
trajectoryType = simulationDetails.trajectoryType;

decisionVectorList_Nodal = readcell( 'decisionVectorList_Nodal.dat','Delimiter',{','} );

if trajectoryType == 'A'
    phase = {'ascent'};
elseif trajectoryType == 'D'
    phase = {'descent'};
elseif trajectoryType == 'AD'
    phase = [{'ascent'} {'descent'}];
end

for p = 1:numel(phase)
    j = 1;
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
    
    evolutions(k).population(i).interpolators.(phase{p}).InputData.normalizedSpecificEnergy  = nan;
    evolutions(k).population(i).interpolators.(phase{p}).Evaluation.normalizedSpecificEnergy = nan;
    
    if length(find(simulationDetails.nodalParameters=='A')) > 0
        j = 2;
        evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
        evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
        
    end
    if length(find(simulationDetails.nodalParameters=='S')) > 0
        j = 3;
        evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
        evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
    end
    if length(find(simulationDetails.nodalParameters=='B')) > 0
        j = 4;
        evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
        evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
    end
    if length(find(simulationDetails.nodalParameters=='T')) > 0
        j = 5;
        evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
        evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
        evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
    end
    
end

%
% elseif trajectoryType == 'D'
%     if length(find(simulationDetails.nodalParameters=='A')) > 0
%         evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.figureTitleContent             = 'Node Interval';
%         evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.figureSaveNameContent          = 'ascentNodeInterval';
%         evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.units                          = '(-)';
%         evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.data.one                       = nan;
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.figureTitleContent            = 'Angle of Attack';
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.figureSaveNameContent         = 'ascentAngleOfAttack';
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.units                         = '(deg)';
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.data.one                      = nan;
%     end
%     if length(find(simulationDetails.nodalParameters=='S')) > 0
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.figureTitleContent          = 'Angle of Sideslip';
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.figureSaveNameContent       = 'ascentAngleOfSideslip';
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.units                       = '(deg)';
%         evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.data.one                    = nan;
%     end
%     if length(find(simulationDetails.nodalParameters=='B')) > 0
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureTitleContent                = 'Bank Angle';
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureSaveNameContent             = 'ascentBankAngle';
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.units                             = '(deg)';
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.data.one                          = nan;
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureTitleContent                = 'Bank Angle';
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureSaveNameContent             = 'ascentBankAngle';
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.units                             = '(deg)';
%         evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.data.one                          = nan;
%     end
%     if length(find(simulationDetails.nodalParameters=='T')) > 0
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureTitleContent          = 'Throttle Setting';
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureSaveNameContent       = 'ascentThrottleSetting';
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.units                       = '(-)';
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.data.one                    = nan;
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureTitleContent          = 'Throttle Setting';
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureSaveNameContent       = 'ascentThrottleSetting';
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.units                       = '(-)';
%         evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.data.one                    = nan;
%     end
%
%     evolutions(k).population(i).interpolators.Descent.InputData.normalizedSpecificEnergy 	=	nan;
%     evolutions(k).population(i).interpolators.Descent.InputData.angleOfAttack        	=	nan;
%     evolutions(k).population(i).interpolators.Descent.InputData.angleOfSideSlip        	=	nan;
%     evolutions(k).population(i).interpolators.Descent.InputData.bankAngle            	=	nan;
%     % evolutions(k).population(i).interpolators.Descent.InputData.thrustElevationAngle 	=	nan;
%     % evolutions(k).population(i).interpolators.Descent.InputData.thrustAzimuthAngle   	=	nan;
%     evolutions(k).population(i).interpolators.Descent.InputData.throttleSetting      	=	nan;
%     %evolutions(k).population(i).interpolators.Descent.InputData.nodeLocation         	=	nan;
%     evolutions(k).population(i).interpolators.Descent.Evaluation.normalizedSpecificEnergy 	=	nan;
%     evolutions(k).population(i).interpolators.Descent.Evaluation.angleOfAttack        	=	nan;
%     evolutions(k).population(i).interpolators.Descent.Evaluation.angleOfSideSlip        	=	nan;
%     evolutions(k).population(i).interpolators.Descent.Evaluation.bankAngle            	=	nan;
%     %evolutions(k).population(i).interpolators.Descent.Evaluation.thrustElevationAngle 	=	nan;
%     %evolutions(k).population(i).interpolators.Descent.Evaluation.thrustAzimuthAngle   	=	nan;
%     evolutions(k).population(i).interpolators.Descent.Evaluation.throttleSetting      	=	nan;
%
%     elseif trajectoryType == 'AD'
%         if length(find(simulationDetails.nodalParameters=='A')) > 0
%             evolutions(k).population(i).decisionVector.parameters.Ascent.nodeInterval.figureTitleContent             = 'Node Interval';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.nodeInterval.figureSaveNameContent          = 'ascentNodeInterval';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.nodeInterval.units                          = '(-)';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.nodeInterval.data.one                       = nan;
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfAttack.figureTitleContent            = 'Angle of Attack';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfAttack.figureSaveNameContent         = 'ascentAngleOfAttack';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfAttack.units                         = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfAttack.data.one                      = nan;
%
%             evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.figureTitleContent             = 'Node Interval';
%             evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.figureSaveNameContent          = 'ascentNodeInterval';
%             evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.units                          = '(-)';
%             evolutions(k).population(i).decisionVector.parameters.Descent.nodeInterval.data.one                       = nan;
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.figureTitleContent            = 'Angle of Attack';
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.figureSaveNameContent         = 'ascentAngleOfAttack';
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.units                         = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfAttack.data.one                      = nan;
%         end
%         if length(find(simulationDetails.nodalParameters=='S')) > 0
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.figureTitleContent          = 'Angle of Sideslip';
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.figureSaveNameContent       = 'ascentAngleOfSideslip';
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.units                       = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Descent.angleOfSideslip.data.one                    = nan;
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfSideslip.figureTitleContent          = 'Angle of Sideslip';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfSideslip.figureSaveNameContent       = 'ascentAngleOfSideslip';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfSideslip.units                       = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.angleOfSideslip.data.one                    = nan;
%         end
%
%         if length(find(simulationDetails.nodalParameters=='B')) > 0
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.figureTitleContent                = 'Bank Angle';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.figureSaveNameContent             = 'ascentBankAngle';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.units                             = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.data.one                          = nan;
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.figureTitleContent                = 'Bank Angle';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.figureSaveNameContent             = 'ascentBankAngle';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.units                             = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.bankAngle.data.one                          = nan;
%
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureTitleContent                = 'Bank Angle';
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureSaveNameContent             = 'ascentBankAngle';
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.units                             = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.data.one                          = nan;
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureTitleContent                = 'Bank Angle';
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.figureSaveNameContent             = 'ascentBankAngle';
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.units                             = '(deg)';
%             evolutions(k).population(i).decisionVector.parameters.Descent.bankAngle.data.one                          = nan;
%         end
%         if length(find(simulationDetails.nodalParameters=='T')) > 0
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.figureTitleContent          = 'Throttle Setting';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.figureSaveNameContent       = 'ascentThrottleSetting';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.units                       = '(-)';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.data.one                    = nan;
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.figureTitleContent          = 'Throttle Setting';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.figureSaveNameContent       = 'ascentThrottleSetting';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.units                       = '(-)';
%             evolutions(k).population(i).decisionVector.parameters.Ascent.throttleSetting.data.one                    = nan;
%
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureTitleContent          = 'Throttle Setting';
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureSaveNameContent       = 'ascentThrottleSetting';
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.units                       = '(-)';
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.data.one                    = nan;
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureTitleContent          = 'Throttle Setting';
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.figureSaveNameContent       = 'ascentThrottleSetting';
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.units                       = '(-)';
%             evolutions(k).population(i).decisionVector.parameters.Descent.throttleSetting.data.one                    = nan;
%         end
%
%         evolutions(k).population(i).interpolators.Ascent.InputData.normalizedSpecificEnergy 	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.InputData.angleOfAttack        	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.InputData.angleOfSideSlip        	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.InputData.bankAngle            	=	nan;
%         %evolutions(k).population(i).interpolators.Ascent.InputData.thrustElevationAngle 	=	nan;
%         %evolutions(k).population(i).interpolators.Ascent.InputData.thrustAzimuthAngle   	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.InputData.throttleSetting      	=	nan;
%         %evolutions(k).population(i).interpolators.Ascent.InputData.nodeLocation         	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.Evaluation.normalizedSpecificEnergy 	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.Evaluation.angleOfAttack        	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.Evaluation.angleOfSideSlip        	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.Evaluation.bankAngle            	=	nan;
%         %evolutions(k).population(i).interpolators.Ascent.Evaluation.thrustElevationAngle 	=	nan;
%         %evolutions(k).population(i).interpolators.Ascent.Evaluation.thrustAzimuthAngle   	=	nan;
%         evolutions(k).population(i).interpolators.Ascent.Evaluation.throttleSetting      	=	nan;
%
%
%
%         evolutions(k).population(i).interpolators.Descent.InputData.normalizedSpecificEnergy 	=	nan;
%         evolutions(k).population(i).interpolators.Descent.InputData.angleOfAttack        	=	nan;
%         evolutions(k).population(i).interpolators.Descent.InputData.angleOfSideSlip        	=	nan;
%         evolutions(k).population(i).interpolators.Descent.InputData.bankAngle            	=	nan;
%         % evolutions(k).population(i).interpolators.Descent.InputData.thrustElevationAngle 	=	nan;
%         % evolutions(k).population(i).interpolators.Descent.InputData.thrustAzimuthAngle   	=	nan;
%         evolutions(k).population(i).interpolators.Descent.InputData.throttleSetting      	=	nan;
%         %evolutions(k).population(i).interpolators.Descent.InputData.nodeLocation         	=	nan;
%         evolutions(k).population(i).interpolators.Descent.Evaluation.normalizedSpecificEnergy 	=	nan;
%         evolutions(k).population(i).interpolators.Descent.Evaluation.angleOfAttack        	=	nan;
%         evolutions(k).population(i).interpolators.Descent.Evaluation.angleOfSideSlip        	=	nan;
%         evolutions(k).population(i).interpolators.Descent.Evaluation.bankAngle            	=	nan;
%         %evolutions(k).population(i).interpolators.Descent.Evaluation.thrustElevationAngle 	=	nan;
%         %evolutions(k).population(i).interpolators.Descent.Evaluation.thrustAzimuthAngle   	=	nan;
%         evolutions(k).population(i).interpolators.Descent.Evaluation.throttleSetting      	=	nan;
%
% end


end

