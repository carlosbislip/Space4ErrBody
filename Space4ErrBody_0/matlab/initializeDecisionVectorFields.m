function [ evolutions ] = initializeDecisionVectorFields( simulationDetails, evolutions, k, i )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

decisionVectorList_Common = readcell( 'decisionVectorList_Common.dat','Delimiter',{','} );

for j = 1:size(decisionVectorList_Common,1)
    evolutions(k).population(i).decisionVector.parameters.common.(decisionVectorList_Common{j,1}).figureSaveNameContent = decisionVectorList_Common{j,1};
    evolutions(k).population(i).decisionVector.parameters.common.(decisionVectorList_Common{j,1}).variableLabel = decisionVectorList_Common{j,2};
    evolutions(k).population(i).decisionVector.parameters.common.(decisionVectorList_Common{j,1}).units = decisionVectorList_Common{j,3};
    evolutions(k).population(i).decisionVector.parameters.common.(decisionVectorList_Common{j,1}).scalingFactor = decisionVectorList_Common{j,4};
    evolutions(k).population(i).decisionVector.parameters.common.(decisionVectorList_Common{j,1}).limits = [decisionVectorList_Common{j,5} decisionVectorList_Common{j,6}];
    evolutions(k).population(i).decisionVector.parameters.common.(decisionVectorList_Common{j,1}).tick = decisionVectorList_Common{j,7};
    evolutions(k).population(i).decisionVector.parameters.common.(decisionVectorList_Common{j,1}).data.one = nan;
end

decisionVectorList_Nodal = readcell( 'decisionVectorList_Nodal.dat','Delimiter',{','} );

trajectoryType = simulationDetails.trajectoryType;
switch (trajectoryType)
    case 'A'
        phase = {'ascent'};
    case 'D'
        phase = {'descent'};
    case 'AD'
        phase = [{'ascent'} {'descent'}];
end


nodalParametersMatrix = readcell( string(strcat('/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/Space4ErrBody_Headers/executables/',simulationDetails.NodalParametersFileName)),'Delimiter',{'\t'} );

ids = [];
for ii = 1:size(nodalParametersMatrix,1)
    if nodalParametersMatrix{ii,1}(1) == '#'
        ids = [ids;ii];
    end
end

nodalParametersMatrix(ids,:) = [];

nodalParametersMatrix = cell2mat(nodalParametersMatrix);

ascentNodalParametersMatrix = nodalParametersMatrix(find(nodalParametersMatrix(:,1) == 0),:);

parameters(1).phase = unique(ascentNodalParametersMatrix(:,2));

descentNodalParametersMatrix = nodalParametersMatrix(find(nodalParametersMatrix(:,1) == 1),:);
parameters(2).phase = unique(descentNodalParametersMatrix(:,2));

possibleNodalParameters = [{'A'} {'S'} {'B'} {'T'}];
for p = 1:numel(phase)
    j = 2;
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
    evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
    
    for kkk = 1:numel(possibleNodalParameters)
        if isempty(find(simulationDetails.nodalParameters==possibleNodalParameters{kkk},1)) == 0 && isempty(find(parameters(p).phase==kkk-1,1)) == 0
            j = kkk + 2;
            
            evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
            evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
            evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
            evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
            evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
            evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
            
        end
        
    end
end

%%
mapping =  [{'inputData'} {'evaluation'} ];

for p = 1:numel(phase)
    for pp = 1
        for j = 1:size(decisionVectorList_Nodal,1)
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).figureSaveNameContent = decisionVectorList_Nodal{j,1};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).value = nan;
        end
    end
end

indices = 1:size(decisionVectorList_Nodal,1);
indices(2) = [];

for p = 1:numel(phase)
    for pp = 2
        for j = indices
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).figureSaveNameContent = decisionVectorList_Nodal{j,1};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
            evolutions(k).population(i).interpolators.(phase{p}).(mapping{pp}).(decisionVectorList_Nodal{j,1}).value = nan;
        end
    end
end



%
%     end
%
%     if isempty(find(simulationDetails.nodalParameters=='A',1)) == 0 && isempty(find(parameters(p).phase==0,1)) == 0
%         j = 2;
%         evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
%         evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
%         
%     end
%     if isempty(find(simulationDetails.nodalParameters=='S',1)) == 0 && isempty(find(parameters(p).phase==1,1)) == 0 
%         j = 3;
%         evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
%         evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
%     end
%     if isempty(find(simulationDetails.nodalParameters=='B',1)) == 0 && isempty(find(parameters(p).phase==2,1)) == 0 
%         j = 4;
%         evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
%         evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
%     end
%     if isempty(find(simulationDetails.nodalParameters=='T',1)) == 0 && isempty(find(parameters(p).phase==3,1)) == 0 
%         j = 5;
%         evolutions(k).population(i).interpolators.(phase{p}).InputData.(decisionVectorList_Nodal{j,1})  = nan;
%         evolutions(k).population(i).interpolators.(phase{p}).Evaluation.(decisionVectorList_Nodal{j,1})	= nan;
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).variableLabel = decisionVectorList_Nodal{j,2};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).units = decisionVectorList_Nodal{j,3};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).scalingFactor = decisionVectorList_Nodal{j,4};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).limits = [decisionVectorList_Nodal{j,5} decisionVectorList_Nodal{j,6}];
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).tick = decisionVectorList_Nodal{j,7};
%         evolutions(k).population(i).decisionVector.parameters.(phase{p}).(decisionVectorList_Nodal{j,1}).data.one = nan;
%     end
%     
% end



%
% evolutions(k).population(i).decisionVector.parameters.Common.initialFlightPathAngle.figureTitleContent           =        'Initial Flight-Path Angle';
% evolutions(k).population(i).decisionVector.parameters.Common.initialFlightPathAngle.figureSaveNameContent        =        'initialFlightPathAngle';
% evolutions(k).population(i).decisionVector.parameters.Common.initialFlightPathAngle.units                        =        '(deg)';
% evolutions(k).population(i).decisionVector.parameters.Common.initialFlightPathAngle.data.one                     =       nan;
% evolutions(k).population(i).decisionVector.parameters.Common.initialLaunchHeadingAngle.figureTitleContent        =        'Initial Launch Heading Angle';
% evolutions(k).population(i).decisionVector.parameters.Common.initialLaunchHeadingAngle.figureSaveNameContent     =        'initialLaunchHeadingAngle';
% evolutions(k).population(i).decisionVector.parameters.Common.initialLaunchHeadingAngle.units                     =        '(deg)';
% evolutions(k).population(i).decisionVector.parameters.Common.initialLaunchHeadingAngle.data.one                  = nan;
% evolutions(k).population(i).decisionVector.parameters.Common.initialVelocity.figureTitleContent                  = 'Initial Velocity';
% evolutions(k).population(i).decisionVector.parameters.Common.initialVelocity.figureSaveNameContent               = 'initialVelocity';
% evolutions(k).population(i).decisionVector.parameters.Common.initialVelocity.units                               = '(m/s)';
% evolutions(k).population(i).decisionVector.parameters.Common.initialVelocity.data.one                            = nan;
% evolutions(k).population(i).decisionVector.parameters.Common.maximumVelocity.figureTitleContent                  =        'Maximum Velocity';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumVelocity.figureSaveNameContent               =        'maximumVelocity';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumVelocity.units                               =        '(m/s)';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumVelocity.data.one                            =       nan;
% evolutions(k).population(i).decisionVector.parameters.Common.maximumHeight.figureTitleContent                    =        'Maximum Height';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumHeight.figureSaveNameContent                 =        'maximumHeight';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumHeight.units                                 =        '(m)';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumHeight.data.one                              =       nan;
% evolutions(k).population(i).decisionVector.parameters.Common.additionalMass.figureTitleContent                   =        'Additional Mass';
% evolutions(k).population(i).decisionVector.parameters.Common.additionalMass.figureSaveNameContent                =        'additionalMass';
% evolutions(k).population(i).decisionVector.parameters.Common.additionalMass.units                                =        '(kg)';
% evolutions(k).population(i).decisionVector.parameters.Common.additionalMass.data.one                             =       nan;
% evolutions(k).population(i).decisionVector.parameters.Common.terminationDistanceRatio.figureTitleContent         =        'Termination Distance Ratio';
% evolutions(k).population(i).decisionVector.parameters.Common.terminationDistanceRatio.figureSaveNameContent      =        'terminationDistanceRatio';
% evolutions(k).population(i).decisionVector.parameters.Common.terminationDistanceRatio.units                      =        '(-)';
% evolutions(k).population(i).decisionVector.parameters.Common.terminationDistanceRatio.data.one                   =       nan;
% evolutions(k).population(i).decisionVector.parameters.Common.finalVelocity.figureTitleContent                    =        'Final Velocity';
% evolutions(k).population(i).decisionVector.parameters.Common.finalVelocity.figureSaveNameContent                 =        'finalVelocity';
% evolutions(k).population(i).decisionVector.parameters.Common.finalVelocity.units                                 =        '(m/s)';
% evolutions(k).population(i).decisionVector.parameters.Common.finalVelocity.data.one                              =       nan;
% evolutions(k).population(i).decisionVector.parameters.Common.skipSuppressionTriggerTime.figureTitleContent       =        'Skip Suppression Trigger Time';
% evolutions(k).population(i).decisionVector.parameters.Common.skipSuppressionTriggerTime.figureSaveNameContent    =        'skipSuppressionTriggerTime';
% evolutions(k).population(i).decisionVector.parameters.Common.skipSuppressionTriggerTime.units                    =        '(s)';
% evolutions(k).population(i).decisionVector.parameters.Common.skipSuppressionTriggerTime.data.one                 =       nan;
% evolutions(k).population(i).decisionVector.parameters.Common.maximumMechanicalLoadAscent.figureTitleContent      =        'Maximum Mechanical Load - Ascent';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumMechanicalLoadAscent.figureSaveNameContent   =        'maximumMechanicalLoadAscent';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumMechanicalLoadAscent.units                   =        '(g_0)';
% evolutions(k).population(i).decisionVector.parameters.Common.maximumMechanicalLoadAscent.data.one                =       nan;


end

