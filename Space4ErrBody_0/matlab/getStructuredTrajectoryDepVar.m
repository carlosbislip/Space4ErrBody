function [ evolutions ] = getStructuredTrajectoryDepVar( evolutions, folderPath, rawData, totalSets, setNumber, validation )
%GET_trajectory Summary of this function goes here
%   Detailed explanation goes here
disp(char(strcat('Getting trajectories for set',{' '},num2str(setNumber),{' '},'/',{' '},num2str(totalSets),':',{' '})))

%% Initialize Fields

% i = rawData.populationSize;
% for k = [ 1 numel(rawData.generationList) ]
%     evolutions(k).evolution                                             = nan;   
%     evolutions(k).population(i).size.collective                         = nan;
%     evolutions(k).population(i).size.printed                            = nan;
%     evolutions(k).population(i).size.nonPrinted                         = nan;
%     evolutions(k).population(i).size.fitness                            = nan;
%     evolutions(k).population(i).indices.nonPrinted                      = nan;
%     evolutions(k).population(i).indices.nonDominatedFront               = nan;
%     evolutions(k).population(i).indices.trajectoryPhaseChange           = nan;
%     evolutions(k).population(i).name                                    = nan;
%     evolutions(k).population(i).printed                                 = nan;
%     
%     evolutions(k).population(i).decisionVector.collective               = nan;
%     evolutions(k).population(i).decisionVector.printed                  = nan;
%     evolutions(k).population(i).decisionVector.nonPrinted               = nan;
%     
%     if( validation == 0 )
%         [ evolutions ] = initializeDecisionVectorFields( rawData, evolutions, k, i );
%         [ evolutions ] = initializeFitnessVectorFields( rawData, evolutions, k, i );
%     end
%     [ evolutions ] = initializeExtremesAndConstraintsFields( evolutions, k, i );
%     [ evolutions ] = initializeDependentVariableFields( evolutions, k, i );
%     
%     
%     evolutions(k).max_tof                     = nan;
%     
% end
% evolutions(1).Common.Bounds.headingError.angularDistanceToGo  = nan;
% evolutions(1).Common.Bounds.headingError.upperBound  = nan;
% evolutions(1).Common.Bounds.headingError.lowerBound  = nan;
% evolutions(1).Common.Bounds.AngleOfAttack.machNumber  = nan;
% evolutions(1).Common.Bounds.AngleOfAttack.upperBound  = nan;
% evolutions(1).Common.Bounds.AngleOfAttack.lowerBound  = nan;
% 


%%

% pp = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(pp)
%     parpool(4);
% end
% 
% strWindowTitle = convertCharsToStrings(char(strcat('Getting individuals for set',{' '},num2str(setNumber),{' '},'/',{' '},num2str(totalSets),':',{' '})));
% ppm_population = ParforProgMon(strWindowTitle, rawData.maxGenerations*rawData.populationSize);
% 
% 
% strWindowTitle = convertCharsToStrings(char(strcat('Getting trajectories - Generational progress of set',{' '},num2str(setNumber),{' '},'/',{' '},num2str(totalSets),':',{' '})));
% ppm_generation = ParforProgMon(strWindowTitle, rawData.maxGenerations);

%%
Folder_prefix = 'OUTPUT*/';
prop_File_prefix = 'propagationHistory/propagationHistory_';
depvar_File_prefix = 'dependentVariables/dependentVariables_';
interp_Ascent_File_prefix = 'evaluatedInterpolatorsAscent/evaluatedInterpolatorsAscent_';
interp_Ascent_LB_File_prefix = 'evaluatedInterpolatorsAscent_LB/evaluatedInterpolatorsAscent_';
interp_Ascent_UB_File_prefix = 'evaluatedInterpolatorsAscent_UB/evaluatedInterpolatorsAscent_';
interp_Descent_File_prefix = 'evaluatedInterpolatorsDescent/evaluatedInterpolatorsDescent_';
interp_Descent_LB_File_prefix = 'evaluatedInterpolatorsDescent_LB/evaluatedInterpolatorsDescent_';
interp_Descent_UB_File_prefix = 'evaluatedInterpolatorsDescent_UB/evaluatedInterpolatorsDescent_';
%monteCarloPopulation_prefix = 'monteCarloPopulation*';

DV_mapped_Ascent_File_prefix = 'map_DV_mapped_Ascent/map_DV_mapped_Ascent_';
DV_mapped_Descent_File_prefix = 'map_DV_mapped_Descent/map_DV_mapped_Descent_';

dependentVariableNames        = fieldnames(evolutions(end).population(end).dependentVariableTimeHistory);
decisionVectorCommonParameterNames = fieldnames(evolutions(end).population(end).decisionVector.parameters.common);
extremesAndContraintsFieldNames = fieldnames(evolutions(end).population(end).extremesAndConstraints);

trajectoryType = rawData.trajectoryType;

%if validation == 0
    nodes = rawData.nodes;
    nodeList = cellstr([ "one";"two";"three";"four";"five";"six";"seven";"eight";"nine";"ten";"eleven";"twelve" ]);
    commonFitnessVectorFieldNames  = fieldnames(evolutions(end).population(end).fitnessVector.common);
    
    ascentPhaseFitnessVectorFieldNames = [];
    descentPhaseFitnessVectorFieldNames = [];
    interpolatorsInputDataFields = [];
    interpolatorsEvaluationFields = [];
    phaseFitnessVectorFieldNames = [];
    
    if isfield(evolutions(end).population(end).fitnessVector,'ascent') || isfield(evolutions(end).population(end).fitnessVector,'descent')
        if trajectoryType == 'A'
            ascentPhaseFitnessVectorFieldNames  = fieldnames(evolutions(end).population(end).fitnessVector.ascent);
            interpolatorsInputDataFields  = fieldnames(evolutions(end).population(end).interpolators.ascent.inputData);
            interpolatorsEvaluationFields = fieldnames(evolutions(end).population(end).interpolators.ascent.evaluation);
            phaseFitnessVectorFieldNames  = fieldnames(evolutions(end).population(end).fitnessVector.ascent);
        elseif trajectoryType == 'D'
            interpolatorsInputDataFields  = fieldnames(evolutions(end).population(end).interpolators.descent.inputData);
            interpolatorsEvaluationFields = fieldnames(evolutions(end).population(end).interpolators.descent.evaluation);
            decisionVectorDescentPhaseParameterNames  = fieldnames(evolutions(end).population(end).decisionVector.parameters.descent);
            phaseFitnessVectorFieldNames  = fieldnames(evolutions(end).population(end).fitnessVector.descent);
        elseif trajectoryType == 'AD'
            interpolatorsInputDataFields  = fieldnames(evolutions(end).population(end).interpolators.ascent.inputData);
            interpolatorsEvaluationFields = fieldnames(evolutions(end).population(end).interpolators.ascent.evaluation);
            decisionVectorAscentPhaseParameterNames  = fieldnames(evolutions(end).population(end).decisionVector.parameters.ascent);
            decisionVectorDescentPhaseParameterNames  = fieldnames(evolutions(end).population(end).decisionVector.parameters.descent);
            phaseFitnessVectorFieldNames  = fieldnames(evolutions(end).population(end).fitnessVector.ascent);
        end
    end
    
    

headingErrorDeadbandBoundsAngularDistanceToGo  = rawData.headingErrorDeadbandBounds(:,1);
headingErrorDeadbandBoundsLowerBound  = rawData.headingErrorDeadbandBounds(:,2);
headingErrorDeadbandBoundsUpperBound  = rawData.headingErrorDeadbandBounds(:,3);
alphaMachBoundsMachNumber  =  rawData.alphaMachBounds(:,1);
alphaMachBoundsLowerBound  =  rawData.alphaMachBounds(:,2);
alphaMachBoundsUpperBound  =  rawData.alphaMachBounds(:,3);

evolutions(1).Common.Bounds.headingErrorDeadBand.angularDistanceToGo  = headingErrorDeadbandBoundsAngularDistanceToGo;
evolutions(1).Common.Bounds.headingErrorDeadBand.lowerBound  = headingErrorDeadbandBoundsLowerBound;
evolutions(1).Common.Bounds.headingErrorDeadBand.upperBound  = headingErrorDeadbandBoundsUpperBound;
evolutions(1).Common.Bounds.AngleOfAttack.machNumber  = alphaMachBoundsMachNumber;
evolutions(1).Common.Bounds.AngleOfAttack.lowerBound  = alphaMachBoundsLowerBound;
evolutions(1).Common.Bounds.AngleOfAttack.upperBound  = alphaMachBoundsUpperBound;

populationDataPerGeneration = rawData.populationDataPerGeneration;
fitnessDataPerGeneration = rawData.fitnessDataPerGeneration;
extremesAndConstraintsDataPerGeneration = rawData.extremesAndConstraintsDataPerGeneration;

fitnessMagnitudeAverageHistory = rawData.fitnessMagnitudeAverageHistory;

for k = 1:numel(rawData.generationList)
    
    evolutions(k).evolution = rawData.generationList(k);

    evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo  = headingErrorDeadbandBoundsAngularDistanceToGo;
    evolutions(k).Common.Bounds.headingErrorDeadBand.lowerBound  = headingErrorDeadbandBoundsLowerBound;
    evolutions(k).Common.Bounds.headingErrorDeadBand.upperBound  = headingErrorDeadbandBoundsUpperBound;
    evolutions(k).Common.Bounds.AngleOfAttack.machNumber  = alphaMachBoundsMachNumber;
    evolutions(k).Common.Bounds.AngleOfAttack.lowerBound  = alphaMachBoundsLowerBound;
    evolutions(k).Common.Bounds.AngleOfAttack.upperBound  = alphaMachBoundsUpperBound;
        
    %pp = 1;
    populationDetailsPerCurrentGeneration = populationDataPerGeneration(:,1:6,k);
    
    populationDataPerCurrentGeneration = populationDataPerGeneration(:,7:end,k);
    fitnessDataPerCurrentGeneration = fitnessDataPerGeneration(:,7:end,k);
    extremesAndConstraintsDataPerCurrentGeneration = extremesAndConstraintsDataPerGeneration(:,7:end,k);
    
    printedIndices       = find([populationDetailsPerCurrentGeneration{:,2}] > 0);
    nonPrintedIndices    = find([populationDetailsPerCurrentGeneration{:,2}] == 0);
    printedPopulation    = populationDataPerCurrentGeneration(printedIndices,:);
    nonPrintedPopulation = populationDataPerCurrentGeneration(nonPrintedIndices,:);
    
    evolutions(k).population(1).size.collective           = rawData.populationSize;
    evolutions(k).population(1).size.printed              = size(printedPopulation,1);
    evolutions(k).population(1).size.nonPrinted           = size(nonPrintedPopulation,1);
    evolutions(k).population(1).indices.printed           = printedIndices;
    evolutions(k).population(1).indices.nonPrinted        = nonPrintedIndices;
    
    rankDataPerCurrentGeneration = populationDataPerGeneration(:,3:4,k);
    
    totalTimeOfFlight = nan(evolutions(k).population(1).size.collective,1);
    fitnessVectorMagnitude = nan(evolutions(k).population(1).size.collective,1);
    evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages = fitnessMagnitudeAverageHistory(k,2:end);
    
    
    for p = 1:rawData.populationSize
        
        evolutions(k).population(p).fitnessVector.RankBySortedVector               = rankDataPerCurrentGeneration{p,1} + 1;
        evolutions(k).population(p).fitnessVector.RankBySortedMagnitude            = rankDataPerCurrentGeneration{p,2} + 1;
             
        evolutions(k).population(p).name          = populationDetailsPerCurrentGeneration{p,1};
        evolutions(k).population(p).printed       = populationDetailsPerCurrentGeneration{p,2};
        evolutions(k).population(p).rankVector    = rankDataPerCurrentGeneration{p,1} + 1;
        evolutions(k).population(p).rankMagnitude = rankDataPerCurrentGeneration{p,2} + 1;
        n = 1;
        
        for i = 1:length(decisionVectorCommonParameterNames)
            evolutions(k).population(p).decisionVector.parameters.common.(decisionVectorCommonParameterNames{i}).data.one = populationDataPerCurrentGeneration{p,n};
            n = n + 1;
        end
        
       if( validation == 0 )
            switch (trajectoryType)
                case 'A'
                    for j = 1:(nodes - 1)
                        evolutions(k).population(p).decisionVector.parameters.ascent.(decisionVectorAscentPhaseParameterNames{1}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                        n = n + 1;
                    end
                    for i = 2:length(decisionVectorAscentPhaseParameterNames)
                        for j = 1:nodes
                            evolutions(k).population(p).decisionVector.parameters.ascent.(decisionVectorAscentPhaseParameterNames{i}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                            n = n + 1;
                        end
                    end
                    
                case 'D'
                    for j = 1:(nodes - 1)
                        evolutions(k).population(p).decisionVector.parameters.descent.(decisionVectorDescentPhaseParameterNames{1}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                        n = n + 1;
                    end
                    for i = 2:length(decisionVectorDescentPhaseParameterNames)
                        for j = 1:nodes
                            evolutions(k).population(p).decisionVector.parameters.descent.(decisionVectorDescentPhaseParameterNames{i}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                            n = n + 1;
                        end
                    end
                    
                case 'AD'
                    for j = 1:(nodes - 1)
                        evolutions(k).population(p).decisionVector.parameters.ascent.(decisionVectorAscentPhaseParameterNames{1}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                        n = n + 1;
                    end
                    for j = 1:(nodes - 1)
                        evolutions(k).population(p).decisionVector.parameters.descent.(decisionVectorDescentPhaseParameterNames{1}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                        n = n + 1;
                    end
                    for i = 2:length(decisionVectorAscentPhaseParameterNames)
                        for j = 1:nodes
                            evolutions(k).population(p).decisionVector.parameters.ascent.(decisionVectorAscentPhaseParameterNames{i}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                            n = n + 1;
                        end
                    end
                    
                    for i = 2:length(decisionVectorDescentPhaseParameterNames)
                        for j = 1:nodes
                            evolutions(k).population(p).decisionVector.parameters.descent.(decisionVectorDescentPhaseParameterNames{i}).data.(nodeList{j}) = populationDataPerCurrentGeneration{p,n};
                            n = n + 1;
                        end
                    end
            end
        end
        
        for i = 1:length(extremesAndContraintsFieldNames)
            evolutions(k).population(p).extremesAndConstraints.(extremesAndContraintsFieldNames{i}).value = extremesAndConstraintsDataPerCurrentGeneration{p,i};
        end
        
        if evolutions(k).population(p).printed == 1
            
            depvar_source            = { convertStringsToChars(strcat( folderPath,"/",depvar_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
            if isfile(depvar_source{:})
                
                fid = fopen(depvar_source{:});
                depvar = dlmread(depvar_source{:},',');
                fclose(fid);
                
                for i = 1:length(dependentVariableNames) - 6
                    evolutions(k).population(p).dependentVariableTimeHistory.(dependentVariableNames{i}).value = depvar(:,i);
                end
                
                evolutions(k).population(p).dependentVariableTimeHistory.acc_x.value      = ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_x.value ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_x.value ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_x.value;
                
                evolutions(k).population(p).dependentVariableTimeHistory.acc_y.value      = ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_y.value ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_y.value ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_y.value;
                
                evolutions(k).population(p).dependentVariableTimeHistory.acc_z.value      = ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_z.value ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_z.value ...
                    + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_z.value;
                
                evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_M.value = sqrt(...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_x.value).^2 ...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_y.value).^2 ...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_z.value).^2);
                
                evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_M.value = sqrt(...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_x.value).^2 ...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_y.value).^2 ...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_z.value).^2);
                
                evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_M.value = sqrt(...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_x.value).^2 ...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_y.value).^2 ...
                    + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_z.value).^2);
                
                I_distanceCoveredRatioAboveTerminationRatio = find( evolutions(k).population(p).dependentVariableTimeHistory.angularDistanceCoveredRatio.value > ...
                    evolutions(k).population(p).decisionVector.parameters.common.terminationDistanceRatio.data.one );
                if isempty(I_distanceCoveredRatioAboveTerminationRatio) == true
                    I_distanceCoveredRatioAboveTerminationRatio = nan;
                end
                
                I_trajectoryPhaseChange =  find( evolutions(k).population(p).dependentVariableTimeHistory.trajectoryPhase.value > 0 );
                
                if isempty(I_trajectoryPhaseChange) == true
                    if isnan(I_distanceCoveredRatioAboveTerminationRatio) == true
                        I_trajectoryPhaseChange = nan;
                    else
                        I_trajectoryPhaseChange = length(evolutions(k).population(p).dependentVariableTimeHistory.trajectoryPhase.value);
                    end
                end
                evolutions(k).population(p).indices.trajectoryPhaseChange = ...
                    [ I_distanceCoveredRatioAboveTerminationRatio(1) I_trajectoryPhaseChange(1) ];
                
                totalTimeOfFlight(p)               = evolutions(k).population(p).dependentVariableTimeHistory.timeOfFlight.value(end);
                evolutions(k).fitness.tof(p)       = evolutions(k).population(p).dependentVariableTimeHistory.timeOfFlight.value(end);
                
                if( validation == 0 )
                    
                    DV_mapped_Ascent_source  = { convertStringsToChars(strcat( folderPath,"/",DV_mapped_Ascent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    DV_mapped_Descent_source = { convertStringsToChars(strcat( folderPath,"/",DV_mapped_Descent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    
                    
                    fid = fopen(DV_mapped_Ascent_source{:});
                    DV_mapped_Ascent = dlmread(DV_mapped_Ascent_source{:},',');
                    fclose(fid);
                    
                    fid = fopen(DV_mapped_Descent_source{:});
                    DV_mapped_Descent = dlmread(DV_mapped_Descent_source{:},',');
                    fclose(fid);
                    
                    for i = 1:length(interpolatorsInputDataFields)
                        evolutions(k).population(p).interpolators.ascent.inputData.(interpolatorsInputDataFields{i}).value  = DV_mapped_Ascent(:,i);
                        evolutions(k).population(p).interpolators.descent.inputData.(interpolatorsInputDataFields{i}).value = DV_mapped_Descent(:,i);
                    end
                    
                    interp_Ascent_source     = { convertStringsToChars(strcat( folderPath,"/",interp_Ascent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    fid = fopen(interp_Ascent_source{:});
                    interp_Ascent = dlmread(interp_Ascent_source{:},',');
                    fclose(fid);
                    
                    interp_Descent_source    = { convertStringsToChars(strcat( folderPath,"/",interp_Descent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    fid = fopen(interp_Descent_source{:});
                    interp_Descent = dlmread(interp_Descent_source{:},',');
                    fclose(fid);
                    
                    %                 interp_Ascent_LB_source     = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Ascent_LB_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    %                 fid = fopen(interp_Ascent_LB_source{:});
                    %                 interp_Ascent_LB = dlmread(interp_Ascent_LB_source{:},',');
                    %                 fclose(fid);
                    %
                    %                 interp_Descent_LB_source    = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Descent_LB_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    %                 fid = fopen(interp_Descent_LB_source{:});
                    %                 interp_Descent_LB = dlmread(interp_Descent_LB_source{:},',');
                    %                 fclose(fid);
                    %
                    %                 interp_Ascent_UB_source     = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Ascent_UB_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    %                 fid = fopen(interp_Ascent_UB_source{:});
                    %                 interp_Ascent_UB = dlmread(interp_Ascent_UB_source{:},',');
                    %                 fclose(fid);
                    %
                    %                 interp_Descent_UB_source    = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Descent_UB_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
                    %                 fid = fopen(interp_Descent_UB_source{:});
                    %                 interp_Descent_UB = dlmread(interp_Descent_UB_source{:},',');
                    %                 fclose(fid);
                    
                    for i = 1:length(interpolatorsEvaluationFields)
                        evolutions(k).population(p).interpolators.ascent.evaluation.(interpolatorsEvaluationFields{i}).value  = interp_Ascent(:,i);
                        evolutions(k).population(p).interpolators.descent.evaluation.(interpolatorsEvaluationFields{i}).value = interp_Descent(:,i);
                    end
                end
            end
        end
        
        for i = 1:length(commonFitnessVectorFieldNames)
            evolutions(k).population(p).fitnessVector.common.(commonFitnessVectorFieldNames{i}).value = fitnessDataPerCurrentGeneration{p,i};
        end
        
        if isfield(evolutions(end).population(end).fitnessVector,'ascent') || isfield(evolutions(end).population(end).fitnessVector,'descent')
            if trajectoryType == 'A'
                for i = 1:length(phaseFitnessVectorFieldNames)
                    evolutions(k).population(p).fitnessVector.ascent.(phaseFitnessVectorFieldNames{i}).value = fitnessDataPerCurrentGeneration{p,i+length(commonFitnessVectorFieldNames)};
                end
                
            end
            if trajectoryType == 'D'
                for i = 1:length(phaseFitnessVectorFieldNames)
                    evolutions(k).population(p).fitnessVector.descent.(phaseFitnessVectorFieldNames{i}).value = fitnessDataPerCurrentGeneration{p,i+length(commonFitnessVectorFieldNames)};
                end
                
            end
            if trajectoryType == 'AD'
                for i = 1:length(phaseFitnessVectorFieldNames)
                    evolutions(k).population(p).fitnessVector.ascent.(phaseFitnessVectorFieldNames{i}).value = fitnessDataPerCurrentGeneration{p,i+length(commonFitnessVectorFieldNames)};
                    evolutions(k).population(p).fitnessVector.descent.(phaseFitnessVectorFieldNames{i}).value = fitnessDataPerCurrentGeneration{p,i+length(commonFitnessVectorFieldNames)+length(phaseFitnessVectorFieldNames)};
                end
            end
        end
        for i = 1:(size(fitnessDataPerCurrentGeneration,2)-length(commonFitnessVectorFieldNames))
            evolutions(k).population(p).fitnessVector.Cost(i,1) = fitnessDataPerCurrentGeneration{p,i};
            %  pop(p).Cost(i,1) = fitnessDataPerCurrentGeneration{p,i};
        end
        
        
        fitnessVectorMagnitude(p) = sqrt(sum(([fitnessDataPerCurrentGeneration{p,:}]).^2, 2));
        evolutions(k).population(p).fitnessVector.commonW.magnitude = fitnessVectorMagnitude(p);
        
        %ppm_population.increment();
        
    end
    
    
    %[pop, idx ] = NonDominatedSorting( pop );
    %evolutions(k).population(1).indices.nonDominatedFront = idx;
    
    ranking = [];
    [~,idxR] = sort([rankDataPerCurrentGeneration{:,1}] + 1,'ascend');
    ranking(:,1) = 1:length([rankDataPerCurrentGeneration{:,1}]);
    ranking(idxR) = ranking;
    
    scaling = [];
    [~,idxS] = sort([rankDataPerCurrentGeneration{:,1}] + 1,'descend');
    scaling(:,1) = 1:length([rankDataPerCurrentGeneration{:,1}]);
    scaling(idxS) = scaling;
    
    for p = 1:evolutions(k).population(1).size.collective
        
        %evolutions(k).population(p).fitnessVector.DominationSet          = pop(p).DominationSet;
        %evolutions(k).population(p).fitnessVector.DominatedCount         = pop(p).DominatedCount;
        %evolutions(k).population(p).fitnessVector.Rank                   = pop(p).Rank;
        evolutions(k).population(p).fitnessVector.commonW.ranking         = ranking(p);
        evolutions(k).population(p).fitnessVector.commonW.scalingFactor   = scaling(p);
    end
    
    evolutions(k).max_tof = max(totalTimeOfFlight);
    
    %ppm_generation.increment();
end



end