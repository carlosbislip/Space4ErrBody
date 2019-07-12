function [ evolutions ] = ...
    ...
    ...
    Get_trajectory(...
    evolutions,...
    Folder_Path_List,...
    pop_Location,...
    pop_path,...
    prop_path,...
    depvar_path,...
    interp_Ascent_path,...
    interp_Descent_path,...
    DV_mapped_Ascent_path,...
    DV_mapped_Descent_path,...
    headingErrorDeadbandBounds_path,...
    alphaMachBounds_path,...
    v_i,gamma_i,pop_i,lon_i_rad,lat_f_deg,lon_f_deg,startEpoch)
%GET_trajectory Summary of this function goes here
%   Detailed explanation goes here
disp('GET_TRAJECTORIES')


if numel(pop_i) > 0
    %populationsize = numel(prop_path)/(numel(pop_i) + 1);
    %a = linspace(populationsize,numel(prop_path),numel(pop_i) + 1);
    %start = [1;a(1:end-1)'+1];
    
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool(4);
    end
    
    strWindowTitle = convertCharsToStrings('Getting trajectory ');
    nEvolutions = (numel(pop_i));
    ppm = ParforProgMon(strWindowTitle, nEvolutions);
else
    %    population = numel(prop_path)/(numel(pop_i) + 1);
    
    %      a = linspace(population,numel(prop_path),numel(pop_i) + 1);
    populationsize =  numel(prop_path);
    
    a = populationsize ;
    start = [1;a(end)];
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool(4);
    end
    nEvolutions = 1;
    strWindowTitle = convertCharsToStrings('Getting trajectory ');
    ppm = ParforProgMon(strWindowTitle, nEvolutions);
    
end

%C = [b a'];
%start = b;
%finish = a';
%evolutions(numel(pop_i)).evolution(population,1).individual = [];


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(4);
end
Folder_prefix = 'OUTPUT*/';

prop_File_prefix = 'propagationHistory/propagationHistory_';
depvar_File_prefix = 'dependentVariables/dependentVariables_';
interp_Ascent_File_prefix = 'evaluatedInterpolatorsAscent/evaluatedInterpolatorsAscent_';
interp_Ascent_LB_File_prefix = 'evaluatedInterpolatorsAscent_LB/evaluatedInterpolatorsAscent_';
interp_Ascent_UB_File_prefix = 'evaluatedInterpolatorsAscent_UB/evaluatedInterpolatorsAscent_';
interp_Descent_File_prefix = 'evaluatedInterpolatorsDescent/evaluatedInterpolatorsDescent_';
interp_Descent_LB_File_prefix = 'evaluatedInterpolatorsDescent_LB/evaluatedInterpolatorsAscent_';
interp_Descent_UB_File_prefix = 'evaluatedInterpolatorsDescent_UB/evaluatedInterpolatorsAscent_';
%monteCarloPopulation_prefix = 'monteCarloPopulation*';

DV_mapped_Ascent_File_prefix = 'map_DV_mapped_Ascent/map_DV_mapped_Ascent_';
DV_mapped_Descent_File_prefix = 'map_DV_mapped_Descent/map_DV_mapped_Descent_';
headingErrorDeadbandBounds_prefix = 'headingErrorDeadBandBounds';
alphaMachBounds_prefix = 'alphaMachBounds';

headingErrorDeadbandBounds_source = { convertStringsToChars(strcat( Folder_Path_List,"/",headingErrorDeadbandBounds_prefix ) ) };
alphaMachBounds_source   = { convertStringsToChars(strcat( Folder_Path_List,"/",alphaMachBounds_prefix ) ) };

fid = fopen(headingErrorDeadbandBounds_source{:});
headingErrorDeadbandBounds = dlmread(headingErrorDeadbandBounds_source{:},',');
fclose(fid);

fid = fopen(alphaMachBounds_source{:});
alphaMachBounds = dlmread(alphaMachBounds_source{:},',');
fclose(fid);


dependentVariableNames        = fieldnames(evolutions(end).population.dependentVariableTimeHistory);
interpolatorsInputDataFields  = fieldnames(evolutions(end).population.interpolators.Ascent.InputData);
interpolatorsEvaluationFields = fieldnames(evolutions(end).population.interpolators.Ascent.Evaluation);

decisionVectorPhaseParameterNames  = fieldnames(evolutions(end).population.decisionVector.parameters.Ascent);
decisionVectorCommonParameterNames = fieldnames(evolutions(end).population.decisionVector.parameters.Common);
interpolatorsInputDataFields       = fieldnames(evolutions(end).population.interpolators.Ascent.InputData);
interpolatorsEvaluationFields      = fieldnames(evolutions(end).population.interpolators.Ascent.Evaluation);

fitnessVectorPhaseNames  = fieldnames(evolutions(end).population.fitnessVector.Ascent);



for k = 1:nEvolutions
    
    evolutions(k).evolution = k-1;
    
    evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo  = headingErrorDeadbandBounds(:,1);
    evolutions(k).Common.Bounds.headingErrorDeadBand.lowerBound  = headingErrorDeadbandBounds(:,2);
    evolutions(k).Common.Bounds.headingErrorDeadBand.upperBound  = headingErrorDeadbandBounds(:,3);
    evolutions(k).Common.Bounds.AngleOfAttack.machNumber  =  alphaMachBounds(:,1);
    evolutions(k).Common.Bounds.AngleOfAttack.lowerBound  =  alphaMachBounds(:,2);
    evolutions(k).Common.Bounds.AngleOfAttack.upperBound  =  alphaMachBounds(:,3);
    
    pp = 1;
    %     prop_source = prop_path(start(k):finish(k),:);
    pop_source           = pop_path(k,:);
    fit_source           = { strrep(pop_source{:},'population','fitness') };
    
    entirePopulation     = loadPopulationFile( pop_source{:} );
    printedIndices       = find([entirePopulation{:,2}] > 0);
    nonPrintedIndices    = find([entirePopulation{:,2}] == 0);
    printedPopulation    = entirePopulation(printedIndices,:);
    nonPrintedPopulation = entirePopulation(nonPrintedIndices,:);
    
    evolutions(k).population.size.collective           = size(entirePopulation,1);
    evolutions(k).population.size.printed              = size(printedPopulation,1);
    evolutions(k).population.size.nonPrinted           = size(nonPrintedPopulation,1);
    evolutions(k).population.indices.printed           = printedIndices;
    evolutions(k).population.indices.nonPrinted        = nonPrintedIndices;
    
    totalTimeOfFlight = nan(evolutions(k).population(1).size.collective,1);
    fitnessVectorMagnitude = totalTimeOfFlight;
    
    % evolutions(k).population.decisionVector.collective     = loadDecisionVectorFromPopulationCellArray( entirePopulation );
    
    entireFitness     = loadFitnessFile( fit_source{:} );
    
    for p = 1:evolutions(k).population.size.collective
        
        evolutions(k).population(p).name    = entirePopulation{p,1};
        evolutions(k).population(p).printed = entirePopulation{p,2};
        n = 1;
        
        for i = 1:length(decisionVectorPhaseParameterNames)
            
            datafields = fieldnames(evolutions(1).population(1).decisionVector.parameters.Ascent.(decisionVectorPhaseParameterNames{i}).data);
            
            for j = 1:(length(datafields))
                evolutions(k).population(p).decisionVector.parameters.Ascent.(decisionVectorPhaseParameterNames{i}).data.(datafields{j}) = entirePopulation{p,n + 3};
                n = n + 1;
            end
            
        end
        
        datafields = fieldnames(evolutions(1).population(1).decisionVector.parameters.Common.(decisionVectorCommonParameterNames{i}).data);
        
        for i = 1:(length(decisionVectorCommonParameterNames) - 3)
            evolutions(k).population(p).decisionVector.parameters.Common.(decisionVectorCommonParameterNames{i}).data.(datafields{1}) = entirePopulation{p,n + 3};
            n = n + 1;
        end
        
        for i = 1:length(decisionVectorPhaseParameterNames)
            
            datafields = fieldnames(evolutions(1).population(1).decisionVector.parameters.Descent.(decisionVectorPhaseParameterNames{i}).data);
            
            for j = 1:length(datafields)
                evolutions(k).population(p).decisionVector.parameters.Descent.(decisionVectorPhaseParameterNames{i}).data.(datafields{j}) = entirePopulation{p,n + 3};
                n = n + 1;
            end
        end
        
        for i = (length(decisionVectorCommonParameterNames) - 2):length(decisionVectorCommonParameterNames)
            evolutions(k).population(p).decisionVector.parameters.Common.(decisionVectorCommonParameterNames{i}).data.one = entirePopulation{p,n + 3};
            n = n + 1;
        end
        
        
        if evolutions(k).population(p).printed == 1
            
            depvar_source            = { convertStringsToChars(strcat( Folder_Path_List,"/",depvar_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
            fid = fopen(depvar_source{:});
            depvar = dlmread(depvar_source{:},',');
            fclose(fid);
            
            for i = 1:(length(dependentVariableNames) - 6)
                evolutions(k).population(p).dependentVariableTimeHistory.(dependentVariableNames{i}) = depvar(:,i);
            end
            
            evolutions(k).population(p).dependentVariableTimeHistory.acc_x      = ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_x ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_x ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_x;
            
            evolutions(k).population(p).dependentVariableTimeHistory.acc_y      = ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_y ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_y ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_y;
            
            evolutions(k).population(p).dependentVariableTimeHistory.acc_z      = ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_z ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_z ...
                + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_z;
            
            evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_M = sqrt(...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_x).^2 ...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_y).^2 ...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_z).^2);
            
            evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_M = sqrt(...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_x).^2 ...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_y).^2 ...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_z).^2);
            
            evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_M = sqrt(...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_x).^2 ...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_y).^2 ...
                + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_z).^2);
            
            
            %distanceToCover = evolutions(k).population(p).dependentVariableTimeHistory.angularDistanceToGo(1);
            %evolutions(k).population(p).dependentVariableTimeHistory.angularDistanceCoveredRatio = ...
            %    evolutions(k).population(p).dependentVariableTimeHistory.angularDistanceTraveled/distanceToCover;
            I_distanceCoveredRatioAboveTerminationRatio = find( evolutions(k).population(p).dependentVariableTimeHistory.angularDistanceCoveredRatio > ...
                evolutions(k).population(p).decisionVector.parameters.Common.terminationDistanceRatio.data.one );
            % XXXXXX = evolutions(k).population(p).dependentVariableTimeHistory.timeOfFlight(I_distanceCoveredRatioAboveTerminationRatio(1));
            I_trajectoryPhaseChange =  find( evolutions(k).population(p).dependentVariableTimeHistory.trajectoryPhase > 0 );
            
            
            evolutions(k).population(p).indices.trajectoryPhaseChange = ...
                [ I_distanceCoveredRatioAboveTerminationRatio(1) I_trajectoryPhaseChange(1) ];
            
            
            %evolutions(k).individuals.gamma_i(p)   = evolutions(k).population(p).individual.flight_path_angle(1);
            %evolutions(k).individuals.chi_i(p)     = evolutions(k).population(p).individual.heading_angle(1);
            %evolutions(k).individuals.lat_f_deg(p) = evolutions(k).population(p).individual.latitude_angle(end);
            %evolutions(k).individuals.lon_f_deg(p) = evolutions(k).population(p).individual.longitude_angle(end);
            totalTimeOfFlight(p) = evolutions(k).population(p).dependentVariableTimeHistory.timeOfFlight(end);
            %evolutions(k).individuals.interp_E_mapped_Ascent(p)  = evolutions(k).population(p).individual.interp_E_mapped_Ascent(end);
            %evolutions(k).individuals.interp_E_mapped_Descent(p)  = evolutions(k).population(p).individual.interp_E_mapped_Descent(end);
            %evolutions(k).fitness.dif_lat(p)   = evolutions(k).population(p).individual.latitude_angle(end) - lat_f_deg;
            %evolutions(k).fitness.dif_lon(p)   = evolutions(k).population(p).individual.longitude_angle(end) - lon_f_deg;
            %evolutions(k).fitness.dif_d_deg(p) = evolutions(k).population(p).individual.distance_to_go(end) - 0.75;
            %evolutions(k).fitness.dif_h(p)     = evolutions(k).population(p).individual.height(end) - 25000;
            evolutions(k).fitness.tof(p)       = evolutions(k).population(p).dependentVariableTimeHistory.timeOfFlight(end);
            
            DV_mapped_Ascent_source  = { convertStringsToChars(strcat( Folder_Path_List,"/",DV_mapped_Ascent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
            DV_mapped_Descent_source = { convertStringsToChars(strcat( Folder_Path_List,"/",DV_mapped_Descent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
            
            
            fid = fopen(DV_mapped_Ascent_source{:});
            DV_mapped_Ascent = dlmread(DV_mapped_Ascent_source{:},',');
            fclose(fid);
            
            fid = fopen(DV_mapped_Descent_source{:});
            DV_mapped_Descent = dlmread(DV_mapped_Descent_source{:},',');
            fclose(fid);
            
            for i = 1:length(interpolatorsInputDataFields)
                evolutions(k).population(p).interpolators.Ascent.InputData.(interpolatorsInputDataFields{i})  = DV_mapped_Ascent(:,i);
                evolutions(k).population(p).interpolators.Descent.InputData.(interpolatorsInputDataFields{i}) = DV_mapped_Descent(:,i);
            end
            
            interp_Ascent_source     = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Ascent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
            interp_Descent_source    = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Descent_File_prefix,evolutions(k).population(p).name,".dat" ) ) };
            
            fid = fopen(interp_Ascent_source{:});
            interp_Ascent = dlmread(interp_Ascent_source{:},',');
            fclose(fid);
            
            fid = fopen(interp_Descent_source{:});
            interp_Descent = dlmread(interp_Descent_source{:},',');
            fclose(fid);
            
            for i = 1:length(interpolatorsEvaluationFields)
                evolutions(k).population(p).interpolators.Ascent.Evaluation.(interpolatorsEvaluationFields{i})  = interp_Ascent(:,i);
                evolutions(k).population(p).interpolators.Descent.Evaluation.(interpolatorsEvaluationFields{i}) = interp_Descent(:,i);
            end
            
        end
        
        
        evolutions(k).population(p).fitnessVector.Common.angularDistanceToGo = entireFitness{p,4};
        
        for i = 1:length(fitnessVectorPhaseNames)
            
            evolutions(k).population(p).fitnessVector.Ascent.(fitnessVectorPhaseNames{i}) = entireFitness{p,i + 4};
            evolutions(k).population(p).fitnessVector.Descent.(fitnessVectorPhaseNames{i}) = entireFitness{p,i + 4 + length(fitnessVectorPhaseNames)};
            
        end
        
        for i = 1:evolutions(k).population(1).size.fitness
            evolutions(k).population(p).fitnessVector.Cost(i,1) = entireFitness{p,i + 3};
            pop(p).Cost(i,1) = entireFitness{p,i + 3};
        end
        
        
        fitnessVectorMagnitude(p) = sqrt(sum(([entireFitness{p,3:end}]).^2, 2));
        evolutions(k).population(p).fitnessVector.Common.magnitude = fitnessVectorMagnitude(p);
        
        
    end
    
    [pop, idx ] = NonDominatedSorting( pop );
    evolutions(k).population(1).indices.nonDominatedFront = idx;
    
    
    [~,idxR] = sort(fitnessVectorMagnitude,'ascend');
    ranking(:,1) = 1:length(fitnessVectorMagnitude);
    ranking(idxR) = ranking;
    
    [~,idxS] = sort(fitnessVectorMagnitude,'descend');
    scaling(:,1) = 1:length(fitnessVectorMagnitude);
    scaling(idxS) = scaling;
    
    for p = 1:evolutions(k).population(1).size.collective
        
        evolutions(k).population(p).fitnessVector.DominationSet          = pop(p).DominationSet;
        evolutions(k).population(p).fitnessVector.DominatedCount         = pop(p).DominatedCount;
        evolutions(k).population(p).fitnessVector.Rank                   = pop(p).Rank;
        evolutions(k).population(p).fitnessVector.Common.ranking         = ranking(p);
        evolutions(k).population(p).fitnessVector.Common.scalingFactor   = scaling(p);        
    end
    
    %evolutions(k).population = population;
    
    evolutions(k).max_tof = max(totalTimeOfFlight);
    
    
    %    decisionVectorPhaseParameterNames = fieldnames(evolutions(1).population.decisionVector.parameters.Ascent);
    %    decisionVectorCommonParameterNames = fieldnames(evolutions(1).population.decisionVector.parameters.Common);
    %
    
    %
    %
    %
    %    if evolutions(k).population.size.printed > 0
    %        evolutions(k).population.decisionVector.printed    = loadDecisionVectorFromPopulationCellArray( printedPopulation );
    %    end
    %     if evolutions(k).population.size.nonPrinted > 0
    %         evolutions(k).population.decisionVector.nonPrinted = loadDecisionVectorFromPopulationCellArray( nonPrintedPopulation );
    %     end
    %
    %
    %     entireFitness     = loadFitnessFile( fit_source{:} );
    %     printedFitness    = entireFitness(printedIndices,:);
    %     nonPrintedFitness = entireFitness(nonPrintedIndices,:);
    %
    %     evolutions(k).population.fitnessVector.collective       = loadFitnessVectorFromFitnessCellArray( entireFitness );
    %     if evolutions(k).population.size.printed > 0
    %        evolutions(k).population.fitnessVector.printed     = loadFitnessVectorFromFitnessCellArray( printedFitness );
    %     end
    %     if evolutions(k).population.size.nonPrinted > 0
    %         evolutions(k).population.fitnessVector.nonPrinted = loadFitnessVectorFromFitnessCellArray( nonPrintedFitness );
    %     end
    %
    %     for p = 1:size(printedPopulation,1)
    %
    %         depvar_source            = { convertStringsToChars(strcat( Folder_Path_List,"/",depvar_File_prefix,printedPopulation{p,1},".dat" ) ) };
    %         interp_Ascent_source     = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Ascent_File_prefix,printedPopulation{p,1},".dat" ) ) };
    %         interp_Descent_source    = { convertStringsToChars(strcat( Folder_Path_List,"/",interp_Descent_File_prefix,printedPopulation{p,1},".dat" ) ) };
    %         DV_mapped_Ascent_source  = { convertStringsToChars(strcat( Folder_Path_List,"/",DV_mapped_Ascent_File_prefix,printedPopulation{p,1},".dat" ) ) };
    %         DV_mapped_Descent_source = { convertStringsToChars(strcat( Folder_Path_List,"/",DV_mapped_Descent_File_prefix,printedPopulation{p,1},".dat" ) ) };
    %
    %         fid = fopen(depvar_source{:});
    %         depvar = dlmread(depvar_source{:},',');
    %         fclose(fid);
    %
    %         fid = fopen(interp_Ascent_source{:});
    %         interp_Ascent = dlmread(interp_Ascent_source{:},',');
    %         fclose(fid);
    %
    %         fid = fopen(interp_Descent_source{:});
    %         interp_Descent = dlmread(interp_Descent_source{:},',');
    %         fclose(fid);
    %
    %         fid = fopen(DV_mapped_Ascent_source{:});
    %         DV_mapped_Ascent = dlmread(DV_mapped_Ascent_source{:},',');
    %         fclose(fid);
    %
    %         fid = fopen(DV_mapped_Descent_source{:});
    %         DV_mapped_Descent = dlmread(DV_mapped_Descent_source{:},',');
    %         fclose(fid);
    %
    %           for i = 1:(length(dependentVariableNames) - 6)
    %             evolutions(k).population(p).dependentVariableTimeHistory.(dependentVariableNames{i}) = depvar(:,i);
    %         end
    %
    %         evolutions(k).population(p).dependentVariableTimeHistory.acc_x      = ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_x ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_x ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_x;
    %
    %         evolutions(k).population(p).dependentVariableTimeHistory.acc_y      = ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_y ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_y ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_y;
    %
    %         evolutions(k).population(p).dependentVariableTimeHistory.acc_z      = ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_z ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_z ...
    %             + evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_z;
    %
    %         evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_M = sqrt(...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_x).^2 ...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_y).^2 ...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_aero_z).^2);
    %
    %         evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_M = sqrt(...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_x).^2 ...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_y).^2 ...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_grav_z).^2);
    %
    %         evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_M = sqrt(...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_x).^2 ...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_y).^2 ...
    %             + (evolutions(k).population(p).dependentVariableTimeHistory.acc_thru_z).^2);
    %
    %
    %
    %         %evolutions(k).individuals.gamma_i(p)   = evolutions(k).population(p).individual.flight_path_angle(1);
    %         %evolutions(k).individuals.chi_i(p)     = evolutions(k).population(p).individual.heading_angle(1);
    %         %evolutions(k).individuals.lat_f_deg(p) = evolutions(k).population(p).individual.latitude_angle(end);
    %         %evolutions(k).individuals.lon_f_deg(p) = evolutions(k).population(p).individual.longitude_angle(end);
    %         evolutions(k).individuals.tof(p)       = evolutions(k).population(p).dependentVariableTimeHistory.timeOfFlight(end);
    %         %evolutions(k).individuals.interp_E_mapped_Ascent(p)  = evolutions(k).population(p).individual.interp_E_mapped_Ascent(end);
    %         %evolutions(k).individuals.interp_E_mapped_Descent(p)  = evolutions(k).population(p).individual.interp_E_mapped_Descent(end);
    %
    %         %evolutions(k).fitness.dif_lat(p)   = evolutions(k).population(p).individual.latitude_angle(end) - lat_f_deg;
    %         %evolutions(k).fitness.dif_lon(p)   = evolutions(k).population(p).individual.longitude_angle(end) - lon_f_deg;
    %         %evolutions(k).fitness.dif_d_deg(p) = evolutions(k).population(p).individual.distance_to_go(end) - 0.75;
    %         %evolutions(k).fitness.dif_h(p)     = evolutions(k).population(p).individual.height(end) - 25000;
    %         evolutions(k).fitness.tof(p)       = evolutions(k).population(p).dependentVariableTimeHistory.timeOfFlight(end);
    %
    %
    %         for i = 1:length(interpolatorsInputDataFields)
    %             evolutions(k).population(p).interpolators.Ascent.InputData.(interpolatorsInputDataFields{i})  = DV_mapped_Ascent(:,i);
    %             evolutions(k).population(p).interpolators.Descent.InputData.(interpolatorsInputDataFields{i}) = DV_mapped_Descent(:,i);
    %         end
    %
    %         for i = 1:length(interpolatorsEvaluationFields)
    %             evolutions(k).population(p).interpolators.Ascent.Evaluation.(interpolatorsEvaluationFields{i})  = interp_Ascent(:,i);
    %             evolutions(k).population(p).interpolators.Descent.Evaluation.(interpolatorsEvaluationFields{i}) = interp_Descent(:,i);
    %         end
    %
    %
    %         pp = pp + 1;
    %         % end
    %     end
    %
    %
    %
    %     if evolutions(k).printedPopulationsize > 0
    %
    %         %     [ aaa , idx1 ] = min(evolutions(k).fitness.dif_d_deg);
    %         %     [ bbb , idx2 ] = min(evolutions(k).fitness.dif_h);
    %         %     [ ccc , idx3 ] = min(evolutions(k).fitness.tof);
    %         [ ddd , idx4 ] = max(evolutions(k).individuals.tof);
    %         %I = [ idx1 idx2 idx3 idx4 ];
    %         I = [ idx4 ];
    %         % criteria = [ {'d_deg'} {'h'} {'tof'} {'max_tof'} ];
    %         criteria = [ {'max_tof'} ];
    %         %  I_12 = intersect(idx1,idx2);
    %         %  I_123 = intersect(I_12,idx3);
    %         %  I_123 = idx1;
    %         evolutions(k).max_tof                     = max(evolutions(k).individuals.tof);
    %         %evolutions(k).max_interp_E_mapped_Ascent  = max(evolutions(k).individuals.interp_E_mapped_Ascent);
    %         %evolutions(k).max_interp_E_mapped_Descent = max(evolutions(k).individuals.interp_E_mapped_Descent);
    %
    %         for i = 1:numel(I)
    %             evolutions(k).best(i).criteria  = criteria(i);
    %             evolutions(k).best(i).index     = I(i);
    %             % evolutions(k).best(i).v_i       = evolutions(k).individuals.v_i(I(i));
    %             %evolutions(k).best(i).gamma_i   = evolutions(k).individuals.gamma_i(I(i));
    %             %evolutions(k).best(i).chi_i     = evolutions(k).individuals.chi_i(I(i));
    %             %evolutions(k).best(i).dif_d_deg = evolutions(k).fitness.dif_d_deg(I(i));
    %             %evolutions(k).best(i).dif_h     = evolutions(k).fitness.dif_h(I(i));
    %             %evolutions(k).best(i).tof       = evolutions(k).fitness.tof(I(i));
    %             evolutions(k).best(i).max_tof   = evolutions(k).individuals.tof(I(i));
    %         end
    %         %}
    %     end
    ppm.increment();
    
    %evolutions(k).filteredIndividuals = filteredIndividuals;
    
end



end