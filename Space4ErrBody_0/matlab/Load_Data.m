function [evolutions,output] = Load_Data(option,p,validation,simulationDetails,Folder_Path_List )
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
        
       
        i = simulationDetails.populationSize;
        for k = [ 1 simulationDetails.maxGenerations ]
            evolutions(k).evolution                                          = nan;
            
            evolutions(k).population(i).size.collective                         = nan;
            evolutions(k).population(i).size.printed                            = nan;
            evolutions(k).population(i).size.nonPrinted                         = nan;
            evolutions(k).population(i).size.fitness                            = nan;
            evolutions(k).population(i).indices.nonPrinted                      = nan;
            evolutions(k).population(i).indices.nonDominatedFront               = nan;
            evolutions(k).population(i).indices.trajectoryPhaseChange           = nan;
            evolutions(k).population(i).name                                    = nan;
            evolutions(k).population(i).printed                                 = nan;
            
            evolutions(k).population(i).decisionVector.collective               = nan;
            evolutions(k).population(i).decisionVector.printed                  = nan;
            evolutions(k).population(i).decisionVector.nonPrinted               = nan;
            
            if( validation == 0 )
                [ evolutions ] = initializeDecisionVectorCommonParameterFields( evolutions, k, i );
                [ evolutions ] = initializeDecisionVectorNodalParameterFields( simulationDetails, evolutions, k, i );
                [ evolutions ] = initializeFitnessVectorFields( simulationDetails, evolutions, k, i );
            end
            [ evolutions ] = initializeExtremesAndConstraintsFields( evolutions, k, i );
            [ evolutions ] = initializeDependentVariableFields( evolutions, k, i );
            
            
            evolutions(k).max_tof                     = nan;
            
        end
        evolutions(1).Common.Bounds.headingError.angularDistanceToGo  = nan;
        evolutions(1).Common.Bounds.headingError.upperBound  = nan;
        evolutions(1).Common.Bounds.headingError.lowerBound  = nan;
        evolutions(1).Common.Bounds.AngleOfAttack.machNumber  = nan;
        evolutions(1).Common.Bounds.AngleOfAttack.upperBound  = nan;
        evolutions(1).Common.Bounds.AngleOfAttack.lowerBound  = nan;
end


end

