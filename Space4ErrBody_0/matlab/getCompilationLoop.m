function [ compilation ] = getCompilationLoop( folderPath, matlabScripts, validation, centralTargetCoordinates, angularDistanceForTermination, workingFolderPath, figurePath)
  % Figure out output name based on what is in folder path
            if contains(convertCharsToStrings(folderPath),'OUTPUT') == 1
                delimiterIndices  = strfind(folderPath,'T_');
                output = char(strcat('simulationCase_',extractAfter(folderPath(delimiterIndices(end):end),'T_'),'.mat'));
            elseif contains(convertCharsToStrings(folderPath),'VALIDATION') == 1
                delimiterIndices  = strfind(folderPath,'N_');
                output = char(strcat('simulationCase_',extractAfter(folderPath(delimiterIndices(end):end),'N_'),'.mat'));
            end
            
            % Import and massage primary data
            rawData = importSimulationData( folderPath );
            
            i = rawData.populationSize;
            for k = [ 1 numel(rawData.generationList) ]
                evolutions(k).evolution                                             = nan;
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
                
                %if( validation == 0 )
                    [ evolutions ] = initializeDecisionVectorFields( rawData, evolutions, k, i );
                    [ evolutions ] = initializeFitnessVectorFields( rawData, evolutions, k, i );
                %end
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
            compilation.evolutions = evolutions;
            %[ evolutions ] = Get_Trajectories( folderPath{p}, rawData, numberOfFolders, p, validation );
            
            compilation.rawData = rawData;
            compilation.mainpath = matlabScripts;
            compilation.case = char(strcat(extractBetween(output,'simulationCase','.mat')));
            compilation.validation = validation;
            compilation.centralTargetCoordinates = centralTargetCoordinates;
            compilation.angularDistanceForTermination = angularDistanceForTermination;
            compilation.workingFolderPath = workingFolderPath;
            compilation.figurePath = figurePath;
            
            
end