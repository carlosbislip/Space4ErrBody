function [ compilation ] = getCompilation( workingFolder, includeDepVar )
%%
online = 0;
matlabScripts = '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/Space4ErrBody_Headers/executables/';
mainpath =      '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/Space4ErrBody_Headers/executables/';
%mainpath = uigetdir;

simulationOutputPath          = strcat(mainpath,'SimulationOutput/');
workingFolderPath             = strcat(simulationOutputPath,workingFolder,'/');
Output_Location               = workingFolderPath;
pop_Location                  = workingFolderPath;
fit_Location                  = workingFolderPath;

Folder_prefix                 = 'OUTPUT*';
figurePath = strcat(workingFolderPath,'figures/');

% Hardcoded coordinates of origin: AMS
lon_i_deg = 4.76416667;
lon_i_rad = deg2rad(lon_i_deg);
% Hardcoded coordinates of final destination: IAD
lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;
validation = 0;

centraTargetCoordinates = [ lat_f_deg, lon_f_deg ];
angularDistanceForTermination = 0.75;
% Hardcoded Starting Epoch
startEpoch = 2458419.95833333;

% Coordinates for Validation case: Re-entry towards Kourou
if contains(pop_Location,'VALIDATION') == 1
    
    lon_i_deg = -106.7;
    lon_i_rad = deg2rad(lon_i_deg);
    lat_f_deg = 5.;
    lon_f_deg = -53;
    validation = 1;
    centraTargetCoordinates = [ lat_f_deg, lon_f_deg ];
    
end

Folders_Containing_Simulations = dir([Output_Location Folder_prefix]);
folderPath = fullfile(Output_Location,{Folders_Containing_Simulations.name}');

%% Scan directories

% option = 1   Try to plot, indiscriminately
% option = 2   Analyze
% option = 3   Analyze and plot

option = 3;


switch option
    
    case 1
        
        textprogressbar('Loading Data.     ')
        %         for p = 1:numel(prop_File_Path_List_prefix)
        %
        %             % Prepare what's needed for the analysis
        %             [evolutions,prop_path,depvar_path,interp_Ascent_path,...
        %                 interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,...
        %                 headingErrorDeadbandBounds_path,...
        %                 alphaMachBounds_path,tof,v_i,gamma_i,chi_i,lat_f,lon_f,pop_path,...
        %                 pop_i,fit_path,fit_i,output] = ...
        %                 ...
        %                 ...
        %                 Path_Prep(option,p,objectiveFunctionCase,...
        %                 prop_File_Path_List_prefix,...
        %                 depvar_File_Path_List_prefix,...
        %                 interp_Ascent_File_Path_List_prefix,...
        %                 interp_Descent_File_Path_List_prefix,...
        %                 DV_mapped_Ascent_File_Path_List_prefix,...
        %                 DV_mapped_Descent_File_Path_List_prefix,...
        %                 headingErrorDeadbandBounds_File_Path_List_prefix,...
        %                 alphaMachBounds_File_Path_List_prefix,...
        %                 Folder_Path_List,pop_file_path_prefix,fit_file_path_prefix);
        %
        %             compilation(p).mainpath = mainpath;
        %             compilation(p).set = char(strcat(extractBetween(output,'evolutions','.mat')));
        %             compilation(p).evolutions = evolutions;
        %             compilation(p).validation = validation;
        %
        %             % textprogressbar(p*100/numel(File_Path_List_prefix))
        %
        %         end
        %
        textprogressbar('    Done. Now plotting/printing.')
        plotsomestuff( compilation , mainpath )
        
    case 2
        
        textprogressbar('Analyzing Simulations.     ')
        %         for p = 1:numel(prop_File_Path_List_prefix)
        % %
        % %             % Prepare what's needed for the analysis
        % %             [evolutions,prop_path,depvar_path,interp_Ascent_path,...
        % %                 interp_Descent_path,DV_mapped_Ascent_path,DV_mapped_Descent_path,...
        % %                 headingErrorDeadbandBounds_path,...
        % %                 alphaMachBounds_path,tof,v_i,gamma_i,chi_i,lat_f,lon_f,pop_path,...
        % %                 pop_i,fit_path,fit_i,output] = ...
        % %                 ...
        % %                 ...
        % %                 Path_Prep(option,p,objectiveFunctionCase,...
        % %                 prop_File_Path_List_prefix,...
        % %                 depvar_File_Path_List_prefix,...
        % %                 interp_Ascent_File_Path_List_prefix,...
        % %                 interp_Descent_File_Path_List_prefix,...
        % %                 DV_mapped_Ascent_File_Path_List_prefix,...
        % %                 DV_mapped_Descent_File_Path_List_prefix,...
        % %                 headingErrorDeadbandBounds_File_Path_List_prefix,...
        % %                 alphaMachBounds_File_Path_List_prefix,...
        % %                 Folder_Path_List,pop_file_path_prefix,fit_file_path_prefix);
        % %             %             if isempty(idx) == 1
        %             %
        %             %                textprogressbar(p*100/numel(File_Path_List_prefix))
        %             %                continue
        %             %
        %             %             end
        %
        %             % Append poplution
        % %             [ evolutions ] = Analyze_Evolution(evolutions,v_i,gamma_i,chi_i,lat_f,lon_f,prop_path,pop_path,pop_i,fit_path,fit_i,lat_f_deg,lon_f_deg);
        %
        %             % Analyze simulations
        % %             [ evolutions ] = ...
        % %                 ...
        % %                 ...
        % %                 Get_Trajectories(...
        % %                 evolutions,...
        % %                 pop_path,...
        % %                 prop_path,...
        % %                 depvar_path,...
        % %                 interp_Ascent_path,...
        % %                 interp_Descent_path,...
        % %                 DV_mapped_Ascent_path,...
        % %                 DV_mapped_Descent_path,...
        % %                 headingErrorDeadbandBounds_path,...
        % %                 alphaMachBounds_path,...
        % %                 v_i,gamma_i,pop_i,lon_i_rad,lat_f_deg,lon_f_deg,startEpoch);
        % %
        % %             save(output,'evolutions');
        % %             %            textprogressbar(p*100/numel(File_Path_List_prefix))
        %
        %         end
        
        % textprogressbar('    Done.')
        
    case 3
        
        pause('on')
        numberOfFolders = numel(folderPath);
        
        if numberOfFolders == 1
            compilation = getCompilationLoop( folderPath{:}, matlabScripts, validation, centraTargetCoordinates, angularDistanceForTermination, workingFolderPath, figurePath);
            if includeDepVar == 1
                compilation.evolutions = getStructuredTrajectoryDepVar( compilation.evolutions, folderPath{:}, compilation.rawData, numberOfFolders, 1, validation );
            end
        else
            parfor p = 1:numberOfFolders
               compilation(p) =  getCompilationLoop( folderPath{p}, matlabScripts, validation, centraTargetCoordinates, angularDistanceForTermination, workingFolderPath, figurePath);
             end
            if includeDepVar == 1
                parfor p = 1:numberOfFolders
                    compilation(p).evolutions = getStructuredTrajectoryDepVar( compilation(p).evolutions, folderPath{p}, compilation(p).rawData, numberOfFolders, p, validation );
                end
            end
            
            poolobj = gcp('nocreate');
            
            delete(poolobj);
        end
        
        %
        %
        %             %for p = 1
        %             % Figure out output name based on what is in folder path
        %             if contains(convertCharsToStrings(folderPath{p}),'OUTPUT') == 1
        %                 delimiterIndices  = strfind(folderPath{p},'T_');
                %                 output = char(strcat('simulationCase_',extractAfter(folderPath{p}(delimiterIndices(end):end),'T_'),'.mat'));
                %             elseif contains(convertCharsToStrings(folderPath{p}),'VALIDATION') == 1
                %                 delimiterIndices  = strfind(folderPath{p},'N_');
                %                 output = char(strcat('simulationCase_',extractAfter(folderPath{p}(delimiterIndices(end):end),'N_'),'.mat'));
                %             end
                %
                %             % Import and massage primary data
                %             rawData = importSimulationData( folderPath{p} );
                %
                %             i = rawData.populationSize;
                %             for k = [ 1 numel(rawData.generationList) ]
                %                 evolutions(k).evolution                                             = nan;
                %                 evolutions(k).population(i).size.collective                         = nan;
                %                 evolutions(k).population(i).size.printed                            = nan;
                %                 evolutions(k).population(i).size.nonPrinted                         = nan;
                %                 evolutions(k).population(i).size.fitness                            = nan;
                %                 evolutions(k).population(i).indices.nonPrinted                      = nan;
                %                 evolutions(k).population(i).indices.nonDominatedFront               = nan;
                %                 evolutions(k).population(i).indices.trajectoryPhaseChange           = nan;
                %                 evolutions(k).population(i).name                                    = nan;
                %                 evolutions(k).population(i).printed                                 = nan;
                %
                %                 evolutions(k).population(i).decisionVector.collective               = nan;
                %                 evolutions(k).population(i).decisionVector.printed                  = nan;
                %                 evolutions(k).population(i).decisionVector.nonPrinted               = nan;
                %
                %                 %if( validation == 0 )
                %                     [ evolutions ] = initializeDecisionVectorFields( rawData, evolutions, k, i );
                %                     [ evolutions ] = initializeFitnessVectorFields( rawData, evolutions, k, i );
                %                 %end
                %                 [ evolutions ] = initializeExtremesAndConstraintsFields( evolutions, k, i );
                %                 [ evolutions ] = initializeDependentVariableFields( evolutions, k, i );
                %
                %
                %                 evolutions(k).max_tof                     = nan;
                %
                %             end
                %             evolutions(1).Common.Bounds.headingError.angularDistanceToGo  = nan;
                %             evolutions(1).Common.Bounds.headingError.upperBound  = nan;
                %             evolutions(1).Common.Bounds.headingError.lowerBound  = nan;
                %             evolutions(1).Common.Bounds.AngleOfAttack.machNumber  = nan;
                %             evolutions(1).Common.Bounds.AngleOfAttack.upperBound  = nan;
                %             evolutions(1).Common.Bounds.AngleOfAttack.lowerBound  = nan;
                %             compilation(p).evolutions = evolutions;
                %             %[ evolutions ] = Get_Trajectories( folderPath{p}, rawData, numberOfFolders, p, validation );
                %
                %             compilation(p).rawData = rawData;
                %             compilation(p).mainpath = matlabScripts;
                %             compilation(p).case = char(strcat(extractBetween(output,'simulationCase','.mat')));
                %             compilation(p).validation = validation;
                %             compilation(p).centraTargetCoordinates = centraTargetCoordinates;
                %             compilation(p).angularDistanceForTermination = angularDistanceForTermination;
                %             compilation(p).workingFolderPath = workingFolderPath;
                %             compilation(p).figurePath = figurePath;
                %
                %
                
                
                %clear evolutions rawData;
                %   savelocation = strcat(workingFolderPath,output);
                %compilationToSave = compilation(p);
                
                %   save(savelocation,'-v7.3', '-struct','compilationToSave')
                % end
                % textprogressbar(p*100/numel(prop_File_Path_List_prefix))
                % ppm.increment();end
                %plotsomestuff( compilation(p) )
       
        
        
    
        %
        %   clear compilation;
        
        
        %    textprogressbar('    Done. Now plotting/printing.')
        
end


end
