function [ simulationData ] = importSimulationData( folderPath )
disp(string(strcat({'     Importing simulaton data'})))


parenthesisIndices = strfind(folderPath,'/');
workingSubFolder = folderPath(parenthesisIndices(end)+1:end);
underscore_indices  = strfind(workingSubFolder,'_');
temp = workingSubFolder;
disp(string(strcat({'     Current folder: '},workingSubFolder)))

[caseNumber, status] = str2num(temp(underscore_indices(1)+1:underscore_indices(2)-1));

if status == 0
    [caseNumber, status] = str2num(temp(underscore_indices(3)+1:underscore_indices(4)-1));
end

if isempty(caseNumber) == 1 && status == false
    [caseNumber, status] = str2num(temp(underscore_indices(4)+1:underscore_indices(5)-1));
end
%caseNumber = str2double(temp(underscore_indices(1)+1:underscore_indices(2)-1));

caseMatrix = readcell( '/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/Space4ErrBody_Headers/executables/inputCaseMatrix1.txt','Delimiter',{'\t'} );
caseDetails = caseMatrix(caseNumber+1,:);
simulationData.caseDetails = caseDetails;
simulationData.NodalParametersFileName = caseDetails{14};
simulationData.seedInitializer = caseDetails{62};
simulationData.populationSize = caseDetails{59};
nodalParameters_underscore_indices  = strfind(caseDetails{14},'_');
nodalParameters_end_of_filename  = strfind(caseDetails{14},'.txt');

disp(string(strcat({'       Determine trajectory type'})))
trajectoryType = caseDetails{13};
disp(string(strcat({'           Trajectory type: '},trajectoryType)))
simulationData.trajectoryType = trajectoryType;
disp(string(strcat({'       Determine number of nodes'})))
simulationData.nodes = str2double(caseDetails{14}(nodalParameters_underscore_indices(2)+1));
disp(string(strcat({'           Number of nodes: '},num2str(simulationData.nodes))))

simulationData.nodalParameters = caseDetails{14}(nodalParameters_underscore_indices(2)+2:nodalParameters_end_of_filename-1);
simulationData.numberOfPhases = numel(trajectoryType);

disp(string(strcat({'       Determine objective function case'})))
simulationData.objectiveFunctionCase = char(caseDetails{64});
disp(string(strcat({'           Objective function case: '},simulationData.objectiveFunctionCase)))

disp(string(strcat({'       Determine optimizing algorithm'})))
optimizerIndex = caseDetails{66}+1;
optimizingAlgorithms = [{'NSGA-II'} {'MOEAD/DE'} {'IHS'}];
simulationData.optimizingAlgorithmIndex = optimizerIndex;
simulationData.optimizingAlgorithm = optimizingAlgorithms{optimizerIndex};
disp(string(strcat({'           Algorithm used:  '},optimizingAlgorithms{optimizerIndex})))

%%
disp(string(strcat({'       Read in Population Data'})))
simulationData.populationSource = strcat(folderPath,'/populationHistory.dat');
populationData = readcell( simulationData.populationSource,'Delimiter',{','} );
simulationData.populationData = populationData;
%%
disp(string(strcat({'       Read in Fitness Data'})))
simulationData.fitnessSource = strcat(folderPath,'/fitnessHistory.dat');
fitnessData = readcell( simulationData.fitnessSource,'Delimiter',{','} );
simulationData.fitnessData = fitnessData;

disp(string(strcat({'       Read in Extremes Data'})))
simulationData.extremesAndConstraintsSource = strcat(folderPath,'/extremesAndConstraintsHistory.dat');
extremesAndConstraintsData = readcell( simulationData.extremesAndConstraintsSource,'Delimiter',{','} );
simulationData.extremesAndConstraintsData = extremesAndConstraintsData;

%% Extract top individuals
topIndividuals = caseDetails{34};
disp(string(strcat({'       Extract Data of Top Individuals'})))
I = find([populationData{:,3}] < topIndividuals);
simulationData.topIndividuals.size = topIndividuals;
simulationData.topIndividuals.population = populationData(I,:);
simulationData.topIndividuals.fitness = fitnessData(I,:);
simulationData.topIndividuals.extremesAndConstraints = extremesAndConstraintsData(I,:);


writecell(simulationData.topIndividuals.population,string(strcat(folderPath,{'/population'},{'_top_'},num2str(topIndividuals),'.dat')))
%writecell(simulationData.topIndividuals.fitness,string(strcat(folderPath,{'/fitnessHistory'},{'_top_'},num2str(topIndividuals),'.dat')))
%writecell(simulationData.topIndividuals.extremesAndConstraints,string(strcat(folderPath,{'/extremesAndConstraintsHistory'},{'_top_'},num2str(topIndividuals),'.dat')))


%%

simulationData.generationList = unique([populationData{:,5}]);
simulationData.maxGenerations = simulationData.generationList(end) + 1;
nlay = numel(simulationData.generationList);
% if optimizerIndex == 0 || optimizerIndex == 1
%     simulationData.maxGenerations = max( [populationData{:,5}] ) + 1;
%     nlay = simulationData.maxGenerations;
% else
%     simulationData.maxGenerations = simulationData.generationList(end) + 1;
%     nlay = 1 + simulationData.maxGenerations/simulationData.populationSize;
% end

[r,c] = size(populationData);
disp(string(strcat({'       Massage Population Data'})))
simulationData.populationDataPerGeneration = permute(reshape((populationData)',[c,r/nlay,nlay]),[2,1,3]);

[r,c] = size(fitnessData);
disp(string(strcat({'       Massage Fitness Data'})))
simulationData.fitnessDataPerGeneration = permute(reshape((fitnessData)',[c,r/nlay,nlay]),[2,1,3]);

[r,c] = size(extremesAndConstraintsData);
disp(string(strcat({'       Massage Extremes Data'})))
simulationData.extremesAndConstraintsDataPerGeneration = permute(reshape((extremesAndConstraintsData)',[c,r/nlay,nlay]),[2,1,3]);

simulationData.fitnessMagnitudeAveragesOutputHistorySource = strcat(folderPath,'/fitnessMagnitudeAveragesOutputHistory.dat' );
fid = fopen(simulationData.fitnessMagnitudeAveragesOutputHistorySource);
disp(string(strcat({'       Read in Fitness Magnitude Data'})))
fitnessMagnitudeAverageHistory = dlmread( simulationData.fitnessMagnitudeAveragesOutputHistorySource,',');
fclose(fid);

simulationData.fitnessMagnitudeAverageHistory = fitnessMagnitudeAverageHistory;

simulationData.headingErrorDeadbandBoundsSource = strcat(folderPath,'/headingErrorDeadbandBounds.dat' );
fid = fopen(simulationData.headingErrorDeadbandBoundsSource);
disp(string(strcat({'       Read in Heading Error Deadband Data'})))
simulationData.headingErrorDeadbandBounds = dlmread( simulationData.headingErrorDeadbandBoundsSource,',');
fclose(fid);

simulationData.alphaMachBoundsSource = strcat(folderPath,'/alphaMachBounds.dat' );
fid = fopen(simulationData.alphaMachBoundsSource);
disp(string(strcat({'       Read in Alpha-Mach Bounds Data'})))
simulationData.alphaMachBounds = dlmread( simulationData.alphaMachBoundsSource,',' );
fclose(fid);

simulationData.targetCircleCoordinatesSource = strcat(folderPath,'/targetCircleCoordinates.dat' );
fid = fopen(simulationData.targetCircleCoordinatesSource);
disp(string(strcat({'       Read in Target Circle Coordinates Data'})))
simulationData.targetCircleCoordinates = dlmread( simulationData.targetCircleCoordinatesSource,',' );
fclose(fid);



end