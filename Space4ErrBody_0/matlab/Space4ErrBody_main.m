
clearvars
close all
clc
dbstop if error

format long
%hold on
%[S_x, S_y, S_z] = sphere(100);
%surf(S_x*r, S_y*r, S_z*r);`
%p = plotearth('Maptype','bluemarble','NEOMap',...
%    '/Users/bislip/Documents/MATLAB/PlotEarth/PlotEarth/BlueMarble.png',...
%    'SampleStep',1,'FullColor',true,'Shape','spherical');


%%
clearvars
close all
clc
dbstop if error

format long

%workingFolder = 'SINGLE_OBJECTIVE_ASCENT';
%workingFolder = 'SINGLE_OBJECTIVE_ASCENT_TOP_INDIVIDUALS_FROM_MO_MC';
workingFolder = 'SINGLE_OBJECTIVE_ASCENT_TOP_INDIVIDUALS_FROM_SO_MC';

includeDepVar = 0;

compilation = getCompilation( workingFolder, includeDepVar );
mkdir(strcat(compilation(1).workingFolderPath,'figures/'));

plot_SO_FitnessComparisons_ForDiscussion( compilation );

plot_SO_ObjectiveComparisons_ForDiscussion( compilation );

%% Monte Carlo Plots
clearvars
close all
clc
dbstop if error
format long

workingFolder = 'MC_SO_ASCENT';
%workingFolder = 'MC_MO_ASCENT';
%workingFolder = 'MC_SO_COUPLED';
workingFolder = 'MC_MO_COUPLED';
workingFolder = 'MC';

includeDepVar = 0;

compilation = getCompilation( workingFolder, includeDepVar );
mkdir(strcat(compilation(1).workingFolderPath,'figures/'));

for p = 1:numel(compilation)
    plot_MC_w_MagnifiedSection( compilation(p), workingFolder )
    plot_MC_BivariateHistogram( compilation(p) )
end


%% Top Individuals with Monte Carlo Underlay
clearvars
close all
clc
dbstop if error
format long

workingFolder = 'MC_SO_ASCENT';
%workingFolder = 'MC_MO_ASCENT';
%workingFolder = 'MC_SO_COUPLED';
workingFolder = 'MC_MO_COUPLED';

includeDepVar = 0;

compilation = getCompilation( workingFolder, includeDepVar );
mkdir(strcat(compilation(1).workingFolderPath,'figures/'));

for p = 1:numel(compilation)
    plotTopIndividuals_w_MC_Underlay( compilation )
end




%% Validation Figures
clearvars
close all
clc
dbstop if error
format long

workingFolder = 'VALIDATION_HORUS';
includeDepVar = 1;

compilation = getCompilation( workingFolder, includeDepVar );
mkdir(strcat(compilation(1).workingFolderPath,'figures/'));


plot_Validation_Figures( compilation );

plot_Validation_Trajectories( compilation );




%%


plotsomestuff( getCompilation( workingFolder, includeDepVar ) );

plotInterpolatorsForDiscussion( compilation )

plotObjectiveComparisonsForSubplots( getCompilation( workingFolder, includeDepVar ) )

plotTimeHistoriesOfConstrainedParameters( getCompilation( workingFolder, includeDepVar ) )

plotFitnessComparisonsForDiscussion( getCompilation( workingFolder, includeDepVar ) )

plotObjectiveComparisonsForDiscussion( compilation )

plotTopIndividualsForDiscussion( compilation )


%plotObjectiveComparisonsForSubplots( compilation );
%plotInterpolators( compilation );

% plotTimeHistoriesOfConstrainedParameters( compilation );

%
plotFitnessComparisonsForDiscussion(  getCompilation( workingFolder, includeDepVar ) );
plotFitnessComparisons(  getCompilation( workingFolder, includeDepVar ) );
%plotObjectiveComparisons( compilation );
%%
% plotTimeHistories( compilation );
%
%
% %%
% plotFitnessComparisons( compilation );
%
% plotObjectiveComparisons( compilation );
%

%{
%%
plotDecisionVector( compilation );

%%
plotDecisionVector( compilation );

%%
plotTrajectories( compilation );

%%
plotInterpolators( compilation );

%%
plotTimeHistories_AeroAngles( compilation );

%%
plotTimeHistories_AeroCoefficients( compilation );

%%
plotTimeHistories_Thermo( compilation );

%%
plotTimeHistories_Mechanical( compilation );

%%
plotTimeHistories_EngineOps( compilation );

%%
plotTimeHistories_SpatialAwareness( compilation );

%%
plotCompoundRelations( compilation );
%}