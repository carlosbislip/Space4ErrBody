function [  ] = plotsomestuff( compilation )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
disp('PLOTSOMESTUFF')


mkdir(strcat(compilation(1).workingFolderPath,'figures/'));


lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;


% Coordinates for Validation case: Re-entry towards Kourou
if compilation(1).validation == 1
    
    lon_i_deg = -106.7;
    lon_i_rad = deg2rad(-22.37);
    lon_i_deg = -22.37;
    lat_f_deg = 5;
    lon_f_deg = -53;
    validation = 1;
end




% 
% 
% compilation(1).rawData.populationDataPerGeneration
% compilation(1).rawData.fitnessDataPerGeneration
% compilation(1).rawData.extremesAndConstraintsDataPerGeneration
% 
% %%
% I_1 = find([compilation(1).rawData.extremesAndConstraintsDataPerGeneration{:,2}]'==1);
% 
% derpExtremes = compilation(1).rawData.extremesAndConstraintsDataPerGeneration(I_1,:);
% derpPopulation = compilation(1).rawData.populationDataPerGeneration(I_1,:);
% derpFitness = compilation(1).rawData.fitnessDataPerGeneration(I_1,:);
% 

%%

%plotObjectiveComparisonsForSubplots( compilation );
%plotInterpolators( compilation );

% plotTimeHistoriesOfConstrainedParameters( compilation );

%
plotFitnessComparisonsForDiscussion( compilation );
plotFitnessComparisons( compilation );
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




%%
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Pareto Front per Evolution:_{ }',' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        
        grid on
        grid minor
        
        [h,ax,BigAx] = gplotmatrix(compilation(p).evolutions(k).population(1).pureFitness);
        
        %
        AllPos   = get(ax(:),'Position');
        AllXLim  = get(ax(:),'XLim');
        AllXTick = get(ax(:),'XTick');
        AllXTickLabel = get(ax(:),'XTickLabel');
        AllYLim  = get(ax(:),'YLim');
        AllYTick = get(ax(:),'YTick');
        
        AllPos   = vertcat(AllPos{:});
        AllXLim  = vertcat(AllXLim{:});
        %AllXTick = vertcat(AllXTick{:});
        AllYLim  = vertcat(AllYLim{:});
        % AllYTick = vertcat(AllYTick{:});
        
        %NewPos  = [AllPos(:,1)+.05 AllPos(:,2)+.05 AllPos(:,3)-.1 AllPos(:,4)-.1];
        NewPos  = [AllPos(:,1)+.005 AllPos(:,2)+.005 AllPos(:,3)-.01 AllPos(:,4)-.01];
        NewXLim = AllXLim;
        NewYLim = AllYLim;
        %NewXTick = AllXTick(:,1);
        nObj = compilation(p).evolutions(1).population(1).size.fitness+1;
        
        column1 = 1:nObj;
        for i = column1
            NewXLim(i,:) = [0 1000];
            NewXTick(i,:) = 0:200:1000;
            
            if i == column1(end) - 1
                %label = cellstr(num2str((0:200:1000)'));
                NewXTickLabel(i,:) = {'0';'200';'400';'600';'800';'1000'};
            end
            
        end
        derp = 1:(nObj):nObj*(nObj);
        row1 = [1 (derp(2:end-1))];
        %row1 = find(NewYLim(:,1) == 0) + 1;
        
        
        for i = row1
            NewYLim(i,:) = [0 1000];
            
        end
        
        
        column = column1;
        row = row1;
        for i = 2:nObj
            column = column1 + nObj*(i - 1);
            for j = column
                NewXLim(j,:) = [0 10];
                NewXTick(j,:) = 0:2:10;
                
                if j == column(end) - 1
                    %label = cellstr(num2str((0:200:1000)'));
                    NewXTickLabel(j,:) = {'0';'2';'4';'6';'8';'10'};
                end
                
                
                
            end
            
            row = row + i - 1;
            for j = row
                NewYLim(j,:) = [0 10];
            end
        end
        
        
        
        %// Update the plot
        %         for kk = 1:numel(ax)
        %             axes(ax(kk))
        %
        %             set(ax(kk),'Position',NewPos(kk,:))
        %             set(ax(kk),'XLim',NewXLim(kk,:))
        %             set(ax(kk),'YLim',NewYLim(kk,:))
        %             set(ax(kk),'XTick',NewXTick(kk,:))
        %             set(ax(kk),'XTickLabel',NewXTickLabel(kk,:))
        %
        %
        %         end
        
        %  set(ax(5),'XTickLabel',NewXTickLabel)
        %
        %          ax(1,1).XLim = ([0  1000]);
        %          ax(1,1).YLim = ([0  1000]);
        %                     ax(1,1).XTick = 0:200:1000;
        %             ax(1,1).YTick = 0:200:1000;
        %
        %         for i = 2:compilation(p).evolutions(1).population(1).size.fitness
        %             ax(2,i).XLim = ([0  10]);
        %             ax(i,2).YLim = ([0  10]);
        %             ax(2,i).XTick = 0:2:10;
        %             ax(i,2).YTick = 0:2:10;
        %         end
        %
        
        
        
        
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'gridPlot_Evolution',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %        close(fig_num);
        
    end
end

%% Printed Distribution per Evolution
for p = 1:numel(compilation)
    
    fig_num = p*100;
    figure(fig_num)
    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
    set (gca,'Fontsize',15)
    title(strcat('Printed Distribution per Evolution:_{ }',' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
    ylim([0 compilation(p).evolutions(1).population(1).size.collective])
    % max_tof = max([compilation(p).evolutions.max_tof]);
    %max_tof =1400;
    %xlim([0 max_tof])
    xlabel('Evolution') % x-axis label
    ylabel('Printed Distribution') % y-axis label
    set(gca,'YTick', 0:compilation(p).evolutions(1).population(1).size.collective/10:compilation(p).evolutions(1).population(1).size.collective);
    set(gca,'XTick', 0:1:numel(compilation(p).evolutions)-1);
    
    hold on
    grid on
    
    for k = 1:numel(compilation(p).evolutions)
        popMatrix(k,:) = [ compilation(p).evolutions(k).population(1).size.printed compilation(p).evolutions(k).population(1).size.nonPrinted ];
        
    end
    
    
    X = linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)) - 1;
    
    if X == 0
        bar(0,sum(popMatrix),'FaceColor',[0.8500, 0.3250, 0.0980])
        bar(0,popMatrix(1),'FaceColor',[0, 0.4470, 0.7410])
        labels = {'Printed Population','Non-Printed Population'};
        legend(labels{end:-1:1});
    else
        bar( X, popMatrix,'stacked')
        legend('Printed Population','Non-Printed Population')
    end
    
    % a = axis;
    % axis([a(1) a(2)-0.8 a(3:4)]);
    
    %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).figurePath,...
        'printedDistribution_Evolution_',...
        '_Set',...
        convertCharsToStrings(compilation(p).set),...
        '.png'),...
        'png');
    %        close(fig_num);
end

%% PARETO!?????
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Pareto Front per Evolution:_{ }',' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 compilation(p).evolutions(1).population(1).size.collective])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        %max_tof =1400;
        %xlim([0 max_tof])
        xlabel('Objective X') % x-axis label
        ylabel('Objective Y') % y-axis label
        zlabel('Objective Z') % y-axis label
        %set(gca,'YTick', 0:compilation(p).evolutions(1).population(1).size.collective/10:compilation(p).evolutions(1).population(1).size.collective);
        %set(gca,'XTick', 0:1:numel(compilation(p).evolutions)-1);
        
        hold on
        grid on
        
        
        
        
        
        
        for i = 1:length(compilation(p).evolutions(k).population(1).indices.nonDominatedFront)
            
            indices = [ compilation(p).evolutions(k).population(1).indices.nonDominatedFront{1,i} ];
            
            for j = 1:length(indices)
                XYZ(j,:) = [ compilation(p).evolutions(k).population(indices(j)).fitnessVector.Common.angularDistanceToGo ...
                    compilation(p).evolutions(k).population(indices(j)).fitnessVector.Ascent.basicDynamics ...
                    compilation(p).evolutions(k).population(indices(j)).fitnessVector.Descent.basicDynamics ];
            end
            scatter3( XYZ(:,1),XYZ(:,2),XYZ(:,3))
            
            
        end
        
        view(3)
        
        
        % a = axis;
        % axis([a(1) a(2)-0.8 a(3:4)]);
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'paretoFront_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %        close(fig_num);
    end
    
end




%}
end
