function [  ] = plot_MC_w_MagnifiedSection( compilation, workingFolder )
disp('   Monte Carlo with Magnified Section')

%%
k = 1;
extremesAndConstraintsFieldNames  = fieldnames(compilation(end).evolutions(end).population(end).extremesAndConstraints);
extremesVectorLength = length(extremesAndConstraintsFieldNames);

X = nan(compilation(end).rawData.populationSize,extremesVectorLength,numel(compilation));

for pp = 1:numel(compilation)
    X(:,:,pp) = cell2mat(compilation(pp).rawData.extremesAndConstraintsDataPerGeneration(:,7:end,k));
    
    
    %
    %     for ii = 1:compilation(p).rawData.populationSize
    %         for j = 1:length(extremesAndConstraintsFieldNames)
    %             X(ii,j,pp) = compilation(pp).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
    %         end
    %     end
    
end
combinations = ...
    [1 8];
%     1 9;...
%     1 10;...
%     1 11;...
%     1 12;...
%     1 21;...
%     1 22;...
%     1 23;...
%     1 24];



's';
%% Point Cloud with zoom Inset

for pp = 1:numel(compilation)
    for kk = 1:size(combinations,1)
        
        zoom = [{''};{'_zoom'};{'_zoom1'}];
        
        for kkk = 1%:numel(zoom)
            fig_num = 100 + 3566000 + kk;
            figure(fig_num)
            
            Xx = combinations(kk,1);
            Xy = combinations(kk,2);
            field_x = extremesAndConstraintsFieldNames{Xx};
            field_x_struct = compilation(end).evolutions(end).population(end).extremesAndConstraints.(field_x);
            variableLabel_x = field_x_struct.variableLabel;
            figureSaveName_x = field_x_struct.figureSaveNameContent;
            units_x = field_x_struct.units;
            limits_x = field_x_struct.limits;
            scalingFactor_x = field_x_struct.scalingFactor;
            tick_x = field_x_struct.tick;
            
            field_y = extremesAndConstraintsFieldNames{Xy};
            field_y_struct = compilation(end).evolutions(end).population(end).extremesAndConstraints.(field_y);
            variableLabel_y = field_y_struct.variableLabel;
            figureSaveName_y = field_y_struct.figureSaveNameContent;
            units_y = field_y_struct.units;
            limits_y = field_y_struct.limits;
            scalingFactor_y = field_y_struct.scalingFactor;
            tick_y = field_y_struct.tick;
            
            %t = tiledlayout(1,3);
            
            m = {'square','+', 'diamond', 'o', '*','square', '<','.', 'x', 'v', '^', '>', 'pentagram'};
            
            fontSize = 15;
            
            %nexttile;
            dscatter(X(:,Xx,pp)/scalingFactor_x,X(:,Xy,pp)/scalingFactor_y,'MARKER',m{1})
            ax1 = gca;
            xlim(ax1,limits_x)
            
            if kkk == 1
                maxValue = 25*ceil(max(max(X(:,Xy,pp)/scalingFactor_y))/25);
                limits_y = [0 maxValue];
                tick_y = maxValue/10;
            elseif kkk == 2
            elseif kkk == 3
                if strcmp(field_y,'maximumBodyFrameMechanicalLoad')
                    limits_y = [0 3];
                    tick_y = 3/10;
                end
            end
            
            
            ylim(ax1,limits_y)
            grid(ax1,'on')
            set(ax1,'Fontsize',fontSize)
            set(ax1,'XTick', round(limits_x(1):tick_x:limits_x(2)));
            set(ax1,'YTick', limits_y(1):tick_y:limits_y(2));
            set(ax1,'TickLabelInterpreter','latex')
            colormap(gca,'jet')
            
            cbh = colorbar;
            cbh.Ticks = caxis( gca ) ;
            cbh.TickLabels = [{'Low'} {'High'}];
            cbh.TickLabelInterpreter = 'latex';
            ylabel(cbh,'Distribution Density','Interpreter','latex','Fontsize',fontSize)
            
            colobarLabelPos = get(get(cbh,'YLabel'),'Position');
            colobarLabelPos(1) = colobarLabelPos(1) - 1;
            set(get(cbh,'YLabel'),'Position',colobarLabelPos)
            
            xlabel(ax1,string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex','Fontsize',fontSize) % x-axis label
            ylabel(ax1,string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex','Fontsize',fontSize) % x-axis label
            
            
            if compilation(pp).rawData.trajectoryType == 'A' && compilation(pp).rawData.objectiveFunctionCase == 'A'
                ax2 = axes('Position',[.15 .7 .2 .2]);
                box(ax2,'on');
                dscatter(X(:,Xx,pp)/scalingFactor_x,X(:,Xy,pp)/scalingFactor_y,'MARKER',m{1})
                grid(ax2, 'on');
                
                limits_x_zoom = [53 55.25];
                limits_y_zoom = [0 2];
                
                xlim(ax2,limits_x_zoom)
                ylim(ax2,limits_y_zoom)
                set(ax2,'Fontsize',fontSize)
                set(ax2,'XTick', limits_x_zoom(1):1.125:limits_x_zoom(2));
                set(ax2,'YTick', limits_y_zoom(1):1:limits_y_zoom(2));
                ax2.YAxisLocation = 'right';
                            set(ax2,'TickLabelInterpreter','latex')

            elseif strcmp(workingFolder,'MC_MO_ASCENT')
            elseif strcmp(workingFolder,'MC_SO_COUPLED')
            elseif compilation(pp).rawData.trajectoryType == 'AD' && compilation(pp).rawData.objectiveFunctionCase == 'I'
                
                if compilation(pp).rawData.nodes == 5
                    
                    ax2 = axes('Position',[.6 .7 .2 .2]);
                    box(ax2,'on');
                    dscatter(X(:,Xx,1)/scalingFactor_x,X(:,Xy,1)/scalingFactor_y,'MARKER',m{1})
                    grid(ax2, 'on');
                    
                    limits_x_zoom = [45 55.25];
                    limits_y_zoom = [0 15];
                    
                    xlim(ax2,limits_x_zoom)
                    ylim(ax2,limits_y_zoom)
                    set(ax2,'Fontsize',fontSize)
                    set(ax2,'XTick', limits_x_zoom(1):5:limits_x_zoom(2));
                    set(ax2,'YTick', limits_y_zoom(1):5:limits_y_zoom(2));
                                set(ax2,'TickLabelInterpreter','latex')

                elseif compilation(pp).rawData.nodes == 6
                    
                     ax2 = axes('Position',[.6 .7 .2 .2]);
                    box(ax2,'on');
                    dscatter(X(:,Xx,1)/scalingFactor_x,X(:,Xy,1)/scalingFactor_y,'MARKER',m{1})
                    grid(ax2, 'on');
                    
                    limits_x_zoom = [45 55.25];
                    limits_y_zoom = [0 15];
                    
                    xlim(ax2,limits_x_zoom)
                    ylim(ax2,limits_y_zoom)
                    set(ax2,'Fontsize',fontSize)
                    set(ax2,'XTick', limits_x_zoom(1):5:limits_x_zoom(2));
                    set(ax2,'YTick', limits_y_zoom(1):5:limits_y_zoom(2));
                                set(ax2,'TickLabelInterpreter','latex')
                    
                    
                elseif compilation(pp).rawData.nodes == 7
                 ax2 = axes('Position',[.6 .7 .2 .2]);
                    box(ax2,'on');
                    dscatter(X(:,Xx,1)/scalingFactor_x,X(:,Xy,1)/scalingFactor_y,'MARKER',m{1})
                    grid(ax2, 'on');
                    
                    limits_x_zoom = [45 55.25];
                    limits_y_zoom = [0 15];
                    
                    xlim(ax2,limits_x_zoom)
                    ylim(ax2,limits_y_zoom)
                    set(ax2,'Fontsize',fontSize)
                    set(ax2,'XTick', limits_x_zoom(1):5:limits_x_zoom(2));
                    set(ax2,'YTick', limits_y_zoom(1):5:limits_y_zoom(2));
                                set(ax2,'TickLabelInterpreter','latex')
                end
                
            end
            
            colormap(gca,'jet')
            
            hold off
            
            seedsAndInitializers = strcat(strcat(num2str(compilation(pp).rawData.nodes),{'_'},num2str(compilation(pp).rawData.seedInitializer)));%,{'_'},...
            figureSaveName = string(strcat(...
                compilation(1).figurePath,...
                'objectiveComparisons_MC_w_Magnified_Section_',...
                strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
                {'_'},seedsAndInitializers,{'_'},compilation(pp).rawData.trajectoryType,{'_'},compilation(pp).rawData.objectiveFunctionCase,'.png'));
            
            saveas(figure(fig_num),figureSaveName,'png');
            pause(10)
            close(fig_num);
        end
    end
end



end