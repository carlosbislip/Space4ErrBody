function [  ] = plotObjectiveComparisonsForSubplots( compilation )
disp('   Objective Comparisons for Subplots')

%%
k = 1;
extremesAndConstraintsFieldNames  = fieldnames(compilation(end).evolutions(end).population(end).extremesAndConstraints);

X = nan(numel(compilation(end).evolutions(end).population(end).size.collective,length(extremesAndConstraintsFieldNames)),numel(compilation));

for pp = 1:numel(compilation)
    for ii = 1:compilation(pp).evolutions(k).population(1).size.collective
        for j = 1:length(extremesAndConstraintsFieldNames)
            X(ii,j,pp) = compilation(pp).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
        end
    end
end

%%

combinations = ...
    [1 8;...
    1 9;...
    1 10;...
    1 11;...
    1 12;...
    1 21;...
    1 22;...
    1 23;...
    1 24];


for kk = 1:size(combinations,1)
    
    zoom = [{''};{'_zoom'};{'_zoom1'}];
    
    for kkk = 1:numel(zoom)
        fig_num = 100 + 3566000 + kk;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,750,420])
        
        Xx = combinations(kk,1);
        Xy = combinations(kk,2);
        field_x = extremesAndConstraintsFieldNames{Xx};
        variableLabel_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).variableLabel;
        figureSaveName_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).figureSaveNameContent;
        units_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).units;
        limits_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).limits;
        tick_x = 2*compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).tick;
        scalingFactor_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).scalingFactor;
        
        field_y = extremesAndConstraintsFieldNames{Xy};
        variableLabel_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).variableLabel;
        figureSaveName_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).figureSaveNameContent;
        units_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).units;
        limits_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).limits;
        tick_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).tick;
        scalingFactor_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).scalingFactor;
        
        t = tiledlayout(1,3);
        
        m = {'square','+', 'diamond', 'o', '*','square', '<','.', 'x', 'v', '^', '>', 'pentagram'};
        
        fontSize = 15;
        
        nexttile;
        dscatter(X(:,Xx,1)/scalingFactor_x,X(:,Xy,1)/scalingFactor_y,'MARKER',m{1})
        ax1 = gca;
        xlim(ax1,limits_x)
        
        if kkk == 1
            maxValue = 25*ceil(max(max(X(:,Xy,:)/scalingFactor_y))/25);
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
        
        nexttile;
        dscatter(X(:,Xx,2)/scalingFactor_x,X(:,Xy,2)/scalingFactor_y,'MARKER',m{1})
        ax2 = gca;
        xlim(ax2,limits_x)
        ylim(ax2,limits_y)
        grid(ax2,'on')
        set(ax2,'Fontsize',fontSize)
        set(ax2,'XTick', round(limits_x(1):tick_x:limits_x(2)));
        set(ax2,'YTick', (limits_y(1):tick_y:limits_y(2)));
        set(gca, 'YTickLabel', [])
        set(ax2,'TickLabelInterpreter','latex')
        colormap(gca,'jet')
        
        nexttile;
        dscatter(X(:,Xx,3)/scalingFactor_x,X(:,Xy,3)/scalingFactor_y,'MARKER',m{1})
        ax3 = gca;
        xlim(ax3,limits_x)
        ylim(ax3,limits_y)
        grid(ax3,'on')
        set(ax3,'Fontsize',fontSize)
        set(ax3,'XTick', round(limits_x(1):tick_x:limits_x(2)));
        set(ax3,'YTick', (limits_y(1):tick_y:limits_y(2)));
        set(gca, 'YTickLabel', [])
        colormap(gca,'jet')
        
        set(ax3,'TickLabelInterpreter','latex')
        
        cbh = colorbar;
        cbh.Ticks = caxis( gca ) ;
        cbh.TickLabels = [{'Low'} {'High'}];
        cbh.TickLabelInterpreter = 'latex';
        ylabel(cbh,'Distribution Density','Interpreter','latex','Fontsize',fontSize)
        
        colobarLabelPos = get(get(cbh,'YLabel'),'Position');
        colobarLabelPos(1) = colobarLabelPos(1) - 3;
        set(get(cbh,'YLabel'),'Position',colobarLabelPos)
        
        linkaxes([ax1 ax2 ax3],'xy');
        
        t.TileSpacing = 'compact';
        xlabel(t,string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex','Fontsize',fontSize) % x-axis label
        ylabel(t,string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex','Fontsize',fontSize) % x-axis label
        
        hold off
        
        seedsAndInitializers = strcat(strcat(num2str(compilation(1).rawData.nodes),{'_'},num2str(compilation(1).rawData.seedInitializer)),{'_'},...
            strcat(num2str(compilation(2).rawData.nodes),{'_'},num2str(compilation(2).rawData.seedInitializer)),{'_'},...
            strcat(num2str(compilation(3).rawData.nodes),{'_'},num2str(compilation(3).rawData.seedInitializer)));
        
        
        
        
        
        
        
        figureSaveName = string(strcat(...
            compilation(1).figurePath,...
            'objectiveComparisons_MC_subplot_',...
            strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
            {'_'},seedsAndInitializers,{'_'},compilation(1).rawData.trajectoryType,zoom{kkk},'.png'));
        
        saveas(figure(fig_num),figureSaveName,'png');
        close(fig_num);
    end
end

%%

combinations = ...
    [1 9 8;...
    1 10 8;...
    1 11 8;...
    1 12 8;...
    1 21 8;...
    1 22 8;...
    1 23 8;...
    1 24 8 ];

if compilation(1).evolutions(1).population(1).size.collective < 10000
    
    
    for kk = 1:size(combinations,1)
        
        fig_num = p*100 + 3566000 + kk;
        figure(fig_num)
        %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
        set (gca,'Fontsize',15)
        
        Xx = combinations(kk,1);
        Xy = combinations(kk,2);
        Xz = combinations(kk,3);
        field_x = extremesAndConstraintsFieldNames{Xx};
        variableLabel_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).variableLabel;
        figureSaveName_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).figureSaveNameContent;
        units_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).units;
        limits_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).limits;
        scalingFactor_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).scalingFactor;
        
        field_y = extremesAndConstraintsFieldNames{Xy};
        variableLabel_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).variableLabel;
        figureSaveName_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).figureSaveNameContent;
        units_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).units;
        limits_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).limits;
        scalingFactor_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).scalingFactor;
        
        field_z = extremesAndConstraintsFieldNames{Xz};
        variableLabel_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).variableLabel;
        figureSaveName_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).figureSaveNameContent;
        units_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).units;
        limits_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).limits;
        scalingFactor_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).scalingFactor;
        
        
        ylim(limits_y)
        xlim(limits_x)
        zlim(limits_z)
        
        xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
        ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
        zlabel(string(strcat(variableLabel_z,{' '},units_z)),'Interpreter','latex') % z-axis label
        set(gca,'TickLabelInterpreter','latex')
        
        
        %title({string(strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' ')));string(strcat(variableLabel_y,{' vs '},variableLabel_x,{' vs '},variableLabel_z))})
        
        hold on
        grid on
        
        m = {'+', 'diamond', 'o', '*','square', '<','.', 'x', 'v', '^', '>', 'pentagram'};
        
        i = 1;
        for k = [1 size(X,3)]
            scatter3(X(:,Xx,k)/scalingFactor_x,X(:,Xy,k)/scalingFactor_y,X(:,Xz,k)/scalingFactor_z,'Marker',m{i})
            i = i + 1;
        end
        
        i = 1;
        for k = [1 size(X,3)]
            
            scatter3(X(:,Xx,k)/scalingFactor_x,limits_y(2)*fillerZ,X(:,Xz,k)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
            scatter3(limits_x(1)*fillerZ,X(:,Xy,k)/scalingFactor_y,X(:,Xz,k)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
            scatter(X(:,Xx,k)/scalingFactor_x,X(:,Xy,k)/scalingFactor_y,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
            
            i = i + 1;
        end
        
        plot3([0 0],[limits_y(2) limits_y(2)],[limits_z(1) limits_z(2)],'Color','k', 'LineStyle','-')
        plot3([0 0],[limits_y(1) limits_y(2)],[0 0],'Color','k', 'LineStyle','-')
        plot3([limits_x(1) limits_x(2)],[limits_y(2) limits_y(2)],[0 0],'Color','k', 'LineStyle','-')
        
        
        hold off
        
        axP = get(gca,'Position');
        set(gca, 'Position', axP)
        % legendtext = [
        %     strcat({'Generation: 0'});
        %     strcat({'Generation: '},num2str(numel(compilation(p).evolutions)))];
        %legendText = [{'Ascent'};{'Violation during ascent'};{'Descent'};{'Violation during descent'}];
        %[~,object_h,~,~] = legendflex(g,legendText,'Interpreter','latex');
        
        %legend(legendtext,'Location','northwest')
        %view([13 49])
        view([50 17])
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'objectiveComparisons_3D_',...
            strcat(figureSaveName_x,'_vs_',figureSaveName_z,'_vs_',figureSaveName_y),...
            '_Case',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
    
    
end




end