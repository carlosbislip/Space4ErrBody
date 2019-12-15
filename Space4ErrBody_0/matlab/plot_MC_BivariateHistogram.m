function [  ] = plot_MC_BivariateHistogram( compilation )
disp('   Monte Carlo Bivariate Histogram')

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

%% Bivariate Histogram


for kk = 1:size(combinations,1)
    
    for kkk = 1%:numel(zoom)
        fig_num = 100 + 3566000 + kk;
        figure(fig_num)
        % set(figure(fig_num),'units','pixels','position',[0,0,750,420])
        
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
        set(ax1,'YTick', limits_y(1):2*tick_y:limits_y(2));
        set(ax1,'TickLabelInterpreter','latex')
        colormap(gca,'jet')
        %[X_edges,Y_edges] = meshgrid(linspace(limits_x(1),limits_x(2),1000),linspace(limits_y(1),limits_y(2),1000));% - limits_x(2))/100:limits_x(2),limits_y(1):tick_y:limits_y(2));
        nbins = 250;
        
        
        
        hist2d(X(:,Xx,1)/scalingFactor_x,X(:,Xy,1)/scalingFactor_y,linspace(min(limits_x),max(limits_x),nbins+1),linspace(min(limits_y),max(limits_y),nbins+1),'bar3');
        pause(1);
        ylim(ax1,limits_y)
        xlim(ax1,limits_x)
        
        limits_z = zlim(ax1);
        
        n = max(floor(log(limits_z(2))./log(10)));
        limits_z(2) = 2.5*10.^(n-1).*ceil(limits_z(2)./(2.5.*10.^(n-1)));
        zlim(ax1,limits_z)
        set(ax1,'ZTick', limits_z(1):limits_z(2)/5:limits_z(2));
        set(ax1,'TickLabelInterpreter','latex')
        
        
        cbh = colorbar;
        cbh.Ticks = caxis( gca ) ;
        cbh.TickLabels = [{'Low'} {'High'}];
        cbh.TickLabelInterpreter = 'latex';
        ylabel(cbh,'Distribution Density','Interpreter','latex','Fontsize',fontSize)
        pause(1);
        
        colobarLabelPos = get(get(cbh,'YLabel'),'Position');
        colobarLabelPos(1) = colobarLabelPos(1) - 1;
        set(get(cbh,'YLabel'),'Position',colobarLabelPos)
        
        xlabel(ax1,string(strcat('$\theta_{\mathrm{ToGo}}$',{' '},units_x)),'Interpreter','latex','Fontsize',fontSize) % x-axis label
        ylabel(ax1,string(strcat('$n_{\mathrm{Body},\mathrm{max}}$',{' '},units_y)),'Interpreter','latex','Fontsize',fontSize) % x-axis label
        zlabel(ax1,string(strcat('Count',{' '},'$(-)$')),'Interpreter','latex','Fontsize',fontSize) % x-axis label
        pause(1);
        
        align_axislabel([], gca)
        % axis equal
        %h = rotate3d(fig_num);
        %set(h,'ActionPostCallback',@align_axislabel)
        %set(h,'ActionPostCallback',@axislabel_rotation)
        %view(3)
        %axesLabelsAlign3D
        %hold off
        
        seedsAndInitializers = strcat(strcat(num2str(compilation(1).rawData.nodes),{'_'},num2str(compilation(1).rawData.seedInitializer)));%,{'_'},...
        figureSaveName = string(strcat(...
            compilation(1).figurePath,...
            'objectiveComparisons_MC_BivariateHistogram_',...
            strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
            {'_'},seedsAndInitializers,{'_'},compilation(1).rawData.trajectoryType,{'_'},compilation(1).rawData.objectiveFunctionCase,'.png'));
        
        saveas(figure(fig_num),figureSaveName,'png');
        close(fig_num);
        
    end
end





end