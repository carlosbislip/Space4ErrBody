function [  ] = plotTopIndividualsForDiscussion( compilation )
disp('   Top Individuals')

%%

totalGenerationsPerRun = nan(numel(compilation),2);
extremesAndConstraintsFieldNames  = fieldnames(compilation(end).evolutions(end).population(end).extremesAndConstraints);
extremesVectorLength = length(extremesAndConstraintsFieldNames);
for p = 1:numel(compilation)
    
    totalGenerationsPerRun(p,:) = [compilation(p).rawData.generationList(end) compilation(p).rawData.optimizingAlgorithmIndex ];
    
    X = nan(compilation(p).rawData.populationSize,extremesVectorLength,numel(compilation(p).rawData.generationList));
    topIDVector = nan(numel(compilation(p).rawData.generationList),1);
    topIndividual = nan(numel(compilation(p).rawData.generationList),extremesVectorLength);
    
   % championVector = nan(1,length(extremesAndConstraintsFieldNames),numel(compilation(p).evolutions));
    %fillerZ = ones(compilation(p).evolutions(1).population(1).size.collective,1);
    maxX = nan( numel(compilation(p).evolutions),length(extremesAndConstraintsFieldNames));
    minX = nan( numel(compilation(p).evolutions),length(extremesAndConstraintsFieldNames));
    
    
     
    for k = 1:numel(compilation(p).rawData.generationList)
        for ii = 1:compilation(p).rawData.populationSize
            X(ii,:,k) = cell2mat(compilation(p).rawData.extremesAndConstraintsDataPerGeneration(ii,7:end,k));
        end
        
        rankingOfIndividuals = [compilation(p).rawData.extremesAndConstraintsDataPerGeneration{:,3,k}];
        topIDVector(k) = find(rankingOfIndividuals == 1 );
        topIndividual(k,:) = cell2mat(compilation(p).rawData.extremesAndConstraintsDataPerGeneration(topIDVector(k),7:end,k));
        
    end
    
    
%     
%     
%     for k = 1:numel(compilation(p).rawData.generationList)
%         for ii = 1:compilation(p).rawData.populationSize
%             if compilation(p).evolutions(k).population(ii).rankVector == 1
%                 for j = 1:length(extremesAndConstraintsFieldNames)
%                     championVector(1,j,k) = compilation(p).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
%                 end
%             end
%             
%             for j = 1:length(extremesAndConstraintsFieldNames)
%                 X(ii,j,k) = compilation(p).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
%             end
%         end
%     end
%     
    for k = 1:numel(compilation(p).evolutions)
        maxX(k,:) = max(X(:,:,k));
        minX(k,:) = min(X(:,:,k));
    end
    
    objectiveComparisonsPerRun(p).topIndividualExtremesVector = topIndividual;
    objectiveComparisonsPerRun(p).allIndividualsExtremesVector = X;
    
end

%%


y01 = [ 1 8 9 10 11 12 21 22 23 24 ];
ylabelOffset = [ 7.25 1.75 9.25 9.25 9.25 9.25 9.25 9.25 9.25 9.25 ];

x0 = 0*ones(numel(y01),1);
x1 = ones(numel(y01),1);

combinations0 = [];%[ x0 y01' ylabelOffset'];
combinations1 = [ x1 y01' ylabelOffset' ];

combinations = [ combinations0;combinations1 ];

m = {'+','o','square','+','o','square','+','o','square','+','o','square','+','o','square','+','o','square'};
I12 = find(totalGenerationsPerRun(:,2)<3);
I3 = find(totalGenerationsPerRun(:,2)==3);
max_evolutions12 = 25*ceil(totalGenerationsPerRun(I12,1)/25);

limits_x0 = zeros(numel(compilation),1);

n = max(floor(log(abs(totalGenerationsPerRun(I3,1)))./log(10)));
max_evolutions3 = 2.5*10.^(n-1).*ceil(totalGenerationsPerRun(I3,1)./(2.5.*10.^(n-1)));

max_evolutions = [[max_evolutions12;max_evolutions3] totalGenerationsPerRun(:,2)];

limits_xAug = [limits_x0 max_evolutions];
limits_x12 = max(limits_xAug(I12,1:2));
nn = 2;
n = n - nn;
IHS_generationCountScalingfactor = 10^nn;

limits_x3 = max(limits_xAug(I3,1:2),[],1)/IHS_generationCountScalingfactor;

optimizingAlgorithms = [{'NSGA-II'} {'MOEAD/DE'} {'IHS'}];

for kk = 1:size(combinations,1)
    
    fig_num = 100 + 3566000 + kk;
    figure(fig_num)
    
    Xx = combinations(kk,1);
    
    if Xx == 0
        
        field_x = 'generation';
        variableLabel_x = 'Generation';
        figureSaveName_x = field_x;
        units_x = '$(-)$';
        scalingFactor_x = 1;
        
    else
        
        field_x = extremesAndConstraintsFieldNames{Xx};
        field_x_struct = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x);
        variableLabel_x = field_x_struct.variableLabel;
        figureSaveName_x = field_x_struct.figureSaveNameContent;
        units_x = field_x_struct.units;
        limits_x12 = field_x_struct.limits;
        limits_x3 = limits_x12;
        scalingFactor_x = field_x_struct.scalingFactor;
        tick_x = field_x_struct.tick;
    end
    
    Xy = combinations(kk,2);
    field_y = extremesAndConstraintsFieldNames{Xy};
    field_y_struct = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y);
    variableLabel_y = field_y_struct.variableLabel;
    figureSaveName_y = field_y_struct.figureSaveNameContent;
    units_y = field_y_struct.units;
    limits_y = field_y_struct.limits;
    scalingFactor_y = field_y_struct.scalingFactor;
    tick_y = field_y_struct.tick;
    
    if strcmp(field_y,'maximumBodyFrameMechanicalLoad')
        limits_y = [0 20];
        tick_y = 4;%(limits_y(2) - limits_y(1))/10;
    end
    hAX=axes;                 % first axes, save handle
    
    hold on
    grid(hAX(1),'on')
    nodes = cell(numel(compilation),1);
    method = cell(numel(compilation),1);
    seed = cell(numel(compilation),1);
    marker = cell(numel(compilation),1);
    colorOrder = get(gca, 'ColorOrder');
    mm = size(colorOrder,1);
    j = 1;
    
    for p = 1:numel(compilation)
        
        ColRow = mod(j,mm);
        if ColRow == 0
            ColRow = mm;
        end
        j = j + 1;
        
        y = objectiveComparisonsPerRun(p).champion(:,Xy)/scalingFactor_y;
        if Xx == 0
            x = linspace(1,size(y,1),size(y,1));
        else
            x = objectiveComparisonsPerRun(p).champion(:,Xx)/scalingFactor_x;
        end
        
        ppp(p) = scatter(hAX(1),x,y(:,1),'Marker',m{p},'MarkerEdgeColor',colorOrder(ColRow,:));
        nodes(p) = {num2str(compilation(p).rawData.nodes)};
        method(p) = {compilation(p).rawData.optimizingAlgorithm};
        seed(p) = {num2str(compilation(p).rawData.seedInitializer)};
        marker(p) = {''};
    end
    pause(1)
    content = [ nodes seed method ];
    griddedData = [ {'Nodes'} {'Seed'} {'Method'} ; content ];
    
    nrow = size(griddedData,1);
    ncol = size(griddedData,2);
    ngrp = nrow*ncol;
    
    for i = 1:ngrp
        hln(i) = line(nan(1,1), nan(1, 1), 'marker','none', 'linestyle', 'none');
    end
    
    allobj = reshape(hln,nrow,ncol);
    mm = size(colorOrder,1);
    j = 1;
    for p = 1:(size(allobj,1)-1)
        
        ColRow = mod(j,mm);
        if ColRow == 0
            ColRow = mm;
        end
        j = j + 1;
        
        allobj(p+1,1) = line(nan(1,1), nan(1,1), 'marker', m{p}, 'linestyle', 'none','Color',colorOrder(ColRow,:));
    end
    pause(1)
    
    alllbl = griddedData;
    legendflex(allobj(:), alllbl(:), 'nrow', nrow, 'xscale', 0.5,'Interpreter','latex');
    
    pos=get(hAX,'position') ;
    pos1=pos(2);              % save the original bottom position
    if Xx == 0
        pos(2)=pos(2)+pos1; pos(4)=pos(4)-pos1;  % raise bottom/reduce height->same overall upper position
        set(hAX,'position',pos)   % and resize first axes
        pos(2)=pos1; pos(4)=0.01; % reset bottom to original and small height
        hAX(2)=axes('position',pos,'color','none');  % and create the second
    end
    pause(1)
    
    limits_x = limits_x12;
    tick_x = limits_x(2)/10;
    xlim(hAX(1),limits_x)
    set(hAX(1),'XTick', round(limits_x(1):tick_x:limits_x(2)));
    set(hAX(1),'TickLabelInterpreter','latex')
    
    set(hAX(1),'Fontsize',15)
    if Xx == 0
        xlabel(hAX(1),string(strcat(optimizingAlgorithms(1),{', '},optimizingAlgorithms(2),{' : '},variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
    else
        xlabel(hAX(1),string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
    end
    
    pause(1)
    
    if Xx == 0
        limits_x = limits_x3;
        tick_x = limits_x(2)/10;
        xlim(hAX(2),limits_x)
        set(hAX(2),'XTick', round(limits_x(1):tick_x:limits_x(2)));
        set(hAX(2),'TickLabelInterpreter','latex')
        xlabel(hAX(2),string(strcat(optimizingAlgorithms(3),{' : '},variableLabel_x,{' $\times$10$^{'},num2str(n),{'}$ $(-)$'})),'Interpreter','latex') % x-axis label
        set(hAX(2),'Fontsize',15)
    end
    
    ylim(hAX(:),limits_y)
    set(hAX(:),'YTick', limits_y(1):tick_y:limits_y(2));
    pause(1)
    
    ylh = ylabel(hAX(1),string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex');
    if Xx == 0
        originalYLabelPosition_y = ylh.Position(2);
        ylh.Position(2) = originalYLabelPosition_y - combinations(kk,3);
    end
    pause(1)
    
    commonFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).fitnessVector.common);
    ascentPhaseFitnessVectorFieldNames = [];
    descentPhaseFitnessVectorFieldNames = [];
    
    if isfield(compilation(p).evolutions(end).population(end).fitnessVector,'ascent') || isfield(compilation(p).evolutions(end).population(end).fitnessVector,'descent')
        if compilation(p).rawData.trajectoryType == 'A'
            ascentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).fitnessVector.ascent);
        elseif compilation(p).rawData.trajectoryType == 'D'
            descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).fitnessVector.descent);
        elseif compilation(p).rawData.trajectoryType == 'AD'
            ascentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).fitnessVector.ascent);
            descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).fitnessVector.descent);
        end
    end
    
    fitnessVectorLength = length(commonFitnessVectorFieldNames) + length(ascentPhaseFitnessVectorFieldNames) + length(descentPhaseFitnessVectorFieldNames);
    
    if fitnessVectorLength == 1
        objectiveTypeForName = 'singleObjective';
    else
        objectiveTypeForName = 'multiObjective';
    end
    
    if compilation(p).rawData.trajectoryType == 'A'
        trajectoryTypeForName = 'ascent';
    elseif compilation(p).rawData.trajectoryType == 'D'
        trajectoryTypeForName = 'descent';
    elseif compilation(p).rawData.trajectoryType == 'AD'
        trajectoryTypeForName = 'coupled';
    end
    
    saveas(...
        figure(fig_num),...
        string(strcat(...
        compilation(p).figurePath,...
        objectiveTypeForName,'_ObjectiveComparisons_InitialPopulation_',figureSaveName_y, {'_vs_'},figureSaveName_x,{'_'},...
        trajectoryTypeForName,'_node_seed_comparison','.png')),...
        'png');
    close(fig_num);
    
end
%%
k = 1;
extremesAndConstraintsFieldNames  = fieldnames(compilation(end).evolutions(end).population(end).extremesAndConstraints);

X = nan(numel(compilation(end).evolutions(end).population(end).size.collective,length(extremesAndConstraintsFieldNames)),numel(compilation));
topIndividuals_extremesAndConstraints = nan(size(compilation(end).rawData.topIndividuals.extremesAndConstraints,1),size(compilation(end).rawData.topIndividuals.extremesAndConstraints(:,7:end),2));
for pp = 1:numel(compilation)
    topIndividuals_extremesAndConstraints(:,:,pp) = cell2mat(compilation(pp).rawData.topIndividuals.extremesAndConstraints(:,7:end));
    for ii = 1:compilation(pp).rawData.populationSize
        for j = 1:length(extremesAndConstraintsFieldNames)
            X(ii,j,pp) = compilation(pp).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
        end
    end
end

%%

XX = nan(numel(compilation(end).rawData.populationSize),2,numel(compilation));
topIndividuals_extremesAndConstraints = XX;%nan(size(compilation(end).rawData.topIndividuals.extremesAndConstraints,1),size(compilation(end).rawData.topIndividuals.extremesAndConstraints(:,7:end),2));
for pp = 1:numel(compilation)
    topIndividuals_extremesAndConstraints(:,1,pp) = cell2mat(compilation(pp).rawData.topIndividuals.extremesAndConstraints(:,6 + combinations(1)));
    topIndividuals_extremesAndConstraints(:,2,pp) = cell2mat(compilation(pp).rawData.topIndividuals.extremesAndConstraints(:,6 + combinations(2)));
    %for ii = 1:compilation(pp).rawData.populationSize
        %for j = 1:length(extremesAndConstraintsFieldNames)
            XX(:,1,pp) = cell2mat(compilation(pp).rawData.populationData(:,6 + combinations(1)));
            XX(:,2,pp) = cell2mat(compilation(pp).rawData.populationData(:,6 + combinations(2)));
        %end
    %end
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

extremesAndConstraintsFieldNames  = fieldnames(compilation(end).evolutions(end).population(end).extremesAndConstraints);

for kk = 1:size(combinations,1)
    
    zoom = [{''};{'_zoom'}];
    
    for kkk = 1:numel(zoom)
        fig_num = 100 + 3566000 + kk;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,750,420])
        
        Xx = combinations(kk,1);
        Xy = combinations(kk,2);
        field_x = extremesAndConstraintsFieldNames{Xx};
        field_x_struct = compilation(end).evolutions(end).population(end).extremesAndConstraints.(field_x);
        variableLabel_x = field_x_struct.variableLabel;
        figureSaveName_x = field_x_struct.figureSaveNameContent;
        units_x = field_x_struct.units;
        limits_x = field_x_struct.limits;
        scalingFactor_x = field_x_struct.scalingFactor;
        tick_x = 2*field_x_struct.tick;
        
        field_y = extremesAndConstraintsFieldNames{Xy};
        field_y_struct = compilation(end).evolutions(end).population(end).extremesAndConstraints.(field_y);
        variableLabel_y = field_y_struct.variableLabel;
        figureSaveName_y = field_y_struct.figureSaveNameContent;
        units_y = field_y_struct.units;
        limits_y = field_y_struct.limits;
        scalingFactor_y = field_y_struct.scalingFactor;
        tick_y = field_y_struct.tick;
        
        XX = nan(compilation(end).rawData.populationSize,2,numel(compilation));
        topIndividuals_extremesAndConstraints = nan(compilation(end).rawData.topIndividuals.size,2,numel(compilation));%nan(size(compilation(end).rawData.topIndividuals.extremesAndConstraints,1),size(compilation(end).rawData.topIndividuals.extremesAndConstraints(:,7:end),2));
        for pp = 1:numel(compilation)
            topIndividuals_extremesAndConstraints(:,1,pp) = cell2mat(compilation(pp).rawData.topIndividuals.extremesAndConstraints(:,6 + Xx));
            topIndividuals_extremesAndConstraints(:,2,pp) = cell2mat(compilation(pp).rawData.topIndividuals.extremesAndConstraints(:,6 + Xy));
            %for ii = 1:compilation(pp).rawData.populationSize
            %for j = 1:length(extremesAndConstraintsFieldNames)
            XX(:,1,pp) = cell2mat(compilation(pp).rawData.populationDataPerGeneration(:,6 + Xx,1));
            XX(:,2,pp) = cell2mat(compilation(pp).rawData.populationDataPerGeneration(:,6 + Xy,1));
            %end
            %end
        end
        
      
    
    m = {'square','+', 'diamond', 'o', '*','square', '<','.', 'x', 'v', '^', '>', 'pentagram'};
    
    fontSize = 15;
    
    colorOrder = get(gca, 'ColorOrder');
    
    if kkk == 1
        maxValue = 25*ceil(max(max(X(:,Xy,:)/scalingFactor_y))/25);
        limits_y = [0 maxValue];
        tick_y = maxValue/10;%(limits_y(2) - limits_y(1))/10;
    end
    
    t = tiledlayout(1,3);
    ax1 = nexttile;
    if Xx == 0
        x1 = linspace(1,size(y,1),size(y,1));
        x2 = linspace(1,size(y,1),size(y,1));
    else
        x1 = X(:,Xx,1)/scalingFactor_x;
        x2 = topIndividuals_extremesAndConstraints(:,Xx,1)/scalingFactor_x;
    end
    y1 = X(:,Xy,1)/scalingFactor_y;
    y2 = topIndividuals_extremesAndConstraints(:,Xy,1)/scalingFactor_y;
    % ax1 = gca;
    hold(ax1,'on')
    dscatter(x1,y1,'MARKER',m{1},'AXES',ax1,'ALPHA',.1)
    scatter(ax1,x2,y2,50,colorOrder(2,:),'Marker','+','MarkerFaceAlpha',1);
    
    hold(ax1,'off')
    
    xlim(ax1,limits_x)
    ylim(ax1,limits_y)
    grid(ax1,'on')
    set(ax1,'Fontsize',fontSize)
    set(ax1,'XTick', round(limits_x(1):tick_x:limits_x(2)));
    set(ax1,'YTick', limits_y(1):tick_y:limits_y(2));
    set(ax1,'TickLabelInterpreter','latex')
    colormap(ax1,'jet')
    
    ax2 = nexttile;
    if Xx == 0
        x1 = linspace(1,size(y,1),size(y,1));
        x2 = linspace(1,size(y,1),size(y,1));
    else
        x1 = X(:,Xx,2)/scalingFactor_x;
        x2 = topIndividuals_extremesAndConstraints(:,Xx,2)/scalingFactor_x;
    end
    y1 = X(:,Xy,2)/scalingFactor_y;
    y2 = topIndividuals_extremesAndConstraints(:,Xy,2)/scalingFactor_y;
    %ax2 = gca;
    hold(ax2,'on')
    dscatter(x1,y1,'MARKER',m{1},'AXES',ax2,'ALPHA',.1)
    scatter(ax2,x2,y2,50,colorOrder(2,:),'Marker','+','MarkerFaceAlpha',.1);
    hold(ax2,'off')
    
    xlim(ax2,limits_x)
    ylim(ax2,limits_y)
    grid(ax2,'on')
    set(ax2,'Fontsize',fontSize)
    set(ax2,'XTick', round(limits_x(1):tick_x:limits_x(2)));
    set(ax2,'YTick', (limits_y(1):tick_y:limits_y(2)));
    set(ax2, 'YTickLabel', [])
    set(ax2,'TickLabelInterpreter','latex')
    colormap(ax2,'jet')
    
    ax3 =  nexttile;
    if Xx == 0
        x1 = linspace(1,size(y,1),size(y,1));
        x2 = linspace(1,size(y,1),size(y,1));
    else
        x1 = X(:,Xx,3)/scalingFactor_x;
        x2 = topIndividuals_extremesAndConstraints(:,Xx,3)/scalingFactor_x;
    end
    y1 = X(:,Xy,3)/scalingFactor_y;
    y2 = topIndividuals_extremesAndConstraints(:,Xy,3)/scalingFactor_y;
    
    hold(ax3,'on')
    dscatter(x1,y1,'MARKER',m{1},'AXES',ax3,'ALPHA',.1)
    dd = scatter(ax3,x2,y2,50,colorOrder(2,:),'Marker','+','MarkerFaceAlpha',.1);
    hold(ax3,'off')
    
    xlim(ax3,limits_x)
    ylim(ax3,limits_y)
    grid(ax3,'on')
    set(ax3,'Fontsize',fontSize)
    set(ax3,'XTick', round(limits_x(1):tick_x:limits_x(2)));
    set(ax3,'YTick', (limits_y(1):tick_y:limits_y(2)));
    set(ax3, 'YTickLabel', [])
    set(ax3,'TickLabelInterpreter','latex')
        legend(dd, {'Top Individuals'},'Interpreter','latex','Location','northoutside');

    
    colormap(ax3,'jet')
    
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
        'topIndividuals_over_MC_subplot_',...
        strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
        {'_'},seedsAndInitializers,{'_'},compilation(1).rawData.trajectoryType,zoom{kkk},'.png'));
    
    saveas(figure(fig_num),figureSaveName,'png');
    close(fig_num);
    end
end




end