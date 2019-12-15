function [  ] = plotObjectiveComparisonsForDiscussion( compilation )
disp('   Objective Comparisons')

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


combinations = [ [x0;x1] [y01';y01'] [ylabelOffset';ylabelOffset'] ];

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

if limits_x12(2) > limits_x3(2)
    
    limits_x3(2) = limits_x12(2);
else
    limits_x12(2) = limits_x3(2);
end

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
        
        y = objectiveComparisonsPerRun(p).topIndividualExtremesVector(:,Xy)/scalingFactor_y;
        if Xx == 0
            x = linspace(1,size(y,1),size(y,1));
        else
            x = objectiveComparisonsPerRun(p).topIndividualExtremesVector(:,Xx)/scalingFactor_x;
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

        legendflex(allobj(:), alllbl(:), 'nrow', nrow,'fontsize',12,'xscale', 0.4,'Interpreter','latex');

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
            descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population.fitnessVector.descent);
        elseif compilation(p).rawData.trajectoryType == 'AD'
            ascentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population.fitnessVector.ascent);
            descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population.fitnessVector.descent);
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
        objectiveTypeForName,'_ObjectiveComparisons_topIndividualExtremesVectors_',figureSaveName_y, {'_vs_'},figureSaveName_x,{'_'},...
        trajectoryTypeForName,'_node_seed_comparison','.png')),...
        'png');
    close(fig_num);
    
end
%%

combinations = ...
    [1 8 1.75;...
    1 9 9.25;...
    1 10 9.25;...
    1 11 9.25;...
    1 12 9.25;...
    1 21 9.25;...
    1 22 9.25;...
    1 23 9.25;...
    1 24 9.25];

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
limits_x3 = max(limits_xAug(I3,1:2),[],1)/100;

optimizingAlgorithms = [{'NSGA-II'} {'MOEAD/DE'} {'IHS'}];

for kk = 1:size(combinations,1)
    
    fig_num = 100 + 3566000 + kk;
    figure(fig_num)
    
    Xx = combinations(kk,1);
    field_x = extremesAndConstraintsFieldNames{Xx};
    field_x_struct = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x);
    variableLabel_x = field_x_struct.variableLabel;
    figureSaveName_x = field_x_struct.figureSaveNameContent;
    units_x = field_x_struct.units;
    limits_x = field_x_struct.limits;
    scalingFactor_x = field_x_struct.scalingFactor;
    tick_x = field_x_struct.tick;
    
    Xy = combinations(kk,2);
    field_y = extremesAndConstraintsFieldNames{Xy};
    field_y_struct = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y);
    variableLabel_y = field_y_struct.variableLabel;
    figureSaveName_y = field_y_struct.figureSaveNameContent;
    units_y = field_y_struct.units;
    limits_y = field_y_struct.limits;
    scalingFactor_y = field_y_struct.scalingFactor;
    tick_y = field_y_struct.tick;
    
    
    
    ylim(limits_y)
    xlim(limits_x)
    xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
    ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
    set(gca,'TickLabelInterpreter','latex')
    
    set(gca,'YTick', round(limits_y(1):tick_y:limits_y(2)));
    
    %title({strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' '));string(strcat(variableLabel_y,{' vs '},variableLabel_x))})
    
    hold on
    grid on
    
    m = {'+', 'diamond', 'o', '*','square', '<','.', 'x', 'v', '^', '>', 'pentagram'};
    
    %XX = reshape(permute(X,[1 3 2]),[],size(X,2));
    q = 1;
    if numel(compilation(p).evolutions) == 1
        
        for i = 1%[1 size(X,3)]
            g(q) = scatter(X(:,Xx,i)/scalingFactor_x,X(:,Xy,i)/scalingFactor_y,'Marker',m{q});
            q = q + 1;
        end
        for i = [1 size(X,3)]
            % g(q) = scatter(topIndividualExtremesVectorVector(i,Xx,:)/scalingFactor_x,topIndividualExtremesVectorVector(i,Xy,:)/scalingFactor_y,'Marker',m{q},'MarkerFaceColor', 'flat');
            %  scatter(topIndividualExtremesVectorMagnitude(i,Xx,:)/scalingFactor_x,topIndividualExtremesVectorMagnitude(i,Xy,:)/scalingFactor_y,'Marker',m{q+1},'MarkerFaceColor', 'flat')
            q = q + 2;
        end
    else
        
        for i = [1 size(X,3)]
            % scatter(XX(i,Xx,:)/scalingFactor_x,XX(i,Xy,:)/scalingFactor_y,'Marker',m{1})
            g(q) = scatter(X(:,Xx,i)/scalingFactor_x,X(:,Xy,i)/scalingFactor_y,'Marker',m{q});
            q = q + 1;
        end
        for i = [1 size(X,3)]
            % g(q) = scatter(topIndividualExtremesVectorVector(i,Xx,:)/scalingFactor_x,topIndividualExtremesVectorVector(i,Xy,:)/scalingFactor_y,'Marker',m{q},'MarkerFaceColor', 'flat');
            %  scatter(topIndividualExtremesVectorMagnitude(i,Xx,:)/scalingFactor_x,topIndividualExtremesVectorMagnitude(i,Xy,:)/scalingFactor_y,'Marker',m{q+1},'MarkerFaceColor', 'flat')
            q = q + 2;
        end
    end
    
    hold off
    
    axP = get(gca,'Position');
    set(gca, 'Position', axP)
    %legendText = [
    %   strcat({'Generation: 0'});
    %   strcat({'Generation: '},num2str(numel(compilation(p).evolutions)));
    %   {'Vector topIndividualExtremesVector - Generation: 0 '};
    %   {'Magnitude topIndividualExtremesVector - Generation: 0 '};
    %   strcat({'Vector topIndividualExtremesVector - '},{'Generation: '},num2str(numel(compilation(p).evolutions)));
    %   strcat({'Magnitude topIndividualExtremesVector - '},{'Generation: '},num2str(numel(compilation(p).evolutions)))];
    %legend(legendtext,'Location','northwest')
    %[~,object_h,~,~] = legendflex(g,legendText,'Interpreter','latex');
    
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).figurePath,...
        'objectiveComparisons_',...
        strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
        '_Case',...
        convertCharsToStrings(compilation(p).case),...
        '.png'),...
        'png');
    close(fig_num);
    
end

%%
for p = 1:numel(compilation)
    
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