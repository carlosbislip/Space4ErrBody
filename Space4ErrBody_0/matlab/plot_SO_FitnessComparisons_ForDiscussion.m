function [  ] = plotFitnessComparisonsForDiscussion( compilation )
disp('   Fitness Comparisons')

%%

totalGenerationsPerRun = nan(numel(compilation),2);
commonFitnessVectorFieldNames  = fieldnames(compilation(end).evolutions(end).population(end).fitnessVector.common);
ascentPhaseFitnessVectorFieldNames = [];
descentPhaseFitnessVectorFieldNames = [];

if isfield(compilation(end).evolutions(end).population(end).fitnessVector,'ascent') || isfield(compilation(end).evolutions(end).population(end).fitnessVector,'descent')
    if compilation(end).rawData.trajectoryType == 'A'
        ascentPhaseFitnessVectorFieldNames  = fieldnames(compilation(end).evolutions(end).population(end).fitnessVector.ascent);
    elseif compilation(end).rawData.trajectoryType == 'D'
        descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(end).evolutions(end).population.fitnessVector.descent);
    elseif compilation(end).rawData.trajectoryType == 'AD'
        ascentPhaseFitnessVectorFieldNames  = fieldnames(compilation(end).evolutions(end).population.fitnessVector.ascent);
        descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(end).evolutions(end).population.fitnessVector.descent);
    end
end

fitnessVectorLength = length(commonFitnessVectorFieldNames)...
    +length(ascentPhaseFitnessVectorFieldNames)...
    +length(descentPhaseFitnessVectorFieldNames);

fitnessVectorNames = [commonFitnessVectorFieldNames ascentPhaseFitnessVectorFieldNames descentPhaseFitnessVectorFieldNames] ;

for p = 1:numel(compilation)
    totalGenerationsPerRun(p,:) = [compilation(p).rawData.generationList(end) compilation(p).rawData.optimizingAlgorithmIndex ];
    topIDVector = nan(numel(compilation(p).rawData.generationList),1);
    topIndividual = nan(numel(compilation(p).rawData.generationList),fitnessVectorLength);
    
    X = nan(compilation(p).rawData.populationSize,fitnessVectorLength,numel(compilation(p).rawData.generationList));
    
    for k = 1:numel(compilation(p).rawData.generationList)
        for ii = 1:compilation(p).rawData.populationSize
            X(ii,:,k) = cell2mat(compilation(p).rawData.fitnessDataPerGeneration(ii,7:7+(fitnessVectorLength-1),k));
        end
        
        rankingOfIndividuals = [compilation(p).rawData.populationDataPerGeneration{:,3,k}];
        topIDVector(k) = find(rankingOfIndividuals == 1 );
        topIndividual(k,:) = cell2mat(compilation(p).rawData.fitnessDataPerGeneration(topIDVector(k),7:7+(fitnessVectorLength-1),k));
        
    end
    
    fitnessComparisons(p).topIndividualFitnessVector = topIndividual;
    fitnessComparisons(p).allIndividualsFintessVector = X;
    
end


%% Fitness Bins - Top Individuals

m = {'+','o','square','+','o','square','+','o','square','+','o','square','+','o','square','+','o','square'};
I12 = find(totalGenerationsPerRun(:,2)<3);
I3 = find(totalGenerationsPerRun(:,2)==3);
max_evolutions12 = 25*ceil(totalGenerationsPerRun(I12,1)/25);

limits_x0 = zeros(p,1);

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

combinations = ...
    [1 1 9.25;...
    1 2 1.75;...
    1 3 9.25;...
    1 4 9.25;...
    1 5 9.25;...
    1 6 9.25;...
    1 7 9.25;...
    1 8 9.25;...
    1 9 9.25;...
    1 10 9.25;...
    1 11 9.25];

for kk = 1:size(combinations,1)
    
    fig_num = 100 + 3466000 + kk;
    figure(fig_num)
    
    Xx = combinations(1,1);
    field_x = 'generation';
    variableLabel_x = 'Generation';
    figureSaveName_x = field_x;
    units_x = '$(-)$';
    scalingFactor_x = 1;
    
    Xy = combinations(kk,2);
    field_y = fitnessVectorNames{Xy};
    if ismember(field_y,convertStringsToChars(string(strcat({'central_target_angular_distance_to_go'},{'heat_load'},{'fuel_mass'}))))
        field_y_struct = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y);
    else
        if isfield(compilation(p).evolutions(end).population(end).fitnessVector,'ascent')
            field_y_struct = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y);
        elseif isfield(compilation(p).evolutions(end).population(end).fitnessVector,'descent')
            field_y_struct = compilation(p).evolutions(end).population(end).fitnessVector.descent.(field_y);
        end
    end
    
    variableLabel_y = field_y_struct.variableLabel;
    figureSaveName_y = field_y_struct.figureSaveNameContent;
    units_y = field_y_struct.units;
    limits_y = field_y_struct.limits;
    scalingFactor_y = field_y_struct.scalingFactor;
    tick_y = field_y_struct.tick;
    
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
        
        y = fitnessComparisons(p).topIndividualsVector(:,Xy);
        x = linspace(1,size(y,1),size(y,1));
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
    
    
    pos=get(hAX,'position') ;
    pos1=pos(2);              % save the original bottom position
    pos(2)=pos(2)+pos1; pos(4)=pos(4)-pos1;  % raise bottom/reduce height->same overall upper position
    set(hAX,'position',pos)   % and resize first axes
    pos(2)=pos1; pos(4)=0.01; % reset bottom to original and small height
    hAX(2)=axes('position',pos,'color','none');  % and create the second
    pause(1)
    
    limits_x = limits_x12;
    tick_x = limits_x(2)/10;
    xlim(hAX(1),limits_x)
    set(hAX(1),'XTick', round(limits_x(1):tick_x:limits_x(2)));
    set(hAX(1),'TickLabelInterpreter','latex')
    
    set(hAX(1),'Fontsize',15)
    xlabel(hAX(1),string(strcat(optimizingAlgorithms(1),{', '},optimizingAlgorithms(2),{' : '},variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
    
    alllbl = griddedData;
    legendflex(allobj(:), alllbl(:), 'nrow', nrow,'fontsize',12,'xscale', 0.4,'Interpreter','latex');
    
    pause(1)
    
    limits_x = limits_x3;
    tick_x = limits_x(2)/10;
    xlim(hAX(2),limits_x)
    set(hAX(2),'XTick', round(limits_x(1):tick_x:limits_x(2)));
    set(hAX(2),'TickLabelInterpreter','latex')
    xlabel(hAX(2),string(strcat(optimizingAlgorithms(3),{' : '},variableLabel_x,{' $\times$10$^{'},num2str(n),{'}$ $(-)$'})),'Interpreter','latex') % x-axis label
    set(hAX(2),'Fontsize',15)
    
    ylim(hAX(:),limits_y)
    set(hAX(:),'YTick', limits_y(1):tick_y:limits_y(2));
    pause(1)
    
    ylh = ylabel(hAX(1),string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex');
    originalYLabelPosition_y = ylh.Position(2);
    ylh.Position(2) = originalYLabelPosition_y - combinations(kk,3);
    
    pause(1)
    
    hold off
    
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
        objectiveTypeForName,'_FitnessComparisons_champions_',figureSaveName_y, {'_vs_'},figureSaveName_x,{'_'},...
        trajectoryTypeForName,'_node_seed_comparison','.png')),...
        'png');
    close(fig_num);
    
end

%%
if isempty(ascentPhaseFitnessVectorFieldNames) == true
    if size(topIndividual,2) > 1
        combinations = nchoosek(1:length(commonFitnessVectorFieldNames),2);
        
        for kk = 1:numel(combinations)
            
            Xx = combinations(kk,1);
            field_x = commonFitnessVectorFieldNames{Xx};
            variableLabel_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).variableLabel;
            figureSaveName_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).figureSaveNameContent;
            units_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).units;
            limits_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).limits;
            scalingFactor_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).scalingFactor;
            
            Xy = combinations(kk,2);
            field_y = commonFitnessVectorFieldNames{Xy};
            variableLabel_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).variableLabel;
            figureSaveName_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).figureSaveNameContent;
            units_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).units;
            limits_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).limits;
            scalingFactor_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).scalingFactor;
            
            Xz = combinations(kk,3);
            field_z = commonFitnessVectorFieldNames{Xz};
            variableLabel_z = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_z).variableLabel;
            figureSaveName_z = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_z).figureSaveNameContent;
            units_z = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_z).units;
            limits_z = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_z).limits;
            scalingFactor_z = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_z).scalingFactor;
            
            
            fig_num = p*100 + 3466000 + kk;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            %set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
            set (gca,'Fontsize',15)
            %title(string(strcat({'Objective Comparisons '},num2str(kk),{' - '},strrep(convertCharsToStrings(compilation(p).case),'_',' '))))
            
            xlim(limits_x)
            ylim(limits_y)
            zlim(limits_z)
            
            xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
            ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
            zlabel(string(strcat(variableLabel_z,{' '},units_z)),'Interpreter','latex') % z-axis label
            set(gca,'TickLabelInterpreter','latex')
            
            hold on
            grid on
            
            i = 1;
            for k = [1 numel(compilation(p).evolutions)]
                
                scatter3(X(:,1,k),1e5*fillerZ,X(:,combinations(kk,2)+1,k),'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                scatter3(0*fillerZ,X(:,combinations(kk,1)+1,k),X(:,combinations(kk,2)+1,k),'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                scatter(X(:,1,k),X(:,combinations(kk,1)+1,k),'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                
                i = i + 1;
            end
            
            plot3([0 0],[1e5 1e5],[1e-5 1e5],'Color','k', 'LineStyle','-')
            plot3([0 0],[1e-5 1e5],[1e-5 1e-5],'Color','k', 'LineStyle','-')
            plot3([0 1000],[1e5 1e5],[1e-5 1e-5],'Color','k', 'LineStyle','-')
            
            hold off
            
            axP = get(gca,'Position');
            set(gca, 'Position', axP)
            legendtext = [{'Initial Population'};{strcat('Generation:_{ }', num2str(numel(compilation(p).evolutions)))}];
            legend(legendtext,'Location','northwest','Interpreter','latex')
            %view([13 49])
            view([50 17])
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'fitnessBinComparisons_',...
                num2str(kk),...
                '_Case',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    else
        if numel(compilation(p).evolutions) > 1
            combinations = [1 1];
            max_evolutions = 1;
            limits_x = [0 max_evolutions];
            Xx = combinations(1,1);
            field_x = 'generation';
            variableLabel_x = 'Generation';
            figureSaveName_x = 'generation';
            units_x = '$(-)$';
            
            max_evolutions = 25*ceil(max((numel(compilation(p).evolutions)-1))/25);
            limits_x = [0 max_evolutions];
            scalingFactor_x = 1;
            
            Xy = combinations(1,2);
            field_y = commonFitnessVectorFieldNames{Xy};
            variableLabel_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).variableLabel;
            figureSaveName_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).figureSaveNameContent;
            units_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).units;
            limits_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).limits;
            scalingFactor_y = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_y).scalingFactor;
            
            fig_num = p*100 + 3466000 + 1;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
            set (gca,'Fontsize',15)
            %title({strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' '));string(strcat(variableLabel_y,{' vs '},{'Generation'}))})
            max_tof = max([compilation(p).evolutions.max_tof]);
            ylim(limits_y)
            
            if numel(compilation(p).evolutions) > 1
                xlim(limits_x)
                set(gca,'XTick', round(0:max_evolutions/10:max_evolutions));
            end
            
            xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
            ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
            set(gca,'TickLabelInterpreter','latex')
            
            hold on
            grid on
            
            %plot(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),topIndividualsVector);
            ax = gca;
            ax.ColorOrderIndex = 1;
            ppp(1) = scatter(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),topIndividual,'Marker',m{1},'DisplayName',string(strcat({'Champion: Vector'})));
            
            % plot(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),topIndividualsMagnitude);
            %ax = gca;
            %ax.ColorOrderIndex = 2;
            %ppp(2) = scatter(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),topIndividualsMagnitude,'Marker',m{2},'DisplayName',string(strcat({'Champion: Magnitude'})));
            
            hold off
            
            axP = get(gca,'Position');
            set(gca, 'Position', axP)
            %legendtext = [{'Initial Population'};{strcat('Generation:_{ }', num2str(numel(compilation(p).evolutions)))}];
            
            
            %legend(ppp(:),'Location','northwest','Interpreter','latex')
            %view([13 49])
            
            saveas(...
                figure(fig_num),...
                string(strcat(...
                compilation(p).figurePath,...
                'fitnessBinComparisons_topIndividuals_',field_y, {'_vs_'},{'generation'},...
                '_Case',...
                convertCharsToStrings(compilation(p).case),...
                '.png')),...
                'png');
            close(fig_num);
            
        end
    end
    
else
    
    if size(ascentPhaseFitnessVectorFieldNames,1) > 1
        
        combinations = nchoosek(1:length(ascentPhaseFitnessVectorFieldNames),2);
        combinations = combinations(find(combinations(:,1) == 1),:);
        XX = X;
        
        X = topIndividual;
        
        for kk = 1:size(combinations,1)
            
            Xx = 1;
            field_x = commonFitnessVectorFieldNames{Xx};
            variableLabel_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).variableLabel;
            figureSaveName_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).figureSaveNameContent;
            units_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).units;
            limits_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).limits;
            scalingFactor_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).scalingFactor;
            
            Xy = combinations(kk,1);
            field_y = ascentPhaseFitnessVectorFieldNames{Xy};
            variableLabel_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).variableLabel;
            figureSaveName_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).figureSaveNameContent;
            units_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).units;
            limits_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).limits;
            scalingFactor_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).scalingFactor;
            
            Xz = combinations(kk,2);
            field_z = ascentPhaseFitnessVectorFieldNames{Xz};
            variableLabel_z = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_z).variableLabel;
            figureSaveName_z = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_z).figureSaveNameContent;
            units_z = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_z).units;
            limits_z = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_z).limits;
            scalingFactor_z = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_z).scalingFactor;
            
            fig_num = p*100 + 3466000 + kk;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            %set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
            set (gca,'Fontsize',15)
            view([50 17])
            
            set(gca,'YScale','log')
            set(gca,'ZScale','log')
            
            xlim(limits_x)
            ylim([1e-5 1e5])
            zlim([1e-5 1e5])
            
            set(gca,'ztick',[1e-5 1e-3 1e-1 1e1 1e3 1e5]);
            set(gca,'ytick',[1e-5 1e-3 1e-1 1e1 1e3 1e5]);
            
            hold on
            grid on
            ax = gca;
            set(ax,'xminorgrid','on','zminorgrid','on','yminorgrid','on')
            
            generationsShown = round(linspace(1,size(XX,3),3));
            XX(find(XX(:,:,:) == 0)) = 1e-5;
            
            fillerZ = ones(size(XX,1),1);
            
            i = 1;
            for k = generationsShown
                scatter3(XX(:,Xx,k)/scalingFactor_x,1e5*fillerZ,XX(:,Xz+1,k)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                scatter3(1e-5*fillerZ,XX(:,Xy+1,k)/scalingFactor_y,XX(:,Xz+1,k)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                scatter(XX(:,Xx,k)/scalingFactor_x,XX(:,Xy+1,k)/scalingFactor_y,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                i = i + 1;
            end
            
            plot3(1e-5*ones(3,1),1e5*ones(3,1),linspace(1e-5,1e5,3),'Color','k', 'LineStyle','-')
            plot3(1e-5*ones(3,1),linspace(1e-5,1e5,3),1e-5*ones(3,1),'Color','k', 'LineStyle','-')
            plot3(linspace(limits_x(1),limits_x(2),3),1e5*ones(3,1),1e-5*ones(3,1),'Color','k', 'LineStyle','-')
            i = 1;
            ppp = [];
            
            ax.ColorOrderIndex = 1;
            
            for k = generationsShown
                ppp(i) = scatter3(XX(:,Xx,k)/scalingFactor_x,XX(:,Xy+1,k)/scalingFactor_y,XX(:,Xz+1,k)/scalingFactor_z,'Marker',m{i},'DisplayName',string(strcat({'Generation: '}, num2str(generationsShown(i) - 1))));
                i = i + 1;
            end
            
            hold off
            
            axP = get(gca,'Position');
            set(gca, 'Position', axP)
            legend(ppp(:),'Location','northwest','Interpreter','latex')
            %title({strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' '));string(strcat(variableLabel_z,{' vs '},variableLabel_x,{' vs '},variableLabel_y))})
            xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
            ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
            zlabel(string(strcat(variableLabel_z,{' '},units_z)),'Interpreter','latex') % z-axis label
            set(gca,'TickLabelInterpreter','latex')
            
            view([50 17])
            
            %                     set(get(gca,'xlabel'),'rotation',-10);
            %                     set(get(h_a,'ylabel'),'rotation',0);
            %                     set(get(h_a,'zlabel'),'rotation',0);
            %
            %
            %                     stored_DataAspectRatio = get(gca,'DataAspectRatio');
            %                     stored_OuterPosition = get(gca,'OuterPosition');
            %                     stored_PlotBoxAspectRatio = get(gca,'PlotBoxAspectRatio');
            %
            %                     align_axislabel_log([], gca);
            %
            %                     set(gca, 'DataAspectRatio', stored_DataAspectRatio)
            %                     set(gca, 'OuterPosition', stored_OuterPosition)
            %                     set(gca, 'PlotBoxAspectRatio', stored_PlotBoxAspectRatio)
            %
            
            
            
            %                     h_a = gca;
            %                     [thx,thy,thz] = axislabel_rotation_angle(h_a);
            %                     set(get(h_a,'xlabel'),'rotation',thx(1));
            %                     set(get(h_a,'ylabel'),'rotation',thy(1));
            %                     set(get(h_a,'zlabel'),'rotation',thz(1));
            
            saveas(...
                figure(fig_num),...
                string(strcat(...
                compilation(p).figurePath,...
                'fitnessBinComparisons_',field_x, {'_vs_'},field_y, {'_vs_'},field_z,...
                '_Case',...
                convertCharsToStrings(compilation(p).case),...
                '.png')),...
                'png');
            close(fig_num);
            
            fig_num = p*100 + 3467000 + kk;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
            set (gca,'Fontsize',15)
            
            xlim(limits_x)
            ylim(limits_y)
            zlim(limits_z)
            
            xlabel(string(strcat(variableLabel_x,{' '},units_x))) % x-axis label
            ylabel(string(strcat(variableLabel_y,{' '},units_y))) % y-axis label
            zlabel(string(strcat(variableLabel_z,{' '},units_z))) % z-axis label
            
            hold on
            grid on
            
            
            fillerZ = ones(size(X,1),1);
            
            scatter3(X(:,Xx)/scalingFactor_x,limits_y(2)*fillerZ,X(:,Xz+1)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{1})
            scatter3(limits_x(1)*fillerZ,X(:,Xy+1)/scalingFactor_y,X(:,Xz+1)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{1})
            scatter(X(:,Xx)/scalingFactor_x,X(:,Xy+1)/scalingFactor_y,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{1})
            
            plot3([0 0],[limits_y(2) limits_y(2)],[limits_z(1) limits_z(2)],'Color','k', 'LineStyle','-')
            plot3([0 0],[limits_y(1) limits_y(2)],[0 0],'Color','k', 'LineStyle','-')
            plot3([limits_x(1) limits_x(2)],[limits_y(2) limits_y(2)],[0 0],'Color','k', 'LineStyle','-')
            
            color = linspace(1,size(topIndividual,1),size(topIndividual,1));
            scale = 25*ones(size(topIndividual,1),1);
            
            scatter3(X(:,Xx)/scalingFactor_x,X(:,Xy+1)/scalingFactor_y,X(:,Xz+1)/scalingFactor_z,scale, color,'Marker',m{1})
            
            caxis([1 size(X,1)])
            cTicks = round(1:size(X,1)/10:(size(X,1) + 1));
            cTicks(end) = cTicks(end)-1;
            for u = 1:length(cTicks)
                cTickLabels(u) = cellstr(num2str(cTicks(u)));
            end
            colorbar('Ticks',cTicks,'TickLabels',cTickLabels);
            
            hold off
            
            %title({strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' '));
            %   string(strcat(variableLabel_z,{' vs '},variableLabel_x,{' vs '},variableLabel_y));
            %   string(strcat({'Generational Champion'}))})
            view([50 17])
            
            saveas(...
                figure(fig_num),...
                string(strcat(...
                compilation(p).figurePath,...
                'fitnessBinComparisons_topIndividuals_',field_x, {'_vs_'},field_y, {'_vs_'},field_z,...
                '_Case',...
                convertCharsToStrings(compilation(p).case),...
                '.png')),...
                'png');
            close(fig_num);
        end
    else
        
        combinations = [1 1];
        
        Xx = combinations(1,1);
        field_x = commonFitnessVectorFieldNames{Xx};
        variableLabel_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).variableLabel;
        figureSaveName_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).figureSaveNameContent;
        units_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).units;
        limits_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).limits;
        scalingFactor_x = compilation(p).evolutions(end).population(end).fitnessVector.common.(field_x).scalingFactor;
        
        Xy = combinations(1,2);
        field_y = ascentPhaseFitnessVectorFieldNames{Xy};
        variableLabel_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).variableLabel;
        figureSaveName_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).figureSaveNameContent;
        units_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).units;
        limits_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).limits;
        scalingFactor_y = compilation(p).evolutions(end).population(end).fitnessVector.ascent.(field_y).scalingFactor;
        
        
        fig_num = p*100 + 3466000 + 1;
        figure(fig_num)
        %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
        set (gca,'Fontsize',15)
        
        set(gca,'YScale','log')
        %set(gca,'XScale','log')
        
        xlim(limits_x)
        ylim([1e-2 1e1])
        
        hold on
        grid on
        
        color = linspace(1,size(topIndividual,1),size(topIndividual,1));
        scale = 25*ones(size(topIndividual,1),1);
        scatter(topIndividual(:,1),topIndividual(:,2),scale,color,'Marker',m{1});
        %[k,av] = convhull(topIndividualsVector);
        %[A,b] = prtp(topIndividualsVector);
        %[membership,member_value]=find_pareto_frontier(topIndividualsVector);
        %            plot(member_value(:,1),member_value(:,2),'r');
        
        % plot(topIndividualsVector(k,1),topIndividualsVector(k,2))
        
        
        
        hold off
        
        caxis([1 size(topIndividual,1)])
        % Create tick markers and labels for the color bar.
        %   The color bar is used to visualize the fitness ranking of each point.
        %   Ranking goes from 1 to the number of individuals, with 1 being the
        %   individual with the lowest magnitude of the Fitness Vector.
        cTicks = round(1:size(topIndividual,1)/10:(size(topIndividual,1) + 1));
        cTicks(end) = cTicks(end)-1;
        for u = 1:length(cTicks)
            cTickLabels(u) = cellstr(num2str(cTicks(u)));
        end
        colorbar('Ticks',cTicks,'TickLabels',cTickLabels);
        
        %title({strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' '));string(strcat(variableLabel_y,{' vs '},variableLabel_x))})
        
        xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
        ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
        
        
        saveas(...
            figure(fig_num),...
            string(strcat(...
            compilation(p).figurePath,...
            'fitnessBinComparisons_',field_y, {'_vs_'},field_x,...
            '_Case',...
            convertCharsToStrings(compilation(p).case),...
            '.png')),...
            'png');
        close(fig_num);
        
    end
end




end

