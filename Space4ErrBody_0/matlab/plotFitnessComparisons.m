function [  ] = plotFitnessComparisons( compilation )
disp('   Fitness Comparisons')

%%


for p = 1:numel(compilation)
    if numel(compilation(p).evolutions) > 1
        
        fitnessMagnitudeAverage                   = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageTop                = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageTopChange          = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageOfAverageTop       = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageOfAverageTopChange = nan(numel(compilation(p).evolutions),1);
        
        for k = 1:numel(compilation(p).evolutions)
            
            fitnessMagnitudeAverage(k)                   = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,1)));
            fitnessMagnitudeAverageTop(k)                = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,2)));
            fitnessMagnitudeAverageTopChange(k)          = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,3)));
            fitnessMagnitudeAverageOfAverageTop(k)       = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,4)));
            fitnessMagnitudeAverageOfAverageTopChange(k) = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,5)));
            
        end
        
        fig_num = p*10000;
        figure(fig_num)
        %set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
        
        semilogy(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),fitnessMagnitudeAverage);
        hold on
        
        title_line1 = strip(strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' ')),'left');
        title_line2 = 'Average Fitness Magnitude per Generation';
        % title({title_line1;title_line2})
        
        if numel(compilation(p).evolutions) > 1
            max_evolutions = 25*ceil(max((numel(compilation(p).evolutions)-1))/25);
            xlim([0 max_evolutions])
            set(gca,'XTick', 0:max_evolutions/10:max_evolutions);
        end
        
        xlabel('Generation (-)','Interpreter','latex') % x-axis label
        ylabel('Average Fitness Magnitude (-)','Interpreter','latex')% y-axis label
        set (gca,'Fontsize',15)
        grid on
        
        semilogy(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),fitnessMagnitudeAverageTop);
        semilogy(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),fitnessMagnitudeAverageTopChange);
        semilogy(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),fitnessMagnitudeAverageOfAverageTop);
        semilogy(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),fitnessMagnitudeAverageOfAverageTopChange);
        
        hold off
        
        legend('Population Average','Top Average','Change of Top Average','Cumulative Average of Top Average','Change of Cumulative Average of Top Average','Location','southwest')
        
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'fitnessMagnitudeAverageHistory',...
            '_Case',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        
    end
end

%%


for p = 1:numel(compilation)
    if numel(compilation(p).evolutions) > 1
        
        fitnessMagnitudeAverage                   = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageTop                = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageTopChange          = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageOfAverageTop       = nan(numel(compilation(p).evolutions),1);
        fitnessMagnitudeAverageOfAverageTopChange = nan(numel(compilation(p).evolutions),1);
        
        for k = 1:numel(compilation(p).evolutions)
            
            fitnessMagnitudeAverage(k)                   = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,1)));
            fitnessMagnitudeAverageTop(k)                = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,2)));
            fitnessMagnitudeAverageTopChange(k)          = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,3)));
            fitnessMagnitudeAverageOfAverageTop(k)       = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,4)));
            fitnessMagnitudeAverageOfAverageTopChange(k) = str2double(string(compilation(p).evolutions(k).population(1).fitnessVector.commonW.magnitudeAverages(:,5)));
            
        end
        
        fig_num = p*10000;
        figure(fig_num)
        %set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
        
        semilogy(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),fitnessMagnitudeAverageTop);
        hold on
        
        title_line1 = strip(strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' ')),'left');
        title_line2 = 'Top Average Fitness Magnitude per Generation';
        % title({title_line1;title_line2})
        
        if numel(compilation(p).evolutions) > 1
            max_evolutions = 25*ceil(max((numel(compilation(p).evolutions)-1))/25);
            xlim([0 max_evolutions])
            set(gca,'XTick', round(0:max_evolutions/10:max_evolutions));
        end
        
        xlabel('Generation (-)','Interpreter','latex') % x-axis label
        ylabel('Average Fitness Magnitude (-)','Interpreter','latex')% y-axis label
        %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        ylim([0 1000])
        set(gca,'YScale','log')
        
        grid on
        
        semilogy(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),fitnessMagnitudeAverageOfAverageTop);
        
        hold off
        
        legend('Top Average','Cumulative Average of Top Average')
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'fitnessMagnitudeAverageTopHistory',...
            '_Case',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        
    end
    
    
end


%%

for p = 1:numel(compilation)
    
    validation = compilation(p).validation;
    trajectoryType = compilation(p).rawData.trajectoryType;
    
    if validation == 0
        nodes = compilation(p).rawData.nodes;
        nodeList = cellstr([ "one";"two";"three";"four";"five";"six";"seven";"eight";"nine";"ten";"eleven";"twelve" ]);
        m = {'+', '<', 'o', '*','square','.', 'diamond','x', 'v', '^', '>', 'pentagram'};
        commonFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).fitnessVector.common);
        ascentPhaseFitnessVectorFieldNames = [];
        descentPhaseFitnessVectorFieldNames = [];
        
        if isfield(compilation(p).evolutions(end).population(end).fitnessVector,'ascent') || isfield(compilation(p).evolutions(end).population(end).fitnessVector,'descent')
            
            if trajectoryType == 'A'
                %interpolatorsInputDataFields  = fieldnames(compilation(p).evolutions(end).population.interpolators.Ascent.InputData);
                %interpolatorsEvaluationFields = fieldnames(compilation(p).evolutions(end).population.interpolators.Ascent.Evaluation);
                %decisionVectorAscentPhaseParameterNames  = fieldnames(compilation(p).evolutions(end).population.decisionVector.parameters.Ascent);
                ascentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).fitnessVector.ascent);
            end
            if trajectoryType == 'D'
                interpolatorsInputDataFields  = fieldnames(compilation(p).evolutions(end).population.interpolators.Descent.InputData);
                interpolatorsEvaluationFields = fieldnames(compilation(p).evolutions(end).population.interpolators.Descent.Evaluation);
                decisionVectorDescentPhaseParameterNames  = fieldnames(compilation(p).evolutions(end).population.decisionVector.parameters.Descent);
                descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population.fitnessVector.descent);
            end
            if trajectoryType == 'AD'
                interpolatorsInputDataFields  = fieldnames(compilation(p).evolutions(end).population.interpolators.Ascent.InputData);
                interpolatorsEvaluationFields = fieldnames(compilation(p).evolutions(end).population.interpolators.Ascent.Evaluation);
                decisionVectorAscentPhaseParameterNames  = fieldnames(compilation(p).evolutions(end).population.decisionVector.parameters.ascent);
                decisionVectorDescentPhaseParameterNames  = fieldnames(compilation(p).evolutions(end).population.decisionVector.parameters.descent);
                ascentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population.fitnessVector.ascent);
                descentPhaseFitnessVectorFieldNames  = fieldnames(compilation(p).evolutions(end).population.fitnessVector.descent);
            end
        end
    end
    
    X = nan(compilation(p).evolutions(1).population(1).size.collective,1 + length(ascentPhaseFitnessVectorFieldNames),numel(compilation(p).evolutions));
    %X2 = nan(length(compilation(p).evolutions(1).population(1).indices.printed),1 + length(fitnessVectorFieldNames),numel(compilation(p).evolutions));
    fillerZ = ones(compilation(p).evolutions(1).population(1).size.collective,1);
    
    
    topIDVector = nan(numel(compilation(p).evolutions),1);
    topIDMagnitude = nan(numel(compilation(p).evolutions),1);
    topIndividualsVector = nan(numel(compilation(p).evolutions),length(commonFitnessVectorFieldNames) + length(ascentPhaseFitnessVectorFieldNames));
    topIndividualsMagnitude = nan(numel(compilation(p).evolutions),length(commonFitnessVectorFieldNames) + length(ascentPhaseFitnessVectorFieldNames));
    
    
    for k = 1:numel(compilation(p).evolutions)
        
        i = 1;
        while (1 ~= compilation(p).evolutions(k).population(i).fitnessVector.RankBySortedVector)
            i = i + 1;
        end
        topIDVector(k) = i;
        
        for ll = 1:length(commonFitnessVectorFieldNames)
            topIndividualsVector(k,ll) = compilation(p).evolutions(k).population(i).fitnessVector.common.(commonFitnessVectorFieldNames{ll}).value;
        end
        for ll = 1:length(ascentPhaseFitnessVectorFieldNames)
            topIndividualsVector(k,length(commonFitnessVectorFieldNames) + ll) = compilation(p).evolutions(k).population(i).fitnessVector.ascent.(ascentPhaseFitnessVectorFieldNames{ll}).value;
        end
        for ll = (length(ascentPhaseFitnessVectorFieldNames)+1):2*length(descentPhaseFitnessVectorFieldNames)
            topIndividualsVector(k,length(commonFitnessVectorFieldNames) + length(ascentPhaseFitnessVectorFieldNames) + ll) = compilation(p).evolutions(k).population(i).fitnessVector.descent.(descentPhaseFitnessVectorFieldNames{ll}).value;
        end
        
        i = 1;
        while (1 ~= compilation(p).evolutions(k).population(i).fitnessVector.RankBySortedMagnitude)
            i = i + 1;
        end
        topIDMagnitude(k) = i;
        
        for ll = 1:length(commonFitnessVectorFieldNames)
            topIndividualsMagnitude(k,ll) = compilation(p).evolutions(k).population(i).fitnessVector.common.(commonFitnessVectorFieldNames{ll}).value;
        end
        for ll = 1:length(ascentPhaseFitnessVectorFieldNames)
            topIndividualsMagnitude(k,length(commonFitnessVectorFieldNames) + ll) = compilation(p).evolutions(k).population(i).fitnessVector.ascent.(ascentPhaseFitnessVectorFieldNames{ll}).value;
        end
        for ll = (length(ascentPhaseFitnessVectorFieldNames)+1):2*length(descentPhaseFitnessVectorFieldNames)
            topIndividualsMagnitude(k,length(commonFitnessVectorFieldNames) + length(ascentPhaseFitnessVectorFieldNames) + ll) = compilation(p).evolutions(k).population(i).fitnessVector.descent.(descentPhaseFitnessVectorFieldNames{ll}).value;
        end
    end
    
    
    if trajectoryType == 'A'
        for k = 1:numel(compilation(p).evolutions)
            for ii = 1:compilation(p).evolutions(k).population(1).size.collective
                X(ii,1,k) = compilation(p).evolutions(k).population(ii).fitnessVector.common.(commonFitnessVectorFieldNames{1}).value;
                for j = 1:length(ascentPhaseFitnessVectorFieldNames)
                    X(ii,j+1,k) = compilation(p).evolutions(k).population(ii).fitnessVector.ascent.(ascentPhaseFitnessVectorFieldNames{j}).value;
                end
            end
            %             for i = 1:length(compilation(p).evolutions(k).population(1).indices.printed)
            %                 X2(i,1)= compilation(p).evolutions(k).population(compilation(p).evolutions(k).population(1).indices.printed(i)).fitnessVector.common.angularDistanceToGo;
            %                 for j = 1:length(fitnessVectorFieldNames)
            %                     X2(i,j+1,k) = compilation(p).evolutions(k).population(compilation(p).evolutions(k).population(1).indices.printed(i)).fitnessVector.Ascent.(fitnessVectorFieldNames{j});
            %                 end
            %             end
        end
        
        if isempty(ascentPhaseFitnessVectorFieldNames) == true
            if size(topIndividualsVector,2) > 1
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
                    ppp(1) = scatter(linspace(1,numel(compilation(p).evolutions),numel(compilation(p).evolutions)),topIndividualsVector,'Marker',m{1},'DisplayName',string(strcat({'Champion: Vector'})));
                    
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
                
                X = topIndividualsVector;
                
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
                    
                    color = linspace(1,size(topIndividualsVector,1),size(topIndividualsVector,1));
                    scale = 25*ones(size(topIndividualsVector,1),1);
                    
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
                
                color = linspace(1,size(topIndividualsVector,1),size(topIndividualsVector,1));
                scale = 25*ones(size(topIndividualsVector,1),1);
                scatter(topIndividualsVector(:,1),topIndividualsVector(:,2),scale,color,'Marker',m{1});
                %[k,av] = convhull(topIndividualsVector);
                %[A,b] = prtp(topIndividualsVector);
                %[membership,member_value]=find_pareto_frontier(topIndividualsVector);
                %            plot(member_value(:,1),member_value(:,2),'r');
                
                % plot(topIndividualsVector(k,1),topIndividualsVector(k,2))
                
                
                
                hold off
                
                caxis([1 size(topIndividualsVector,1)])
                % Create tick markers and labels for the color bar.
                %   The color bar is used to visualize the fitness ranking of each point.
                %   Ranking goes from 1 to the number of individuals, with 1 being the
                %   individual with the lowest magnitude of the Fitness Vector.
                cTicks = round(1:size(topIndividualsVector,1)/10:(size(topIndividualsVector,1) + 1));
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
    elseif trajectoryType == 'D'
        for i = 1:length(ascentPhaseFitnessVectorFieldNames)
            %    evolutions(k).population(p).fitnessVector.Descent.(fitnessVectorFieldNames{i}) = pureFitness{p,i};
        end
        
    elseif trajectoryType == 'AD'
        for i = 1:length(ascentPhaseFitnessVectorFieldNames)
            %        evolutions(k).population(p).fitnessVector.Ascent.(fitnessVectorFieldNames{i}) = pureFitness{p,i};
            %        evolutions(k).population(p).fitnessVector.Descent.(fitnessVectorFieldNames{i}) = pureFitness{p,i +length(fitnessVectorFieldNames)};
        end
    end
    
end



end

