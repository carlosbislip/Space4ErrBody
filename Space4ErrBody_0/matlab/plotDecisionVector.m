function [  ] = plotDecisionVector( compilation )
%% Decision Vector Plotting

% Extract the names of the primary kinds of parameters.
%    They are organized as:
%           Ascent
%           Common
%           Descent
trajectoryPhaseNames = fieldnames(compilation(end).evolutions(end).population(end).decisionVector.parameters);


% Create tick markers and labels for the color bar.
%   The color bar is used to visualize the fitness ranking of each point.
%   Ranking goes from 1 to the number of individuals, with 1 being the
%   individual with the lowest magnitude of the Fitness Vector.
cTicks = 0:compilation(1).rawData.populationSize/10:compilation(1).rawData.populationSize;
%cTicks(end) = cTicks(end)-1;
for u = 1:length(cTicks)
    cTickLabels(u) = cellstr(num2str(cTicks(u)));
end


%%

% Plot the figures
for p = 1:numel(compilation)
    
    % Loop through all
    %for m = 1:length(trajectoryPhaseNames)
    for m = 1
        
        trajectoryPhase = trajectoryPhaseNames{m};
        variableNames = fieldnames(compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase));
        nodes =  compilation(1).rawData(1).nodes;
        nodeList = cellstr([ "one";"two";"three";"four";"five";"six";"seven";"eight";"nine";"ten";"eleven";"twelve" ]);
        
        
        % Loop through each relevant type of parameter
        for l = 1:length(variableNames)
            variable              = variableNames{l};
            variableSubFields     = fieldnames(compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable));
            dataFields            = fieldnames(compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable).data);
            variableLabel         = compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable).variableLabel;
            figureSaveNameContent = compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable).figureSaveNameContent;
            units                 = compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable).units;
            limits                = compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable).limits;
            scalingFactor         = compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable).scalingFactor;
            tick                  = compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable).tick;
            
            
            
            
            field          = variableNames{l};
            field_struct   = compilation(1).evolutions(end).population(end).decisionVector.parameters.(trajectoryPhase).(variable);
            variableLabel  = field_struct.variableLabel;
            figureSaveName = field_struct.figureSaveNameContent;
            units          = field_struct.units;
            limits         = field_struct.limits;
            scalingFactor  = field_struct.scalingFactor;
            tick           = field_struct.tick;
            
            
            
            % Individual Generational History per Parameter
            % Create a single figure for each nodal parameter through all
            % evolutions
            
            % Loop through nodal parameters
            %for j = 1:(length(variableSubFields)-3)
            for jj = 1:length(dataFields)
                
                fig_num = p*10000000 + jj*1;
                figure(fig_num)
                set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
                set (gca,'Fontsize',15)
                
                xlabel(string(strcat('Generation $(-)$')),'Interpreter','latex') % x-axis label
                
                if numel(dataFields) == 1
                    ylabel(string(strcat(variableLabel,{' '},units)),'Interpreter','latex')% y-axis label
                else
                    ylabel(string(strcat(trajectoryPhase,{' : '},variableLabel,{' '}, num2str(jj),{' '},units)),'Interpreter','latex')% y-axis label
                end
                
                if numel(compilation(p).evolutions) > 1
                    max_evolutions = 25*ceil(compilation(p).rawData.generationList(end)/25);
                    xlim([0 max_evolutions])
                    set(gca,'XTick', 0:max_evolutions/10:max_evolutions);
                else
                    max_evolutions = 1;
                    xlim([0 2])
                    set(gca,'XTick', 0:1:max_evolutions);
                    
                end
                
                ylim(limits)
                set(gca,'YTick', limits(1):tick:limits(2));
                set(gca,'TickLabelInterpreter','latex')
                
                hold on
                grid on
                
                % Loop through Generations
                for k = 1:numel(compilation(p).rawData.generationList)
                    
                    % Plot all NON-PRINTED individuals
                    indices = compilation(p).evolutions(k).population(1).indices.nonPrinted;
                    Y = nan(compilation(p).rawData.populationSize,1);
                    X = k*ones(compilation(p).rawData.populationSize,1);
                    color = Y;
                    scaling = Y;
                    
                    for iii = 1:numel(indices)
                        Y(indices(iii))       = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).data.(dataFields{jj});
                        %color(indices(iii))   = nan;
                        scaling(indices(iii)) = (color(indices(iii))/compilation(p).rawData.populationSize)*100;
                    end
                    scatter(X(indices),Y(indices),'+');
                    
                    % Plot all PRINTED individuals
                    indices = compilation(p).evolutions(k).population(1).indices.printed;
                    for iii = 1:numel(indices)
                        Y(indices(iii))       = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).data.(dataFields{jj});
                        color(indices(iii))   = compilation(p).evolutions(k).population(indices(iii)).rankVector;
                        scaling(indices(iii)) = (color(indices(iii))/compilation(p).rawData.populationSize)*100;
                    end
                    x = X(indices);
                    y = Y(indices);
                    c = color(indices);
                    scatter(x,y,70,c,'filled');
                    
                end
                
                colormap(gca,'jet')
                cbh = colorbar;
                cbh.Ticks = caxis( gca ) ;
                cbh.TickLabels = [{'Low'} {'High'}];
                cbh.TickLabelInterpreter = 'latex';
                
                ylabel(cbh,'Rank','Interpreter','latex','Fontsize',15)
                colobarLabelPos = get(get(cbh,'YLabel'),'Position');
                colobarLabelPos(1) = colobarLabelPos(1);% - 3;
                set(get(cbh,'YLabel'),'Position',colobarLabelPos)
                
                if strcmp(trajectoryPhase,'common')
                    saveNameMainTerm = strcat( compilation(p).figurePath, 'decisionVector_',figureSaveNameContent,'_vs_generation',...
                        '_Case', convertCharsToStrings(compilation(p).case), '.png');
                else
                    saveNameMainTerm = strcat( compilation(p).figurePath, 'decisionVector_',figureSaveNameContent,'_',num2str(jj),'_vs_generation',...
                        '_Case', convertCharsToStrings(compilation(p).case), '.png');
                end
                
                hold off
                saveas(figure(fig_num),saveNameMainTerm,'png');
                close(fig_num);
            end
            
            
            %%{
            % Individual Parameter per Generation
            % Create a single figure for all Individuals of each nodal
            % parameter at each evolution
            
            % Loop through Generations
            for k = 1:numel(compilation(p).rawData.generationList)
                
                % Loop through nodal parameters (all dataFields of all variables)
                %for j = 1:(length(variableSubFields)-3)
                for jj = 1:length(dataFields)
                    
                    fig_num = p*10000000 + k*100000 + jj*1;
                    figure(fig_num)
                    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
                    set (gca,'Fontsize',15)
                    %title(strcat(trajectoryPhase,':_{ }',variableLabel,'_{ }', num2str(jj), ' Individuals - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
                    xlim([0 compilation(p).rawData.populationSize])
                    
                    xlabel(string(strcat('Individual $(-)$')),'Interpreter','latex') % x-axis label
                    
                    if numel(dataFields) == 1
                        ylabel(string(strcat(variableLabel,{' '},units)),'Interpreter','latex')% y-axis label
                    else
                        ylabel(string(strcat(trajectoryPhase,{' : '},variableLabel,{' '}, num2str(jj),{' '},units)),'Interpreter','latex')% y-axis label
                    end
                    
                ylim(limits)
                set(gca,'YTick', limits(1):tick:limits(2));
                set(gca,'TickLabelInterpreter','latex')
                hold on
                    grid on
                    
                    % Plot all NON-PRINTED individuals
                    indices = compilation(p).evolutions(k).population(1).indices.nonPrinted;
                    Y = nan(numel(indices),1);
                    color = Y;
                    scaling = Y;
                    for iii = 1:compilation(p).evolutions(k).population(1).size.nonPrinted
                        Y(indices(iii)) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).data.(dataFields{jj});
                        color(indices(iii))   = compilation(p).evolutions(k).population(indices(iii)).rankVector;
                        scaling(indices(iii)) = (color(indices(iii))/compilation(p).rawData.populationSize)*100;
                    end
                    %scatter(indices,Y(indices),scaling(indices),color(indices));
                    scatter(indices,Y(indices),20,'s','filled','MarkerFaceAlpha',.1);
                    
                    % Plot all PRINTED individuals
                    indices = compilation(p).evolutions(k).population(1).indices.printed;
                    Y = nan(numel(indices),1);
                    scaling = Y;
                    for iii = 1:compilation(p).evolutions(k).population(1).size.printed
                        Y(indices(iii)) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).data.(dataFields{jj});
                        color(indices(iii))   = compilation(p).evolutions(k).population(indices(iii)).rankVector;
                        scaling(indices(iii)) = (color(indices(iii))/compilation(p).rawData.populationSize)*100;
                    end
                    
                    y = Y(indices);
                    c = color(indices);
                    scatter(indices,y,70,c,'filled');
                    
                    %scatter(indices,Y(indices),scaling(indices),'r','filled');
                    
                    %caxis([1 compilation(p).rawData.populationSize])
                    %colorbar('Ticks',cTicks,'TickLabels',cTickLabels,'Interpreter','latex');
                    caxis([0 numel(indices)])
                    colormap(gca,'jet')
                    cbh = colorbar;
                    cbh.Ticks = [0 numel(indices)] ;
                    cbh.TickLabels = [{'High'} {'Low'}];
                    cbh.TickLabelInterpreter = 'latex';
                    
                    ylabel(cbh,'Rank','Interpreter','latex','Fontsize',15)
                    colobarLabelPos = get(get(cbh,'YLabel'),'Position');
                    colobarLabelPos(1) = colobarLabelPos(1);% - 3;
                    set(get(cbh,'YLabel'),'Position',colobarLabelPos)
                    
                    hold off
                    
                    if strcmp(trajectoryPhase,'common')
                        saveNameMainTerm = strcat( compilation(p).figurePath, 'decisionVector_',figureSaveNameContent,{'_vs_'},'individual_perGeneration_',num2str(compilation(p).rawData.generationList(k)),...
                            '_Case', convertCharsToStrings(compilation(p).case), '.png');
                    else
                        saveNameMainTerm = strcat( compilation(p).figurePath, 'decisionVector_',figureSaveNameContent,'_',num2str(jj),{'_vs_'},'individual_perGeneration_',num2str(compilation(p).rawData.generationList(k)),...
                            '_Case', convertCharsToStrings(compilation(p).case), '.png');
                    end
                    
                    hold off
                    saveas(figure(fig_num),saveNameMainTerm,'png');
                    close(fig_num);
                    
                    
                    saveas(...
                        figure(fig_num),...
                        strcat(...
                        compilation(p).figurePath,...
                        figureSaveNameContent,num2str(jj),'Individuals_Generation',...
                        num2str(k - 1),...
                        '_Case',...
                        convertCharsToStrings(compilation(p).case),...
                        '.png'),...
                        'png');
                    close(fig_num);
                end
                %end
            end
            
            
            % This section only applies to the ASCENT and DESCENT phases
            if m ~= 2
                % Parameter progression through nodes per Generation
                % Create a single figure for all individuals at each node
                % for each evolution
                
                % Loop through evolutions
                for k = 1:numel(compilation(p).evolutions)
                    
                    if compilation(p).evolutions(k).population(1).size.printed > 0
                        printed = nan(length(dataFields),compilation(p).evolutions(k).population(1).size.printed);
                        scalingPrinted = printed;
                    end
                    
                    fig_num = p*100 + k*1 + 100000;
                    figure(fig_num)
                    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
                    set (gca,'Fontsize',15)
                    title(strcat(trajectoryPhase,':_{ }',variableLabel,' at Node - Generation:_{ }',num2str(k-1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
                    %xlim([1 length(dataFields)])
                    xlabel('Node (-)') % x-axis label
                    ylabel(strcat(trajectoryPhase,':_{ }',variableLabel,'_{ }',units) )% y-axis label
                    set(gca,'XTick', 1:1:length(dataFields));
                    hold on
                    grid on
                    
                    % Loop through each node (dataField)
                    %for j = 1:(length(variableSubFields)-3)
                    for jj = 1:length(dataFields)
                        
                        % Plot all NON-PRINTED individuals
                        indices = compilation(p).evolutions(k).population(1).indices.nonPrinted;
                        Y = nan(compilation(p).evolutions(k).population(1).size.nonPrinted,1);
                        color = Y;
                        scaling = Y;
                        for iii = 1:compilation(p).evolutions(k).population(1).size.nonPrinted
                            Y(iii) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).data.(dataFields{jj});
                            color(iii) = compilation(p).evolutions(k).population(indices(iii)).fitnessVector.common.ranking;
                            scaling(iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.scalingFactor;
                        end
                        scatter(jj*ones(compilation(p).evolutions(k).population(1).size.nonPrinted,1),Y,scaling,color);
                        
                        indices = compilation(p).evolutions(k).population(1).indices.printed;
                        for iii = 1:compilation(p).evolutions(k).population(1).size.printed
                            printed(jj,iii) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).data.(dataFields{jj});
                            scalingPrinted(jj,iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.common.scalingFactor;
                        end
                        
                    end
                    %end
                    
                    % Plot PRINTED individuals at node with connecting lines
                    if compilation(p).evolutions(k).population(1).size.printed > 0
                        
                        for kkk = 1:compilation(p).evolutions(k).population(1).size.printed
                            
                            h = plot(linspace(1,length(dataFields),length(dataFields)),printed(:,kkk));
                            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                            scatter(linspace(1,length(dataFields),length(dataFields)),printed(:,kkk),scalingPrinted(:,kkk),'filled');
                            
                        end
                        
                    end
                    
                    caxis([1 compilation(p).evolutions(k).population(1).size.collective])
                    colorbar('Ticks',cTicks,'TickLabels',cTickLabels);
                    
                    hold off
                    saveas(...
                        figure(fig_num),...
                        strcat(...
                        compilation(p).figurePath,...
                        figureSaveNameContent,'IndividualsAtNode_Generation',...
                        num2str(k - 1),...
                        '_Case',...
                        convertCharsToStrings(compilation(p).case),...
                        '.png'),...
                        'png');
                    close(fig_num);
                end
                
            end
            %}
        end
    end
end

end