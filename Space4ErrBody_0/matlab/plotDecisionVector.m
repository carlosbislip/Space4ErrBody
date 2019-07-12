function [  ] = plotDecisionVector( compilation )
%% Decision Vector Plotting

% Extract the names of the primary kinds of parameters.
%    They are organized as: 
%           Ascent
%           Common
%           Descent 
trajectoryPhaseNames = fieldnames(compilation(1).evolutions(1).population(1).decisionVector.parameters);


% Create tick markers and labels for the color bar.
%   The color bar is used to visualize the fitness ranking of each point.
%   Ranking goes from 1 to the number of individuals, with 1 being the
%   individual with the lowest magnitude of the Fitness Vector.
cTicks = 1:compilation(1).evolutions(1).population(1).size.collective/10:(compilation(1).evolutions(1).population(1).size.collective + 1);
cTicks(end) = cTicks(end)-1;
for u = 1:length(cTicks)
    cTickLabels(u) = cellstr(num2str(cTicks(u)));
end


% Plot the figures
for p = 1:numel(compilation)
    
    % Loop through Ascent and Descent ONLY
    for m = 1:length(trajectoryPhaseNames)
        
        trajectoryPhase = trajectoryPhaseNames{m};
        variableNames = fieldnames(compilation(1).evolutions(1).population(1).decisionVector.parameters.(trajectoryPhase));
   
        % Loop through each relevant type of parameter
        for l = 1:length(variableNames)
            variable              = variableNames{l};
            variableSubFields     = fieldnames(compilation(1).evolutions(1).population(1).decisionVector.parameters.(trajectoryPhase).(variable));
            dataFields            = fieldnames(compilation(1).evolutions(1).population(1).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{4}));
            figureTitleContent    = compilation(1).evolutions(1).population(1).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{1});
            figureSaveNameContent = compilation(1).evolutions(1).population(1).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{2});
            units                 = compilation(1).evolutions(1).population(1).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{3});
            
            % Individual Evolutionary History per Parameter
            % Create a single figure for each nodal parameter through all
            % evolutions
            
            % Loop through nodal parameters
            for j = 1:(length(variableSubFields)-3)
                for jj = 1:length(dataFields)
                    
                    fig_num = p*10000000 + j*1;
                    figure(fig_num)
                    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
                    set (gca,'Fontsize',15)
                    title(strcat(trajectoryPhase,':_{ }',figureTitleContent,'_{ }', num2str(jj), ' per Evolution:_{ }',' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
                    
                    if numel(compilation(p).evolutions) > 1
                        xlim([0 (numel(compilation(p).evolutions)-1)])
                    end
                    xlabel('Evolution (-)') % x-axis label
                    ylabel(strcat(trajectoryPhase,':_{ }',figureTitleContent,'_{ }', num2str(jj),'_{ }',units) )% y-axis label
                    set(gca,'XTick', 0:1:(numel(compilation(p).evolutions)-1));
                    hold on
                    grid on
                    
                    % Loop through Evolutions
                    for k = 1:numel(compilation(p).evolutions)
                        
                        % Plot all NON-PRINTED individuals
                        indices = compilation(p).evolutions(k).population(1).indices.nonPrinted;
                        Y = nan(compilation(p).evolutions(k).population(1).size.nonPrinted,1);
                        color = Y;
                        scaling = Y;
                        for iii = 1:compilation(p).evolutions(k).population(1).size.nonPrinted
                            Y(iii)       = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{ j + 3 }).(dataFields{jj});
                            color(iii)   = compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.ranking;
                            scaling(iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.scalingFactor;
                        end
                        scatter((k - 1 )*ones(compilation(p).evolutions(k).population(1).size.nonPrinted,1),Y,scaling,color);
                        
                        % Plot all PRINTED individuals
                        if compilation(p).evolutions(k).population(1).size.printed > 0
                            
                            indices = compilation(p).evolutions(k).population(1).indices.printed;
                            Y = nan(compilation(p).evolutions(k).population(1).size.printed,1);
                            scaling = Y;
                            for iii = 1:compilation(p).evolutions(k).population(1).size.printed
                                Y(iii)       = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{ j + 3 }).(dataFields{jj});
                                scaling(iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.scalingFactor;
                            end
                            
                            scatter((k - 1 )*ones(compilation(p).evolutions(k).population(1).size.printed,1),Y,scaling,'r','filled');
                            
                        end
                        
                    end
                    
                    caxis([1 compilation(p).evolutions(k).population(1).size.collective])
                    colorbar('Ticks',cTicks,'TickLabels',cTickLabels);
                    
                    hold off
                    saveas(...
                        figure(fig_num),...
                        strcat(...
                        compilation(p).mainpath,...
                        '/figures/',figureSaveNameContent,num2str(jj),'PerEvolution',...
                        '_Set',...
                        convertCharsToStrings(compilation(p).set),...
                        '.png'),...
                        'png');
                    close(fig_num);
                end
            end
            
            % Individual Parameter per Evolution
            % Create a single figure for all Individuals of each nodal
            % parameter at each evolution
            
            % Loop through Evolutions
            for k = 1:numel(compilation(p).evolutions)
                
                % Loop through nodal parameters (all dataFields of all variables)
                for j = 1:(length(variableSubFields)-3)
                    for jj = 1:length(dataFields)
                        
                        fig_num = p*10000000 + k*100000 + j*1;
                        figure(fig_num)
                        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
                        set (gca,'Fontsize',15)
                        title(strcat(trajectoryPhase,':_{ }',figureTitleContent,'_{ }', num2str(jj), ' Individuals - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
                        xlim([0 compilation(p).evolutions(k).population(1).size.collective])
                        xlabel('Individual (-)') % x-axis label
                        ylabel(strcat(trajectoryPhase,':_{ }',figureTitleContent,'_{ }', num2str(jj),'_{ }',units) )% y-axis label
                        hold on
                        grid on
                        
                        % Plot all NON-PRINTED individuals
                        indices = compilation(p).evolutions(k).population(1).indices.nonPrinted;
                        Y = nan(compilation(p).evolutions(k).population(1).size.nonPrinted,1);
                        color = Y;
                        scaling = Y;
                        for iii = 1:compilation(p).evolutions(k).population(1).size.nonPrinted
                            Y(iii) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{ j + 3 }).(dataFields{jj});
                            color(iii) = compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.ranking;
                            scaling(iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.scalingFactor;
                        end
                        scatter(indices,Y,scaling,color);
                        
                        % Plot all PRINTED individuals
                        if compilation(p).evolutions(k).population(1).size.printed > 0
                            
                            indices = compilation(p).evolutions(k).population(1).indices.printed;
                            Y = nan(compilation(p).evolutions(k).population(1).size.printed,1);
                            scaling = Y;
                            for iii = 1:compilation(p).evolutions(k).population(1).size.printed
                                Y(iii) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{ j + 3 }).(dataFields{jj});
                                scaling(iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.scalingFactor;
                            end
                            
                            scatter(indices,Y,scaling,'r','filled');
                            
                        end
                        
                        caxis([1 compilation(p).evolutions(k).population(1).size.collective])
                        colorbar('Ticks',cTicks,'TickLabels',cTickLabels);
                        
                        hold off
                        saveas(...
                            figure(fig_num),...
                            strcat(...
                            compilation(p).mainpath,...
                            '/figures/',figureSaveNameContent,num2str(jj),'Individuals_Evolution',...
                            num2str(k - 1),...
                            '_Set',...
                            convertCharsToStrings(compilation(p).set),...
                            '.png'),...
                            'png');
                        close(fig_num);
                    end
                end
            end
            
            
            % This section only applies to the ASCENT and DESCENT phases
            if m ~= 2
                % Parameter progression through nodes per Evolution
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
                    title(strcat(trajectoryPhase,':_{ }',figureTitleContent,' at Node - Evolution:_{ }',num2str(k-1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
                    xlim([1 length(dataFields)])
                    xlabel('Node (-)') % x-axis label
                    ylabel(strcat(trajectoryPhase,':_{ }',figureTitleContent,'_{ }',units) )% y-axis label
                    set(gca,'XTick', 1:1:length(dataFields));
                    hold on
                    grid on
                    
                    % Loop through each node (dataField)
                    for j = 1:(length(variableSubFields)-3)
                        for jj = 1:length(dataFields)
                            
                            % Plot all NON-PRINTED individuals
                            indices = compilation(p).evolutions(k).population(1).indices.nonPrinted;
                            Y = nan(compilation(p).evolutions(k).population(1).size.nonPrinted,1);
                            color = Y;
                            scaling = Y;
                            for iii = 1:compilation(p).evolutions(k).population(1).size.nonPrinted
                                Y(iii) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{ j + 3 }).(dataFields{jj});
                                color(iii) = compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.ranking;
                                scaling(iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.scalingFactor;
                            end
                            scatter(jj*ones(compilation(p).evolutions(k).population(1).size.nonPrinted,1),Y,scaling,color);
                            
                            indices = compilation(p).evolutions(k).population(1).indices.printed;
                            for iii = 1:compilation(p).evolutions(k).population(1).size.printed
                                printed(jj,iii) = compilation(p).evolutions(k).population(indices(iii)).decisionVector.parameters.(trajectoryPhase).(variable).(variableSubFields{ j + 3 }).(dataFields{jj});
                                scalingPrinted(jj,iii) = (100/compilation(p).evolutions(k).population(1).size.collective)*compilation(p).evolutions(k).population(indices(iii)).fitnessVector.Common.scalingFactor;
                            end
                            
                        end
                    end
                    
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
                        compilation(p).mainpath,...
                        '/figures/',figureSaveNameContent,'IndividualsAtNode_Evolution',...
                        num2str(k - 1),...
                        '_Set',...
                        convertCharsToStrings(compilation(p).set),...
                        '.png'),...
                        'png');
                    close(fig_num);
                end
                
            end
        end
    end
end

end