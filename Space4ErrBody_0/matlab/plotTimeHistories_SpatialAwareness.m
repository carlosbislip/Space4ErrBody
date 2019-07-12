function [  ] = plotTimeHistories_SpatialAwareness( compilation )

%% Time History: Cumulative Cartesian Distance Travelled - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 833000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat(' Cumulative Cartesian Distance Travelled through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            ylim([0 8000])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel(' Cumulative Cartesian Distance Travelled (km)') % y-axis label
            set(gca,'YTick', 0:500:8000);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                h = stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.cumulativeCartesianDistanceTravelled/1e3);
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.cumulativeCartesianDistanceTravelled(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/1e3,'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.cumulativeCartesianDistanceTravelled(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/1e3,'s');
            end
            %plot([0 4000],(25)*[1 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryCumulativeDistanceTravelled_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).set),...
                '.png'),...
                'png');
            %    close(fig_num);
            
        end
    end
end

%% Time History: Height - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*200 + 823000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            ylim([10 130])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 10:10:130);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            plot([0 max_tof],(10)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                if  compilation(1).validation == 1
                    stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    ylim([20 130])
                    set(gca,'YTick', 20:10:130);
                else
                    h = stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/1e3,'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/1e3,'s');
               end
                
            end
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryHeight_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).set),...
                '.png'),...
                'png');
            %    close(fig_num);
            
        end
    end
end

%% Time History: Distance Traveled - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 723000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Distance Traveled through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            max_tof = max([compilation(p).evolutions.max_tof]);
            ylim([0 60])
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Distance Traveled (deg)') % y-axis label
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceTraveled);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryAngularDistanceTraveled_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).set),...
                '.png'),...
                'png');
            %      close(fig_num);
        end
    end
end

%% Time History: Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 724000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Distance To Go through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            ylim([0 60])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Distance To Go (deg)') % y-axis label
            set(gca,'YTick', 0:10:60);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            plot([0 max_tof],(.75)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                if  compilation(1).validation == 1
                    stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    plot([0 1400],(.75)*[1 1],'k','LineWidth',2)
                    
                else
                    h = stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo);
                 set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
 end
                
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryAngularDistanceToGo_v_T_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).set),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Time History: Distance Covered Ratio - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 725000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Distance Covered Ratio - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Distance Covered Ratio (-)') % y-axis label
        set(gca,'YTick', 0:.1:1);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],(.4)*[1 1],'k','LineWidth',2)
        plot([0 max_tof],(.1)*[1 1],'k','LineWidth',2)
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            %distanceToCover = compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(1);
            %distanceCoveredRatio = compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceTraveled/distanceToCover;
            %h = stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
            %   distanceCoveredRatio);
            h = stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceCoveredRatio);
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
            
            terminationDistanceRatio = compilation(p).evolutions(k).population(ii).decisionVector.parameters.Common.terminationDistanceRatio.data.one;
            I = find( compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceCoveredRatio > terminationDistanceRatio );
            scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(I(1)),terminationDistanceRatio);
            'e';
            
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/timeHistoryDistanceCoveredRatio_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end


%% Time History: Heading Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 783000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 360])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading Angle (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngle);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/timeHistoryHeadingAngle_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Time History: Heading to Target - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 727000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading to Target through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 360])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading to Target (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingToTarget);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/timeHistoryHeadingToTarget_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
        
    end
end


%% Time History: Heading Error - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 734000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Error through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading Error (deg)') % y-axis label
        set(gca,'YTick', -30:5:30);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingError);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/timeHistoryHeadingError_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end
end