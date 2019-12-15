function [  ] = plotTimeHistories_SpatialAwareness( compilation )

%% Time History: Cumulative Cartesian Distance Travelled - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 833000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat(' Cumulative Cartesian Distance Travelled through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
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
            plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'cumulativeCartesianDistanceTravelled'}, 1e3 )

         end
        %plot([0 4000],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryCumulativeDistanceTravelled_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %    close(fig_num);
        
    end
end

%% Time History: Height - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            
            fig_num = p*200 + 823000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([10 200])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:20:200);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(20)*[1 1],'k','LineWidth',2)
            
            if  compilation(1).validation == 1
                for ii = compilation(p).evolutions(k).population(1).indices.printed
                    
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    ylim([20 130])
                    set(gca,'YTick', 20:10:130);
                end
            else
                plotFieldVsField( compilation, p, k, compilation(p).evolutions(k).population(1).indices.printed , {'timeOfFlight'}, 1, {'height'}, 1e3 )
                
            end
                
                hold off
                
                saveas(...
                    figure(fig_num),...
                    strcat(...
                compilation(p).figurePath,...
                'timeHistoryHeight_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Time History: Distance Traveled - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 723000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Distance Traveled through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            max_tof = max([compilation(p).evolutions.max_tof]);
            ylim([0 60])
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Distance Traveled (deg)') % y-axis label
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                        plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'angularDistanceTraveled'}, 1 )
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryAngularDistanceTraveled_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %      close(fig_num);
        end
    end
end

%% Time History: Central Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 724000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Central Target Angular Distance To Go through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
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
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    plot([0 1400],(.75)*[1 1],'k','LineWidth',2)
                    
                else
                                     plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'centralTargetAngularDistanceToGo'}, 1 )

                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryCentralTargetAngularDistanceToGo_v_T_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

%% Time History: Dynamic Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 724000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Dynamic Target Angular Distance To Go through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
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
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                plot([0 1400],(.75)*[1 1],'k','LineWidth',2)
                
            else
                                                   plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'dynamicTargetAngularDistanceToGo'}, 1 )

            end
        end
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryDynamicTargetAngularDistanceToGo_v_T_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %  close(fig_num);
    end
end


%% Time History: Distance Covered Ratio - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 725000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Distance Covered Ratio - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
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
            plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'angularDistanceCoveredRatio'}, 1 )
        end
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryDistanceCoveredRatio_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end


%% Time History: Heading Angle - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 783000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Angle through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
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
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngle);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryHeadingAngle_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Heading to Central Target - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 727000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading to Central through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([0 360])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading to Central Target (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToCentralTarget);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryHeadingToCentralTarget_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        
    end
end


%% Time History: Intersection Case - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 927000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Intersection Case through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([0 20])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Intersection Case (-)') % y-axis label
        set(gca,'YTick', 0:1:20);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.intersectionCase);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryIntersectionCase_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        
    end
end


%% Time History: Heading to Dynamic Target - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 727000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading to Dynamic Target through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([0 360])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading to Dynamic Target (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToDynamicTarget);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryHeadingToDynamicTarget_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        
    end
end

%% Time History: Dynamic Target Heading Bounds - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)

                if compilation(p).evolutions(k).population(1).indices.printed > 0

        fig_num = p*100 + 727000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Dynamic Target Heading Angle Bounds through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([0 360])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Dynamic Target Heading Angle Bounds (deg)') % y-axis label
        set(gca,'YTick', 0:30:360);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight.value,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetHeadingAngleBounds1.value);
            set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)));
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight.value,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetHeadingAngleBounds2.value);
            
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight.value,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngle.value);
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight.value,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToDynamicTarget.value);
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight.value,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToDynamicTargetError.value);
            
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryDynamicTargetHeadingAngleBounds_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        
    end
    end
end

%% Time History: Heading to Dynamic Target Error - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 734000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading to Dynamic Target Error through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([-30 30])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Heading to Dynamic Target Error (deg)') % y-axis label
        set(gca,'YTick', -30:5:30);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToDynamicTargetError);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryHeadingToDynamicTargetError_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end

%% Time History: Dynamic Target: Latitude - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)

                if compilation(p).evolutions(k).population(1).indices.printed > 0

                    fig_num = p*100 + 735000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Dynamic Target: Latitude through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        %yLimits = [ (compilation(1).centraTargetCoordinates(1) - compilation(p).angularDistanceForTermination), (compilation(1).centraTargetCoordinates(1) + compilation(p).angularDistanceForTermination)];
        %ylim(yLimits)
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Dynamic Target: Longitude (deg)') % y-axis label
        %set(gca,'YTick', yLimits(1):.2:yLimits(2));
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            %for ii = 19
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight.value,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetLatitude.value);
        end
        
        plot([0 max_tof],compilation(1).centraTargetCoordinates(1)*[1 1],'k','LineWidth',2)
        plot([0 max_tof],(compilation(1).centraTargetCoordinates(1) + compilation(p).angularDistanceForTermination)*[1 1],'k','LineWidth',2)
        plot([0 max_tof],(compilation(1).centraTargetCoordinates(1) - compilation(p).angularDistanceForTermination)*[1 1],'k','LineWidth',2)
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryDynamicTargetLatitude_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
    end
end


%% Time History: Dynamic Target: Longitude - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 736000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Dynamic Target: Longitude through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        %  yLimits = [ (compilation(1).centraTargetCoordinates(2) - compilation(p).angularDistanceForTermination), (compilation(1).centraTargetCoordinates(2) + compilation(p).angularDistanceForTermination)];
        %ylim(yLimits)
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Dynamic Target: Longitude (deg)') % y-axis label
        %set(gca,'YTick', yLimits(1):.2:yLimits(2));
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetLongitude);
            %set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)));
            
            %plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
            %   compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.intersectionCase);
            
        end
        
        plot([0 max_tof],compilation(1).centraTargetCoordinates(2)*[1 1],'k','LineWidth',2)
        plot([0 max_tof],-76.491504280000001*[1 1],'k','LineWidth',2)
        plot([0 max_tof],-78.420162379999994*[1 1],'k','LineWidth',2)
        
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryDynamicTargetLongitude_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Vertex Coordinates: Latitude - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 735000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Vertex Coordinates: Latitude through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        ylim([0 90])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Vertex Coordinates: Longitude (deg)') % y-axis label
        set(gca,'YTick', 0:15:90);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.vertexLatitude);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryVertexLatitude_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end


%% Time History: Vertex Coordinates: Longitude - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 735000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Vertex Coordinates: Longitude through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        ylim([-180 180])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Vertex Coordinates: Longitude (deg)') % y-axis label
        set(gca,'YTick', -180:30:180);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.vertexLongitude);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryVertexLongitude_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end
%% Time History: Flight Corridor - Dynamic Pressure - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 735000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Time History: Flight Corridor - Dynamic Pressure - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([10 200])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:20:200);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.entryCorridorBoundary_DynamicPressure/1e3);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryFlightCorridorDynamicPressure_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end
%% Time History: Flight Corridor - Skip Suppression - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 735000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Time History: Flight Corridor - Skip Suppression - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([10 200])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:20:200);
        set(gca,'XTick', 0:200:max_tof);
        
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_SkipSuppression/1e3);
        end
        
        %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
        
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryFlightCorridorDynamicPressure_Generation_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end




end