function [  ] = plotTimeHistories_Mechanical( compilation )


%% Time History: Total Body G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3957000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Total Body G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            %ylim([0 5])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Total Body G-load (g_0)') % y-axis label
            %set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
           % plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
%             for ii = compilation(p).evolutions(k).population(1).indices.printed
%                 h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
%                     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude);
%                 set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
%                 
%                 if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
%                     scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
%                         compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
%                     set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
%                 end
%                 
%                 if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
%                     scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
%                         compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
%                 end
%             end
            
            
             ii = compilation(p).evolutions(k).population(1).indices.printed;
                
            plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'bodyFrameTotalGLoadMagnitude'}, 1 )
                
            
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryBodyFixedTotalGLoad_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %  close(fig_num);
        end
    end
end


%% Time History: Body Frame X-Component Acceleration - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*400 + 6977000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Body Frame X-Component Acceleration through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-10 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Body Frame X-Component Acceleration (m/s^2)') % y-axis label
            %set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            field = {'bodyFrameTotalAcceleration_x'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryBodyFrameAccelerationX_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end

%% Time History: Passenger Frame X-Component Acceleration - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*400 + 3977000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Passenger Frame X-Component Acceleration through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-10 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Passenger Frame X-Component Acceleration (m/s^2)') % y-axis label
            %set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            
            field = {'passengerFrameTotalAcceleration_x'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryPassengerFrameAccelerationX_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end

%% Time History: Passenger Frame Z-Component Acceleration - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*400 + 3277000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Passenger Frame Z-Component Acceleration through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-10 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Passenger Frame Z-Component Acceleration (m/s^2)') % y-axis label
            %set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            field = {'passengerFrameTotalAcceleration_z'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryPassengerFrameAccelerationZ_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end

%% Time History: Passenger Frame X-Component Jerk - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*400 + 5967000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Passenger Frame X-Component Jerk through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-10 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Passenger Frame X-Component Jerk (m/s^3)') % y-axis label
            %set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            field = {'passengerFrameJerk_x_calc'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryPassengerFrameJerkX_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end

%% Time History: Passenger Frame Y-Component Jerk - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*400 + 3967000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Passenger Frame Y-Component Jerk through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-10 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Passenger Frame Y-Component Jerk (m/s^3)') % y-axis label
            %set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            
            field = {'passengerFrameJerk_y_calc'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryPassengerFrameJerkY_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end

%% Time History: Passenger Frame Z-Component Jerk - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*400 + 3967000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Passenger Frame Z-Component Jerk through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-10 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Passenger Frame Z-Component Jerk (m/s^3)') % y-axis label
            %set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            
            field = {'passengerFrameJerk_z_calc'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryPassengerFrameJerkZ_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end

%% Time History: Total Body Z-Component G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*200 + 3457000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Total Body Z-Component G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-5 5])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Total Body Z-Component G-load (g)') % y-axis label
            set(gca,'YTick', -5:1:10);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            field = {'bodyFrameTotalGLoad_z'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryBodyFixedTotalGLoadZ_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Time History: Total Passenger Z-Component G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*200 + 3458000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Total Passenger Z-Component G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-5 5])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Total Passenger Z-Component G-load (g)') % y-axis label
            set(gca,'YTick', -5:1:10);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[0 1],'k','LineWidth',2)
            
            field = {'passengerFrameTotalGLoad_z'};
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryPassengerFrameTotalGLoadZ_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Time History: Aero G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3457000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Aero G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 5])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Aero G-load (g)') % y-axis label
            set(gca,'YTick', 0:.2:5);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.aerodynamicFrameAerodynamicGLoad,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    ylim([0 2])
                    set(gca,'YTick', 0:.2:2);
                else
                                field = {'aerodynamicFrameAerodynamicGLoad'};
                    plotSpecificFieldOfIndividual( compilation, p, k, ii, field )
                end
            end
            
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryAerodynamicGLoad_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            % close(fig_num);
        end
    end
end

%% Time History: Lift Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3485100 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Lift Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 15])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Lift Force Magnitude (10^3 kN)') % y-axis label
            set(gca,'YTick', 0:1:15);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.currentLiftForce/1e6);
            end
            %legend('x-dir','y-dir','z-dir')
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryCurrentLiftForce_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %   close(fig_num);
        end
    end
end

%% Time History: Drag Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3475100 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Drag Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 3])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Drag Force Magnitude (10^3 kN)') % y-axis label
            set(gca,'YTick', 0:1:3);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.currentDragForce/1e6);
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.currentDragForce(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/1e6,'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                end
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.currentDragForce(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/1e6,'s');
                end
            end
            %legend('x-dir','y-dir','z-dir')
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryCurrentDragForce_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %   close(fig_num);
        end
    end
end

%% Time History: Thrust Acc. Components - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3465000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust Components through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([0 5])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Thrust Acc. Components (m/s^2)') % y-axis label
            %set(gca,'YTick', 0:1:10);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            
            for ii = 1:numel(compilation(p).evolutions(k).trajectories)
                if ( compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.distance_to_go(end) < 30 )
                    % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_x,'k');
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_y,'r');
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_z,'b');
                end
            end
            legend('x-dir','y-dir','z-dir')
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryThrustAccelerationComponents_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

%% Time History: Thrust Acc. Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3465100 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust Acc. Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([0 5])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Thrust Acc. Magnitude (m/s^2)') % y-axis label
            %set(gca,'YTick', 0:1:10);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_M);
            end
            %legend('x-dir','y-dir','z-dir')
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryThrustAccelerationMagnitude_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %   close(fig_num);
        end
    end
end

%% Time History: Thrust Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3465100 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 15])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Thrust Force Magnitude (10^3 kN)') % y-axis label
            set(gca,'YTick', 0:1:15);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.thrustMagnitude/1e6);
            end
            %legend('x-dir','y-dir','z-dir')
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryThrustMagnitude_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %   close(fig_num);
        end
    end
end

%% Time History: Thrust g-Load - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3465100 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust g-Load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Thrust g-Load (g)') % y-axis label
            set(gca,'YTick', 0:1:10);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            for ii = 1:numel(compilation(p).evolutions(k).trajectories)
                %   plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                %      compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.totalThrustGLoad);
            end
            %legend('x-dir','y-dir','z-dir')
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryThrustGLoad_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %   close(fig_num);
        end
    end
end

%% Time History: Gravity Force Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 3466300 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Gravity Force Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Gravity Force Magnitude (N)') % y-axis label
            set(gca,'YTick', 0:1:10);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            for ii = 1:numel(compilation(p).evolutions(k).trajectories)
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    ((compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.localGravity_1).^2 + (compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.localGravity_2).^2).^(1/2));
                
            end
            %legend('x-dir','y-dir','z-dir')
            
            %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryGravityMagnitude_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %   close(fig_num);
        end
    end
end

%% Time History: Dynamic Pressure - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 654180 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Dynamic Pressure through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 35])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Dynamic Pressure (kPa)') % y-axis label
            set(gca,'YTick', 0:5:35);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],25*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicPressure/1e3,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    %ylim([0 1])
                    %set(gca,'YTick', 0:.2:2);
                else
                    h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicPressure/1e3);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicPressure(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/1e3,'x');
                        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    end
                    
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicPressure(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/1e3,'s');
                    end
                end
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryDynamicPressure_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %      close(fig_num);
        end
    end
end

%% Time History: Bending Moment (Q-alpha) - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 674180 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Bending Moment (Q-alpha) through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 6])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Bending Moment, aka (Q-alpha) (kPa-rad)') % y-axis label
            set(gca,'YTick', 0:1:6);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(4363.3231299858)*[1 1]/1e3,'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bendingMoment/1e3);
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bendingMoment(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/1e3,'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                end
                
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bendingMoment(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/1e3,'s');
                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryBendingMoment_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %      close(fig_num);
        end
    end
end

end