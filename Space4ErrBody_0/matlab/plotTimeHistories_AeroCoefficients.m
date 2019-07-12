function [  ] = plotTimeHistories_AeroCoefficients( compilation )

%% Time History: Pitch Moment Coefficient - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 725000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Cm through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            ylim([-.01 .01])
            max_tof = max([compilation(p).evolutions.max_tof]);
            %max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Time (s)') % x-axis label
            ylabel('C_m (-)') % y-axis label
            set(gca,'YTick', -.01:0.001:.01);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],(0)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if  compilation(1).validation == 1
                    stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.aerodynamicCoefficient_Cm,'k','LineWidth',2);
                else
                    stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.aerodynamicCoefficient_Cm);
                end
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryMomentCoefficientPitch_Evolution',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).set),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end

%% Time History: Cm Increment due to Control Surfaces - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 724000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Cm Increment due to Control Surfaces through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
            ylim([-.1 .03])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Time (s)') % x-axis label
            ylabel('\Delta C_m (-)') % y-axis label
            set(gca,'YTick', -.1:0.01:.03);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 8000],(0)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if  compilation(1).validation == 1
                    stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.increment_Cm_bodyflap_dif,'k','LineWidth',2);
                else
                    h = stairs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.increment_Cm_bodyflap_dif);
                      set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.increment_Cm_bodyflap_dif(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.increment_Cm_bodyflap_dif(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
     
                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/timeHistoryMomentCoefficientIncrementPitch_Evolution',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).set),...
                '.png'),...
                'png');
            %        close(fig_num);
        end
    end
end






end