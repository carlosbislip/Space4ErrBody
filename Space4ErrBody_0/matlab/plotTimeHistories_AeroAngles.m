function [  ] = plotTimeHistories_AeroAngles( compilation )


%% Time History: Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 1000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Angle of Attack  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 50])
            max_tof = max([compilation(p).evolutions.max_tof]);
            %max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Commanded Angle of Attack (deg)') % y-axis label
            set(gca,'YTick', 0:5:50);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            if  compilation(1).validation == 1
                
                for ii = compilation(p).evolutions(k).population(1).indices.printed
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    ylim([10 45])
                    set(gca,'YTick', 0:5:50);
                end
            else
                ii = compilation(p).evolutions(k).population(1).indices.printed;
                
                plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'angleOfAttack'}, 1 )
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryAngleOfAttack_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

%% Time History: Evaluated Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 2000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Angle of Attack through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([0 50])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:5:50);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluatedAngleOfAttack);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/timeHistoryEvaluatedAngleOfAttack_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Commanded Angle of Attack - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Angle of Attack  through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([0 50])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', 0:5:50);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            if  compilation(1).validation == 1
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commandedAngleOfAttack,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([10 45])
                set(gca,'YTick', 0:5:50);
                
            else
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commandedAngleOfAttack);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            '/figures/timeHistoryCommandedAngleOfAttack_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
                if compilation(p).evolutions(k).population(1).indices.printed > 0
        
        fig_num = p*100 + 4000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        % max_tof = 1400;
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        if  compilation(1).validation == 1
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bankAngle,'k','LineWidth',2);
                plot([0 1400],[0 0],'k','LineWidth',2)
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
                ylim([-80 100])
                set(gca,'YTick', -80:20:100);
                
            end
            
        else
            ii = compilation(p).evolutions(k).population(1).indices.printed;
            
            plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'bankAngle'}, 1 )
        end
        
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryBankAngle_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        %  close(fig_num);
    end
    end
end

%% Time History: Skip Suppresion Limit
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 5000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Skip Suppression Limit through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(90*[-1 1])
            max_tof = max([compilation(p).evolutions.max_tof]);
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Skip Suppression Limit (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            plot([0 max_tof],(0)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    abs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.skipSuppressionBankAngleLimit));
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.skipSuppressionBankAngleLimit(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.skipSuppressionBankAngleLimit(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistorySkipSuppressionBankAngleLimit_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
        
    end
end



%% Time History: Evaluated Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        %for k = 1
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 6000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Evaluated Bank Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Evaluated Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluatedBankAngle);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            '/figures/timeHistoryEvaluatedBankAngle_Evolution',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Commanded Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 7000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Commanded Bank Angle through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([-90 90])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        grid on
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commandedBankAngle);
        end
        
        plot([0 max_tof],[0 0],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'timeHistoryCommandedBankAngle_Evolution',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
       % close(fig_num);
    end
end

%% Time History: Bank Reversal Trigger - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + 8000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Reversal Trigger - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        ylim([0 1.5])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Commanded Bank Reversal Trigger (-)') % y-axis label
        set(gca,'YTick', 0:1:2);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        % for ii = numel(compilation(p).evolutions(k).trajectories):numel(compilation(p).evolutions(k).trajectories)
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger);
        end
        
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            '/figures/timeHistoryBankAngleReversalTrigger_Evolution',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Time History: Flight-Path Angle - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 9000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Flight-Path Angle through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-15 10])
            max_tof = max([compilation(p).evolutions.max_tof]);
            %max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Flight-Path Angle (deg)') % y-axis label
            %set(gca,'YTick', -15:3:3);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],[0 0],'k','LineWidth',2)
            
            if  compilation(1).validation == 1
                
                for ii = compilation(p).evolutions(k).population(1).indices.printed
                    
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngle,'k','LineWidth',2);
                    plot([0 1400],[0 0],'k','LineWidth',2)
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    ylim([-10 0])
                    set(gca,'YTick', -10:1:0);
                end                
            else
                ii = compilation(p).evolutions(k).population(1).indices.printed;
                
                plotFieldVsField( compilation, p, k, ii, {'timeOfFlight'}, 1, {'flightPathAngle'}, 1 )
            end
            
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryFlightPathAngle_Generation',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Time History: Flight-Path Angle Rate - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 10000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Flight-Path Angle Rate through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-1 1])
            max_tof = max([compilation(p).evolutions.max_tof]);
            %max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Flight-Path Angle Rate (deg/s)') % y-axis label
            set(gca,'YTick', -1:.25:1);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngleRate,'k','LineWidth',2);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngleRate,'b','LineWidth',2);
                    plot([0 1400],[0 0],'k','LineWidth',2)
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                    ylim([-2 2])
                    set(gca,'YTick', -2:1:2);
                    
                else
                    h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngleRate);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngleRate(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    end
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngleRate(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
                    end
                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryFlightPathAngleRate_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            % close(fig_num);
        end
    end
end



%% Time History: BodyFlap Deflection Angle - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 11000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('BodyFlap Deflection Angle through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-20 30])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof =1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('BodyFlap Deflection Angle (deg)') % y-axis label
            set(gca,'YTick', -20:5:30);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyflapDeflectionAngle,'k','LineWidth',2);
                    xlim([0 1400])
                    plot([0 1400],[0 0],'k','LineWidth',2)
                    set(gca,'XTick', 0:200:1400);
                    ylim([-20 20])
                    set(gca,'YTick', -20:5:20);
                else
                    h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyflapDeflectionAngle);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyflapDeflectionAngle(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    end
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyflapDeflectionAngle(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
                    end
                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryBodyFlapDeflectionAngle_Generation',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

%% Time History: Elevon Deflection Angle - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            fig_num = p*100 + 12000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Elevon Deflection Angle through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-40 40])
            max_tof = max([compilation(p).evolutions.max_tof]);
            %max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Elevon Deflection Angle (deg)') % y-axis label
            set(gca,'YTick', -40:10:40);
            set(gca,'XTick', 0:200:max_tof);
            
            hold on
            grid on
            
            plot([0 max_tof],[0 0],'k','LineWidth',2)
            
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.elevonDeflectionAngle,'k','LineWidth',2);
                    xlim([0 1400])
                    plot([0 1400],[0 0],'k','LineWidth',2)
                    set(gca,'XTick', 0:200:1400);
                    ylim([-40 40])
                    set(gca,'YTick', -40:5:40);
                else
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.elevonDeflectionAngle);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.elevonDeflectionAngle(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.elevonDeflectionAngle(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'timeHistoryElevonDeflectionAngle_Generation',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

end