function [  ] = plotCompoundRelations( compilation )


for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 100000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height vs. Freestream Temperature - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 130])
            xlim([150 400])
            xlabel('Temperature (K)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:10:130);
            set(gca,'XTick', 150:25:400);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.intersectionCase,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                figurePath,...
                'freestreamTemperature_v_height_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Time History: Height & Commanded Throttle Setting & Total g-Load - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + 823000 + k*1;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height through Time - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            
            ylim([0 1])
            max_tof = max([compilation(p).evolutions.max_tof]);
            % max_tof = 1400;
            xlim([0 max_tof])
            xlabel('Propagation Time (s)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:.1:1);
            set(gca,'XTick', 0:200:max_tof);
            hold on
            grid on
            
            %plot([0 max_tof],(25)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/max(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height));
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commandedThrottleSetting);
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude/max(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude));
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack/max(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack));
                
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                figurePath,...
                'height_v_T_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Height vs. Freestream Temperature - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 100000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height vs. Freestream Temperature - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 130])
            xlim([150 400])
            xlabel('Temperature (K)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:10:130);
            set(gca,'XTick', 150:25:400);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.freestreamTemperature,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                figurePath,...
                'freestreamTemperature_v_height_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Angle of Attack vs. Mach Number - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 100000;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            %title(strcat('Angle of Attack vs. Mach Number - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 50])
            xlim([0 30])
            xlabel('Mach Number $(-)$','Interpreter','latex') % x-axis label
            ylabel('Angle of Attack $($deg$)$','Interpreter','latex') % y-axis label
            set(gca,'YTick', 0:10:50);
            set(gca,'XTick', 0:2.5:30);
            set(gca,'TickLabelInterpreter','latex')
            
            hold on
            grid on
            
            plot(compilation(p).evolutions(k).Common.Bounds.AngleOfAttack.machNumber,compilation(p).evolutions(k).Common.Bounds.AngleOfAttack.upperBound,'k','LineWidth',2);
            plot(compilation(p).evolutions(k).Common.Bounds.AngleOfAttack.machNumber,compilation(p).evolutions(k).Common.Bounds.AngleOfAttack.lowerBound,'k','LineWidth',2);
            if  compilation(1).validation == 1
                
                for ii = compilation(p).evolutions(k).population(1).indices.printed
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.machNumber,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack,'k','LineWidth',2);
                end
            else
                %ii = compilation(p).evolutions(k).population(1).indices.printed;
                
                id = [1544 86] + 1;
                p = 1;
                ii = id(p);
                
                ax = gca;
                ax.ColorOrderIndex = 1;
                % plotFieldVsField( compilation, p, k, ii, 'machNumber', 1, 'angleOfAttack', 1, 2 )
                p = 2;
                ii = id(p);
                ax = gca;
                ax.ColorOrderIndex = 2;
                % plotFieldVsField( compilation, p, k, ii, 'machNumber', 1, 'angleOfAttack', 1, 2 )
                
            end
            
            
            field_x = 'machNumber';
            field_x_struct = compilation(end).evolutions(end).population(end).dependentVariableTimeHistory.(field_x);
            variableLabel_x = field_x_struct.variableLabel;
            figureSaveName_x = field_x_struct.figureSaveNameContent;
            units_x = field_x_struct.units;
            limits_x12 = field_x_struct.limits;
            limits_x3 = limits_x12;
            scalingFactor_x = field_x_struct.scalingFactor;
            tick_x = field_x_struct.tick;
            
            
            
            field_y = 'angleOfAttack';
            field_y_struct = compilation(end).evolutions(end).population(end).dependentVariableTimeHistory.(field_y);
            variableLabel_y = field_y_struct.variableLabel;
            figureSaveName_y = field_y_struct.figureSaveNameContent;
            units_y = field_y_struct.units;
            limits_y = field_y_struct.limits;
            scalingFactor_y = field_y_struct.scalingFactor;
            tick_y = field_y_struct.tick;
            
            
            
            p = 1;
            selectIndividual = id(p);
            
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1)) == false
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_y;
            else
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2)) == false
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_y;
            else
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
            y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
            
            j = 1;
            ax = gca;
            ax.ColorOrderIndex = 1;
            h(1) = plot(x,y,'LineWidth',2);
            
            scatter(phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            
            ax = gca;
            ax.ColorOrderIndex = 2;
           % j = j + 1;
            
            %ascentRange = 1:compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2);
            %descentRange = compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2):numel(x);
            
           % h(1) = plot(x,y,'LineWidth',2);
            
            p = 2;
            selectIndividual = id(p);
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1)) == false
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_y;
            else
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2)) == false
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_y;
            else
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
            y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
            
            %ascentRange = 1:compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2);
            %descentRange = compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2):numel(x);
           
            ax = gca;
            ax.ColorOrderIndex = 2;
            h(2) = plot(x,y,'LineStyle','--','LineWidth',2);
           
            h(3) = scatter( phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            
            legend(h(1:3),[{'Initial Individual'};{'Re-evaluated Individual'};{'Phase Change'}],'Interpreter','latex','Location','southeast')
            
            
            hold off
            saveas(...
                figure(fig_num),...
                string(strcat(...
                compilation(p).figurePath,...
                figureSaveName_y, {'_vs_'},figureSaveName_x,'_generation_',...
                num2str(k - 1),'_case',convertCharsToStrings(compilation(p).case),'.png')),...
                'png');
            
            %     close(fig_num);
            
        end
    end
end



%% Total Body G-load vs. Central Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 101000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Total Body G-load vs. Central Target Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 10])
            xlim([0 60])
            xlabel('Central Target Angular Distance To Go (deg)') % x-axis label
            ylabel('Total Body G-load (g_0)') % y-axis label
            set(gca,'YTick', 0:2:10);
            set(gca,'XTick', 0:10:60);
            
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            ii = compilation(p).evolutions(k).population(1).indices.printed;
            plotFieldVsField( compilation, p, k, ii, 'centralTargetAngularDistanceToGo', 1, 'bodyFrameTotalGLoadMagnitude', 1 )
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bodyFrameTotalGLoadMagnitude_v_centralTargetAngularDistanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

%% Total Body G-load vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 111000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Total Body G-load vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 30])
            %xlim([0 60])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Total Body G-load (g_0)') % y-axis label
            set(gca,'YTick', 0:5:30);
            %set(gca,'XTick', 0:10:60);
            
            hold on
            grid on
            
            plot([0 3],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude);
                %  plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluated_bank_angle);
                % plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle);
                %
                %                 plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                %                     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commanded_bank_angle);
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bodyFrameTotalGLoadMagnitude_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Skip Suppression Limit vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 101000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Skip Suppression Limit vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(rad2deg([-pi/2 pi/2]))
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Skip Suppression Limit (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 15
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.skipSuppressionBankAngleLimit);
                    %  plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluated_bank_angle);
                    % plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle);
                    %
                    %                 plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    %                     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commanded_bank_angle);
                    
                end
            end
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'skip_suppression_limit_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Skip Suppression Limit vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 101000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Skip Suppression Limit vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(rad2deg([-pi/2 pi/2]))
            %xlim([0 1])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Skip Suppression Limit (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            %set(gca,'XTick', 0:.1:1);
            
            hold on
            grid on
            
            plot([0 3],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.skipSuppressionBankAngleLimit);
                %  plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluated_bank_angle);
                % plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle);
                %
                %                 plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                %                     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commanded_bank_angle);
                
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'skip_suppression_limit_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %    close(fig_num);
        end
    end
end

%% Total Passenger Z-Component G-load vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 102000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Total Passenger Z-Component G-load vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-5 5])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Total Passenger Z-Component G-load (g)') % y-axis label
            set(gca,'YTick', -5:1:5);
            set(gca,'XTick', 0:10:60);
            hold on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                %   if sum(abs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heading_error) > 30) == 0
                %       if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 30
                %          if sum(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger) > 2
                %             if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy(end) < 0.3
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.passenger_fixed_total_g_load_z);
                %     end
                % end
                %    end
                %end
            end
            
            plot([0 60],(0)*[1 1],'k','LineWidth',2)
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'total_passenger_g_load_z_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Specific Energy vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 103000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Specific Energy vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([0 1.5])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Normalized Specific Energy (-)') % y-axis label
            %set(gca,'YTick', 0:.1:1.5);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            plot([0 60],(1)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.specificEnergy);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'specificEnergy_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Normalized Specific Energy vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 103000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Normalized Specific Energy vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 1.5])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Normalized Specific Energy (-)') % y-axis label
            set(gca,'YTick', 0:.1:1.5);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            plot([0 60],(1)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/normalizedSpecificEnergy_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Thrust Magnitude vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 104000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust Magnitude vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([0 1100])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Thrust Magnitude (kN)') % y-axis label
            %set(gca,'YTick', 0:100:1100);
            set(gca,'XTick', 0:10:60);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.thrustMagnitude/1e3);
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'thrustMagnitude_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Thrust Magnitude vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 104000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust Magnitude vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([0 1100])
            %xlim([0 1])
            xlabel('Normalized Specific Energy ()') % x-axis label
            ylabel('Thrust Magnitude (kN)') % y-axis label
            %set(gca,'YTick', 0:100:1100);
            set(gca,'XTick', 0:.1:1);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.thrustMagnitude/1e3);
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'thrustMagnitude_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Thrust Acceleration Components vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 105000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust Acceleration Components vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-40 40])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Thrust Acceleration Components (m/s^2)') % y-axis label
            set(gca,'YTick', -40:10:40);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_x,'k');
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_y,'r');
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_z,'b');
                
            end
            
            legend('x-dir','y-dir','z-dir')
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'thrustAccelerationsComponents_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Thrust Acceleration Components vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 106000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Thrust Acceleration Magnitude vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 40])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Thrust Acceleration Magnitude (m/s^2)') % y-axis label
            set(gca,'YTick', 0:5:40);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.acc_thru_M);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'thrustAccelerationMagnitude_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Height vs. Density - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 107000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height vs. Density - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 100])
            xlim([0 .05])
            xlabel('Density (kg/m^3)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:10:100);
            set(gca,'XTick', 0:.01:.05);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.localDensity,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'height_v_localDensity_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Density vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 108000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Density vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            % ylim([0 800])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Density (kg/m^3)') % y-axis label
            % set(gca,'YTick', 0:100:800);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                %   if ( compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 20 )
                %     if( sum(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger) < 6 )
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.localDensity);
                % end
                %     end
            end
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'density_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Height vs. Dynamic Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 109000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height vs. Dynamic Target Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 200])
            xlim([0 60])
            xlabel('Angular Distance To Go (deg)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:20:200);
            set(gca,'XTick', 0:10:60);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
                
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/1e3,'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                end
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/1e3,'s');
                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'height_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

%% Mass vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 110000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Mass vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([0 150])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Mass (10^3 kg)') % y-axis label
            %set(gca,'YTick', 0:15:150);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass/1e3);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'mass_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Mass vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 110000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Mass vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            % ylim([0 150])
            xlim([0 1])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Mass (10^3 kg)') % y-axis label
            %set(gca,'YTick', 0:15:150);
            set(gca,'XTick', 0:.1:1);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass/1e3);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'mass_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Mass Rate vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 111000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Mass Rate vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 250])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Mass Rate (kg/s)') % y-axis label
            set(gca,'YTick', 0:25:250);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    abs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.massRate) );
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'massRate_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Mass vs. Height - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 112000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Mass vs. Height - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            % ylim([0 150])
            % max_tof = max([compilation(p).evolutions.max_tof]);
            % xlim([0 8000])
            xlabel('Height (km)') % x-axis label
            ylabel('Mass (kg)') % y-axis label
            % set(gca,'YTick', 0:10:150);
            % set(gca,'XTick', 0:500:8000);
            hold on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'mass_v_height_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Height vs. Evaluated Throttle - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 113000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Evaluated Throttle Setting vs. Height  - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            % ylim([0 150])
            % xlim([0 8000])
            ylabel('Height (km)') % x-axis label
            xlabel('Evaluated Throttle Setting (-)') % y-axis label
            % set(gca,'YTick', 0:10:150);
            % set(gca,'XTick', 0:500:8000);
            hold on
            grid on
            %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluatedThrottleSetting,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'evaluatedThrottle_v_height_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Evaluated Throttle Setting vs. Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 114000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Evaluated Throttle Setting vs. Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([180 360])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Evaluated Throttle Setting (-)') % y-axis label
            %set(gca,'YTick', 180:30:360);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluatedThrottleSetting);
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'evaluatedThrottleSetting_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Commanded Throttle Setting vs. Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 115000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Commanded Throttle Setting vs. Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([180 360])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Commanded Throttle Setting (-)') % y-axis label
            %set(gca,'YTick', 180:30:360);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commandedThrottleSetting);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'commandedThrottleSetting_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Commanded Throttle Setting vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 115000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Commanded Throttle Setting vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            %ylim([180 360])
            xlim([0 1])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Commanded Throttle Setting (-)') % y-axis label
            %set(gca,'YTick', 180:30:360);
            set(gca,'XTick', 0:.1:1);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if ( compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 30 )
                    %    if( sum(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger) < 6 )
                    
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.commandedThrottleSetting);
                    %     end
                end
                
                
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'commandedThrottleSetting_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Height vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 118000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 130])
            % xlim([0 1])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:10:130);
            %set(gca,'XTick', 0:.1:1);
            hold on
            grid on
            
            %plot([0 3],(25)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
            end
            
            %plot(V,compilation(1).evolutions(1).trajectories(1).individual.height/1e3)
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'height_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Airspeed vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 116000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Airspeed vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 8000])
            xlim([0 1])
            ylabel('Airspeed (m/s)') % x-axis label
            xlabel('Normalized Specific Energy (-)') % y-axis label
            set(gca,'YTick', 0:500:8000);
            set(gca,'XTick', 0:.1:1);
            hold on
            grid on
            
            %  plot([0 1],(25)*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed);
            end
            
            %plot(V,compilation(1).evolutions(1).trajectories(1).individual.height/1e3)
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'airspeed_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Height vs. Airspeed - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 116000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height vs. Airspeed - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 130])
            xlim([0 8000])
            xlabel('Airspeed (m/s)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:10:130);
            set(gca,'XTick', 0:500:8000);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plotrange = 1:compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2);
                if  compilation(1).validation == 1
                    a = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(plotrange)/1e3,'k','LineWidth',2);
                    b = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_SkipSuppression(plotrange)/1e3);
                    c = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_DynamicPressure(plotrange)/1e3);
                    d = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_MechanicalLoad_descent(plotrange)/1e3);
                    e = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_ThermalLoad(plotrange)/1e3);
                    f = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_BendingMoment(plotrange)/1e3);
                    legend([a b c d e f],'Height','Equilibrium Glide','Dynamic Pressure','Mechanical Load','Thermal Load','Bending Moment','Location','northwest')
                    
                else
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(plotrange)/1e3,'k','LineWidth',2);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_SkipSuppression(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_DynamicPressure(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_MechanicalLoad_descent(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_MechanicalLoad_ascent(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_ThermalLoad(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_BendingMoment(plotrange)/1e3);
                    
                    legend('Height','Skip Suppression','Dynamic Pressure','Mechanical Load - Descent','Mechanical Load - Ascent','Thermal Load','Bending Moment')
                end
                
                clear plotrange;
            end
            %plot([0 8000],(10)*[1 1],'k','LineWidth',2)
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'height_v_V_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Flight Corridor: Descent - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 116000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Flight Corridor: Descent - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 130])
            xlim([0 8000])
            xlabel('Airspeed (m/s)') % x-axis label
            ylabel('Height (km)') % y-axis label
            set(gca,'YTick', 0:10:130);
            set(gca,'XTick', 0:500:8000);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plotrange = (compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2) + 1):length(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed);
                if  compilation(1).validation == 1
                    a = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(plotrange)/1e3,'k','LineWidth',2);
                    b = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_SkipSuppression(plotrange)/1e3);
                    c = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_DynamicPressure(plotrange)/1e3);
                    d = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_MechanicalLoad_descent(plotrange)/1e3);
                    e = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_ThermalLoad(plotrange)/1e3);
                    f = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_BendingMoment(plotrange)/1e3);
                    legend([a b c d e f],'Height','Equilibrium Glide','Dynamic Pressure','Mechanical Load','Thermal Load','Bending Moment','Location','northwest')
                    
                else
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height(plotrange)/1e3,'k','LineWidth',2);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_SkipSuppression(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_DynamicPressure(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_MechanicalLoad_descent(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_MechanicalLoad_ascent(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_ThermalLoad(plotrange)/1e3);
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed(plotrange),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightCorridorBoundary_BendingMoment(plotrange)/1e3);
                    
                    legend('Height','Skip Suppression','Dynamic Pressure','Mechanical Load - Descent','Mechanical Load - Ascent','Thermal Load','Bending Moment')
                end
                
                clear plotrange;
            end
            %plot([0 8000],(10)*[1 1],'k','LineWidth',2)
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'flightCorridor_Descent_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Heading Angle vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 117000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Heading Angle vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([180 360])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Heading Angle (deg)') % y-axis label
            set(gca,'YTick', 180:30:360);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                %if ( compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 20 )
                %   if( sum(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger) < 6 )
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngle);
                %  end
                % end
                
                
            end
            
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'headingAngle_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Heading Error vs. Dynamic Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 118000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Heading Error vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(30*[-1 1])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Heading Error (deg)') % y-axis label
            set(gca,'YTick', -30:10:30);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot(compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo,compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.lowerBound,'k','LineWidth',2);
            plot(compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo,compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.upperBound,'k','LineWidth',2);
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingError,'k','LineWidth',2);
                    
                else
                    h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo.value,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToDynamicTargetError.value);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo.value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToDynamicTargetError.value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    end
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo.value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToDynamicTargetError.value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
                    end
                    
                end
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'dynamicTargetHeadingError_v_angularDistanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Heading Error to Central Target vs. Central Target Angular Low Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 119000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Heading Error vs. Low Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(30*[-1 1])
            xlim([.7 1.5])
            xlabel('Low Distance To Go (deg)') % x-axis label
            ylabel('Heading Error (deg)') % y-axis label
            set(gca,'YTick', -30:10:30);
            set(gca,'XTick', .7:.1:1.5);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot(compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo,compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.lowerBound,'k','LineWidth',2);
            plot(compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo,compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.upperBound,'k','LineWidth',2);
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToCentralTargetError,'k','LineWidth',2);
                    
                else
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToCentralTargetError);
                end
            end
            
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).FigurePath,...
                'centralTargetHeadingError_v_centralTargetAngularLowDistanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Heading Error to Central Target vs. Central Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 118000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Heading Error vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(30*[-1 1])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Heading Error (deg)') % y-axis label
            set(gca,'YTick', -30:10:30);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot(compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo,compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.lowerBound,'k','LineWidth',2);
            plot(compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.angularDistanceToGo,compilation(p).evolutions(k).Common.Bounds.headingErrorDeadBand.upperBound,'k','LineWidth',2);
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                if  compilation(1).validation == 1
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToCentralTargetError,'k','LineWidth',2);
                    
                else
                    h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToCentralTargetError);
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToCentralTargetError(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
                        set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                    end
                    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                        scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.headingAngleToCentralTargetError(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
                    end
                    
                end
                
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'centralTargetHeadingError_v_angularDistanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Angle of Attack vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 119000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Angle of Attack vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 50])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Angle of Attack (deg)') % y-axis label
            set(gca,'YTick', 0:10:50);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                %if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 15
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack);
                %end
            end%
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'angleOfAttack_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Angle of Attack vs. Distance Covered Ratio - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 119000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Angle of Attack vs. Distance Covered Ratio - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 50])
            xlim([0 1])
            xlabel('Distance Covered Ratio (-)') % x-axis label
            ylabel('Angle of Attack (deg)') % y-axis label
            set(gca,'YTick', 0:10:50);
            set(gca,'XTick', 0:.1:1);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                
                distanceToCover = compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(1);
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceTraveled/distanceToCover,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'angleOfAttack_v_distanceCoveredRatio_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Angle of Attack vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 119000;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            %title(strcat('Angle of Attack vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 50])
            
            xlabel('Normalized Specific Energy $(-)$','Interpreter','latex') % x-axis label
            ylabel('Angle of Attack $($deg$)$','Interpreter','latex') % y-axis label
            set(gca,'YTick', 0:10:50);
            % set(gca, 'XDir','reverse')
            set(gca,'TickLabelInterpreter','latex')
            
            hold on
            grid on
            
            %for ii = compilation(p).evolutions(k).population(1).indices.printed
            id = [1544 86] + 1;
            p = 1;
            % ii = id(p);
            %                 h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy.value,...
            %                     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack.value,'LineWidth',2);
            %
            %                        p = 2;
            % ii = id(p);
            %                      h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy.value,...
            %                     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angleOfAttack.value,'LineWidth',2);
            %
            
            field_x = 'normalizedSpecificEnergy';
            field_x_struct = compilation(end).evolutions(end).population(end).dependentVariableTimeHistory.(field_x);
            variableLabel_x = field_x_struct.variableLabel;
            figureSaveName_x = field_x_struct.figureSaveNameContent;
            units_x = field_x_struct.units;
            limits_x12 = field_x_struct.limits;
            limits_x3 = limits_x12;
            scalingFactor_x = field_x_struct.scalingFactor;
            tick_x = field_x_struct.tick;
            
            
            
            field_y = 'angleOfAttack';
            field_y_struct = compilation(end).evolutions(end).population(end).dependentVariableTimeHistory.(field_y);
            variableLabel_y = field_y_struct.variableLabel;
            figureSaveName_y = field_y_struct.figureSaveNameContent;
            units_y = field_y_struct.units;
            limits_y = field_y_struct.limits;
            scalingFactor_y = field_y_struct.scalingFactor;
            tick_y = field_y_struct.tick;
            
            
            
            p = 1;
            selectIndividual = id(p);
            
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1)) == false
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_y;
            else
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2)) == false
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_y;
            else
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            x1 = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
            y1 = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
            
            j = 1;
            ax = gca;
            ax.ColorOrderIndex = j;
            h(1) = plot(x1,y1,'LineWidth',2);
            
            %scatter( phaseChange_x(1),phaseChange_y(1),70,'kx')
            scatter( phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            
            j = j + 1;
            ax = gca;
            ax.ColorOrderIndex = j;
            
            
            %ascentRange = 1:compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2);
            %descentRange = compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2):numel(x);
            
            p = 2;
            selectIndividual = id(p);
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1)) == false
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_y;
            else
                phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2)) == false
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_y;
            else
                phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
                phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
            end
            
            x2 = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
            y2 = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
            
            %ascentRange = 1:compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2);
            %descentRange = compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2):numel(x);
            
            h(2) = plot(x2,y2,'LineStyle','--','LineWidth',2);
            
            
            ax = gca;
            %h(3) = scatter( phaseChange_x(1),phaseChange_y(1),70,'kx');
            
            h(3) = scatter( phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            j = j + 1;
            
            
            maxValue = .25*ceil(max([x1;x2])/.25);
            
            limits_x = [0 maxValue];
            xlim(limits_x)
            set(gca,'XTick', limits_x(1):.25:limits_x(2));
            
            legend(h(:),[{'Initial Individual'};{'Re-evaluated Individual'};{'Phase Chase'}],'Interpreter','latex','Location','northeast')
            
            
            hold off
            
            
            saveas(...
                figure(fig_num),...
                string(strcat(...
                compilation(p).figurePath,...
                figureSaveName_y, {'_vs_'},figureSaveName_x,'_generation_',...
                num2str(k - 1),'_case',convertCharsToStrings(compilation(p).case),'.png')),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Bank Angle vs. Central Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 220000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Bank Angle vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(90*[-1 1])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Bank Angle (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.centralTargetAngularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bankAngle);
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bankAngle_v_centralTargetAngularDistanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Bank Angle vs. Dynamic Target Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 120000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Bank Angle vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(90*[-1 1])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Bank Angle (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicTargetAngularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bankAngle);
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bankAngle_v_dynamicTargetAngularDistanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Bank Angle vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 120000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Bank Angle vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(90*[-1 1])
            % xlim([0 1])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Bank Angle (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            % set(gca,'XTick', 0:.1:1);
            % set(gca, 'XDir','reverse')
            hold on
            grid on
            
            plot([0 3],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                %if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 15
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bankAngle);
                % end
            end
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bankAngle_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Bank Angle vs. Airspeed - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 120000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Bank Angle vs. Airspeed - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(90*[-1 1])
            xlim([0 8000])
            xlabel('Airspeed (m/s)') % x-axis label
            ylabel('Bank Angle (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            % set(gca,'XTick', 0:.1:1);
            % set(gca, 'XDir','reverse')
            set(gca,'XTick', 0:500:8000);
            hold on
            grid on
            
            %plot([0 3],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                %if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 15
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bankAngle);
                % end
            end
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bankAngle_v_airspeed_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Evaluated Bank Angle vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 121000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Evaluated Bank Angle vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(90*[-1 1])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Evaluated Bank Angle (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 15
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluatedBankAngle);
                end
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'evaluatedBankAngle_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Evaluated Bank Angle vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 121000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Evaluated Bank Angle vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim(90*[-1 1])
            % xlim([0 3])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Evaluated Bank Angle (deg)') % y-axis label
            set(gca,'YTick', -90:15:90);
            % set(gca,'XTick', 0:.1:3);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot([0 3],[0 0],'k','LineWidth',2)
            
            %for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.evaluatedBankAngle);
            end
            
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'evaluatedBankAngle_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Bank Reversal Trigger vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 122000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Bank Reversal Trigger vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 2])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Bank Reversal Trigger (-)') % y-axis label
            set(gca,'YTick', 0:1:1);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger);
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bankReversalTrigger_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Flight-Path Angle vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 123000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Flight-Path Angle vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-15 10])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Flight-Path Angle (deg)') % y-axis label
            set(gca,'YTick', -15:3:10);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 15
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngle);
                end
            end
            hold off
            saveas(...
                figure(fig_num),...
                sstrcat(...
                compilation(p).figurePath,...
                'flightPathAngle_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Flight-Path Angle vs. Normalized Specific Energy - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 123000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Flight-Path Angle vs. Normalized Specific Energy - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-15 10])
            %xlim([0 60])
            xlabel('Normalized Specific Energy (-)') % x-axis label
            ylabel('Flight-Path Angle (deg)') % y-axis label
            set(gca,'YTick', -15:3:10);
            %set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            
            plot([0 3],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.normalizedSpecificEnergy,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngle);
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'flightPathAngle_v_normalizedSpecificEnergy_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Flight-Path Angle Rate vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 124000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Flight-Path Angle Rate vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-1 1])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Flight-Path Angle Rate (deg/s)') % y-axis label
            set(gca,'YTick', -1:.1:1);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 15
                    plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flightPathAngleRate);
                end
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'flightPathAngleRate_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Dynamic Pressure vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 125000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Dynamic Pressure vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 10000]/1e3)
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Dynamic Pressure (kPa)') % y-axis label
            % set(gca,'YTick', -90:15:90);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            hold on
            grid on
            plot([0 60],[100 100]/1e3,'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.dynamicPressure/1e3);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'dynamicPressure_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Bending Moment vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 126000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Bending Moment vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 10])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Bending Moment (kPa-rad)') % y-axis label
            % set(gca,'YTick', -90:15:90);
            set(gca,'XTick', 0:10:60);
            % set(gca, 'XDir','reverse')
            
            hold on
            grid on
            plot([0 60],[5014 5014]/1e3,'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bendingMoment/1e3);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bendingMoment_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Pitch Moment Coefficient vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 127000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Cm vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-.001 .001])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('C_m (-)') % y-axis label
            set(gca,'YTick', -.001:0.0001:.001);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.aerodynamicCoefficient_Cm );
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'Cm_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% BodyFlap Deflection Angle vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 128000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('BodyFlap Deflection Angle vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-30 30])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('BodyFlap Deflection Angle (deg)') % y-axis label
            set(gca,'YTick', -30:5:30);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyflapDeflectionAngle);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'bodyflap_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Elevon Deflection Angle vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 129000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Elevon Deflection Angle vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([-45 45])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Elevon Deflection Angle (deg)') % y-axis label
            set(gca,'YTick', -45:10:45);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            plot([0 60],[0 0],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.elevonDeflectionAngle);
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'elevon_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Height vs. Tauber Heat Flux at Leading Edge - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 131000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Height vs. Tauber Heat Flux at Leading Edge - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 150])
            % xlim([0 8000])
            xlabel('Tauber Heat Flux at Leading Edge (W/m^2)') % x-axis label
            ylabel('Height (km)') % y-axis label
            % set(gca,'YTick', 0:10:150);
            % set(gca,'XTick', 0:500:8000);
            
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heat_flux_tauber_leadingedge,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3);
            end
            
            plot([0 8000],(25)*[1 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'height_v_heatFluxTauber_leadingedge_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            % close(fig_num);
        end
    end
end

%% Tauber Heat Flux at Leading Edge vs. Airspeed - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 132000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Tauber Heat Flux at Leading Edge vs. Airspeed - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 1000])
            xlim([0 8000])
            xlabel('Airspeed (m/s)') % x-axis label
            ylabel('Tauber Heat Flux at Leading Edge (W/m^2)') % y-axis label
            set(gca,'YTick', 0:100:1000);
            set(gca,'XTick', 0:500:8000);
            hold on
            
            
            %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heat_flux_tauber_leadingedge/1e3);
            end
            
            plot([0 8000],(25)*[1 1],'k','LineWidth',2)
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'heatFluxTauber_leadingedge_v_V_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            % close(fig_num);
        end
    end
end

%% Tauber Heat Flux at Leading Edge vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 133000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Tauber Heat Rate at Leading Edge vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 1000])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Tauber Heat Rate at Leading Edge (kW/m^2)') % y-axis label
            set(gca,'YTick', 0:100:1000);
            set(gca,'XTick', 0:10:60);
            hold on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heat_flux_tauber_leadingedge/1e3 );
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'heatRateTauber_leadingedge_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Chapman Heat Flux at Nose vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 134000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Chapman Heat Flux at Nose vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 800])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Chapman Heat Flux at Nose (kW/m^2)') % y-axis label
            set(gca,'YTick', 0:50:800);
            set(gca,'XTick', 0:10:60);
            
            hold on
            grid on
            
            plot([0 60],530*[1 1],'k','LineWidth',2)
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heatFluxChapmanNose/1e3 );
                set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heatFluxChapmanNose(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/1e3,'x');
                    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
                end
                
                if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
                    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
                        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heatFluxChapmanNose(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/1e3,'s');
                end
            end
            
            hold off
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'heatFluxChapman_nose_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% Chapman Eq. Wall Temp. at Nose vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 135000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('Chapman Eq. Wall Temp. at Nose vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 2000])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('Chapman Eq. Wall Temp. at Nose (K)') % y-axis label
            set(gca,'YTick', 0:200:2000);
            set(gca,'XTick', 0:10:60);
            hold on
            grid on
            
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.wall_temperature_chapman );
            end
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'wallTempChapman_nose_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end

%% TUDAT Heat Rate at Nose vs. Angular Distance To Go - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).population(1).indices.printed > 0
            
            fig_num = p*100 + k*1 + 136000;
            figure(fig_num)
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set (gca,'Fontsize',15)
            title(strcat('TUDAT Heat Rate at Nose vs. Angular Distance To Go - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            ylim([0 1000])
            xlim([0 60])
            xlabel('Distance To Go (deg)') % x-axis label
            ylabel('TUDAT Heat Rate at Nose (kW/m^2)') % y-axis label
            set(gca,'YTick', 0:100:1000);
            set(gca,'XTick', 0:10:60);
            hold on
            
            %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            for ii = compilation(p).evolutions(k).population(1).indices.printed
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heat_rate_TUDAT_nose/1e3 );
                
                
                
                plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                    (compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heat_rate_TUDAT_nose - 0.5*(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.localDensity).*(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.airspeed).^3)/1e3 );
                
            end
            
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'heatRateTUDAT_nose_v_distanceToGo_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end
end


%% Height vs. Groundtrack - per Generation
lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;

if compilation(1).validation == 1
    lon_i_deg = -106.7;
    lon_i_rad = deg2rad(-22.37);
    lon_i_deg = -22.37;
    lat_f_deg = 5;
    lon_f_deg = -53;
    validation = 1;
end



th = 0:pi/50:2*pi;
xunit = .75 * cos(th) + lon_f_deg;
yunit = .75 * sin(th) + lat_f_deg;
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        if compilation(p).evolutions(k).printedPopulationSize > 0
            %for k = 1
            fig_num = p*100 + k*1 + 137000;
            figure(fig_num)
            hold on
            title(strcat('Trajectories - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
            set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            xlabel('\tau (deg)') % x-axis label
            ylabel('\delta (deg)') % y-axis label
            zlabel('Height (km)') % z-axis label
            img = imread('img.jpg');
            imagesc([-180 180], [-90 90], (flipud(img)));
            set (gca,'Fontsize',20)
            set(gca,'XTick', -90:15:30);
            set(gca,'YTick', 0:15:90);
            set(gca,'ZTick', 0:15:150);
            xlim([-90 30])
            ylim([0 90])
            zlim([0 150])
            
            
            if compilation(1).validation == 1
                xlim([-82 -42])
                ylim([-8 12])
                
            end
            for ii =1:numel(compilation(p).evolutions(k).trajectories)
                
                color_line3(...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.longitude_angle,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.latitude_angle,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight);
                plot(...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.longitude_angle,...
                    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.latitude_angle,'g')
                %                 scatter3(...
                %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
                %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
                %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
                %                 scatter3(...
                %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
                %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
                %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
                
            end
            
            plot(lon_f_deg,lat_f_deg,'MarkerSize',20)
            plot(lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2)
            plot(180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2)
            plot(xunit, yunit,'k','LineWidth',2);
            scatter(lon_f_deg,lat_f_deg,100,'r','x')
            axP = get(gca,'Position');
            %  legend(legendtext,'Location','southeastoutside')
            set(gca, 'Position', axP)
            view([13 49])
            grid on
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'height_v_Groundtrack_Generation_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).case)',...
                '.png'),...
                'png');
        end
    end
end

%% Height vs. Angular Distance To Go vs. Mass - per Generation
for p = 1:numel(compilation)
    
    
    % for k = 1:numel(compilation(p).evolutions)
    fig_num = p*100 + k*1 + 138000;
    figure(fig_num)
    hold on
    title(strcat('Height vs. Angular Distance To Go vs. Mass - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
    xlabel('Distance To Go (deg)') % x-axis label
    ylabel('Mass (kg))') % y-axis label
    zlabel('Height (km)') % z-axis label
    % img = imread('img.jpg');
    %imagesc([-180 180], [-90 90], (flipud(img)));
    set (gca,'Fontsize',20)
    %set(gca,'XTick', -90:15:30);
    %set(gca,'YTick', 0:15:90);
    %set(gca,'ZTick', 0:15:150);
    %xlim([-90 30])
    %ylim([0 90])
    zlim([0 175])
    
    
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1        %     for k = 1:numel(compilation(p).evolutions)
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            if sum(abs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heading_error) > 30) == 0
                if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 30
                    if sum(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger) > 2
                        
                        %color_line3(...
                        %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.longitude_angle,...
                        %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.latitude_angle,...
                        %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,...
                        %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass);
                        plot3(...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass,...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3)
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
                    end
                end
            end
        end
    end
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([13 49])
    grid on
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).mainpath,...
        '/figures/Height_v_distanceToGo_v_mass_Set',...
        convertCharsToStrings(compilation(p).case)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end



%% Height vs. Angular Distance To Go vs. Bank Angle - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = 1
        fig_num = p*100 + k*1 + 139000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Angular Distance To Go vs. Bank Angle - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Bank Angle (deg))') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        set(gca,'XTick', 0:10:60);
        % set(gca,'YTick', 0:15:90);
        set(gca,'ZTick', 0:25:175);
        xlim([0 60])
        %ylim([0 90])
        zlim([0 175])
        hold on
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            if sum(abs(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.heading_error) > 30) == 0
                if compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo(end) < 30
                    if sum(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle_reversal_trigger) > 2
                        
                        plot3(...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle,...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3)
                        %plot3(...
                        %   compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        %  compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.skip_suppression_limit,...
                        % compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3)
                        
                        % plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bank_angle,'g');
                        % plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.skip_suppression_limit,'k');
                        
                        % plot3(...
                        %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                        %     ones(numel(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo))*100,...
                        %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,'g')
                        
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
                        %                 scatter3(...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
                        %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
                        % end
                    end
                end
            end
        end
        
        axP = get(gca,'Position');
        %  legend(legendtext,'Location','southeastoutside')
        set(gca, 'Position', axP)
        view(3)
        %view([24 26])
        view([50 49])
        
        grid on
        hold off
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'height_v_distanceToGo_v_bankAngle_Set',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).case)',...
            '.png'),...
            'png');
        % close(fig_num);
        %  end
    end
    
end

%% Height vs. Angular Distance To Go vs. Flight-Path Angle - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 140000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Angular Distance To Go vs. Flight-Path Angle - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Flight-Path Angle (deg))') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        %set(gca,'XTick', -90:15:30);
        set(gca,'YTick', -90:15:90);
        %set(gca,'ZTick', 0:15:150);
        xlim([0 60])
        %ylim([0 90])
        zlim([0 100])
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            
            %color_line3(...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.longitude_angle,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.latitude_angle,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass);
            plot3(...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flight_path_angle,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3)
            
            
            %plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
            %   compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.flight_path_angle,'g');
            
            
            % plot3(...
            %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
            %     ones(numel(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo))*100,...
            %     compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,'g')
            
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
            %end
            
            
            
        end
        
        axP = get(gca,'Position');
        %  legend(legendtext,'Location','southeastoutside')
        set(gca, 'Position', axP)
        view([24 26])
        grid on
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'height_v_distanceToGo_v_flightpathAngle_Set',...
            convertCharsToStrings(compilation(p).case)',...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Height vs. Angular Distance To Go vs. Total Body G-load - per Generation
for p = 1:numel(compilation)
    
    
    %for k = 1:numel(compilation(p).evolutions)
    for k = numel(compilation(p).evolutions)
        fig_num = p*100 + k*1 + 141000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Angular Distance To Go vs. Total Body G-load - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Total Body G-load (g)') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        set(gca,'XTick', 0:5:60);
        %set(gca,'YTick', 0:15:90);
        set(gca,'ZTick', 0:15:165);
        xlim([0 60])
        %ylim([0 90])
        zlim([0 165])
        
        
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            
            %color_line3(...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.longitude_angle,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.latitude_angle,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass);
            plot3(...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3)
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.bodyFrameTotalGLoadMagnitude);
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
        end
    end
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([13 49])
    grid on
    hold off
    
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).figurePath,...
        'height_v_distanceToGo_v_total_body_g_load_Set',...
        convertCharsToStrings(compilation(p).case)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end

%% Height vs. Angular Distance To Go vs. Equilibrium Wall Temperature - per Generation
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        
        fig_num = p*100 + k*1 + 142000;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Angular Distance To Go vs. Eq. Wall Temp. - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Eq. Wall Temp. (K)') % y-axis label
        zlabel('Height (km)') % z-axis label
        % img = imread('img.jpg');
        %imagesc([-180 180], [-90 90], (flipud(img)));
        set (gca,'Fontsize',20)
        set(gca,'XTick', 0:10:60);
        set(gca,'YTick', 0:100:1800);
        %set(gca,'ZTick', 0:15:150);
        %xlim([0 60])
        %ylim([0 90])
        zlim([0 175])
        
        grid on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = compilation(p).evolutions(k).population(1).indices.printed
            %      for ii = 1
            %color_line3(...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.longitude_angle,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.latitude_angle,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3,...
            %    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.mass);
            plot3(...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.wall_temperature_chapman,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.height/1e3)
            plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.angularDistanceToGo,...
                compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.wall_temperature_chapman);
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(1),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(1),'filled')
            %                 scatter3(...
            %                     compilation(p).evolutions(i).trajectories(k).individual.x_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.y_R(end),...
            %                     compilation(p).evolutions(i).trajectories(k).individual.z_R(end),'filled')
        end
    end
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([13 49])
    
    grid on
    hold off
    
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).figurePath,...
        'height_v_distanceToGo_v_total_body_g_load_Set',...
        convertCharsToStrings(compilation(p).case)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end

end
