function [  ] = plotCompoundRelations( compilation )

%% Skip Suppression Limit vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    % for k = 1:numel(compilation(p).evolutions)
    for k = 1
        fig_num = p*100 + 794000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Skip Suppression Limit vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(rad2deg([-pi/2 pi/2]))
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Skip Suppression Limit (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        hold on
        plot([0 60],[0 0],'k','LineWidth',2)
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
       % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
             for ii = 194
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 45 )
                %       if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit);
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.evaluated_bank_angle);
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
                %
                %                 stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                %                     compilation(p).evolutions(k).trajectories(ii).individual.commanded_bank_angle);
                
                
                
                % end
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/skip_suppression_limit_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Total Passenger Z-Component G-load vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 894000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Passenger Z-Component G-load vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-5 10])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Total Passenger Z-Component G-load (g)') % y-axis label
        set(gca,'YTick', -5:1:10);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
                %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.passenger_fixed_total_g_load_z);
                %  end
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_passenger_g_load_z_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Thrust Magnitude vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 894000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Magnitude vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1100])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Thrust Magnitude (kN)') % y-axis label
        set(gca,'YTick', 0:100:1100);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
       % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
        for ii = 194
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 45 )
                %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.thrustMagnitude/1e3);
                %  end
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/thrustMagnitude_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Height vs. Density - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 794000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Density - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 100])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 .05])
        xlabel('Density (kg/m^3)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:100);
        set(gca,'XTick', 0:.01:.05);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %   if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
            %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.localDensity,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
            % end
            %     end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_localDensity_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Density vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 794000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Density vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 800])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Density (kg/m^3)') % y-axis label
        % set(gca,'YTick', 0:100:800);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %   if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
            %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.localDensity);
            % end
            %     end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/density_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Height vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
   %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 974000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 800])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Angular Distance To Go (deg)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:100:800);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
       % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
        for ii = 194
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 45 )
                %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
            end
            %     end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Mass vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 785000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Mass vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Mass (kg)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            % if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
            %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass);
            %    end
            % end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/mass_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Mass vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 785000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Mass vs. Angular Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Mass (kg)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            % if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
            %     if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass);
            %    end
            % end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/mass_v_distanceToGo_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Mass vs. Height - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 784000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Mass vs. Height - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        % xlim([0 8000])
        xlabel('Height (km)') % x-axis label
        ylabel('Mass (kg)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        % set(gca,'XTick', 0:500:8000);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
                if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                        compilation(p).evolutions(k).trajectories(ii).individual.mass);
                end
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/mass_v_height_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Throttle vs. Height - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 774000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle vs. Height  - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        % xlim([0 8000])
        xlabel('Height (-)') % x-axis label
        ylabel('Throttle (-)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        % set(gca,'XTick', 0:500:8000);
        hold on
        
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.evaluated_throttle_setting,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/throttle_v_height_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Throttle Setting vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 814000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Throttle Setting vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([180 360])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Throttle Setting (-)') % y-axis label
        %set(gca,'YTick', 180:30:360);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
           if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30 )
            %    if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        compilation(p).evolutions(k).trajectories(ii).individual.evaluated_throttle_setting);
           %     end
            end
            
            
        end
        
        %                 for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        %             %for ii = 1:numel(compilation(p).evolutions(k).trajectories)
        %            plot(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
        %                compilation(p).evolutions(k).trajectories(ii).individual.heading_to_target);
        %
        %         end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/evaluatedThrottleSetting_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end






%% Height vs. Velocity - per Evolution
for p = 1:numel(compilation)
    
    %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1    
        fig_num = p*100 + 764000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Velocity - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 8000])
        xlabel('Velocity (m/s)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:150);
        set(gca,'XTick', 0:500:8000);
        hold on
        
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
            %    if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                        compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
             %   end
            end
        end
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_V_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Heading Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
  %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 854000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([180 360])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Heading Angle (deg)') % y-axis label
        set(gca,'YTick', 180:30:360);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
        %   for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            %if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 20 )
             %   if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        compilation(p).evolutions(k).trajectories(ii).individual.heading_angle);
              %  end
           % end
            
            
        end
        
        %                 for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        %             %for ii = 1:numel(compilation(p).evolutions(k).trajectories)
        %            plot(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
        %                compilation(p).evolutions(k).trajectories(ii).individual.heading_to_target);
        %
        %         end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/headingAngle_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Heading Error vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
     %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 754000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heading Error vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(30*[-1 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Heading Error (deg)') % y-axis label
        set(gca,'YTick', -30:10:30);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
             if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
            %    if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_error);
            % stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger);
            %   end
            end
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_distance,...
                compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_LB,'k','LineWidth',2);
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_distance,...
                compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_UP,'k','LineWidth',2);
            
            %stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
            %   compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
            
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/headingError_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end




%% Angle of Attack vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
     %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 544000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Angle of Attack vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Angle of Attack (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
              %  if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        compilation(p).evolutions(k).trajectories(ii).individual.angle_of_attack);
                %end
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/angleOfAttack_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end
%% Bank Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
     %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Bank Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
              %  if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
                %end
            end
        end
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankAngle_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Bank Reversal Trigger vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1  
        fig_num = p*100 + 744000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Bank Reversal Trigger vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        % ylim(30*[-1 1])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Bank Reversal Trigger (-)') % y-axis label
        %set(gca,'YTick', -30:10:30);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
       %for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
             stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger);
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bankReversalTrigger_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Flight-Path Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
      %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1 
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Flight-Path Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim(90*[-1 1])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Flight-Path Angle (deg)') % y-axis label
        set(gca,'YTick', -90:15:90);
        set(gca,'XTick', 0:10:60);
        % set(gca, 'XDir','reverse')
        
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-20):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
               % if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        compilation(p).evolutions(k).trajectories(ii).individual.flight_path_angle);
              %  end
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/flightPathAngle_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% Pitch Moment Coefficient vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
      %   for k = 1:numel(compilation(p).evolutions)
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1 
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Cm vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-.1 .03])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('C_m (-)') % y-axis label
        set(gca,'YTick', -.1:0.01:.03);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
               % if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                    stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                        compilation(p).evolutions(k).trajectories(ii).individual.aero_moment_coefficient_C_m );
               % end
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Cm_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end


%% BodyFlap Deflection Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = 1:numel(compilation(p).evolutions)
        % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
         for k = 1
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('BodyFlap Deflection Angle vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-30 30])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('BodyFlap Deflection Angle (deg)') % y-axis label
        set(gca,'YTick', -30:10:30);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 60 )
                % if( sum(compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger) < 6 )
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bodyflap_deflection);
                % end
            end
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/bodyflap_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end






%% Total Body G-load vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Body G-load vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Total Body G-load (g)') % y-axis label
        set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag );
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_body_g_load_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Height vs. Tauber Heat Flux at Leading Edge - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 724100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Tauber Heat Flux at Leading Edge - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        % xlim([0 8000])
        xlabel('Tauber Heat Flux at Leading Edge (W/m^2)') % x-axis label
        ylabel('Height (km)') % y-axis label
        % set(gca,'YTick', 0:10:150);
        % set(gca,'XTick', 0:500:8000);
        hold on
        
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_tauber_leadingedge,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_heatFluxTauber_leadingedge_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Tauber Heat Flux at Leading Edge vs. Velocity - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 724200 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Tauber Heat Flux at Leading Edge vs. Velocity - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 8000])
        xlabel('Velocity (m/s)') % x-axis label
        ylabel('Tauber Heat Flux at Leading Edge (W/m^2)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:500:8000);
        hold on
        
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_tauber_leadingedge/1e3);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heatFluxTauber_leadingedge_v_V_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Tauber Heat Flux at Leading Edge vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Tauber Heat Rate at Leading Edge vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Tauber Heat Rate at Leading Edge (kW/m^2)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_tauber_leadingedge/1e3 );
        end
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heatRateTauber_leadingedge_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end

%% Chapman Heat Flux at Nose vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    % for k = 1:numel(compilation(p).evolutions)
    for k = 1
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Chapman Heat Flux at Nose vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Chapman Heat Flux at Nose (kW/m^2)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
                
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.heat_flux_chapman_nose/1e3 );
            end
        end 
            hold off
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).mainpath,...
                '/figures/heatFluxChapman_nose_v_distance_Evolution_',...
                num2str(k - 1),...
                '_Set',...
                convertCharsToStrings(compilation(p).set),...
                '.png'),...
                'png');
            %     close(fig_num);
        end
    end

%% TUDAT Heat Rate at Nose vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 954000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('TUDAT Heat Rate at Nose vs. Distance To Go - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 1000])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 60])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('TUDAT Heat Rate at Nose (kW/m^2)') % y-axis label
        set(gca,'YTick', 0:100:1000);
        set(gca,'XTick', 0:10:60);
        hold on
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heat_rate_TUDAT_nose/1e3 );
            
            
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                (compilation(p).evolutions(k).trajectories(ii).individual.heat_rate_TUDAT_nose - 0.5*(compilation(p).evolutions(k).trajectories(ii).individual.localDensity).*(compilation(p).evolutions(k).trajectories(ii).individual.airspeed).^3)/1e3 );
            
            
            
        end
        
        
        
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/heatRateTUDAT_nose_v_distance_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %     close(fig_num);
    end
end



%% Time History: Height vs. Groundtrack - per Evolution
lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;
th = 0:pi/50:2*pi;
xunit = .75 * cos(th) + lon_f_deg;
yunit = .75 * sin(th) + lat_f_deg;
for p = 1:numel(compilation)
    
    
   % for k = 1:numel(compilation(p).evolutions)
    for k = 1
        fig_num = p*100 + 724300 + k*1;
        figure(fig_num)
        hold on
        title(strcat('Trajectories - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
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
        
        
        
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %   for k = 1:numel(compilation(p).evolutions)
        %    for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
       % for ii =1:numel(compilation(p).evolutions(k).trajectories)
        for ii = 194
                    if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 45 )

            color_line3(...
                compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                compilation(p).evolutions(k).trajectories(ii).individual.time_vector);
            stairs(...
                compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
                compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,'g')
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
        %  end
        
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
            compilation(p).mainpath,...
            '/figures/Height_v_Groundtrack_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set)',...
            '.png'),...
            'png');
    end
end




%% Height vs. Distance To Go vs. Mass - per Evolution
for p = 1:numel(compilation)
    
    
    % for k = 1:numel(compilation(p).evolutions)
    fig_num = p*100 + 724400;% + k*1;
    figure(fig_num)
    hold on
    title(strcat('Height vs. Distance To Go vs. Mass - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
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
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
                            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
    
            %color_line3(...
            %    compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.mass);
            plot3(...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
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
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end



%% Height vs. Distance To Go vs. Bank Angle - per Evolution
for p = 1:numel(compilation)
    
    
    % for k = 1:numel(compilation(p).evolutions)
    fig_num = p*100 + 724400;% + k*1;
    figure(fig_num)
    hold on
    title(strcat('Height vs. Distance To Go vs. Bank Angle - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
    set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
    xlabel('Distance To Go (deg)') % x-axis label
    ylabel('Bank Angle (deg))') % y-axis label
    zlabel('Height (km)') % z-axis label
    % img = imread('img.jpg');
    %imagesc([-180 180], [-90 90], (flipud(img)));
    set (gca,'Fontsize',20)
    %set(gca,'XTick', -90:15:30);
    %set(gca,'YTick', 0:15:90);
    %set(gca,'ZTick', 0:15:150);
    xlim([38 60])
    %ylim([0 90])
    zlim([0 100])
    
    
    
    % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    %     for k = 1:numel(compilation(p).evolutions)
    for k = 1
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 40 )
                
                %color_line3(...
                %    compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
                %    compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
                %    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                %    compilation(p).evolutions(k).trajectories(ii).individual.mass);
                plot3(...
                    compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
                  plot3(...
                    compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
                
                stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.bank_angle,'g');
                 stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    compilation(p).evolutions(k).trajectories(ii).individual.skip_suppression_limit,'k');
                
                plot3(...
                    compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                    ones(numel(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go))*100,...
                    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,'g')
                
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
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([24 26])
    grid on
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        compilation(p).mainpath,...
        '/figures/Height_v_distanceToGo_v_bankAngle_Set',...
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end

%% Height vs. Distance To Go vs. Total Body G-load - per Evolution
for p = 1:numel(compilation)
    
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 724400;% + k*1;
        figure(fig_num)
        hold on
        title(strcat('Height vs. Distance To Go vs. Total Body G-load - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        xlabel('Distance To Go (deg)') % x-axis label
        ylabel('Total Body G-load (g)') % y-axis label
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
        
        
        
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
            %color_line3(...
            %    compilation(p).evolutions(k).trajectories(ii).individual.longitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.latitude_angle,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
            %    compilation(p).evolutions(k).trajectories(ii).individual.mass);
            plot3(...
                compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag);
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
        compilation(p).mainpath,...
        '/figures/Height_v_distanceToGo_v_total_body_g_load_Set',...
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end




end