function [  ] = plotCompoundRelations( compilation )



%% Height vs. Angular Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 794000 + k*1;
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
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
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
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
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
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass);
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
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                compilation(p).evolutions(k).trajectories(ii).individual.mass);
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
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.height/1e3,...
                compilation(p).evolutions(k).trajectories(ii).individual.evaluated_throttle_setting);
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

%% Height vs. Velocity - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
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
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
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
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
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
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_angle);
            
            
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
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
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
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.heading_error);
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger);
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_distance,...
                compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_LB,'k','LineWidth',2);
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_distance,...
                compilation(p).evolutions(k).trajectories(ii).individual.headingErrorDeadBand_UP,'k','LineWidth',2);
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
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




%% Bank Angle vs. Distance To Go - per Evolution
for p = 1:numel(compilation)
    
    %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
    for k = 1:numel(compilation(p).evolutions)
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
            
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle);
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
    
    for k = 1:numel(compilation(p).evolutions)
        %for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
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
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go,...
                compilation(p).evolutions(k).trajectories(ii).individual.bank_angle_reversal_trigger);
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

%% Height vs. Heating Rate at Leading Edge - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 724100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Height vs. Heating Rate at Leading Edge - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        % xlim([0 8000])
        xlabel('Heating Rate at Leading Edge (W/m^2)') % x-axis label
        ylabel('Height (km)') % y-axis label
        set(gca,'YTick', 0:10:150);
        % set(gca,'XTick', 0:500:8000);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.q_dot_LE,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/height_v_q_dot_LE_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end

%% Heating Rate at Leading Edge vs. Velocity - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 724200 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Heating Rate at Leading Edge vs. Velocity - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 150])
        % max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 8000])
        xlabel('Velocity (m/s)') % x-axis label
        ylabel('Heating Rate at Leading Edge (W/m^2)') % y-axis label
        %set(gca,'YTick', 0:10:150);
        set(gca,'XTick', 0:500:8000);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
            stairs(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                compilation(p).evolutions(k).trajectories(ii).individual.q_dot_LE);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/q_dot_LE_v_V_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end


%% Time History: Height vs. Groundtrack - per Evolution
lat_f_deg = 38.9444444444444;
lon_f_deg = -77.4558333333;
th = 0:pi/50:2*pi;
xunit = .75 * cos(th) + lon_f_deg;
yunit = .75 * sin(th) + lat_f_deg;
for p = 1:numel(compilation)
    
    
    for k = 1:numel(compilation(p).evolutions)
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
        for ii =1:numel(compilation(p).evolutions(k).trajectories)
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
    title(strcat('Trajectories - Final Evolution of Set ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
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
    
    
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        %     for k = 1:numel(compilation(p).evolutions)
        %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            
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
    
    axP = get(gca,'Position');
    %  legend(legendtext,'Location','southeastoutside')
    set(gca, 'Position', axP)
    view([13 49])
    grid on
    hold off
    saveas(...
        figure(fig_num),...
        strcat(...
        compilations(1).evolutions.mainpath,...
        '/figures/Height_v_distanceToGo_v_mass_Set',...
        convertCharsToStrings(compilation(p).set)',...
        '.png'),...
        'png');
    % close(fig_num);
    %  end
end





end