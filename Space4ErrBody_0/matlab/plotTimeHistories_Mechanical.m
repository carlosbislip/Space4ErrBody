function [  ] = plotTimeHistories_Mechanical( compilation )


%% Time History: Total Body G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3457000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Body G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Total Body G-load (g)') % y-axis label
        set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
         %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
          stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_mag);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_body_g_load_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %        close(fig_num);
    end
end



%% Time History: Total Body Z-Component G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*200 + 3457000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Body Z-Component G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-5 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Total Body Z-Component G-load (g)') % y-axis label
        set(gca,'YTick', -5:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
         %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
          stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.body_fixed_total_g_load_z);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_body_g_load_z_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end


%% Time History: Total Passenger Z-Component G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*200 + 3458000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Total Passenger Z-Component G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([-5 10])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Total Passenger Z-Component G-load (g)') % y-axis label
        set(gca,'YTick', -5:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
         %for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
          stairs(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.passenger_fixed_total_g_load_z);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/total_passenger_g_load_z_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %    close(fig_num);
    end
end



%% Time History: Aero G-load - per Evolution
for p = 1:numel(compilation)
    
    for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 3457000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Aero G-load through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 5])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Aero G-load (g)') % y-axis label
        set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)):numel(compilation(p).evolutions(k).trajectories)
            % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.total_aero_g_load);
        end
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/aero_g_load_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end

%% Time History: Thrust Acc. Components - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3465000 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Components through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 5])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Thrust Acc. Components (m/s^2)') % y-axis label
        %set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:50:max_tof);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)):numel(compilation(p).evolutions(k).trajectories)
            % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_x,'k');
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_y,'r');
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_z,'b');
        end
        legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Thrust_acc_components_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end



%% Time History: Thrust Acc. Magnitude - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 3465100 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Thrust Acc. Magnitude through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        %ylim([0 5])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Thrust Acc. Magnitude (m/s^2)') % y-axis label
        %set(gca,'YTick', 0:1:10);
        set(gca,'XTick', 0:50:max_tof);
        hold on
        
        for ii = (numel(compilation(p).evolutions(k).trajectories)):numel(compilation(p).evolutions(k).trajectories)
            % for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            if ( compilation(p).evolutions(k).trajectories(ii).individual.distance_to_go(end) < 30 )
                plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                    compilation(p).evolutions(k).trajectories(ii).individual.acc_thru_M);
            end
        end
        %legend('x-dir','y-dir','z-dir')
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/Thrust_acc_magnitude_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
     %   close(fig_num);
    end
end




%% Time History: Dynamic Pressure - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
   % for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
        fig_num = p*100 + 654180 + k*1;
        figure(fig_num)
        set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        title(strcat('Dynamic Pressure through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
        ylim([0 100])
        max_tof = max([compilation(p).evolutions.max_tof]);
        xlim([0 max_tof])
        xlabel('Propagation Time (s)') % x-axis label
        ylabel('Dynamic Pressure (kPa)') % y-axis label
        %set(gca,'YTick', -90:30:90);
        set(gca,'XTick', 0:200:max_tof);
        hold on
        
       % for ii = (numel(compilation(p).evolutions(k).trajectories)-10):numel(compilation(p).evolutions(k).trajectories)
             for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.time_vector,...
                compilation(p).evolutions(k).trajectories(ii).individual.dynamic_pressure/1e3);
        end
        
        %plot([0 max_tof],(6371 + 25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).mainpath,...
            '/figures/dynamic_pressure_v_T_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        %      close(fig_num);
    end
end



end