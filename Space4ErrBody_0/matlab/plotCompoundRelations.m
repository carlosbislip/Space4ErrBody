function [  ] = plotCompoundRelations( compilation, mainpath )




%% Height vs. Velocity - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
        fig_num = p*100 + 724000 + k*1;
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
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/height_v_V_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        close(fig_num);
    end
end


%% Height vs. Heating Rate at Leading Edge - per Evolution
for p = 1:numel(compilation)
    
    for k = 1:numel(compilation(p).evolutions)
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
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.q_dot_LE,...
                compilation(p).evolutions(k).trajectories(ii).individual.height/1e3);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
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
    
    for k = 1:numel(compilation(p).evolutions)
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
        
        for ii = 1:numel(compilation(p).evolutions(k).trajectories)
            plot(compilation(p).evolutions(k).trajectories(ii).individual.airspeed,...
                compilation(p).evolutions(k).trajectories(ii).individual.q_dot_LE);
        end
        
        plot([0 8000],(25)*[1 1],'k','LineWidth',2)
        hold off
        saveas(...
            figure(fig_num),...
            strcat(...
            mainpath,...
            '/figures/q_dot_LE_v_V_Evolution_',...
            num2str(k - 1),...
            '_Set',...
            convertCharsToStrings(compilation(p).set),...
            '.png'),...
            'png');
        % close(fig_num);
    end
end






end