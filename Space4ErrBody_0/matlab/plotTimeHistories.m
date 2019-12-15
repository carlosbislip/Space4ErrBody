function [  ] = plotTimeHistories( compilation )

%% Time History of Dependent Variables

dependentVariableFieldNames  = fieldnames(compilation(1).evolutions(end).population(end).dependentVariableTimeHistory);

field_x = dependentVariableFieldNames{2};
variableLabel_x = compilation(1).evolutions(end).population(end).dependentVariableTimeHistory.(field_x).variableLabel;
figureSaveName_x = compilation(1).evolutions(end).population(end).dependentVariableTimeHistory.(field_x).figureSaveNameContent;
units_x = compilation(1).evolutions(end).population(end).dependentVariableTimeHistory.(field_x).units;
limits_x = compilation(1).evolutions(end).population(end).dependentVariableTimeHistory.(field_x).limits;
scalingFactor_x = compilation(1).evolutions(end).population(end).dependentVariableTimeHistory.(field_x).scalingFactor;
id = [1544 86] + 1;

for p = 1:numel(compilation)
    
    title_line1 = strip(strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' ')),'left');
    max_tof = 50 * ceil(max([compilation(p).evolutions.max_tof])/50);
    
    
    for j = 1:length(dependentVariableFieldNames)
        %for j = 16
        
        %field_y = dependentVariableFieldNames{j};
                    field_y = 'angleOfAttack';
variableLabel_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).variableLabel;
        figureSaveName_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).figureSaveNameContent;
        units_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).units;
        limits_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).limits;
        scalingFactor_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).scalingFactor;
        tick_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).tick;
        
        title_line3 = string(strcat(variableLabel_y,{' vs '},variableLabel_x));
        
        for k = 1:numel(compilation(p).evolutions)
            
            
            if ~isnan(compilation(p).evolutions(k).max_tof)
                
                fig_num = p*100 + 3466000 + k*1;
                figure(fig_num)
                %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
                set (gca,'Fontsize',15)
                
                xlim([0 max_tof])
                set(gca,'XTick', 0:max_tof/10:max_tof);
                
                %title_line2 = string(strcat({'Generation '},num2str(k - 1)));
                %title({title_line1;title_line2;title_line3})
                
                ylim(limits_y)
                xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
                ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
                
                if norm(double(isinf(limits_y))) == 0
                    set(gca,'YTick', limits_y(1):tick_y:limits_y(2));
                end
                
                hold on
                grid on
                
                if  compilation(1).validation == 1
                    for ii = compilation(p).evolutions(k).population(1).indices.printed
                        
                        plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(field_x).value,...
                            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(field_y).value,'k','LineWidth',2);
                        xlim([0 1400])
                        set(gca,'XTick', 0:200:1400);
                    end
                else
                    %plotFieldVsField( compilation, p, k, compilation(p).evolutions(k).population(1).indices.printed, field_x, scalingFactor_x, field_y, scalingFactor_y, 2 );
                   % plotFieldVsField( compilation, 1, k, id(1), field_x, scalingFactor_x, field_y, scalingFactor_y, 2 );
                    %plotFieldVsField( compilation, 2, k, id(2), field_x, scalingFactor_x, field_y, scalingFactor_y, 2 );
                end
                    set(gca,'TickLabelInterpreter','latex')
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
            
            %scatter( phaseChange_x(1),phaseChange_y(1),70,'kx')
            scatter( phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            
            ax = gca;
            ax.ColorOrderIndex = 2;
            j = j + 1;
            
        
            
            %ascentRange = 1:compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2);
            %descentRange = compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2):numel(x);
            
            %h(1) = plot(x,y,'LineWidth',2);
            
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
            
            h(2) = plot(x,y,'LineStyle','--','LineWidth',2);
            
            
            ax = gca;
            %h(3) = scatter( phaseChange_x(1),phaseChange_y(1),70,'kx');
   
            h(3) = scatter( phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            j = j + 1;
            
            
            legend(h(:),[{'Initial Individual'};{'Re-evaluated Individual'};{'Phase Change'}],'Interpreter','latex','Location','northeast')
            
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                hold off
                
                saveas(...
                    figure(fig_num),...
                    strcat(...
                    compilation(p).figurePath,...
                    'dependentVariables_',...
                    strcat(figureSaveName_x,'_vs_',figureSaveName_y),'_generation_',...
                    num2str(k - 1),...
                    '_case',...
                    convertCharsToStrings(compilation(p).case),...
                    '.png'),...
                    'png');
                close(fig_num);
                
            end
        end
    end
end



end