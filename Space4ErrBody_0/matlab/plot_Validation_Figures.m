function [  ] = plot_Validation_Figures( compilation )

%% Figures for Validation: Time Histories

dependentVariableFieldNames  = [ {'timeOfFlight'} {'angleOfAttack'} {'bankAngle'} {'flightPathAngle'} {'dynamicPressure'} {'machNumber'}...
    {'headingAngleToCentralTargetError'} {'bodyFrameTotalGLoadMagnitude'} {'centralTargetAngularDistanceToGo'} {'bodyflapDeflectionAngle'} {'heatFluxChapmanNose'}...
    {'airspeed'} {'height'} ];

numberOfFields = numel(dependentVariableFieldNames);


combinationsX = [ ones(1,numberOfFields-1) (numberOfFields-1) 9 9 ];
combinationsY = [ 2:numberOfFields numberOfFields 7 7];
combinations = [combinationsX' combinationsY'];

originaldepvarlist = fieldnames(compilation(end).evolutions(end).population(end).dependentVariableTimeHistory);

max_tof = 1400;%50 * ceil(max([compilation(p).evolutions.max_tof])/50);

for kk = 1:size(combinations,1)
    
    Xx = combinations(kk,1);
    Xy = combinations(kk,2);
    
    field_x = dependentVariableFieldNames{Xx};
    field_x_struct = compilation(end).evolutions(end).population(end).dependentVariableTimeHistory.(field_x);
    variableLabel_x = field_x_struct.variableLabel;
    figureSaveName_x = field_x_struct.figureSaveNameContent;
    units_x = field_x_struct.units;
    limits_x = field_x_struct.limits;
    scalingFactor_x = field_x_struct.scalingFactor;
    tick_x = field_x_struct.tick;
    
    field_y = dependentVariableFieldNames{Xy};
    field_y_struct = compilation(end).evolutions(end).population(end).dependentVariableTimeHistory.(field_y);
    variableLabel_y = field_y_struct.variableLabel;
    figureSaveName_y = field_y_struct.figureSaveNameContent;
    units_y = field_y_struct.units;
    limits_y = field_y_struct.limits;
    scalingFactor_y = field_y_struct.scalingFactor;
    tick_y = field_y_struct.tick;
    
    fig_num = 100 + 3466000;
    figure(fig_num)
    
    set (gca,'Fontsize',15)
    
    
    if strcmp(field_x,'timeOfFlight')
        limits_x = [0 max_tof];
        xlim(limits_x)
        set(gca,'XTick', 0:200:max_tof);
    elseif strcmp(field_x,'centralTargetAngularDistanceToGo')
        xlim(limits_x)
        set(gca,'XTick',  limits_x(1):tick_x:limits_x(2));
    elseif strcmp(field_x,'airspeed')
        %xlim(limits_x)
        % set(gca,'XTick',  limits_x(1):tick_x:limits_x(2));
    end
    
    ylim(limits_y)
    xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
    ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
    
    if norm(double(isinf(limits_y))) == 0
        set(gca,'YTick', limits_y(1):tick_y:limits_y(2));
    end
    
    set(gca,'TickLabelInterpreter','latex')
    
    
    hold on
    grid on
    
    x = compilation(1).evolutions(1).population(1).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
    y = compilation(1).evolutions(1).population(1).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
    
    h(1) = plot(x,y,'k','LineWidth',2);
    
    
    if strcmp(field_x,'centralTargetAngularDistanceToGo')
        x2 = compilation(1).rawData.headingErrorDeadbandBounds(:,1);
        y21 = compilation(1).rawData.headingErrorDeadbandBounds(:,2);
        y22 = compilation(1).rawData.headingErrorDeadbandBounds(:,3);
        
        plot(x2,y21,'k','LineWidth',2);
        plot(x2,y22,'k','LineWidth',2);
        plot(limits_x,[0 0],':k','LineWidth',1);
        
        
        if kk == size(combinations,1)
            limits_x = [0 5];
            
            set(gca,'XTick',  limits_x(1):1:limits_x(2));
            xlim(limits_x)
            
        end
        
        
        
        
        
        
    elseif strcmp(field_y,'bankAngle')
        plot(limits_x,[0 0],':k','LineWidth',1);
    elseif strcmp(field_y,'flightPathAngle')
        plot(limits_x,[0 0],':k','LineWidth',1);
        limits_y = [-10 0];
        ylim(limits_y)
        set(gca,'YTick', limits_y(1):1:limits_y(2));
        
    elseif strcmp(field_y,'dynamicPressure')
        limits_y = [0 10];
        ylim(limits_y)
        set(gca,'YTick', limits_y(1):limits_y(2)/10:limits_y(2));
        
    elseif strcmp(field_y,'machNumber')
        limits_y = [0 30];
        ylim(limits_y)
        set(gca,'YTick', limits_y(1):5:limits_y(2));
        
    elseif strcmp(field_y,'bodyFrameTotalGLoadMagnitude')
        limits_y = [0 2];
        ylim(limits_y)
        set(gca,'YTick', limits_y(1):limits_y(2)/10:limits_y(2));
        
    elseif strcmp(field_y,'bodyflapDeflectionAngle')
        plot(limits_x,[0 0],':k','LineWidth',1);
        
        set(gca,'YTick', limits_y(1):5:limits_y(2));
        
    elseif strcmp(field_y,'heatFluxChapmanNose')
        plot(limits_x,[530 530],':k','LineWidth',1);
        
        limits_y = [0 700];
        ylim(limits_y)
    elseif strcmp(field_y,'height')
        if strcmp(field_x,'timeOfFlight')
            limits_y = [0 140];
            ylim(limits_y)
        end
        if strcmp(field_x,'airspeed')
            limits_y = [0 140];
            ylim(limits_y)
            
            
            label{1} = variableLabel_y;
            flightCorridorNames = [{'flightCorridorBoundary_DynamicPressure'} {'flightCorridorBoundary_BendingMoment'}...
                {'flightCorridorBoundary_ThermalLoad'} {'flightCorridorBoundary_SkipSuppression'}...
                {'flightCorridorBoundary_MechanicalLoad_ascent'} {'flightCorridorBoundary_MechanicalLoad_descent'}];
            
            
            for i = 1:numel(flightCorridorNames)
                
                % label{1+i} = compilation(end).evolutions(end).population(end).dependentVariableTimeHistory.(flightCorridorNames{i}).variableLabel;
                
                % y = compilation(1).evolutions(1).population(1).dependentVariableTimeHistory.(flightCorridorNames{i}).value/scalingFactor_y;
                
                % h(1+i) = plot(x,y,'LineWidth',1);
                
                
            end
            %  legendflex(h(:), label(:),'xscale', 0.5,'Interpreter','latex');
            
        end
        
    end
    
    saveName = string(strcat(...
        compilation.figurePath,...
        'validation_',strcat(figureSaveName_x,'_vs_',figureSaveName_y),'.png'));
  
    if strcmp(field_x,'centralTargetAngularDistanceToGo') && kk == size(combinations,1)
        saveName = string(strcat(...
            compilation.figurePath,...
            'validation_',strcat(figureSaveName_x,'_vs_',figureSaveName_y),'_refined','.png'));
    end
    
    hold off
    
    saveas(...
        figure(fig_num),...
        saveName,'png');
    close(fig_num);
    
end







end