function [  ] = plotTimeHistoriesOfConstrainedParameters( compilation )

dependentVariableFieldNames = [{'bodyFrameTotalGLoadMagnitude'} {'heatFluxChapmanNose'} {'flightPathAngle'} {'integratedHeatLoad'} {'centralTargetAngularDistanceToGo'}];
    field_x = 'timeOfFlight';

    
%% Time History of Dependent Variables

%compilation = compilations(2);

for p = 2%:numel(compilation)

 field_x_struct = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_x);
            variableLabel_x = field_x_struct.variableLabel;
            figureSaveName_x = field_x_struct.figureSaveNameContent;
            units_x = field_x_struct.units;
            limits_x = field_x_struct.limits;
            scalingFactor_x = field_x_struct.scalingFactor;
            tick_x = field_x_struct.tick;

    
    % title_line1 = strip(strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' ')),'left');
    %max_tof = 50 * ceil(max([compilation(p).evolutions.max_tof])/50);
    
    idx_of_printed_individuals = compilation(p).evolutions(1).population(1).indices.printed;
   % colors =[ 0.4940 0.1840 0.5560 ; 0 0.4470 0.7410 ];
    colors =[ 1 0 0 ; 0 0.4470 0.7410 ];
    linetypes = [{'-'};{'--'}];
    %[0 0.4470 0.7410]
    
    combo = [ 1 4;...
        1 5;...
        3 5;...
        4 4;...
        5 4 ];
    %combo = [ 4 4 ];
    
    for j = 1:size(combo,1)
        
        
        field_y = dependentVariableFieldNames{combo(j,1)};
        field_y_struct = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y);
        variableLabel_y = field_y_struct.variableLabel;
        figureSaveName_y = field_y_struct.figureSaveNameContent;
        units_y = field_y_struct.units;
        limits_y = field_y_struct.limits;
        scalingFactor_y = field_y_struct.scalingFactor;
        tick_y = field_y_struct.tick;
        
        if strcmp(field_y,'bodyFrameTotalGLoadMagnitude')
            limits_y = [0 5];
            tick_y = 1;
        end
        if strcmp(field_y,'heatFluxChapmanNose')
            limits_y = [0 2000];
            tick_y = 250;
        end
        if strcmp(field_y,'integratedHeatLoad')
            limits_y = [0 5.0e2];
            tick_y = 1.0e2;
        end
        
        
        %for k = 5%1:length(idx_of_printed_individuals)
        k = combo(j,2);
        fig_num = p*100 + 3466000 + k*1;
        figure(fig_num)
        set (gca,'Fontsize',15)
        
        ylim(limits_y)
        xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
        ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
        set(gca,'TickLabelInterpreter','latex')
        
        if norm(double(isinf(limits_y))) == 0
            set(gca,'YTick', limits_y(1):tick_y:limits_y(2));
        end
        
        hold on
        grid on
        %set(gca,'yminorgrid','on')
        
        selectIndividual = idx_of_printed_individuals(k);
        
        
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
        
        if  compilation(1).validation == 1
            for ii = idx_of_printed_individuals
                
                plot(compilation(p).evolutions(1).population(ii).dependentVariableTimeHistory.(field_x).value,...
                    compilation(p).evolutions(1).population(ii).dependentVariableTimeHistory.(field_y).value,'k','LineWidth',2);
                xlim([0 1400])
                set(gca,'XTick', 0:200:1400);
            end
        else
            x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
            y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
            
            ascentRange = 1:compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2);
            descentRange = compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2):numel(x);
            
            h(1) = plot(x(ascentRange),y(ascentRange),linetypes{1},'Color',colors(1,:),'LineWidth',2);
            h(2) = plot(x(descentRange),y(descentRange),linetypes{2},'Color',colors(2,:),'LineWidth',2);
            
        end
        
        
        if combo(j,1) < 4
            switch true
                case strcmp(field_y,'bodyFrameTotalGLoadMagnitude')
                    constraint = [ 1 1.5 2.2 3];
                case strcmp(field_y,'heatFluxChapmanNose')
                    constraint = [ 400 530 750 ];
                case strcmp(field_y,'flightPathAngle')
                    constraint = 0;
            end
            for i = 1:numel(constraint)
                for pp = 1:2
                    
                    if pp == 1
                        phaseIndices = ascentRange;
                    else
                        phaseIndices = descentRange;
                    end
                    
                    x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(phaseIndices)/scalingFactor_x;
                    tof(i,pp) = (max(x) - min(x));
                    
                    y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(phaseIndices)/scalingFactor_y;
                    y_constraint = constraint(i)*ones(numel(compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(phaseIndices)/scalingFactor_y),1);
                    
                    if strcmp(field_y,'flightPathAngle') && pp == 1
                        y = -y;
                    end
                    
                    h(numel(h) + 1) = plot(x,y_constraint,'--k','LineWidth',1);
                    [xi,yi] = polyxpoly(x,y,x,y_constraint);
                    x = [x;xi];
                    [x,index] = sortrows(x);
                    
                    y = [y;yi];
                    y = y(index);
                    y_constraint = [y_constraint;yi];
                    y_constraint = y_constraint(index);
                    
                    idx=y>=y_constraint;
                    maxViolation(i,pp) = 0;
                    violationIntegral(i,pp) = 0;
                    violationIntegralAverage(i,pp) = 0;
                    %                         totalViolationsAscent(i,:) = [ 0 0 0 0 ];
                    %                         totalViolationsDescent(i,:) = [ 0 0 0 0] ;
                    maxY = max(y(idx));
                    maxX = max(x(idx));
                    if ismember(1,idx)
                        
                        if pp == 1
                            maxViolation(i,pp) = (max(abs(y(idx))) - max(abs(constraint(i))));
                            violationIntegral(i,pp) = sum(abs(y(idx)) - abs(constraint(i)))*(x(2) - x(1));
                            violationIntegralAverage(i,pp) = violationIntegral(i,pp)/tof(i,pp);
                            totalViolations(i,pp) = (maxViolation(i,pp) + violationIntegralAverage(i,pp));
                            totalViolationsAscent(i,1) = constraint(i);
                            totalViolationsAscent(i,2) = maxViolation(i,pp);
                            totalViolationsAscent(i,3) = violationIntegralAverage(i,pp);
                            totalViolationsAscent(i,4) = totalViolations(i,pp);
                            totalViolationsAscent(i,5) = totalViolations(i,pp)/constraint(i);
                            
                        else
                            maxViolation(i,pp) = (max(abs(y(idx))) - max(abs(constraint(i))));
                            violationIntegral(i,pp) = sum(abs(y(idx)) - abs(constraint(i)))*(x(2) - x(1));
                            violationIntegralAverage(i,pp) = violationIntegral(i,pp)/tof(i,pp);
                            totalViolations(i,pp) = (maxViolation(i,pp) + violationIntegralAverage(i,pp));
                            totalViolationsDescent(i,1) = constraint(i);
                            totalViolationsDescent(i,2) = maxViolation(i,pp);
                            totalViolationsDescent(i,3) = violationIntegralAverage(i,pp);
                            totalViolationsDescent(i,4) = totalViolations(i,pp);
                            totalViolationsDescent(i,5) = totalViolations(i,pp)/constraint(i);
                        end
                        %                             totalViolationsstruct(i,pp).name = field_y;
                        %                             totalViolationsstruct(i,pp).value = constraint(i);
                        %                             totalViolationsstruct(i,pp).violations = [totalViolations(i,pp)  maxViolation(i,pp) violationIntegralAverage(i,pp)];
                        
                    end
                    idx_dif = [nan;diff(idx)];
                    I1 = find(idx_dif == 1 );
                    if(idx(1)==1)
                        I1 = [1;I1];
                    end
                    I2 = find(idx_dif == -1 );
                    if(idx(end)==1)
                        I2 = [I2;(numel(idx)+1)];
                    end
                    
                    I = [(I1) (I2-1)];
                    
                    if strcmp(field_y,'flightPathAngle') && pp == 1
                        %y = -y;
                        upper = y_constraint;
                        lower = -y;
                    else
                        upper = y;
                        lower = y_constraint;
                    end
                    
                    for ii = 1:size(I,1)
                        range = I(ii,1):I(ii,2);
                        
                        [h(numel(h) + 1),msg]=jbfill(transpose(x(range)),transpose(upper(range)),transpose(lower(range)),colors(pp,:),colors(pp,:)/2,0,.2);
                    end
                    
                end
                %                         legendText = [{'Calculated Value'};{'Phase Change'};{'Constraint'};{'Region Applicable to Ascent'};{'Region Applicable to Descent'};{'5'};{'Constraint'};{''};{''}];
                %                         legend([h(1) h(3) h(4) h(5) h(7)],legendText,'Location','northeast')
                %
            end
            
            
            patches = find( isgraphics(h,'patch') );
            
            if isempty(patches) == 0
                
                trajectoryPatches(1).indices = find(ismember(vec2mat([h(patches).FaceColor],3),colors(1,:),'rows'));
                trajectoryPatches(2).indices = find(ismember(vec2mat([h(patches).FaceColor],3),colors(2,:),'rows'));
                
                hatchAngle = [135 45];
                hatchType = [{'single'} {'single'}];
                
                for kkk = 1:numel(trajectoryPatches)
                    for jkjk = 1:numel(trajectoryPatches(kkk).indices)
                        hatchfill2(h(patches(trajectoryPatches(kkk).indices(jkjk))),hatchType{kkk},'HatchAngle',hatchAngle(kkk),'HatchDensity',40,'HatchColor','k','HatchLineWidth',.5);
                    end
                end
                
            end
            
            g(1) = plot(NaN,NaN,'Color',colors(1,:),'LineWidth',2);
            g(2) = fill([nan nan],[nan nan],colors(1,:));
            g(3) = plot(NaN,NaN,'Color',colors(2,:),'LineWidth',2);
            g(4) = fill([nan nan],[nan nan],colors(2,:));
            legendText = [{'Ascent'};{'Violation during ascent'};{'Descent'};{'Violation during descent'}];
            [~,object_h,~,~] = legendflex(g,legendText,'Interpreter','latex');
            patches = find( isgraphics(object_h,'patch') );
            
            for kkk = 1:numel(trajectoryPatches)
                hatchfill2(object_h(patches(kkk)),'single','HatchAngle',hatchAngle(kkk),'HatchDensity',10,'HatchColor','k','HatchLineWidth',.5);
            end
            
            scatter(phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            
            switch true
                case strcmp(field_y,'bodyFrameTotalGLoadMagnitude')
                    cleanUnits = strcat('$',erase(erase(units_y,'('),')'),'$');
                    positioning = [ {'SW'} {'NW'} {'NW'} {'NW'} ];
                    
                    if k == 4
                        
                        for ppp = 1:length(constraint)
                            finalNumberAsText = strcat(sprintf('%.1f', constraint(ppp)),'$');
                            labelpoints(1400,constraint(ppp), strcat('$x_{n_{T,',num2str(ppp),'}} =',finalNumberAsText,{' '},cleanUnits), positioning{ppp}, 'FontSize', 15,'interpreter','latex','buffer',0.05)
                        end
                        labelpoints(420,phaseChange_y(2), [{'$\leftarrow$'}], 'E', 'FontSize', 15,'interpreter','latex','buffer',0.1)
                        labelpoints(420,.06, [{'Phase Change'}], 'E', 'FontSize', 15,'interpreter','latex','buffer',0.45)
                        
                    elseif k == 5
                        for ppp = 1:length(constraint)
                            finalNumberAsText = strcat(sprintf('%.1f', constraint(ppp)),'$');
                            labelpoints(800,constraint(ppp), strcat('$x_{n_{T,',num2str(ppp),'}} =',finalNumberAsText,{' '},cleanUnits), positioning{ppp},'FontSize', 15,'interpreter','latex','buffer',0.05)
                        end
                        labelpoints(phaseChange_x(2),phaseChange_y(2), [{'Phase'};{'Change'}],'NW', 'FontSize', 15,'interpreter','latex','stacked', 'down','buffer',0.1)
                        
                    end
                case strcmp(field_y,'heatFluxChapmanNose')
                    cleanUnits = strcat('$',erase(erase(units_y,'('),')'),'$');
                    
                    for ppp = 1:length(constraint)
                        finalNumberAsText = strcat(sprintf('%.1f', constraint(ppp)),'$');
                        labelpoints(0,constraint(ppp), strcat('$x_{\dot{q}_{',num2str(ppp),'}} =',finalNumberAsText,{' '},cleanUnits), 'NE', 'FontSize', 15,'interpreter','latex')
                    end
                    labelpoints(phaseChange_x(2),phaseChange_y(2), [{'Phase Change'}],'SE', 'FontSize', 15,'interpreter','latex')
                case strcmp(field_y,'flightPathAngle')
                    constraint = 0;
                    labelpoints(phaseChange_x(2),phaseChange_y(2), [{'Phase Change'}],'SE', 'FontSize', 15,'interpreter','latex')
            end
            
            
            
            %labelpoints(420,phaseChange_y(2), [{'$\leftarrow$'}], 'E', 'FontSize', 15,'interpreter','latex','buffer',0.1)
            %labelpoints(420,.06, [{'Phase Change'}], 'E', 'FontSize', 15,'interpreter','latex','buffer',0.45)
            
        end
        
        if combo(j,1) >= 4
            
            switch true
                case strcmp(field_y,'integratedHeatLoad')
                    % ytickformat('%.1f')
                    %set(gca,'yminorgrid','on')
                    
                    x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
                    y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
                    y_constraint = y(end)*ones(numel(compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y),1);
                    plot(x,y_constraint,'--k','LineWidth',1);
                    labelpoints(phaseChange_x(2),phaseChange_y(2), {'Phase Change'}, 'SE', 'FontSize', 15,'interpreter','latex')
                    scatter(0,y_constraint(1),'k');
                    
                    cleanUnits = strcat(erase(erase(units_y,'('),')'));
                    finalNumberAsText = strcat(num2str(fix((y_constraint(1)/(10^2)*100))/100),{'\times10^2$'});
                    labelpoints(0,y_constraint(1), strcat({'$\dot{Q}_{f} ='},finalNumberAsText,{' '},cleanUnits), 'NE', 'FontSize', 15,'interpreter','latex','buffer',0.25)
                    
                    plot(x(end)*[1 1],[0 y(end)],'--k','LineWidth',1);
                    scatter(x(end),0,'k');
                    
                    cleanUnits = '$ s';
                    %finalNumberAsText = strcat(num2str(fix((y_constraint(1)/(10^5)*100))/100),{'\times10^5$'});
                    labelpoints(x(end),0, strcat({'$t_{tof} ='},num2str(x(end)),{' '},cleanUnits), 'NW', 'FontSize', 15,'interpreter','latex','buffer',0.25)
                    
                    constraint = [ 400 530 750 ] *1000;
                    cost = y(end)*1000;%/((x(end) - x(1))*constraint(2));
                    cost = y(end)*1000./((x(end) - x(1))*constraint);
                    %colors = [ 0.5843    0.8157    0.9882];
                    0.525551613070314;
                    
                case strcmp(field_y,'centralTargetAngularDistanceToGo')
                    %set(gca,'yminorgrid','on')
                    
                    x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
                    y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
                    y_constraint = y(end)*ones(numel(compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y),1);
                    plot(x,y_constraint,'--k','LineWidth',1);
                    
                    labelpoints(phaseChange_x(2),phaseChange_y(2), {'Phase Change'}, 'NE', 'FontSize', 15,'interpreter','latex')
                    scatter(0,y_constraint(1),'k');
                    
                    cleanUnits = strcat(erase(erase(units_y,'('),')'));
                    finalNumberAsText = strcat(num2str(fix((y_constraint(1)/(10^0)*100))/100),{'$'});
                    labelpoints(0,y_constraint(1), strcat({'$\theta_{ToGo,f} ='},finalNumberAsText,{' '},cleanUnits), 'NE', 'FontSize', 15,'interpreter','latex','buffer',0.25)
                    
                    initialNumberAsText = strcat(num2str(fix((y(1)/(10^0)*100))/100),{'$'});
                    labelpoints(0,y(1), strcat({'$\theta_{ToGo,0} ='},initialNumberAsText,{' '},cleanUnits), 'NE', 'FontSize', 15,'interpreter','latex','buffer',0.25)
                    scatter(0,y(1),'k');
                    cost = 100 * y(end)/y(1);
                    
                    
            end
            
            legendText = [{'Ascent'};{'Descent'}];
            legendflex(h(1:2),legendText,'Interpreter','latex');
            scatter(phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            
            
        end
        
        
        hold off
        
        max_tof = 50 * ceil(max(x)/50);
        xlim([0 max_tof])
        set(gca,'XTick', 0:max_tof/10:max_tof);
        
        if j < 4
            saveNameMainTerm = strcat('fitnessComputationExamples_Individual_',num2str(k),'_','constraintViolation_');
        elseif j >= 4
            saveNameMainTerm = strcat('fitnessComputationExamples_Individual_',num2str(k),'_','cost_');
        end
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            saveNameMainTerm,...
            strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
            '_Case',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        %end
    end
end

%% e

for p = 2%1:numel(compilation)
    
    idx_of_printed_individuals = compilation(p).evolutions(1).population(1).indices.printed;
    
    for j = 2
        
        field_y = dependentVariableFieldNames{j};
        variableLabel_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).variableLabel;
        figureSaveName_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).figureSaveNameContent;
        units_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).units;
        limits_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).limits;
        scalingFactor_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).scalingFactor;
        tick_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).tick;
        
        if strcmp(field_y,'heatFluxChapmanNose')
            limits_y = [0 2000];
            tick_y = 250;
        end
        
        for k = 4;%1:length(idx_of_printed_individuals)
            
            fig_num = p*100 + 3466000 + k*1;
            figure(fig_num)
            
            set (gca,'Fontsize',15)
            
            ylim(limits_y)
            xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
            ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
            set(gca,'TickLabelInterpreter','latex')
            
            if norm(double(isinf(limits_y))) == 0
                set(gca,'YTick', limits_y(1):tick_y:limits_y(2));
            end
            
            hold on
            grid on
            %set(gca,'yminorgrid','on')
            
            selectIndividual = idx_of_printed_individuals(k);
            %colors =[ 0.4940 0.1840 0.5560 ; 0.5843 0.8157 0.9882 ];
            %    colors =[ 0.4940 0.1840 0.5560 ; 0 0.4470 0.7410 ];
    colors =[ 1 0 0 ; 0 0.4470 0.7410 ];
    linetypes = [{'-'};{'--'}];
            
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
                y_constraint = 0*ones(numel(compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y),1);
                
                %h(numel(h) + 1) = plot(x,y_constraint,'--k','LineWidth',1);
                shadedColor = [ 0.5843 0.8157 0.9882 ];
                
                [h(numel(h) + 1),msg]=jbfill(transpose(x),transpose(y),transpose(y_constraint),shadedColor,shadedColor/2,0,.5);
                
                %[h(numel(h) + 1),msg]=jbfill(transpose(x(range)),transpose(upper(range)),transpose(lower(range)),colors(pp,:),colors(pp,:)/2,0,.2);
                
                
            if  compilation(1).validation == 1
                for ii = idx_of_printed_individuals
                    
                    plot(compilation(p).evolutions(1).population(ii).dependentVariableTimeHistory.(field_x).value,...
                        compilation(p).evolutions(1).population(ii).dependentVariableTimeHistory.(field_y).value,'k','LineWidth',2);
                    xlim([0 1400])
                    set(gca,'XTick', 0:200:1400);
                end
            else
                
                x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
                y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
                
                ascentRange = 1:compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2);
                descentRange = compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2):numel(x);
                
            h(1) = plot(x(ascentRange),y(ascentRange),linetypes{1},'Color',colors(1,:),'LineWidth',2);
            h(2) = plot(x(descentRange),y(descentRange),linetypes{2},'Color',colors(2,:),'LineWidth',2);
            
            end
            
           legendText = [{'Ascent'};{'Descent'}];
            legendflex(h(1:2),legendText,'Interpreter','latex');
            scatter(phaseChange_x(2),phaseChange_y(2),70,'filled','ks');
            
            labelpoints(phaseChange_x(2),phaseChange_y(2), {'Phase Change'}, 'SE', 'FontSize', 15,'interpreter','latex')
            
            hold off
            
            max_tof = 50 * ceil(max(x)/50);
            xlim([0 max_tof])
            set(gca,'XTick', 0:max_tof/10:max_tof);
            
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'fitnessComputationExamples_Individual_',num2str(k),'_','areaUnderCurve_',...
                strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
                '_Case',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
    end
end

%%
% for p = 1:numel(compilation)
%
%     idx_of_printed_individuals = compilation(p).evolutions(1).population(1).indices.printed;
%
%     for j = 4
%
%         field_y = dependentVariableFieldNames{j};
%         variableLabel_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).variableLabel;
%         figureSaveName_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).figureSaveNameContent;
%         units_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).units;
%         limits_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).limits;
%         scalingFactor_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).scalingFactor;
%         tick_y = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y).tick;
%
%         limits_y = [0 5.0e5];
%         tick_y = 1.0e5;
%
%
%         for k = 1:length(idx_of_printed_individuals)
%
%             fig_num = p*100 + 3466000 + k*1;
%             figure(fig_num)
%
%             set (gca,'Fontsize',15)
%
%             ylim(limits_y)
%             xlabel(string(strcat(variableLabel_x,{' '},units_x))) % x-axis label
%             ylabel(string(strcat(variableLabel_y,{' '},units_y))) % y-axis label
%
%             if norm(double(isinf(limits_y))) == 0
%                 set(gca,'YTick', limits_y(1):tick_y:limits_y(2));
%             end
%
%             hold on
%             grid on
%
%             selectIndividual = idx_of_printed_individuals(k);
%
%
%             if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1)) == false
%                 phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_x;
%                 phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(1))/scalingFactor_y;
%             else
%                 phaseChange_x(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
%                 phaseChange_y(1) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
%             end
%
%             if isnan(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2)) == false
%                 phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_x;
%                 phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(compilation(p).evolutions(1).population(selectIndividual).indices.trajectoryPhaseChange(2))/scalingFactor_y;
%             else
%                 phaseChange_x(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value(end)/scalingFactor_x;
%                 phaseChange_y(2) = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value(end)/scalingFactor_y;
%             end
%
%             if  compilation(1).validation == 1
%                 for ii = idx_of_printed_individuals
%
%                     plot(compilation(p).evolutions(1).population(ii).dependentVariableTimeHistory.(field_x).value,...
%                         compilation(p).evolutions(1).population(ii).dependentVariableTimeHistory.(field_y).value,'k','LineWidth',2);
%                     xlim([0 1400])
%                     set(gca,'XTick', 0:200:1400);
%                 end
%             else
%                 h = plotFieldVsField( compilation, p, 1, selectIndividual, field_x, scalingFactor_x, field_y, scalingFactor_y, 2 );
%             end
%             labelpoints(phaseChange_x(2),phaseChange_y(2), {'Phase Change'}, 'SE', 'FontSize', 15,'interpreter','latex')
%
%
%             % ytickformat('%.1f')
%             set(gca,'yminorgrid','on')
%
%             x = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_x).value/scalingFactor_x;
%             y = compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y;
%             y_constraint = y(end)*ones(numel(compilation(p).evolutions(1).population(selectIndividual).dependentVariableTimeHistory.(field_y).value/scalingFactor_y),1);
%             %plot(x,y_constraint,'--k','LineWidth',1);
%             scatter(0,y_constraint(1),'k');
%
%             cleanUnits = strcat('$',erase(erase(units_y,'('),')'),'$');
%             finalNumberAsText = strcat(num2str(fix((y_constraint(1)/(10^5)*100))/100),{'\times10^5$'});
%             labelpoints(0,y_constraint(1), strcat({'$\dot{Q}_{f} ='},finalNumberAsText,{' '},cleanUnits), 'NE', 'FontSize', 15,'interpreter','latex')
%
%             %cost = y(end)*1000/((x(end) - x(1))*constraint(2));
%             colors = [ 0.5843 0.8157 0.9882 ];
%             0.525551613070314;
%
%             hold off
%
%             max_tof = 50 * ceil(max(x)/50);
%             xlim([0 max_tof])
%             set(gca,'XTick', 0:max_tof/10:max_tof);
%
%             saveas(...
%                 figure(fig_num),...
%                 strcat(...
%                 compilation(p).figurePath,...
%                 'fitnessComputationExamples_Individual_',num2str(k),'_','constraintViolation',...
%                 strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
%                 '_Case',...
%                 convertCharsToStrings(compilation(p).case),...
%                 '.png'),...
%                 'png');
%             close(fig_num);
%         end
%     end
% end
%






end