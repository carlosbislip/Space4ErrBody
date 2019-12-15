function [  ] = plotInterpolatorsForDiscussion( compilation )


interpolatorsEvaluationFields  = fieldnames(compilation(end).evolutions(end).population(end).interpolators.ascent.evaluation);
interpolatorInputDataFields  = fieldnames(compilation(end).evolutions(end).population(end).interpolators.ascent.inputData);
phases = fieldnames(compilation(end).evolutions(end).population(end).interpolators);

%% Interpolators: Angle of Attack

    k = 1;%:numel(compilation(p).evolutions)
%ii = compilation(p).evolutions(k).population(1).indices.printed(2);
p = [ 1 1 1 1 1 2 2 2 2 2 ];
ii1 = [ 1 33 51 64 81 ];
ii2 = [ 8 37 44 52 82 ];

p = [ 1 2 ];
ii1 = [ 51 ];
ii2 = [ 8 ];

ii = [ ii1, ii2 ];

combo = [ p' ii' ];

for p = 1:size(combo,1)
    
    pp = combo(p,1);
    ii = combo(p,2);
    
    %if compilation(p).evolutions(k).population(1).indices.printed > 0
    for q = 1:numel(phases)
        
        field_x = interpolatorsEvaluationFields{1};
        variableLabel_x = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_x).variableLabel;
        figureSaveName_x = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_x).figureSaveNameContent;
        units_x = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_x).units;
        limits_x = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_x).limits;
        scalingFactor_x = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_x).scalingFactor;
        tick_x = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_x).tick;
        
        
        fig_num = p*100 + 3459000;
        figure(fig_num)
        %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        set (gca,'Fontsize',15)
        %title(strcat('Interpolated AoA for Ascent - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
        %ylim([0 50])
        %xlim([0 1])
        
        
        field_y = interpolatorsEvaluationFields{2};
        variableLabel_y = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_y).variableLabel;
        figureSaveName_y = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_y).figureSaveNameContent;
        units_y = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_y).units;
        limits_y = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_y).limits;
        scalingFactor_y = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_y).scalingFactor;
        tick_y = compilation(pp).evolutions(end).population(end).interpolators.(phases{q}).evaluation.(field_y).tick;
        
        xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
        ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
        
        ylim(limits_y)
        xlim(limits_x)
        
        set(gca,'XTick', limits_x(1):tick_x:limits_x(2));
        set(gca,'YTick', limits_y(1):tick_y:limits_y(2));
        set(gca,'TickLabelInterpreter','latex')
        
        hold on
        grid on
        
        
        x = compilation(pp).evolutions(k).population(ii).interpolators.(phases{q}).evaluation.(field_x).value;
        y = compilation(pp).evolutions(k).population(ii).interpolators.(phases{q}).evaluation.(field_y).value;
        
        upperBound = 20*ones(numel(x),1);
        lowerBound = 10*ones(numel(x),1);
       
        plot(x,y,'LineWidth',2);
         plot(x,upperBound,'--k','LineWidth',1);
        plot(x,lowerBound,'--k','LineWidth',1);
        %ax = gca;
        %ax.ColorOrderIndex = 1;
        
        x_scatter = compilation(pp).evolutions(k).population(ii).interpolators.(phases{q}).inputData.(field_x).value;
        y_scatter = compilation(pp).evolutions(k).population(ii).interpolators.(phases{q}).inputData.(field_y).value;
        %set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
        scatter( x_scatter, y_scatter, 70,'filled','ks');
        
        if size(x_scatter,1) == 5
            positioning = [ {'S'} {'N'} {'S'} {'S'} {'NW'} ];
            buffer = [ 0.1 0.1 0.1 0.1 0.1 0.1 0.0 ];
        else
            positioning = [ {'S'} {'N'} {'SW'} {'E'} {'W'} {'SE'} {'NW'} ];
            buffer = [ 0.1 0.1 0.1 0.25 0.35 0.1 0.1 ];
        end
        
        for ppp = 1:size(x_scatter,1)
            %finalNumberAsText = strcat(sprintf('%.1f', constraint(ppp)),'$');
            labelpoints(x_scatter(ppp),y_scatter(ppp), string(strcat({'cn '},num2str(ppp))), positioning{ppp}, 'FontSize', 15,'interpreter','latex', 'buffer',buffer(ppp))
        end
        
        
        labelpoints(x(end),upperBound(end), strcat({'$\alpha_{ub}$ $=$ '},num2str(upperBound(end)),{' deg'}),'NW','FontSize', 15,'interpreter','latex','buffer',0.05)
        labelpoints(x(end),lowerBound(end), strcat({'$\alpha_{lb}$ $=$ '},num2str(lowerBound(end)),{' deg'}),'SW','FontSize', 15,'interpreter','latex','buffer',0.05)
        
        
        %plot([0 max_tof],(10)*[0 1],'k','LineWidth',2)
        hold off
        figureSaveName = string(strcat(...
            compilation(pp).figurePath,...
            'interpolatorExamples_Individual_',num2str(ii),'_',(phases{q}),'_',...
            strcat(figureSaveName_y,'_vs_',figureSaveName_x),...
            '_',num2str(k - 1),'_Case',...
            convertCharsToStrings(compilation(pp).case),'.png'));
        
        saveas(figure(fig_num),figureSaveName,'png');
        close(fig_num);
        
    end
    
    %end
end
end