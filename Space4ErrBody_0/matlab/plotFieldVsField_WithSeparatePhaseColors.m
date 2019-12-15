function [ h ] = plotFieldVsField_WithSeparatePhaseColors( compilation, p, k, i, X_field, X_scale, Y_field, Y_scale, colors, widthOfLine )

%j = 1;
%for ii = iii
    %if max(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.('bodyFrameTotalGLoadMagnitude').value) < 5
    x = compilation(p).evolutions(k).population(i).dependentVariableTimeHistory.(X_field).value/X_scale;
    y = compilation(p).evolutions(k).population(i).dependentVariableTimeHistory.(Y_field).value/Y_scale;
    
    ascentRange = 1:compilation(p).evolutions(k).population(i).indices.trajectoryPhaseChange(2);
    descentRange = compilation(p).evolutions(k).population(i).indices.trajectoryPhaseChange(2):numel(x);
    
    h(1) = plot(x(ascentRange),y(ascentRange),'Color',colors(1,:),'LineWidth',widthOfLine);
    h(2) = plot(x(descentRange),y(descentRange),'Color',colors(2,:),'LineWidth',widthOfLine);
    
    h(3) = scatter(x(ascentRange(end)),y(ascentRange(end)),70,'filled','ks');
    
    %  j = j + 1;
    %  end
end