function [ h ] = plotFieldVsField( compilation, p, k, iii, X_field, X_scale, Y_field, Y_scale, widthOfLine )

j = 1;
for ii = iii
    %if max(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.('bodyFrameTotalGLoadMagnitude').value) < 5
    
    h(1) = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(X_field).value/X_scale,...
        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(Y_field).value/Y_scale,'LineWidth',widthOfLine);
    ax = gca;
    ax.ColorOrderIndex = j;
    
    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
        h(2) = scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(X_field).value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/X_scale,...
            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(Y_field).value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1))/Y_scale,20,'x');
        ax.ColorOrderIndex = j;
    else
        h(2) = scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(X_field).value(end)/X_scale,...
            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(Y_field).value(end)/Y_scale,20,'x');
        ax.ColorOrderIndex = j;
    end
    
    if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
        h(3) = scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(X_field).value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/X_scale,...
            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(Y_field).value(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2))/Y_scale,20,'s');
    else
        h(3) = scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(X_field).value(end)/X_scale,...
            compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(Y_field).value(end)/Y_scale,20,'s');
        
    end
    j = j + 1;
   %  end
end