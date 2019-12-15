function [  ] = plotSpecificFieldOfIndividual( compilation, p, k, ii, field )


h = plot(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight,...
    compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(field{1}));
set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))

if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)) == false
    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),...
        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(field{1})(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(1)),'x');
    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
else
    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(end),...
        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(field{1})(end),'x');
    set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
end

if isnan(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)) == false
    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),...
        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(field{1})(compilation(p).evolutions(k).population(ii).indices.trajectoryPhaseChange(2)),'s');
else
    scatter(compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.timeOfFlight(end),...
        compilation(p).evolutions(k).population(ii).dependentVariableTimeHistory.(field{1})(end),'s');
    
end


end