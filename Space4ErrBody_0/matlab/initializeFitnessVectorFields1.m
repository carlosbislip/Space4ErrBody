function [ evolutions ] = initializeFitnessVectorFields( simulationDetails, evolutions, k, i )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

trajectoryType = simulationDetails.trajectoryType;
objectiveFunctionCase = simulationDetails.objectiveFunctionCase;

if trajectoryType == 'A'
    if objectiveFunctionCase == 'A'
        evolutions(k).population(i).size.fitness = 2;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'B'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energy   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'C'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.monotonicEnergy   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'D'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.monotonicEnergyAndFlightPath   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'E'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'F'
        evolutions(k).population(i).size.fitness = 4;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'G'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'H'
        evolutions(k).population(i).size.fitness = 4;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsZ   = nan;
    elseif objectiveFunctionCase == 'I'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerJerkXYZ   = nan;
    elseif objectiveFunctionCase == 'J'
        evolutions(k).population(i).size.fitness = 4;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelZJerkXYZ   = nan;
    elseif objectiveFunctionCase == 'K'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.allDynamics   = nan;
    elseif objectiveFunctionCase == 'L'
        evolutions(k).population(i).size.fitness = 6;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerJerkXYZ   = nan;
    elseif objectiveFunctionCase == 'M'
        evolutions(k).population(i).size.fitness = 16;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.maximumNormalizedEnergyCost   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.maximumNormalizedEnergyPenalty   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.monotonicEnergy   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.flightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.heatLoadCost   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.heatLoadPenalty   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.dynamicPressure   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.mechanicalLoad   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.bendingMoment   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsPosZ   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsNegZ = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerJerkX   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerJerkY   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerJerkZ   = nan;
    end
elseif trajectoryType == 'D'
    if objectiveFunctionCase == 'A'
        evolutions(k).population(i).size.fitness = 2;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'B'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energy   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'C'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.monotonicEnergy   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'D'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.monotonicEnergyAndFlightPath   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'E'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'F'
        evolutions(k).population(i).size.fitness = 4;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'G'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'H'
        evolutions(k).population(i).size.fitness = 4;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelerationsZ   = nan;
    elseif objectiveFunctionCase == 'I'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerJerkXYZ   = nan;
    elseif objectiveFunctionCase == 'J'
        evolutions(k).population(i).size.fitness = 4;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelZJerkXYZ   = nan;
    elseif objectiveFunctionCase == 'K'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.allDynamics   = nan;
    elseif objectiveFunctionCase == 'L'
        evolutions(k).population(i).size.fitness = 6;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerJerkXYZ   = nan;
    end
    
elseif trajectoryType == 'AD'
    if objectiveFunctionCase == 'A'
        evolutions(k).population(i).size.fitness = 3;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'B'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energy   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energy   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'C'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.monotonicEnergy   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.monotonicEnergy   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'D'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.monotonicEnergyAndFlightPath   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.monotonicEnergyAndFlightPath   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'E'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'F'
        evolutions(k).population(i).size.fitness = 7;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'G'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
    elseif objectiveFunctionCase == 'H'
        evolutions(k).population(i).size.fitness = 7;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelerationsZ   = nan;
    elseif objectiveFunctionCase == 'I'
        evolutions(k).population(i).size.fitness = 9;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerJerkXYZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerJerkXYZ   = nan;
    elseif objectiveFunctionCase == 'J'
        evolutions(k).population(i).size.fitness = 7;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelZJerkXYZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelZJerkXYZ   = nan;
    elseif objectiveFunctionCase == 'K'
        evolutions(k).population(i).size.fitness = 5;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.allDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyMassFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.allDynamics   = nan;
    elseif objectiveFunctionCase == 'L'
        evolutions(k).population(i).size.fitness = 11;
        evolutions(k).population(i).fitnessVector.Common.angularDistanceToGo   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Ascent.passengerJerkXYZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.energyAndFlightPathAngle   = nan;
        evolutions(k).population(i).fitnessVector.Descent.fuelMassCost   = nan;
        evolutions(k).population(i).fitnessVector.Descent.basicDynamics   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerAccelerationsZ   = nan;
        evolutions(k).population(i).fitnessVector.Descent.passengerJerkXYZ   = nan;
    end
end
evolutions(k).population(i).fitnessVector.Cost                             = nan;
evolutions(k).population(i).fitnessVector.DominationSet                    = nan;
evolutions(k).population(i).fitnessVector.DominatedCount                   = nan;
evolutions(k).population(i).fitnessVector.Rank                             = nan;
evolutions(k).population(i).fitnessVector.RankBySortedVector               = nan;
evolutions(k).population(i).fitnessVector.RankBySortedMagnitude            = nan;
evolutions(k).population(i).fitnessVector.Common.magnitude                 = nan;
evolutions(k).population(i).fitnessVector.Common.ranking                   = nan;
evolutions(k).population(i).fitnessVector.Common.scalingFactor             = nan;
evolutions(k).population(i).fitnessVector.Common.magnitudeAverages         = nan;


end

