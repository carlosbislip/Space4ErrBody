function [ fitness ] = loadFitnessVectorFromFitnessCellArray( fitnessData )

fitnessVector = cell2mat(fitnessData(:,3:end));

fitness.distanceToGo                    = fitnessVector(:,1);
fitness.AscentEnergyAndFlightPathAngle  = fitnessVector(:,2);
fitness.AscentBasicDynamics             = fitnessVector(:,3);
fitness.AscentPassengerAccelerationsZ   = fitnessVector(:,4);
fitness.AscentPassengerJerkXYZ          = fitnessVector(:,5);
fitness.DescentEnergyAndFlightPathAngle = fitnessVector(:,6);
fitness.DescentBasicDynamics            = fitnessVector(:,7);
fitness.DescentPassengerAccelerationsZ  = fitnessVector(:,8);
fitness.DescentPassengerJerkXYZ         = fitnessVector(:,9);
fitness.DescentPassengerJerkXYZ         = fitnessVector(:,10);
fitness.DescentPassengerJerkXYZ         = fitnessVector(:,11);
fitness.Magnitude                       = sqrt(sum(fitnessVector.^2, 2));



            evolutions(k).population.fitnessVector.Common.distanceToGo              = fitnessVector(:,1);
            evolutions(k).population.fitnessVector.Ascent.energyAndFlightPathAngle  = fitnessVector(:,2);
            evolutions(k).population.fitnessVector.Ascent.fuelMassCost              = fitnessVector(:,3);
            evolutions(k).population.fitnessVector.Ascent.basicDynamics             = fitnessVector(:,4);
            evolutions(k).population.fitnessVector.Ascent.passengerAccelerationsZ   = fitnessVector(:,5);
            evolutions(k).population.fitnessVector.Ascent.passengerJerkXYZ          = fitnessVector(:,6);
            evolutions(k).population.fitnessVector.Descent.energyAndFlightPathAngle = fitnessVector(:,7);
            evolutions(k).population.fitnessVector.Descent.basicDynamics            = fitnessVector(:,8);
            evolutions(k).population.fitnessVector.Descent.fuelMassCost             = fitnessVector(:,9);
            evolutions(k).population.fitnessVector.Descent.passengerAccelerationsZ  = fitnessVector(:,10);
            evolutions(k).population.fitnessVector.Descent.passengerJerkXYZ         = fitnessVector(:,11);
            evolutions(k).population.fitnessVector.Common.magnitude                  = sqrt(sum(fitnessVector.^2, 2));
            evolutions(k).population.fitnessVector.Common.ranking                  = sqrt(sum(fitnessVector.^2, 2));



[~,p] = sort(fitness.Magnitude,'ascend');
r(:,1) = 1:length(fitness.Magnitude);
r(p) = r;

fitness.Ranking         = r;

end