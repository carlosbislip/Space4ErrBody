function [ population ] = loadDecisionVectorFromPopulationCellArray( populationData )

decisionVector = cell2mat(populationData(:,3:end));

population.decisionVector.parameters.Ascent.nodeInterval.data.one                    = decisionVector(:,1);
population.decisionVector.parameters.Ascent.nodeInterval.data.two                    = decisionVector(:,2);
population.decisionVector.parameters.Ascent.nodeInterval.data.three                  = decisionVector(:,3);
population.decisionVector.parameters.Ascent.nodeInterval.data.four                   = decisionVector(:,4);
population.decisionVector.parameters.Ascent.nodeInterval.data.five                   = decisionVector(:,5);
population.decisionVector.parameters.Ascent.nodeInterval.data.six                    = decisionVector(:,6);
population.decisionVector.parameters.Ascent.nodeInterval.data.seven                  = decisionVector(:,7);
population.decisionVector.parameters.Ascent.nodeInterval.data.eight                  = decisionVector(:,8);
population.decisionVector.parameters.Ascent.nodeInterval.data.nine                   = decisionVector(:,9);

population.decisionVector.parameters.Ascent.angleOfAttack.data.one                    = decisionVector(:,10);
population.decisionVector.parameters.Ascent.angleOfAttack.data.two                    = decisionVector(:,11);
population.decisionVector.parameters.Ascent.angleOfAttack.data.three                  = decisionVector(:,12);
population.decisionVector.parameters.Ascent.angleOfAttack.data.four                   = decisionVector(:,13);
population.decisionVector.parameters.Ascent.angleOfAttack.data.five                   = decisionVector(:,14);
population.decisionVector.parameters.Ascent.angleOfAttack.data.six                    = decisionVector(:,15);
population.decisionVector.parameters.Ascent.angleOfAttack.data.seven                  = decisionVector(:,16);
population.decisionVector.parameters.Ascent.angleOfAttack.data.eight                  = decisionVector(:,17);
population.decisionVector.parameters.Ascent.angleOfAttack.data.nine                   = decisionVector(:,18);
population.decisionVector.parameters.Ascent.angleOfAttack.data.ten                    = decisionVector(:,19);

population.decisionVector.parameters.Ascent.bankAngle.data.one                    = decisionVector(:,20);
population.decisionVector.parameters.Ascent.bankAngle.data.two                    = decisionVector(:,21);
population.decisionVector.parameters.Ascent.bankAngle.data.three                  = decisionVector(:,22);
population.decisionVector.parameters.Ascent.bankAngle.data.four                   = decisionVector(:,23);
population.decisionVector.parameters.Ascent.bankAngle.data.five                   = decisionVector(:,24);
population.decisionVector.parameters.Ascent.bankAngle.data.six                    = decisionVector(:,25);
population.decisionVector.parameters.Ascent.bankAngle.data.seven                  = decisionVector(:,26);
population.decisionVector.parameters.Ascent.bankAngle.data.eight                  = decisionVector(:,27);
population.decisionVector.parameters.Ascent.bankAngle.data.nine                   = decisionVector(:,28);
population.decisionVector.parameters.Ascent.bankAngle.data.ten                    = decisionVector(:,29);

population.decisionVector.parameters.Ascent.thrustElevationAngle.data.one                    = decisionVector(:,30);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.two                    = decisionVector(:,31);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.three                  = decisionVector(:,32);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.four                   = decisionVector(:,33);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.five                   = decisionVector(:,34);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.six                    = decisionVector(:,35);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.seven                  = decisionVector(:,36);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.eight                  = decisionVector(:,37);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.nine                   = decisionVector(:,38);
population.decisionVector.parameters.Ascent.thrustElevationAngle.data.ten                    = decisionVector(:,39);

population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.one                    = decisionVector(:,40);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.two                    = decisionVector(:,41);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.three                  = decisionVector(:,42);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.four                   = decisionVector(:,43);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.five                   = decisionVector(:,44);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.six                    = decisionVector(:,45);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.seven                  = decisionVector(:,46);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.eight                  = decisionVector(:,47);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.nine                   = decisionVector(:,48);
population.decisionVector.parameters.Ascent.thrustAzimuthAngle.data.ten                    = decisionVector(:,49);

population.decisionVector.parameters.Ascent.throttleSetting.data.one                    = decisionVector(:,50);
population.decisionVector.parameters.Ascent.throttleSetting.data.two                    = decisionVector(:,51);
population.decisionVector.parameters.Ascent.throttleSetting.data.three                  = decisionVector(:,52);
population.decisionVector.parameters.Ascent.throttleSetting.data.four                   = decisionVector(:,53);
population.decisionVector.parameters.Ascent.throttleSetting.data.five                   = decisionVector(:,54);
population.decisionVector.parameters.Ascent.throttleSetting.data.six                    = decisionVector(:,55);
population.decisionVector.parameters.Ascent.throttleSetting.data.seven                  = decisionVector(:,56);
population.decisionVector.parameters.Ascent.throttleSetting.data.eight                  = decisionVector(:,57);
population.decisionVector.parameters.Ascent.throttleSetting.data.nine                   = decisionVector(:,58);
population.decisionVector.parameters.Ascent.throttleSetting.data.ten                    = decisionVector(:,59);

population.decisionVector.parameters.Common.initialFlightPathAngle.data.one                    = decisionVector(:,60);
population.decisionVector.parameters.Common.initialLaunchHeadingAngle.data.one                  = decisionVector(:,61);
population.decisionVector.parameters.Common.initialVelocity.data.one                            = decisionVector(:,62);
population.decisionVector.parameters.Common.maximumVelocity.data.one                            = decisionVector(:,63);
population.decisionVector.parameters.Common.maximumHeight.data.one                              = decisionVector(:,64);
population.decisionVector.parameters.Common.additionalMass.data.one                             = decisionVector(:,65);
population.decisionVector.parameters.Common.terminationDistanceRatio.data.one                   = decisionVector(:,66);

population.decisionVector.parameters.Descent.nodeInterval.data.one                   = decisionVector(:,67);
population.decisionVector.parameters.Descent.nodeInterval.data.two                   = decisionVector(:,68);
population.decisionVector.parameters.Descent.nodeInterval.data.three                 = decisionVector(:,69);
population.decisionVector.parameters.Descent.nodeInterval.data.four                  = decisionVector(:,70);
population.decisionVector.parameters.Descent.nodeInterval.data.five                  = decisionVector(:,71);
population.decisionVector.parameters.Descent.nodeInterval.data.six                   = decisionVector(:,72);
population.decisionVector.parameters.Descent.nodeInterval.data.seven                 = decisionVector(:,73);
population.decisionVector.parameters.Descent.nodeInterval.data.eight                 = decisionVector(:,74);
population.decisionVector.parameters.Descent.nodeInterval.data.nine                  = decisionVector(:,75);

population.decisionVector.parameters.Descent.angleOfAttack.data.one                    = decisionVector(:,76);
population.decisionVector.parameters.Descent.angleOfAttack.data.two                    = decisionVector(:,77);
population.decisionVector.parameters.Descent.angleOfAttack.data.three                  = decisionVector(:,78);
population.decisionVector.parameters.Descent.angleOfAttack.data.four                   = decisionVector(:,79);
population.decisionVector.parameters.Descent.angleOfAttack.data.five                   = decisionVector(:,80);
population.decisionVector.parameters.Descent.angleOfAttack.data.six                    = decisionVector(:,81);
population.decisionVector.parameters.Descent.angleOfAttack.data.seven                  = decisionVector(:,82);
population.decisionVector.parameters.Descent.angleOfAttack.data.eight                  = decisionVector(:,83);
population.decisionVector.parameters.Descent.angleOfAttack.data.nine                   = decisionVector(:,84);
population.decisionVector.parameters.Descent.angleOfAttack.data.ten                    = decisionVector(:,85);


population.decisionVector.parameters.Descent.bankAngle.data.one                    = decisionVector(:,86);
population.decisionVector.parameters.Descent.bankAngle.data.two                    = decisionVector(:,87);
population.decisionVector.parameters.Descent.bankAngle.data.three                  = decisionVector(:,88);
population.decisionVector.parameters.Descent.bankAngle.data.four                   = decisionVector(:,89);
population.decisionVector.parameters.Descent.bankAngle.data.five                   = decisionVector(:,90);
population.decisionVector.parameters.Descent.bankAngle.data.six                    = decisionVector(:,91);
population.decisionVector.parameters.Descent.bankAngle.data.seven                  = decisionVector(:,92);
population.decisionVector.parameters.Descent.bankAngle.data.eight                  = decisionVector(:,93);
population.decisionVector.parameters.Descent.bankAngle.data.nine                   = decisionVector(:,94);
population.decisionVector.parameters.Descent.bankAngle.data.ten                    = decisionVector(:,95);

population.decisionVector.parameters.Descent.thrustElevationAngle.data.one                    = decisionVector(:,96);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.two                    = decisionVector(:,97);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.three                  = decisionVector(:,98);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.four                   = decisionVector(:,99);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.five                   = decisionVector(:,100);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.six                    = decisionVector(:,101);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.seven                  = decisionVector(:,102);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.eight                  = decisionVector(:,103);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.nine                   = decisionVector(:,104);
population.decisionVector.parameters.Descent.thrustElevationAngle.data.ten                    = decisionVector(:,105);


population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.one                    = decisionVector(:,106);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.two                    = decisionVector(:,107);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.three                  = decisionVector(:,108);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.four                   = decisionVector(:,109);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.five                   = decisionVector(:,110);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.six                    = decisionVector(:,111);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.seven                  = decisionVector(:,112);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.eight                  = decisionVector(:,113);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.nine                   = decisionVector(:,114);
population.decisionVector.parameters.Descent.thrustAzimuthAngle.data.ten                    = decisionVector(:,115);


population.decisionVector.parameters.Descent.throttleSetting.data.one                    = decisionVector(:,116);
population.decisionVector.parameters.Descent.throttleSetting.data.two                    = decisionVector(:,117);
population.decisionVector.parameters.Descent.throttleSetting.data.three                  = decisionVector(:,118);
population.decisionVector.parameters.Descent.throttleSetting.data.four                   = decisionVector(:,119);
population.decisionVector.parameters.Descent.throttleSetting.data.five                   = decisionVector(:,120);
population.decisionVector.parameters.Descent.throttleSetting.data.six                    = decisionVector(:,121);
population.decisionVector.parameters.Descent.throttleSetting.data.seven                  = decisionVector(:,122);
population.decisionVector.parameters.Descent.throttleSetting.data.eight                  = decisionVector(:,123);
population.decisionVector.parameters.Descent.throttleSetting.data.nine                   = decisionVector(:,124);
population.decisionVector.parameters.Descent.throttleSetting.data.ten                    = decisionVector(:,125);
population.decisionVector.parameters.Common.finalVelocity.data.one                          = decisionVector(:,126);
population.decisionVector.parameters.Common.skipSuppressionTriggerTime.data.one                  = decisionVector(:,127);
population.decisionVector.parameters.Common.maximumMechanicalLoad.data.one                       = decisionVector(:,128);






end