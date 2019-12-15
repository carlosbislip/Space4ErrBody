function [ evolutions ] = initializeFitnessVectorFields( simulationDetails, evolutions, k, i )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

trajectoryType = simulationDetails.trajectoryType;
objectiveFunctionCase = simulationDetails.objectiveFunctionCase;


switch (trajectoryType)
    case 'A'
        phase = {'ascent'};
    case 'D'
        phase = {'descent'};
    case 'AD'
        phase = [{'ascent'} {'descent'}];
end



fitnessVectorList = readcell( 'fitnessVectorList.dat','Delimiter',{','} );

for j = 1:size(fitnessVectorList,1)
    structuredList.(fitnessVectorList{j,1}).figureSaveNameContent = fitnessVectorList{j,1};
    structuredList.(fitnessVectorList{j,1}).variableLabel = fitnessVectorList{j,2};
    structuredList.(fitnessVectorList{j,1}).units = fitnessVectorList{j,3};
    structuredList.(fitnessVectorList{j,1}).scalingFactor = fitnessVectorList{j,4};
    structuredList.(fitnessVectorList{j,1}).limits = [fitnessVectorList{j,5} fitnessVectorList{j,6}];
    structuredList.(fitnessVectorList{j,1}).tick = fitnessVectorList{j,7};
    structuredList.(fitnessVectorList{j,1}).value = nan;
end

evolutions(k).population(i).fitnessVector.common.(fitnessVectorList{1,1})   = structuredList.(fitnessVectorList{1,1});

switch true
    case ismember(objectiveFunctionCase,convertStringsToChars(string(strcat({'H'},{'M'}))))
        evolutions(k).population(i).fitnessVector.common.(fitnessVectorList{2,1})   = structuredList.(fitnessVectorList{2,1});
    case ismember(objectiveFunctionCase,convertStringsToChars(string(strcat({'I'},{'N'}))))
        evolutions(k).population(i).fitnessVector.common.(fitnessVectorList{2,1})   = structuredList.(fitnessVectorList{2,1});
        evolutions(k).population(i).fitnessVector.common.(fitnessVectorList{3,1})   = structuredList.(fitnessVectorList{3,1});
end



for p = 1:numel(phase)
    switch true
        case ismember(objectiveFunctionCase,'B')
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{4,1})   = structuredList.(fitnessVectorList{4,1});
        case ismember(objectiveFunctionCase,'C')
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{4,1})   = structuredList.(fitnessVectorList{4,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{5,1})   = structuredList.(fitnessVectorList{5,1});
        case ismember(objectiveFunctionCase,'D')
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{4,1})   = structuredList.(fitnessVectorList{4,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{5,1})   = structuredList.(fitnessVectorList{5,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{6,1})   = structuredList.(fitnessVectorList{6,1});
        case ismember(objectiveFunctionCase,'E')
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{4,1})   = structuredList.(fitnessVectorList{4,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{5,1})   = structuredList.(fitnessVectorList{5,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{6,1})   = structuredList.(fitnessVectorList{6,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{7,1})   = structuredList.(fitnessVectorList{7,1});
        case ismember(objectiveFunctionCase,'F')
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{4,1})   = structuredList.(fitnessVectorList{4,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{5,1})   = structuredList.(fitnessVectorList{5,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{6,1})   = structuredList.(fitnessVectorList{6,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{7,1})   = structuredList.(fitnessVectorList{7,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{8,1})   = structuredList.(fitnessVectorList{8,1});
        case ismember(objectiveFunctionCase,convertStringsToChars(string(strcat({'G'},{'H'},{'I'}))))
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{4,1})   = structuredList.(fitnessVectorList{4,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{5,1})   = structuredList.(fitnessVectorList{5,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{6,1})   = structuredList.(fitnessVectorList{6,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{7,1})   = structuredList.(fitnessVectorList{7,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{8,1})   = structuredList.(fitnessVectorList{8,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{9,1})   = structuredList.(fitnessVectorList{9,1});
        case ismember(objectiveFunctionCase,'J')
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{10,1})   = structuredList.(fitnessVectorList{10,1});
        case ismember(objectiveFunctionCase,'K')
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{10,1})   = structuredList.(fitnessVectorList{10,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{8,1})   = structuredList.(fitnessVectorList{8,1});
        case ismember(objectiveFunctionCase,convertStringsToChars(string(strcat({'L'},{'M'},{'N'}))))
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{10,1})   = structuredList.(fitnessVectorList{10,1});
            evolutions(k).population(i).fitnessVector.(phase{p}).(fitnessVectorList{11,1})   = structuredList.(fitnessVectorList{11,1});
    end
end



evolutions(k).population(i).fitnessVector.Cost                             = nan;
evolutions(k).population(i).fitnessVector.DominationSet                    = nan;
evolutions(k).population(i).fitnessVector.DominatedCount                   = nan;
evolutions(k).population(i).fitnessVector.Rank                             = nan;
evolutions(k).population(i).fitnessVector.RankBySortedVector               = nan;
evolutions(k).population(i).fitnessVector.RankBySortedMagnitude            = nan;
evolutions(k).population(i).fitnessVector.commonW.magnitude                 = nan;
evolutions(k).population(i).fitnessVector.commonW.ranking                   = nan;
evolutions(k).population(i).fitnessVector.commonW.scalingFactor             = nan;
evolutions(k).population(i).fitnessVector.commonW.magnitudeAverages         = nan;




end

