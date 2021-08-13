
function [prob, occCumSum] = getOccProb(occCumSum, newProb, n)
% [prob, occCumSum] = getOccProb(occCumSum, newProb, n);
% 
% Update occupation probability for a data point given the previous
% occupation probabilty and the known number of instances generating the
% current probability
%
% Argin :
% occCumSum --> mat with sum occupancy for each pixel
% newProb --> mat with new occupation values for n
% n --> total number of datapoints thus far
%
% Argout: 
% prob --> mat probability for time point n
% occCumSum --> mat with total occupancy value at n
%
% ES Dickison
% Yale University 
% Created Aug 2021
occCumSum = (occCumSum + newProb);
prob = occCumSum / n;

end