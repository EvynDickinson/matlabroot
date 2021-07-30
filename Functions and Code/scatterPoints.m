
function x = scatterPoints(ROI, nSamples)
% x = scatterPoints(ROI,nSamples)
% ROI =  [minVal,maxVal];
% Randomly draw nSamples numbers between the minVal and maxVal
% to use for x points in a scatter plot
% 
% ES Dickinson
% University of Washington, 2020

% nSamples= 45;


minVal = ROI(1);
maxVal = ROI(2);

rangeVal = abs(diff([minVal,maxVal]));
a = minVal:rangeVal/(nSamples-1):maxVal;
loc = randperm(nSamples);

x = a(loc);

end