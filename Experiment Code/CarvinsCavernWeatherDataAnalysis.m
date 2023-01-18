%% Load Data

clear; clc
rootdir = 'G:\My Drive\Jeanne Lab\DATA\Carvins Cove Virginia Meterological Dataset/';
% sheetName = 'CCR_Met_final_2021';
dataPath = [rootdir 'Temperature Data.xlsx'];
data = readmatrix(dataPath,'sheet',sheetName);

%% 
fig = figure;

plot(data(:,4))



ROI = 1:600; % 60 minutes of temperature data
fig = figure;
plot(data(ROI,4))

