

% FILE: INP599_modelgeneration.m
%   
% DESCRIPTION:  
%               
%
% REQUIREMENTS: 
% BUGS: ---
% NOTES:
%
% AUTHOR: Evyn Dickinson | esd39 
% COMPANY: Yale University
% VERSION: 0.1.1
% CREATED: 04-25-22

clear; clc

% set path to directory housing the program code:
dirName = 'INP599_ESD39_finalproject';
if strcmp(getenv('COMPUTERNAME'),'EVYNPC') | strcmp(getenv('COMPUTERNAME'),'ACADIA')
    directoryPath = ['G:\My Drive\INP 599 Stats and Data Analysis\' dirName];
else
    directoryPath = uigetdir('',['Select directory containing''' dirName '']);
end
addpath(directoryPath);
addpath([directoryPath '\Formating Scripts'])
figPath = [directoryPath '\figures\'];


%% Load & Partition Preformatted Data
rng(1234) % for reproduceability

% Load data
load([directoryPath '\data']);

% clean and split into predictors and response variables
regressionData(any(isnan(regressionData),2),:) = []; % remove data that has nans in any of the variables
x = regressionData(:,1:length(predictors)-1);
y = regressionData(:,length(predictors));

% Partition data
l = size(x,1); % length of data
fullIdx = randi(l,[l,1]); % random interation of numbers
trainIdx = fullIdx(1:round(0.6*l)); % training data index
testIdx = fullIdx(round(0.6*l)+1:end); % test data index

trainX = x(trainIdx,:);
trainY = y(trainIdx);
testX = x(testIdx,:);
testY = y(testIdx);
 
initial_vars = who; initial_vars{end+1} = 'initial_vars'; % names of variables to 'restart'
initial_vars{end+1} = 'models';

%% FIGURE: summary example of the data

fig = figure; hold on
    % test data
    [warming, cooling, tempList] = heatingvscooling(testX,testY);
    plot(tempList,warming,'color', Color('red'),'linewidth', 2)
    plot(tempList,cooling,'color', Color('dodgerblue'),'linewidth', 2)
    % train data
    [warming, cooling, tempList] = heatingvscooling(trainX,trainY);
    plot(tempList,warming,'color', Color('darkred'),'linewidth', 2)
    plot(tempList,cooling,'color', Color('blue'),'linewidth', 2)

    xlabel('temperature (\circC)')
    ylabel('distance from food (mm)')
    title('raw data')
    formatFig(fig);
    legend({'test warming', 'test cooling','train warming', 'train cooling'},'box', 'off')

save_figure(fig, [figPath 'Raw data test train division'], '-png');


%% Model 1 
% simple GLM with no interactions or regularization
clearvars('-except',initial_vars{:})

mdlName = 'simple GLM';

mdl = fitglm(testX, testY,'VarNames', predictors);  % train model
y_hat = predict(mdl,testX);                         % predict on withheld test data
mse = sum((testY-y_hat).^2)/length(testY);          % model fit measure (mse)

% store data for later comparison
models(1).name = mdlName;
models(1).mdl = mdl;
models(1).y_hat = y_hat;
models(1).mse = mse;


% Plot example of model predictions:
fig = figure; hold on
    % test data
    [warming, cooling, tempList] = heatingvscooling(testX,testY);
    plot(tempList,warming,'color', Color('red'),'linewidth', 2)
    plot(tempList,cooling,'color', Color('dodgerblue'),'linewidth', 2)
    % prediction data
    [warming, cooling, tempList] = heatingvscooling(testX,y_hat);
    plot(tempList,warming,'color', Color('darkred'),'linewidth', 2)
    plot(tempList,cooling,'color', Color('blue'),'linewidth', 2)

    xlabel('temperature (\circC)')
    ylabel('distance from food (mm)')
    title(mdlName)
    formatFig(fig);
    legend({'test warming', 'test cooling','model warming', 'model cooling'},'box', 'off')

save_figure(fig, [figPath mdlName ' preditions'], '-png');


%% Model 2
% simple GLM with no interactions or regularization
clearvars('-except',initial_vars{:})

mdlName = 'GLM with interactions';

mdl = fitglm(testX, testY,'interactions','VarNames', predictors); % train model
y_hat = predict(mdl,testX);                         % predict on withheld test data
mse = sum((testY-y_hat).^2)/length(testY);          % model fit measure (mse)

% store data for later comparison
models(2).name = mdlName;
models(2).mdl = mdl;
models(2).y_hat = y_hat;
models(2).mse = mse;

% Plot example of model predictions:
fig = figure; hold on
    % test data
    [warming, cooling, tempList] = heatingvscooling(testX,testY);
    plot(tempList,warming,'color', Color('red'),'linewidth', 2)
    plot(tempList,cooling,'color', Color('dodgerblue'),'linewidth', 2)
    % prediction data
    [warming, cooling, tempList] = heatingvscooling(testX,y_hat);
    plot(tempList,warming,'color', Color('darkred'),'linewidth', 2)
    plot(tempList,cooling,'color', Color('blue'),'linewidth', 2)

    xlabel('temperature (\circC)')
    ylabel('distance from food (mm)')
    title(mdlName)
    formatFig(fig);
    legend({'test warming', 'test cooling','model warming', 'model cooling'},'box', 'off')

save_figure(fig, [figPath mdlName ' preditions'], '-png');







