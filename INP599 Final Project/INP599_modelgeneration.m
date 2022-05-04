

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


%% Model 3
% L1 regularization GLM
clearvars('-except',initial_vars{:})

mdlName = 'L1 regularization';

% pull structures off GPU
X = gather(trainX);
Y = gather(trainY); 

% -------------------------------- lamda optimization  ----------------------------------
% Test a set of lamda values to find the best model version
lambda_val = [];
nLambda = 1000;
lambda_val(:,1) = logspace(-4,2,nLambda);

tic
% test various values of lamda on the training data to see the best model fits
[mdl,lassoData] = deal([]);
for ii = 1:nLambda
    % pick a lambda value
    mdl(ii).lambda = lambda_val(ii); % goes up by an order of magnitude each time
    
    % fit linear regression model with LASSO regularization
    [mdl(ii).B,mdl(ii).fitInfo] = lasso(X,Y,'lambda',mdl(ii).lambda);
    mdl(ii).train_MSE = mdl(ii).fitInfo.MSE;
    
    % fit the test data for each lambda model
    data = mdl(ii).B.*testX';
    mdl(ii).y_hat = sum(data,1);
    mdl(ii).test_MSE = sum((testY'-mdl(ii).y_hat).^2)/length(testY);  % model fit measure (mse)

    % plot data
    lassoData.lambda(ii) = mdl(ii).lambda;
    lassoData.train_MSE(ii) = mdl(ii).train_MSE;
    lassoData.test_MSE(ii) = mdl(ii).test_MSE;
    lassoData.B(:,ii) = mdl(ii).B;
end
toc

% 'best' model fit: (aka best mix of bias and variance)
[~,loc] = min(lassoData.test_MSE);
lamda = lambda_val(loc);

% FIGURE: lambda values and prediction MSE 
LW = 1.5;
CList = Color('midnightblue','powderblue',20);
fig = figure; 
hold on
    plot(lassoData.lambda,lassoData.train_MSE,'color', Color('teal'),'LineWidth',LW)
    plot(lassoData.lambda,lassoData.test_MSE,'color', Color('orange'),'LineWidth',LW)
    set(gca,'XScale','log')
    v_line(lamda,'darkgrey',':',2)
    % labels and formatting
    title(mdlName)
    xlabel('\lambda')
    ylabel('MSE')
    set(gca,'fontsize', 15)
    legend({'Train Data', 'Test Data'},'FontSize',13,'Color','k','Box','off','Location','northwest')
formatFig(fig);

save_figure(fig, [figPath mdlName ' lamda parameters'], '-png');
% ---------------------------------------------------------------------------------------

% Get full model with optimized lambda:
y_hat = mdl(loc).y_hat';
mse = gather(mdl(loc).test_MSE);

% store data for later comparison
models(3).name = mdlName;
models(3).mdl = mdl(loc);
models(3).y_hat = y_hat;
models(3).mse = mse;

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

%% Model 4 
% L1 regularization with interactions GLM
clearvars('-except',initial_vars{:})

mdlName = 'L1 with interactions';

% pull structures off GPU
X = gather(trainX);
Y = gather(trainY); 

% build in interactions terms
D = x2fx(X,'interaction');
D(:,1) = []; % No constant term
testD = x2fx(gather(testX),'interaction');
testD(:,1) = []; 


% -------------------------------- lamda optimization  ----------------------------------
% Test a set of lamda values to find the best model version
lambda_val = [];
nLambda = 1000;
lambda_val(:,1) = logspace(-4,2,nLambda);

tic
% test various values of lamda on the training data to see the best model fits
[mdl,lassoData] = deal([]);
for ii = 1:nLambda
    % pick a lambda value
    mdl(ii).lambda = lambda_val(ii); % goes up by an order of magnitude each time
    
    % fit linear regression model with LASSO regularization
    [mdl(ii).B,mdl(ii).fitInfo] = lasso(D,Y,'lambda',mdl(ii).lambda);
    mdl(ii).train_MSE = mdl(ii).fitInfo.MSE;
    
    % fit the test data for each lambda model
    data = mdl(ii).B.*testD';
    mdl(ii).y_hat = sum(data,1);
    mdl(ii).test_MSE = sum((testY'-mdl(ii).y_hat).^2)/length(testY);  % model fit measure (mse)

    % plot data
    lassoData.lambda(ii) = mdl(ii).lambda;
    lassoData.train_MSE(ii) = mdl(ii).train_MSE;
    lassoData.test_MSE(ii) = mdl(ii).test_MSE;
    lassoData.B(:,ii) = mdl(ii).B;
end
toc

% 'best' model fit: (aka best mix of bias and variance)
[~,loc] = min(lassoData.test_MSE);
lamda = lambda_val(loc);

% FIGURE: lambda values and prediction MSE 
LW = 1.5;
CList = Color('midnightblue','powderblue',20);
fig = figure; 
hold on
    plot(lassoData.lambda,lassoData.train_MSE,'color', Color('teal'),'LineWidth',LW)
    plot(lassoData.lambda,lassoData.test_MSE,'color', Color('orange'),'LineWidth',LW)
    set(gca,'XScale','log')
    v_line(lamda,'darkgrey',':',2)
    % labels and formatting
    title(mdlName)
    xlabel('\lambda')
    ylabel('MSE')
    set(gca,'fontsize', 15)
    legend({'Train Data', 'Test Data'},'FontSize',13,'Color','k','Box','off','Location','northwest')
formatFig(fig);

save_figure(fig, [figPath mdlName ' lamda parameters'], '-png');
% ---------------------------------------------------------------------------------------

% Get full model with optimized lambda:
y_hat = mdl(loc).y_hat';
mse = gather(mdl(loc).test_MSE);

% store data for later comparison
models(4).name = mdlName;
models(4).mdl = mdl(loc);
models(4).y_hat = y_hat;
models(4).mse = mse;

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


%% SAVE working data

save('G:\My Drive\INP 599 Stats and Data Analysis\INP599_ESD39_finalproject\model data.mat','models')
save('G:\My Drive\INP 599 Stats and Data Analysis\INP599_ESD39_finalproject\working data.mat','-v7.3');


% LOAD working data:
load('G:\My Drive\INP 599 Stats and Data Analysis\INP599_ESD39_finalproject\working data.mat','-v7.3');


