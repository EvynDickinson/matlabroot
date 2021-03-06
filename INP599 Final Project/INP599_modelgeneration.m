

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
    % All data
    [warming, cooling, tempList] = heatingvscooling(regressionData(:,1:11),regressionData(:,12));
    plot(tempList,warming,'color', Color('black'),'linewidth', 2)
    plot(tempList,cooling,'color', Color('darkgrey'),'linewidth', 2)
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
    legend({'all data warming', 'all data cooling','test warming', 'test cooling','train warming', 'train cooling'},'box', 'off')

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
nLambda = 500;
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
    if rem(ii,25)==0
        disp(ii)
    end
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


%% Model 5
% L2 regularization GLM
clearvars('-except',initial_vars{:})

mdlName = 'L2 regularization';

% pull structures off GPU
X = gather(trainX);
Y = gather(trainY); 
tX = gather(testX);
tY = gather(testY);

% -------------------------------- lamda optimization  ----------------------------------
% Test a set of lamda values to find the best model version
lambda_val = [];
nLambda = 500;
% lambda_val(:,1) = logspace(-2,12,nLambda);
lambda_val(:,1) = logspace(6,10,nLambda);


tic
% test various values of lamda on the training data to see the best model fits
[mdl,ridgeData] = deal([]);
for ii = 1:nLambda
    % pick a lambda value
    mdl(ii).lambda = lambda_val(ii); % goes up by an order of magnitude each time
    
    % fit linear regression model with LASSO regularization
    mdl(ii).B = ridge(Y,X,mdl(ii).lambda);

    % find Y-train-hat and MSE 
    data = mdl(ii).B.*X';
    mdl(ii).train_y_hat = sum(data,1);
    mdl(ii).train_MSE = (1/length(Y))*sum((Y-mdl(ii).train_y_hat').^2);

    % fit the test data for each lambda model
    data = mdl(ii).B.*tX';
    mdl(ii).y_hat = sum(data,1);
    mdl(ii).MSE = (1/length(tY))*sum((tY-mdl(ii).y_hat').^2);

    % plot data
    ridgeData.lambda(ii) = mdl(ii).lambda;
    ridgeData.train_MSE(ii) = mdl(ii).train_MSE;
    ridgeData.test_MSE(ii) = mdl(ii).MSE;
    ridgeData.B(:,ii) = mdl(ii).B;

    if rem(ii,25)==0
        disp(ii)
    end
end
toc

% 'best' model fit: (aka best mix of bias and variance)
[~,loc] = min(ridgeData.test_MSE);
lamda = lambda_val(loc);

% FIGURE: lambda values and prediction MSE 
LW = 1.5;
CList = Color('midnightblue','powderblue',20);
fig = figure; 
hold on
    plot(ridgeData.lambda,ridgeData.train_MSE,'color', Color('teal'),'LineWidth',LW)
    plot(ridgeData.lambda,ridgeData.test_MSE,'color', Color('orange'),'LineWidth',LW)
    set(gca,'XScale','log')
    set(gca,'YScale', 'log')
    v_line(lamda,'darkgrey',':',2)
    % labels and formatting
    title(mdlName)
    xlabel('\lambda')
    ylabel('MSE')
    set(gca,'fontsize', 15)
    legend({'Train Data', 'Test Data'},'FontSize',13,'Color','k','Box','off','Location','southwest')
formatFig(fig);

save_figure(fig, [figPath mdlName ' lamda parameters'], '-png');
% ---------------------------------------------------------------------------------------

% Get full model with optimized lambda:
y_hat = mdl(loc).y_hat';
mse = gather(mdl(loc).MSE);

% store data for later comparison
models(5).name = mdlName;
models(5).mdl = mdl(loc);
models(5).y_hat = y_hat;
models(5).mse = mse;

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


%% Model 6 
% L2 regularization with interactions GLM
clearvars('-except',initial_vars{:})

mdlName = 'L2 with interactions';

% pull structures off GPU
X = gather(trainX);
Y = gather(trainY); 
tX = gather(testX);
tY = gather(testY);

% build in interactions terms
D = x2fx(X,'interaction');
testD = x2fx(tX,'interaction');


% -------------------------------- lamda optimization  ----------------------------------
% Test a set of lamda values to find the best model version
lambda_val = [];
nLambda = 250;
lambda_val(:,1) = logspace(5,11,nLambda);

tic
% test various values of lamda on the training data to see the best model fits
[mdl,ridgeData] = deal([]);
for ii = 1:nLambda
    % pick a lambda value
    mdl(ii).lambda = lambda_val(ii); % goes up by an order of magnitude each time
    
    % fit linear regression model with RIDGE regularization
    mdl(ii).B = ridge(Y,D,mdl(ii).lambda);

    % find Y-train-hat and MSE 
    data = mdl(ii).B.*D';
    mdl(ii).train_y_hat = sum(data,1);
    mdl(ii).train_MSE = (1/length(Y))*sum((Y-mdl(ii).train_y_hat').^2);

    % fit the test data for each lambda model
    data = mdl(ii).B.*testD';
    mdl(ii).y_hat = sum(data,1);
    mdl(ii).MSE = (1/length(tY))*sum((tY-mdl(ii).y_hat').^2);

    % plot data
    ridgeData.lambda(ii) = mdl(ii).lambda;
    ridgeData.train_MSE(ii) = mdl(ii).train_MSE;
    ridgeData.test_MSE(ii) = mdl(ii).MSE;
    ridgeData.B(:,ii) = mdl(ii).B;
    if rem(ii,25)==0
        disp(ii)
    end
end
toc

% 'best' model fit: (aka best mix of bias and variance)
[~,loc] = min(ridgeData.test_MSE);
lamda = lambda_val(loc);

% FIGURE: lambda values and prediction MSE 
LW = 1.5;
CList = Color('midnightblue','powderblue',20);
fig = figure; 
hold on
    plot(ridgeData.lambda,ridgeData.train_MSE,'color', Color('teal'),'LineWidth',LW)
    plot(ridgeData.lambda,ridgeData.test_MSE,'color', Color('orange'),'LineWidth',LW)
    set(gca,'XScale','log','YScale','log')
    
    v_line(lamda,'darkgrey',':',2)
    % labels and formatting
    title(mdlName)
    xlabel('\lambda')
    ylabel('MSE')
    set(gca,'fontsize', 15)
    legend({'Train Data', 'Test Data'},'FontSize',13,'Color','k','Box','off','Location','southwest')
formatFig(fig);

save_figure(fig, [figPath mdlName ' lamda parameters'], '-png');
% ---------------------------------------------------------------------------------------

% Get full model with optimized lambda:
y_hat = mdl(loc).y_hat';
mse = gather(mdl(loc).MSE);

% store data for later comparison
models(6).name = mdlName;
models(6).mdl = mdl(loc);
models(6).y_hat = y_hat;
models(6).mse = mse;

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

%% Coefficent comparisons
nModel = 6;   % number of models included  
nParams = 67; % largest number of model features

% ---- Pull the coefficients for each model ----
B = struct;
% simple GLM
B(1).coeff = gather(models(1).mdl.Coefficients.Estimate);
B(1).predictors = ['intercept',predictors(1:end-1)];
% simple GLM with interactions
B(2).coeff = gather(models(2).mdl.Coefficients.Estimate);
B(2).predictors = models(2).mdl.CoefficientNames;  
% L1 GLM
B(3).coeff = [nan; models(3).mdl.B];
B(3).predictors = ['intercept',predictors(1:end-1)];
% L1 GLM with interactions
B(4).coeff = [nan; models(4).mdl.B];
B(4).predictors = B(2).predictors; 
% L2 GLM
B(5).coeff = [nan; models(5).mdl.B];
B(5).predictors = ['intercept',predictors(1:end-1)];
% L2 GLM with interactions
B(6).coeff = models(6).mdl.B;
B(6).predictors = B(2).predictors; % doesn't include an intercept

mdlWeight = nan(nParams,nModel);
for ii = 1:nModel
    mdlWeight(1:length(B(ii).coeff),ii) = B(ii).coeff;
end

% Model Color Assignment
CList = {'Indigo','DarkViolet','Orange','Gold','Blue','DodgerBlue'};
for ii = 1:nModel
    B(ii).color = Color(CList{ii});
    mdlNames{ii} = models(ii).name;
    mse(ii,1) = gather(models(ii).mse);
end

tbl = table(mdlNames',mse,'VariableNames', {'Model','MSE'});
fprintf('\n\n Model MSEs\n\n')
disp(tbl)

% FIGURES: plot the coefficients for each model
for ii = 1:nModel
    fig = figure; set(fig,'pos', [1950 764 989 498])
        hold on
        plotY = mdlWeight(:,ii);
        loc = isnan(plotY);
        plotY(loc) = [];
        plotX = 1:nParams;
        plotX(loc) = [];
        
        h = bar(plotX,plotY);
        set(h,'BarWidth',0.8,'FaceColor',B(ii).color, 'EdgeColor','k','FaceAlpha',1)
%         set(gca,'YScale','log')
        xlabel('Feature')
        ylabel('Coefficient (feature weight)')
        title(models(ii).name)
        formatFig(fig);
    
    save_figure(fig, [figPath models(ii).name ' coefficient weights not log'], '-png');
end

%% Plot all the coefficients across the models together


% FIGURE: plot the coefficients compared across models
fig = figure; set(fig, 'pos',[-1003 377 940 915])
subplot(2,1,1)
hold on; roi = 1:2;
    for ii = roi
        scatter(1:nParams,mdlWeight(:,ii),50,B(ii).color,'filled')
    end
    xlabel('Feature')
    ylabel('Coefficient (feature weight)')
    legend(mdlNames(roi),'box','off')
    
subplot(2,1,2); roi = 3:nModel;
    hold on
    for ii = roi
        scatter(1:nParams,mdlWeight(:,ii),50,B(ii).color,'filled')
    end
    xlabel('Feature')
    ylabel('Coefficient (feature weight)')
    legend(mdlNames(roi),'box','off')
    
formatFig(fig, false, [2,1]);
save_figure(fig, [figPath 'All models coefficient weights'], '-png');


% Make tables showing the largest features and what they represent for each model

for ii = 1:2
    % organize the coefficients to find the largest 
    sortedFeatures = [abs(mdlWeight(:,ii)),mdlWeight(:,ii),(1:nParams)'];
    sortedFeatures(isnan(sortedFeatures(:,1)),:) = [];
    % add p-values into the structure
    sortedFeatures(:,4) = gather(models(ii).mdl.Coefficients.pValue);
    Features = sortrows(sortedFeatures,'descend');
    % pull the names of the largest cofficients
    tbl = table(B(2).predictors(Features(:,3))', Features(:,2), Features(:,4),'VariableNames', {'Feature', 'Coefficient','P-Value'});
    fprintf(['\n\nModel: ' mdlNames{ii} '\n\n'])
    disp(tbl)
end


for ii = 3:6
    % organize the coefficients to find the largest 
    sortedFeatures = [abs(mdlWeight(:,ii)),mdlWeight(:,ii),(1:nParams)'];
    sortedFeatures(isnan(sortedFeatures(:,1)),:) = [];
    Features = sortrows(sortedFeatures,'descend');
    % pull the names of the largest cofficients
    tbl = table(B(2).predictors(Features(:,3))', Features(:,2), 'VariableNames', {'Feature', 'Coefficient'});
    fprintf(['\n\nModel: ' mdlNames{ii} '\n\n'])
    disp(tbl)
end
    
    
    
%% Comparision of model MSE

% FIGURE: plot the first twelve feature coefficients across models
buff = linspace(-0.3,0.3,nModel);

fig = figure; set(fig, 'pos',[1960 664 473 601])
    hold on
    for ii = 1:nModel
        plot([ii,ii],[0,models(ii).mse],'linewidth',1.5,'Color',B(ii).color)
        scatter(ii,models(ii).mse,75,B(ii).color,'filled')
    end
    xlim([0,nModel+1])
    ax = gca;
    ax.XTick = 1:nModel;
    ax.XTickLabel = mdlNames;
    ax.XTickLabelRotation = 30;
    ylabel('MSE')
    formatFig(fig);

save_figure(fig, [figPath 'All models MSE'], '-png');



%% SAVE working data

save('G:\My Drive\INP 599 Stats and Data Analysis\INP599_ESD39_finalproject\model data.mat','models');
save('G:\My Drive\INP 599 Stats and Data Analysis\INP599_ESD39_finalproject\working data.mat','-v7.3');


% LOAD working data:
load('G:\My Drive\INP 599 Stats and Data Analysis\INP599_ESD39_finalproject\working data.mat');


