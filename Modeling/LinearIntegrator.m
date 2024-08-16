
%% Prep and load data
clearvars('-except',initial_vars{:})

% 1) load data from step 4.2 : QuadStep4.2.m

exp = 1;
grouped(exp).name

% pull temperature and avg distance
x = grouped(exp).temp; % stimulus data (temperature time series)
y = grouped(exp).dist.all; % response data (w/in trial avg distance to food) 

% account for nans...
loc = (any(isnan(y),2));
y(loc,:) = [];
x(loc,:) = [];

% 2) organize data for the model
trainRatio = 0.7; % how much of the data to use for training e.g. 70%
trainIndex = round(trainRatio * length(x));

% split data for training and testing
X_train = x(1:trainIndex);
Y_train = y(1:trainIndex, :);

X_test = x(trainIndex+1:end);
Y_test = y(trainIndex+1:end, :);


%% Simple linear integrator model

% Define the model function with time delay (same as before)
LIM = @(params, X) ...
    params(1) * exp(-X / params(2)) - params(3) * exp(-(X - params(5)) / params(4)) .* heaviside(X - params(5));

% Initial parameter guesses [A1, tau1, A2, tau2, delta_t]
initialParams = [1, 1, 1, 1, 100];  % Initial guesses including time delay

% Objective function: sum of squared errors across all trials in the training set
objectiveFunc = @(params) sum(sum((LIM(params, X_train) - Y_train).^2, 2));

% Fit the model using fminsearch (simple unconstrained optimization)
tic
options = optimset('Display', 'iter');  % Show iteration output
fitParams = fminsearch(objectiveFunc, initialParams, options);
toc


objectiveFunc = @(params) ...
    sum(sum((modelFunc(params, X_train) - Y_train).^2, 2));

% For debugging, add a line like this:
fprintf('params: %f %f %f %f %f\n', params);
fprintf('Objective function value: %f\n', objectiveFunc(initialParams));




% Generate predictions for the training set and test set using the fitted parameters
Y_train_pred = LIM(fitParams, X_train);
Y_test_pred = LIM(fitParams, X_test);

% Plot the observed vs. predicted for the training set
figure;
subplot(2,1,1);
plot(X_train, Y_train, 'b', 'DisplayName', 'Observed Training Data');
hold on;
plot(X_train, Y_train_pred, 'r', 'DisplayName', 'Model Prediction');
legend;
xlabel('Time (s)');
ylabel('Response');
title('Training Data: Observed vs. Model Prediction');

% Plot the observed vs. predicted for the test set
subplot(2,1,2);
plot(X_test, Y_test, 'b', 'DisplayName', 'Observed Test Data');
hold on;
plot(X_test, Y_test_pred, 'r', 'DisplayName', 'Model Prediction');
legend;
xlabel('Time (s)');
ylabel('Response');
title('Test Data: Observed vs. Model Prediction');

% Calculate performance metrics for the training and test sets
trainError = sum(sum((Y_train - Y_train_pred).^2));
testError = sum(sum((Y_test - Y_test_pred).^2));

fprintf('Training set error (Sum of squared errors): %f\n', trainError);
fprintf('Test set error (Sum of squared errors): %f\n', testError);

%% Attemped linear regression


% Example data (replace with your actual stimulus and response data)
X = (0:0.01:10)';  % Example stimulus (column vector)
Y = [sin(2 * pi * 0.5 * X) + 0.1 * randn(size(X)), ...
     sin(2 * pi * 0.5 * X) + 0.15 * randn(size(X))];  % Example response (2 trials)

% Number of time lags to consider
numLags = 50;

% Create the time-lagged predictor matrix for each trial
numTrials = size(Y, 2);
X_lagged = zeros(length(X) - numLags, numLags);
for lag = 1:numLags
    X_lagged(:, lag) = X(numLags-lag+1:end-lag);  % Lagged stimulus data
end

% Adjust Y to match the dimensions of X_lagged
Y_adjusted = Y(numLags+1:end, :);

% Flatten Y_adjusted and X_lagged for all trials
Y_flat = Y_adjusted(:);
X_flat = repmat(X_lagged, numTrials, 1);

% Perform linear regression
b = X_flat \ Y_flat;  % Regression coefficients

% Predict the response using the model
Y_pred_flat = X_flat * b;

% Reshape Y_pred_flat back into the original trial structure
Y_pred = reshape(Y_pred_flat, [], numTrials);

% Plot the observed vs. predicted for the first trial
figure;
plot(X(numLags+1:end), Y_adjusted(:, 1), 'b', 'DisplayName', 'Observed Response');
hold on;
plot(X(numLags+1:end), Y_pred(:, 1), 'r--', 'DisplayName', 'Predicted Response');
legend;
xlabel('Time (s)');
ylabel('Response');
title('Observed vs. Predicted Response (Trial 1)');

% Optionally, calculate performance metrics like R-squared for each trial
rSquared = 1 - sum((Y_adjusted - Y_pred).^2) ./ sum((Y_adjusted - mean(Y_adjusted)).^2);
fprintf('R-squared for each trial: %f\n', rSquared);
