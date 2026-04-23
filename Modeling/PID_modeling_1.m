
clear 
close all
clc

%% PID Modeling Data

%% Load from parquet style flattened data: 

dataFile = 'D:\Evyn Lab Data\Data structures\PID Modeling Data\Berlin Caviar all data.parquet';
tic
T = parquetread(dataFile);
toc

%% Extract meta information from the dataset

% need a set of identifiers for each type of experiment
temp_protocols = unique(T.TempProtocol);
trial_IDs = unique(T.TrialID);
nTrials = length(trial_IDs);
nProtocols = length(temp_protocols);


%% Create temperature transformations that we think would be likely to generate the response data: 

% =========================================================
%  STEP 1 — Precompute transforms
% 
%  The temperature signal x is identical across all trials, so we compute
%  the I and D transforms once and reuse them in every CV fold.
%  Each column of X_I / X_D corresponds to one window size in t_list.
%  First t-1 rows of each column are NaN (window not yet full).
% ==========================================================



tic
X_I = nan(n_samples, t_num);
X_D = nan(n_samples, t_num);

for ti = 1:t_num
    X_I(:, ti) = computeIntegral(x_data,   t_list(ti), dt);
    X_D(:, ti) = computeDerivative(x_data, t_list(ti), dt);
end
x_P = x_data;   % proportional: raw signal, no windowing needed
fprintf('\nTemperature transformations done\n')
toc





































