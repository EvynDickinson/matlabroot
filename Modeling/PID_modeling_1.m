
clear 
close all
clc

%% PID Modeling Data

%% Load from parquet style flattened data: 

dataDir = 'D:\Evyn Lab Data\Data structures\PID Modeling Data\';
dataFile = [dataDir 'Berlin Caviar all.parquet'];

tic
T = parquetread(dataFile);
toc

%% Extract meta information from the dataset
figDir = createFolder([dataDir, 'Figures\']);
blkbgd = true; % black background figures?
foreColor = formattingColors(blkbgd);

% Get unique trial IDs in order of appearance
[~, first_idx] = unique(T.TrialID, 'stable'); % start index for each trial
% Sort to get them in order
trial_ids = T.TrialID(sort(first_idx));
nTrials = length(trial_ids);

% Get unique temp
protocols = unique(T.TempProtocol, 'stable');
nProtocols = length(protocols);

clear first_idx
initial_vars = who;
initial_vars = add_var(initial_vars, 'initial_vars');

%% FIGURE: quick overview figure of the full set of data
kolor = Color('vaporwavepurple');
LW = 0.5;

r = 4;
c = 1;

fig = getfig;
    subplot(r,c,1)
    plot(T.Temperature, 'color', kolor, 'LineWidth',LW)
    ylabel('temp \circC')
    subplot(r,c,2)
    plot(T.InnerFoodQuad, 'color', kolor, 'LineWidth',LW)
    ylabel('food')
    subplot(r,c,3)
    plot(T.EscapeRing, 'color', kolor, 'LineWidth',LW)
    ylabel('escape')
    subplot(r,c,4)
    plot(T.Sleep, 'color', kolor, 'LineWidth',LW)
    ylabel('sleep')
formatFig(fig, blkbgd, [r,c]);
for ii = 1:r
    subplot(r,c,ii)
    set(gca, XColor='none')
end

clear r c kolor LW

%% Create temperature transformations that we think would be likely to generate the response data: 
% clearvars('-except',initial_vars{:})

% =========================================================
%  STEP 1 — Precompute transforms --
%  The timecourse data for a given temperature protocol is already
%  frame-matched across the trials (since these are being pulled from the
%  post QuadStep4_2.m data organization 
%  This means that we only need to load up one representative sample for
%  each experiment type to generate temperature transforms and then can
%  load it back into the full dataset to make the full temp-behavior table
% 
%  The temperature signal x is identical across all trials, so we compute
%  the I and D transforms once and reuse them in every CV fold.
%  Each column of X_I / X_D corresponds to one window size in t_list.
%  First t-1 rows of each column are NaN (window not yet full).
% ==========================================================

% pull out experiment indexes
initial_vars = add_var(initial_vars, 'protocol_info');
initial_vars = add_var(initial_vars, 'protocol_groups');

% find the start and stop points for each of the trials -- (confirm that
% all within a temp type are consistent) 

% Find start/stop indices for each trial
protocol_info = struct();
for trial = 1:nTrials
    idx = find(strcmp(T.TrialID, trial_ids{trial}));
    protocol_info(trial).trial_id  = trial_ids{trial};
    protocol_info(trial).trial_idx = [idx(1), idx(end)];
    protocol_info(trial).protocol = T.TempProtocol{idx(1)};
end

protocol_groups = struct();
for p = 1:nProtocols
    proto_mask = strcmp({protocol_info.protocol}, protocols{p});
    
    protocol_groups(p).name  = protocols{p};
    protocol_groups(p).trial_ids  = {protocol_info(proto_mask).trial_id};
    protocol_groups(p).trial_idxs = vertcat(protocol_info(proto_mask).trial_idx); % start and stop idx of each protocol
end

%% Data sanity check
% check that all the trials within a temp protocol have the same number of
% data points and overlapping temperature protocols: 
clearvars('-except',initial_vars{:})


% Check that all trials within each temp protocol have the same length
for p = 1:nProtocols
    lengths = diff(protocol_groups(p).trial_idxs, 1, 2) + 1;  % +1 for inclusive length
    
    if ~all(lengths == lengths(1))
        warning('Protocol "%s": trials have unequal lengths [%s] — results may be unreliable.', ...
            protocol_groups(p).name, num2str(lengths'));
            bad = find(lengths ~= lengths(1));
            warning('... offending trials: %s', num2str(bad'));
    else
        fprintf('Protocol "%s": %d trials x %d samples — OK\n', ...
            protocol_groups(p).name, length(lengths), lengths(1));
    end
end

% Visual alignment off all the temperature trial data: 
[r,c] = subplot_numbers(nProtocols);
fig = getfig;
for p = 1:nProtocols
    subplot(r,c,p); hold on
    % plot the time data for each protocol within this protocol type
    % extracted from the T table

end

% Visual alignment of all the temperature trial data:
[r,c] = subplot_numbers(nProtocols);

fig = getfig('',0,[2419 1225]);
for p = 1:nProtocols
    subplot(r,c,p); hold on
    
    idxs = protocol_groups(p).trial_idxs;
    for trial = 1:size(idxs, 1)
        roi = idxs(trial,1) : idxs(trial,2);
        plot(T.Time(roi), T.Temperature(roi), 'Color', Color('Vaporwaveblue'), 'LineWidth',2);
    end
    
    title(protocol_groups(p).name, 'Interpreter', 'none');
    xlabel('time (min)'); ylabel('temp (\circC)');
end
formatFig(fig, blkbgd, [r,c]);

for p = 1:nProtocols
    subplot(r,c,p)
    ylim([13, 37])
    xlim([0, 1000])
end

save_figure(fig, [figDir, 'All temp protocols raw']);

%% Extract temperature data for each exp type
clearvars('-except',initial_vars{:})

% pull the temperature protocol for each experiment type
for p = 1:nProtocols
    roi = protocol_groups(p).trial_idxs(1,1):protocol_groups(p).trial_idxs(1,2);
    protocol_groups(p).temp = T.Temperature(roi);
end

%% Generate temperature transformations for the temperature data

% plot the behavior data vs temperature data to get an idea of the general
% pattern that we are trying to emulate:

% plot out an example of the food attraction behavior for each type of
% trial based on temperature bin (0.5C bins)

sSpan = 30; % 10 second smoothing bin for temperature direction
tempDir = createFolder([figDir, 'temp protocols\']);


% Temperature bins at 0.5C intervals
bin_edges = floor(min(T.Temperature)*2)/2 : 0.5 : ceil(max(T.Temperature)*2)/2; % snap to nearest 0.5
bin_centers = bin_edges(1:end-1) + 0.25;

for p = 1: nProtocols
    temp = protocol_groups(p).temp;
    temp = smooth(temp, sSpan, 'moving'); % smoothed temperature timecourse
    % Temperature direction (logical)
    dT = diff(temp);
    is_increasing = [false; dT > 0]; % shift to match original length
    is_decreasing = [false; dT < 0];

    % bin by temp interval
    bin_idx = discretize(T.Temperature, bin_edges);
    protocol_groups(p).is_increasing = is_increasing;
    protocol_groups(p).is_decreasing = is_decreasing;
    protocol_groups(p).temp_bin = bin_idx;


    % Plot out the temperature heating vs cooling stats: 
    x = 1:length(temp); 
    sz = 35; r = 3; c = 1;
    fig = getfig; 
    subplot(r,c,1); hold on
    scatter(x(is_decreasing), temp(is_decreasing), sz, 'b');
    scatter(x(is_increasing), temp(is_increasing),sz, 'r');
    title(protocol_groups(p).name, 'Interpreter','none')
    subplot(r,c,2)
    scatter(x(is_decreasing), temp(is_decreasing), sz, 'b');
    subplot(r,c,3)
    scatter(x(is_increasing), temp(is_increasing), sz, 'r');
    formatFig(fig, blkbgd, [r,c]);
    for ii = 1:r
        subplot(r,c,ii)
        set(gca, XColor='none')
    end
    save_figure(fig, [tempDir, protocol_groups(p).name, ' heating vs cooling']);

end


    








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





































