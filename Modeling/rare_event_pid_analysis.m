%% rare_event_pid_analysis.m
%
% Logistic regression comparison of three temperature signal transforms
% (Proportional, Integrative, Derivative) as predictors of a rare binary event.
%
% Framework:
%   - Sweeps over t_num window sizes from 1 sec to max_t minutes
%   - LOO cross-validation; log-loss on held-out trial selects best window
%     for I and D per fold (replaces BIC — see Step 2 comments for rationale)
%   - Full-data model fit using modal best window across folds
%   - Evaluation: AUC-ROC, Precision-Recall, Log-loss, Odds Ratios
%
% Window selection rationale (important):
%   BIC penalizes model parameters, not feature construction complexity.
%   Since every window size produces the same 2-parameter model, BIC cannot
%   distinguish a short noisy window from a long smooth one — it will always
%   prefer the smoother predictor regardless of biological relevance.
%   Log-loss computed on the held-out test trial measures true out-of-sample
%   generalization, so it penalizes windows that are only benefiting from
%   noise smoothing rather than capturing real signal.
%
% To use with real data: replace the SYNTHETIC DATA section with your own
% data loading, ensuring x_data (n_samples x 1) and y_data (n_samples x nTrials)
% are defined, along with num.fps and data.time.
%
% Local copies of computeIntegral / computeDerivative are defined at the
% bottom of this file. Remove them if the originals are already on your path
% and replace the *Local suffix calls above with the original function names.
%
% Requires: Statistics and Machine Learning Toolbox
% ES DICKINSON, YALE, 2026

% clear; clc; close all;
% rng(42);   % reproducibility

%% =========================================================
%  ORGANIZE DATA
%  Required variables:
%    x_data    — n_samples x 1, smoothed temperature signal
%    y_data    — n_samples x nTrials, binary event matrix (~1% rate)
%    num.fps   — scalar, frames per second
%    data.time — n_samples x 1, time vector in minutes
% ==========================================================
clearvars('-except',initial_var{:})

saveDir = createFolder([figDir, 'PID Model/']);

dt = 1/num.fps; % sampling frequency
t_num = 10; % number of samples across time span to check for fit
max_t = 60; % max duration in the past to sample for window time (min)
t_list = round(linspace(num.fps, max_t*num.fps*60, t_num)); % from 1 second to 60 minutes
t_list_time = t_list/num.fps/60; % time windows in minutes 

fprintf('Window sizes to sweep (min): ');
fprintf('%.1f  ', t_list / num.fps / 60);
fprintf('\n');

best = struct();
results = struct();

data_type = 'OutterRing';

nTrials = num.trials*2; % treat each fly independently 
fields = {'P', 'I', 'D'}; % model types

% Pull the data together for this data type
x_data = smooth(data.temperature(:,1), 15, 'moving')-25;
y_data = [squeeze(data.(data_type)(:,M,:)), squeeze(data.(data_type)(:,F,:))]; 
trial_idx = (1:nTrials) .* ones(size(y_data)); % identifying number for each data point

n_samples = size(x_data,1);

fprintf('Data event rate: %.2f%%\n', 100 * mean(y_data(:)));


% num.fps    = 30;          % 1 fps for fast demo; set to real fps (e.g. 30) for actual data
% num.trials = 20;         % trials per group (total trials = num.trials * 2)
% nTrials    = num.trials * 2;
% dt         = 1 / num.fps;
% n_samples  = 630000;
% data.time  = (0:n_samples-1)' / num.fps / 60;   % time axis in minutes
% 
% % Temperature: slow drift + two oscillation frequencies + noise,
% % then smoothed and zero-centered to mimic a real experimental signal
% t_vec    = (1:n_samples)';
% temp_raw = 2   * sin(2*pi*t_vec / 800) + ...
%            1   * sin(2*pi*t_vec / 200) + ...
%            0.5 * randn(n_samples, 1);
% x_data   = smooth(temp_raw, 15, 'moving') - mean(temp_raw);
% 
% % Generate rare binary events driven by the integrative signal at a
% % mid-range window. The proportional term adds a small secondary contribution.
% % This ground truth (I wins) lets us verify the analysis recovers the right answer.
% % Change log_odds to test scenarios where P or D is the true driver instead.
% win_true             = round(n_samples / 10);
% x_I_true             = computeIntegralLocal(x_data, win_true, dt);
% x_I_true(isnan(x_I_true)) = 0;
% 
% % Normalize predictors to zero mean, unit variance before computing log-odds.
% % x_I_true is an unbounded integral that grows large, so using it raw causes
% % the intercept to be overwhelmed and the event rate to blow up.
% % With normalized inputs, the intercept directly sets the baseline log-odds:
% %   intercept = log(p / (1-p)), so -4.6 targets ~1% event rate.
% x_I_norm = (x_I_true - mean(x_I_true)) / std(x_I_true);
% x_P_norm = (x_data   - mean(x_data))   / std(x_data);
% log_odds             = -4.6 + 0.5 * x_I_norm + 0.2 * x_P_norm;
% p_event              = 1 ./ (1 + exp(-log_odds));
% y_data               = double(rand(n_samples, nTrials) < repmat(p_event, 1, nTrials));

% fprintf('Synthetic event rate: %.2f%%\n', 100 * mean(y_data(:)));

% saveDir = 'D:\Evyn Lab Data\PID Modeling\';      % output folder path for saving figures, e.g. './PID Model/'
% blkbgd  = true;   % set true if using dark background figure formatting

% % =========================================================
% %  PID ANALYSIS SETUP
% % ==========================================================
% 
% t_num = 10;        % number of window sizes to evaluate across the sweep
% max_t = 60;        % maximum window duration to consider (minutes)
% min_t = 
% 
% % Window sizes in samples, log-spaced from 1 second up to max_t minutes.
% % Spaced evenly on the sample axis — adjust to linspace or logspace depending
% % on whether you expect the relevant timescale to be in a particular range.
% t_list = round(linspace(num.fps, max_t * num.fps * 60, t_num));
% 
% fields  = {'P', 'I', 'D'};   % model types, used to loop over results struct
% CV      = struct();
% results = struct();

%% =========================================================
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
    X_I(:, ti) = computeIntegralLocal(x_data,   t_list(ti), dt);
    X_D(:, ti) = computeDerivativeLocal(x_data, t_list(ti), dt);
end
x_P = x_data;   % proportional: raw signal, no windowing needed
fprintf('\nTemperature transformations done\n')
toc


%% =========================================================
%  FIGURES — Visualize transforms across window sizes
%
%  Use this to sanity-check the transforms before modeling.
%  Things to look for:
%    I: longer windows should be smoother and track slow temperature trends
%    D: longer windows should suppress fast noise; check for runaway values
%       at short windows (noisy derivative). ylim is clipped to [1st, 99th]
%       percentile to prevent a few outliers from collapsing the y-axis.
% ==========================================================

cmap = parula(t_num);   % replace with Color('dodgerblue','whitesmoke',t_num) if available
r = 1; c = 3;
fig_signals = figure('Name','PID Signal Inputs','Position',[100 100 1067 528]);
subplot(r,c,1); hold on
    plot(data.time, x_P, 'Color', cmap(round(t_num/2),:), 'LineWidth', 2);
    xlabel('time (min)'); title('Proportional (P)'); grid on;

subplot(r,c,2); hold on
    for ii = 1:t_num
        plot(data.time, X_I(:,ii), 'Color', cmap(ii,:), 'LineWidth', 1.5);
    end
    xlabel('time (min)'); title('Integrative (I)'); grid on;

subplot(r,c,3); hold on
    for ii = 1:t_num
        plot(data.time, X_D(:,ii), 'Color', cmap(ii,:), 'LineWidth', 1.5);
    end
    xlabel('time (min)'); title('Derivative (D)');
    ylim(prctile(X_D(~isnan(X_D)), [1 99])); grid on;

sgtitle('PID Signal Inputs Across Window Sizes', Color=foreColor); % main figure title
cb = colorbar(subplot(r,c,3), 'southoutside');
colormap(cmap);
cb.Ticks      = [0 1];
cb.TickLabels = {sprintf('%.1f min', t_list_time(1)), ...
                            sprintf('%.1f min', t_list_time(end))};
cb.Color = foreColor;
cb.Label.String = 'Window size';
cb.Label.Color = foreColor;
formatFig(fig_signals, blkbgd, [r,c]);

save_figure(fig_signals, [saveDir, 'temperature transformations']);

%% =========================================================
%  STEP 2 — LOO cross-validation with log-loss window selection
%
%  Structure: hold out one trial at a time. Train on all remaining trials.
%  For I and D, evaluate every window size on the HELD-OUT trial using
%  log-loss, then select the window with the lowest test log-loss.
%
%  WHY LOG-LOSS INSTEAD OF BIC:
%    BIC is a training-set metric that penalizes the number of fitted
%    parameters. Since every window size produces the same 2-parameter
%    model (intercept + slope), BIC cannot penalize window complexity at
%    all — it will always favor longer, smoother windows that reduce noise
%    regardless of whether they carry real signal. Log-loss computed on
%    the held-out test trial measures genuine out-of-sample generalization,
%    so a longer window only wins if it actually predicts the test events
%    better, not just because it is smoother.
%
%  WHY LOG-LOSS INSTEAD OF AUC FOR SELECTION:
%    With ~1% event rates, a single held-out trial may contain only a
%    handful of positive events. AUC computed on 5-10 positives is very
%    noisy — small fluctuations in those few events can swing the result
%    substantially. Log-loss uses every timepoint in the test trial (both
%    0s and 1s), so it produces a much more stable estimate per fold.
%    AUC and average precision are still used for the final reported
%    performance in Step 4, where predictions are pooled across all folds.
%
%  FOR REAL DATA: swap `for` → `parfor` and call parpool() first.
%  parfor cannot slice into structs, so we extract to plain *_const arrays.
% ==========================================================

% initialize the parallel pool for analysis
if isempty(gcp('nocreate'))
    parpool()
end

% Pre-allocate output cells and arrays
P_y_pred   = cell(nTrials, 1);
I_y_pred   = cell(nTrials, 1);
D_y_pred   = cell(nTrials, 1);
I_t_chosen = nan(nTrials, 1);   % best window (in samples) selected per fold
D_t_chosen = nan(nTrials, 1);

n_train = nTrials - 1;

% Plain-variable copies for parfor broadcast
x_P_const    = x_P;
X_I_const    = X_I;
X_D_const    = X_D;
t_list_const = t_list;
t_num_const  = t_num;

% Log-loss helper: clips predictions away from 0/1 to prevent log(0).
% Operates on column vectors of true labels and predicted probabilities.
logloss_fn = @(y_true, y_pred) ...
    -mean(y_true .* log(y_pred + eps) + (1 - y_true) .* log(1 - y_pred + eps));

fprintf('Running LOO CV (%d folds)...\n', nTrials);
tic

parfor held_out = 1:nTrials   % <-- swap to parfor for large datasets

    % Split into training and test sets.
    % Training: all trials except the held-out one, stacked into a long vector.
    % Test: the single held-out trial (n_samples x 1).
    train_trials = setdiff(1:nTrials, held_out);
    y_train      = y_data(:, train_trials);
    y_train      = y_train(:);          % n_samples*(nTrials-1) x 1
    y_test       = y_data(:, held_out); % n_samples x 1

    % Replicate the (single) x signal to match the stacked y_train length.
    % x is the same across all trials, so we just tile it nTrials-1 times.
    x_P_train = repmat(x_P_const, n_train, 1);
    x_P_test  = x_P_const;

    % ---------------------------------------------------------
    % Proportional — no window to select, fit once per fold
    % ---------------------------------------------------------
    mdl = fitglm(x_P_train, y_train, 'Distribution','binomial','Link','logit');
    P_y_pred{held_out} = predict(mdl, x_P_test);

    % ---------------------------------------------------------
    % Integral — log-loss window selection on held-out trial
    % ---------------------------------------------------------
    best_ll_I = Inf;
    best_ti_I = 1;

    for ti = 1:t_num_const
        x_I_col   = X_I_const(:, ti);
        x_I_train = repmat(x_I_col, n_train, 1);
        x_I_test  = x_I_col;

        % NaN mask: the first (t-1) rows of each window are NaN.
        % We drop NaNs from both train and test before fitting/predicting.
        valid_tr = ~isnan(x_I_train);
        valid_te = ~isnan(x_I_test);

        if sum(valid_tr) < 10, continue; end   % skip if window leaves too few training points

        % Fit on training data
        mdl_ti = fitglm(x_I_train(valid_tr), y_train(valid_tr), ...
                        'Distribution','binomial','Link','logit');

        % Predict on held-out test trial and score with log-loss.
        % We only score on the valid (non-NaN) portion of the test trial.
        y_pred_te        = nan(length(y_test), 1);
        y_pred_te(valid_te) = predict(mdl_ti, x_I_test(valid_te));
        valid_score      = ~isnan(y_pred_te);
        if sum(valid_score) < 10, continue; end

        ll = logloss_fn(y_test(valid_score), y_pred_te(valid_score));

        if ll < best_ll_I
            best_ll_I = ll;
            best_ti_I = ti;
        end
    end
    I_t_chosen(held_out) = t_list_const(best_ti_I);

    % Refit using the selected window and store test predictions for Step 4
    x_I_train = repmat(X_I_const(:, best_ti_I), n_train, 1);
    x_I_test  = X_I_const(:, best_ti_I);
    valid_tr  = ~isnan(x_I_train);
    valid_te  = ~isnan(x_I_test);
    mdl = fitglm(x_I_train(valid_tr), y_train(valid_tr), ...
                 'Distribution','binomial','Link','logit');
    y_pred_I           = nan(length(y_test), 1);
    y_pred_I(valid_te) = predict(mdl, x_I_test(valid_te));
    I_y_pred{held_out} = y_pred_I;

    % ---------------------------------------------------------
    % Derivative — log-loss window selection on held-out trial
    % ---------------------------------------------------------
    best_ll_D = Inf;
    best_ti_D = 1;

    for ti = 1:t_num_const
        x_D_col   = X_D_const(:, ti);
        x_D_train = repmat(x_D_col, n_train, 1);
        x_D_test  = x_D_col;

        valid_tr = ~isnan(x_D_train);
        valid_te = ~isnan(x_D_test);

        if sum(valid_tr) < 10, continue; end

        mdl_ti = fitglm(x_D_train(valid_tr), y_train(valid_tr), ...
                        'Distribution','binomial','Link','logit');

        y_pred_te           = nan(length(y_test), 1);
        y_pred_te(valid_te) = predict(mdl_ti, x_D_test(valid_te));
        valid_score         = ~isnan(y_pred_te);
        if sum(valid_score) < 10, continue; end

        ll = logloss_fn(y_test(valid_score), y_pred_te(valid_score));

        if ll < best_ll_D
            best_ll_D = ll;
            best_ti_D = ti;
        end
    end
    D_t_chosen(held_out) = t_list_const(best_ti_D);

    % Refit using the selected window and store test predictions for Step 4
    x_D_train = repmat(X_D_const(:, best_ti_D), n_train, 1);
    x_D_test  = X_D_const(:, best_ti_D);
    valid_tr  = ~isnan(x_D_train);
    valid_te  = ~isnan(x_D_test);
    mdl = fitglm(x_D_train(valid_tr), y_train(valid_tr), ...
                 'Distribution','binomial','Link','logit');
    y_pred_D           = nan(length(y_test), 1);
    y_pred_D(valid_te) = predict(mdl, x_D_test(valid_te));
    D_y_pred{held_out} = y_pred_D;

    fprintf('  Fold %2d / %d complete\n', held_out, nTrials);
end
toc

% Collect per-fold predictions and window choices into CV struct
CV.P.y_pred   = P_y_pred;
CV.I.y_pred   = I_y_pred;
CV.D.y_pred   = D_y_pred;
CV.I.t_chosen = I_t_chosen;   % in samples; divide by num.fps/60 for minutes
CV.D.t_chosen = D_t_chosen;

%% =========================================================
%  STEP 3 — Full-data fit using modal best window
%
%  After CV, we have a distribution of selected windows across folds.
%  The modal window is used for the final interpretable model fit on
%  all data combined. This is the window we report in figures and tables.
%
%  The distribution of t_chosen is itself informative:
%    - Tight distribution (low variance): data strongly support one timescale
%    - Spread distribution (high variance): the optimal timescale is
%      inconsistent across trials, suggesting the signal may be weak or
%      that the event is not well-described by a single fixed window
% ==========================================================

t_best_I  = mode(CV.I.t_chosen);
t_best_D  = mode(CV.D.t_chosen);
ti_best_I = find(t_list == t_best_I);
ti_best_D = find(t_list == t_best_D);

x_I_best = X_I(:, ti_best_I);
x_D_best = X_D(:, ti_best_D);

y_all   = y_data(:);
x_P_all = repmat(x_P,      nTrials, 1);
x_I_all = repmat(x_I_best, nTrials, 1);
x_D_all = repmat(x_D_best, nTrials, 1);

valid_I = ~isnan(x_I_all);
valid_D = ~isnan(x_D_all);

% Fit full-data models on raw (unscaled) predictors.
% These are used for BIC comparison and as the base for the z-scored OR models below.
results.P.mdl = fitglm(x_P_all,          y_all,          'Distribution','binomial','Link','logit');
results.I.mdl = fitglm(x_I_all(valid_I), y_all(valid_I), 'Distribution','binomial','Link','logit');
results.D.mdl = fitglm(x_D_all(valid_D), y_all(valid_D), 'Distribution','binomial','Link','logit');

results.P.BIC    = results.P.mdl.ModelCriterion.BIC;
results.I.BIC    = results.I.mdl.ModelCriterion.BIC;
results.D.BIC    = results.D.mdl.ModelCriterion.BIC;
results.I.t_best = t_best_I;
results.D.t_best = t_best_D;

fprintf('\nBest window — I: %.1f min | D: %.1f min\n', ...
    t_best_I / num.fps / 60, t_best_D / num.fps / 60);

% --- Window selection stability figure ---
% A tight histogram = consistent optimal timescale across animals/trials.
% A flat or bimodal histogram = timescale is poorly determined by the data.
fig_tstab = figure('Name','Window Selection Stability','Position',[100 540 700 320]);
subplot(1,2,1);
    histogram(CV.I.t_chosen / num.fps / 60, t_num, 'FaceColor','#1f77b4');
    xlabel('Selected window (min)'); title('Integral — window stability');
    ylabel('Folds'); grid on;
subplot(1,2,2);
    histogram(CV.D.t_chosen / num.fps / 60, t_num, 'FaceColor','#2ca02c');
    xlabel('Selected window (min)'); title('Derivative — window stability');
    ylabel('Folds'); grid on;
sgtitle('Log-loss-selected window size across LOO folds');
set(fig_tstab, 'color', 'white')


%% =========================================================
%  STEP 4 — CV performance metrics + evaluation figures
%
%  Predictions from all held-out folds are pooled into a single vector
%  before computing metrics. This gives a global out-of-sample evaluation
%  that uses every trial exactly once as a test case.
%
%  Metrics used and why:
%    AUC-ROC: threshold-free rank discrimination. Asks "does the model
%      score events above non-events?" Good overall summary but can be
%      optimistic when the negative class dominates (common with rare events).
%    Average Precision (PR-AUC): focuses on the event class. More sensitive
%      to how well the model actually finds events vs. how it ranks everything.
%      Preferred for rare event problems — a random classifier scores
%      equal to the event rate (shown as dashed baseline), not 0.5.
%    Log-loss: measures calibration of predicted probabilities. Lower = better.
%      Useful for comparing models but harder to interpret in isolation.
%    BIC: reported for completeness from the full-data models, but note it
%      is a training-set metric and should not be used to rank models here.
% ==========================================================

colors = [0.122 0.471 0.706;   % blue   — P
          0.839 0.373 0.102;   % orange — I
          0.173 0.627 0.173];  % green  — D

for m = 1:numel(fields)
    f = fields{m};

    % Pool held-out predictions across all folds
    y_pred_all = vertcat(CV.(f).y_pred{:});
    y_true_all = y_data(:);

    % Drop any timepoints where prediction is NaN (window warmup period)
    valid      = ~isnan(y_pred_all);
    y_pred_v   = y_pred_all(valid);
    y_true_v   = y_true_all(valid);

    % ROC curve and AUC
    [fpr, tpr, ~, auc_val] = perfcurve(y_true_v, y_pred_v, 1);

    % Precision-Recall curve and Average Precision (PR-AUC)
    [rec, prec, ~, ap_val] = perfcurve(y_true_v, y_pred_v, 1, ...
                                       'XCrit','reca', 'YCrit','prec');

    % Log-loss across pooled predictions
    ll_val = logloss_fn(y_true_v, y_pred_v);

    % Store everything for figures
    CV.(f).AUC     = auc_val;
    CV.(f).AP      = ap_val;
    CV.(f).logloss = ll_val;
    CV.(f).fpr     = fpr;
    CV.(f).tpr     = tpr;
    CV.(f).rec     = rec;
    CV.(f).prec    = prec;

    fprintf('%s — CV AUC: %.3f | Avg Prec: %.3f | Log-loss: %.4f | BIC: %.1f', ...
        f, auc_val, ap_val, ll_val, results.(f).BIC);
    if ismember(f, {'I','D'})
        fprintf(' | t_best: %.1f min', results.(f).t_best / num.fps / 60);
    end
    fprintf('\n');
end

% --- ROC curves ---
fig_roc = figure('Name','ROC Curves','Position',[820 540 540 460]);
    hold on;
    for m = 1:numel(fields)
        f = fields{m};
        plot(CV.(f).fpr, CV.(f).tpr, 'Color', colors(m,:), 'LineWidth', 2.5, ...
            'DisplayName', sprintf('%s  (AUC = %.3f)', f, CV.(f).AUC));
    end
    plot([0 1],[0 1],'k:','LineWidth',1,'HandleVisibility','off');
    xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title('ROC Curves — LOO Cross-Validation');
    legend('Location','southeast'); grid on; xlim([0 1]); ylim([0 1]);

% --- Precision-Recall curves ---
% The dashed baseline equals the event rate — a random classifier sits here.
% Any useful model should be clearly above this line.
fig_pr = figure('Name','Precision-Recall Curves','Position',[1380 540 540 460]);
    hold on;
    for m = 1:numel(fields)
        f = fields{m};
        plot(CV.(f).rec, CV.(f).prec, 'Color', colors(m,:), 'LineWidth', 2.5, ...
            'DisplayName', sprintf('%s  (AP = %.3f)', f, CV.(f).AP));
    end
    yline(mean(y_data(:)), 'k:', 'LineWidth', 1, 'HandleVisibility','off');
    xlabel('Recall'); ylabel('Precision');
    title('Precision-Recall Curves — LOO Cross-Validation');
    legend('Location','northeast'); grid on; xlim([0 1]); ylim([0 1]);

% --- BIC comparison (full-data models) ---
% Note: BIC here is from the full-data models, not CV. It reflects
% in-sample fit quality and should be interpreted with caution — a lower
% BIC does not guarantee better generalization to new trials.
fig_bic = figure('Name','BIC Comparison','Position',[100 80 440 340]);
    bic_vals = [results.P.BIC, results.I.BIC, results.D.BIC];
    b = bar(categorical({'P','I','D'}), bic_vals, 'FaceColor','flat');
    b.CData = colors;
    ylabel('BIC (lower = better in-sample fit)');
    title('Full-Data Model BIC by Transform'); grid on;

% --- CV AUC and Average Precision ---
fig_metrics = figure('Name','CV Metrics','Position',[560 80 700 340]);
subplot(1,2,1);
    b1 = bar(categorical({'P','I','D'}), [CV.P.AUC, CV.I.AUC, CV.D.AUC], 'FaceColor','flat');
    b1.CData = colors;
    ylabel('AUC-ROC'); title('CV AUC-ROC');
    yline(0.5,'k:','LineWidth',1); ylim([0.5 1]); grid on;
subplot(1,2,2);
    b2 = bar(categorical({'P','I','D'}), [CV.P.AP, CV.I.AP, CV.D.AP], 'FaceColor','flat');
    b2.CData = colors;
    ylabel('Average Precision'); title('CV Avg Precision (PR-AUC)');
    yline(mean(y_data(:)),'k:','LineWidth',1); grid on;
sgtitle('Cross-Validated Performance by Transform');

% --- Odds ratios with 95% CI (full-data models, z-scored predictors) ---
%
% WHY Z-SCORE BEFORE FITTING:
%   P is in °C, I is in °C·seconds, D is in °C/second.
%   A raw odds ratio for I reflects the change in odds per 1 °C·second,
%   which is not comparable to 1 °C for P. Z-scoring each predictor
%   (subtract mean, divide by SD) puts all three on the same scale: the
%   OR now reflects the change in odds per 1 standard deviation increase
%   in each transform, making them directly comparable.
%
%   Z-scoring is done NaN-safe (mean/std computed on valid values only)
%   and only affects the OR model — all other metrics above use raw predictors.
zscore_valid = @(v) (v - mean(v(~isnan(v)))) / std(v(~isnan(v)));

x_P_z = zscore_valid(x_P_all);
x_I_z = zscore_valid(x_I_all);
x_D_z = zscore_valid(x_D_all);

results.P.mdl_z = fitglm(x_P_z,          y_all,          'Distribution','binomial','Link','logit');
results.I.mdl_z = fitglm(x_I_z(valid_I), y_all(valid_I), 'Distribution','binomial','Link','logit');
results.D.mdl_z = fitglm(x_D_z(valid_D), y_all(valid_D), 'Distribution','binomial','Link','logit');

fig_or = figure('Name','Odds Ratios (z-scored)','Position',[1280 80 540 320]);
hold on;
plot_order = {'D','I','P'};   % plotted bottom to top
for m = 1:3
    f    = plot_order{m};
    cidx = find(strcmp(fields, f));
    coef = results.(f).mdl_z.Coefficients.Estimate(2);
    ci   = results.(f).mdl_z.coefCI;
    or_v  = exp(coef);
    or_lo = exp(ci(2,1));
    or_hi = exp(ci(2,2));
    errorbar(or_v, m, or_v - or_lo, or_hi - or_v, ...
        'horizontal', 'o', ...
        'Color',           colors(cidx,:), ...
        'MarkerFaceColor', colors(cidx,:), ...
        'MarkerSize', 8, 'LineWidth', 2);
end
xline(1,'k--','LineWidth',1);   % OR = 1 is the null (no effect)
yticks(1:3); yticklabels({'D (Derivative)','I (Integral)','P (Proportional)'});
xlabel('Odds Ratio (per 1 SD increase in predictor, 95% CI)');
title('Odds Ratios — Full-Data Models (z-scored)');
grid on;

fprintf('\nDone.\n');

%% =========================================================
%  LOCAL FUNCTION DEFINITIONS
%
%  Exact copies of computeIntegral / computeDerivative.
%  Remove these and replace *Local calls above with the original
%  function names if they are already on your MATLAB path.
% ==========================================================

function x_I = computeIntegralLocal(x, t, dt)
% Sliding window trapezoidal integral. See computeIntegral.m for full docs.
    kernel          = ones(t, 1);
    kernel([1 end]) = 0.5;              % trapz correction: half-weight endpoints
    x_conv          = conv(x, kernel, 'full') * dt;
    n               = length(x);
    x_I             = [nan(t-1, 1); x_conv(t : t+n-t)];
end

function x_D = computeDerivativeLocal(x, t, dt)
% Sliding window least-squares derivative. See computeDerivative.m for full docs.
    ramp   = (-(t-1)/2 : (t-1)/2)';         % centered linear ramp = OLS slope weights
    w      = ramp / (sum(ramp.^2) * dt);    % normalize to give output in units of x/sec
    x_conv = conv(x, flipud(w), 'full');
    n      = length(x);
    x_D    = [nan(t-1, 1); x_conv(t : t+n-t)];
end
