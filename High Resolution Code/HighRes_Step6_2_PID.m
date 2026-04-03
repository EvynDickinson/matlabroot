

%% PID Analysis Setup
clearvars('-except',initial_var{:})

saveDir = createFolder([figDir, 'PID Model/']);

dt = 1/num.fps; % sampling frequency
t_num = 10; % number of samples across time span to check for fit
max_t = 60; % max duration in the past to sample for window time (min)
t_list = round(linspace(num.fps, max_t*num.fps*60, t_num)); % from 1 second to 60 minutes

best = struct();
results = struct();

data_type = 'OutterRing';

nTrials = num.trials*2; % treat each fly independently 
fields = {'P', 'I', 'D'}; % model types

% Pull the data together for this data type
x_data = smooth(data.temperature(:,1), 15, 'moving');
y_data = [squeeze(data.(data_type)(:,M,:)), squeeze(data.(data_type)(:,F,:))]; 
trial_idx = (1:nTrials) .* ones(size(y_data)); % identifying number for each data point

n_samples = size(x_data,1);

%% ANALYSIS: Pull the data together across trials:  

% % Visual check of stability of t (after step 3) 
% figure;
% subplot(1,2,1); histogram(CV.I.t_chosen); xlabel('t chosen (I)'); title('Integral t stability');
% subplot(1,2,2); histogram(CV.D.t_chosen); xlabel('t chosen (D)'); title('Derivative t stability');

% =========================================================
% STEP 1 — precompute transforms ONCE (x is same across trials)
% =========================================================
% note: takes like 1/3 a sec per iteration (e.g. n = 10 takes 3 seconds to run)
X_I = nan(n_samples, t_num);
X_D = nan(n_samples, t_num);
for ti = 1:t_num
    X_I(:, ti) = computeIntegral(x_data,   t_list(ti), dt);
    X_D(:, ti) = computeDerivative(x_data, t_list(ti), dt);
end
x_P = x_data;

%% FIGURES: visualize the different time integrations / derivatives for the different
% t values: 
cList = Color('dodgerblue', 'whitesmoke', t_num);
fig = getfig('', 1, [ 954 676]); 
    subplot(1,3,1); hold on
        ax = gca; ax.ColorOrder = cList; 
        plot(data.time, x_P, 'linewidth', 2)
        title('proportional')
    subplot(1,3,2); hold on
        ax = gca; ax.ColorOrder = cList; 
        plot(data.time, X_I, 'linewidth', 2)
        title('integrative')
    subplot(1,3,3); hold on
        ax = gca; ax.ColorOrder = cList;
        plot(data.time, X_D, 'linewidth', 2)
        title('derivative')
        ylim([-0.02, 0.02])
formatFig(fig, blkbgd, [1,3]);
for ii = 1:3;   subplot(1,3,ii); xlabel('time (min)'); end
save_figure(fig, [saveDir, 'signal inputs']);

% derivative model only: 
[r, c] = subplot_numbers(t_num,5);
fig = figure; 
for ii = 1:t_num
    subplot(r,c,ii)
    plot(data.time, X_D(:,ii), LineWidth=2, Color=cList(ii,:))
end
formatFig(fig, blkbgd, [r,c]);
for ii = 1:t_num; subplot(r,c,ii); set(gca, XColor='none'); end

save_figure(fig, [saveDir, 'derivative signal inputs']);

% visualize the time lag periods
fig = getfig('', 1, [408 468]); hold on
    t_list_time = t_list/num.fps/60;
    plot(t_list_time, Color='m', LineStyle=':')
    for ii = 1:num_t
        scatter(ii, t_list_time(ii), 35, cList(ii,:), 'filled')
    end
    formatFig(fig, blkbgd);
    xlabel('t index')
    ylabel('t window (min)')
    save_figure(fig, [saveDir, 'time lag t']);

%% ANALYSIS SINGLE CHANNEL PROCESSING: 
% =========================================================
% STEP 2 — LOO cross-validation
% =========================================================
CV.P.y_pred   = cell(nTrials, 1);
CV.I.y_pred   = cell(nTrials, 1);
CV.D.y_pred   = cell(nTrials, 1);
CV.I.t_chosen = nan(nTrials, 1);
CV.D.t_chosen = nan(nTrials, 1);

n_train = nTrials - 1;

tic
for held_out = 1:nTrials % approx 4 mins per loop 

    train_trials = setdiff(1:nTrials, held_out);
    y_train = y_data(:,train_trials);
    y_train = y_train(:); % reshape to single vector
    y_test  = y_data(:,held_out);
    
    % same temp data for all the trials for train and test, since it's
    % uniform across trials
    x_P_train = repmat(x_P, n_train, 1);
    x_P_test  = x_P;

    % ---------------------------------------------------------
    % Proportional
    % ---------------------------------------------------------
    mdl = fitglm(x_P_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');
    CV.P.y_pred{held_out} = predict(mdl, x_P_test);

    % ---------------------------------------------------------
    % Integration  — select best t on training set
    % ---------------------------------------------------------
    best_bic_I = Inf;
    best_ti_I  = 1;
    % select the best t for this data set: (takes a few minutes depending
    % on the number of t windows ~15seconds/iteration)
    for ti = 1:t_num
        x_I_col   = X_I(:, ti);
        x_I_train = repmat(x_I_col, n_train, 1);
        valid     = ~isnan(x_I_train);
        if sum(valid) < 10, continue; end
        mdl = fitglm(x_I_train(valid), y_train(valid), ...
                     'Distribution', 'binomial', 'Link', 'logit');
        if mdl.ModelCriterion.BIC < best_bic_I
            best_bic_I = mdl.ModelCriterion.BIC;
            best_ti_I  = ti;
        end
    end

    % save the best fit version
    CV.I.t_chosen(held_out) = t_list(best_ti_I);

    x_I_train = repmat(X_I(:, best_ti_I), n_train, 1);
    x_I_test  = X_I(:, best_ti_I);
    valid_tr  = ~isnan(x_I_train);
    valid_te  = ~isnan(x_I_test);

    % run the model for this hold out data set with the
    % best t for this group
    mdl = fitglm(x_I_train(valid_tr), y_train(valid_tr), ...
                 'Distribution', 'binomial', 'Link', 'logit');

    y_pred_I = nan(length(y_test), 1);
    y_pred_I(valid_te) = predict(mdl, x_I_test(valid_te));
    CV.I.y_pred{held_out} = y_pred_I;

    % ---------------------------------------------------------
    % Derivative — select best t on training set
    % ---------------------------------------------------------
    best_bic_D = Inf;
    best_ti_D  = 1;
    for ti = 1:t_num
        x_D_col   = X_D(:, ti);
        x_D_train = repmat(x_D_col, n_train, 1);
        valid     = ~isnan(x_D_train);
        if sum(valid) < 10, continue; end
        mdl = fitglm(x_D_train(valid), y_train(valid), ...
                     'Distribution', 'binomial', 'Link', 'logit');
        if mdl.ModelCriterion.BIC < best_bic_D
            best_bic_D = mdl.ModelCriterion.BIC;
            best_ti_D  = ti;
        end
    end
    CV.D.t_chosen(held_out) = t_list(best_ti_D);

    x_D_train = repmat(X_D(:, best_ti_D), n_train, 1);
    x_D_test  = X_D(:, best_ti_D);
    valid_tr  = ~isnan(x_D_train);
    valid_te  = ~isnan(x_D_test);

    mdl = fitglm(x_D_train(valid_tr), y_train(valid_tr), ...
                 'Distribution', 'binomial', 'Link', 'logit');

    y_pred_D = nan(length(y_test), 1);
    y_pred_D(valid_te) = predict(mdl, x_D_test(valid_te));
    CV.D.y_pred{held_out} = y_pred_D;
 
end

%% ANALYSIS : PARALLEL STRUCTURES
% =========================================================
% STEP 2 — LOO cross-validation (parfor)
% =========================================================

% parfor cannot slice into structs directly — extract to plain variables
P_y_pred   = cell(nTrials, 1);
I_y_pred   = cell(nTrials, 1);
D_y_pred   = cell(nTrials, 1);
I_t_chosen = nan(nTrials, 1);
D_t_chosen = nan(nTrials, 1);

n_train = nTrials - 1;

% broadcast these to workers as local copies — parfor requires them to be
% explicitly named variables rather than struct fields for clean slicing
x_P_const = x_P;
X_I_const = X_I;
X_D_const = X_D;
t_list_const = t_list;
t_num_const  = t_num;

tic
parfor held_out = 1:nTrials

    train_trials = setdiff(1:nTrials, held_out);
    y_train      = y_data(:, train_trials);
    y_train      = y_train(:);
    y_test       = y_data(:, held_out);

    x_P_train = repmat(x_P_const, n_train, 1);
    x_P_test  = x_P_const;

    % ---------------------------------------------------------
    % Proportional
    % ---------------------------------------------------------
    mdl = fitglm(x_P_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');
    P_y_pred{held_out} = predict(mdl, x_P_test);

    % ---------------------------------------------------------
    % Integral — select best t on training set
    % ---------------------------------------------------------
    best_bic_I = Inf;
    best_ti_I  = 1;

    for ti = 1:t_num_const
        x_I_col   = X_I_const(:, ti);
        x_I_train = repmat(x_I_col, n_train, 1);
        valid     = ~isnan(x_I_train);
        if sum(valid) < 10, continue; end
        mdl = fitglm(x_I_train(valid), y_train(valid), ...
                     'Distribution', 'binomial', 'Link', 'logit');
        if mdl.ModelCriterion.BIC < best_bic_I
            best_bic_I = mdl.ModelCriterion.BIC;
            best_ti_I  = ti;
        end
    end

    I_t_chosen(held_out) = t_list_const(best_ti_I);

    x_I_train = repmat(X_I_const(:, best_ti_I), n_train, 1);
    x_I_test  = X_I_const(:, best_ti_I);
    valid_tr  = ~isnan(x_I_train);
    valid_te  = ~isnan(x_I_test);

    mdl = fitglm(x_I_train(valid_tr), y_train(valid_tr), ...
                 'Distribution', 'binomial', 'Link', 'logit');

    y_pred_I           = nan(length(y_test), 1);
    y_pred_I(valid_te) = predict(mdl, x_I_test(valid_te));
    I_y_pred{held_out} = y_pred_I;

    % ---------------------------------------------------------
    % Derivative — select best t on training set
    % ---------------------------------------------------------
    best_bic_D = Inf;
    best_ti_D  = 1;

    for ti = 1:t_num_const
        x_D_col   = X_D_const(:, ti);
        x_D_train = repmat(x_D_col, n_train, 1);
        valid     = ~isnan(x_D_train);
        if sum(valid) < 10, continue; end
        mdl = fitglm(x_D_train(valid), y_train(valid), ...
                     'Distribution', 'binomial', 'Link', 'logit');
        if mdl.ModelCriterion.BIC < best_bic_D
            best_bic_D = mdl.ModelCriterion.BIC;
            best_ti_D  = ti;
        end
    end

    D_t_chosen(held_out) = t_list_const(best_ti_D);

    x_D_train = repmat(X_D_const(:, best_ti_D), n_train, 1);
    x_D_test  = X_D_const(:, best_ti_D);
    valid_tr  = ~isnan(x_D_train);
    valid_te  = ~isnan(x_D_test);

    mdl = fitglm(x_D_train(valid_tr), y_train(valid_tr), ...
                 'Distribution', 'binomial', 'Link', 'logit');

    y_pred_D           = nan(length(y_test), 1);
    y_pred_D(valid_te) = predict(mdl, x_D_test(valid_te));
    D_y_pred{held_out} = y_pred_D;

end
toc

% collect results back into CV struct after parfor completes
CV.P.y_pred   = P_y_pred;
CV.I.y_pred   = I_y_pred;
CV.D.y_pred   = D_y_pred;
CV.I.t_chosen = I_t_chosen;
CV.D.t_chosen = D_t_chosen;

%%
% =========================================================
% STEP 3 — full-data fit using modal t
% =========================================================
t_best_I  = mode(CV.I.t_chosen);
t_best_D  = mode(CV.D.t_chosen);
ti_best_I = find(t_list == t_best_I);
ti_best_D = find(t_list == t_best_D);

x_I_best  = X_I(:, ti_best_I);
x_D_best  = X_D(:, ti_best_D);

y_all     = vertcat(y_data{:});
x_P_all   = repmat(x_P,      nTrials, 1);
x_I_all   = repmat(x_I_best, nTrials, 1);
x_D_all   = repmat(x_D_best, nTrials, 1);

valid_I   = ~isnan(x_I_all);
valid_D   = ~isnan(x_D_all);

results.P.mdl    = fitglm(x_P_all,          y_all,          'Distribution', 'binomial', 'Link', 'logit');
results.I.mdl    = fitglm(x_I_all(valid_I), y_all(valid_I), 'Distribution', 'binomial', 'Link', 'logit');
results.D.mdl    = fitglm(x_D_all(valid_D), y_all(valid_D), 'Distribution', 'binomial', 'Link', 'logit');
results.P.BIC    = results.P.mdl.ModelCriterion.BIC;
results.I.BIC    = results.I.mdl.ModelCriterion.BIC;
results.D.BIC    = results.D.mdl.ModelCriterion.BIC;
results.I.t_best = t_best_I;
results.D.t_best = t_best_D;

% =========================================================
% STEP 4 — CV performance metrics
% =========================================================
for m = 1:numel(fields)
    f = fields{m};

    y_pred_all = vertcat(CV.(f).y_pred{:});
    y_true_all = vertcat(y_data{:});
    valid      = ~isnan(y_pred_all);

    y_pred_all = y_pred_all(valid);
    y_true_all = y_true_all(valid);

    [~, ~, ~, CV.(f).AUC] = perfcurve(y_true_all, y_pred_all, 1);
    CV.(f).logloss = -mean(y_true_all .* log(y_pred_all + eps) + ...
                          (1 - y_true_all) .* log(1 - y_pred_all + eps));

    fprintf('%s — CV AUC: %.3f | Log-loss: %.3f', f, CV.(f).AUC, CV.(f).logloss);
    if ismember(f, {'I', 'D'})
        fprintf(' | t_best: %d (modal across folds)', results.(f).t_best);
    end
    fprintf('\n');
end


%%  Extras 
% 
% parfor ti = 1:t_num
%     X_I(:, ti) = computeIntegral(x_data{1},   t_list(ti), dt);
%     X_D(:, ti) = computeDerivative(x_data{1}, t_list(ti), dt);
% end