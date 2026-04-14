

%% PID Analysis Setup
clearvars('-except',initial_var{:})

saveDir = createFolder([figDir, 'PID Model/']);

dt = 1/num.fps; % sampling frequency
t_num = 10; % number of samples across time span to check for fit
max_t = 90; % max duration in the past to sample for window time (min)
t_list = round(linspace(num.fps, max_t*num.fps*60, t_num)); % from 1 second to 60 minutes
t_list_time = t_list/num.fps/60;


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

%% ANALYSIS: Pull the data together across trials:  

% =========================================================
% STEP 1 — precompute transforms (x is same across trials)
% =========================================================
% note: takes a few seconds to run 
tic
X_I = nan(n_samples, t_num);
X_D = nan(n_samples, t_num);
for ti = 1:t_num
    X_I(:, ti) = computeIntegral(abs(x_data),   t_list(ti), dt);
    X_D(:, ti) = computeDerivative(x_data, t_list(ti), dt);
end
x_P = x_data;
toc

%% Comparison of integral temperature transformations

% X_I_shifted = nan(n_samples, t_num);
X_I_abs = nan(n_samples, t_num);
for ti = 1:t_num
     X_I_abs(:, ti) = computeIntegral(abs(x_data),   t_list(ti), dt);
end



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
    plot(t_list_time, Color='m', LineStyle=':')
    for ii = 1:t_num
        scatter(ii, t_list_time(ii), 35, cList(ii,:), 'filled')
    end
    formatFig(fig, blkbgd);
    xlabel('t index')
    ylabel('t window (min)')
    save_figure(fig, [saveDir, 'time lag t']);



%% ANALYSIS : PARALLEL STRUCTURES
% % take ~90 mins to run 40 trials with 10 time lag iterations
% % =========================================================
% % STEP 2 — LOO cross-validation (parfor)
% % =========================================================
% 
% % turn on and create a parallel processing pool
% parpool();
% % 16 workers on the mac 
% % 16 workers on ACADIA
% 
% 
% % parfor cannot slice into structs directly — extract to plain variables
% P_y_pred   = cell(nTrials, 1);
% I_y_pred   = cell(nTrials, 1);
% D_y_pred   = cell(nTrials, 1);
% I_t_chosen = nan(nTrials, 1);
% D_t_chosen = nan(nTrials, 1);
% 
% n_train = nTrials - 1;
% 
% % broadcast these to workers as local copies — parfor requires them to be
% % explicitly named variables rather than struct fields for clean slicing
% x_P_const = x_P;
% X_I_const = X_I;
% X_D_const = X_D;
% t_list_const = t_list;
% t_num_const  = t_num;
% 
% tic
% parfor held_out = 1:nTrials
% 
%     train_trials = setdiff(1:nTrials, held_out);
%     y_train      = y_data(:, train_trials);
%     y_train      = y_train(:);
%     y_test       = y_data(:, held_out);
% 
%     x_P_train = repmat(x_P_const, n_train, 1);
%     x_P_test  = x_P_const;
% 
%     % ---------------------------------------------------------
%     % Proportional
%     % ---------------------------------------------------------
%     mdl = fitglm(x_P_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');
%     P_y_pred{held_out} = predict(mdl, x_P_test);
% 
%     % ---------------------------------------------------------
%     % Integral — select best t on training set
%     % ---------------------------------------------------------
%     best_bic_I = Inf;
%     best_ti_I  = 1;
% 
%     for ti = 1:t_num_const
%         x_I_col   = X_I_const(:, ti);
%         x_I_train = repmat(x_I_col, n_train, 1);
%         valid     = ~isnan(x_I_train);
%         if sum(valid) < 10, continue; end
%         mdl = fitglm(x_I_train(valid), y_train(valid), ...
%                      'Distribution', 'binomial', 'Link', 'logit');
%         if mdl.ModelCriterion.BIC < best_bic_I
%             best_bic_I = mdl.ModelCriterion.BIC;
%             best_ti_I  = ti;
%         end
%     end
% 
%     I_t_chosen(held_out) = t_list_const(best_ti_I);
% 
%     x_I_train = repmat(X_I_const(:, best_ti_I), n_train, 1);
%     x_I_test  = X_I_const(:, best_ti_I);
%     valid_tr  = ~isnan(x_I_train);
%     valid_te  = ~isnan(x_I_test);
% 
%     mdl = fitglm(x_I_train(valid_tr), y_train(valid_tr), ...
%                  'Distribution', 'binomial', 'Link', 'logit');
% 
%     y_pred_I           = nan(length(y_test), 1);
%     y_pred_I(valid_te) = predict(mdl, x_I_test(valid_te));
%     I_y_pred{held_out} = y_pred_I;
% 
%     % ---------------------------------------------------------
%     % Derivative — select best t on training set
%     % ---------------------------------------------------------
%     best_bic_D = Inf;
%     best_ti_D  = 1;
% 
%     for ti = 1:t_num_const
%         x_D_col   = X_D_const(:, ti);
%         x_D_train = repmat(x_D_col, n_train, 1);
%         valid     = ~isnan(x_D_train);
%         if sum(valid) < 10, continue; end
%         mdl = fitglm(x_D_train(valid), y_train(valid), ...
%                      'Distribution', 'binomial', 'Link', 'logit');
%         if mdl.ModelCriterion.BIC < best_bic_D
%             best_bic_D = mdl.ModelCriterion.BIC;
%             best_ti_D  = ti;
%         end
%     end
% 
%     D_t_chosen(held_out) = t_list_const(best_ti_D);
% 
%     x_D_train = repmat(X_D_const(:, best_ti_D), n_train, 1);
%     x_D_test  = X_D_const(:, best_ti_D);
%     valid_tr  = ~isnan(x_D_train);
%     valid_te  = ~isnan(x_D_test);
% 
%     mdl = fitglm(x_D_train(valid_tr), y_train(valid_tr), ...
%                  'Distribution', 'binomial', 'Link', 'logit');
% 
%     y_pred_D           = nan(length(y_test), 1);
%     y_pred_D(valid_te) = predict(mdl, x_D_test(valid_te));
%     D_y_pred{held_out} = y_pred_D;
% 
% end
% toc
% 
% % collect results back into CV struct after parfor completes
% CV.P.y_pred   = P_y_pred;
% CV.I.y_pred   = I_y_pred;
% CV.D.y_pred   = D_y_pred;
% CV.I.t_chosen = I_t_chosen;
% CV.D.t_chosen = D_t_chosen;

%%
% % =========================================================
% % STEP 3 — full-data fit using modal t
% % =========================================================
% t_best_I  = mode(CV.I.t_chosen);
% t_best_D  = mode(CV.D.t_chosen);
% ti_best_I = find(t_list == t_best_I);
% ti_best_D = find(t_list == t_best_D);
% 
% x_I_best  = X_I(:, ti_best_I);
% x_D_best  = X_D(:, ti_best_D);
% 
% y_all     = y_data(:);
% x_P_all   = repmat(x_P,      nTrials, 1);
% x_I_all   = repmat(x_I_best, nTrials, 1);
% x_D_all   = repmat(x_D_best, nTrials, 1);
% 
% valid_I   = ~isnan(x_I_all);
% valid_D   = ~isnan(x_D_all);
% 
% results.P.mdl    = fitglm(x_P_all,          y_all,          'Distribution', 'binomial', 'Link', 'logit');
% results.I.mdl    = fitglm(x_I_all(valid_I), y_all(valid_I), 'Distribution', 'binomial', 'Link', 'logit');
% results.D.mdl    = fitglm(x_D_all(valid_D), y_all(valid_D), 'Distribution', 'binomial', 'Link', 'logit');
% results.P.BIC    = results.P.mdl.ModelCriterion.BIC;
% results.I.BIC    = results.I.mdl.ModelCriterion.BIC;
% results.D.BIC    = results.D.mdl.ModelCriterion.BIC;
% results.I.t_best = t_best_I;
% results.D.t_best = t_best_D;
% 
% % Visual check of stability of t (after step 3) 
% figure;
% subplot(1,2,1); histogram(CV.I.t_chosen); xlabel('t chosen (I)'); title('Integral t stability');
% subplot(1,2,2); histogram(CV.D.t_chosen); xlabel('t chosen (D)'); title('Derivative t stability');
% 
% % =========================================================
% % STEP 4 — CV performance metrics
% % =========================================================
% for m = 1:numel(fields)
%     f = fields{m};
% 
%     y_pred_all = vertcat(CV.(f).y_pred{:});
%     y_true_all = y_data(:);
%     valid      = ~isnan(y_pred_all);
% 
%     y_pred_all = y_pred_all(valid);
%     y_true_all = y_true_all(valid);
% 
%     [~, ~, ~, CV.(f).AUC] = perfcurve(y_true_all, y_pred_all, 1);
%     CV.(f).logloss = -mean(y_true_all .* log(y_pred_all + eps) + ...
%                           (1 - y_true_all) .* log(1 - y_pred_all + eps));
% 
%     fprintf('%s — CV AUC: %.3f | Log-loss: %.3f', f, CV.(f).AUC, CV.(f).logloss);
%     if ismember(f, {'I', 'D'})
%         fprintf(' | t_best: %d (modal across folds)', results.(f).t_best);
%     end
%     fprintf('\n');
% end


%%  Extras 
% 
% parfor ti = 1:t_num
%     X_I(:, ti) = computeIntegral(x_data{1},   t_list(ti), dt);
%     X_D(:, ti) = computeDerivative(x_data{1}, t_list(ti), dt);
% end



%% 
% why is the selected time bin for t always the same ?
% check the model fits and process


%% Quick visuals: new tuning curves based on the different input parameters

% how do the different integration windows line up with the behavior
% viusally?

% need to find a way to bin the predictor variable somehow ...
% could just plot it as x vs y
tt = 10; 
y_avg = mean(y_data,2, 'omitnan');

typeList = {'proportional', 'integrative', 'derivative'};
roiList = {'h_idx', 'c_idx'};
cList = {'red', 'dodgerblue'};
r = 1; c = 2;
FA = 0.3;
SZ = 35; % plot size
x_all = [x_P, X_I(:,tt), X_D(:,tt)];

for type  = 1:3

    x = x_all(:,type) ;
    typeStr = typeList{type};
    
    fig = getfig('', 1,[ 954 543]); 
    for ii = 1:2 % heating and cooling

        subplot(r,c,ii); hold on

        roi = data.tempbin.(roiList{ii});
        Y = y_avg(roi);
        X = x(roi);
        [X, Y] = rmnan(X, Y, 2, true);

        % raw data
        scatter(X, Y, SZ, Color(cList{ii}), "filled", 'MarkerFaceAlpha', FA)
        % linear line of best fit
        coefficients = polyfit(X, Y, 1);
        xFit = linspace(min(X), max(X), 100); % Create smooth x-values
        yFit = polyval(coefficients, xFit);   % Calculate corresponding y-values
        plot(xFit, yFit, Color=foreColor, LineWidth=2);
        rho = corr(X, Y);
        xlabel(['Temp (' typeStr ')'])
        ylabel(data_type)
        title([typeStr ' rho: ' num2str(rho)])

    end
    formatFig(fig, blkbgd, [r,c]);
    matchSubplotAxes(fig);

end


%% ANALYSIS: H & C separated values 
% what transformation of temperature gives the best fit?

data_type = 'jump';
y_data = [squeeze(data.(data_type)(:,M,:)), squeeze(data.(data_type)(:,F,:))]; 
roiList = {'h_idx', 'c_idx'};

% initialize variables: 
results = struct; % where the model fits will go


% PROPORTIONAL MODEL :
x = x_P;
rho = nan(nTrials, 2);
for type = 1:2 % heating and cooling
    roi = data.tempbin.(roiList{type}); % data points to include
    for trial = 1:nTrials
        y = y_data(roi,trial);
        x = x_P(roi);
        [X,Y] = rmnan(x, y, 2, true);
        rho(trial,type) = corr(X, Y);
    end
    results(1).rho = rho; % save the correlation coeffecients
end
fig = getfig('',1,[308 504]); hold on
    scatter(1:nTrials,rho(:,1), 35, 'r', 'filled', 'MarkerFaceAlpha',0.8)
    scatter(1:nTrials,rho(:,2), 35, 'b', 'filled', 'MarkerFaceAlpha',0.8)
    xlim([-5, nTrials+5])
    h_line(0, 'grey');
    title(data_type)
    ylabel('correlation coefficient')
    formatFig(fig, blkbgd);
    set(gca, XColor = 'none')
save_figure(fig, [saveDir, data_type ' proportional coeff scatter']);



% INTEGRAL MODEL : 
rho = nan(nTrials, t_num,2);
for type = 1:2 % heating then cooling
    roi = logical(replaceNaN(data.tempbin.(roiList{type}),false));
    for tt = 1:t_num
        for trial = 1:nTrials
            x = X_I(roi,tt);
            y = y_data(roi,trial);
            [X,Y] = rmnan(x,y,2, true);
            rho(trial,tt,type) = corr(X, Y);
        end
    end
end
results(2).rho = rho;
fig = getfig('',1,[1040 650]); hold on
    xC = t_list_time-1;
    xH = t_list_time+1;
    scatter(repmat(xH, nTrials,1), rho(:,:,1), 35, 'r','filled', XJitter='density',XJitterWidth=0.1,MarkerFaceAlpha=0.5);
    scatter(repmat(xC, nTrials,1), rho(:,:,2), 35, 'b', 'filled', XJitter='density',XJitterWidth=0.1,MarkerFaceAlpha=0.5);
    scatter(xH, median(rho(:,:,1),1, 'omitnan'), 50, foreColor, 'filled')
    scatter(xC, median(rho(:,:,2),1, 'omitnan'), 50, foreColor, 'filled')
    h_line(0)
    xlabel('time integration window (min)')
    ylabel('correlation coefficent (pearsons)')
    title([data_type ' Integrative Model'])
formatFig(fig, blkbgd);
xlim([-5,max(t_list_time)+5])
save_figure(fig, [saveDir, data_type ' integral coeff scatter']);


% DERIVATIVE MODEL : 
rho = nan(nTrials, t_num,2);
for type = 1:2 % heating then cooling
    roi = logical(replaceNaN(data.tempbin.(roiList{type}),false));
    for tt = 1:t_num
        for trial = 1:nTrials
            x = X_D(roi,tt);
            y = y_data(roi,trial);
            [X,Y] = rmnan(x,y,2, true);
            rho(trial,tt,type) = corr(X, Y);
        end
    end
end
results(3).rho = rho;
fig = getfig('',1,[1040 650]); hold on
    xC = t_list_time-1;
    xH = t_list_time+1;
    scatter(repmat(xH, nTrials,1), rho(:,:,1), 35, 'r','filled', XJitter='density',XJitterWidth=0.1,MarkerFaceAlpha=0.5);
    scatter(repmat(xC, nTrials,1), rho(:,:,2), 35, 'b', 'filled', XJitter='density',XJitterWidth=0.1,MarkerFaceAlpha=0.5);
    scatter(xH, mode(rho(:,:,1),1, 'omitnan'), 50, foreColor, 'filled')
    scatter(xC, mode(rho(:,:,2),1, 'omitnan'), 50, foreColor, 'filled')
    h_line(0)
    title([data_type ' Derivative Model'])
    xlabel('time integration window (min)')
    ylabel('correlation coefficent (pearsons)')
formatFig(fig, blkbgd);
xlim([-5,max(t_list_time)+5])
save_figure(fig, [saveDir, data_type ' derivative coeff scatter']);


%% FIGURE: replot the temp tuning curve data simply by the transformed temperture data: 
nbins = 20;
D_best = 2;
I_best = 2;

x_sources = {x_P, X_D(:,D_best), X_I(:,D_best)};
x_labels   = {'Proportional', 'Derivative', 'Integrative'};

h_idx = data.tempbin.h_idx;
c_idx = data.tempbin.c_idx;
idx_sets = {h_idx, c_idx};

kolor = Color('magenta', 'dodgerblue', 2);
fig = figure;

for xx = 1:3
    x = x_sources{xx};
    [c, edges] = discretize(x, nbins);
    temp_ax = edges(1:end-1) + diff(edges);

    % bin and sort data: 
    y_plot = struct();
    [y_plot(1).all, y_plot(2).all] = deal(nan(nTrials, nbins));
    [y_plot(1).avg, y_plot(2).avg, y_plot(1).sem, y_plot(2).sem] = deal(nan(nbins,1));

    for tt = 1:2
        for ii = 1:nbins
            roi = (c == ii) & idx_sets{tt};
            y = y_data(roi, :)*100;
            y_plot(tt).all(:,ii) = mean(y, 1, 'omitnan');
        end
        y_plot(tt).avg = mean(y_plot(tt).all, 1, 'omitnan')';
        y_plot(tt).sem = std( y_plot(tt).all, 0, 1, 'omitnan')' / sqrt(nTrials);
    end
    % plot the data for this transformation type:
    subplot(1, 3, xx); hold on;
    title(x_labels{xx});
    for tt = 1:2
        plot(temp_ax, y_plot(tt).avg, 'color', kolor(tt,:));
        plot_error_fills(true, temp_ax, y_plot(tt).avg, y_plot(tt).sem, kolor(tt,:), '-png');
    end
    ylabel('flies in escape ring (%)')
end
formatFig(fig, blkbgd, [1,3]);
save_figure(fig, [saveDir, 'time lag t']);



%% Analysis: Plot out the best models of the data based on the best time lags 
for ii = 2:3
    rho = results(ii).rho; % find all the correlation coefficients
    avg = squeeze(median(rho,1)); % find the average of the heating and cooling correlations
    [~,idx] = (max(mean(avg,2))); % find the strongest correlation for the model base
    results(ii).rho_idx = idx; % save the location of the correlation

    % smooth bin the data using this time lag: 
    for type = 1:2 % heating and cooling (same time lag for both? -- doesn't have to be...)
        roi = logical(replaceNaN(data.tempbin.(roiList{type}),false));
        switch ii
            case 1
            case 2 % integral
                x_sel = X_I(:,idx);
            case 3 % derivative
                x_sel = X_D(:,idx);
        end

        % find the population correlation
        x_temp = x_sel(roi);
        x = repmat(x_temp,[1,nTrials]);
        [X,Y] = rmnan(x(:),y(:),2, true);
        p_corr = corr(X, Y);

        % pull up the predictor and response variable: 
        [x_idx, Edges] = discretize(x_sel,nbins);
        x = x_idx(roi,:);
        y = y_data(roi,:);
        
        % pull the avg response value for each of the trials in each bin: 
        for bin_i = 1:nbins
            x_roi = % WORKING HERE: PLOT THIS BY BIN TYPE WITH ERROR
        end







