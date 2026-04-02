

%% PID Analysis Setup
clearvars('-except',initial_var{:})

dt = 1/num.fps; % sampling frequency
num_t = 25; % number of samples across time span to check for fit
max_t = 60; % max duration in the past to sample for window time (min)
t_list = linspace(1, max_t*num.fps*60, num_t);

best = struct();
results = struct();


%% Proportional Fit 
% need to update this to run for all the different trials (and for both
% sexes) and then to hold out some of the data randomly 

for trial = 1:num.trials
    x = data.temperature(:,trial);
    y = data.OutterRing(:, M, trial);

    mdl_P = fitglm(x, y, 'Distribution', 'binomial', 'Link', 'logit');

    results.P.mdl = mdl_P;
    results.P.BIC = mdl_P.ModelCriterion.BIC;


end




%% 
% --- Setup ---
dt = 1;
t_candidates = 5:5:100;
results = struct();

% =========================================================
% 1. PROPORTIONAL
% =========================================================
mdl_P = fitglm(x, y, 'Distribution', 'binomial', 'Link', 'logit');
results.P.mdl = mdl_P;
results.P.BIC = mdl_P.ModelCriterion.BIC;

% =========================================================
% 2. INTEGRAL
% =========================================================
best_bic_I = Inf;
for t = t_candidates
    x_I   = computeIntegral(x, t, dt);
    valid = ~isnan(x_I);
    mdl   = fitglm(x_I(valid), y(valid), 'Distribution', 'binomial', 'Link', 'logit');
    bic   = mdl.ModelCriterion.BIC;
    if bic < best_bic_I
        best_bic_I       = bic;
        results.I.mdl    = mdl;
        results.I.BIC    = bic;
        results.I.t_best = t;
        results.I.x_I    = x_I;
    end
end

% =========================================================
% 3. DERIVATIVE
% =========================================================
best_bic_D = Inf;
for t = t_candidates
    x_D   = computeDerivative(x, t, dt);
    valid = ~isnan(x_D);
    mdl   = fitglm(x_D(valid), y(valid), 'Distribution', 'binomial', 'Link', 'logit');
    bic   = mdl.ModelCriterion.BIC;
    if bic < best_bic_D
        best_bic_D       = bic;
        results.D.mdl    = mdl;
        results.D.BIC    = bic;
        results.D.t_best = t;
        results.D.x_D    = x_D;
    end
end

%% Functions

function x_I = computeIntegral(x, t, dt)
    % Trailing trapezoidal integral over window of t samples
    x_I = nan(size(x));
    for i = t:length(x)
        x_I(i) = trapz(x(i-t+1:i)) * dt;
    end
end

function x_D = computeDerivative(x, t, dt)
    % Slope from a linear fit over trailing window of t samples
    x_D = nan(size(x));
    tt = (0:t-1)' * dt;
    for i = t:length(x)
        p = polyfit(tt, x(i-t+1:i), 1);
        x_D(i) = p(1);           % slope = derivative estimate
    end
end

function auc = getAUC(mdl, y_true)
    y_pred = mdl.Fitted.Probability;
    [~, ~, ~, auc] = perfcurve(y_true, y_pred, 1);
end

% A Few Design Decisions to Think About
% Choosing t: AUC is a reasonable default criterion, but if you care about calibration (actual probabilities) rather than just ranking, use log-likelihood or BIC instead. You can grab it from mdl.ModelCriterion.BIC.
% Edge handling: The first t samples will be NaN for I and D since the window isn't full. The code above excludes them with valid. Make sure y is aligned with x so the exclusion is consistent.
% Derivative method: polyfit over the window gives a robust slope estimate that's less sensitive to noise than a simple two-point difference — good since you said x isn't guaranteed to be clean. If speed becomes an issue with large sweeps, you can replace it with (x(i) - x(i-t+1)) / ((t-1)*dt) as a faster approximation.
% Comparing across models: Since P has no t, and I/D have optimized t, their AUCs aren't directly comparable without correction for the search. Consider holding out a test set or using cross-validation for a fairer comparison.