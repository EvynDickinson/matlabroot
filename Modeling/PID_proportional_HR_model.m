
%% Proportional ~exponential~ Model

% Load in data from the high-resolution experiments and fit the temperature
% data to an exponential function and then get a measure of the fit for the
% data

% This is best used to model 'jump' escape data from the high resolution
% experiments

% Load in data via HighRes_Step6_1


%% Testing area
% 
% 
% tic
% X_I = nan(n_samples, t_num);
% X_D = nan(n_samples, t_num);
% 
% for ti = 1:t_num
%     X_I(:, ti) = computeIntegralLocal(x_data,   t_list(ti), dt);
%     X_D(:, ti) = computeDerivativeLocal(x_data, t_list(ti), dt);
% end
% x_P = x_data;   % proportional: raw signal, no windowing needed
% fprintf('\nTemperature transformations done\n')
% toc
% 

% tP = 25; % temperature of preference
% 
% % what do different transformations of temperature look like
% x_P = data.temp; % untransformed data
% x_P_centered = data.temp -tP; 
% x_P_nonlin_mag = (data.temp-tP).^2;
% 
% r = 1; c = 3;
% fig = figure; 
% subplot(r,c,1)
% plot(data.time, x_P)
% title('temperature')
% subplot(r,c,2)
% plot(data.time, x_P_centered)
% subplot(r,c,3)
% plot(data.time, x_P_nonlin_mag)


% what do the z-score versions of these inputs look like?


%% ANALYSIS & FIGURE: exponential fit of data to temperature
clearvars('-except',initial_var{:})
saveDir = createFolder([figDir, 'PID Model/']);


% select the type of data to model: 
data_type = 'jump';

x_data = smooth(data.temp, floor(num.fps), 'moving'); % half second smoothing for the input stim
y_data = [squeeze(data.(data_type)(:,M,:)), squeeze(data.(data_type)(:,F,:))]; % extract the output data

% reshape the data into a single vector
x = repmat(x_data, [num.trials*2, 1]);
x = x(:);
y = y_data(:);

% Find the locations for heating and cooling data (for later color
% plotting etc but need to match these locations with nan removals)
h_idx = repmat(data.tempbin.h_idx, [num.trials*2, 1]);
h_idx = h_idx(:);
c_idx = repmat(data.tempbin.c_idx, [num.trials*2, 1]);
c_idx = c_idx(:);

% remove nan data from all locations synchronously: 
loc = isnan(y);
x(loc) = [];
y(loc) = [];
c_idx(loc) = [];
h_idx(loc) = [];

% FIT DATA MODEL (simple exponential of the behavior)
fprintf('\nMaking model fit...')
[f, gof] = fit(x, y, 'exp1');
fprintf('done \n')

% PLOT THE DATA VS THE MODEL FIT
% Define bins **this might need to flex depending on temp protocol etc
nBins = 20;
edges = linspace(min(x), max(x), nBins + 1);
binCenters = (edges(1:end-1) + edges(2:end)) / 2; % what is the avg temp

% Assign each observation to a bin
binIdx = discretize(x, edges);  % returns 1..nBins, NaN outside edges

% Heating and cooling masks
y_h = y .* h_idx;
y_c = y .* c_idx;

% Compute stats per bin
binMean_h = accumarray(binIdx(h_idx), y(h_idx), [nBins 1], @mean,  NaN)';
binSEM_h  = accumarray(binIdx(h_idx), y(h_idx), [nBins 1], @(v) std(v)/sqrt(num.trials*2), NaN)';
binMean_c = accumarray(binIdx(c_idx), y(c_idx), [nBins 1], @mean,  NaN)';
binSEM_c  = accumarray(binIdx(c_idx), y(c_idx), [nBins 1], @(v) std(v)/sqrt(num.trials*2), NaN)';

% Evaluate fit over smooth x range
xFit = linspace(min(x), max(x), 300)'; % 300 data points to test for a smooth fit
yFit = feval(f, xFit);
yFit_bins = feval(f, binCenters)'; % evaluation at the center of the bins

% ASSESSING MODEL FITS: 
% weighted binned R2:
binN_h = histcounts(x(h_idx), edges);
ss_res_wh = sum(binN_h .* (binMean_h - yFit_bins).^2, 'omitnan');
ss_tot_wh = sum(binN_h .* (binMean_h - mean(binMean_h, 'omitnan')).^2, 'omitnan');
r2_weighted_h = 1 - ss_res_wh / ss_tot_wh;

binN_c = histcounts(x(c_idx), edges);
ss_res_wc = sum(binN_c .* (binMean_c - yFit_bins).^2, 'omitnan');
ss_tot_wc = sum(binN_c .* (binMean_c - mean(binMean_c, 'omitnan')).^2, 'omitnan');
r2_weighted_c = 1 - ss_res_wc / ss_tot_wc;

% For heating
yhat_h = feval(f, x(h_idx)); % predicted values
y_h    = y(h_idx);
ss_res_h = sum((y_h - yhat_h).^2, 'omitnan');
ss_tot_h = sum((y_h - mean(y_h, 'omitnan')).^2, 'omitnan');
raw_r2_h = 1 - ss_res_h / ss_tot_h;

% For cooling
yhat_c = feval(f, x(c_idx));
y_c    = y(c_idx);
ss_res_c = sum((y_c - yhat_c).^2, 'omitnan');
ss_tot_c = sum((y_c - mean(y_c, 'omitnan')).^2, 'omitnan');
raw_r2_c = 1 - ss_res_c / ss_tot_c;

fprintf('R² heating: %.3f,  R² cooling: %.3f\n', raw_r2_h, raw_r2_c);


% Get model-predicted probabilities (clamp to avoid log(0))
p_hat = feval(f, x);
p_hat = max(min(p_hat, 1-1e-10), 1e-10);  % clamp to (0,1)

% Bernoulli log-likelihood
ll = sum(y .* log(p_hat) + (1 - y) .* log(1 - p_hat), 'omitnan');
ll_per_obs = ll / numel(y);
fprintf('Log-likelihood per observation: %.6f\n', ll_per_obs);

% Null model log-likelihood (intercept only — just the base rate)
p_null = mean(y, 'omitnan');
ll_null = sum(y .* log(p_null) + (1 - y) .* log(1 - p_null), 'omitnan');
ll_null_per_obs = ll_null / numel(y);

% McFadden's pseudo-R² (0–1 scale, analogous to R²)
pseudo_r2 = 1 - (ll / ll_null);

fprintf('Log-likelihood per observation (model): %.6f\n', ll_per_obs);
fprintf('Log-likelihood per observation (null):  %.6f\n', ll_null_per_obs);
fprintf('Difference:  %.6f\n', ll_per_obs - ll_null_per_obs);
fprintf('McFadden pseudo-R²:    %.4f\n', pseudo_r2);


% Evaluate fit at bin centers
yFit_bins = feval(f, binCenters(:))';  % force row vector

% Heating
ss_res_h = nansum((binMean_h - yFit_bins).^2);
ss_tot_h = nansum((binMean_h - nanmean(binMean_h)).^2);
r2_h = 1 - ss_res_h / ss_tot_h;

% Cooling
ss_res_c = nansum((binMean_c - yFit_bins).^2);
ss_tot_c = nansum((binMean_c - nanmean(binMean_c)).^2);
r2_c = 1 - ss_res_c / ss_tot_c;

fprintf('Heating  — R²: %.4f,  Weighted R²: %.4f\n', r2_h, r2_weighted_h);
fprintf('Cooling  — R²: %.4f,  Weighted R²: %.4f\n', r2_c, r2_weighted_c);

%% FIGURE: Plot the temperature vs binned behavior RAW VS FIT
foreColor = formattingColors(blkbgd);
% plot parameters
LW = 2; % line width
CZ = 5; % error bar cap size
MZ = 7; % marker size 
hColor = Color('red');
cColor = Color('dodgerblue');

% legend captions: 
fit_str =  sprintf('Fit: %.3g · e^{%.3g·x}', f.a, f.b);
c_str = sprintf('Cooling: R^2=%.2g', r2_weighted_c);
h_str = sprintf('Heating: R^2=%.2g', r2_weighted_h);

% Plot raw binned data vs fit model data
fig = getfig('', 1,[554 635]); 
    hold on
    % cooling data
    errorbar(binCenters, binMean_c, binSEM_c, 'o', 'MarkerSize', MZ,...
        'MarkerFaceColor', cColor, 'MarkerEdgeColor', 'none', ...
        'Color', cColor,  'LineWidth', LW, 'CapSize', CZ);
    % heating data
    errorbar(binCenters, binMean_h, binSEM_h, 'o','MarkerSize', MZ,...
        'MarkerFaceColor', hColor,  'MarkerEdgeColor', 'none', ...
        'Color', hColor, 'LineWidth', LW, 'CapSize', CZ);
    % model fit
    plot(xFit, yFit,  'Color', foreColor, 'LineWidth', LW);
    % formatting
    xlabel('temperature (\circC)');
    ylabel([data_type ' (fraction of flies)']);
    title(sprintf('Binned %s vs Temperature', data_type));
    legend(c_str, h_str, fit_str, ...
        'Location', 'northwest', 'box', 'off', 'TextColor', foreColor);
     formatFig(fig, blkbgd);
    ylim([-0.001, max(ylim)])

    % set(gca, 'YScale', 'log')
  
save_figure(fig, [saveDir, 'exponential proportional fit ', data_type]);
   

%% FIGURE: plot the measures of fit metrics 
% Organize metrics into a table-style figure
metrics = {
    'Raw R² (heating)',        raw_r2_h;
    'Raw R² (cooling)',        raw_r2_c;
    'Binned R² (heating)',     r2_h;
    'Binned R² (cooling)',     r2_c;
    'Weighted R² (heating)',   r2_weighted_h;
    'Weighted R² (cooling)',   r2_weighted_c;
    'McFadden pseudo-R²',      pseudo_r2;
    'Log-likelihood',          ll_per_obs;
    'Null log-likelihood',     ll_null_per_obs;
};

% Split into two groups: R²-like (0-1) and log-likelihoods
r2_labels = metrics(1:7, 1);
r2_vals   = cell2mat(metrics(1:7, 2));
ll_labels = metrics(8:9, 1);
ll_vals   = cell2mat(metrics(8:9, 2));

% Colors: heating = red, cooling = blue, combined/other = grey
hColor = Color('red');
cColor = Color('dodgerblue');
colors = [
    hColor;   % raw R² heating
    cColor;   % raw R² cooling
    hColor;   % binned R² heating
    cColor;   % binned R² cooling
    hColor;   % weighted R² heating
    Color('grey');   % weighted R² cooling
    Color('lightgrey')];    % McFadden


fig = getfig('',1);
% --- Top panel: R²-like metrics ---
ax1 = subplot(2,1,1); hold on;
for i = 1:numel(r2_vals)
    bar(i, r2_vals(i), 'FaceColor', colors(i,:), 'EdgeColor', 'none');
end
yline(0, 'k-', 'LineWidth', 0.8);
set(ax1, 'XTick', 1:numel(r2_labels), 'XTickLabel', r2_labels, ...
    'XTickLabelRotation', 30, 'TickLabelInterpreter', 'none');
ylabel('R² value');
title(sprintf('Model fit metrics — %s', data_type));
box off;
for i = 1:numel(r2_vals)
    if r2_vals(i) >= 0
        ypos = r2_vals(i) + 0.02;
        valign = 'bottom';
    else
        ypos = 0.02;  % just above zero line regardless of how negative the bar is
        valign = 'bottom';
    end
    text(i, ypos, sprintf('%.4f', r2_vals(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   valign, ...
        'FontSize', 8, ...
        'Color',    colors(i,:));
end
ylim([min(r2_vals)*1.2, max(r2_vals)*1.2]);

% --- Bottom panel: log-likelihoods ---
ax2 = subplot(2,1,2); hold on;
bar(1, ll_vals(1),  'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
bar(2, ll_vals(2),  'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
set(ax2, 'XTick', 1:2, 'XTickLabel', ll_labels, ...
    'XTickLabelRotation', 30, 'TickLabelInterpreter', 'none');
ylabel('Log-likelihood per observation');
box off;

tight_layout = @() set(fig, 'Units', 'normalized');
tight_layout();
formatFig(fig, blkbgd, [2,1]);
save_figure(fig, [saveDir, 'exponential proportional fit model fit ', data_type]);

%% FIGURE: Plot the temp vs transformed temp metric

T0 = 25;  % reference temperature (preference point)

x_range  = linspace(min(x), max(x), 300)';
y_prop   = exp(f.b .* (x_range - T0));  % fold-change in behavioral drive relative to T0

fig = getfig('', 1, [400 500]);
hold on;
yline(1, '--', 'Color', foreColor, 'LineWidth', 1, 'Alpha', 0.5);  % reference line at T0
plot(x_range, y_prop, 'Color', foreColor, 'LineWidth', LW); %Color('vaporwavepurple')
scatter(T0, 1, 60, foreColor);  % mark T0

xlabel('temperature (\circC)');
ylabel({'relative behavioral drive'; '(fold change from 25°C)'});
title({'"proportional" thermal gain'; 'function for jumping'});
formatFig(fig, blkbgd);
ylim([-2, max(ylim)])

save_figure(fig, [saveDir, 'proportional thermal gain function ']);

%% FIGURE: temp and transformed temp overlaid
% how does the transformed temperature look over time for this given temp protocol?

% Compute drive signal over time
drive = exp(f.b .* (x_data - T0));  % fold-change drive at each timepoint
LW = 2;

fig = getfig('', 1);
    hold on;
    % plot cooling and heating segments separately
    plot(data.time, drive, 'Color', foreColor , 'linewidth', LW);
    yline(1, '--', 'Color', Color('vaporwavepink'), 'LineWidth', 1);  % T0 reference
    xlabel('time (min)');
    ylabel({'relative behavioral drive'; '(fold change from 25°C)'});title('Thermal drive signal over time');
    formatFig(fig, blkbgd);
    ylim([-2, max(ylim)])
    
    yyaxis right
    plot(data.time, x_data, 'Color', Color('grey'), 'LineWidth', LW, 'linestyle', '-');
    ylabel('temperature (\circC)');
    set(gca, 'ycolor', Color('grey'))

    yyaxis left 
    set(gca, 'ycolor', foreColor)

save_figure(fig, [saveDir, 'proportional thermal drive signal over time ', data_type]);

%% FIGURE:

% Compute drive signal over time
drive = exp(f.b .* (x_data - T0));  % fold-change drive at each timepoint
LW = 2;

fig = getfig('', 1,[1000 500]);

% --- Subplot 1: raw temperature over time ---
subplot(1,3,1); hold on;
plot(data.time, x_data, 'Color', foreColor, 'LineWidth', LW);
yline(T0, '--', 'Color', Color('vaporwavepink'), 'LineWidth', 1);
xlabel('time (min)');
ylabel('temperature (\circC)');
title({'Raw temperature'; ' '});
xlim([-10, 400])


% --- Subplot 2: drive signal over time ---
subplot(1,3,2); hold on;
plot(data.time, drive, 'Color', foreColor, 'LineWidth', LW);
yline(1, '--', 'Color', Color('vaporwavepink'), 'LineWidth', 1);
xlabel('time (min)');
ylabel({'relative behavioral drive'; '(fold change from 25°C)'});
title({'Proportional thermal drive'; ' '});
ylim([-2, max(ylim)]);
xlim([-10, 400])

% --- Subplot 3: gain function (temp vs drive) ---
subplot(1,3,3); hold on;
x_range = linspace(min(x), max(x), 300)';
y_prop  = exp(f.b .* (x_range - T0));
yline(1, '--', 'Color', Color('vaporwavepink'), 'LineWidth', 1);
plot(x_range, y_prop, 'Color', foreColor, 'LineWidth', LW);
scatter(T0, 1, 60, Color('vaporwavepink'), 'filled');
xlabel('temperature (\circC)');
ylabel({'relative behavioral drive'; '(fold change from 25°C)'});
title({'Thermal gain function'; ' '});
ylim([-2, max(ylim)]);
xlim([10, 40])

formatFig(fig, blkbgd, [1,3]);
save_figure(fig, [saveDir, 'proportional thermal drive signal ', data_type]);



%% FIGURE: data + model summary
[foreColor, backColor] = formattingColors(blkbgd);

baselineColor = Color('WongOrange');
LW = 2;
smooth_win = floor(num.fps) * 10; % smooth window for display
time_xlim = [-10, 400];

% set up figure aligments
r = 5;
c = 4;
sb(3).idx = 1:c:r*c;               % left column  (1,5,9,13,17)
sb(1).idx = [3,4];                  % top — skips col 2 (blank)
sb(2).idx = [7,8,11,12,15,16,19,20]; % remaining rows — skips col 2


fig = getfig('', false, [1238 876]);

% --- Subplot 1: raw temperature over time ---
ax1 = subplot(r,c,sb(1).idx); hold on
plot(data.time, x_data, 'Color', foreColor, 'LineWidth', LW)
yline(T0, '--', 'Color', baselineColor, 'LineWidth', 1)
xlabel('time (min)')
ylabel('temperature (\circC)')
% title({'Temperature protocol'; ' '})
xlim(time_xlim)
set(ax1, 'ytick', 15:10:35)
ylim([15, 35])

% --- Subplot 2: raw data timecourse + model predicted drive ---
ax2 = subplot(r,c,sb(2).idx); hold on

% raw data — average across trials and sexes for display
y_avg = mean(y_data, 2, 'omitnan');
y_avg_smooth = smooth(y_avg, smooth_win, 'moving');
% model predicted drive over time
drive = exp(f.b .* (x_data - T0));
% normalize both to [0 1] for overlay on same axis
norm = @(v) (v - min(v)) ./ (max(v) - min(v));
plot(data.time, norm(y_avg_smooth), 'Color', foreColor, 'LineWidth', LW)
plot(data.time, norm(drive), 'Color', Color('WongGreen'), ...
    'LineWidth', LW, 'LineStyle', '-')
l = legend({'real data', 'model'}, 'Color', backColor, 'TextColor', foreColor, 'Box', 'off');
% % update legend: 
% l.Color = backColor;
% l.TextColor = foreColor;
% l.Box = 'off';
% yline(0, '--', 'Color', foreColor, 'LineWidth', 0.5, 'Alpha', 0.3)
xlabel('time (min)')
ylabel('jumps (min-max normalization)')
% title({'Data vs predicted drive'; ' '})
xlim(time_xlim)
ylim([0, 1])

% --- Subplot 3: thermal gain function ---
ax3 = subplot(r,c,sb(3).idx); hold on
x_range = linspace(min(x), max(x), 300)';
y_prop = exp(f.b .* (x_range - T0));
yline(1, '--', 'Color', baselineColor, 'LineWidth', 1)
plot(x_range, y_prop, 'Color', foreColor, 'LineWidth', LW)
scatter(T0, 1, 60, baselineColor, 'filled')
xlabel('temperature (\circC)')
ylabel('relative behavioral drive (fold change from 25\circC)')
% title({'Thermal gain function'; ' '})
ylim([-2, max(ylim)])
xlim([10, 40])

formatFig(fig, blkbgd, [r,c], sb);
set(ax1, 'XColor', 'none')
set(ax2, 'XColor', 'none')
linkaxes([ax1,ax2],'X')


% Add scale bar to the bottom subplot (ax4)
axes(ax2)
xl = xlim;
yl = ylim;

% --- Scale bar parameters ---
scaleBar_duration = 100;        % length in x-axis units (e.g., min)
scaleBar_x_start = xl(1) + 0.03 * diff(xl);   % left-aligned with small margin
scaleBar_x_end   = scaleBar_x_start + scaleBar_duration;
scaleBar_y = yl(1) - 0.04 * diff(yl);       % just below the data

% Draw the scale bar line
annotation_ax = ax2;
line(annotation_ax, ...
    [scaleBar_x_start, scaleBar_x_end], ...
    [scaleBar_y, scaleBar_y], ...
    'Color', foreColor, 'LineWidth', 2, 'Clipping', 'off', 'HandleVisibility','off')

% Add label centered above the bar
scale_label_str = sprintf('%d min', scaleBar_duration);
text(ax2, scaleBar_x_start, scaleBar_y - 0.04*diff(yl), ...
    scale_label_str, ...
    'Color', foreColor, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 18, ...
    'Clipping', 'off', 'Tag', 'ScaleBar')
addlistener(ax2, 'XLim', 'PostSet', @(~,~) updateScaleBar(ax4, foreColor, scaleBar_duration, scale_label_str));

% save figure
save_figure(fig, [saveDir, 'data and model fit ', data_type]);

