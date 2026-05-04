
saveDir = createFolder([figDir, 'Paper Figures/']);
initial_var = add_var(initial_var, 'saveDir');

%% FIGURE: normalized 0 to max temp tuning curves for the five essential behaviors
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing', 'jump'};
[colors, param_names] = getPaperColors(paramList); % colors determined by universal color scheme

ntemps = length(data.tempbin.temps);
nParams = length(paramList);
TC = struct;
for i = 1:nParams
    [TC.(paramList{i}).avg, TC.(paramList{i}).std] = deal(nan(ntemps,2)); % cooling then warming for columns
end

% preformat the z-score for the data that will be extracted and used: 
Zdata = [];
Zdata(:,1) = (sum(data.CI,2)./num.trials); % courtship index
for i = 2:nParams
    Zdata(:,i) = sum(sum(data.(paramList{i}),2),3)./(num.trials*2);
end
Zdata(isnan(Zdata)) = 0;


for t = 1:ntemps
    for type = 1:2 % cooling then warming
        switch type 
            case 1 
                ROI = data.tempbin.cooling(:,t);
            case 2
                ROI = data.tempbin.warming(:,t);
        end
        
        for i = 1:nParams % courtship, sleep, food, escape
            y = Zdata(ROI,i);
            y_avg = mean(y,'omitnan');
            TC.(paramList{i}).avg(t,type) = y_avg;
            if ~isempty(y)
                y_err = std(y, 0,1,'omitnan');
                TC.(paramList{i}).std(t,type) = y_err;
            end
        end
    end
end

% ------------ Plot --------------
plotErr = false;
sSpan = 8;
r = 1;
c = 2; 
LW = 3;

x = data.tempbin.temps';

fig = getfig('', 0,[ 868 806]);
for type = 1:2
    subplot(r,c,type); hold on
    for param = 1: nParams
        raw_y = smooth(TC.(paramList{param}).avg(:,type), sSpan, 'moving');
        scaleY = rescale(raw_y);
        kolor = colors(param, :);
        plot(x, scaleY, 'color', kolor, 'linewidth',LW)
    end
end

% --------- Formatting -------------
formatFig(fig, blkbgd, [r,c]);
matchAxis(fig, true);
subplot(r,c,2)
set(gca, 'ycolor', 'none')
xlabel('temperature (\circC)')
title('warming','color', foreColor,'fontname', 'Arial','FontAngle','italic')
subplot(r,c,1)
ylabel('flies (norm %)')
set(gca, 'xdir', 'reverse')
xlabel('temperature (\circC)')
title('cooling','color', foreColor,'fontname','Arial','FontAngle','italic')

% TODO: update this to work for other temp protocols
switch groupName
    case {'Berlin LTS caviar', 'Berlin LTS 35-15 No Food', 'TrpA1-Gal4 x UAS-Kir2.1 LTS Caviar'}
        for i = 1:2
            subplot(r,c,i)
            set(gca, 'xtick', 15:5:35)
            xlim([13, 37])
            temp_lims = [14.5, 35.5];
        end
end


% update y axis to fit temp region indicator
ylims = [-0.05, 1.05];
for ii = 1:2
    subplot(r,c,ii)
    ylim(ylims)
    set(gca, 'ytick', 0:0.25:1)
end

set(findall(fig, 'Type', 'axes'), 'FontSize', 22)
set(findall(fig, 'Type', 'text'), 'FontSize', 22)

% add time arrows 
for ii = 1:2
    subplot(r,c,ii)
    addTimeArrow(gca, foreColor)
end

% plot a line for the 'safe' vs 'threat' zones: 
y = rangeLine(fig, 0, false);
LW = 5;
x_less = [temp_lims(1), 25];
x_more = [25, temp_lims(2)];
subplot(r,c,2)
    plot(x_more, [y,y], 'Color', Color('darkred'), 'LineWidth',LW, 'HandleVisibility','off') % threat
    plot(x_less, [y,y], 'Color', Color('grey'), 'LineWidth',LW, 'HandleVisibility','off') % safe
subplot(r,c,1)
    plot(x_less, [y,y], 'Color', Color('darkred'), 'LineWidth',LW, 'HandleVisibility','off') % threat
    plot(x_more, [y,y], 'Color', Color('grey'), 'LineWidth',LW, 'HandleVisibility','off') % safe


% add arrows for the time at half-max? (TODO 5.1)

save_figure(fig, [saveDir 'normalized behavior for timing temp tuning curve' figcolor(blkbgd)],fig_type);

%% ANALYSIS & FIGURES: T-50 for behaviors trial independent and flex temp bins
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors
T50_save = createFolder([saveDir, 'T50/']);

paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing', 'jump'};
[colors, param_names] = getPaperColors(paramList); % colors determined by universal color scheme
nParams = length(paramList);

temp_bins = 14:2:36; % adjust this as needed
ntemps    = length(temp_bins) - 1; % number of bins
bin_cents = (temp_bins(1:end-1) + temp_bins(2:end)) / 2;  % bin centers for plotting/interpolation
threshold = 0.5;

T50 = struct;
temp_type = {'cooling', 'warming'};
num.sexes = 2;

% initialize empty structure
% CI: trials only; all others: trials x sexes
for ii = 1:nParams
    if strcmpi(paramList{ii}, 'CI')
        nCols = num.trials;
    else
        nCols = num.trials * num.sexes;
    end
    for type = 1:2
        field = temp_type{type};
        T50.(paramList{ii}).(field) = nan(nCols, 1);
        T50.(paramList{ii}).([field '_norm']) = nan(ntemps, nCols);
        T50.(paramList{ii}).([field '_raw']) = nan(ntemps, nCols);
    end
end


for trial = 1:num.trials
    for p = 1:nParams
        if strcmpi(paramList{p}, 'CI')
            raw_all = {data.CI(:, trial)};
            col_idx = {trial};
        else
            raw_all = {data.(paramList{p})(:, 1, trial), ...
                       data.(paramList{p})(:, 2, trial)};
            col_idx = {(trial-1)*num.sexes + 1, ...
                       (trial-1)*num.sexes + 2};
        end

        for s = 1:length(raw_all)
            raw = raw_all{s};
            col = col_idx{s};

            for type = 1:2
                field = temp_type{type};

                % select the correct direction index
                switch type
                    case 1 % cooling
                        dir_idx = data.tempbin.c_idx;
                    case 2 % warming
                        dir_idx = data.tempbin.h_idx;
                end

                % bin the data by temperature
                y_by_temp = nan(ntemps, 1);
                for t = 1:ntemps
                    % frames within this direction AND within this temp bin
                    bin_idx = dir_idx & ...
                              data.temp >= temp_bins(t) & ...
                              data.temp <  temp_bins(t+1);
                    y_by_temp(t) = mean(raw(bin_idx), 'omitnan');
                end

                % normalize
                trial_max = max(y_by_temp, [], 'omitnan');
                trial_min = min(y_by_temp, [], 'omitnan');
                y_norm = (y_by_temp - trial_min) ./ (trial_max - trial_min);

                % save tuning curves
                T50.(paramList{p}).([field '_norm'])(:, col) = y_norm;
                T50.(paramList{p}).([field '_raw'])(:, col)  = y_by_temp;

                % directional threshold crossing
                switch type
                    case 1 % cooling: search high temp -> low temp
                        search_norm  = y_norm(end:-1:1);
                        search_temps = bin_cents(end:-1:1);
                    case 2 % warming: search low temp -> high temp
                        search_norm  = y_norm;
                        search_temps = bin_cents;
                end

                crossing_idx = find(search_norm >= threshold, 1, 'first');
                if ~isempty(crossing_idx)
                    if crossing_idx > 1
                        x1 = search_temps(crossing_idx - 1);
                        x2 = search_temps(crossing_idx);
                        y1 = search_norm(crossing_idx - 1);
                        y2 = search_norm(crossing_idx);
                        try
                            T50.(paramList{p}).(field)(col) = interp1([y1 y2], [x1 x2], threshold);
                        catch
                            T50.(paramList{p}).(field)(col) = search_temps(crossing_idx);
                        end
                    else
                        T50.(paramList{p}).(field)(col) = search_temps(crossing_idx);
                    end
                end
            end
        end
    end
end

% summary stats across trials x sexes
T50_summary = struct;
for ii = 1:nParams
    for type = 1:2
        field = temp_type{type};
        vals = T50.(paramList{ii}).(field);
        T50_summary.(paramList{ii}).(field).mean = mean(vals, 'omitnan');
        T50_summary.(paramList{ii}).(field).sem = std(vals, 0, 'omitnan') ./ sqrt(sum(~isnan(vals)));
        T50_summary.(paramList{ii}).(field).all = vals;
    end
end

% -------------------------------------------------------------------------------
switch groupName
    case {'Berlin LTS 35-15 No Food', 'Berlin LTS caviar', 'TrpA1-Gal4 x UAS-Kir2.1 LTS Caviar'}
        % temp_lims = [14.5, 35.5];
        temp_ticks = 15:5:35;
end


% Plot the different temp tuning curves — one figure per parameter
r = 1; c = 3;
SZ = 75;
x_pad = 0.5;
temp_lims = [min(bin_cents) - x_pad, max(bin_cents) + x_pad];
occ_lims = [-0.01, 1.01];
occ_ticks = 0:0.2:1;
threat_color = Color('darkred');
safe_color = Color('grey');

if strcmp(questdlg('display behavior-by-behavior figures?'), 'Yes')

    for pp = 1:nParams
        fig = getfig('', false);
    
        faded = colors(pp,:) + (1 - colors(pp,:)) * 0.65;  % blend toward white
    
        % --- COOLING tuning curves ---
        ax1 = subplot(r, c, 1); hold on
        plot(bin_cents, T50.(paramList{pp}).cooling_norm, ...
            'Color', faded, 'LineWidth', 0.8)
        avg_norm = mean(T50.(paramList{pp}).cooling_norm, 2, 'omitnan');
        plot(bin_cents, avg_norm, 'Color', colors(pp,:), 'LineWidth', 2.5)
    
        % mean across all trials then normalize
        avg_raw  = mean(T50.(paramList{pp}).cooling_raw, 2, 'omitnan');
        raw_max  = max(avg_raw, [], 'omitnan');
        raw_min  = min(avg_raw, [], 'omitnan');
        avg_norm = (avg_raw - raw_min) ./ (raw_max - raw_min);
        plot(bin_cents, avg_norm, 'Color', colors(pp,:), 'LineWidth', 2.5, 'linestyle', ':')
    
        h_line(threshold, foreColor, '--', 1)
        xlabel('Temperature (\circC)')
        ylabel('Normalized Occupancy')
        title({param_names{pp}; 'during cooling'})
        xlim(temp_lims)
        ylim(occ_lims)
        set(ax1, 'ytick', occ_ticks, 'xtick', temp_ticks, 'XDir', 'reverse')
    
        % --- WARMING tuning curves ---
        ax2 = subplot(r, c, 2); hold on
        plot(bin_cents, T50.(paramList{pp}).warming_norm, ...
            'Color', faded, 'LineWidth', 0.8)
        avg_norm = mean(T50.(paramList{pp}).warming_norm, 2, 'omitnan');
        plot(bin_cents, avg_norm, 'Color', colors(pp,:), 'LineWidth', 2.5)
        
        % mean across all trials then normalize
        avg_raw  = mean(T50.(paramList{pp}).warming_raw, 2, 'omitnan');
        raw_max  = max(avg_raw, [], 'omitnan');
        raw_min  = min(avg_raw, [], 'omitnan');
        avg_norm = (avg_raw - raw_min) ./ (raw_max - raw_min);
        plot(bin_cents, avg_norm, 'Color', colors(pp,:), 'LineWidth', 2.5, 'linestyle', ':')
    
        h_line(threshold, foreColor, '--', 1)
        xlabel('Temperature (\circC)')
        ylabel('Normalized Occupancy')
        title({param_names{pp}; 'during warming'})
        xlim(temp_lims)
        ylim(occ_lims)
        set(ax2, 'ytick', occ_ticks, 'xtick', temp_ticks)
    
        % --- T50 scatter ---
        ax3 = subplot(r, c, 3); hold on
        all_vals = [];
        for ii = 1:2
            vals = T50.(paramList{pp}).(temp_type{ii});
            all_vals = [all_vals; vals(:)];
            scatter(ii, vals, SZ, colors(pp,:), 'filled', ...
                'XJitter', 'density', 'XJitterWidth', 0.2, ...
                'MarkerFaceAlpha', 0.6)
            m   = T50_summary.(paramList{pp}).(temp_type{ii}).mean;
            sem = T50_summary.(paramList{pp}).(temp_type{ii}).sem;
            errorbar(ii, m, sem, ...
                'Color', foreColor, 'LineWidth', 1.5, 'CapSize', 6)
            scatter(ii, m, SZ+25, foreColor, 'filled', 'square')
        end
        xlim([0.5, 2.5])
        y_range = [min(all_vals, [], 'omitnan'), max(all_vals, [], 'omitnan')];
        y_pad   = 0.05 * diff(y_range);
        ylim([y_range(1) - y_pad, y_range(2) + y_pad])
        xticks([1, 2])
        xticklabels(temp_type)
        ylabel('Temperature at 50% max (\circC)')
        title('T_{50}')
    
        % Formatting
        linkaxes([ax1, ax2], 'xy')
        
        formatFig(fig, blkbgd, [r, c]);
        
        % add arrows and temp regime markers: 
        addTimeArrow(ax1, foreColor);
        addTimeArrow(ax2, foreColor);
        axes(ax1) % activate tuning curve axes
        y = 1.005;
        xrange = xlim;
        x_less = [xrange(1), 25];
        x_more = [25, xrange(2)];
        axes(ax2)
            plot(x_more, [y,y], 'Color', threat_color, 'LineWidth',3, 'HandleVisibility','off') % threat
            plot(x_less, [y,y], 'Color', safe_color, 'LineWidth',3, 'HandleVisibility','off') % safe
        axes(ax1)
            plot(x_less, [y,y], 'Color', threat_color, 'LineWidth',3, 'HandleVisibility','off') % threat
            plot(x_more, [y,y], 'Color', safe_color, 'LineWidth',3, 'HandleVisibility','off') % safe
    
    
        save_figure(fig, [T50_save 'T50 ' param_names{pp} figcolor(blkbgd)], '-pdf', true);
    end
end
% ------------------------------------------------------------------------------------------

%% FIGURE: T50 summary scatter — all params, cooling vs warming subplots
% temp on x-axis, behaviors on y-axis, ordered by mean T50
r = 1; c = 2;
SZ = 150;
fig = getfig('', false,[1311 758]); %  1078 531
ax_handles = gobjects(1, 2);

for ii = 1:2 % cooling then warming
    % sort behaviors by mean T50
    mean_T50 = nan(1, nParams);
    for i = 1:nParams
        mean_T50(i) = T50_summary.(paramList{i}).(temp_type{ii}).mean;
    end

    if ii == 1 % cooling: earliest in time = highest temp = top of plot
        [~, sort_idx] = sort(mean_T50, 'ascend');   % ascend so rank 1 = lowest T50 = bottom
        rank_order = nParams:-1:1;                   % flip rank assignment so highest T50 plots at top
    else       % warming: earliest in time = lowest temp = bottom of plot
        [~, sort_idx] = sort(mean_T50, 'descend');
        rank_order = 1:nParams;
    end

    ax_handles(ii) = subplot(r, c, ii); hold on
    all_vals = [];

    for rank = 1:nParams
        % y = rank_y(rank); % smaller y value distances for plotting
        i = sort_idx(rank);
        vals = T50.(paramList{i}).(temp_type{ii});
        all_vals = [all_vals; vals(:)];
        scatter(vals, rank, SZ, colors(i,:), 'filled', ...
            'YJitter', 'density', 'YJitterWidth', 0.2, ...
            'MarkerFaceAlpha', 0.8)
        m   = T50_summary.(paramList{i}).(temp_type{ii}).mean;
        sem = T50_summary.(paramList{i}).(temp_type{ii}).sem;
        errorbar(m, rank, sem, 'horizontal', ...
            'Color', foreColor, 'LineWidth', 2, 'CapSize', 6)
        scatter(m, rank, SZ+50, foreColor, 'filled', 'square')
    end

    % x-axis: temperature
    xlim(temp_lims)
    set(ax_handles(ii), 'xtick', temp_ticks)
    xlabel({' ' ; 'Temperature at 50% max (\circC)'})

    % y-axis: behaviors in sorted order
    ylim([0.5, nParams + 0.5])
    yticks(1:nParams)
    yticklabels(param_names(sort_idx))

    % flip x for cooling
    if ii == 1
       set(ax_handles(ii), XDir='reverse')
    end

    title(temp_type{ii})

    % threat/safe color bar at top of axes
    y = nParams + 0.45;   % just inside top of y-axis
    xrange = temp_lims;
    x_less = [xrange(1), 25];
    x_more = [25, xrange(2)];
    switch ii 
        case 2  % warming
        plot(x_less, [y, y], 'Color', safe_color, 'LineWidth', 3, 'HandleVisibility', 'off', 'Clipping', 'off')
        plot(x_more, [y, y], 'Color', threat_color, 'LineWidth', 3, 'HandleVisibility', 'off', 'Clipping', 'off')
        case 1 % cooling 
        plot(x_more, [y, y], 'Color', safe_color, 'LineWidth', 3, 'HandleVisibility', 'off', 'Clipping', 'off')
        plot(x_less, [y, y], 'Color', threat_color, 'LineWidth', 3, 'HandleVisibility', 'off', 'Clipping', 'off')
    end
end

% sgtitle('T_{50} across behaviors', 'Color', foreColor, 'FontSize', 20)
formatFig(fig, blkbgd, [r, c]);

set(findall(fig, 'Type', 'axes'), 'FontSize', 20)        % tick labels
set(findall(fig, 'Type', 'text'), 'FontSize', 20)        % all text/titles/labels

% add time arrow after formatFig
for ii = 1:2
    addTimeArrow(ax_handles(ii), foreColor, -0.09);
end

save_figure(fig, [T50_save 'T50 all params' figcolor(blkbgd)]);


%% FIGURE: TEMP TUNING CURVE FOR SELECTED PARAMETER

clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors
plot_threat_zone = true;

% Select the type of information to plot: 
[title_str,pName,y_dir,y_lab,nullD,scaler,dType,~,sexSep,ylimits] = ... 
    PlotParamSelectionHR('Heating and Cooling');
[colors, param_names] = getPaperColors({pName}); % colors determined by universal color scheme

switch pName
    case 'jump'
        ylimits = [0, 0.5];
end

plot_err = true;

if isempty(title_str)
    return
end
 
% ------------- Plotting parameters -------------

% select temp protocol specific plotting features
autoYLim = false; % Y LIMITS
if any(isnan(ylimits)) % in case new data is added without specs for axis limits
    autoYLim = true;
end

% TODO: update this as new protocols are added to the pipeline
if contains(groupName,'LTS') % X LIMITS
    xlimits = [13, 37];
    autoXLim = false;
    xPos = [15, 25];  % for the shaded threat region on the plot
    data_cut_off = [15, 35];
else
    autoXLim = true;
end

r = 1; % figure rows
c = 2; % heating and cooling separated plots
LW = 2; % plotting line width
FA = 0.35; % SEM shading face alpha level

% ------------- DATA AND PLOTTING -------------

% pull larger data type group
yy = data.(pName);
if sexSep % data separated per fly or per group of flies
    y_all = [squeeze(yy(:,M,:)), squeeze(yy(:,F,:))];
else 
    y_all = yy;
end
x  = data.tempbin.temps; % temp bins
nTemps = length(x); % number of temperature bins
types  = {'cooling', 'warming'};

% Extract and Plot data:
fig = getfig('',0,[998 882]); % short and fat: [1230 637]
for ii = 1:2
    subplot(r,c,ii); hold on
        Idx = data.tempbin.(types{ii});
        rawY = nan([nTemps, 2]); % first col = avg, second = sem

        % extract the tuning information across the temp bins
        for tt = 1:nTemps 
            % skip the temp bin if it's not one included in the protocol
            % (e.g. for cases where the temp slightly overshoots or
            % undershoots the value and then we get a 'read' for something
            % like the 14.5 bin when really there are a small number of
            % data points due to temp overshoot
            if x(tt)<data_cut_off(1) ||  x(tt)>data_cut_off(2)
                continue
            end
            % all the cooling data across the flies that fits this temp bin
            y = yy(Idx(:,tt),:); 

            % fill the temp bin data into the appropriate structure
            raw = mean(y,'omitnan').*scaler; % find cooling fly data 
            rawY(tt,1) = mean(raw,'omitnan');
            rawY(tt,2) = sem(raw);
        end

        % plot the data: 
         plot_error_fills(plot_err, x, rawY(:,1), rawY(:,2), colors, fig_type, FA);
         plot(x,rawY(:,1),'color', colors, 'LineWidth', LW)
end
         
% ------------ formatting ------------
formatFig(fig, blkbgd,[r,c]);
matchSubplotAxes(fig); % match x across both subplots and y across subplots
for ii = 1:2
    subplot(r,c,ii) 
    title(types{ii},'color', foreColor,'FontAngle','italic')
    xlabel('temperature (\circC)')
    if ~autoXLim;  xlim(xlimits);  end
    if ~autoYLim;  ylim(ylimits);  end
    set(gca, 'ydir', y_dir)
    % subplot specific adjustments
    if ii==1 % cooling
        set(gca, 'XDir','reverse')
        ylabel([param_names{:} ' % flies'])
    end
    if ii==2 % warming
        set(gca, 'YColor', 'none')
    end
    h_line(nullD, 'gray', '--',1)
end

% add shaded area for 'threat' temp region
set(findall(fig, 'Type', 'axes'), 'FontSize', 22)
set(findall(fig, 'Type', 'text'), 'FontSize', 22)

subplot(r,c,1) 
set(gca, 'ytick', 0:0.1:0.5)

% add time arrows 
for ii = 1:2
    subplot(r,c,ii)
    addTimeArrow(gca, foreColor)
end

switch groupName
    case {'Berlin LTS caviar', 'Berlin LTS 35-15 No Food', 'TrpA1-Gal4 x UAS-Kir2.1 LTS Caviar'}
        for i = 1:2
            subplot(r,c,i)
            set(gca, 'xtick', 15:5:35)
            xlim([13, 37])
            temp_lims = [14.5, 35.5];
        end
end

% plot a line for the 'safe' vs 'threat' zones: 
y = rangeLine(fig, 0, false);
LW = 5;
x_less = [temp_lims(1), 25];
x_more = [25, temp_lims(2)];
subplot(r,c,2)
    plot(x_more, [y,y], 'Color', Color('darkred'), 'LineWidth',LW, 'HandleVisibility','off') % threat
    plot(x_less, [y,y], 'Color', Color('grey'), 'LineWidth',LW, 'HandleVisibility','off') % safe
subplot(r,c,1)
    plot(x_less, [y,y], 'Color', Color('darkred'), 'LineWidth',LW, 'HandleVisibility','off') % threat
    plot(x_more, [y,y], 'Color', Color('grey'), 'LineWidth',LW, 'HandleVisibility','off') % safe

% Save the Figure
save_figure(fig, [saveDir title_str ' tuning curve' figcolor(blkbgd)]);