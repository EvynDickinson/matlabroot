
saveDir = createFolder([figDir, 'Paper Figures/']);
initial_var = add_var(initial_var, 'saveDir');

%% FIGURE: normalized 0 to max temp tuning curves for the four essential behaviors
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing', 'jump'};
colors = getPaperColors(paramList); % colors determined by universal color scheme

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

x = data.tempbin.temps';

fig = getfig('', 0,[ 868 806]);
for type = 1:2
    subplot(r,c,type); hold on
    for param = 1: nParams
        raw_y = smooth(TC.(paramList{param}).avg(:,type), sSpan, 'moving');
        scaleY = rescale(raw_y);
        kolor = colors(param, :);
        plot(x, scaleY, 'color', kolor, 'linewidth',2)
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
    case {'Berlin LTS caviar', 'Berlin LTS 35-15 No Food'}
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

save_figure(fig, [saveDir 'normalized behavior for timing temp tuning curve'],fig_type);

%% FIGURE: scatter plot of normalized half-max timepoint for each behavior
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing', 'jump'};
colors = getPaperColors(paramList); % colors determined by universal color scheme

ntemps = length(data.tempbin.temps);
nParams = length(paramList);
TC = struct;
for i = 1:nParams
    [TC.(paramList{i}).avg, TC.(paramList{i}).std] = deal(nan(ntemps,2)); % cooling then warming for columns
end

% preformat the z-score for the data that will be extracted and used: 
% this matrix is by time

% TODO: need to zscore within an individual fly (based on temp bins perhaps), rather than across the flies, 
% so we can have an individualized measure of behavior order timing


Zdata = [];
Zdata(:,1) = (sum(data.CI,2)./num.trials); % courtship index
for i = 2:nParams
    Zdata(:,i) = sum(sum(data.(paramList{i}),2),3)./(num.trials*2);
end
Zdata(isnan(Zdata)) = 0; 

% sort the data into temperature groups rather than by time alone
for t = 1:ntemps
    for type = 1:2 % cooling then warming
        switch type 
            case 1 
                ROI = data.tempbin.cooling(:,t);
            case 2
                ROI = data.tempbin.warming(:,t);
        end
        
        for i = 1:nParams % courtship, sleep, food, escape, jumps, etc
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

% --------------------------------------------------------

% For each parameter and each trial, find the temperature at which
% occupancy crosses 50% of that trial's max — separately for cooling and warming
threshold = 0.5; % 50% of max
T50 = struct;

for i = 1:nParams
    T50.(paramList{i}).cooling = nan(num.trials, 1);
    T50.(paramList{i}).warming = nan(num.trials, 1);
end

for trial = 1:num.trials
    for i = 1:nParams

        % get raw data for this trial
        if i == 1 % CI
            raw = data.CI(:, trial);
        else
            raw = sum(data.(paramList{i})(:,:,trial), 2); % sum across sides/sexes if needed
        end

        for type = 1:2 % 1 = cooling, 2 = warming
            % get frame indices and corresponding temps for this ramp
            switch type
                case 1
                    field = 'cooling';
                case 2
                    field = 'warming';
            end

            % collect data across all temp bins for this ramp
            y_by_temp = nan(ntemps, 1);
            for t = 1:ntemps
                ROI = data.tempbin.(field)(:, t);
                y_by_temp(t) = mean(raw(ROI), 'omitnan');
            end

            % normalize to this trial's max on this ramp
            trial_max = max(y_by_temp, [], 'omitnan');
            trial_min = min(y_by_temp, [], 'omitnan');
            y_norm = (y_by_temp - trial_min) ./ (trial_max - trial_min);

            % find first crossing of threshold
            crossing_idx = find(y_norm >= threshold, 1, 'first');
            if ~isempty(crossing_idx)
                % interpolate for sub-bin precision
                if crossing_idx > 1
                    x1 = data.tempbin.temps(crossing_idx - 1);
                    x2 = data.tempbin.temps(crossing_idx);
                    y1 = y_norm(crossing_idx - 1);
                    y2 = y_norm(crossing_idx);
                    T50.(paramList{i}).(field)(trial) = interp1([y1 y2], [x1 x2], threshold);
                else
                    T50.(paramList{i}).(field)(trial) = data.tempbin.temps(crossing_idx);
                end
            end
        end
    end
end

% summary stats across trials
T50_summary = struct;
for i = 1:nParams
    for field = {'cooling', 'warming'}
        vals = T50.(paramList{i}).(field{:});
        T50_summary.(paramList{i}).(field{:}).mean = mean(vals, 'omitnan');
        T50_summary.(paramList{i}).(field{:}).sem  = std(vals, 0, 'omitnan') ./ sqrt(sum(~isnan(vals)));
        T50_summary.(paramList{i}).(field{:}).all  = vals;
    end
end










