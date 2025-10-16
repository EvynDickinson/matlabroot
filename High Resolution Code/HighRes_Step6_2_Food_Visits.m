

initial_var{end+1} = 'FV';
initial_var{end+1} = 'maxTime';
initial_var{end+1} = 'fps';



%% Extract the periods for food visits
clearvars('-except',initial_var{:})

FV = struct; % food visits (this will later be added to the DATA structure)
fps =  fly(M).fps;

% how many frames can be skipped before a sustained fly on food period is ended
frameDropAllowance = ceil((1/3) * fps);  % 1/3 of a second

% find instances of the flies on the food
% test for a single fly first: 
for trial = 1:num.trials
    for sex = 1:2
    
        % find periods of sustained time on the food aka, when there are long
        % periods of adjactent frames with flies on the food 
        maxTime = 873000; %588000; % drops the last up and down to have even sampling and account for food time-dependency
        frames = find(data.FlyOnFood(:,sex,trial));
        frames(frames>maxTime) = [];
        onFood = diff(frames)<=frameDropAllowance; % allow for a small gap in frames (dropped etc) [logical about 'frames']
        
        idx = find(onFood==0); % locations where a running list of flies on food frames exceed the max skip allowance
        frame_loc_stop = [frames(idx); frames(end)];
        frame_loc_start = [frames(1); frames(idx+1)];
        
        onfoodROI = [frame_loc_start, frame_loc_stop]; % frame indexes of when the fly started on the food and left the food
        nOnFood = size(onfoodROI,1); % how many times is the fly on the food total in the experiment
        onFoodDuration = (diff(onfoodROI,1,2))/fps; % duration of time (s) fly spent on the food

        % save data into the FV struct
        FV(trial,sex).ROI = onfoodROI;
        FV(trial,sex).nROI = nOnFood;
        FV(trial,sex).duration = onFoodDuration;

        % WORKING HERE 
        % **** subdivide the protocol into the four different types of temp regimes:
        % temps approaching 25 from each direction **** 

        % subdivide by heating/cooling temperature regime first (then subdivide later)
        % if start of food visit is during a temp regime it counts probably
        % could add a buffer to the turn points later, but this will suffice for now
        temp_regimes = {'cooling', 'warming', 'hold'};
        temp_types = {'hot', 'cold'};
       
        hot_temp = find(data.temperature(:,trial)>25); % when is the temp warmer than neutral
        cold_temp = find(data.temperature(:,trial)<25); % when is the temp cooler than neutral
        
        for tt = 1:3 % cooling, warming, holds
            type_str = temp_regimes{tt};
            % find which food visits start in each temp regime
            r_locs =  ismember(frame_loc_start, find(data.(type_str)));  % which locations match the temp change profile
            h_locs = ismember(onfoodROI(r_locs,1), hot_temp); %logical for hot temps + this trial type
            c_locs = ismember(onfoodROI(r_locs,1), cold_temp); %logical for cold temps + this trial type

            % find overlap between heating/cooling and each temp type
            if tt<=2 % dont do this for holds...
                str = ['hot_and_' type_str];
                locs = h_locs;
                FV(trial,sex).(str).ROI = onfoodROI(locs,:);
                FV(trial,sex).(str).nROI = size(FV(trial,sex).(str).ROI,1);
                FV(trial,sex).(str).duration = onFoodDuration(locs);
                str = ['cold_and_' type_str];
                locs = c_locs;
                FV(trial,sex).(str).ROI = onfoodROI(locs,:);
                FV(trial,sex).(str).nROI = size(FV(trial,sex).(str).ROI,1);
                FV(trial,sex).(str).duration = onFoodDuration(locs);
            else
                % sort into groups for the temp regime types: 
                locs = r_locs;
                FV(trial,sex).(type_str).ROI = onfoodROI(locs,:);
                FV(trial,sex).(type_str).nROI = size(FV(trial,sex).(type_str).ROI,1);
                FV(trial,sex).(type_str).duration = onFoodDuration(locs);
            end
        end
    end
end

%% TODO: group trends across the flies for the full experiment to see how they trend
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors


% Giant histogram across the two different fly sexes
fig = getfig('',1); hold on
% histogram
for sex = 1:2
    plotData = [];
    for trial = 1:num.trials
        plotData = [plotData; FV(trial, sex).duration];
    end
    max(plotData)
    histogram(plotData,'FaceColor',data.color(sex,:),'BinEdges',0:2:120,'FaceAlpha',0.6)
end
% formatting the figure
set(gca, 'YScale','log')
xlabel('food visit duration (s)')
ylabel('instances (#)')
formatFig(fig, blkbgd);

save_figure(fig, [figDir 'food visit duration histogram'],fig_type)

% Cumulative visit duration distribution -- find the 95% number and mean
% for each of the flies in the trial
% does the 95% number or the mean change across temperature regimes?
% scatter plot of the 95% cumulative distribution time for food visits --
% could ask if this changes (and the mean) for each temp regime...
r = 1;
c = 2;
offset = 0.3;
LW = 1;
sz = 50;
fig = getfig('',1);

for sex = 1:2
    kolor = data.color(sex,:);
    plotData = [];
    for trial = 1:num.trials
        y = FV(trial, sex).duration;
        [f,x] = ecdf(y);
        % find closest to 95%
        all_idx = find(f<=0.95);
        loc = all_idx(end); % max bin closest but less than 95%
        visit_dur = x(loc); % duration of visits that encapsulates 95% of visits
        % find mean visit duration
        mean_dur = mean(y,'omitnan');
        plotData = [plotData; visit_dur, mean_dur];
    end
    mean_dur = mean(plotData(:,2),'omitnan');
    mean_CDF = mean(plotData(:,1),'omitnan');
    x = shuffle_data(linspace(sex-offset+0.1, sex+offset-0.1, size(plotData,1)));
    
    % plot 95% cumulative dist. data
    subplot(r, c, 1)
    hold on
    scatter(x, plotData(:,1),sz,kolor,'filled')
    plot([sex-offset, sex+offset], [mean_CDF, mean_CDF],"Color",kolor, 'linewidth', LW)
    % plot avg visit duration data
    subplot(r, c, 2)
    hold on
    scatter(x, plotData(:,2),sz,kolor,'filled')
    plot([sex-offset, sex+offset], [mean_dur, mean_dur],"Color",kolor, 'linewidth', LW)
end
% formatting
formatFig(fig, blkbgd, [r c]);
subplot(r, c, 1)
    title('95% CDF','color', foreColor)
    ylabel('food visit duration (s) capturing 95%')
    set(gca, 'xcolor', 'none')
    xlim([0, 2.5])
subplot(r, c, 2)
    title('Mean','color', foreColor)
    ylabel('avg food visit duration (s)')
    set(gca, 'xcolor', 'none')
    xlim([0, 2.5])

save_figure(fig, [figDir 'food visit duration and CDF'],fig_type)

% check significance?
% quick t-test

%% Duration of food visits for each temperature regime...hot&warming, cold&warming etc.
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

r = 1;
c = 2;
offset = 0.3;
LW = 1;
sz = 50;
FA = 0.6; %scatter face alpha
fig = getfig('',1);
x_locs = [1, 2, 4, 5, 7]; % locations for the 5 different conditions (H&C, C&C, C&W, H&W, H)
type_str = {'hot_and_warming','hot_and_cooling', 'cold_and_cooling', 'cold_and_warming',  'hold'};
% type_str = {'hot_and_cooling', 'cold_and_cooling', 'cold_and_warming', 'hot_and_warming', 'hold'};
xlims = [min(x_locs)-1, max(x_locs)+1];

for t = 1:length(x_locs) % (H&C, C&C, C&W, H&W, H)
    str = type_str{t};
    for sex = 1:2
        kolor = data.color(sex,:);
        plotData = [];
        for trial = 1:num.trials
            y = FV(trial, sex).(str).duration;
            if ~isempty(y)
                [f,x] = ecdf(y);
                % find closest to 95%
                all_idx = find(f<=0.95);
                loc = all_idx(end); % max bin closest but less than 95%
                visit_dur = x(loc); % duration of visits that encapsulates 95% of visits
                % find mean visit duration
                mean_dur = mean(y,'omitnan');
                plotData = [plotData; visit_dur, mean_dur];
            end
        end
        mean_dur = mean(plotData(:,2),'omitnan');
        mean_CDF = mean(plotData(:,1),'omitnan');
        x = x_locs(t);
        x_all = shuffle_data(linspace(x-offset+0.1, x+offset-0.1, size(plotData,1)));
        
        % plot 95% cumulative dist. data
        subplot(r, c, 1)
        hold on
        scatter(x_all, plotData(:,1),sz,kolor,'filled', 'MarkerFaceAlpha', FA)
        plot([x-offset, x+offset], [mean_CDF, mean_CDF],"Color",kolor, 'linewidth', LW)
        % plot avg visit duration data
        subplot(r, c, 2)
        hold on
        scatter(x_all, plotData(:,2),sz, kolor,'filled', 'MarkerFaceAlpha', FA)
        plot([x-offset, x+offset], [mean_dur, mean_dur],"Color",kolor, 'linewidth', LW)
    end
end



% formatting
formatFig(fig, blkbgd, [r c]);
subplot(r, c, 1)
    title('95% CDF','color', foreColor)
    ylabel('food visit duration (s) capturing 95%')
    set(gca, 'XTick', x_locs, 'XTickLabel', strrep(type_str,'_',' '))
    xlim(xlims)
subplot(r, c, 2)
    title('Mean','color', foreColor)
    ylabel('avg food visit duration (s)')
    set(gca, 'XTick', x_locs, 'XTickLabel', strrep(type_str,'_',' '))
    xlim(xlims)

save_figure(fig, [figDir 'food visit by temp regime duration and CDF'],fig_type)

% check significance?
% quick t-test

% Plot the temperature time course for the plotted data above
fig = figure; 
    plot(data.temperature(:,1), 'color', 'w', 'linewidth', 3)
    formatFig(fig, true);
    xlim([-20000, maxTime])
    set(gca, 'xcolor', 'none')
    h_line(25, 'grey', '--', 1)
    ylabel('temperature (\circC)')
    set(gca, 'FontSize', 25)
    save_figure(fig, [figDir 'food visit temp regime temp plot'],fig_type)


%% How does visit duration change over time? [only for LTS currently]
foreColor = formattingColors(blkbgd); % get background colors

r = 5;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2:r;

x = []; y = [];
fig = getfig('',1);
subplot(r,c,sb(1).idx)
    plot(data.time, data.temp, 'color', foreColor, 'linewidth', 1.5)
    xlims = xlim;
    ylabel('(\circC)')

subplot(r,c,sb(2).idx)
    hold on
    for sex = 1:2
        for trial = 1:num.trials
            x = [x; FV(trial, sex).ROI(:,1)];
            y = [y; FV(trial, sex).duration];
        end
    end

    % plot the points and the rolling avg of the duration
    [X,I] = sort(x);
    Y = y(I);
    X = X./(fps*60);
    scatter(X, Y, 35, Color('grey'),'filled', 'MarkerFaceAlpha', 0.6)
    plot(X, smooth(Y, 500, 'moving'), 'color', foreColor, 'linewidth', 1.5)

    set(gca, 'yscale', 'log')
    ylabel('food visit duration (s)')
    xlabel('time (min)')
    r1 = corrcoef(x,y);
    xlims = [xlims,xlim];
    xlims = [min(xlims), max(xlims)];
    xlim(xlims);

formatFig(fig, blkbgd,[r,c], sb);
subplot(r,c,sb(1).idx)
    set(gca, 'xcolor', 'none')
    xlim(xlims);
    ylim([15,35])
    set(gca, 'YTick', [15,25,35])
    title(['duration-to-time pearson corr coef: ' num2str(r1(1,2))],'color', foreColor)
subplot(r,c,sb(2).idx)
    v_line(data.warming_idx./(fps*60), 'red')
    v_line(data.cooling_idx./(fps*60), 'dodgerblue')

save_figure(fig, [figDir 'food visit duration over time'],fig_type)


% Food visit duration simply by temperature not time
temp = data.temp(x);
dur = y;
% remove temps just at 25C (eg. hold data)
hCut = 25.25; % threshold above 25 for cutoff
lCut = 24.75; % threshold below 25 for cutoff
loc = temp<hCut & temp>lCut;
temp(loc) = [];
dur(loc) = [];

cool_idx = temp<=lCut;
hot_idx = temp>=hCut;
all_hot = dur(hot_idx);
all_cold = dur(cool_idx);
hot = mean(all_hot,'omitnan');
cold = mean(all_cold,'omitnan');

fig = getfig('', 1, [384 500]);
    hold on
    scatter(temp, dur, 35, Color('grey'), "filled", 'MarkerFaceAlpha', 0.6)
    % plot hot and cold avg
    plot([15, lCut], [cold, cold], 'color', Color('dodgerblue'), 'linewidth', 2)
    plot([hCut, 35], [hot, hot], 'color', Color('red'), 'linewidth', 2)
    set(gca, 'yscale', 'log')
    formatFig(fig, blkbgd);
    xlabel('temperature (\circC)')
    ylabel('duration of food visit (s)')
    xlim([12.5, 37.5])
    set(gca, 'XTick', 15:5:35)
    [h, p] = ttest2(all_hot, all_cold);
    disp(['Hot vs Cold duration P-Value ' num2str(p)])
    title(['hot vs cold: p=' num2str(p)], 'color', foreColor)

save_figure(fig, [figDir 'food visit duration over temp'],fig_type,true,false);


% plot the lines of best fit for the above 25 and below 25 data...
rmse1 = []; rmse2 = [];
for fit_val = 0:6

    cX = temp(cool_idx);
    [p1,S1] = polyfit(cX,all_cold,fit_val);
    x1 = 15:lCut;
    y1 = polyval(p1,x1);
    plot(x1,y1,'color', foreColor,'LineStyle','--')
    
    hX = temp(hot_idx);
    [p2,S2] = polyfit(hX,all_hot,fit_val);
    x2 = hCut:35;
    y2 = polyval(p2,x2);
    plot(x2,y2,'color', foreColor,'LineStyle','--')

    % % find best fitting for both sides: 
    % y1_fit = polyval(p1, cX);
    % rmse1 = [rmse1, sqrt(mean((all_cold - y1_fit).^2))];% Calculate Root Mean Square Error (RMSE)
    % 
    % y2_fit = polyval(p2, hX);
    % rmse2 = [rmse2, sqrt(mean((all_hot - y2_fit).^2))];

    % save_figure(fig, [figDir 'food visit duration over temp with best fit ' num2str(fit_val)],fig_type,1,0);
end



% 1. Create some sample data with noise
x = 0:0.1:4;
y = 1.5*x.^2 + 2*x + 1 + randn(size(x))*5;

% 2. Try different polynomial degrees
degrees = 1:5; % Test degrees from 1 to 5
figure;
hold on;
plot(x, y, 'o'); % Plot original data

best_n = -1;
min_rmse = inf;

for n = degrees
    % Perform the polynomial fit
    p = polyfit(x, y, n);
    
    % Evaluate the polynomial over the same x-range
    y_fit = polyval(p, x);
    
    % Calculate Root Mean Square Error (RMSE)
    rmse = sqrt(mean((y - y_fit).^2));
    
    % Track the best degree and RMSE
    if rmse < min_rmse
        min_rmse = rmse;
        best_n = n;
    end
    
    % Plot each polynomial fit
    plot(x, y_fit, 'DisplayName', ['Degree ', num2str(n), ', RMSE: ', num2str(rmse, 2)]);
end

hold off;
legend show;
title('Polynomial Fits of Different Degrees');
xlabel('x');
ylabel('y');
disp(['The best degree based on lowest RMSE on the training data is: ', num2str(best_n)]);












%% 



% 
% % Giant histogram across the two different fly sexes
% fig = getfig('',1); hold on
% % histogram
% for sex = 1:2
%     plotData = [];
%     for trial = 1:num.trials
%         plotData = [plotData; FV(trial, sex).duration];
%     end
%     max(plotData)
%     histogram(plotData,'FaceColor',data.color(sex,:),'BinEdges',0:2:120,'FaceAlpha',0.6)
% end
% % formatting the figure
% set(gca, 'YScale','log')
% xlabel('food visit duration (s)')
% ylabel('instances (#)')
% formatFig(fig, blkbgd);
% 
% save_figure(fig, [figDir 'food visit duration histogram'],fig_type)




%% TODO: Plot out the trajectories of the flies from before to after they visit the food
% define the time periods that want to be plotted before and after the food visit
clearvars('-except',initial_var{:})

fps =  fly(M).fps;

preD = 15; % pre period to plot in seconds
postD = preD; % post period to plot in seconds;
%convert to frames
preD = preD * fps; % 
postD = postD * fps; %


% figure; histogram(onFoodDuration/fps)


%% TODO: which of these are within certain temperature regimes? 

%% TODO: pull out the trajectory before and after the fly is on food to see
% what region they came from/go to


%% TODO: what is the inter-eating-bout frequency? How does this change or
% not change across temperature?












