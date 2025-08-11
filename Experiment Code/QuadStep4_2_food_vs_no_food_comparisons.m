

%% Food vs no food trial comparisons 

% Make food to no food alignment list: 
switch expGroup
    case 'Berlin LTS 15-35 plate comparisons'
        foodPairs = [1,3; 2,4];  
    case 'Berlin F LRR 25-17 plate comparisons'
        foodPairs = [1,3; 2,4];  
    case 'Berlin temperature holds'
        foodPairs = [1,2; 3,4; 5,6; 7,8; 9,10; 11,12; 13,14; 15,16];
end
initial_vars{end+1} = 'foodPairs';

%% FIGURE: Time Course for single parameter -- select your metric
% TODO add flies on food to this section
clearvars('-except',initial_vars{:})

plot_err = true;
autoLim = false;
xlim_auto = false; % change the time range for the x axis
% time_limits = [100,700]; % time limits if manual control over x-axis range
time_limits = [0, 400];
nMax =  num.exp; 
plot_high_null = true; % plot the low or high null occupancy for empty trials
foreColor = formattingColors(blkbgd); %get background colors

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,quad_regions] = PlotParamSelection(true);
switch questdlg('Plot error?','','True','False', 'Cancel','True')
    case 'True'
        plot_err = true;
    case 'False'
        plot_err = false;
    case 'Cancel'
        return
    case ''
        return
end
if isempty(title_str)
    return
end
fig_dir = [saveDir, dir_end];
% set figure folder
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; % dependent var timecourse
sb(3).idx = 3:c:r*c; %dependent var temp tuning curve

lw = 0.5;
hLW = 1.5;
buff = 0.5; % tuning curve buffer
sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true);
for i = num.exp:-1:1
    % general parameters:
    x = grouped(i).time;
    if any(foodPairs(:,2)==i) % if this trial is an empty trial
        FP_idx = foodPairs(foodPairs(:,2)==i,1); % find it's paired food trial for pulling color
        ls = '-';
        highlight = true;
        nP = 2; % plot both high and low occ
        LW = lw;
    else
        FP_idx = i;
        ls = '-';
        highlight = false;
        nP = 1; % plot only food occupancy
        LW = lw + hLW;
    end
    kolor = grouped(FP_idx).color; 

    % temperature time course 
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor,'LineStyle',ls)

   % selected parameter time course
    subplot(r,c,sb(2).idx); hold on
    for nn = 1:nP % for number of subfields (eg. 2 low/high or 1 food)
        
        % determine if this is a food or null trial
        if data(i).emptytrial % empty trial, so use high occupancy null instead of random food assignment
           switch nn
               case 1 % high occupancy null
                    subfield = 'high'; 
                    hColor = foreColor;
               case 2 % low occupancy null
                   subfield = 'low'; 
                   hColor = foreColor;
           end
        else 
            subfield = 'food';
        end
        % pull the subfield data structure from the grouped data set
        if quad_regions % sub regions (requires '.food' or '.low' extension etc)
            yy = grouped(i).(pName).(subfield);
        else
            yy = grouped(i).(pName); % no subregions in the metric (e.g., ring)
        end
        % plot the data lines based on the type (e.g. single trial or avg)
        switch dType
            case 1 % single trial lines
                for trial = 1:num.trial(i)
                    y = smooth(yy.all(:,trial),sSpan, 'moving')*scaler;
                    if highlight
                        plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle',ls)
                    end
                    plot(x,y,'LineWidth',LW,'Color',kolor,'LineStyle',ls)
                end
            case {2, 3} % avg line (for full timecourse -- this gets separated by H/
                    y = smooth(yy.avg,sSpan, 'moving')*scaler;
                    if highlight
                        plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle',ls)
                    end
                    plot(x,y,'LineWidth',LW,'Color',kolor, 'LineStyle',ls)
        end
    end

    % temp vs dependent variable tuning curve
    subplot(r,c,sb(3).idx); hold on
     for nn = 1:nP % for number of subfields (eg. 2 low/high or 1 food)
        % determine if this is a food or null trial
        if data(i).emptytrial % empty trial, so use high occupancy null instead of random food assignment
           if nn==1
                subfield = 'high'; 
                hColor = foreColor;
           elseif nn==2
               subfield = 'low'; 
               hColor = foreColor;
           end
        else 
            subfield = 'food';
        end
        % pull the subfield data structure from the grouped data set
        if quad_regions % sub regions (requires '.food' or '.low' extension etc)
            yy = grouped(i).(pName).(subfield);
        else
            yy = grouped(i).(pName); % no subregions in the metric (e.g., ring)
        end

         switch dType
             case 1 % single trial lines
                for trial = 1:num.trial(i)
                    x = yy.temps;
                    rawY = [grouped(i).(pName).increasing.raw(:,trial),grouped(i).(pName).decreasing.raw(:,trial)];
                    y = mean(rawY,2,'omitnan')*scaler;
                    if highlight
                        plot(x,y,'LineWidth',LW+hLW,'Color',hColor)
                    end
                    plot(x,y,'color',kolor,'linewidth',LW + buff)
                end
    
             case 2 % avg lines (combined heating and cooling)
                x = yy.temps;
                rawY = [yy.increasing.raw,yy.decreasing.raw];
                y = mean(rawY,2,'omitnan')*scaler;
                y_err = (std(rawY,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle',ls)
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'LineStyle',ls)
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    
             case 3 % separated heating and cooling
                x = yy.temps;
                YC = yy.decreasing.raw;
                YH = yy.increasing.raw;
                % cooling
                y = mean(YC,2,'omitnan')*scaler;
                y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle','--')
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'linestyle', '--')
                % heating
                y = mean(YH,2,'omitnan')*scaler;
                y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle','-')
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'linestyle', '-')          
         end
         dataString{i} = grouped(i).name;
     end
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",'none')
curr_lims = xlim;

% distance
subplot(r,c,sb(2).idx)
ylabel(y_lab)
xlabel('time (min)')
set(gca,'ydir',y_dir)
curr_lims = [curr_lims; xlim];

% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel(y_lab)
xlabel('temp (\circC)')
% if ~autoLim
%     ylim(dt_lim)
% end
h_line(nullD,'grey',':',2) %36.2
set(gca,'ydir',y_dir)

% align xlimits on the two timecourse plots
if xlim_auto
    time_limits = [min(curr_lims(:)), max(curr_lims(:))]; %#ok<UNRCH>
end
subplot(r,c,sb(1).idx)
set(gca, 'xlim', time_limits)
subplot(r,c,sb(2).idx)
set(gca, 'xlim', time_limits)

% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

% null string notation (high or low occupancy)
% if plot_high_null && any([data(:).emptytrial])
%     null_str = 'null high';
% elseif ~plot_high_null && any([data(:).emptytrial])
%     null_str = 'null low';
% else
%     null_str = '';
% end
%  % ' ' null_str
% save figure
save_figure(fig,[fig_dir 'Timecourse summary ' title_str],fig_type);


%% FIGURE: testing how many of the empty trials have spatial biases?
clearvars('-except',initial_vars{:})

fig = getfig('',1); hold on
    high_occ = [];
    for i = 1:size(foodPairs,1)
        exp = foodPairs(i,2);
        high_occ = [high_occ; grouped(exp).occ_idx(:,2)]; 
    end
    h = histogram(high_occ);
    h.FaceColor = Color('grey');
    h.FaceAlpha = 0.8;
    % formatting
    formatFig(fig, blkbgd);
    xlabel('Well location')
    ylabel('exp count')
    set(gca, 'xtick', 1:4)
    title([expGroup ' empty trials'], 'color', 'w')
save_figure(fig,[figDir, 'null well distribution'],fig_type);

fig = getfig('',1);    hold on
   % both
    high_occ = [];
    for i = 1:size(foodPairs,1)
        exp = foodPairs(i,1);
        high_occ = [high_occ; grouped(exp).occ_idx(:,2)]; 
    end
    h = histogram(high_occ);
    h.FaceColor = Color('grey');
    h.FaceAlpha = 0.8;
    % % formatting
    formatFig(fig, blkbgd);
    xlabel('Well location')
    ylabel('exp count')
    set(gca, 'xtick', 1:4)
    title([expGroup ' food trials'], 'color', 'w')
save_figure(fig,[figDir, 'food well distribution'],fig_type);

% Comparison of food wells to highest occupancy well over the full experiment: 

fig = getfig('',1); hold on
    % percent of food trials with the highest occ being the food well
    trial_per = nan([size(foodPairs,1),1]);
    for i = 1:size(foodPairs,1) %food trials only
        exp = foodPairs(i,1);
        high_occ = grouped(exp).occ_idx(:,2);
        foodWell = data(exp).T.foodLoc;
        alignedN = sum((high_occ-foodWell)==0);
        alignedPercent = (alignedN/size(high_occ,1))*100;
        trial_per(i) = alignedPercent; % percent of aligned trials this group
    end
    bar(trial_per, 'FaceColor',Color('grey'))
    trials = {expNames{foodPairs(:,1)}};
    trials = strrep(trials, '_', ' ');
    set(gca, 'xtick', 1:size(foodPairs,1),'xticklabel',trials,'XTickLabelRotation',30)
    ylabel('Percent trials with food and highest occupancy aligned')
    formatFig(fig, blkbgd);
save_figure(fig,[figDir, 'food well to high occ alignment'],fig_type);


















%% LTS STATIC TEMP HOLDS occupancy temp tuning scatter plots

%% FIGURE: [STATIC TRIALS ONLY] Scatter plot for single parameter -- select your metric
clearvars('-except',initial_vars{:})

% plot_err = true;
autoLim = false;
xlim_auto = false; % change the time range for the x axis
temp_lims = [12, 35];
nMax =  num.exp; 
plot_high_null = true; % plot the low or high null occupancy for empty trials
foreColor = formattingColors(blkbgd); %get background colors

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,quad_regions] = PlotParamSelection(false);
switch questdlg('Plot error?','','True','False', 'Cancel','True')
    case 'True'
        plot_err = true;
    case 'False'
        plot_err = false;
    case 'Cancel'
        return
    case ''
        return
end
if isempty(title_str)
    return
end
fig_dir = [saveDir, dir_end];
% set figure folder
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

% Fictive or real time? 
% TODO: (8.10) update this to work for both fictive time and time duration
FT = false;
if ~FT 
    timeROI = 1:data(1).fps * 3600 * 5; % time region to average over...5 hours
end

food_color = Color('gold');
empty_color = Color('black');
empty_type = 2; % high occupancy null data


% Plotting Parameters: 
LW = 1.5; % avg line width
eLW = 0.75; % error bar width
SZ = 35; % scatter point size 
buff = 0.75; % scatter plot buffer size
dataString = cell([1,num.exp]);
r = 1; % rows
c = 2; % columns

% FIGURE:
fig = getfig('',true);

for i = 1:num.exp
    % get temperature for the experiment group: 
    x = strrep(data(i).temp_protocol,'Hold','');
    x = str2double(strrep(x,'C',''));
    
    % determine which data to grab depending on food or no food: 
    switch data(i).emptytrial
        case true  % empty trial 
            kolor = empty_color;
            subfield = 'high'; 
        case false % food present
            kolor = food_color;
            subfield = 'food';
    end

    % pull the subfield data structure from the grouped data set
    if quad_regions % sub regions (requires '.food' or '.low' extension etc)
        yy = grouped(i).(pName).(subfield);
    else
        yy = grouped(i).(pName); % no subregions in the metric (e.g., ring)
    end
    % MAKE SWITCH HERE FOR FICTIVE OR NOT...

    for type = 1:2 % heating and cooling
        subplot(r,c,type)
        hold on
        if ~FT
            hROI = timeROI;
            cROI = timeROI;
        end
        switch type
            case 1 % heating
                ROI = hROI;
            case 2 % cooling
                ROI = cROI;
        end

        % heating data: 
        y = mean(yy.all(hROI,:),1,'omitnan');
        y_avg = mean(y);
        y_err = std(y, 0,2)/sqrt(length(y));
        X = linspace(x-buff, x+buff, length(y));
        scatter(X, y, SZ, kolor, 'filled')
        plot([x-buff, x+buff],[y_avg, y_avg], 'color', kolor, 'linewidth', LW)
        if plot_err
            errorbar(x,y_avg,y_err,'color', kolor, 'linewidth', eLW)
        end
    end
end

formatFig(fig, blkbgd,[r,c]);
subplot(r,c,1)
title('warming')
xlim(temp_lims)
subplot(r,c,2)
title('cooling')
xlim(temp_lims)






















    if any(foodPairs(:,2)==i) % if this trial is an empty trial
        FP_idx = foodPairs(foodPairs(:,2)==i,1); % find it's paired food trial for pulling color
        ls = '-';
        highlight = true;
        nP = 2; % plot both high and low occ
        LW = lw;
    else
        FP_idx = i;
        ls = '-';
        highlight = false;
        nP = 1; % plot only food occupancy
        LW = lw + hLW;
    end
    kolor = grouped(FP_idx).color; 

    % temperature time course 
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor,'LineStyle',ls)

   % selected parameter time course
    subplot(r,c,sb(2).idx); hold on
    for nn = 1:nP % for number of subfields (eg. 2 low/high or 1 food)
        
        % determine if this is a food or null trial
        if data(i).emptytrial % empty trial, so use high occupancy null instead of random food assignment
           switch nn
               case 1 % high occupancy null
                    subfield = 'high'; 
                    hColor = foreColor;
               case 2 % low occupancy null
                   subfield = 'low'; 
                   hColor = foreColor;
           end
        else 
            subfield = 'food';
        end
        % pull the subfield data structure from the grouped data set
        if quad_regions % sub regions (requires '.food' or '.low' extension etc)
            yy = grouped(i).(pName).(subfield);
        else
            yy = grouped(i).(pName); % no subregions in the metric (e.g., ring)
        end
        % plot the data lines based on the type (e.g. single trial or avg)
        switch dType
            case 1 % single trial lines
                for trial = 1:num.trial(i)
                    y = smooth(yy.all(:,trial),sSpan, 'moving')*scaler;
                    if highlight
                        plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle',ls)
                    end
                    plot(x,y,'LineWidth',LW,'Color',kolor,'LineStyle',ls)
                end
            case {2, 3} % avg line (for full timecourse -- this gets separated by H/
                    y = smooth(yy.avg,sSpan, 'moving')*scaler;
                    if highlight
                        plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle',ls)
                    end
                    plot(x,y,'LineWidth',LW,'Color',kolor, 'LineStyle',ls)
        end
    end

    % temp vs dependent variable tuning curve
    subplot(r,c,sb(3).idx); hold on
     for nn = 1:nP % for number of subfields (eg. 2 low/high or 1 food)
        % determine if this is a food or null trial
        if data(i).emptytrial % empty trial, so use high occupancy null instead of random food assignment
           if nn==1
                subfield = 'high'; 
                hColor = foreColor;
           elseif nn==2
               subfield = 'low'; 
               hColor = foreColor;
           end
        else 
            subfield = 'food';
        end
        % pull the subfield data structure from the grouped data set
        if quad_regions % sub regions (requires '.food' or '.low' extension etc)
            yy = grouped(i).(pName).(subfield);
        else
            yy = grouped(i).(pName); % no subregions in the metric (e.g., ring)
        end

         switch dType
             case 1 % single trial lines
                for trial = 1:num.trial(i)
                    x = yy.temps;
                    rawY = [grouped(i).(pName).increasing.raw(:,trial),grouped(i).(pName).decreasing.raw(:,trial)];
                    y = mean(rawY,2,'omitnan')*scaler;
                    if highlight
                        plot(x,y,'LineWidth',LW+hLW,'Color',hColor)
                    end
                    plot(x,y,'color',kolor,'linewidth',LW + buff)
                end
    
             case 2 % avg lines (combined heating and cooling)
                x = yy.temps;
                rawY = [yy.increasing.raw,yy.decreasing.raw];
                y = mean(rawY,2,'omitnan')*scaler;
                y_err = (std(rawY,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle',ls)
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'LineStyle',ls)
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    
             case 3 % separated heating and cooling
                x = yy.temps;
                YC = yy.decreasing.raw;
                YH = yy.increasing.raw;
                % cooling
                y = mean(YC,2,'omitnan')*scaler;
                y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle','--')
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'linestyle', '--')
                % heating
                y = mean(YH,2,'omitnan')*scaler;
                y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle','-')
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'linestyle', '-')          
         end
         dataString{i} = grouped(i).name;
     end
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",'none')
curr_lims = xlim;

% distance
subplot(r,c,sb(2).idx)
ylabel(y_lab)
xlabel('time (min)')
set(gca,'ydir',y_dir)
curr_lims = [curr_lims; xlim];

% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel(y_lab)
xlabel('temp (\circC)')
% if ~autoLim
%     ylim(dt_lim)
% end
h_line(nullD,'grey',':',2) %36.2
set(gca,'ydir',y_dir)

% align xlimits on the two timecourse plots
if xlim_auto
    time_limits = [min(curr_lims(:)), max(curr_lims(:))]; %#ok<UNRCH>
end
subplot(r,c,sb(1).idx)
set(gca, 'xlim', time_limits)
subplot(r,c,sb(2).idx)
set(gca, 'xlim', time_limits)

% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

% null string notation (high or low occupancy)
% if plot_high_null && any([data(:).emptytrial])
%     null_str = 'null high';
% elseif ~plot_high_null && any([data(:).emptytrial])
%     null_str = 'null low';
% else
%     null_str = '';
% end
%  % ' ' null_str
% save figure
save_figure(fig,[fig_dir 'Timecourse summary ' title_str],fig_type);



%% clearvars('-except',initial_vars{:})
[title_str, param,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);
foreColor = formattingColors(blkbgd);
r = 1;
c = 2;
autoLims = true;

timeROI = 1:data(1).fps * 3600 * 5; % time region to average over...5 hours
% add in fictive temp option here later (8/11)







plotData = [];
for exp = 1:num.exp
    if ext
        y_all = mean(grouped(exp).(param).food.all(timeROI,:),1,'omitnan');
    else
        y_all = mean(grouped(exp).(param).all(timeROI,:),1,'omitnan');
    end
    plotData = autoCat(plotData, y_all',false);
end

temp_list = [15 17 20 23 25 27 30 33]; % temps that we have temp hold data for...
y_avg = mean(plotData, 1,'omitnan');
y_sem = std(plotData,0,1,'omitnan')./sqrt(num.trial);

sz = 35;
buff = 0.4;
if autoLims
    ylims = [];
else
    ylims = [0, 90];
end
xlims = [14, 36];

kolor = foreColor;
LW = 2;

fig = getfig('',1);
% cooling cooling orientation 
subplot(r,c,1); hold on
    for exp = 1:num.exp
        temp = temp_list(exp);
        y = plotData(:,exp);
        y(isnan(y)) = [];
        x = shuffle_data(linspace(temp-buff, temp+buff, length(y)));
        scatter(x,y,sz, Color('grey'), 'filled')
        scatter(temp, y_avg(exp), 70,foreColor,'filled', 'square')
        errorbar(temp, y_avg(exp), y_sem(exp),'Color', foreColor, 'LineWidth',LW)
    end
    plot(temp_list, y_avg,'color', foreColor, 'linewidth', LW, 'linestyle', ':')
    plot_error_fills(true, temp_list, y_avg, y_sem, foreColor,fig_type);
    %formatting
    set(gca,'xdir', 'reverse')
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('static')
    if autoLims
        ylims = [ylims; ylim];
    end

% warming orientation
subplot(r,c,2); hold on
    for exp = 1:num.exp
        temp = temp_list(exp);
        y = plotData(:,exp);
        y(isnan(y)) = [];
        x = shuffle_data(linspace(temp-buff, temp+buff, length(y)));
        scatter(x,y,sz, Color('grey'), 'filled')
        scatter(temp, y_avg(exp), 70,foreColor,'filled', 'square')
        errorbar(temp, y_avg(exp), y_sem(exp),'Color', foreColor, 'LineWidth',LW)
    end
    plot(temp_list, y_avg,'color', foreColor, 'linewidth', LW, 'linestyle', ':')
    plot_error_fills(true, temp_list, y_avg, y_sem, foreColor,fig_type);
    %formatting
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('static')
    if autoLims
        ylims = [ylims; ylim];
    end
formatFig(fig, blkbgd,[r c]);
for i = 1:2
    subplot(r,c,i); hold on
    if autoLims
        ylim([min(ylims(:,1)),max(ylims(:,2))]);
    else 
        ylim(ylims)
    end
end

save_figure(fig,[figDir, param ' occ temp tuning curve scatter'],fig_type);


% figure;
% hold on
% for exp = 1:num.exp
%     plot(grouped(exp).(param).food.avg(timeROI))
% end
