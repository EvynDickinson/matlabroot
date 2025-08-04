

%% Food vs no food trial comparisons 

% Make food to no food alignment list: 
switch expGroup
    case 'Berlin LTS 15-35 plate comparisons'
        foodPairs = [1,3; 2,4];  
    case 'Berlin F LRR 25-17 plate comparisons'
        foodPairs = [1,3; 2,4];  
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

%% Occupancy 'null' distribution for no food trials
clearvars('-except',initial_vars{:})

% Plot out the quadrant data:

np = 2; %null-pair idx
exp = 1; % active trial
% find the quads with the highest and lowest occupancy over the course of
% the experiment:

sSpan = 180;

dummy = [];
plotData = [];
for i = 1:4
    a = grouped(np).fullquad.(quadOrder{i}).all;
    dummy(i,:) = sum(a,1,'omitnan');
    plotData(:,:,i) = a;
end
% find min and max occupancy quadrants: 
[~, lowerIDX] = min(dummy);
[~, upperIDX] = max(dummy);

[minOcc,maxOcc] = deal([]);
for i = 1:num.trial(np)
    minOcc(:,i) = squeeze(plotData(:,i,lowerIDX(i)));
    maxOcc(:,i) = squeeze(plotData(:,i,upperIDX(i)));
end
    y_err = smooth(mean((minOcc-maxOcc)./2,2,'omitnan'),sSpan,'moving');
    y_avg = smooth(mean([minOcc,maxOcc],2,'omitnan'),sSpan, 'moving');

fig = getfig('',1); 
    hold on
    % plot the null distribution data
    y_err = smooth(mean((minOcc-maxOcc)./2,2,'omitnan'),sSpan,'moving');
    y_avg = smooth(mean([minOcc,maxOcc],2,'omitnan'),sSpan, 'moving');
    kolor = grouped(np).color;
    time = grouped(np).time;
    y1 = smooth(mean(minOcc,2),sSpan, 'moving');
    y2 = smooth(mean(maxOcc,2),sSpan, 'moving');
    plot_error_fills(true, time, y_avg,y_err,kolor,fig_type,0.5);
    % plot(time,y1,'color',kolor)
    % plot(time,y2,'color',kolor)
    % plot the paired food trial on top:
    x = grouped(exp).time;
    kolor = grouped(exp).color;
    y = grouped(exp).fullquad.food.avg;
    y_err = grouped(exp).fullquad.food.std./sqrt(num.trial(exp));
    plot_error_fills(true, x, y, y_err, kolor, fig_type);
    plot(x,y,'color',kolor,'linewidth', 1)
    % formatting
    xlabel('time (min)')
    ylabel('food quadrant occupancy (%)')
    xlim([0 700])
    formatFig(fig,false);
    
save_figure(fig,[figDir, 'full quad occ over time'],fig_type);






















