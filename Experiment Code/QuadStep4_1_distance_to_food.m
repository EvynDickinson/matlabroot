


%% FIGURE: Basic over-lap of time-trials and temperature protocols w/ SPEED
clearvars('-except',initial_vars{:})
plot_err = true;

autoLim = true;
% Y limit ranges
speed_lim = [0,10]; %speed
dist_lim = [5, 30]; %distance
dt_lim = [12, 34];      %distance-temp
% dt_lim = [10, 30];        %distance-temp
nMax =  num.exp;%
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = [10,11,13,14]; %speed timecourse
sb(4).idx = 3:3:15; %binned distance alignment

LW = 0.75;
sSpan = 360; % 2 minute smoothing 
dataString = cell([1,num.exp]);

% FIGURE:  
fig = getfig('',true);
for i = 1:nMax
%     i = expOrder(ii);
    x = grouped(i).time;   
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',1,'Color',kolor)

    %distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).dist.avg,'moving',sSpan);
        % for jj = 1:50
        %     y = smooth(y,'moving',sSpan);
        % end
%         y_err = smooth(grouped(i).dist.err,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        if ~autoLim
            ylim(dist_lim)
        end

%     %speed
    subplot(r,c,sb(3).idx); hold on
        y = smooth(grouped(i).speed.avg,'moving',sSpan);
        % for jj = 1:50
        %     y = smooth(y,'moving',sSpan);
        % end
        %         y_err = smooth(grouped(i).speed.err,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        if ~autoLim
            ylim(speed_lim)
        end

    %temp dependent distance
    subplot(r,c,sb(4).idx); hold on
        x = grouped(i).dist.distavgbytemp(:,1);
        y = grouped(i).dist.distavgbytemp(:,2);
        y_err = grouped(i).dist.distavgbytemp_err(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.4);
        plot(x,y,'color',kolor,'linewidth',LW+1)

        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",backColor)
% distance
subplot(r,c,sb(2).idx)
ylabel('proximity to food (mm)')
set(gca,"XColor",backColor)
set(gca,'ydir','reverse')
% speed
subplot(r,c,sb(3).idx)
ylabel('speed (mm/s)')
xlabel('time (min)')
% temp-distance relationship
subplot(r,c,sb(4).idx)
ylabel('proximity to food (mm)')

xlabel('temp (\circC)')
if ~autoLim
    ylim(dt_lim)
end
h_line(18.1,'grey',':',1) %36.2
set(gca,'ydir','reverse')
%
legend(dataString,'textcolor', foreColor, 'location', 'northwest', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' timecourse summary'],fig_type);

%% FIGURE: Basic over-lap of time-trials and temperature protocols NO SPEED

clearvars('-except',initial_vars{:})
plot_err = true;
autoLim = true;
% Y limit ranges
dist_lim = [5,35];       %distance
dt_lim = [10,32];        %distance-temp
auto_time = true;      % automated time axis limits
time_lim = [0,400];     %time limit (x-axis)
nMax =  num.exp;%
[~,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = 3:c:r*c; %binned distance alignment

LW = 0.75;
sSpan = 360;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true);
for i = 1:nMax
%     i = expOrder(ii);
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor)

    %distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).dist.avg,'moving',sSpan);
%         y_err = smooth(grouped(i).dist.err,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        if ~autoLim
            ylim(dist_lim)
        end

    %temp dependent distance
    subplot(r,c,sb(3).idx); hold on
        x = grouped(i).dist.distavgbytemp(:,1);
        y = grouped(i).dist.distavgbytemp(:,2);
        y_err = grouped(i).dist.distavgbytemp_err(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];

        plot(x,y,'color',kolor,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",backColor)
if ~auto_time
    xlim(time_lim)
end
% distance
subplot(r,c,sb(2).idx)
ylabel('proximity to food (mm)')
xlabel('time (min)')
set(gca,'ydir','reverse')
if ~auto_time
    xlim(time_lim)
end
% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel('proximity to food (mm)')
xlabel('temp (\circC)')
if ~autoLim
    ylim(dt_lim)
end
h_line(18.1,'grey',':',1) %36.2
set(gca,'ydir','reverse')
%
% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

% % get the avg distance over the whole experiment
% roi = [1000 159935];
% [dist_mean, dist_err] = deal([]);
% for i = 1:num.exp
%     dist_all = mean(grouped(i).dist.all(roi(1):roi(2),:),1,'omitnan');
%     dist_err(i) = std(dist_all);
%     dist_mean(i) = mean(dist_all);
% end
%
% e = errorbar([17,20,23,25],dist_mean,dist_err);
% e.Marker  = 'o';
% e.Color = 'r';
% e.MarkerSize = 10;

% save figure
save_figure(fig,[saveDir expGroup ' timecourse summary no speed food only'],fig_type);

%% FIGURE: WORKING highlight specific trials within the grouped data:

%1) select the group
i = 1;
trial_list = [17:20];
colorList = {'white'};


for idx = 1:length(trial_list)
     % trial specific parameters: 
    trial = trial_list(idx);
    disp([data(i).T.Date(trial) data(i).T.Arena(trial)])

    for jj = 1:2
        if jj ==1
            kolor = Color('white');
        else 
           kolor =  grouped(i).color;
        end
        % plot distance
        subplot(r,c,sb(2).idx); hold on
        x = grouped(i).time;
        y = smooth(grouped(i).dist.all(:,trial),sSpan, 'moving');
        plot(x,y,'LineWidth',LW+2,'Color',kolor)
        if ~autoLim
            ylim(dist_lim)
        end
            
        % plot temp dependent distance
        subplot(r,c,sb(3).idx); hold on
        x = grouped(i).dist.tempList(trial,:);
        y = grouped(i).dist.tempBinned(trial,:);
        loc = isnan(y);
        x(loc) = [];  y(loc) = [];
        plot(x,y,'color',kolor,'linewidth',3)
        h = warndlg('wait for identificaiton');
        uiwait(h)
    end
end

%2) select the trials
%3) overlay the data


%% FIGURE & STATS: Hysteresis for each genotype / trial
clearvars('-except',initial_vars{:})
LW = 0.75;
buff = 0.2;
SZ = 50;
r = 1; %rows
c = 3; %columns
plot_err = false;
plotSig = true; %plot significance stars
[foreColor,backColor] = formattingColors(blkbgd);

% FIGURE:
fig = getfig('',true);
% Hystersis
subplot(r,c,1)
hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    %increasing
    x = grouped(i).increasing.temps;
    y = grouped(i).increasing.avg;
    y_err = grouped(i).increasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans
    y(loc) = []; x(loc) = []; y_err(loc) = [];
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
    plot(x,y,'LineWidth',LW+0.5,'Color',kolor,'linestyle','-')
    %decreasing
    x = grouped(i).decreasing.temps;
    y = grouped(i).decreasing.avg;
    y_err = grouped(i).decreasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans
    y(loc) = []; x(loc) = []; y_err(loc) = [];

    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
    plot(x,y,'LineWidth',LW+.5,'Color',kolor,'linestyle','--','HandleVisibility','off');

    % Names and Colors of included data
    dataString{i} = grouped(i).name;
end
subplot(r,c,1)
ylabel('proximity to food (mm)')
xlabel('temp (\circC)')
set(gca, 'ydir', 'reverse')

% Pull difference in distance heating-cooling
subplot(r,c,2)
hold on
for i = 1:num.exp
    x = repmat(grouped(i).decreasing.temps,[1,num.trial(i)]);
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    kolor = grouped(i).color;
%     plot(x,y,'color',kolor,'LineWidth',LW);
    plot(mean(x,2),mean(y,2),'color',kolor,'LineWidth',2)
end
h_line(0,foreColor,':',1)
xlabel('temp (\circC)')
ylabel('distance difference (mm)')

% Cumulative difference in proximity
subplot(r,c,3)
hold on
for ii = 1:num.exp
    i = expOrder(ii);

    kolor = grouped(i).color;
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    plotY = sum(y,1,'omitnan');
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
    scatter(x,plotY,SZ,kolor,"filled","o")
    plot([ii-buff,ii+buff],[mean(plotY),mean(plotY)],'color',foreColor,'LineWidth',2)
end
xlim([0.5,num.exp+0.5])
h_line(0,foreColor,':',1)
ylabel('cumulative difference (mm)')

formatFig(fig,blkbgd,[r,c]);
set(gca,'XTick',[],'xcolor',backColor)
xlabel('Group','color',foreColor)

% STATS: are the means of any groups different from zero?
[p, mlt, id] = deal([]);
for ii = 1:num.exp
    i = expOrder(ii);
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    plotY = sum(y,1,'omitnan');
    [~,p(ii)] = ttest(plotY);
    group_name{ii} = expNames{i};
    %multicompare
    mlt = autoCat(mlt, plotY',false);
    id = autoCat(id,i*ones(length(plotY),1),false);
end
%Bonferonni correction:
alpha = 0.05;
m = num.exp;
p_limit = alpha/m;
h = p<=p_limit;
stats_tbl = table(group_name',h',p','VariableNames',{'group','significant','p value'});
disp(stats_tbl)

% add significance stars to the figure:
if plotSig
    y_pos = rangeLine(fig,1);
    subplot(r,c,3); hold on
    for ii = 1:num.exp
        if h(ii)
            scatter(ii,y_pos,100,foreColor,'*')
        end
    end
end

% Multicompare across the groups for significance
% STATS:
% TODO -- update all other stats to reflect this vv
% determine which groups differ from each other
if num.exp>1
    [~,~,stats] = anova1(mlt(:),id(:),'off');
    alpha = 0.05; %significance level
    [c,~,~,~] = multcompare(stats,alpha,'off');
    % bonferonni multiple comparisons correction
    m = size(c,1); %number of hypotheses
    sigThreshold = alpha/m;
    %find p-values that fall under the threshold
    significantHypotheses = c(:,6)<=sigThreshold;
    fprintf('\n\nPosition hysteresis cross group comparison statistics\n\n')
    [Group1,Group2,P_Value] = deal([]);
    idx = 0;
    for i = 1:length(significantHypotheses)
        if significantHypotheses(i)
            idx = idx+1;
            Group1{idx,1} = expNames{c(i,1)};
            Group2{idx,1} = expNames{c(i,2)};
            P_Value(idx,1) = c(i,6);
        end
    end
    sig_comp = table(Group1,Group2,P_Value);
    disp(sig_comp)
end

% save figure
save_figure(fig,[saveDir expGroup ' hysteresis summary'],fig_type);

%% FIGURE: ramp by ramp hysteresis comparision
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);

% FIGURE: plot the heating vs cooling plots for each of the four ramps
% acros each of the different trials
buff = 0.3;
LW = 0.75;

fig = getfig('',true,[431, 497]); hold on
    for ii = 1:num.exp
        i = expOrder(ii);
        kolor = grouped(i).color;
        nRamps = size(mat(i).hysteresis,2);
        x = repmat(linspace(ii-buff,ii+buff,nRamps),[num.trial(i),1]);
        plot(x',mat(i).cumHist','color',kolor,'linewidth',LW,...
             'Marker','o', 'MarkerFaceColor',kolor)
    end

h_line(0,foreColor,':',1)
formatFig(fig,blkbgd);
set(gca,'XTick',1:num.exp)
xlabel('Group')
ylabel('Cumulative dist difference (cool-heat)')

% save figure
save_figure(fig,[saveDir expGroup ' ramp by ramp cumulative hysteresis'],fig_type);

%% ANALYSIS AND FIGURES: Event-aligned comparisons
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
autoLim = true;
autoSave = true;
ylimits = [-20,15]; %for manual y-limit selection
sSpan = 180;
plot_err = true;

temp_limits = [nan,nan];
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).temp_protocol);
    if tp.nHold>0 % account for protocols with no holding period
        sections = {'decreasing','increasing','holding'};
    else
        sections = {'decreasing','increasing'};
    end
    for ss = 1:length(sections)
        switch sections{ss}
            case 'increasing'
                tpBin = 'up';
                nrr = tp.nUp;
            case 'decreasing'
                tpBin = 'down';
                nrr = tp.nDown;
            case 'holding'
                tpBin = 'hold';
                nrr = tp.nHold;
        end
        [temp,temperature] = deal([]);
        duration = min(tp.(tpBin)(:,2)-tp.(tpBin)(:,1)); %get shortest ramp period
        for rr = 1:nrr % for each ramp...
            ROI = tp.(tpBin)(rr,1):tp.(tpBin)(rr,1)+duration;
            temp(:,:,rr) = grouped(i).dist.all(ROI,:);
            temperature(:,rr) = grouped(i).temp(ROI);
        end
        temp_norm = temp-mean(temp(1:10,:,:),'omitnan'); %normalize to zero distance
        temp_avg = mean(temp_norm,3);

        % add to the grouped data
        grouped(i).aligned.([sections{ss} '_avg']) = temp_avg;
        grouped(i).aligned.([sections{ss} '_norm']) = temp_norm;
        grouped(i).aligned.([sections{ss} '_all']) = temp;
        grouped(i).aligned.([sections{ss} '_SEM']) = std(temp_avg,0,2,'omitnan')/sqrt(num.trial(i));
        grouped(i).aligned.([sections{ss} '_MEAN']) = mean(temp_avg,2,'omitnan');
        grouped(i).aligned.([sections{ss} '_temperature']) = mean(temperature,2,'omitnan');
    end
    temp_limits(1) = min([temp_limits(1),tp.threshLow]);
    temp_limits(2) = max([temp_limits(2),tp.threshHigh]);
end

% FIGURES: SINGLE EXPERIMENT COMPARISON
if blkbgd==true
    s_color = {'dodgerblue','red','white'};
else
    s_color = {'dodgerblue','red','black'};
end

for i = 1:num.exp
    dispROI = 15;
    if strcmpi(data(i).temp_protocol,'linear_ramp_XF_25-17') %shorter than 15mins
        dispROI = 7.5;
    end
    duration = ceil(dispROI*f2m);
    x = linspace(0,dispROI,duration+1);
    LW = 1.5;

    fig = getfig('',true,[428 832]);
    hold on
    for ss = 1:length(sections)
        y = grouped(i).aligned.([sections{ss} '_MEAN'])(1:duration+1);
        y_err = grouped(i).aligned.([sections{ss} '_SEM'])(1:duration+1);

        y = smooth(y,sSpan,'moving');
        y_err = smooth(y_err, sSpan,'moving');
        plot_error_fills(plot_err, x, y, y_err, Color(s_color{ss}),  fig_type, 0.4);
        plot(x, y,'color',Color(s_color{ss}),'LineWidth',LW)
    end
    xlabel('time (min)')
    set(gca,'ydir', 'reverse')
    ylabel('proximity to food (mm)')
    if ~autoLim
        ylim(ylimits)
    end
    formatFig(fig,blkbgd);
    title([grouped(i).name],'Color',foreColor,'FontSize',12,'FontName','times')

    save_figure(fig,[saveDir grouped(i).name...
                ' event aligned distance -duration ' num2str(dispROI) ' min'],...
                fig_type,autoSave);
end

% Select duration of time-since-event to plot
durList = {'all','5','10','15','20','30','50','other'};
selIdx = listdlg('ListString', durList,'PromptString','Display duration (min)?','ListSize',[150,150]);
if isempty(selIdx); return; end

switch durList{selIdx}
    case 'all'
        dispROI = 0;
    case 'other'
        dispROI = str2double(cell2mat(inputdlg('Visible duration?')));
    case {'5','10','15','20','30','50'}
        dispROI = str2double(durList{selIdx});
end

% FIGURE: CROSS EXPERIMENT COMPARISION (WITHIN HEATING,COOLING,HOLDING)
% Parameters
r = 5; c = 3;
sb(1).idx = 1; sb(2).idx = 2; sb(3).idx = 3; % temperature ramps
sb(4).idx = 4:3:(r*c); sb(5).idx = 5:3:(r*c);  sb(6).idx = 6:3:(r*c);
LW = 1.5;
SEM_shading = false;
sSpan = 1;

% Plotting
fig = getfig('',false,[1265 542]);
for i = 1:num.exp
    for ss = 1:length(sections)
        % Get the appropriate time and ROI:
        y = grouped(i).aligned.([sections{ss} '_temperature']);
        maxTime = length(y)/f2m;
        if dispROI==0 || dispROI>maxTime
            ROI = 1:length(y);
            duration = length(y);
            timeLength = duration/f2m;
            x = linspace(0,timeLength,duration);
        else
            duration = ceil(dispROI*f2m);
            ROI = (1:duration+1);
            timeLength = dispROI;
            x = linspace(0,timeLength,duration+1);
        end
     % -- temp ramp --
        subplot(r,c,sb(ss).idx); hold on
        y = grouped(i).aligned.([sections{ss} '_temperature'])(ROI);

        plot(x,y,'color', grouped(i).color,'linewidth',LW)

    % -- proximity --
        subplot(r,c,sb(ss+3).idx); hold on
        y = grouped(i).aligned.([sections{ss} '_MEAN'])(ROI);
        y_err = grouped(i).aligned.([sections{ss} '_SEM'])(ROI);
        plot_error_fills(SEM_shading, x, y, y_err, grouped(i).color,  fig_type, 0.4);
        plot(x, smooth(y,'moving',sSpan),'color',grouped(i).color,'LineWidth',LW)
    end
end

% Formatting
y_lim = [];
temp_limits = [floor(temp_limits(1)),ceil(temp_limits(2))];
formatFig(fig,blkbgd,[r,c],sb);
for ss = 1:length(sections)
    % Temperature boxes
    subplot(r,c,sb(ss).idx) %
        set(gca,'xcolor',backColor)
        ylim(temp_limits)
        set(gca,'ytick',temp_limits)
        ylabel('\circC')
    % Distance plots
    subplot(r,c,sb(ss+3).idx)
        set(gca,'ydir','reverse')
        ylabel('proximity to food (mm)')
        xlabel('time (min)')
        h_line(0,foreColor,':')
    % find axes info:
    y_lim = autoCat(y_lim,ylim);
end
% Set Y-Limits for distance data
if autoLim
    ylimits = [min(min(y_lim)),max(max(y_lim))];
end
for ss = 1:length(sections)
    subplot(r,c,sb(ss+3).idx)
    ylim(ylimits)
%     ylim([-15, 14])
%     ylim([-5, 5])
%     xlim([0 150])
%     subplot(r,c,sb(ss).idx)
%     xlim([0 150])
end

save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' ...
            num2str(sSpan) ' duration ' num2str(dispROI) ' min'],fig_type);
%
% save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' num2str(sSpan) ' duration ' num2str(dispROI) ' min 2'],'-pdf');

% % FIGURE: CROSS EXPERIMENT COMPARISION (WITHIN HEATING,COOLING,HOLDING -- NO TEMP PLOTS)
% r = 1;
% c = 3;
%
% dispROI = 50;
% duration = ceil(dispROI*3*60);
% x = linspace(0,dispROI,duration+1);
% LW = 1.5;
% SEM_shading = false;
% sSpan = 1;
%
%
% fig = figure; set(fig,'pos',[1932 690 1050 438])
% for ss = 1:length(sections) %
%     subplot(r,c,ss); hold on
%     for i = 1:num.exp
%         y = grouped(i).aligned.([sections{ss} '_MEAN'])(1:duration+1);
%         plot_error_fills(SEM_shading, x, y, y_err, grouped(i).color,  fig_type, 0.4);
%         plot(x, smooth(y,'moving',sSpan),'color',grouped(i).color,'LineWidth',LW)
%     end
%     xlabel('time (min)')
%     ylabel('distance from food (mm)')
%     title(sections{ss})
%     ylim(ylimits)
% end
% formatFig(fig,true,[r,c]);

% save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' ...
%             num2str(sSpan) ' duration ' num2str(dispROI) ' min'],'-png',true);

%% FIGURE: Event aligned with preceding data before event
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);

switch questdlg('Short or long display range?','','Short (15 min)','Long (50 min)','Cancel','Short (15 min)')
    case 'Short (15 min)'
%         ylimits = [-8,4];
        ylimits = [-5,2];
        dispROI = 10;
        pre_disp = 10;
    case 'Long (50 min)'
        ylimits = [-15,15];
        dispROI = 50;
        pre_disp = 20;
end


% SELECT THE TIME ROIs
timeROI = 60; % how many minutes to look at behavior after each event
preROI = pre_disp;  % how many minutes to look at behavior before each event
post_duration = ceil(timeROI*3*60);
pre_duration = ceil(preROI*3*60);
total_duration = pre_duration+post_duration;

sections = {'decreasing','increasing','holding'};
temp_limits = [nan,nan];
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).temp_protocol);
    for ss = 1:length(sections)
        switch sections{ss}
            case 'increasing'
                tpBin = 'up';
                nrr = tp.nUp; %number of events for this category
            case 'decreasing'
                tpBin = 'down';
                nrr = tp.nDown;
            case 'holding'
                tpBin = 'hold';
                nrr = tp.nHold;
        end
        [temp,temperature] = deal([]);
        for rr = 1:nrr
            ROI = tp.(tpBin)(rr,1)-pre_duration:tp.(tpBin)(rr,1)+post_duration;
            if any(ROI<=0)
                loc = (ROI<=0);
                ROI(loc) = [];
                xx = find(loc);
                MT = nan(xx(end),num.trial(i));
                temp(:,:,rr) = [MT ; grouped(i).dist.all(ROI,:)];
                temperature(:,rr) = [MT(:,1); grouped(i).temp(ROI)];
            else
                temp(:,:,rr) = grouped(i).dist.all(ROI,:);
                temperature(:,rr) = grouped(i).temp(ROI);
            end

        end

        temp_norm = temp-mean(temp(pre_duration-5:pre_duration+5,:,:),'omitnan'); %normalize to zero distance
        temp_avg = mean(temp_norm,3,'omitnan');

        % add to the grouped data
        grouped(i).ext_aligned.([sections{ss} '_avg']) = temp_avg;
        grouped(i).ext_aligned.([sections{ss} '_norm']) = temp_norm;
        grouped(i).ext_aligned.([sections{ss} '_all']) = temp;
        grouped(i).ext_aligned.([sections{ss} '_SEM']) = std(temp_avg,0,2,'omitnan')/sqrt(num.trial(i));
        grouped(i).ext_aligned.([sections{ss} '_MEAN']) = mean(temp_avg,2,'omitnan');
        grouped(i).ext_aligned.([sections{ss} '_temperature']) = mean(temperature,2,'omitnan');
    end
    temp_limits(1) = min([temp_limits(1),tp.threshLow]);
    temp_limits(2) = max([temp_limits(2),tp.threshHigh]);
    grouped(i).ext_aligned.postTime = timeROI;
    grouped(i).ext_aligned.preTime = preROI ;
    grouped(i).ext_aligned.post_duration = post_duration;
    grouped(i).ext_aligned.pre_duration = pre_duration;
end


% FIGURE: CROSS EXPERIMENT COMPARISION (WITHIN HEATING,COOLING,HOLDING)
r = 5;
c = 3;
sb(1).idx = 1; sb(2).idx = 2; sb(3).idx = 3; % temperature ramps
sb(4).idx = 4:3:(r*c); sb(5).idx = 5:3:(r*c);  sb(6).idx = 6:3:(r*c);


post_dur = ceil(dispROI*3*60);
pre_dur = ceil(pre_disp*3*60);
x = [linspace(-pre_disp,0,pre_dur), linspace(0,dispROI,post_dur+1)];
LW = 1.5;
SEM_shading = false;
sSpan = 30;
xlimit = [-pre_disp,dispROI];

fig = getfig('',true);
for ss = 1:length(sections) %
    % temp ramp
    subplot(r,c,sb(ss).idx); hold on
    for i = 1:num.exp
        y = grouped(i).ext_aligned.([sections{ss} '_temperature'])(1:length(x));
        plot(x,y,'color', grouped(i).color,'linewidth',LW)
    end
    ylabel('\circC')
    ylim(temp_limits) % TODO 1/4[tp.threshLow,tp.threshHigh]
    v_line(0,foreColor,':',1)
    xlim(xlimit)
    % event-aligned plot
    subplot(r,c,sb(ss+3).idx); hold on
    for i = 1:num.exp
        y = grouped(i).ext_aligned.([sections{ss} '_MEAN'])(1:length(x));
        plot_error_fills(SEM_shading, x, y, y_err, grouped(i).color,  fig_type, 0.4);
        plot(x, smooth(y,sSpan,'moving'),'color',grouped(i).color,'LineWidth',LW)
    end
    xlabel('time (min)')
    ylabel('proximity to food (mm)')
%     title(sections{ss})
    set(gca,'ydir','reverse')
    ylim(ylimits)
    v_line(0,foreColor,':',1)
    xlim(xlimit)
end
formatFig(fig,blkbgd,[r,c],sb);
for ss = 1:length(sections) %
    subplot(r,c,sb(ss).idx)
    set(gca,'xcolor',backColor)
end

save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' ...
            num2str(sSpan) ' duration ' num2str(dispROI) ' min with pretime'],fig_type,false);

clearvars('-except',initial_vars{:})
autoSave = true;
[foreColor,backColor] = formattingColors(blkbgd);
SZ = 40;
types = {'hold','down','up'};

for i = 1:num.exp
  fig = getfig('',false,[1368 557]);
     for trial = 1:num.trial(i)
        wells = mat(i).position(trial).wells;
        wellLoc = mat(i).position(trial).wellLoc;
        kolor = mat(i).position(trial).cMap;

        for tt = 1:3
            subplot(1,3,tt); hold on
            x = mat(i).position(trial).(types{tt}).data(:,1);
            y = mat(i).position(trial).(types{tt}).data(:,2);

            % Make food well the origin
            x_offset = wells(1,wellLoc);
            y_offset = wells(2,wellLoc);
            wells_x = wells(1,:)-x_offset;
            wells_y = wells(2,:)-y_offset;
            X = x-x_offset;
            Y = y-y_offset;

            % Rotate to correct orientation
            switch wellLoc
                case 1
                    plotData(:,1) = Y;
                    plotData(:,2) = -X;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 2
                    plotData(:,1) = X;
                    plotData(:,2) = -Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 3
                    plotData(:,1) = -Y;
                    plotData(:,2) = X;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
                case 4
                    plotData(:,1) = X;
                    plotData(:,2) = Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
            end
            % PLOT
            scatter(WELLS(1:4,1),WELLS(1:4,2),SZ,foreColor,'filled')
            scatter(WELLS(wellLoc,1),WELLS(wellLoc,2),SZ,'green','filled')
            scatter(plotData(:,1),plotData(:,2),15,kolor,'filled')
        end
     end
     % Formatting
     formatFig(fig,blkbgd,[1,3]);
     for tt = 1:3
        subplot(1,3,tt)
        viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color',foreColor);
        axis square;
        axis equal;
        title(types{tt},'color',foreColor)
        set(gca,'XColor',backColor,'YColor',backColor);
        if ~strcmpi(getenv('COMPUTERNAME'),'EVYNPC')
            % set color bar information
            colorData = uint8(round(kolor.*255)); % convert color map to uint8
            colormap(colorData);
            c = colorbar('color',foreColor);
            set(c,'Ticks',[0,1],'TickLabels',[mat(i).position(trial).tempBins(1),mat(i).position(trial).tempBins(end)]);
            c.Label.String = 'Temperature (\circC)';
            c.Label.VerticalAlignment = "bottom";
        end
     end

%     title([data(i).ExpGroup],'color',foreColor)
    save_figure(fig,[saveDir expGroup grouped(i).name ' rate divided COM position'],fig_type,autoSave);

end

%% FIGURE: Normalized increasing/decreasing temp responses
% TODO: add option to plot the heating | cooling on different subplots
% TODO  1/26 invert the zscore ...
clearvars('-except',initial_vars{:})
autoSave = true;
[foreColor,backColor] = formattingColors(blkbgd);
sep_h_c = true; %separate heating and cooling ramps
plot_err = false;
LW = 1.5;
sSpan = 180;

% set up figure aligments
if sep_h_c
    r = 5; %rows
    c = 11; %columns
    sb(1).idx = 1:6; %temp timecourse
    sb(2).idx = [12:(12+5),23:(23+5),34:(34+5),45:(45+5)]; %distance from food timecourse
    sb(3).idx = [8:c:r*c,9:c:r*c]; %binned distance alignment cooling
    sb(4).idx = [10:c:r*c,11:c:r*c]; %binned distance alignment heating
    sb(5).idx = 7:c:r*c;%buffer space
else
    r = 5; %rows
    c = 3; %columns
    sb(1).idx = 1:2; %temp timecourse
    sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse
    sb(3).idx = 3:3:15; %binned distance alignment
end

% FIGURE:
fig = getfig('normalized repsonses',true);
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',LW,'Color',kolor)

    %zscore distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).dist.zscore.avg,'moving',sSpan);
%         y_err = smooth(grouped(i).speed.err,'moving',sSpan);
%         plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW,'Color',kolor)

    %cooling dependent distance
    subplot(r,c,sb(3).idx); hold on

        %decreasing
        x = grouped(i).decreasing.temps;
        y = grouped(i).decreasing.zscore.avg;
        y_err = grouped(i).decreasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];

        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW,'Color',kolor,'linestyle','--','HandleVisibility','off');

    if sep_h_c
%         subplot(r,c,sb(5).idx);
        %heating dependent distance
        subplot(r,c,sb(4).idx); hold on
    end
        %increasing
        x = grouped(i).increasing.temps;
        y = grouped(i).increasing.zscore.avg;
        y_err = grouped(i).increasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];

        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW,'Color',kolor,'linestyle','-')

        % Names and Colors of included data
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",backColor)
% distance
subplot(r,c,sb(2).idx)
ylabel('distance zscore')
xlabel('time (min)')

% temp rate
subplot(r,c,sb(3).idx)
ylimits = ylim;
ylabel('distance (zscore)')
xlabel('temp (\circC)')
if sep_h_c
    subplot(r,c,sb(4).idx)
    xlabel('temp (\circC)')
    set(gca,'ycolor',backColor)

    %set same y scale range
    ylimits_2 = ylim;
    ymax = max([ylimits(2),ylimits_2(2)]);
    ymin = min([ylimits(1),ylimits_2(1)]);
    ylim([ymin, ymax])
    subplot(r,c,sb(3).idx)
    ylim([ymin, ymax])
end
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' zscore timecourse summary'],fig_type,autoSave);


% FIGURE: Separated heating & cooling
fig = getfig('',true,[1032 704]);
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;
    %cooling dependent distance
    subplot(1,2,1); hold on
        x = grouped(i).decreasing.temps;
        y = grouped(i).decreasing.zscore.avg;
        y_err = grouped(i).decreasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW,'Color',kolor,'HandleVisibility','off');
        title('Cooling')
        xlabel('Temperature (\circC)')
        ylabel('Distance z-score')
    %heating dependent distance
    subplot(1,2,2); hold on
        %increasing
        x = grouped(i).increasing.temps;
        y = grouped(i).increasing.zscore.avg;
        y_err = grouped(i).increasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];

        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        title('Heating')
        xlabel('Temperature (\circC)')
        ylabel('Distance z-score')
        % Names and Colors of included data
        dataString{i} = grouped(i).name;
end
% FORMATING AND LABELS
formatFig(fig,blkbgd,[1,2]);
fig = matchAxis(fig,true);
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 6)

save_figure(fig,[saveDir expGroup ' zscore heating and cooling overlay'],fig_type, autoSave);

%% FIGURE & STATS: Distance traveled over temperature variation
clearvars('-except',initial_vars{:})
autoAxis = true;
[~,backColor] = formattingColors(blkbgd);

% set up figure aligments | parameters
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse
sb(3).idx = 3:3:15; %binned distance alignment

LW = 2;
SZ = 50;
sSpan = 180;
buff = 0.2;

sizeUnit = 85; %60
figSize = [sizeUnit+(sizeUnit*num.exp), 650]; %590

% FIGURE:
fig = getfig('',true,figSize); hold on
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
    y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
%     datastats.all = autoCat(datastats.all,y,false);
%     datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
    scatter(x,y,SZ,kolor,'filled')
    plot([ii-(2.5*buff),ii+(2.5*buff)],[mean(y),mean(y)],'color', kolor,'linewidth',LW)
end
set(gca,'xtick',0:num.exp+1,'XTickLabel',[],'XColor',backColor)
xlim([0,num.exp+1])
if ~autoAxis
    ylim([0,22])
end
% xlabel('Fly Strain')
ylabel('temp-driven spatial range (mm)')
formatFig(fig,blkbgd);
set(gca,'XColor',backColor)
save_figure(fig,[saveDir expGroup ' distance modulation by temperature'],fig_type);

% STATS:
[datastats.all, datastats.id] = deal([]);
for i = 1:num.exp
    y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
    datastats.all = autoCat(datastats.all,y,false);
    datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
end

% determine which groups differ from each other
[~,~,stats] = anova1(datastats.all,datastats.id,'off');
alpha = 0.05; %significance level
[c,~,~,~] = multcompare(stats,alpha,'off');
% bonferonni multiple comparisons correction
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nTemp-driven position range statistics\n\n')
[Group1,Group2,P_Value] = deal([]);
idx = 0;
for i = 1:length(significantHypotheses)
    if significantHypotheses(i)
        idx = idx+1;
        Group1{idx,1} = expNames{c(i,1)};
        Group2{idx,1} = expNames{c(i,2)};
        P_Value(idx,1) = c(i,6);
    end
end
sig_comp = table(Group1,Group2,P_Value);
disp(sig_comp)

%% FIGURE & STATS: Min and max distance from the food (for a temperature bin)
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
y_limits = [5,35];
buff = 0.2;
r = 2;
c = 2;
SZ = 30;
LW = 1.5;

% GATHER DATA
for tt = 1:2
  switch tt
     case 1
        type = 'increasing';
     case 2
        type = 'decreasing';
  end

    % Pull the avg min|max distance
    [min_dist,max_dist]  = deal(nan([max(num.trial),num.exp]));
    for i = 1:num.exp
        for trial = 1:num.trial(i)
            y = grouped(i).(type).all(:,trial);
            % Minimum distance
            min_dist(trial,i) = min(y); %I is the temp index
            % Maximum distance
            max_dist(trial,i) = max(y);
        end
        grouped(i).(type).minDist.dist = min_dist(:,i);
        grouped(i).(type).maxDist.dist = max_dist(:,i);
    end
end


fig = getfig('',true,[631 718]);
for ii = 1:num.exp
    i = expOrder(ii); %plot in correct order
    kolor = grouped(i).color;

    % Increasing temp rate -> minimum distance to food
    subplot(r,c,1); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).increasing.minDist.dist;
        y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)

    % Increasing temp rate -> maximum distance to food
    subplot(r,c,2); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).increasing.maxDist.dist;
         y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)
    % Decreasing temp rate -> minimum distance to food
    subplot(r,c,3); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).decreasing.minDist.dist;
         y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)
    % Decreasing temp rate -> maximum distance to food
    subplot(r,c,4); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).decreasing.maxDist.dist;
         y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)
end
%labels and formatting
formatFig(fig,blkbgd,[r,c]);
for ii = 1:4
    subplot(r,c,ii)
    set(gca,'xcolor', backColor,'ydir','reverse')
    xlim([0,num.exp+1])
    ylim(y_limits)
    ylabel('food proximity (mm)')
    switch ii
        case 1
            title('Min | Heating','color',foreColor)
        case 2
            title('Max | Heating','color',foreColor)
        case 3
            title('Min | Cooling','color',foreColor)
        case 4
            title('Max | Cooling','color',foreColor)
    end
end
save_figure(fig,[saveDir expGroup ' food proximity modulation by temperature'],fig_type);


% STATS: are the distances different from each other for each condition (min/max | heat/cool)?
for ll = 1:2
    switch ll
        case 1
            distGroup = 'minDist';
        case 2
            distGroup = 'maxDist';
    end
    for tt = 1:2
        switch tt
         case 1
            type = 'increasing';
         case 2
            type = 'decreasing';
        end

        [datastats.all, datastats.id] = deal([]);
        for i = 1:num.exp
            y = grouped(i).(type).(distGroup).dist; y(isnan(y)) = [];
            datastats.all = autoCat(datastats.all,y',false);
            datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
        end
        % % determine which groups differ from each other
        [~,~,stats] = anova1(datastats.all,datastats.id,'off'); close all
        alpha = 0.05; %significance level
        [c,~,~,~] = multcompare(stats,alpha,'off');
        % bonferonni multiple comparisons correction
        m = size(c,1); %number of hypotheses
        sigThreshold = alpha/m;
        %find p-values that fall under the threshold
        significantHypotheses = c(:,6)<=sigThreshold;
        fprintf(['\n\n' type ' temp rate | ' distGroup ' comparison\n\n'])
        [Group1,Group2,P_Value] = deal([]);
        idx = 0;
        if any(significantHypotheses)
            for i = 1:length(significantHypotheses)
                if significantHypotheses(i)
                    idx = idx+1;
                    Group1{idx,1} = expNames{c(i,1)};
                    Group2{idx,1} = expNames{c(i,2)};
                    P_Value(idx,1) = c(i,6);
                end
            end
            sig_comp = table(Group1,Group2,P_Value);
            disp(sig_comp)
        else
            disp('no significance')
        end

    end
end

%% ANALYSIS AND FIGURE: temp range at minimum and max distances
% TODO: statistical comparisions within groups 1/2/23
clearvars('-except',initial_vars{:})
[~,backColor] = formattingColors(blkbgd);
autoLim = true;
yLimits = [14,28];

figUnit = 70;
figSize = [2*(figUnit+(num.exp*figUnit)), 590];

for tt = 1:2
  switch tt
     case 1
        type = 'increasing';
     case 2
        type = 'decreasing';
  end

    % Pull the avg min|max distance and temperature index for each trial
    [min_dist,max_dist,I_min,I_max]  = deal(nan([max(num.trial),num.exp]));
    for i = 1:num.exp
        for trial = 1:num.trial(i)
            y = grouped(i).(type).all(:,trial);
            % Minimum distance
            [min_dist(trial,i), I_min(trial,i)] = min(y); %I is the temp index
            % Maximum distance
            [max_dist(trial,i), I_max(trial,i)] = max(y);
        end
    end

    % For each trial -- find which temperatures have a mean distance that falls
    % within the SE of the distance for the 'min|max' distance found above
    [minimumTempRangeStart,minimumTempRangeEnd,mimimumTemp] = deal(nan([max(num.trial),num.exp]));
    [maxTempRangeStart,maxTempRangeEnd,maxTemp] = deal(nan([max(num.trial),num.exp]));
    for i = 1:num.exp
        for trial = 1:num.trial(i)
            switch type
                case 'increasing'
                    rateIDX = find(data(i).G(trial).TR.rates>0); % increasing temperature rates
                case 'decreasing'
                    rateIDX = find(data(i).G(trial).TR.rates<0); % decreasing temperature rates
            end
            if ~(length(rateIDX)==1)
                warndlg('Too many or few increasing temperature rates detected')
            end
            tempList = data(i).G(trial).TR.temps;
            y = data(i).G(trial).TR.data(:,3); % distance from food for this whole trial
            positionAvg = grouped(i).(type).all(:,trial); % avg distance at temps for trial

            % ------ minimum temp -------
            % find the SE for the distance at the minumum temperature:
            loc = (data(i).G(trial).TR.tempIdx==I_min(trial,i)) &...
                  (data(i).G(trial).TR.rateIdx==rateIDX);
            testData = y(loc);
            SE = std(testData,'omitnan');
            boundLimit = min_dist(trial,i)+SE; %since this is minumum value, the limit would only exist above
            tempLocWithinBounds = positionAvg<=boundLimit;

            % find other temps that have their distance within 1 STD of the min food distance
            upTemp = tempLocWithinBounds(I_min(trial,i):end);
            upLoc = find(diff(upTemp)<0);
            if ~isempty(upLoc)
                tempLocWithinBounds(I_min(trial,i)+upLoc(1):end) = false;
            end
            downTemp = tempLocWithinBounds(1:I_min(trial,i));
            downLoc = find(diff(downTemp)>0);
            if ~isempty(downLoc)
                tempLocWithinBounds(1:downLoc(end)) = false;
            end
            tempswithinBounds = tempList(tempLocWithinBounds);

            % Save the temp ranges for minumum distance:
            minimumTempRangeStart(trial,i) = tempswithinBounds(1);
            minimumTempRangeEnd(trial,i) = tempswithinBounds(end);
            mimimumTemp(trial,i) = tempList(I_min(trial,i));

             % ------ maximum temp -------
            % find the SE for the distance at the minumum temperature:
            loc = (data(i).G(trial).TR.tempIdx==I_max(trial,i)) & ...
                  (data(i).G(trial).TR.rateIdx==rateIDX);
            testData = y(loc);
            SE = std(testData,'omitnan');
            boundLimit = max_dist(trial,i)-SE; %since this is maximum value, the limit would only exist below
            tempLocWithinBounds = positionAvg>=boundLimit;

            % find other temps that have their distance within 1 STD of the max food distance
            upTemp = tempLocWithinBounds(I_max(trial,i):end);
            upLoc = find(diff(upTemp)<0);
            if ~isempty(upLoc)
                tempLocWithinBounds(I_max(trial,i)+upLoc(1):end) = false;
            end
            downTemp = tempLocWithinBounds(1:I_max(trial,i));
            downLoc = find(diff(downTemp)>0);
            if ~isempty(downLoc)
                tempLocWithinBounds(1:downLoc(end)) = false;
            end
            tempswithinBounds = tempList(tempLocWithinBounds);

            % Save the temp ranges for minumum distance:
            maxTempRangeStart(trial,i) = tempswithinBounds(1);
            maxTempRangeEnd(trial,i) = tempswithinBounds(end);
            maxTemp(trial,i) = tempList(I_max(trial,i));
        end
        grouped(i).(type).minDist.TempRange = [minimumTempRangeStart(:,i),minimumTempRangeEnd(:,i)];
        grouped(i).(type).minDist.Temp = mimimumTemp(:,i);
        grouped(i).(type).maxDist.TempRange = [maxTempRangeStart(:,i),maxTempRangeEnd(:,i)];
        grouped(i).(type).maxDist.Temp = maxTemp(:,i);
    end
end

% FIGURE:
buff = 0.15;
r = 1;
c = 2;
fig = getfig('',true,figSize);
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    % MAX DISTANCE TEMP RANGE
    subplot(r,c,1); hold on
        % -- increasing ---
        y1 = mean(grouped(i).increasing.maxDist.TempRange(:,1),'omitnan');
        y2 = mean(grouped(i).increasing.maxDist.TempRange(:,2),'omitnan');
        plot([ii+buff,ii+buff],[y1,y2],'color',kolor,'LineWidth',1.5)
        scatter(ii+buff,mean(grouped(i).increasing.maxDist.Temp,'omitnan'),50,kolor,'filled')
        % -- decreasing ---
        y1 = mean(grouped(i).decreasing.maxDist.TempRange(:,1),'omitnan');
        y2 = mean(grouped(i).decreasing.maxDist.TempRange(:,2),'omitnan');
        plot([ii-buff,ii-buff],[y1,y2],'color',kolor,'LineWidth',1.5,'LineStyle',':')
        scatter(ii-buff,mean(grouped(i).decreasing.maxDist.Temp,'omitnan'),50,kolor,'filled')

    % MIN DISTANCE TEMP RANGE
    subplot(r,c,2); hold on
        % -- increasing ---
        y1 = mean(grouped(i).increasing.minDist.TempRange(:,1),'omitnan');
        y2 = mean(grouped(i).increasing.minDist.TempRange(:,2),'omitnan');
        plot([ii+buff,ii+buff],[y1,y2],'color',kolor,'LineWidth',1.5)
        scatter(ii+buff,mean(grouped(i).increasing.minDist.Temp,'omitnan'),50,kolor,'filled')
        % -- decreasing ---
        y1 = mean(grouped(i).decreasing.minDist.TempRange(:,1),'omitnan');
        y2 = mean(grouped(i).decreasing.minDist.TempRange(:,2),'omitnan');
        plot([ii-buff,ii-buff],[y1,y2],'color',kolor,'LineWidth',1.5,'linestyle',':')
        scatter(ii-buff,mean(grouped(i).decreasing.minDist.Temp,'omitnan'),50,kolor,'filled')
end
formatFig(fig,blkbgd,[r,c]);
subplot(r,c,2)
    xlim([0,num.exp+1]);
    if ~autoLim
        ylim(yLimits)
    end
    set(gca,'xcolor', backColor) % 'xtick',0:num.exp+1,'XTickLabel',[])
    xlabel(' ')
    ylabel('temp (\circC) when closest to food')
subplot(r,c,1)
    xlim([0,num.exp+1]);
    if ~autoLim
        ylim(yLimits)
    end
    set(gca,'xcolor', backColor)%  set(gca,'xtick',0:num.exp+1,'XTickLabel',[])
    xlabel(' ')
    ylabel('temp (\circC) when furthest from food')
if autoLim
    matchAxis(fig,true);
end

save_figure(fig,[saveDir expGroup ' min and max ranges'],fig_type);

%% FIGURE: Temperature - distance correlation
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
corr_coef = [];
buff = 0.2;
SZ = 50;
LW = 1.5;

pixWidth = 60; % additional pixel size for image for each extra experiment group
figSize = [pixWidth + (pixWidth*num.exp),590];

% get correlation data
for i = 1:num.exp
    pooledData = [];
    % get speed / distance information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp;
        y = grouped(i).dist.all(:,trial);
        temp = [];
        temp = autoCat(temp,x,false);%temp for each trial within the exp.
        temp = autoCat(temp,y,false);%dist for each trial within the exp
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-distance correlation
        rho = corr(temp);
        corr_coef(i).all(trial) = rho(1,2);
        % save data for pooled comparison
        pooledData = [pooledData; temp];
    end
    % Pooled speed-distance correlation
    rho = corr(pooledData);
    corr_coef(i).group = rho(1,2);
end

% correlation coefficients
fig = getfig('',true,figSize);
hold on
 for ii = 1:num.exp
   i = expOrder(ii);
   kolor = grouped(i).color;
   xlow = ii-buff-0.1;
   xhigh = ii+buff+0.1;
   x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
   y = corr_coef(i).all;
   y_avg = mean(corr_coef(i).all);
   scatter(x,y,SZ,kolor,'filled')
   plot([xlow,xhigh],[corr_coef(i).group,corr_coef(i).group],'color',kolor,'linestyle',':','linewidth',LW)
   plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',LW)
 end
 xlim([0.5,num.exp+.5])
 ylabel('temp-distance corr. coef.')
 h_line(0,foreColor,':',1)
 formatFig(fig,blkbgd);
 set(gca,'xcolor',backColor)

% save figure
save_figure(fig,[saveDir expGroup ' temp distance correlation'],fig_type);

%% FIGURE: Temp-distance
% correlation ONLY during ramps

clearvars('-except',initial_vars{:})

[foreColor,backColor] = formattingColors(blkbgd);
ylimits = [-0.8,0.1];
autoLim = true;
corr_coef = [];
buff = 0.2;
SZ = 50;
LW = 1.5;

pixWidth = 60; % additional pixel size for image for each extra experiment group
figSize = [pixWidth + (pixWidth*num.exp),590];

% get correlation data
plotData = struct;
for exp = 1:num.exp
    i = expOrder(exp);
    disp(expNames{i})

    pooled_temp = [];
    pooled_dist = [];
    % get speed / distance information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp;
        y = grouped(i).dist.all(:,trial);
        pooled_temp = autoCat(pooled_temp,x,false);
        pooled_dist = autoCat(pooled_dist,y,false);
    end
    % screen out control periods (recovery holds and start/end of exp)
    tp = getTempTurnPoints(data(i).temp_protocol);
    loc = [tp.UpROI, tp.DownROI];
    loc = sort(loc);

    temp = pooled_temp(loc,:);
    dist = pooled_dist(loc,:);
    [rho,pval] = corr(temp,dist);

    plotData(exp).rho = rho(logical(eye(num.trial(i))));
    plotData(exp).pval = pval(logical(eye(num.trial(i))));
    plotData(exp).groupName = expNames(i);
    plotData(exp).color =  grouped(i).color;
end

% correlation coefficients
fig = getfig('',true,figSize);
hold on
 for i = 1:num.exp
   disp(plotData(i).groupName)
   kolor = plotData(i).color;
   xlow = i-buff-0.1;
   xhigh = i+buff+0.1;
   y = plotData(i).rho;
   y_avg = mean(plotData(i).rho);
   x = shuffle_data(linspace(i-buff,i+buff,length(y)));

   scatter(x,y,SZ,kolor,'filled')
   plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',LW)
 end
 xlim([0.5,num.exp+.5])
 ylabel('temp-distance corr. coef.')
 h_line(0,foreColor,':',1)
 formatFig(fig,blkbgd);
 set(gca,'xcolor',backColor)
 if ~autoLim
     ylim(ylimits)
 end

% save figure
save([saveDir expGroup ' temp distance correlation ramps only'],'plotData');
save_figure(fig,[saveDir expGroup ' temp distance correlation ramps only'],'-png',true,false);
save_figure(fig,[saveDir expGroup ' temp distance correlation ramps only'],'-pdf',true,true);

%% FIGURE: surface map of distance-temp tuning curve
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
autoSave = false;
autoDist = true; % automatically determine distance axis limits
autoColor = true; % automatically determine color limits for z (distance) axis
autoTemp = true; % automatically set temperature limits
buff = 0.3;
if blkbgd
    gridalpha = 0.3;
else
    gridalpha = 0.5;
end

x_limits = [0,num.exp+1];
y_limits = [14, 23]; % Temperature limits
z_limits = [14, 28];  % Distance limits

% get strain names
strain_label = cell([1,num.exp]);
for exp = 1:num.exp
    i = expOrder(exp);
    temp = strsplit(expNames{i},' ');
    strain_label{exp} = temp{1};
end

fig = getfig('',true);
set(fig,'color',backColor);
colorLimits = [];
for exp = 1:num.exp
    i = expOrder(exp);
%     temp = strsplit(expNames{i},' ');
%     strain_label{exp} = temp{1};
    z = grouped(i).dist.distavgbytemp(:,2); % distance
    loc = isnan(z); %no-data areas (temp bins outside perview of this exp)
    y = grouped(i).dist.distavgbytemp(:,1); % temperature
    % remove nans
    z(loc) = [];
    y(loc) = [];

    % Adjust temp and distance vectors for a surface plot
    Y = repmat(y,[1,2]);
    Z = repmat(z,[1,2]);

    % Build x-location map for this experiment group
    x = [exp-buff,exp+buff];
    X = repmat(x,[length(y),1]);

    surf(X,Y,Z,'edgecolor','none','facecolor','interp');

    % save the surf shape into a new folder...
    plotData = [];
    plotData.X = X; plotData.Y = Y; plotData.Z = Z;
    save([saveDir expNames{i} ' distance tuning curve surf map data.mat'],'plotData','-mat');

    % next experimental group
    colorLimits = autoCat(colorLimits,[min(z) max(z)]);
    hold on
end

% FORMATTING
colormap(flipud(parula))
set(gca,'zdir','reverse')
% xlabel('Strain')
ylabel('temperature (\circC)')
zlabel('proximity to food (mm)')
set(gca,'color',backColor,'xcolor',foreColor,'ycolor',foreColor,'zcolor',foreColor)
set(gca,'gridAlpha',gridalpha)
set(gca,'linewidth',2,'fontsize',12,'fontname','arial')

colorLimits = [min(min(colorLimits)),max(max(colorLimits))];
if autoColor
    set(gca,'CLim',colorLimits)
end
if ~autoDist
    zlim(z_limits)
end
if ~autoTemp
    ylim(y_limits)
end
xlim(x_limits)
set(gca,'XTick',1:num.exp,'xticklabel',strain_label)


% % save figure
save_figure(fig,[saveDir expGroup ' distance tuning curve surf map'],fig_type,autoSave);

%% FIGURE: Flies on food analysis
clearvars('-except',initial_vars{:})
plot_err = true;
[foreColor,backColor] = formattingColors(blkbgd);
num_lim = [0,5];
num_temp_lim = [0,1.4];
autoLim = true; %keep automatically determining the y limits
xlim_auto = true; % change the time range for the x axis
time_limits = [50,365]; % time limits if manual control over x-axis range
well_radius = 3; % 5 mm diameter of the physical well -- give 0.5mm buffer zone outside well
well_rad = well_radius * pix2mm; %convert mm to pixels

% Calculate number of flies on food for each frame
plotData = [];
for i = 1:num.exp
    N = [];
    for trial = 1:num.trial(i)
        onFood = [];
        well_loc = data(i).T.foodLoc(trial);
        well_center = data(i).data(trial).data.wellcenters(:,well_loc);
        x = data(i).data(trial).data.x_loc;
        y = data(i).data(trial).data.y_loc;
        for frame = 1:size(x,1)
            distances = sqrt((x(frame,:)-well_center(1)).^2+((y(frame,:)-well_center(2)).^2));
            onFood(frame,1) = sum(distances<well_rad);
        end
        N = autoCat(N,onFood,false);
    end
    plotData(i).N = N;
    disp(['Done exp ' num2str(i)])
end

% Cluster the flies on food by temperature?
for i = 1:num.exp
    temps = unique(data(i).G(1).TR.temps);
    rateIdx = data(i).G(1).TR.rateIdx;
    tempIdx = data(i).G(1).TR.tempIdx;
    % find rate index
    heatRate = find(data(i).G(1).TR.rates>0);
    coolRate = find(data(i).G(1).TR.rates<0);
    try
        holdRate = find(data(i).G(1).TR.rates==0);
        ntypes = 3;
    catch
        ntypes = 2;
    end

    for temp = 1:length(temps)
        for type = 1:ntypes
            switch type
                case 1 %heating
                    g_name = 'increasing';
                    idxSex = heatRate;
                case 2 %cooling
                    g_name = 'decreasing';
                    idxSex = coolRate;
                case 3
                    g_name = 'holding';
                    idxSex = holdRate;
            end
            % increasing rates:
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            plotData(i).(g_name)(temp,1) = mean(mean(plotData(i).N(loc,:),1,'omitnan'),'omitnan'); %avg
            plotData(i).(g_name)(temp,2) = std(mean(plotData(i).N(loc,:),1,'omitnan'),'omitnan');%./num.trial(i); %err
        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        plotData(i).temp_all(temp,1) = mean(mean(plotData(i).N(loc,:),1,'omitnan'),'omitnan'); %avg
        plotData(i).temp_all(temp,2) = std(mean(plotData(i).N(loc,:),1,'omitnan'),'omitnan')./num.trial(i);% %err
    end
    plotData(i).temps = temps;
end
disp('All finished')

% % Correlation between flies on food and distance:
% for i = 1:num.exp
%     N = plotData(i).N;
%     dist = grouped(i).dist.all;
%     % filter for only ramps?
%     tp = getTempTurnPoints(data(i).T.TempProtocol{1});
%     loc = [tp.UpROI,tp.DownROI];
%     loc = sort(loc');
%     N(~loc,:) = [];
%     dist(~loc,:) = [];
%     loc2 = any(isnan(N),2) | any(isnan(dist),2);
%     N(loc2,:) = [];
%     dist(loc2,:) = [];
%
%
%     for trial = 1:num.trial(i)
%         rho(trial) = corr(N(:,trial),dist(:,trial));
%     end
% end

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = 3:c:r*c; %binned distance alignment

LW = 0.75;
sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true);
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor)

    %number of flies on food over time
    subplot(r,c,sb(2).idx); hold on
        y_err = std(plotData(i).N,0,2,'omitnan');
        y_err = smooth(y_err,sSpan,'moving');
        y = mean(plotData(i).N,2,'omitnan');
        y = smooth(y,sSpan,'moving');
        % plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'color',kolor,'linewidth',LW)
        if ~autoLim
            ylim(num_lim)
        end

    %temp dependent distance
    subplot(r,c,sb(3).idx); hold on
        x = plotData(i).temps;
        y = plotData(i).temp_all(:,1);
        y_err = plotData(i).temp_all(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];

        plot(x,y,'color',kolor,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",backColor)

% flies on food
subplot(r,c,sb(2).idx)
ylabel('# flies on food')
xlabel('time (min)')
set(gca,"XColor",foreColor)
set(gca,'xlim',[50,365])
% temp-flies on food relationship
subplot(r,c,sb(3).idx)
ylabel('# flies on food')
xlabel('temp (\circC)')
if ~autoLim
    ylim(num_temp_lim)
end
% change x limits
if ~xlim_auto
    subplot(r,c,sb(1).idx)
    set(gca,'xlim',time_limits)
    subplot(r,c,sb(2).idx)
    set(gca,'xlim',time_limits)
end
% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

save_figure(fig,[saveDir expGroup ' flies on food'],fig_type);


% FIGURE: Flies on food by heating / cooling
plot_err = true;
equalLim = true;
LW = 1.5;
r = 1;
c = 2;
nMax = num.exp;

% FIGURE:
fig = getfig('',true);
% AVG
subplot(r,c,1)
hold on
for i = 1:nMax
    kolor = grouped(i).color;
    x = plotData(i).temps;
    y = plotData(i).temp_all(:,1);
    y_err = plotData(i).temp_all(:,2);
    loc = isnan(y)|isnan(y_err);
    x(loc) = [];
    y(loc) = [];
    y_err(loc) = [];

    plot(x,y,'color',kolor,'linewidth',LW+1)
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
end
xlabel('Temperature (\circC)')
ylabel('flies on food (#)')
% SEP HEAT / COOL
subplot(r,c,2)
hold on
for i = 1:nMax
    kolor = grouped(i).color;
    for type = 1:2 %increasing | decreasing
        switch type
            case 1
                section_type = 'increasing';
                l_style = '-';
            case 2
                section_type = 'decreasing';
                l_style = '--';
        end
        x = plotData(i).temps;
        y = plotData(i).(section_type)(:,1);
        y_err = plotData(i).(section_type)(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        % plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',l_style)
    end
end
xlabel('Temperature (\circC)')
ylabel('flies on food (#)')

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
if equalLim
    fig = matchAxis(fig,true);
end
% ylim(num_temp_lim)

% save figure
save_figure(fig,[saveDir 'Flies on food heating and cooling'],fig_type);

%% FIGURES: eccentricity across time
clearvars('-except',initial_vars{:})

% Pull data together
for i = 1:num.exp
    ecent = [];
    for trial = 1:num.trial(i)
        y = data(i).data(trial).data.occupancy.eccentricity(:,1);
        ecent = autoCat(ecent,y,false);
    end
    grouped(i).ecent.all = ecent;
end

% Cluster the sleeping flies by temperature
for i = 1:num.exp
    temps = unique(data(i).G(1).TR.temps);
    rateIdx = data(i).G(1).TR.rateIdx;
    tempIdx = data(i).G(1).TR.tempIdx;
    % find rate index
    heatRate = find(data(i).G(1).TR.rates>0);
    coolRate = find(data(i).G(1).TR.rates<0);
    try
        holdRate = find(data(i).G(1).TR.rates==0);
        ntypes = 3;
    catch
        ntypes = 2;
    end
    for temp = 1:length(temps)
        for type = 1:ntypes
            switch type
                case 1 %heating
                    g_name = 'increasing';
                    idxSex = heatRate;
                case 2 %cooling
                    g_name = 'decreasing';
                    idxSex = coolRate;
                case 3 %holding
                    g_name = 'holding';
                    idxSex = holdRate;
            end
            % Divided by heating / cooling
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            grouped(i).ecent.(g_name)(temp,1) = mean(mean(grouped(i).ecent.all(loc,:),2,'omitnan'),'omitnan'); %avg
            grouped(i).ecent.(g_name)(temp,2) = std(mean(grouped(i).ecent.all(loc,:),1,'omitnan'),'omitnan');%./num.trial(i); %err
        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        grouped(i).ecent.temp_all(temp,1) = mean(mean(grouped(i).ecent.all(loc,:),2,'omitnan'),'omitnan'); %avg
        grouped(i).ecent.temp_all(temp,2) = std(mean(grouped(i).ecent.all(loc,:),1,'omitnan'),'omitnan')./num.trial(i);% %err
    end
    grouped(i).ecent.temps = temps;
end

plot_err = false;
equalLim = true;
LW = 1.5;
r = 1;
c = 2;
nMax = num.exp;

% FIGURE:
fig = getfig('',true);
% AVG
subplot(r,c,1)
hold on
for i = 1:nMax
    kolor = grouped(i).color;
    x = grouped(i).ecent.temps;
    y = grouped(i).ecent.temp_all(:,1);
    y_err = grouped(i).ecent.temp_all(:,2);
    loc = isnan(y)|isnan(y_err);
    x(loc) = [];
    y(loc) = [];
    y_err(loc) = [];

    plot(x,y,'color',kolor,'linewidth',LW+1)
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
end
xlabel('Temperature (\circC)')
ylabel('eccentricity (mm)')
% SEP HEAT / COOL
subplot(r,c,2)
hold on
for i = [4,7,8] %1:nMax
    kolor = grouped(i).color;
    for type = 1:2 %increasing | decreasing
        switch type
            case 1
                section_type = 'increasing';
                l_style = '-';
            case 2
                section_type = 'decreasing';
                l_style = '--';
        end
        x = grouped(i).ecent.temps;
        y = grouped(i).ecent.(section_type)(:,1);
        y_err = grouped(i).ecent.(section_type)(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        if plot_err && ~strcmpi(fig_type,'-pdf')
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        elseif plot_err && strcmpi(fig_type,'-pdf')
                plot(x,y-y_err,'color',kolor, 'linewidth', 0.5)
                plot(x,y+y_err,'color',kolor, 'linewidth', 0.5)
        end
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',l_style)
    end
end
xlabel('Temperature (\circC)')
ylabel('eccentricity (mm)')

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
if equalLim
    fig = matchAxis(fig,true);
else
    subplot(r,c,1)
    ylim([16 30])
    subplot(r,c,2)
    ylim([16 30])
end
%

% ylim(num_temp_lim)


% save figure
save_figure(fig,[saveDir  'Eccentricity by temp'],fig_type);

%% FIGURE: 3D temperature modulation of behavior
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);
SZ = 75;

%for each trial -- plot it in 3D space based on:
% X) distance range
% Y) minimum distance to food (either during heating/cooling)
% Z) hysteresis

fig = getfig('',true);
for i = 1:num.exp
    kolor = grouped(i).color;
    x = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
    y = min([grouped(i).increasing.minDist.dist';grouped(i).decreasing.minDist.dist']);
    y(isnan(y)) = [];
    z = sum(grouped(i).decreasing.all-grouped(i).increasing.all,'omitnan');

    scatter3(x,y,z,SZ,kolor,'filled')
    hold on
end
% Labels and formatting
formatFig(fig,blkbgd);
xlabel('range (mm)')
ylabel('min distance (mm)')
zlabel('hysteresis (mm)')
set(gca,'zcolor',foreColor)

save_figure(fig,[saveDir expGroup ' 3D space modulation by temperature'],fig_type);

%% FIGURE: avg distance between wells using the pix2mm conversion
[foreColor,backColor] = formattingColors(blkbgd);

d = [];
for i = 1:num.exp
  for trial = 1:num.trial(i)
    % distance between well 1 and 3
    well_1 = data(i).data(trial).data.wellcenters(:,1);
    well_3 = data(i).data(trial).data.wellcenters(:,3);
    d = [d;sum((well_1-well_3).^2).^0.5];
    % distance between well 2 and 4
    well_1 = data(i).data(trial).data.wellcenters(:,2);
    well_3 = data(i).data(trial).data.wellcenters(:,4);
    d = [d;sum((well_1-well_3).^2).^0.5];
  end
end

well_distance = d./pix2mm;
% mean(well_distance)

fig = getfig('',true,[390 573]);
violin(well_distance,'xlabel',{''},'facecolor',...
       Color('teal'),'edgecolor',foreColor,'alp',1,...
       'bw',0.35,'mc',foreColor,'medc','r--')
ylabel('well distance (mm)','FontSize',14)
formatFig(fig,blkbgd);
h_line(36.2,'gold','--',1.5)
ylim([35,40])
set(gca,'xcolor',backColor,'ytick',35:1:40)
save_figure(fig,[saveDir expGroup ' physcial well distance selection'],fig_type);

%% FIGURE & STATS: Hysteresis within a specific time range
clearvars('-except',initial_vars{:})
ROI = [16,18]; % hysteresis range %TODO: add autochoices here
plot_err = false;
plotSig = true; %plot significance stars
[foreColor,backColor] = formattingColors(blkbgd);
LW = 1;
buff = 0.2;
SZ = 50;
alpha = 0.05;
r = 1; %rows
c = 3; %columns

% FIGURE:
fig = getfig('',true);
% Hystersis
subplot(r,c,1)
hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    %increasing
    x = grouped(i).increasing.temps;
    y = grouped(i).increasing.avg;
    y_err = grouped(i).increasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans
    y(loc) = []; x(loc) = []; y_err(loc) = [];

    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'LineWidth',LW+0.5,'Color',kolor,'linestyle','-')
    %decreasing
    x = grouped(i).decreasing.temps;
    y = grouped(i).decreasing.avg;
    y_err = grouped(i).decreasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans
    y(loc) = []; x(loc) = []; y_err(loc) = [];

    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'LineWidth',LW+.5,'Color',kolor,'linestyle','--','HandleVisibility','off');

    % Names and Colors of included data
    dataString{i} = grouped(i).name;
end
subplot(r,c,1)
set(gca,'ydir','reverse')
ylabel('proximity to food (mm)')
xlabel('temp (\circC)')
xlim([ROI(1),ROI(2)])

% Pull difference in distance heating-cooling
subplot(r,c,2)
hold on
for i = 1:num.exp
    x = repmat(grouped(i).decreasing.temps,[1,num.trial(i)]);
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    kolor = grouped(i).color;
%     plot(x,y,'color',kolor,'LineWidth',LW);
    plot(mean(x,2),mean(y,2),'color',kolor,'LineWidth',2)
end
h_line(0,foreColor,':',1)
xlim([ROI(1),ROI(2)])
xlabel('temp (\circC)')
ylabel('distance difference (mm)')

% Cumulative difference in proximity
subplot(r,c,3)
hold on
%TODO: update for smaller range or range that doesn't match exactly the temperature bins
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    %find temp ROI region:
    loc = [];
    loc(1) = find(grouped(i).decreasing.temps==ROI(1));
    loc(2) = find(grouped(i).decreasing.temps==ROI(2));
    y1 = grouped(i).decreasing.all(loc(1):loc(2),:);
    y2 = grouped(i).increasing.all(loc(1):loc(2),:);
    y = y1-y2;
    plotY = sum(y,1,'omitnan');
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
    scatter(x,plotY,SZ,kolor,"filled","o")
    plot([ii-buff,ii+buff],[mean(plotY),mean(plotY)],'color',foreColor,'LineWidth',2)
end
xlim([0.5,num.exp+0.5])
h_line(0,foreColor,':',1)
ylabel('cumulatice difference (mm)')

formatFig(fig,true,[r,c]);
set(gca,'XTick',[],'xcolor',backColor)
xlabel('Group','color',foreColor)

% STATS: are the means of any groups different from zero?
[p, mlt, id] = deal([]);
for ii = 1:num.exp
    i = expOrder(ii);
    loc = [];
    loc(1) = find(grouped(i).decreasing.temps==ROI(1));
    loc(2) = find(grouped(i).decreasing.temps==ROI(2));
    y1 = grouped(i).decreasing.all(loc(1):loc(2),:);
    y2 = grouped(i).increasing.all(loc(1):loc(2),:);
    y = y1-y2;
    plotY = sum(y,1,'omitnan');
    [~,p(ii)] = ttest(plotY);
    group_name{ii} = expNames{i};
    %multicompare
    mlt = autoCat(mlt, plotY',false);
    id = autoCat(id,i*ones(length(plotY),1),false);
end
%Bonferonni correction:
m = num.exp;
p_limit = alpha/m;
h = p<=p_limit;
stats_tbl = table(group_name',h',p','VariableNames',{'group','significant','p value'});
disp(stats_tbl)

% add significance stars to the figure:
if plotSig
    y_pos = rangeLine(fig,1);
    subplot(r,c,3); hold on
    for ii = 1:num.exp
        if h(ii)
            scatter(ii,y_pos,100,foreColor,'*')
        end
    end
end


% STATS:
% determine which groups differ from each other
[~,~,stats] = anova1(mlt(:),id(:),'off');
[c,~,~,~] = multcompare(stats,alpha,'off');
% bonferonni multiple comparisons correction
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nPosition hysteresis cross group comparison statistics\n\n')
[Group1,Group2,P_Value] = deal([]);
idx = 0;
for i = 1:length(significantHypotheses)
    if significantHypotheses(i)
        idx = idx+1;
        Group1{idx,1} = expNames{c(i,1)};
        Group2{idx,1} = expNames{c(i,2)};
        P_Value(idx,1) = c(i,6);
    end
end
sig_comp = table(Group1,Group2,P_Value);
disp(sig_comp)

% save figure
save_figure(fig,[saveDir 'Hysteresis summary ' num2str(ROI(1)) ' to ' num2str(ROI(2)) ' deg'],fig_type);

%% FIGURE: [temp shift experiments] align distance by events not temp
% TODO: do this for combined and separated heating/cooling distances
clearvars('-except',initial_vars{:})
plot_err = true;
equalLim = true;
LW = 1.5;
r = 1;
c = 2;
autoLim = false;
ylimits = [10 28];
nMax = num.exp;
xmax = [];

% FIGURE:
fig = getfig('Temp Shift Event Aligned',true); %set(fig,'color','w',"Position",[1934 468 1061 590])
% AVG
subplot(r,c,1)
hold on
for i = 1:nMax
    kolor = grouped(i).color;
    y = grouped(i).dist.distavgbytemp(:,2);
    x = 1:length(y);

    y_err = grouped(i).dist.distavgbytemp_err(:,2);
    loc = isnan(y)|isnan(y_err);
    x(loc) = [];
    y(loc) = [];
    y_err(loc) = [];

    plot(x,y,'color',kolor,'linewidth',LW+1)
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
    dataString{i} = grouped(i).name;
    xmax = autoCat(xmax,x(end));
end
xmax = max(xmax);
xlim([0,xmax+1])
set(gca,'xtick',[1,xmax],'xticklabels',{'Min','Max'},'ydir', 'reverse')
xlabel('Temperature')
ylabel('proxiity to food (mm)')
% SEP HEAT / COOL
subplot(r,c,2)
hold on
xmax = [];
for i = 1:nMax
    kolor = grouped(i).color;
    %increasing
    y = grouped(i).increasing.avg;
    x = 1:length(y);
    y_err = grouped(i).increasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans
    y(loc) = []; x(loc) = []; y_err(loc) = [];
    % plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
    plot(x,y,'LineWidth',LW+0.5,'Color',kolor,'linestyle','-')
    %decreasing
    y = grouped(i).decreasing.avg;
    x = 1:length(y);
    y_err = grouped(i).decreasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans
    y(loc) = []; x(loc) = []; y_err(loc) = [];

    % plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
    plot(x,y,'LineWidth',LW+.5,'Color',kolor,'linestyle','--','HandleVisibility','off');

    % Names and Colors of included data
    dataString{i} = grouped(i).name;
    xmax = autoCat(xmax,x(end));
end
xmax = max(xmax);
xlim([0,xmax+1])
set(gca,'xtick',[1,xmax],'xticklabels',{'Min','Max'},'ydir', 'reverse')
xlabel('Temperature')
ylabel('proxiity to food (mm)')
% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
if equalLim
    fig = matchAxis(fig,true);
end
if ~autoLim
    subplot(r,c,1)
    ylim(ylimits)
    subplot(r,c,2)
    ylim(ylimits)
end

% save figure
save_figure(fig,[saveDir 'Min max temp aligned proximity to food'],fig_type);

%% FIGURE: [WT comp] latutide organized
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);
corr_coef = [];
buff = 0.2;
SZ = 50;
LW = 1.5;

lat_list = {'Swedish', 'Berlin', 'Oregon','Canton','Malawi', 'Zimbabwe'};
latitudes = [60.1282,  52.5200,    43.8041,   40.7989, -13.2543, -19.0154];

% get correlation data
for i = 1:num.exp
    pooledData = [];
    % get speed / distance information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp;
        y = grouped(i).dist.all(:,trial);
        temp = [];
        temp = autoCat(temp,x,false);%temp for each trial within the exp.
        temp = autoCat(temp,y,false);%dist for each trial within the exp
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-distance correlation
        rho = corr(temp);
        corr_coef(i).all(trial) = rho(1,2);
        % save data for pooled comparison
        pooledData = [pooledData; temp];
    end
    % Pooled speed-distance correlation
    rho = corr(pooledData);
    corr_coef(i).group = rho(1,2);
end


% correlation coefficients
fig = getfig('',true,[453 590]); hold on
hold on
 for ii = 1:num.exp
   i = find(contains(expNames,lat_list(ii)));
   kolor = grouped(i).color;
   lat = abs(latitudes(ii));
   xlow = lat-2;
   xhigh = lat+2;
%    x = shuffle_data(linspace(lat-buff,lat+buff,num.trial(i)));
   x = lat*ones(1,num.trial(i));
   y = corr_coef(i).all;
   y_avg = mean(corr_coef(i).all);
   scatter(x,y,SZ,kolor,'filled')
%    plot([xlow,xhigh],[corr_coef(i).group,corr_coef(i).group],'color',kolor,'linestyle',':','linewidth',LW)
%    plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',1.5)
 end
 xlim([-5,65])
 ylabel('temperature-distance correlation')
 h_line(0,foreColor,':',1)
 formatFig(fig,blkbgd);
 set(gca,'xcolor',foreColor)
 xlabel('latitude (\circ)')
% save figure
save_figure(fig,[saveDir expGroup ' temp distance correlation by latitude 2'],fig_type);

%% FIGURE: [temperature rate experiments] surf plot tuning curve
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
blkbgd = true;
autoSave = false; % auto save figure?
autoDist = true; % automatically determine distance axis limits
autoColor = true; % automatically determine color limits for z (distance) axis
autoTemp = true; % automatically set temperature limits
buff = 0.3;
LW = 1.5;

if blkbgd
    gridalpha = 0.3;
else
    gridalpha = 0.5;
end

x_limits = [0,num.exp+1];
y_limits = [14, 23]; % Temperature limits
z_limits = [14, 28];  % Distance limits

% get data to form 3-d mesh grid...
strain_label = cell([1,num.exp]);
[X,Y,Z,plotData] = deal([]);
for exp = 1:num.exp
    i = expOrder(exp);
    z = grouped(i).dist.distavgbytemp(:,2); % distance
    loc = isnan(z); %no-data areas (temp bins outside perview of this exp)
    y = grouped(i).dist.distavgbytemp(:,1); % temperature
    % remove nans
    z(loc) = [];
    y(loc) = [];
    % get temp-rate data
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    tempRates = tp.rates;
    if any(tempRates(tempRates==0))
        tempRates(tempRates==0) = [];
    end
    if ~any(diff(abs(tempRates))==0)
        warndlg('Cannot run this with different rate protocols')
        return
    end
    TR = abs(tempRates(1));
    x = TR*ones(size(y));
    % Build plotting surface mesh
    X = autoCat(X,x,false);
    Y = autoCat(Y,y,false);
    Z = autoCat(Z,z,false);
    plotData = autoCat(plotData,[x,y,z]);
end

% Plot figure:
fig = getfig('',true);
set(fig,'color',backColor);
colorLimits = [];
surf(X,Y,Z,'edgecolor','none','facecolor','interp');

% FORMATTING
colormap(flipud(parula))
set(gca,'zdir','reverse')
xlabel('Temp rate (\circC/min)')
ylabel('temperature (\circC)')
zlabel('proximity to food (mm)')
set(gca,'color',backColor,'xcolor',foreColor,'ycolor',foreColor,'zcolor',foreColor)
set(gca,'gridAlpha',gridalpha)
set(gca,'linewidth',2,'fontsize',18,'fontname','arial')

% % % save figure
% save_figure(fig,[saveDir expGroup ' temp rate distance tuning curve'],fig_type,true);


% add actual plot lines
hold on
for i = 1:num.exp
    plot3(X(:,i),Y(:,i),Z(:,i),'color',foreColor,'linewidth',LW)
end

c = colorbar;
c.FontSize = 15;
c.Color = foreColor;
c.Label.String = 'Proximity to food (mm)';
set(gca,'xgrid','off','ygrid','off','zgrid','off')



%%

% Plot the number of flies in the arena over time... --> how many dropped
% flies are there?? (dying flies)
lostFlies  = [];
for i = 1:num.exp
MT = [];
    for trial = 1:num.trial(i)
        MT = autoCat(MT, data(i).data(trial).data.T.flyCount,false);
    end
    lostFlies(i).data = MT;
end


sSpan = 10;
fig = getfig('',1); hold on
for i = 1:num.exp
    x = grouped(i).time;
    plot(x, smooth(mean(lostFlies(i).data, 2, 'omitnan'),sSpan, 'moving'),'color', grouped(i).color)
end
xlabel('time (min)')
ylabel('# flies')

formatFig(fig,true);



%% FIGURE: SINGLE TRIAL LINES over-lap of time-trials and temperature protocols NO SPEED
clearvars('-except',initial_vars{:})
plot_err = true;
autoLim = true;
% Y limit ranges
dist_lim = [10,35];       %distance
dt_lim = [14, 32];        %distance-temp
xlim_auto = true; % change the time range for the x axis
time_limits = [50,365]; % time limits if manual control over x-axis range
nMax =  num.exp;%
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = 3:c:r*c; %binned distance alignment

LW = 0.75;
sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true);
for i = num.exp:-1:1
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor)

    %distance
    subplot(r,c,sb(2).idx); hold on
        for trial = 1:num.trial(i)
            y = smooth(grouped(i).dist.all(:,trial),sSpan, 'moving');
            plot(x,y,'LineWidth',LW,'Color',kolor)
            if ~autoLim
                ylim(dist_lim)
            end
        end

    %temp dependent distance
    subplot(r,c,sb(3).idx); hold on
        for trial = 1:num.trial(i)
            x = grouped(i).dist.tempList(trial,:);
            y = grouped(i).dist.tempBinned(trial,:);
            loc = isnan(y);
            x(loc) = [];  y(loc) = [];
            plot(x,y,'color',kolor,'linewidth',1.25)
        end
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",backColor)
% distance
subplot(r,c,sb(2).idx)
ylabel('proximity to food (mm)')
xlabel('time (min)')
set(gca,'ydir','reverse')
% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel('proximity to food (mm)')
xlabel('temp (\circC)')
if ~autoLim
    ylim(dt_lim)
end
h_line(18.1,'grey',':',1) %36.2
set(gca,'ydir','reverse')

if ~xlim_auto
    subplot(r,c,sb(1).idx)
    set(gca, 'xlim', time_limits)
    subplot(r,c,sb(2).idx)
    set(gca, 'xlim', time_limits)
end



% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' timecourse summary all trial lines'],fig_type);
