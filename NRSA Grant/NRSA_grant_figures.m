
%% Quick hypothesis testing:
clear

% folder = getCloudPath;
figDir = [getCloudPath, 'Electrophysiology Modeling\'];

% raw start data
T = 15:5:35;
PN_base = [0.1, 0.5, 0.7, 0.8, 0.85];
TRN_base = [1, 0.25, 0, 0.25, 1];

% upsample the data
N = 100;
temp = linspace(15,35,N);
PN = interp1(T,PN_base, temp,'spline');
TRN = interp1(T,TRN_base, temp,'spline');

% determine the effect of subtractive inhibition with rectification
PN_TRN = PN - TRN;
PN_TRN(PN_TRN<0) = 0; % rectify

% Plot base figure:
LW = 2;
fig = figure; hold on
    plot(temp, PN, 'color', 'k','linewidth', LW) % upsampled PN activity
    plot(temp, TRN, 'color', 'r','linewidth', LW) % upsampled  TRN activity
    % PN with TRN activity
    plot(temp, PN_TRN, 'color',Color('dodgerblue'),'linewidth', LW) % upsampled
% Labels
xlabel('Temperature (\circC)')
ylabel('Relative activity (a.u.)')
formatFig(fig, false);

save_figure(fig,[figDir 'Strong TRN inhibition'],'-png',false,false);


% Plot all lines figure:
LW = 0.5;
n = 30;
gains = linspace(0.1,1,n); % test a range of TRN activity levels
CList = Color('Cyan', 'Darkblue', n);
ZList = Color('yellow', 'darkred',n);
Z = [];
fig = figure; hold on
    % TRN activity
    for i = 1:n
        Z(:,i) = TRN*gains(i);
        plot(temp, Z(:,i), 'color', ZList(i,:),'linewidth', LW) % upsampled PN activity
    end
    % new PN activity
    for i = 1:n
        Y = PN' - Z(:,i);
        Y(Y<0) = 0;% rectify
        plot(temp, Y, 'color',CList(i,:),'linewidth', LW) % upsampled
    end
    plot(temp, PN, 'color', 'k','linewidth', 2) % upsampled PN activity

% Labels
xlabel('Temperature (\circC)')
ylabel('Relative activity (a.u.)')
formatFig(fig, false);
save_figure(fig,[figDir 'Range of TRN inhibition'],'-png',false,false);

%% Quick hypothesis testing: Linear PN activity change
clear

% folder = getCloudPath;
figDir = [getCloudPath, 'Electrophysiology Modeling\'];

% raw start data
T = 15:5:35;
PN_base = [0.1, 0.5, 0.7, 0.8, 0.85];
TRN_base = [1, 0.25, 0, 0.25, 1];

% upsample the data
N = 100;
temp = linspace(15,35,N);
PN = interp1(T,PN_base, temp,'spline');
TRN = interp1(T,TRN_base, temp,'spline');

% determine the effect of subtractive inhibition with rectification
PN_TRN = PN - TRN;
PN_TRN(PN_TRN<0) = 0; % rectify

% Plot base figure:
LW = 2;
fig = figure; hold on
    plot(temp, PN, 'color', 'k','linewidth', LW) % upsampled PN activity
    plot(temp, TRN, 'color', 'r','linewidth', LW) % upsampled  TRN activity
    % PN with TRN activity
    plot(temp, PN_TRN, 'color',Color('dodgerblue'),'linewidth', LW) % upsampled
% Labels
xlabel('Temperature (\circC)')
ylabel('Relative activity (a.u.)')
formatFig(fig, false);

save_figure(fig,[figDir 'Strong TRN inhibition'],'-png',false,false);


% Plot all lines figure:
LW = 0.5;
n = 30;
gains = linspace(0.1,1,n); % test a range of TRN activity levels
CList = Color('Cyan', 'Darkblue', n);
ZList = Color('yellow', 'darkred',n);
Z = [];
fig = figure; hold on
    % TRN activity
    for i = 1:n
        Z(:,i) = TRN*gains(i);
        plot(temp, Z(:,i), 'color', ZList(i,:),'linewidth', LW) % upsampled PN activity
    end
    % new PN activity
    for i = 1:n
        Y = PN' - Z(:,i);
        Y(Y<0) = 0;% rectify
        plot(temp, Y, 'color',CList(i,:),'linewidth', LW) % upsampled
    end
    plot(temp, PN, 'color', 'k','linewidth', 2) % upsampled PN activity

% Labels
xlabel('Temperature (\circC)')
ylabel('Relative activity (a.u.)')
formatFig(fig, false);
save_figure(fig,[figDir 'Range of TRN inhibition'],'-png',false,false);



%% FIGURE: heating and cooling separated vertical temp colored OCCUPATION PROBABILITY
% load data from QuadStep4.2 first (specifically, data with waxed antenna / MP vs
% control caviar data ('Berlin F LRR 25-17 sensory components')
clearvars('-except',initial_vars{:})
% blkbgd = true;  fig_type = '-png'; 
fig_type = '-pdf'; blkbgd = false;
buff = 0.1;
sz = 50;
autoLim = true;
y_lim = [0,1];
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

dataString = cell([1,num.exp]);

% Find the max temp and min temp of all the experiments 
temp_min = 16; 
temp_max = 26;
temp_bin = 0.5;
buff = 0.25;
sz = 30;
LW = 3;

expList = [2, 2, 4, 4];
tempList = [17, 25, 17, 25];
cList = {'dodgerblue', 'red', 'dodgerblue', 'red'};

% FIGURE:
fig = getfig('',true,[565 649]); hold on
testData = []; % save data here for statistical comparisons
for i = 1:length(expList)
    kolor = Color(cList{i});
    exp = expList(i); 
    x = shuffle_data(i + linspace(-buff,buff,num.trial(exp)));
    tempLoc = find(grouped(exp).occ.temps==tempList(i));
    y = mean([grouped(exp).occ.increasing.raw(tempLoc,:);...
           grouped(exp).occ.decreasing.raw(tempLoc,:)]).*100; % mean per trial of warming and cooling
    % y = grouped(exp).occ.decreasing.raw(tempLoc,:).*100; % mean per trial of warming and cooling
    testData = autoCat(testData, y', false);
    y_avg = mean(y, 'omitnan');
    scatter(x,y,sz,kolor,'filled')
    plot([min(x)-0.1,max(x)+0.1],[y_avg,y_avg],'Color',kolor,'linewidth', LW)
end

% FORMATING AND LABELS
formatFig(fig,blkbgd);
set(gca, 'TickDir', 'out') 
ylabel('flies in food quadrant (%)')

save_figure(fig,[saveDir 'NRSA grant occupancy comparison waxed vs intact'],fig_type);


% STATS
datastats.all = testData;
datastats.id = {'intact cold', 'intact warm', 'blocked cold', 'blocked warm'};

% determine which groups differ from each other
[~,~,stats] = anova1(datastats.all,datastats.id,'off');
alpha = 0.05; %significance level
[c,~,~,~] = multcompare(stats,alpha,'off');

% bonferonni multiple comparisons correction
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nAverage occupancy statistics\n\n')
[Group1,Group2, P_Value] = deal([]);
idx = 0;
if ~any(significantHypotheses)
    disp('No statistical differences in avg speed between groups')
else
    for i = 1:length(significantHypotheses)
        if significantHypotheses(i)
            idx = idx+1;
            Group1{idx,1} = datastats.id{c(i,1)};
            Group2{idx,1} = datastats.id{c(i,2)};
            P_Value(idx,1) = c(i,6);
        end
    end
    sig_comp = table(Group1,Group2,P_Value);
    disp(sig_comp)
end




%% 
%% FIGURE: scattered OCCUPATION PROBABILITY for selected temperatures
% load data from QuadStep4.2 first (specifically, data with waxed antenna / MP vs
% control caviar data ('Berlin F LRR 25-17 sensory components')
clearvars('-except',initial_vars{:})
% blkbgd = true;  fig_type = '-png'; 
fig_type = '-pdf'; blkbgd = false;
buff = 0.1;
sz = 50;
autoLim = true;
y_lim = [0,1];
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

dataString = cell([1,num.exp]);

% Find the max temp and min temp of all the experiments 
temp_min = 16; 
temp_max = 26;
temp_bin = 0.5;
buff = 0.25;
sz = 30;
LW = 3;

expList = [2,  2,  4, 4];
tempList = [17,  25, 17, 25];
cList = {'dodgerblue',  'red', 'dodgerblue',  'red'};

% FIGURE:
fig = getfig('',true,[565 649]); hold on
testData = []; % save data here for statistical comparisons
for i = 1:length(expList)
    kolor = Color(cList{i});
    exp = expList(i); 
    x = shuffle_data(i + linspace(-buff,buff,num.trial(exp)));
    tempLoc = find(grouped(exp).occ.temps==tempList(i));
    y = mean([grouped(exp).occ.increasing.raw(tempLoc,:);...
           grouped(exp).occ.decreasing.raw(tempLoc,:)]).*100; % mean per trial of warming and cooling
    % y = grouped(exp).occ.decreasing.raw(tempLoc,:).*100; % mean per trial of warming and cooling
    testData = autoCat(testData, y', false);
    y_avg = mean(y, 'omitnan');
    scatter(x,y,sz,kolor,'filled')
    plot([min(x)-0.1,max(x)+0.1],[y_avg,y_avg],'Color',kolor,'linewidth', LW)
end

% FORMATING AND LABELS
formatFig(fig,blkbgd);
set(gca, 'TickDir', 'out') 
ylabel('flies in food quadrant (%)')
ylim([0,100])

save_figure(fig,[saveDir 'NRSA grant occupancy comparison waxed vs intact 3-temp'],fig_type);

% STATS   
datastats.all = testData;
datastats.id = {'intact 17',  'intact 25', 'blocked 17', 'blocked 25'};

% determine which groups differ from each other
[~,~,stats] = anova1(datastats.all,datastats.id,'off');
alpha = 0.05; %significance level
[c,~,~,~] = multcompare(stats,alpha,'off');

% bonferonni multiple comparisons correction
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nAverage occupancy statistics\n\n')
[Group1,Group2, P_Value] = deal([]);
idx = 0;
if ~any(significantHypotheses)
    disp('No statistical differences in avg speed between groups')
else
    for i = 1:length(significantHypotheses)
        if significantHypotheses(i)
            idx = idx+1;
            Group1{idx,1} = datastats.id{c(i,1)};
            Group2{idx,1} = datastats.id{c(i,2)};
            P_Value(idx,1) = c(i,6);
        end
    end
    sig_comp = table(Group1,Group2,P_Value);
    disp(sig_comp)
end

%% FIGURE:  OCCUPATION timecourse standard figure

clearvars('-except',initial_vars{:})

combined_dt = true; 
% Fast temp protocols:
% x_limit = [17,25];
% stepSize = 2;
% % giant ramp protocols: 
x_limit = [10,35];
stepSize = 5;

plot_err = true;
autoLim = false;
% Y limit ranges
dist_lim = [0, 60];       %distance
dt_lim = [0, 60];        %distance-temp
auto_time = true;      % automated time axis limits
time_lim = [0,400];     %time limit (x-axis)
show_exp = 1:2; %[2,4];
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
for i = show_exp % 1:nMax
%     i = expOrder(ii);
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor)

    %distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).occ.avg,'moving',sSpan)*100;
%         y_err = smooth(grouped(i).dist.err,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        if ~autoLim
            ylim(dist_lim)
        end

    %temp dependent distance
    subplot(r,c,sb(3).idx); hold on
    if combined_dt
            x = grouped(i).occ.temps; 
            y1 = [grouped(i).occ.increasing.raw, grouped(i).occ.decreasing.raw];
            y =  mean(y1,2,'omitnan').*100;
            y_err = std(y1,0,2,'omitnan').*100;
            loc = isnan(y)|isnan(y_err);
            x(loc) = [];
            y(loc) = [];
            y_err(loc) = [];
            plot(x,y,'color',kolor,'linewidth',LW+1)
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
            dataString{i} = grouped(i).name;
    else
        for type = 1:2
            x = grouped(i).occ.temps; 
            switch type
                case 1
                    y1 = grouped(i).occ.increasing.raw;
                    LS = '-';
                case 2 
                    y1 = grouped(i).occ.decreasing.raw;
                    LS = '--';
            end
            y =  mean(y1,2,'omitnan').*100;
            y_err = std(y1,0,2,'omitnan').*100;
            loc = isnan(y)|isnan(y_err);
            x(loc) = [];
            y(loc) = [];
            y_err(loc) = [];
    
            plot(x,y,'color',kolor,'linewidth',LW+1,'LineStyle',LS)
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
            dataString{i} = grouped(i).name;
        end
    end
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
ylabel('flies in food region (%)')
xlabel('time (min)')
set(gca,'ydir','reverse')
if ~auto_time
    xlim(time_lim)
end
% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel('flies in food region (%)')
xlabel('temp (\circC)')
if ~autoLim
    ylim(dt_lim)
end
h_line(14.4,'grey',':',1) %36.2
xlim(x_limit)
% ylim([10,80])
set(gca, 'tickdir', 'out')
set(gca, 'xtick', x_limit(1):stepSize:x_limit(2))


% save figure
if combined_dt
    fig_title = ' H C combined';
else 
    fig_title = ' H C separated';
end

save_figure(fig,[saveDir expGroup ' occupation percent' fig_title ' timecourse summary'],fig_type);



%% FIGURE & STATS: scatter plot speed for selected temperature points
% load data from QuadStep4.2 first (specifically, data with waxed antenna / MP vs
% control caviar data ('Berlin F LRR 25-17 sensory components')
clearvars('-except',initial_vars{:})
% blkbgd = true;  fig_type = '-png'; 
fig_type = '-pdf'; blkbgd = false;
buff = 0.1;
sz = 50;
autoLim = true;
y_lim = [0,1];
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

dataString = cell([1,num.exp]);

% Find the max temp and min temp of all the experiments 
temp_min = 16; 
temp_max = 26;
temp_bin = 0.5;
buff = 0.25;
sz = 30;
LW = 3;

expList = [2, 2,  4, 4];
tempList = [17,  25, 17, 25];
cList = {'dodgerblue', 'red', 'dodgerblue', 'red'};

% FIGURE:
fig = getfig('',true,[565 649]); hold on
testData = []; % save data here for statistical comparisons
for i = 1:length(expList)
    kolor = Color(cList{i});
    exp = expList(i); 
    x = shuffle_data(i + linspace(-buff,buff,num.trial(exp)));
    tempLoc = find(grouped(exp).speed.temps==tempList(i));
    y = mean([grouped(exp).speed.increasing.raw(tempLoc,:);...
           grouped(exp).speed.decreasing.raw(tempLoc,:)]); % mean per trial of warming and cooling
    % y = grouped(exp).occ.decreasing.raw(tempLoc,:); % mean per trial of warming and cooling
    testData = autoCat(testData, y', false);
    y_avg = mean(y, 'omitnan');
    scatter(x,y,sz,kolor,'filled')
    plot([min(x)-0.1,max(x)+0.1],[y_avg,y_avg],'Color',foreColor,'linewidth', LW)
end

% FORMATING AND LABELS
formatFig(fig,blkbgd);
set(gca, 'TickDir', 'out') 
ylabel('avg fly speed (mm/s)')
% ylim([0,100])

save_figure(fig,[saveDir 'NRSA grant speed comparison waxed vs intact 3-temp'],fig_type);


% STATS   
datastats.all  = testData;
datastats.id = {'intact 17',  'intact 25', 'blocked 17',  'blocked 25'};

% determine which groups differ from each other
[~,~,stats] = anova1(datastats.all,datastats.id,'off');
alpha = 0.05; %significance level
[c,~,~,~] = multcompare(stats,alpha,'off');

% bonferonni multiple comparisons correction
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nAverage occupancy statistics\n\n')
[Group1,Group2, P_Value] = deal([]);
idx = 0;
if ~any(significantHypotheses)
    disp('No statistical differences in avg speed between groups')
else
    for i = 1:length(significantHypotheses)
        if significantHypotheses(i)
            idx = idx+1;
            Group1{idx,1} = datastats.id{c(i,1)};
            Group2{idx,1} = datastats.id{c(i,2)};
            P_Value(idx,1) = c(i,6);
        end
    end
    sig_comp = table(Group1,Group2,P_Value);
    disp(sig_comp)
end



%% FIGURE: Temperature -- occupancy correlation
% correlation ONLY during ramps

clearvars('-except',initial_vars{:})

[foreColor,backColor] = formattingColors(blkbgd);
ylimits = [-0.8,0.8];
autoLim = false;
corr_coef = [];
buff = 0.2;
SZ = 50;
LW = 1.5;

pixWidth = 60; % additional figure pixel size for image for each extra experiment group
figSize = [pixWidth + (pixWidth*num.exp),590];

% get correlation data
plotData = struct;
for exp = 1:num.exp
    i = expOrder(exp);
    disp(expNames{i})
    if data(i).hold_exp % skip hold trials since there is only one temperature point...
        [plotData(exp).rho, plotData(exp).pval]  = deal(nan);
        plotData(exp).groupName = expNames(i);
        plotData(exp).color =  grouped(i).color;
        continue
    end
    pooled_temp = [];
    pooled_dist = [];
    % get speed / occupancy information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp;
        well_loc = data(i).T.foodLoc(trial);
        y = data(i).data(trial).occupancy.occ(:,well_loc);
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
 ylabel('Pearson correlation of temperature and occupation')
 h_line(0,foreColor,':',1)
 formatFig(fig,blkbgd);
 set(gca,'xcolor',backColor)
 if ~autoLim
     ylim(ylimits)
 end

 % Perform the Mann-Whitney U test (Wilcoxon rank-sum test)
[p, h, stats] = ranksum(plotData(1).rho, plotData(3).rho); % caviar (1) vs waxed (3)
% Display the results
fprintf('p-value: %.4f\n', p);
if h == 0
    fprintf('The null hypothesis cannot be rejected: The two groups are not significantly different.\n');
else
    fprintf('The null hypothesis is rejected: The two groups are significantly different.\n');
end

% save figure
save([saveDir expGroup ' temp occupation correlation ramps only'],'plotData');
save_figure(fig,[saveDir expGroup ' temp occupation correlation ramps only'],'-png',true,false);
save_figure(fig,[saveDir expGroup ' temp occupation correlation ramps only'],'-pdf',true,true);

%% FIGURE: Temperature -- speed correlation
% correlation ONLY during ramps

clearvars('-except',initial_vars{:})

[foreColor,backColor] = formattingColors(blkbgd);
ylimits = [-0.8,0.8];
autoLim = false;
corr_coef = [];
buff = 0.2;
SZ = 50;
LW = 1.5;

pixWidth = 60; % additional figure pixel size for image for each extra experiment group
figSize = [pixWidth + (pixWidth*num.exp),590];

% get correlation data
plotData = struct;
for exp = 1:num.exp
    i = expOrder(exp);
    disp(expNames{i})
    if data(i).hold_exp % skip hold trials since there is only one temperature point...
        [plotData(exp).rho, plotData(exp).pval]  = deal(nan);
        plotData(exp).groupName = expNames(i);
        plotData(exp).color =  grouped(i).color;
        continue
    end
    pooled_temp  = [];
    % get speed  information
    for trial = 1:num.trial(i)
         x = data(i).data(trial).occupancy.temp;
         pooled_temp = autoCat(pooled_temp,x,false);
    end
    pooled_dist = grouped(i).speed.all;

    % screen out control periods (recovery holds and start/end of exp)
    tp = getTempTurnPoints(data(i).temp_protocol);
    loc = [tp.UpROI, tp.DownROI];
    loc = sort(loc);
    
    temp = pooled_temp(loc,:);
    dist = pooled_dist(loc,:);
    nanLoc = any(isnan(dist),2);
    temp(nanLoc,:) = [];
    dist(nanLoc,:) = [];
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
 ylabel('Pearson correlation of temperature and speed')
 h_line(0,foreColor,':',1)
 formatFig(fig,blkbgd);
 set(gca,'xcolor',backColor)
 if ~autoLim
     ylim(ylimits)
 end

% Perform the Mann-Whitney U test (Wilcoxon rank-sum test)
[p, h, stats] = ranksum(plotData(1).rho, plotData(3).rho); % caviar (1) vs waxed (3)
% Display the results
fprintf('p-value: %.4f\n', p);
if h == 0
    fprintf('The null hypothesis cannot be rejected: The two groups are not significantly different.\n');
else
    fprintf('The null hypothesis is rejected: The two groups are significantly different.\n');
end


% save figure
save([saveDir expGroup ' temp speed correlation ramps only'],'plotData');
save_figure(fig,[saveDir expGroup ' temp speed correlation ramps only'],'-png',true,false);
save_figure(fig,[saveDir expGroup ' temp speed correlation ramps only'],'-pdf',true,true);



%%  statistical comparison on two correlations: 


%% Determine the area contained within the 

pix2mm = 12.8; %conversion from pixels to mm for these videos
radii = 165; %well surround regions
r = 435; % radius of the arena

% what percent of the arena is taken but the circle?
full_arena_area = pi*r^2;
special_circle_area = pi*radii^2;
(special_circle_area/full_arena_area)*100


%% Antoine Equation
% https://webbook.nist.gov/chemistry/name-ser/

clear

ko = 272.15; % kelvin offset
t = 15:35;
T = t+ko;

odors = {'2-butanone', 'ethyl acetate', 'acetic acid','3-octanol','ethyl butyrate','methyl acetate'};
coeffc = [3.9894	1150.207	-63.904;...
                4.22809	1245.702	-55.189;...
                4.68206	1642.54	-39.764;...
                4.8465	1663.322	-97.47;...
                4.33187	1509.443	-45.284;...
                4.20364	1164.426	-52.69] ;

CList = Color('teal', 'orange', length(odors));
fold_change = [];
c= 3;
r = 1;
sb(1).idx = 1:2;
sb(2).idx = 3;

fig = getfig('',1);
subplot(r,c,sb(1).idx)
hold on
for i = 1:length(odors)
    A = coeffc(i, 1);
    B = coeffc(i,2);
    C = coeffc(i,3);

    P = 10.^(A-(B./(T + C))); %pressure in BAR
    fold_change(i) = (P(end)-P(1))/P(1);
    plot(t, P, 'DisplayName', odors{i},'linewidth', 1.5,'color', CList(i,:));
end
xlabel('Temperature (\circC)');
ylabel('Vapor Pressure (bar)');
hold off;
set(gca, 'tickdir', 'out')
legend show;
legend('box', 'off','location', 'northwest','FontSize',10)
 
subplot(r,c,sb(2).idx)
scatter(1:length(odors), fold_change,100,CList,'filled')
xlim([0,length(odors)+1])
ax = gca;
set(ax, 'xtick', 1:length(odors),'xticklabel', odors)
ylabel('fold change from 15\circC to 35\circC')
formatFig(fig, false,[r,c],sb);

