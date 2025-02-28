
%% 

% test if the data comes from a normal distribution:
h = kstest(x);


%% FIGURES & STATS: Thermal Threat / Temperature Escape Behavior
clearvars('-except',initial_vars{:})
fig_dir = createFolder([saveDir, 'Stats/']);
[foreColor,~] = formattingColors(blkbgd);

% temperature region of interest: 
%TODO: make this more dynamic so that it works for different temperture protocols
temps = [16.5,18]; % temperatures to bin
SD = []; % stats data structure
[SD.w_occ,SD.c_occ,SD.w_speed,SD.c_speed] = deal(nan([max(num.trial),num.exp])); % empty mat for the trial average data to fill
SD.speed = nan([max(num.trial),num.exp]); % empty mat for the trial average data to fill

for ii = 1 : num.exp
    i = expOrder(ii); % set this in experiment order
    % Find the selected temperature regime time points: 
    x = grouped(i).ring.temps;
    idx = find(x>=temps(1) & x<=temps(2)); % index for pulling the address of time points within testing temp band
    cool_frames = [];  warm_frames = []; % set empty matrex for the frame list to fill
    for t = 1:length(idx)
        coolIdx = find(grouped(i).position.temp_rates<0); % warming rate index
        warmIdx = find(grouped(i).position.temp_rates>0); % cooling rate index
        % compile the frame numbers that occur during the cooling and warming
        % periods within the above selected temperature range
        cool_frames = [cool_frames; grouped(i).position.loc(coolIdx,idx(t)).frames]; 
        warm_frames = [warm_frames; grouped(i).position.loc(warmIdx,idx(t)).frames];
    end
    
    % load the raw ring occupancy and speed data: 
    rawOcc = grouped(i).ring.all; % ring occupancy data for all trials
    rawSpeed = grouped(i).speed.all; % speed data for all trials
    
    % pull out the trial average for the heating and cooling periods that will be statistically compared and plotted: 
    SD.w_occ(1:num.trial(i),ii) = mean(rawOcc(warm_frames,:),1,'omitnan');
    SD.c_occ(1:num.trial(i),ii) = mean(rawOcc(cool_frames,:),1,'omitnan');
    SD.w_speed(1:num.trial(i),ii) = mean(rawSpeed(warm_frames,:),1,'omitnan');
    SD.c_speed(1:num.trial(i),ii) = mean(rawSpeed(cool_frames,:),1,'omitnan');
    SD.name{ii} = grouped(i).name;
end

% statistics: are there differences betwen heating and cooling in the occupancy data for the selected region?
[group_name, p, ~, ~] = deal([]);
for ii = 1:num.exp
    [~,p(ii)] = ttest(SD.c_occ(:,ii), SD.w_occ(:,ii));
    group_name{ii} = SD.name{ii};
    % %multicompare
    % mlt = autoCat(mlt, plotY',false);
    % id = autoCat(id,ii*ones(length(plotY),1),false);
end
%Bonferonni correction:
alpha = 0.05;
m = num.exp;
p_limit = alpha/m;
h = p<=p_limit;
disp('Differences between selected temp region in ring occupancy:')
stats_tbl = table(group_name',h',p','VariableNames',{'group','significant','p value'});
disp(stats_tbl)

% FIGURE:
buff = 0.1;
buff2 = 0.2;
sz = 35;
lw = 2;
fig = getfig('',true,[711 680]); hold on
for i = 1:num.exp
    % warming
    xx = i-buff2;
    kolor = Color('red');
    rawY = SD.w_occ(:,i);
    rawY(isnan(rawY)) = [];
    x = shuffle_data(linspace(xx-buff, xx+buff,length(rawY)));
    scatter(x,rawY,sz,kolor, "filled")
    plot([xx-buff2,xx+buff2],[mean(rawY),mean(rawY)],'Color',kolor,'LineWidth',lw)
    xx = i+buff2;
    kolor = Color('dodgerblue');
    rawY = SD.c_occ(:,i);
    rawY(isnan(rawY)) = [];
    x = shuffle_data(linspace(xx-buff, xx+buff,length(rawY)));
    scatter(x,rawY,sz,kolor, "filled")
    plot([xx-buff2,xx+buff2],[mean(rawY),mean(rawY)],'Color',kolor,'LineWidth',lw)
end
% stats
y = rangeLine(fig, 0, false);
y2 = rangeLine(fig, 3, false);
for i = 1:num.exp
    if h(i)
        plot([i-buff2, i+buff2],[y,y],'linewidth', 1.5, 'Color',foreColor)
        scatter(i,y2,100,foreColor,'*')
    end
end
% formatting
formatFig(fig, blkbgd);
ylabel('flies in outer ring (%)')
set(gca, 'XTick',1:num.exp, 'XTickLabel',strrep(SD.name,'_',' '))
title_str = [num2str(temps(1)) ' to ' num2str(temps(2))];
title(title_str)
save_figure(fig,[fig_dir 'outer ring occupancy from ' title_str],fig_type);


% timecourse for the same temperature data:
scaler = 1; plot_err = 1;
fig = getfig('',true,[480 680]); hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    x = grouped(i).ring.temps;
    YC = grouped(i).ring.decreasing.raw;
    YH = grouped(i).ring.increasing.raw;
    % cooling
    y = mean(YC,2,'omitnan')*scaler;
    y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '--')
    % heating
    y = mean(YH,2,'omitnan')*scaler;
    y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '-')        
end
y = rangeLine(fig,2,true);
plot(temps, [y,y],'color', foreColor,'linewidth', 1.5)
formatFig(fig, blkbgd);
ylabel('flies in outer ring (%)')
xlabel('temperature (\circC)')

save_figure(fig,[fig_dir 'outer ring occupancy tuning curve'],fig_type);


% STATISTICS: 
% Differences between heating and cooling? ANOVA of absolute area between
% the curves -- asking if it differs significantly from zero and then
% asking if it differs between the different experiment conditions: 

% Is there a difference in the percent of flies in the outer rim during
% heating and cooling for each of the conditions? 

% How do the COOLING temperature dependent lines differ? Slope based analysis
temp = []; occ = []; exp = [];
%format the data for the statistical test
for i = 1:num.exp
    x = repmat((grouped(i).ring.temps'),[num.trial(i),1]);
    temp = [temp; x];
    y = grouped(i).ring.decreasing.raw(:);
    occ  = [occ; y];
    exp = [exp; i*ones(size(y))];
end
% run the ancova test
[~,~,~,stats] = aoctool(temp,occ,exp);
%run the multicomparison to look at differences in slope
c = multcompare(stats,"Estimate", "slope","CriticalValueType","bonferroni","Display","on");
% display the statistical differences: 
thresh = 0.05;
cPlot = []; pPlot = [];
for i = 1:size(c,1)
    cPlot(c(i,1),c(i,2)) = c(i,6);
    cPlot(c(i,2),c(i,1)) = c(i,6);
    if c(i,6)<=thresh
        pPlot = [pPlot; c(i,1),c(i,2); c(i,2),c(i,1)];
    end
end
fig = getfig('',1,[705 649]);
% imagesc(cPlot);
% axis square equal
xlim([0.5,num.exp+0.5])
ylim([0.5,num.exp+0.5])
hold on
scatter(pPlot(:,1), pPlot(:,2), 100, foreColor, '*')
v_line(0.5:1:num.exp+0.5,'grey','-',0.5)
h_line(0.5:1:num.exp+0.5,'grey','-',0.5)
formatFig(fig, blkbgd);
names = strrep({grouped(:).name},'_',' ');
set(gca,'XTick',1:1:num.exp,'YTick',1:1:num.exp,'XTickLabel',names,'YTickLabel',names)
set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
save_figure(fig,[fig_dir 'outer ring occupancy cooling slope stats'],fig_type);

% How do the COOLING temperature dependent lines differ? Intercept based analysis
temp = []; occ = []; exp = [];
% format the data for the statistical test
for tt = 1:num.exp
    i = expOrder(tt);
    x = repmat((grouped(i).ring.temps'),[num.trial(i),1]);
    temp = [temp; x];
    y = grouped(i).ring.decreasing.raw(:);
    occ  = [occ; y];
    exp = [exp; tt*ones(size(y))];
end
%run the ancova test
[~,~,~,stats] = aoctool(temp,occ,exp);
%run the multicomparison to look at differences in slope
c = multcompare(stats,"Estimate", "intercept","CriticalValueType","bonferroni","Display","on");
% display the statistical differences: 
thresh = 0.05;
cPlot = []; pPlot = [];
for i = 1:size(c,1)
    cPlot(c(i,1),c(i,2)) = c(i,6);
    cPlot(c(i,2),c(i,1)) = c(i,6);
    if c(i,6)<=thresh
        pPlot = [pPlot; c(i,1),c(i,2); c(i,2),c(i,1)];
    end
end
fig = getfig('',1,[705 649]);
% imagesc(cPlot);
% axis square equal
lims = [0.5,num.exp+0.5];
xlim(lims)
ylim(lims)
hold on
scatter(pPlot(:,1), pPlot(:,2), 100, foreColor, '*')
v_line(0.5:1:num.exp+0.5,'grey','-',0.5)
h_line(0.5:1:num.exp+0.5,'grey','-',0.5)
formatFig(fig, blkbgd);
names = strrep({grouped(expOrder).name},'_',' ');
set(gca,'XTick',1:1:num.exp,'YTick',1:1:num.exp,'XTickLabel',names,'YTickLabel',names)
set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
plot(lims, lims)
save_figure(fig,[fig_dir 'outer ring occupancy cooling intercept stats'],fig_type);

% How do the WARMING temperature dependent lines differ? Slope based analysis
temp = []; occ = []; exp = [];
%format the data for the statistical test
for i = 1:num.exp
    x = repmat((grouped(i).ring.temps'),[num.trial(i),1]);
    temp = [temp; x];
    y = grouped(i).ring.increasing.raw(:);
    occ  = [occ; y];
    exp = [exp; i*ones(size(y))];
end
%run the ancova test
[~,~,~,stats] = aoctool(temp,occ,exp);
%run the multicomparison to look at differences in slope
c = multcompare(stats,"Estimate", "slope","CriticalValueType","bonferroni","Display","on");
% display the statistical differences: 
thresh = 0.05;
cPlot = []; pPlot = [];
for i = 1:size(c,1)
    cPlot(c(i,1),c(i,2)) = c(i,6);
    cPlot(c(i,2),c(i,1)) = c(i,6);
    if c(i,6)<=thresh
        pPlot = [pPlot; c(i,1),c(i,2); c(i,2),c(i,1)];
    end
end
fig = getfig('',1,[705 649]);
% imagesc(cPlot);
% axis square equal
xlim([0.5,num.exp+0.5])
ylim([0.5,num.exp+0.5])
hold on
scatter(pPlot(:,1), pPlot(:,2), 100, foreColor, '*')
v_line(0.5:1:num.exp+0.5,'grey','-',0.5)
h_line(0.5:1:num.exp+0.5,'grey','-',0.5)
formatFig(fig, blkbgd);
names = strrep({grouped(:).name},'_',' ');
set(gca,'XTick',1:1:num.exp,'YTick',1:1:num.exp,'XTickLabel',names,'YTickLabel',names)
set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
save_figure(fig,[fig_dir 'outer ring occupancy warming slope stats'],fig_type);

% How do the WARMING temperature dependent lines differ? Intercept based analysis
temp = []; occ = []; exp = [];
%format the data for the statistical test
for i = 1:num.exp
    x = repmat((grouped(i).ring.temps'),[num.trial(i),1]);
    temp = [temp; x];
    y = grouped(i).ring.increasing.raw(:);
    occ  = [occ; y];
    exp = [exp; i*ones(size(y))];
end
%run the ancova test
[~,~,~,stats] = aoctool(temp,occ,exp);
%run the multicomparison to look at differences in slope
c = multcompare(stats,"Estimate", "intercept","CriticalValueType","bonferroni","Display","on");
% display the statistical differences: 
thresh = 0.05;
cPlot = []; pPlot = [];
for i = 1:size(c,1)
    cPlot(c(i,1),c(i,2)) = c(i,6);
    cPlot(c(i,2),c(i,1)) = c(i,6);
    if c(i,6)<=thresh
        pPlot = [pPlot; c(i,1),c(i,2); c(i,2),c(i,1)];
    end
end
fig = getfig('',1,[705 649]);
% imagesc(cPlot);
% axis square equal
xlim([0.5,num.exp+0.5])
ylim([0.5,num.exp+0.5])
hold on
scatter(pPlot(:,1), pPlot(:,2), 100, foreColor, '*')
v_line(0.5:1:num.exp+0.5,'grey','-',0.5)
h_line(0.5:1:num.exp+0.5,'grey','-',0.5)
formatFig(fig, blkbgd);
names = strrep({grouped(:).name},'_',' ');
set(gca,'XTick',1:1:num.exp,'YTick',1:1:num.exp,'XTickLabel',names,'YTickLabel',names)
set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
save_figure(fig,[fig_dir 'outer ring occupancy warming intercept stats'],fig_type);



%% FIGURE & STATS: Hysteresis for each genotype / trial
clearvars('-except',initial_vars{:})
fig_dir = createFolder([saveDir, 'Stats/']);

SD = struct;
for ii = 1:num.exp
    i = expOrder(ii);
    SD(ii).name = grouped(i).name;
    SD(ii).kolor = grouped(i).color;
    SD(ii).p = [];
    y1 = grouped(i).ring.decreasing.raw;
    y2 = grouped(i).ring.increasing.raw;
    for t = 1:size(y1,1) % for each temperature bin
       [~,SD(ii).p(t)] = ttest(y1(t,:), y2(t,:));
    end
    % Bonferonni correction:
    alpha = 0.05;
    m = size(y1,1); % number of temperature comparisions
    p_limit = alpha/m;
    SD(ii).h = SD(ii).p<=p_limit;
    % cumulative difference between heating and cooling
end

r = 10;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2:r;

% timecourse for the same temperature data:
scaler = 1; plot_err = 1; y_lim = [];
fig = getfig('',true,[480 680]); 
%significance plot
subplot(r,c,sb(1).idx)
hold on
for i = 1:num.exp
    x = grouped(i).ring.temps;
    y = i *ones(size(x));
    y(~SD(i).h) = nan;
    plot(x,y,'color', SD(i).kolor,'linewidth', 1.5)
end
y_lim = [y_lim; ylim];
% plot temp tuning curve
subplot(r,c,sb(2).idx)
hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    x = grouped(i).ring.temps;
    YC = grouped(i).ring.decreasing.raw;
    YH = grouped(i).ring.increasing.raw;
    % cooling
    y = mean(YC,2,'omitnan')*scaler;
    y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '--')
    % heating
    y = mean(YH,2,'omitnan')*scaler;
    y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '-')        
end
y_lim = [y_lim; ylim];
ylabel('flies in outer ring (%)')
xlabel('temperature (\circC)')

formatFig(fig, blkbgd,[r,c],sb);
subplot(r,c,sb(2).idx)
ylim = [min(y_lim(:,1)) max(y_lim(:,2))]; 
subplot(r,c,sb(1).idx)
ylim = [min(y_lim(:,1)) max(y_lim(:,2))]; 
set(gca,'XColor','none', 'YColor','none')

save_figure(fig,[fig_dir 'outer ring occupancy tuning curve with running heat cool stats'],fig_type);


% FIGURE: Cumulative difference in outer ring occupancy
LW = 0.75;
buff = 0.2;
SZ = 50;
r = 1; %rows
c = 3; %columns
plot_err = false;
plotSig = true; %plot significance stars
[foreColor,~] = formattingColors(blkbgd);

fig = getfig('',true,[480 680]); 
hold on
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    y = (grouped(i).ring.decreasing.raw-grouped(i).ring.increasing.raw);
    plotY = sum(y,1,'omitnan');
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
    scatter(x,plotY,SZ,kolor,"filled","o")
    plot([ii-buff,ii+buff],[mean(plotY),mean(plotY)],'color',foreColor,'LineWidth',2)
end
xlim([0.5,num.exp+0.5])
h_line(0,foreColor,':',1)
ylabel('cumulative difference (mm)')
formatFig(fig,blkbgd);
set(gca, 'xcolor', 'none')

% STATS: are the means of any groups different from zero?
[p, mlt, id] = deal([]);
for ii = 1:num.exp
    i = expOrder(ii);
    y = (grouped(i).ring.decreasing.raw-grouped(i).ring.increasing.raw);
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
    for ii = 1:num.exp
        if h(ii)
            scatter(ii,y_pos,100,foreColor,'*')
        end
    end
end

% Multicompare across the groups for significance
% STATS:
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
save_figure(fig,[fig_dir expGroup ' outer ring occupancy cumulative hysteresis summary'],fig_type);



%% TODO: incorporate the dynamic structure from below to compare different types of experimental data
plot_err = true;
autoLim = true;
xlim_auto = true; % change the time range for the x axis
time_limits = [0,900]; % time limits if manual control over x-axis range
nMax =  num.exp;%
[~,backColor] = formattingColors(blkbgd); %get background colors

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end] = PlotParamSelection(true);
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

   % selected parameter time course
    subplot(r,c,sb(2).idx); hold on
        switch dType
            case 1 % single trial lines
                for trial = 1:num.trial(i)
                    y = smooth(grouped(i).(pName).all(:,trial),sSpan, 'moving')*scaler;
                    plot(x,y,'LineWidth',LW,'Color',kolor)
                end
            case {2, 3} % avg line
                y = smooth(grouped(i).(pName).avg,sSpan, 'moving')*scaler;
                plot(x,y,'LineWidth',LW,'Color',kolor)
        end

    %temp vs dependent variable tuning curve
    subplot(r,c,sb(3).idx); hold on
    
     switch dType
         case 1 % single trial lines
            for trial = 1:num.trial(i)
                if strcmp(pName, 'dist')
                    x = grouped(i).(pName).distavgbytemp(:,1);
                    rawY = [grouped(i).increasing.all(:,trial),grouped(i).decreasing.all(:,trial)];
                else
                    x = grouped(i).(pName).temps;
                    rawY = [grouped(i).(pName).increasing.raw(:,trial),grouped(i).(pName).decreasing.raw(:,trial)];
                end
                y = mean(rawY,2,'omitnan')*scaler;
                plot(x,y,'color',kolor,'linewidth',1.25)
            end

         case 2 % avg lines (combined heating and cooling)
            if strcmp(pName, 'dist')
                x = grouped(i).(pName).distavgbytemp(:,1);
                rawY = [grouped(i).increasing.all,grouped(i).decreasing.all];
            else
                x = grouped(i).(pName).temps;
                rawY = [grouped(i).(pName).increasing.raw,grouped(i).(pName).decreasing.raw];
            end
            y = mean(rawY,2,'omitnan')*scaler;
            y_err = std(rawY,0,2,'omitnan')*scaler;
            plot(x,y,'color',kolor,'linewidth',1.25)
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);

         case 3 % separated heating and cooling
            if strcmp(pName, 'dist')
                x = grouped(i).(pName).distavgbytemp(:,1);
                YC = grouped(i).decreasing.all;
                YH = grouped(i).increasing.all;
            else
                x = grouped(i).(pName).temps;
                YC = grouped(i).(pName).decreasing.raw;
                YH = grouped(i).(pName).increasing.raw;
            end
            % cooling
            y = mean(YC,2,'omitnan')*scaler;
            y_err = std(YC,0,2,'omitnan')*scaler;
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
            plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '--')
            % heating
            y = mean(YH,2,'omitnan')*scaler;
            y_err = std(YH,0,2,'omitnan')*scaler;
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
            plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '-')          
     end
     dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",'none')

% distance
subplot(r,c,sb(2).idx)
ylabel(y_lab)
xlabel('time (min)')
set(gca,'ydir',y_dir)
% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel(y_lab)
xlabel('temp (\circC)')
% if ~autoLim
%     ylim(dt_lim)
% end
h_line(nullD,'grey',':',2) %36.2
set(gca,'ydir',y_dir)

if ~xlim_auto
    subplot(r,c,sb(1).idx)
    set(gca, 'xlim', time_limits)
    subplot(r,c,sb(2).idx)
    set(gca, 'xlim', time_limits)
    % ylim([0,100])
    % ylim([20,100])
end

% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[fig_dir 'Timecourse summary ' title_str],fig_type);


%% FIGURE & STATS: Linear regression for each trial for cooling outer ring
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd); %get background colors
fig_dir = createFolder([saveDir, 'Stats/']);

% for each trial, pull out all the data within the cooling regions and fit
% a linear model to those data points (not binned): 
stats = struct;
for i = 1:num.exp
    [stats(i).R2, stats(i).slp, stats(i).RMSE, stats(i).intercept] = deal(nan(num.trial(i),1));
    tp = getTempTurnPoints(data(i).temp_protocol);
    raw_y = grouped(i).ring.all; % pull all the raw data for this parameter
    roi = tp.DownROI; % find the region for all cooling points
    for trial = 1:num.trial(i)
        temp = grouped(i).temp(roi);
        y = raw_y(roi,trial);
        mdl = fitlm(temp,y);
        stats(i).R2(trial) = mdl.Rsquared.Ordinary;
        stats(i).slp(trial) = table2array(mdl.Coefficients(2,1));
        stats(i).intercept(trial) = table2array(mdl.Coefficients(1,1));
        stats(i).RMSE(trial) = mdl.RMSE;
    end
    disp(['Finished ' num2str(i)])
end

% PLOT LINEAR FIT PARAMETERS
r = 1;
c = 3;
buff = 0.2;
buff2 = 0.3;
fig = getfig('',1); 
% slope
subplot(r,c,1)
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    x = shuffle_data(linspace(exp-buff,exp+buff,num.trial(i)));
    y = stats(i).slp;
    scatter(x,y,35, grouped(i).color,'filled')
    y_avg = mean(y);
    plot([exp-buff2,exp+buff2],[y_avg, y_avg],'Color', foreColor,'linewidth', 1.5)
end
h_line(0,'grey','--')
ylabel('Slope of temp - outer ring occupancy')
% R2
subplot(r,c,2)
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    x = shuffle_data(linspace(exp-buff,exp+buff,num.trial(i)));
    y = stats(i).R2;
    scatter(x,y,35, grouped(i).color,'filled')
    y_avg = mean(y);
    plot([exp-buff2,exp+buff2],[y_avg, y_avg],'Color', foreColor,'linewidth', 1.5)
end
% h_line(0,'grey','--')
ylabel('R2 value')
% RMSE value
subplot(r,c,3)
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    x = shuffle_data(linspace(exp-buff,exp+buff,num.trial(i)));
    y = stats(i).RMSE;
    scatter(x,y,35, grouped(i).color,'filled')
    y_avg = mean(y);
    plot([exp-buff2,exp+buff2],[y_avg, y_avg],'Color', foreColor,'linewidth', 1.5)
end
% h_line(0,'grey','--')
ylabel('RMSE')
formatFig(fig,blkbgd,[r,c]);
for i = 1:c
    subplot(r,c,i)
    set(gca, 'xcolor', 'none')
end
save_figure(fig,[fig_dir 'outer ring individual fit cooling slope fits'],fig_type);

% STATS: n-way anova of slope values for each group
slopes = [];
for i = 1:num.exp
    exp = expOrder(i);
    slopes = [slopes; i*ones(num.trial(exp),1), stats(exp).slp];
end

aov = anova(slopes(:,1), slopes(:,2));
c = multcompare(aov,"CriticalValueType","bonferroni");
c = table2array(c);
names = strrep({grouped(expOrder).name},'_',' ');
fig = getMultiCompSignTable(c, 1:num.exp, blkbgd, 0.05, names);
title('outer ring individual fit cooling slope stats')

save_figure(fig,[fig_dir 'outer ring individual fit cooling slope stats'],fig_type);




%% FIGURES & STATS: relative and absolute time in specific regions during warming vs cooling (temp independent) 
clearvars('-except',initial_vars{:})
fig_dir = createFolder([saveDir, 'Stats/']);
[foreColor,~] = formattingColors(blkbgd); %get background colors
[title_str, pName,~,~,~,scaler,~,~] = PlotParamSelection(false,true);

% ------------ DATA EXTRACTION ------------ 
% for each trial, pull out all the data within the cooling regions and
% warming regions for comparision:
pd = struct;
[pd.all, pd.c_avg,  pd.h_avg] = deal(nan(max(num.trial),num.exp));
[pd.avg, pd.sem] = deal(nan(num.exp,1));
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).temp_protocol);
    raw_y = grouped(i).(pName).all; % pull all the raw data for this parameter (this is the percent of flies in this region...)
    raw_y = (raw_y.*data(i).T.NumFlies')./100;  % total number of flies in the region
    roi_c = tp.DownROI; % find the region for all cooling points
    roi_w = tp.UpROI; % find the region for all cooling points
    c_count = sum(raw_y(roi_c,:))./data(i).T.NumFlies'; % avg frames per fly in the region (cooling)
    h_count = sum(raw_y(roi_w,:))./data(i).T.NumFlies'; % avg frames per fly in the region (warming)
    c_avg = ((c_count./length(roi_c)))*100; % each fly spent this % of time in the region
    h_avg = ((h_count./length(roi_w)))*100; % each fly spent this % of time in the region
    y = h_avg - c_avg; 
    pd.c_avg(1:num.trial(i),i) = c_avg;
    pd.h_avg(1:num.trial(i),i) = h_avg;
    pd.all(1:num.trial(i),i) = y;
    pd.avg(i) = mean(y,'omitnan');
    pd.sem(i) = std(y)/sqrt(num.trial(i));
end
% note: PD structure is not is experiment order, it is in grouped order.

% ------------ STATISTICAL ANALYSES ------------ 

% ---------------------------------------------------------
% ANOVA with Interactions:
% Comparison between peak % time in region during warming & cooling across conditions
y_cool = pd.c_avg(:);
y_warm = pd.h_avg(:);
exp = repmat({grouped(:).name},[size(pd.c_avg,1),1]);
exp = reshape(exp, [numel(exp), 1]);
c_regime = repmat({'cooling'},size(exp));
w_regime = repmat({'warming'},size(exp)); 
% concatenate the vectors
temp_regime = [c_regime; w_regime];
experiments = [exp; exp];
y = [y_cool; y_warm];
% remove nans
loc = isnan(y);
temp_regime(loc) = [];
y(loc) = [];
experiments(loc) = [];
% test for interaction between temperature regime (heat vs cool) and experiment groups
tbl = table(experiments,temp_regime,y,VariableNames=["experiments" "temp_regime" "Y"]);
aovInteraction = anova(tbl,"Y ~ experiments + temp_regime + experiments:temp_regime");
disp('ANOVA with interaction test for data:')
disp(aovInteraction)

% ---------------------------------------------------------

% One-way anova of % time spent in the region during warming | cooling across
% experiment groups
names = strrep({grouped(:).name},'_',' ');
for tt = 1:2
    switch tt
        case 1
            y_type = 'cooling';
            y_raw = y_cool;
        case 2 
            y_type = 'warming';
            y_raw = y_warm;
    end
    [~,~,stats] = anova1(y_raw, exp, 'off'); 
    [c, ~,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni",'Display','off');
    cold_stats = array2table(c, "VariableNames", ["Group", "Control Group", "Lower Limit","Difference","Upper Limit","P-value"]);
    cold_stats.("Group") = gnames(cold_stats.("Group"));
    cold_stats.("Control Group") = gnames(cold_stats.("Control Group"));
    fprintf(['\n \n ' y_type '\n'])
    disp(cold_stats)
    
    % PLOT anova mulitple comparisons matrix
    fig = getMultiCompSignTable(c, expOrder, blkbgd, 0.05, names,true);
    title(['Time in ' title_str ' during ' y_type],'Color',foreColor)
    save_figure(fig,[fig_dir title_str ' differences in time during ' y_type],fig_type);
end

% ---------------------------------------------------------

% One-way anova in the warming/cooling assymetry behavior between exp groups
names = strrep({grouped(:).name},'_',' ');
exp = repmat({grouped(:).name},[size(pd.all,1),1]);
y = pd.all(:); %pd is not expOrder but grouped order
exp = reshape(exp,[numel(y),1]);
loc = isnan(y);
y(loc) = [];
exp(loc) = [];
[~,~,stats] = anova1(y, exp, 'off'); 
[c, ~,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni",'Display','off');
group_stats = array2table(c, "VariableNames", ["Group", "Control Group", "Lower Limit","Difference","Upper Limit","P-value"]);
group_stats.("Group") = gnames(group_stats.("Group"));
group_stats.("Control Group") = gnames(group_stats.("Control Group"));
fprintf('\n \n Grouped data: \n')
disp(group_stats)

% PLOT anova mulitple comparisons matrix
fig = getMultiCompSignTable(c, expOrder, blkbgd, 0.05, names,true);
title(['Assymetry btwn W&C in ' title_str],'Color',foreColor)
save_figure(fig,[fig_dir title_str ' asymetry across W and C percent time'],fig_type);

% ---------------------------------------------------------

% Stats: are there significant differences during heating and cooling? (for each genotype)
alpha = 0.05/num.exp; %bonferonni MCC
h = ttest((pd.c_avg), (pd.h_avg),'Alpha',alpha);
disp('Groups with significant differences between heating and cooling:')
disp({grouped(logical(h)).name}')

% ---------------------------------------------------------

% Group FIGURE: 
r = 1; c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;
buff = 0.2; buff2 = 0.3;
fig = getfig('',1); 
% heating vs cooling
subplot(r,c,sb(1).idx)
hold on
x_ticks = [];
for i = 1:num.exp
    exp = expOrder(i);
    x1 = (i-buff2)*ones(num.trial(exp),1);
    y1 = rmnan(pd.c_avg(:,exp))*scaler;
    x2 = (i+buff2)*ones(num.trial(exp),1);
    y2 = rmnan(pd.h_avg(:,exp))*scaler;
    x = [i-buff2; i+buff2];
    y = [mean(y1); mean(y2)];
    plot([x1,x2]',[y1,y2]', 'color', grouped(exp).color,'linewidth', 1)
    plot(x,y, 'color', foreColor,'linewidth', 2.5)
    x_ticks = [x_ticks; x];
    % y_avg = [pd.avg(exp), pd.avg(exp)];
    % % errorbar
    % errorbar(i,y_avg(1), pd.sem(exp),'Color',foreColor,'linewidth', 1.5,'Marker','square','MarkerFaceColor',foreColor)
    % % plot([i-buff2,i+buff2],y_avg,'Color', foreColor,'linewidth', 1.5)
end
xlim([1-(buff2*2),num.exp+(2*buff2)])
set(gca, 'xtick', x_ticks, 'XTickLabel', repmat({'C','H'},[1, num.exp]))
ylabel(['Avg time spent in ' title_str ' (%)'])

% difference:
subplot(r,c,sb(2).idx)
hold on
for i = 1:num.exp
    exp = expOrder(i);
    x = shuffle_data(linspace(i-buff,i+buff,num.trial(exp)));
    y = rmnan(pd.all(:,exp))*scaler;
    scatter(x,y,35, grouped(exp).color,'filled')
    y_avg = [pd.avg(exp), pd.avg(exp)];
    y_avg = y_avg*scaler;
    % errorbar
    errorbar(i,y_avg(1), pd.sem(exp)*scaler,'Color',foreColor,'linewidth', 1.5,'Marker','square','MarkerFaceColor',foreColor)
    % plot([i-buff2,i+buff2],y_avg,'Color', foreColor,'linewidth', 1.5)
end
xlim([0,num.exp+1])
h_line(0,'grey', '--')
ylabel(['Difference in time spent in ' title_str ' btwn C vs H (%)'])

formatFig(fig, blkbgd,[r,c], sb);
subplot(r,c,sb(2).idx)
set(gca, 'xcolor', 'none')

% plot the statistics
H = h(expOrder);
subplot(r,c,sb(1).idx)
y = rangeLine(fig,2,true);
x = 1:num.exp;
Y = y*ones(size(x));
x(~H) = [];
Y(~H) = [];
scatter(x,Y,150,foreColor, 'filled', 'pentagram')

save_figure(fig,[fig_dir 'Percent time spent in ' title_str ' H & C'],fig_type);

% --------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------
% ------------------------------- Absolute time comparisons --------------------------------------------

% Show as total time spent in a region, rather than the proportion of time (mostly
% interesting for the different temperature rate experiments) 
plotData = []; [c_avg, h_avg] = deal(nan(max(num.trial), num.exp));
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).temp_protocol);
    % find the total time in heating or cooling:
    dur_c = (length(tp.DownROI)/fps)/60; % total time cooling (min)
    dur_h = (length(tp.UpROI)/fps)/60; % total time cooling (min)
    % get the translated time: 
    plotData(i).c = ((rmnan(pd.c_avg(:,i))*dur_c)/100)*scaler; % minutes in the region during cooling
    plotData(i).h = ((rmnan(pd.h_avg(:,i))*dur_h)/100)*scaler; % minutes in the region during warming
    plotData(i).diff = plotData(i).h-plotData(i).c;
    plotData(i).sem = std(plotData(i).diff)./sqrt(num.trial(i));
    c_avg(1:num.trial(i),i) = plotData(i).c;
    h_avg(1:num.trial(i),i) = plotData(i).h;
end

% STATISTICAL ANALYSIS: 
% Stats: are there significant differences during heating and cooling? (for each genotype)
% **SHOULD** be identical to the data above
alpha = 0.05/num.exp; %bonferonni MCC
h = ttest((c_avg), (h_avg),'Alpha',alpha);

% FIGURE: 
r = 1; c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;
buff = 0.2; buff2 = 0.3;
fig = getfig('',1); 
% heating vs cooling
subplot(r,c,sb(1).idx)
hold on
x_ticks = [];
for i = 1:num.exp
    exp = expOrder(i);
    x1 = (i-buff2)*ones(num.trial(exp),1);
    y1 = plotData(exp).c;
    x2 = (i+buff2)*ones(num.trial(exp),1);
    y2 = plotData(exp).h;
    x = [i-buff2; i+buff2];
    y = [mean(y1); mean(y2)];
    plot([x1,x2]',[y1,y2]', 'color', grouped(exp).color,'linewidth', 1)
    plot(x,y, 'color', foreColor,'linewidth', 2.5)
    x_ticks = [x_ticks; x];
    % y_avg = [pd.avg(exp), pd.avg(exp)];
    % % errorbar
    % errorbar(i,y_avg(1), pd.sem(exp),'Color',foreColor,'linewidth', 1.5,'Marker','square','MarkerFaceColor',foreColor)
    % % plot([i-buff2,i+buff2],y_avg,'Color', foreColor,'linewidth', 1.5)
end
xlim([1-(buff2*2),num.exp+(2*buff2)])
set(gca, 'xtick', x_ticks, 'XTickLabel', repmat({'C','H'},[1, num.exp]))
ylabel(['Avg time spent in ' title_str ' (min)'])

% difference:
subplot(r,c,sb(2).idx)
hold on
for i = 1:num.exp
    exp = expOrder(i);
    x = shuffle_data(linspace(i-buff,i+buff,num.trial(exp)));
    y = plotData(exp).diff;
    scatter(x,y,35, grouped(exp).color,'filled')
    y_avg = [mean(y), mean(y)];
    % errorbar
    errorbar(i,y_avg(1), plotData(exp).sem,'Color',foreColor,'linewidth', 1.5,'Marker','square','MarkerFaceColor',foreColor)
    % plot([i-buff2,i+buff2],y_avg,'Color', foreColor,'linewidth', 1.5)
end
xlim([0,num.exp+1])
h_line(0,'grey', '--')
ylabel('Time (min) spent in heating (+) vs cooling (-)')

formatFig(fig, blkbgd,[r,c], sb);
subplot(r,c,sb(2).idx)
set(gca, 'xcolor', 'none')

% plot the statistics: TODO
H = h(expOrder);
subplot(r,c,sb(1).idx)
y = rangeLine(fig,2,true);
x = 1:num.exp;
Y = y*ones(size(x));
x(~H) = [];
Y(~H) = [];
scatter(x,Y,150,foreColor, 'filled', 'pentagram')
save_figure(fig,[fig_dir 'Time spent in ' title_str ' during H & C'],fig_type);

% ---------------------------------------------------------

% One-way anova in the warming/cooling assymetry behavior between exp groups
y = []; exp = [];
for i = 1:num.exp
    y = [y; plotData(i).diff];
    exp = [exp; repmat({grouped(i).name},size(plotData(i).diff))];
end
[p,~,stats] = anova1(y, exp, 'off'); 
[c, m,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni",'Display','off');
group_stats = array2table(c, "VariableNames", ["Group", "Control Group", "Lower Limit","Difference","Upper Limit","P-value"]);
group_stats.("Group") = gnames(group_stats.("Group"));
group_stats.("Control Group") = gnames(group_stats.("Control Group"));
fprintf('\n \n Grouped data: \n')
disp(group_stats)

% PLOT anova mulitple comparisons matrix
fig = getMultiCompSignTable(c, expOrder, blkbgd, 0.05, names,true);
title(['time assymetry btwn W&C in ' title_str],'Color',foreColor)
save_figure(fig,[fig_dir title_str ' asymetry across W and C absolute time'],fig_type);


%% FIGURE + STATS: dynamic variable linearity comparision between experiment types
% fully grouped models -- not run on individual trials but fit across all trials' data
clearvars('-except',initial_vars{:})
fig_dir = createFolder([saveDir, 'Stats/']);
[foreColor,~] = formattingColors(blkbgd); %get background colors
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end] = PlotParamSelection(false);
a = questdlg('Increasing or decreasing?','','Warming', 'Cooling', 'Cancel', 'Cooling');
switch a
    case 'Warming'
        d_type = 'increasing';
        l_type = '-';
        x_dir = 'normal';
        roi_name = 'UpROI';
    case 'Cooling'
        d_type = 'decreasing';
        l_type = '--';
        x_dir = 'reverse';
        roi_name = 'DownROI';
    case {'Cancel',''}
        disp('Canceled selection')
        return
end

% switch questdlg('Analyze the slope or intercept?','','Slope', 'Intercept', 'Cancel', 'Slope')
%     case 'Slope'
%         a_type = 'slope';
%         b_type = 'slp';
%     case 'Intercept'
%         a_type = 'intercept';
%         b_type = 'intercept';
%     case {'Cancel',''}
%         disp('Canceled selection')
%         return
% end

% PLOT RAW DATA BEING COMPARED:
plot_err = 1;
fig = getfig('',true,[480 680]); hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    x = grouped(i).(pName).temps;
    YC = grouped(i).(pName).(d_type).raw;
    y = mean(YC,2,'omitnan')*scaler;
    y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '--')
end
% y = rangeLine(fig,2,true);
% plot(temps, [y,y],'color', foreColor,'linewidth', 1.5)
formatFig(fig, blkbgd);
ylabel(y_lab)
xlabel('temperature (\circC)')
set(gca, 'YDir',y_dir,'XDir', x_dir)
save_figure(fig,[fig_dir title_str ' ' d_type ' temp tuning curve'],fig_type);

% How do the selected temperature dependent lines differ?
temp = []; y_data = []; exp = [];
%format the data for the statistical test
for tt = 1:num.exp
    i = expOrder(tt);
    x = repmat((grouped(i).(pName).temps'),[num.trial(i),1]);
    temp = [temp; x];
    y = grouped(i).(pName).(d_type).raw(:);
    y_data  = [y_data; y];
    exp = [exp; tt*ones(size(y))];
end
% run the ancova test
[~,~,~,stats] = aoctool(temp,y_data,exp,0.05,'','','', "off");
for ii = 1:2
    switch ii
        case 1
            param_type = 'slope';
        case 2 
            param_type = 'intercept';
    end
    %run the multicomparison to look at differences in slope
    c = multcompare(stats,"Estimate", param_type,"CriticalValueType","bonferroni","Display","off");
    names = strrep({grouped(expOrder).name},'_',' ');
    fig = getMultiCompSignTable(c, 1:num.exp, blkbgd, 0.05, names);
    title([title_str ' ' d_type ' temp ' param_type ' stats'],'color', foreColor)
    save_figure(fig,[fig_dir title_str ' ' d_type ' temp ' param_type ' stats'],fig_type);
end

% -------------------------------------------------------------------------------------------------------------
% -----------------------------    Individual fits for the data     -----------------------------------
% -------------------------------------------------------------------------------------------------------------

% for each trial, pull out all the data within the cooling regions and fit
% a linear model to those data points (not binned): 
stats = struct;
for i = 1:num.exp
    [stats(i).R2, stats(i).slp, stats(i).RMSE, stats(i).intercept] = deal(nan(num.trial(i),1));
    tp = getTempTurnPoints(data(i).temp_protocol);
    raw_y = grouped(i).(pName).all; % pull all the raw data for this parameter
    roi = tp.(roi_name); % find the region for all cooling points
    for trial = 1:num.trial(i)
        temp = grouped(i).temp(roi);
        y = raw_y(roi,trial);
        mdl = fitlm(temp,y);
        stats(i).R2(trial) = mdl.Rsquared.Ordinary;
        stats(i).slp(trial) = table2array(mdl.Coefficients(2,1));
        stats(i).intercept(trial) = table2array(mdl.Coefficients(1,1));
        stats(i).RMSE(trial) = mdl.RMSE;
    end
    disp(['Finished ' num2str(i)])
end

% PLOT LINEAR FIT PARAMETERS
r = 1;
c = 4;
buff = 0.2;
buff2 = 0.3;
fig = getfig('',1); 
% slope
subplot(r,c,1)
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    x = shuffle_data(linspace(exp-buff,exp+buff,num.trial(i)));
    y = stats(i).slp.*scaler;
    scatter(x,y,35, grouped(i).color,'filled')
    y_avg = mean(y);
    plot([exp-buff2,exp+buff2],[y_avg, y_avg],'Color', foreColor,'linewidth', 1.5)
end
h_line(0,'grey','--')
ylabel([ 'Slope of temp - ' title_str ' occupancy'])
% intercept
subplot(r,c,2)
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    x = shuffle_data(linspace(exp-buff,exp+buff,num.trial(i)));
    y = stats(i).intercept.*scaler;
    scatter(x,y,35, grouped(i).color,'filled')
    y_avg = mean(y);
    plot([exp-buff2,exp+buff2],[y_avg, y_avg],'Color', foreColor,'linewidth', 1.5)
end
h_line(0,'grey','--')
ylabel(['Intercept of temp - ' title_str ' occupancy'])
% R2
subplot(r,c,3)
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    x = shuffle_data(linspace(exp-buff,exp+buff,num.trial(i)));
    y = stats(i).R2;
    scatter(x,y,35, grouped(i).color,'filled')
    y_avg = mean(y);
    plot([exp-buff2,exp+buff2],[y_avg, y_avg],'Color', foreColor,'linewidth', 1.5)
end
ylabel('R2 value')
% RMSE value
subplot(r,c,4)
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    x = shuffle_data(linspace(exp-buff,exp+buff,num.trial(i)));
    y = stats(i).RMSE;
    scatter(x,y,35, grouped(i).color,'filled')
    y_avg = mean(y);
    plot([exp-buff2,exp+buff2],[y_avg, y_avg],'Color', foreColor,'linewidth', 1.5)
end
ylabel('RMSE')
formatFig(fig,blkbgd,[r,c]);
for i = 1:c
    subplot(r,c,i)
    set(gca, 'xcolor', 'none')
end
save_figure(fig,[fig_dir title_str ' individual fits'],fig_type);

% STATS: n-way anova of slope values for each group
for ii = 1:2
    switch ii
        case 1 
            param_str = 'slp';
            b = 'slope';
        case 2
            param_str = 'intercept';
            b = 'intercept';
    end
    y_raw = [];
    for i = 1:num.exp
        exp = expOrder(i);
        y_raw = [y_raw; i*ones(num.trial(exp),1), stats(exp).(param_str)];
    end
    
    aov = anova(y_raw(:,1), y_raw(:,2));
    c = multcompare(aov,"CriticalValueType","bonferroni");
    c = table2array(c);
    names = strrep({grouped(expOrder).name},'_',' ');
    fig = getMultiCompSignTable(c, 1:num.exp, blkbgd, 0.05, names);
    title([title_str ' individual fit ' a ' ' b ' fits'],'color', foreColor)

    save_figure(fig,[fig_dir title_str ' individual fit ' a ' ' b ' fit stats'],fig_type);
end



%%

%% FIGURE: Time Course for single parameter -- select your metric
% TODO (2/26) update this to be able to plot single heating and cooling lines for all the trials 
% clearvars('-except',initial_vars{:})
% 
% plot_err = true;
% autoLim = true;
% xlim_auto = true; % change the time range for the x axis
% time_limits = [0,900]; % time limits if manual control over x-axis range
% nMax =  num.exp;%
% [~,backColor] = formattingColors(blkbgd); %get background colors
% 
% % Select the type of information to plot: 
% [title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end] = PlotParamSelection(true);
% switch questdlg('Plot error?','','True','False', 'Cancel','True')
%     case 'True'
%         plot_err = true;
%     case 'False'
%         plot_err = false;
%     case 'Cancel'
%         return
%     case ''
%         return
% end
% if isempty(title_str)
%     return
% end
% fig_dir = [saveDir, dir_end];
% % set figure folder
% if ~exist(fig_dir, 'dir')
%     mkdir(fig_dir)
% end
% 
% % set up figure aligments
% r = 5; %rows
% c = 3; %columns
% sb(1).idx = [1,2]; %temp timecourse
% sb(2).idx = [4,5,7,8,10,11,13,14]; % dependent var timecourse
% sb(3).idx = 3:c:r*c; %dependent var temp tuning curve
% 
% LW = 0.75;
% sSpan = 180;
% dataString = cell([1,num.exp]);
% 
% % FIGURE:
% fig = getfig('',true);
% for i = num.exp:-1:1
%     x = grouped(i).time;
%     kolor = grouped(i).color;
% 
%     %temp
%     subplot(r,c,sb(1).idx); hold on
%         y = grouped(i).temp;
%         plot(x,y,'LineWidth',2,'Color',kolor)
% 
%    % selected parameter time course
%     subplot(r,c,sb(2).idx); hold on
%         switch dType
%             case 1 % single trial lines
%                 for trial = 1:num.trial(i)
%                     y = smooth(grouped(i).(pName).all(:,trial),sSpan, 'moving')*scaler;
%                     plot(x,y,'LineWidth',LW,'Color',kolor)
%                 end
%             case {2, 3} % avg line
%                 y = smooth(grouped(i).(pName).avg,sSpan, 'moving')*scaler;
%                 plot(x,y,'LineWidth',LW,'Color',kolor)
%         end
% 
%     %temp vs dependent variable tuning curve
%     subplot(r,c,sb(3).idx); hold on
% 
%      switch dType
%          case 1 % single trial lines
%             for trial = 1:num.trial(i)
%                 if strcmp(pName, 'dist')
%                     x = grouped(i).(pName).distavgbytemp(:,1);
%                     rawY = [grouped(i).increasing.all(:,trial),grouped(i).decreasing.all(:,trial)];
%                 else
%                     x = grouped(i).(pName).temps;
%                     rawY = [grouped(i).(pName).increasing.raw(:,trial),grouped(i).(pName).decreasing.raw(:,trial)];
%                 end
%                 y = mean(rawY,2,'omitnan')*scaler;
%                 plot(x,y,'color',kolor,'linewidth',1.25)
%             end
% 
%          case 2 % avg lines (combined heating and cooling)
%             if strcmp(pName, 'dist')
%                 x = grouped(i).(pName).distavgbytemp(:,1);
%                 rawY = [grouped(i).increasing.all,grouped(i).decreasing.all];
%             else
%                 x = grouped(i).(pName).temps;
%                 rawY = [grouped(i).(pName).increasing.raw,grouped(i).(pName).decreasing.raw];
%             end
%             y = mean(rawY,2,'omitnan')*scaler;
%             y_err = (std(rawY,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
%             plot(x,y,'color',kolor,'linewidth',1.25)
%             plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
% 
%          case 3 % separated heating and cooling
%             if strcmp(pName, 'dist')
%                 x = grouped(i).(pName).distavgbytemp(:,1);
%                 YC = grouped(i).decreasing.all;
%                 YH = grouped(i).increasing.all;
%             else
%                 x = grouped(i).(pName).temps;
%                 YC = grouped(i).(pName).decreasing.raw;
%                 YH = grouped(i).(pName).increasing.raw;
%             end
%             % cooling
%             y = mean(YC,2,'omitnan')*scaler;
%             y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
%             plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
%             plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '--')
%             % heating
%             y = mean(YH,2,'omitnan')*scaler;
%             y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
%             plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
%             plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '-')          
%      end
%      dataString{i} = grouped(i).name;
% end
% 
% % FORMATING AND LABELS
% formatFig(fig,blkbgd,[r,c],sb);
% % temp
% subplot(r,c,sb(1).idx)
% ylabel('\circC')
% set(gca,"XColor",backColor)
% 
% % distance
% subplot(r,c,sb(2).idx)
% ylabel(y_lab)
% xlabel('time (min)')
% set(gca,'ydir',y_dir)
% % temp-distance relationship
% subplot(r,c,sb(3).idx)
% ylabel(y_lab)
% xlabel('temp (\circC)')
% % if ~autoLim
% %     ylim(dt_lim)
% % end
% h_line(nullD,'grey',':',2) %36.2
% set(gca,'ydir',y_dir)
% 
% if ~xlim_auto
%     subplot(r,c,sb(1).idx)
%     set(gca, 'xlim', time_limits)
%     subplot(r,c,sb(2).idx)
%     set(gca, 'xlim', time_limits)
%     % ylim([0,100])
%     % ylim([20,100])
% end
% 
% % legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)
% 
% % save figure
% save_figure(fig,[fig_dir 'Timecourse summary ' title_str],fig_type);

%% 3D visualization of parameters in both a temp and time domain
clearvars('-except',initial_vars{:})
fig_dir = createFolder([saveDir, 'Stats/']);
[foreColor,~] = formattingColors(blkbgd); %get background colors

% [title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end] = PlotParamSelection(false);


%% TODO: update this to work for all other trial types: currently plots the avg occupancy for the coldest two degrees of the fast temp ramp
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);

% SHOW THE FLIES THAT ARE IN THE FOOD QUADRANT DURING 'FICTIVE' COLDEST 2C (FOR F LRR)
buff = 0.2;
fig = getfig('',1,[508 680]);
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    % find time points: 
    tp = getTempTurnPoints(data(i).temp_protocol);
    offset = fps*60*abs(tp.rates(1));
    roi1 = tp.down(:,2);
    roi2 = (roi1-720);
    roi = [];
    for ii = 1:length(roi1)
        roi = [roi; (roi2(ii):roi1(ii))'];
    end
    % plot time points: 
     y = mean(grouped(i).quadrant.all(roi,:),1,'omitnan');
     x = shuffle_data(linspace(exp-buff, exp+buff, num.trial(i)));
     scatter(x,y,50,grouped(i).color,'filled')
     plot([exp-0.35, exp+0.35],[mean(y), mean(y)],'color', foreColor,'linewidth', 2)
end
formatFig(fig, blkbgd);
h_line(25,'grey', '--', 1.5)
set(gca, 'xcolor', 'none')
ylim([0,100])
ylabel('Food quadrant occupancy (%)')
save_figure(fig, [saveDir, 'quad occ during coldest 2C cooling no food'],fig_type);


% SHOW THE FLIES THAT ARE IN THE FOOD QUADRANT DURING 'FICTIVE' COLDEST 2C (FOR F LRR)
buff = 0.2;
fig = getfig('',1,[508 680]);
hold on
for exp = 1:num.exp
    i = expOrder(exp);
    % find time points: 
    tp = getTempTurnPoints(data(i).temp_protocol);
    offset = fps*60*abs(tp.rates(1));
    roi1 = tp.down(:,2);
    roi2 = (roi1-720);
    roi = [];
    for ii = 1:length(roi1)
        roi = [roi; (roi2(ii):roi1(ii))'];
    end
    % plot time points: 
     y = mean(grouped(i).ring.all(roi,:),1,'omitnan');
     x = shuffle_data(linspace(exp-buff, exp+buff, num.trial(i)));
     scatter(x,y,50,grouped(i).color,'filled')
     plot([exp-0.35, exp+0.35],[mean(y), mean(y)],'color', foreColor,'linewidth', 2)
end
formatFig(fig, blkbgd);
h_line(25,'grey', '--', 1.5)
set(gca, 'xcolor', 'none')
ylim([0,100])
ylabel('Outer ring occupancy (%)')
save_figure(fig, [saveDir, 'outer ring occ during coldest 2C cooling no food'],fig_type);


%% TODO = this works for the temp hold trials but  needs updating or saving & dynamics

[foreColor,~] = formattingColors(blkbgd);

tp = getTempTurnPoints('linear_ramp_F_25-17');
tRange = 17:19;
roi1 = tp.down(:,2);
roi2 = (roi1-720);
roi = [];
for i = 1:length(roi1)
    roi = [roi; (roi2(i):roi1(i))'];
end

% SHOW THE FLIES THAT ARE IN THE FOOD QUADRANT DURING 'FICTIVE' COLDEST 2C (FOR F LRR)
cList = repmat({'darkorange','grey'},[1,num.exp/2]);
buff = 0.2;
fig = getfig('',1,[508 680]);
hold on
for i = 2:2:num.exp
     y = mean(grouped(i).quadrant.all(roi,:),1,'omitnan');
     x = shuffle_data(linspace(i-buff, i+buff, num.trial(i)));
     scatter(x,y,50,Color(cList{i}),'filled')
     plot([i-0.35, i+0.35],[mean(y), mean(y)],'color', foreColor,'linewidth', 2)
end
formatFig(fig, blkbgd);
h_line(25,'grey', '--', 1.5)
set(gca, 'xcolor', 'none')
ylim([0,100])
ylabel('Food quadrant occupancy (%)')
save_figure(fig, [saveDir, 'quad occ during coldest 2C cooling no food'],fig_type,1,0);
for i = 1:2:num.exp
     y = mean(grouped(i).quadrant.all(roi,:),1,'omitnan');
     x = shuffle_data(linspace(i-buff, i+buff, num.trial(i)));
     scatter(x,y,50,Color(cList{i}),'filled')
     plot([i-0.35, i+0.35],[mean(y), mean(y)],'color', foreColor,'linewidth', 2)
end
save_figure(fig, [saveDir, 'quad occ during coldest 2C cooling food'],fig_type,1,1);

% SHOW THE FLIES THAT ARE IN THE FOOD QUADRANT DURING 'FICTIVE' COLDEST 2C (FOR F LRR)
cList = repmat({'darkorange','grey'},[1,num.exp/2]);
buff = 0.2;
buff2 = 0.4;
fig = getfig('',1,[508 680]);
hold on
for i = 2:2:num.exp
     y = mean(grouped(i).ring.all(roi,:),1,'omitnan');
     x = shuffle_data(linspace(i-buff, i+buff, num.trial(i)));
     scatter(x,y,50,Color(cList{i}),'filled')
     plot([i-buff2, i+buff2],[mean(y), mean(y)],'color', foreColor,'linewidth', 2)
end
ylim([0,100])
formatFig(fig, blkbgd);
h_line(25,'grey', '--', 1.5)
set(gca, 'xcolor', 'none')
ylabel('Outer ring occupancy (%)')
save_figure(fig, [saveDir, 'Outer ring occ during coldest 2C cooling no food'],fig_type,1,0);
for i = 1:2:num.exp
     y = mean(grouped(i).ring.all(roi,:),1,'omitnan');
     x = shuffle_data(linspace(i-buff, i+buff, num.trial(i)));
     scatter(x,y,50,Color(cList{i}),'filled')
     plot([i-buff2, i+buff2],[mean(y), mean(y)],'color', foreColor,'linewidth', 2)
end
save_figure(fig, [saveDir, 'Outer ring occ during coldest 2C cooling with food'],fig_type,1,1);
























