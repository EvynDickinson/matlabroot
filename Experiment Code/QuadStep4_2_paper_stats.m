


%% FIGURE: Thermal Threat / Temperature Escape Behavior
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
        coolIdx = find(grouped(i).position.temp_rates<0); %warming rate index
        warmIdx = find(grouped(i).position.temp_rates>0); %cooling rate index
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
    cooling
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
set(gca,'XTick',0:1:num.exp,'YTick',0:1:num.exp,'XTickLabel',names,'YTickLabel',names)
set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
save_figure(fig,[fig_dir 'outer ring occupancy cooling slope stats'],fig_type);


% How do the COOLING temperature dependent lines differ? Intercept based analysis
temp = []; occ = []; exp = [];
%format the data for the statistical test
for i = 1:num.exp
    x = repmat((grouped(i).ring.temps'),[num.trial(i),1]);
    temp = [temp; x];
    y = grouped(i).ring.decreasing.raw(:);
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
set(gca,'XTick',0:1:num.exp,'YTick',0:1:num.exp,'XTickLabel',names,'YTickLabel',names)
set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
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
set(gca,'XTick',0:1:num.exp,'YTick',0:1:num.exp,'XTickLabel',names,'YTickLabel',names)
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
set(gca,'XTick',0:1:num.exp,'YTick',0:1:num.exp,'XTickLabel',names,'YTickLabel',names)
set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
save_figure(fig,[fig_dir 'outer ring occupancy warming intercept stats'],fig_type);

%% 

% test if the data comes from a normal distribution:
h = kstest(x)


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
scaler = 1; plot_err = 1;
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
ylabel('flies in outer ring (%)')
xlabel('temperature (\circC)')

formatFig(fig, blkbgd,[r,c],sb);
subplot(r,c,sb(1).idx)
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
[foreColor,backColor] = formattingColors(blkbgd);

fig = getfig('',true,[480 680]); 
hold on
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    y = abs(grouped(i).ring.decreasing.raw-grouped(i).ring.increasing.raw);
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
    y = abs(grouped(i).ring.decreasing.raw-grouped(i).ring.increasing.raw);
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
set(gca,"XColor",backColor)

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



