

% SPEED DATA COMPARISON


%% FIGURE: speed histogram
clearvars('-except',initial_vars{:})
autoSave = true;
[foreColor,~] = formattingColors(blkbgd);

% Simple histogram of speed across full experiment
bins = 0:0.5:15;
% xlimit = [-0.5,7];
xlimit = [-.5,20];

[nrows, ncols] = subplot_numbers(num.exp, 3);
if num.exp==3
    nrows = 1;
    ncols = 3;
end
if nrows==1
    figPos = [953 439];
else
    figPos = [953 1021];
end

fig = getfig('',true,figPos);
for ii = 1:num.exp
    i = expOrder(ii);
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,true);
    end
    % histogram for each group
    subplot(nrows,ncols,ii)
    histogram(temp,bins,'FaceColor',grouped(i).color,'facealpha',1,'EdgeColor',foreColor)
    v_line(mean(temp,'omitnan'),foreColor,'--',1.5)
%     title(grouped(i).name)
    xlabel('speed (mm/s)')
    ylabel('count')
    xlim(xlimit)
end
formatFig(fig,blkbgd,[nrows,ncols]);
save_figure(fig,[saveDir expGroup ' speed histograms'],fig_type,autoSave);

%% FIGURE AND STATS: average experiment speed
% Comparision of avg speed and range across full experiments
clearvars('-except',initial_vars{:})
autoSave = false;
plot_err = true;
[~,backColor] = formattingColors(blkbgd);
buff = 0.2;
SZ = 50;
LW = 1.5;
r = 1;
c = 3;
sb(1).idx = 1:2; %avg speed per temp
sb(2).idx = 3;   %avg speed


fig = getfig('',true);
subplot(r,c,sb(1).idx); hold on
    %pull temp avg data
    [speed,speed_err,temp] = deal([]);
    for i = 1:num.exp
        temp(:,i) = data(i).G(1).TR.temps;
        tempIdx = data(i).G(1).TR.tempIdx;
        for t = 1:size(temp,1)-1
            loc = (tempIdx==t);
            raw = mean(grouped(i).speed.all(loc,:),'omitnan');
            y = mean(raw);
            y_err = std(raw)./sqrt(num.trial(i));
            speed(t,i) = y;
            speed_err(t,i) = y_err;
        end
    end
    for i = 1:num.exp
        kolor = grouped(i).color;
        loc = ~isnan(speed(:,i));
        x = temp(loc,i);
        y = speed(loc,i);
        y_err = speed_err(loc,i);

       plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
        plot(x,y,'color',kolor,'linewidth', LW)
    end
    xlabel('Temperature (\circC)')
    ylabel('Speed (mm/s)')

subplot(r,c,sb(2).idx);hold on
for ii = 1:num.exp
    i = expOrder(ii);
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,false);
    end
    y = mean(temp,1,'omitnan');
    y_avg = mean(mean(temp,'omitnan'));
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
    scatter(x,y,SZ,grouped(i).color,'filled')
    plot([ii-buff,ii+buff],[y_avg,y_avg],'color', grouped(i).color,'linewidth',LW)
end
xlim([0,num.exp+(buff*2)])
ylabel('avg speed (mm/s)')
formatFig(fig,blkbgd,[r,c],sb);
set(gca,'xcolor',backColor)
save_figure(fig,[saveDir expGroup ' avg speed'],fig_type,autoSave);

% AVG SPEED STATS:
[datastats.all, datastats.id] = deal([]);
for i = 1:num.exp
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,false);
    end
    y = mean(temp,1,'omitnan');
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
fprintf('\n\nAverage speed statistics\n\n')
[Group1,Group2,P_Value] = deal([]);
idx = 0;
if ~any(significantHypotheses)
    disp('No statistical differences in avg speed between groups')
else
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

%% FIGURE: cumulative distribution function of speed
clearvars('-except',initial_vars{:})
autoSave = true;
autoLim = false;
xlimit = [0,15];

LW = 1.5;

% xlimit = [0,20];

fig = figure; hold on
for i = 1:num.exp
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,true);
    end
    h = cdfplot(temp);
    set(h,'color',grouped(i).color,'linewidth',LW)
end
if ~autoLim
    xlim(xlimit)
end
xlabel('Speed (mm/s)')
ylabel('Empirical CDF')
title('')
formatFig(fig, blkbgd);
set(gca,'xgrid','off','ygrid','off')
set(gca,'ytick',0:0.2:1)
save_figure(fig,[saveDir expGroup ' speed CDF'],fig_type,autoSave);

%% FIGURE & STATS: speed hysteresis across groups
clearvars('-except',initial_vars{:})
plotSig = true;
plot_err = false;
LW = 1;
buff = 0.2;
alpha = 0.05;
[foreColor,backColor] = formattingColors(blkbgd);

SZ = 50;
r = 1; %rows
c = 3; %columns

[groupName, stats] = deal([]);

fig = getfig('',true);
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    x = data(i).G(1).TR.temps(1:end-1);
    [h_speed,c_speed] = deal([]);
    rateIdx = data(i).G(1).TR.rateIdx;

    % heating locations
    idx = find(data(i).G(1).TR.rates>0);
    heatloc = rateIdx==idx;
    h_style = '-';
    % cooling locations
    idx = find(data(i).G(1).TR.rates<0);
    coolloc = rateIdx==idx;
    c_style = '--';

    %pull temp avg data
    tempIdx = data(i).G(1).TR.tempIdx;
    speed_data = grouped(i).speed.all;
    for t = 1:length(x)
        temp_loc = (tempIdx==t);
        % heating
        loc = temp_loc & heatloc;
        raw = mean(speed_data(loc,:),'omitnan');
        h_speed(t,:) =  raw;
        % cooling
        loc = temp_loc & coolloc;
        raw = mean(speed_data(loc,:),'omitnan');
        c_speed(t,:) =  raw;
    end

    % plot heat and cooling lines
    subplot(r,c,1); hold on
        %heating
        heat_y = mean(h_speed,2,'omitnan');
        heat_err = std(h_speed, 0,2, 'omitnan');
        plot_error_fills(plot_err, x, heat_y, heat_err, kolor,  fig_type, 0.35);
        plot(x,heat_y,'color',kolor,'linewidth',LW,'linestyle',h_style)
        % cooling
        cool_y = mean(c_speed,2,'omitnan');
        cool_err = std(c_speed, 0,2, 'omitnan');
        plot_error_fills(plot_err, x, cool_y, cool_err, kolor,  fig_type, 0.35);
        plot(x,cool_y,'color',kolor,'linewidth',LW,'linestyle',c_style)

    % plot speed hysteresis for each temp bin
    subplot(r,c,2); hold on
        y = cool_y-heat_y;
        plot(x,y,'color',kolor,'linewidth',LW)

    % plot cumulative hysteresis for each trial
    subplot(r,c,3); hold on
        y = sum((c_speed-h_speed),1,'omitnan');
        y_avg = mean(y);
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[y_avg,y_avg],'color',kolor,'linewidth',1.3)
        % save cumulative difference for stats
        stats = autoCat(stats,y,true);
        groupName{ii} = grouped(i).name;
end

%formatting and labels
formatFig(fig,blkbgd,[r,c]);
subplot(r,c,1)
    xlabel('Temp (\circC)')
    ylabel('Speed (mm/s)')
subplot(r,c,2)
    xlabel('Temp (\circC)')
    ylabel('Speed diff (cool-heat)')
    h_line(0,foreColor,':',1)
subplot(r,c,3)
    xlim([0.5,num.exp+.5])
    h_line(0,foreColor,':',1)
    ylabel('Cumulative speed difference')
    set(gca,'xcolor',backColor)


% SIG HYSTERESIS STATS: are the means of any groups different from zero?
p = [];
for ii = 1:num.exp
    [~,p(ii)] = ttest(stats(ii,:));
end
%Bonferonni correction:
m = num.exp;
p_limit = alpha/m;
h = p<=p_limit;
stats_tbl = table(groupName',h',p','VariableNames',{'group','significant','p value'});
disp(stats_tbl)

% add significance stars to the figure:
if plotSig
    y_pos = rangeLine(fig,10,false);
    subplot(r,c,3); hold on
    for ii = 1:num.exp
        if h(ii)
            scatter(ii,y_pos,100,foreColor,'*')
        end
    end
end

% SIG DIFF BETWEEN GROUPS STATS:
%format for anova
stats_id = repmat((1:num.exp)',[1,size(stats,2)]);
stats_id = reshape(stats_id,[numel(stats),1]);
x = reshape(stats,[numel(stats),1]);
loc = isnan(x);
stats_id(loc) = [];
x(loc) = [];

% determine which groups differ from each other
[~,~,data_stats] = anova1(x,stats_id,'off');
[c,~,~,~] = multcompare(data_stats,"Display",'off');

% bonferonni multiple comparisons correction
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nSpeed hysteresis comparison statistics:\n\n')
[Group1,Group2,P_Value] = deal([]);
idx = 0;
if any(significantHypotheses)
    for i = 1:length(significantHypotheses)
        if significantHypotheses(i)
            idx = idx+1;
            Group1{idx,1} = groupName{c(i,1)};
            Group2{idx,1} = groupName{c(i,2)};
            P_Value(idx,1) = c(i,6);
        end
    end
    sig_comp = table(Group1,Group2,P_Value);
    disp(sig_comp)
else
    disp('no significant differences between groups')
end


% save figure
save_figure(fig,[saveDir expGroup ' speed hysteresis summary'],fig_type);

%% FIGURE: Speed - distance to food correlation
clearvars('-except',initial_vars{:})
plot_err = true;
corr_coef = [];
num_points = 300;
buff = 0.2;
speedLim = [0,15];
SZ = 50;
LW = 1.5;
[foreColor,backColor] = formattingColors(blkbgd);
r = 1;
c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;


% get correlation data
for i = 1:num.exp
    % get speed / distance information
    speed = grouped(i).speed.all;   %speed for each trial within the exp.
    fly_dist = grouped(i).dist.all; %dist for each trial within the exp

    for trial = 1:num.trial(i)
        temp = [];
        temp(:,1) = speed(:,trial);
        temp(:,2) = fly_dist(:,trial);
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-distance correlation
        rho = corr(temp);
        corr_coef(i).all(trial) = rho(1,2);
    end
 end

% PLOT FIGURE
fig = getfig('',true);
% linear fit of speed-distance to food|source
subplot(r,c,sb(1).idx); hold on
    for i = 1:num.exp

        % get speed / distance information
        speed = grouped(i).speed.all;   %speed for each trial within the exp.
        fly_dist = grouped(i).dist.all; %dist for each trial within the exp
        % reshape & remove NaNs
        y = [];
        y(:,1) = speed(:);
        y(:,2) = fly_dist(:);
        loc = any(isnan(y),2);
        y(loc,:) = [];

        % speed-distance correlation
        rho = corr(y);
        corr_coef(i).group = rho(1,2);

        % sort by speed
        sorted_data = sortrows(y,1);
        idx = randi(length(y),[num_points,1]);
        X = sorted_data(idx,1);
        Y = sorted_data(idx,2);

        % Plot random selection of points & linear regression line of best fit
        kolor = grouped(i).color;
        scatter(X,Y,10,kolor)
        fitX = sorted_data(:,1);
        [p,S] = polyfit(fitX,sorted_data(:,2),1);
        [y_fit,delta] = polyval(p,fitX,S);
        plot(fitX,y_fit,'color',kolor,'linestyle','-','linewidth',1)
        if plot_err
            plot(fitX,y_fit+2*delta,'color',kolor,'linestyle','--','linewidth',1)
            plot(fitX,y_fit-2*delta,'color',kolor,'linestyle','--','linewidth',1)
        end
    end
    xlim(speedLim)
    xlabel('Speed (mm/s)')
    ylabel('Distance from source (mm)')

% correlation coefficients
subplot(r,c,sb(2).idx); hold on
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
     ylabel('correlation coefficient')
     h_line(0,foreColor,':',1)
formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(2).idx)
set(gca,'xcolor',backColor)

% save figure
save_figure(fig,[saveDir expGroup ' speed distance correlation'],fig_type);

%% TODO : FIGURE: Speed - eccentricity correlation
clearvars('-except',initial_vars{:})
plot_err = true;
corr_coef = [];
num_points = 300;
buff = 0.2;
speedLim = [0,15];
SZ = 50;
LW = 1.5;
[foreColor,backColor] = formattingColors(blkbgd);
r = 1;
c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;

% get correlation data
for i = 1:num.exp
    % get speed / distance information
    speed = grouped(i).speed.all;   %speed for each trial within the exp.
    fly_dist = grouped(i).ecent.all; %dist for each trial within the exp

    for trial = 1:num.trial(i)
        temp = [];
        temp(:,1) = speed(:,trial);
        temp(:,2) = fly_dist(:,trial);
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-distance correlation
        rho = corr(temp);
        corr_coef(i).all(trial) = rho(1,2);
    end
 end

% PLOT FIGURE
fig = getfig('',true);
% linear fit of speed-distance to food|source
subplot(r,c,sb(1).idx); hold on
    for i = 1:num.exp

        % get speed / distance information
        speed = grouped(i).speed.all;   %speed for each trial within the exp.
        fly_dist = grouped(i).ecent.all; %dist for each trial within the exp
        % reshape & remove NaNs
        y = [];
        y(:,1) = speed(:);
        y(:,2) = fly_dist(:);
        loc = any(isnan(y),2);
        y(loc,:) = [];

        % speed-distance correlation
        rho = corr(y);
        corr_coef(i).group = rho(1,2);

        % sort by speed
        sorted_data = sortrows(y,1);
        idx = randi(length(y),[num_points,1]);
        X = sorted_data(idx,1);
        Y = sorted_data(idx,2);

        % Plot random selection of points & linear regression line of best fit
        kolor = grouped(i).color;
        scatter(X,Y,10,kolor)
        fitX = sorted_data(:,1);
        [p,S] = polyfit(fitX,sorted_data(:,2),1);
        [y_fit,delta] = polyval(p,fitX,S);
        plot(fitX,y_fit,'color',kolor,'linestyle','-','linewidth',1)
        if plot_err
            plot(fitX,y_fit+2*delta,'color',kolor,'linestyle','--','linewidth',1)
            plot(fitX,y_fit-2*delta,'color',kolor,'linestyle','--','linewidth',1)
        end
    end
    xlim(speedLim)
    xlabel('Speed (mm/s)')
    ylabel('Distance from center (mm)')

% correlation coefficients
subplot(r,c,sb(2).idx); hold on
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
     ylabel('correlation coefficient')
     h_line(0,foreColor,':',1)
formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(2).idx)
set(gca,'xcolor',backColor)

% save figure
save_figure(fig,[saveDir expGroup ' speed eccentricity correlation'],fig_type);

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

% % Perform the Mann-Whitney U test (Wilcoxon rank-sum test)
% [p, h, stats] = ranksum(plotData(1).rho, plotData(3).rho); % caviar (1) vs waxed (3)
% % Display the results
% fprintf('p-value: %.4f\n', p);
% if h == 0
%     fprintf('The null hypothesis cannot be rejected: The two groups are not significantly different.\n');
% else
%     fprintf('The null hypothesis is rejected: The two groups are significantly different.\n');
% end


% save figure
save([saveDir expGroup ' temp speed correlation ramps only'],'plotData');
save_figure(fig,[saveDir expGroup ' temp speed correlation ramps only'],'-png',true,false);
save_figure(fig,[saveDir expGroup ' temp speed correlation ramps only'],'-pdf',true,true);

%% FIGURE: Speed - outter ring eccentricity correlation [for FULL experiment]
clearvars('-except',initial_vars{:})
plot_err = true;
corr_coef = [];
num_points = 300;
buff = 0.2;
speedLim = [0,15];
SZ = 50;
LW = 1.5;
[foreColor,backColor] = formattingColors(blkbgd);
r = 1;
c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;

% get correlation data
for i = 1:num.exp
    % get speed / distance information
    speed = grouped(i).speed.all;   %speed for each trial within the exp.
    fly_dist = grouped(i).ring.all; %dist for each trial within the exp

    for trial = 1:num.trial(i)
        temp = [];
        temp(:,1) = speed(:,trial);
        temp(:,2) = fly_dist(:,trial);
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-occupancy correlation
        rho = corr(temp);
        corr_coef(i).all(trial) = rho(1,2);
    end
 end

% PLOT FIGURE
fig = getfig('',true);
% linear fit of speed-distance to food|source
subplot(r,c,sb(1).idx); hold on
    for i = 1:num.exp

        % get speed / distance information
        speed = grouped(i).speed.all;   %speed for each trial within the exp.
        fly_dist = grouped(i).ring.all; %dist for each trial within the exp
        % reshape & remove NaNs
        y = [];
        y(:,1) = speed(:);
        y(:,2) = fly_dist(:);
        loc = any(isnan(y),2);
        y(loc,:) = [];

        % speed-distance correlation
        rho = corr(y);
        corr_coef(i).group = rho(1,2);

        % sort by speed
        sorted_data = sortrows(y,1);
        idx = randi(length(y),[num_points,1]);
        X = sorted_data(idx,1);
        Y = sorted_data(idx,2);

        % Plot random selection of points & linear regression line of best fit
        kolor = grouped(i).color;
        scatter(X,Y,10,kolor)
        fitX = sorted_data(:,1);
        [p,S] = polyfit(fitX,sorted_data(:,2),1);
        [y_fit,delta] = polyval(p,fitX,S);
        plot(fitX,y_fit,'color',kolor,'linestyle','-','linewidth',1)
        if plot_err
            plot(fitX,y_fit+2*delta,'color',kolor,'linestyle','--','linewidth',1)
            plot(fitX,y_fit-2*delta,'color',kolor,'linestyle','--','linewidth',1)
        end
    end
    xlim(speedLim)
    xlabel('Speed (mm/s)')
    ylabel('Ring occupancy (%)')

% correlation coefficients
subplot(r,c,sb(2).idx); hold on
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
     ylabel('correlation coefficient')
     h_line(0,foreColor,':',1)
formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(2).idx)
set(gca,'xcolor',backColor)

% save figure
save_figure(fig,[saveDir expGroup ' speed ring occupancy correlation'],fig_type);


%% FIGURE: Speed - outter ring eccentricity correlation [for RAMPS only]
clearvars('-except',initial_vars{:})
plot_err = true;
corr_coef = [];
num_points = 300;
buff = 0.2;
speedLim = [0,15];
SZ = 50;
LW = 1.5;
[foreColor,backColor] = formattingColors(blkbgd);
r = 1;
c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;

% Select the temperature protocol for ramp comparison (if not all the protocols are the same)
protocolList = unique({data(:).temp_protocol});
if length(protocolList)>1
    idx = listdlg('PromptString', 'Select the temperature protocol for time ROI selection', 'ListString',['Respective protocols'; protocolList(:)]);
    if idx>1
        t_proto = cell(num.exp,1);
        t_proto(:) = protocolList(idx-1);
    else
        t_proto = {data(:).temp_protocol};
    end
else
     t_proto = {data(:).temp_protocol};
end
        
% get correlation data
for i = 1:num.exp
    tp = getTempTurnPoints(t_proto{i});
    ROI = [tp.UpROI, tp.DownROI];

    % get speed / distance information
    speed = grouped(i).speed.all(ROI,:);     %speed for each trial within the exp.
    fly_dist = grouped(i).ring.all(ROI,:);         %dist for each trial within the exp

    for trial = 1:num.trial(i)
        temp = [];
        temp(:,1) = speed(:,trial);
        temp(:,2) = fly_dist(:,trial);
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-occupancy correlation
        rho = corr(temp);
        corr_coef(i).all(trial) = rho(1,2);
    end
 end

pixWidth = 60; % additional figure pixel size for image for each extra experiment group
figSize = [pixWidth + (pixWidth*num.exp),590];

% PLOT FIGURE correlation coefficients
fig = getfig('',true);
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
           plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',LW,'HandleVisibility','off')
     end
     xlim([0.5,num.exp+.5])
     ylabel('speed-ecentricity correlation')
     h_line(0,foreColor,':',1)
formatFig(fig,blkbgd);
legend({grouped(:).name},'fontsize', 8,'box', 'off')

% save figure
save_figure(fig,[saveDir expGroup ' speed-eccentricity ramp only correlation'],fig_type);


