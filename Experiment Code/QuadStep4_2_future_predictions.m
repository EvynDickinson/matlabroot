
% load existing model data: 

load([saveDir 'temp occupation correlation ramps only.mat'])


%% Predictions of CURRENT behavior based on PAST conditions
clearvars('-except',initial_vars{:})
models = struct;
figDir = [saveDir, 'models/'];
if ~exist(figDir,'dir')
    mkdir(figDir)
end

trial = 1;

exp = 1;

raw_dist = grouped(exp).dist.all; % distance to food for all trials in the group
temp = grouped(exp).temp;
temprate = nan(size(temp));
y = data(exp).G(trial).TR.data(:,4); % temp rate over time for the experiment
temprate(1:length(y)) = y;

% set up the range of testing points in the past to measure the 'predicted temp'
delta_t = [0, 15]; % range of times into the past to test (min)
delta_t = delta_t*60; % convert to seconds from minutes
skip_size = 15; % frames to skip between test times
delta_T = (delta_t(1)*fps):skip_size:delta_t(2)*fps;
if delta_T(1)==0 % get rid of a zero delay?
    delta_T(1) = [];
end

% set up the range of temp rate integration regions
delta_r = [0,2]; % period of temp rate integration that the fly might have used (min)
delta_r = delta_r*60; % convert to seconds from minutes
skip_size = 10; % seconds
delta_R = (delta_r(1)*fps):skip_size:delta_r(2)*fps;
if delta_R(1)==0 % get rid of a zero delay?
    delta_R(1) = [];
end

nRates = length(delta_R);
nTemps = length(delta_T);

% ********************************************************************
% ************** !!!!! CAUTION EXTREMELY SLOW !!!!! **************
% what was the temp and temp rate at dt?
for jj = 1:nRates
    tic
    dr = delta_R(jj); % set this as the current temprate integration time region

    for ii = 1:nTemps
        dt = delta_T(ii); % what is the time point in the past from which the temp might be predicted?

        MT  = nan(size(temp)); % allocate an empty data vector for the temp information
        start_offset = dt + dr; % points at the start to omit because we don't have data in the past for them
        Ti = (1:length(MT)-start_offset)'; % 'start time points' index number

        temp_past = temp(Ti); % points for the temp offset in the past
        temprate_past = [];
        for i = 1:length(Ti)
            temprate_past(i,1) = mean(temprate(Ti(i):Ti(i)+dr),'omitnan');
        end
        temp_prediction = MT;
        change_in_temp = temprate_past*(dt/60/fps); % convert from deg/min to deg
        temp_prediction(start_offset+1:end) =  temp(Ti) + change_in_temp;

        % save temperature prediction information
        models(ii,jj).dt = dt;
        models(ii,jj).dr = dr;
        models(ii,jj).predicted_temp = temp_prediction;
    end
    toc
    disp([num2str(jj) '/' num2str(length(delta_R)) ' temp rate integration time'])
end
% ********************************************************************
% ********************************************************************
save([saveDir 'models'],'models','-v7.3')
initial_vars{end+1} = 'models';

% Plot the predicted temperatures:

% test a small integration window: 
rr = 1;
kolor = Color('white', 'blue', nTemps);
time = grouped(exp).time;

fig = getfig('',1);
hold on
for tt = 1:nTemps
    plot(time, models(tt,rr).predicted_temp,'Color', kolor(tt,:))
end
plot(time, grouped(exp).temp,'color', 'k','LineWidth',1)
xlabel('time (min)')
ylabel('temp (\circC)')
title(['model predicted temperatures | rate integration: ' num2str(delta_R(rr)) ' seconds'],'color', 'k')
formatFig(fig);
save_figure(fig, [figDir grouped(exp).name ' predicted temp lines'], '-png');

% test a small integration window: 
kolor = Color('white', 'blue', nTemps);
time = grouped(exp).time;
r = 6; 
c = 6;
fig = getfig('',1);
for rr = 1:nRates
    subplot(r,c,rr)
    hold on
    for tt = 1:5:nTemps
        plot(time, models(tt,rr).predicted_temp,'Color', kolor(tt,:))
    end
    plot(time, grouped(exp).temp,'color', 'k','LineWidth',1)
    % xlabel('time (min)')
 
    title([num2str(delta_R(rr)) ' s'],'color', 'k')
end
% add labels
leftedge = 1:r:r*c;
bottomrow = (r*c)-c+1:r*c;
for rr = leftedge
    subplot(r,c,rr)
    ylabel('temp (\circC)')
end
for rr = bottomrow
    subplot(r,c,rr)
    xlabel('time (min)')
end
for rr = 1:r*c
    subplot(r,c,rr)
    if ~any(rr==bottomrow)
        set(gca,'xcolor','none')
    end
    if ~any(rr==leftedge)
        set(gca, 'ycolor', 'none')
    end
    ylim([10,30])
end

save_figure(fig,[figDir grouped(trial).name ' model temperature predictions'],'-pdf',1,0);
save_figure(fig,[figDir grouped(trial).name ' model temperature predictions'],'-png');

%% ANALYSIS: correlation for measure of linearity

% Screen out the recovery period: 
% pre / post buffer: 16 mins
buff = 16;
buff = buff*fps;
tP = getTempTurnPoints(data(exp).temp_protocol);
roi = tP.down(:,1)-buff : tP.up(:,2)-buff; % selected region around temp ramps

% CAUTION: SLOW
correlations = nan([nRates,nTemps,num.trial(exp)]);
for rr = 1:nRates
    for tt = 1:nTemps
        A = models(tt,rr).predicted_temp(roi);
        for trial = 1:num.trial(exp)
            B = [A, grouped(exp).dist.all(roi,trial)];
            R = corrcoef(B,'Rows','pairwise');
            correlations(rr,tt,trial) = R(1,2);
        end
    end
    disp([num2str(rr) '/' num2str(nRates)])
end
% save([saveDir 'correlations'],'correlations','-v7.3')

% ground truth correlation: 
g_corr = nan([num.trial(exp),1]);
 for trial = 1:num.trial(exp)
    B = [grouped(exp).temp, grouped(exp).dist.all(:,trial)];
    R = corrcoef(B,'Rows','pairwise');
    g_corr(trial) = R(1,2);
end


% plot the correlations for each of the delays
sz = 10;
avg_sz = sz + 10;
kolor = Color('grey', 'blue', nTemps);
r = 6; 
c = 6;
fig = getfig('',1);
for rr = 1:nRates
    subplot(r,c,rr)
    hold on
    for trial = 1:num.trial(exp)
        C = Color('grey');
        x = delta_T;
        y = squeeze(correlations(rr,:,trial));
        plot(x,y,'color', C, 'linewidth',0.5)
    end
    % avg across trials
    y_avg = mean(correlations(rr,:,:),3);
    plot(x,y_avg,'color', 'k', 'linewidth',1)
    % plot the original data
    x = zeros([1,num.trial(exp)]);
    y = g_corr;
    scatter(x,y,sz,'r')
    scatter(x(1), mean(y),avg_sz, 'k','filled')
    h_line(mean(g_corr),'r','-')
    [~,idx] = min(y_avg);
    v_line(delta_T(idx),'green','-')
    title([num2str(delta_R(rr)) 's | min ' num2str(delta_T(idx))],'color', 'k')
end
% add labels
leftedge = 1:r:r*c;
bottomrow = (r*c)-c+1:r*c;
for rr = leftedge
    subplot(r,c,rr)
    ylabel('corr coeff')
end
for rr = bottomrow
    subplot(r,c,rr)
    xlabel('prediction delay (sec)')
end
for rr = 1:r*c
    subplot(r,c,rr)
    if ~any(rr==bottomrow)
        set(gca,'xcolor','none')
    end
    if ~any(rr==leftedge)
        set(gca, 'ycolor', 'none')
    end
    ylim([-1, 0])
    xlim([-50,3600])
end
save_figure(fig,[figDir grouped(exp).name ' model correlations for all times'],'-png');


%  Plot the strongest correlation for each type/trial...
kolor = Color('blue','yellow',nRates);
min_corr = [];
past_time = [];
for rr = 1:nRates
     y_avg = mean(correlations(rr,:,:),3);
     [~,idx] = min(y_avg);
     min_corr(rr) = y_avg(idx);
     past_time(rr) = delta_T(idx);
end
corr_diff =  mean(g_corr) - min_corr;
[~, max_idx] = max(corr_diff); % 'best' model prediction time
delay_t = past_time(max_idx);
integration_t = delta_R(max_idx);

fig = getfig('',1);
hold on
scatter(past_time,corr_diff,50,kolor,"filled")
h_line(0,'k',':')
% plot the best model
scatter(past_time(max_idx),corr_diff(max_idx),100,'r')
title(['Highest correlation: ' num2str(delay_t) 's prior with ' num2str(integration_t) 's temp rate integration window'])
xlabel('time in past for current temp prediction (s)')
ylabel('Difference in correlation from baseline')
formatFig(fig);

save_figure(fig,[figDir grouped(exp).name ' model performance'],'-png')

% plot the behavior of the best model vs the real data...
tempRange = 10:0.1:30;
nBins = length(tempRange);

% find best model numbers: 
T_idx = find(delay_t==delta_T);
R_idx = max_idx;
x = models(T_idx,R_idx).predicted_temp; % predicted temp
y = grouped(exp).dist.all; %distance to food
tempBins = discretize(x,tempRange); % assign each temp point to a 0.5degC temp bin
%predicted version
[y_all,  y_err] = deal(nan(size(tempRange)));
for i = 1:nBins
    idx  = find(tempBins==i);
    if ~isempty(idx)
        y_all(i) = mean(mean(y(idx,:),'omitnan'),'omitnan');
        y_err(i) = std(mean(y(idx,:),'omitnan'),0,'omitnan');
    end
end
%real version
tempBins = discretize(grouped(exp).temp,tempRange); % assign each temp point to a 0.5degC temp bin
[yR_all,  yR_err] = deal(nan(size(tempRange)));
for i = 1:nBins
    idx  = find(tempBins==i);
    if ~isempty(idx)
        yR_all(i) = mean(mean(y(idx,:),'omitnan'),'omitnan');
        yR_err(i) = std(mean(y(idx,:),'omitnan'),0,'omitnan');
    end
end

fig = getfig('',1);
    hold on
    scatter(tempRange,y_all,50,'r','filled')
    scatter(tempRange,yR_all,50,'k','filled')
    xlabel('temp (\circC)')
    ylabel('distance to food (mm)')
    formatFig(fig); 
save_figure(fig,[figDir grouped(exp).name ' model vs real temp-dist'],'-png')









%%  Bin behavior by temperature...

bintime = 1; % seconds for averaging temp and temperature rate of change
tempRange = 10:0.5:30;
nBins = length(tempRange);

% calculate the instant temp rate: 
tic
for tt = 1:nTemps
    for rr = 1:nRates
        y = models(tt,rr).predicted_temp;
        tempBins = discretize(y,tempRange); % assign each temp point to a 0.5degC temp bin
        tempIdx = [];
        for n = 1:nBins
            tempIdx(n).idx = find(tempBins==n);
        end
        model(tt,rr).tempIdx = tempIdx;
    end
end
toc

% find the avg distance and error for each temp bin now... 
distances = struct;
y = grouped(exp).dist.all;
for rr = 1:nRates
    [distances(rr).avg, distances(rr).err] = deal(nan(nBins,nTemps));
    for tt = 1:nTemps  
          for n = 1:nBins
              idx = model(tt,rr).tempIdx(n).idx;
              if ~isempty(idx)
                  y_trials = mean(y(idx,:),1,'omitnan');
                  distances(rr).avg(n,tt) = mean(y_trials,'omitnan');
                  distances(rr).err(n,tt) = std(y_trials,0,'omitnan');
              end
          end
    end
end




kolor = Color('white', 'blue', nTemps);
x = tempRange;
r = 5; 
c = 6;
fig = getfig('',1);
for rr = 1:nRates
    subplot(r,c,rr)
    hold on
    for tt = nTemps:-10:1
        plot(x, distances(rr).avg(:,tt),'Color', kolor(tt,:))
    end
    plot(grouped(exp).dist.distavgbytemp(:,1), grouped(exp).dist.distavgbytemp(:,2),'color', 'r','LineWidth',1)
    % xlabel('time (min)')
    title([num2str(delta_R(rr)) ' s'],'color', 'k')
end
% add labels
leftedge = 1:r+1:r*c;
bottomrow = (r*c)-c+1:r*c;
for rr = leftedge
    subplot(r,c,rr)
    ylabel('dist (mm)')
end
for rr = bottomrow
    subplot(r,c,rr)
    xlabel('temp (\circC)')
end
for rr = 1:r*c
    subplot(r,c,rr)
    if ~any(rr==bottomrow)
        set(gca,'xcolor','none')
    end
    if ~any(rr==leftedge)
        set(gca, 'ycolor', 'none')
    end
    ylim([10,28])
end




























































































%% ======================================================
% ======================================================
% ======================================================
%% Predictions of CURRENT behavior based on FUTURE conditions
clearvars('-except',initial_vars{:})
models = struct;
initial_vars{end+1} = 'M';
figDir = [saveDir, 'models/'];
if ~exist(figDir,'dir')
    mkdir(figDir)
end

% set up the range of testing points in the 'future' to measure the 'predicted temp'
delta_t = [0, 45]; % range of times into the future to test (min)
skip_size = 1; % how many per minute?
delta_T = (delta_t(1)):skip_size:delta_t(2);
if delta_T(1)==0 % get rid of a zero delay?
    delta_T(1) = [];
end

% ********************************************************************
% ************** !!!!! CAUTION EXTREMELY SLOW !!!!! **************
M = struct;
for exp = 1:num.exp

    raw_dist = grouped(exp).dist.all; % distance to food for all trials in the group
    temp = grouped(exp).temp;
    temprate = nan(size(temp));
    y = data(exp).G(1).TR.data(:,4); % temp rate over time for the experiment
    temprate(1:length(y)) = y;
    
    % run this based on instantaneous rate of change only
    skip_size = 0.5;
    delta_r = [0,5.5];
    delta_R = (delta_r(1)):skip_size:delta_r(2);
    
    nRates = length(delta_R);
    nTemps = length(delta_T);

    % what was the temp and temp rate at dt?
    tic
    for jj = 1:nRates        
        dr = delta_R(jj); % set this as the current temprate integration time region
        for ii = 1:nTemps
            dt = delta_T(ii); % what is the time point (in frame#) from which the temp might be predicted?
            dtDT = smooth(temprate,(delta_R(jj)*fps*60)+1,'moving');
            dummy = temp + (dt.*dtDT);
            % save temperature prediction information
            models(ii,jj).dt = dt;
            models(ii,jj).dr = dr;
            models(ii,jj).predicted_temp = dummy;
        end
    end
    % save data into structure: 
    M(exp).models = models;
    M(exp).delta_R = delta_R;
    M(exp).delta_T = delta_T;
    M(exp).nTemps = nTemps;
    M(exp).nRates = nRates;
    toc
    disp([num2str(exp) '/' num2str(num.exp) ' temp rate integration time'])
end

% find the minimum and maximum predicted temps: 
for exp = 1:num.exp
    minT = []; maxT = [];
    for tt = 1:nTemps
        for rr = 1:nRates
            y = M(exp).models(tt,rr).predicted_temp;
            minT = [minT, min(y)];
            maxT = [maxT, max(y)];
        end
    end
    minT = min(minT);
    maxT = max(maxT);
    M(exp).minT = minT;
    M(exp).maxT = maxT;
end

% ********************************************************************
% ********************************************************************

% test a small integration window: 
for exp = 1:num.exp
    rr = 1;
    kolor = Color('white', 'blue', nTemps);
    time = grouped(exp).time;
    
    fig = getfig('',1);
    hold on
    for tt = 1:nTemps
        plot(time,  M(exp).models(tt,rr).predicted_temp,'Color', kolor(tt,:))
    end
    plot(time, grouped(exp).temp,'color', 'k','LineWidth',1)
    xlabel('time (min)')
    ylabel('temp (\circC)')
    title({grouped(exp).name; ['model predicted temperatures | rate integration: ' num2str(delta_R(rr)) ' minutes']},'color', 'k')
    formatFig(fig);
    save_figure(fig, [figDir grouped(exp).name ' predicted temp lines ' num2str(delta_R(rr)) ' int window'], '-png');
end

%%  FIGURES: visualize the fictive temperature predictions
% save([saveDir 'models'],'models','-v7.3')

for exp = 1:num.exp

    % test a small integration window: 
    kolor = Color('white', 'blue', nTemps);
    time = grouped(exp).time;
    % [nrows, ncols] = subplot_numbers(nRates*nTemps);
    ylims = [];
    r = 3; 
    c = 4;
    fig = getfig('',1);
    for rr = 1:nRates
        subplot(r,c,rr)
        hold on
        for tt = 1:5:nTemps
            plot(time, M(exp).models(tt,rr).predicted_temp,'Color', kolor(tt,:))
        end
        plot(time, grouped(exp).temp,'color', 'k','LineWidth',1)
        title([num2str(delta_R(rr)) ' min'],'color', 'k')
        ylims = [ylims, ylim];

    end
    % add labels
    leftedge = 1:c:r*c;
    bottomrow = (r*c)-c+1:r*c;
    for rr = leftedge
        subplot(r,c,rr)
        ylabel('temp (\circC)')
    end
    for rr = bottomrow
        subplot(r,c,rr)
        xlabel('time (min)')
    end
    for rr = 1:r*c
        subplot(r,c,rr)
        if ~any(rr==bottomrow)
            set(gca,'xcolor','none')
        end
        if ~any(rr==leftedge)
            set(gca, 'ycolor', 'none')
        end
        ylim([min(ylims) max(ylims)])
    end
    
    save_figure(fig,[figDir grouped(exp).name ' model temperature predictions'],'-pdf',1,0);
    save_figure(fig,[figDir grouped(exp).name ' model temperature predictions'],'-png');
end

%% FIGURE: plot side by side the distance vs predicted temp compared to distance vs actual temp
clearvars('-except',initial_vars{:})
figDir = [saveDir, 'models/'];
if ~exist(figDir,'dir')
    mkdir(figDir)
end

rr = 1; % zero lag for starters
skip_size = 0.25; % how many degrees to bin...
LW = 1;

for exp = 1:num.exp
    plotData = [];
    % general parameters for this experiment: 
    tBins = floor(M(exp).minT):skip_size:ceil(M(exp).maxT);
    nBins = length(tBins)-1;
    M(exp).binned.diff = nan([M(exp).nTemps+1,1]);
    tP = getTempTurnPoints(data(exp).temp_protocol);
    % heating and cooling separated average distance to food by predicted temp
    for tt = 0:M(exp).nTemps
        y_cooling = grouped(exp).dist.all(tP.DownROI,:);
        y_heating = grouped(exp).dist.all(tP.UpROI,:);
        if tt > 0
            x_cooling = M(exp).models(tt,rr).predicted_temp(tP.DownROI);
            x_heating = M(exp).models(tt,rr).predicted_temp(tP.UpROI);
        else
            x_cooling = grouped(exp).temp(tP.DownROI);
            x_heating = grouped(exp).temp(tP.UpROI);
        end
        c_idx = discretize(x_cooling,tBins);
        h_idx = discretize(x_heating,tBins);
        [y_cooling_avg, y_cooling_err,y_heating_avg, y_heating_err] = deal(nan(num.trial(exp),1));
        for bin = 1:nBins
            % cooling
            loc = c_idx==bin; % location of temperature points that fit the temp range
            y_cooling_avg(bin) = mean(mean(y_cooling(loc,:),'omitnan'),'omitnan');
            y_cooling_err(bin) = std(mean(y_cooling(loc,:),'omitnan'),'omitnan');
            % heating
            loc = h_idx==bin; % location of temperature points that fit the temp range
            y_heating_avg(bin) = mean(mean(y_heating(loc,:),'omitnan'),'omitnan');
            y_heating_err(bin) = std(mean(y_heating(loc,:),'omitnan'),'omitnan');
        end
        t_idx = tt+1;
        plotData(t_idx).y_cooling = [y_cooling_avg,y_cooling_err];
        plotData(t_idx).y_heating = [y_heating_avg,y_heating_err];
        % measure of 'fit'
        a = median(abs(y_cooling_avg-y_heating_avg),'omitnan'); % total difference between heating & cooling
        a = a*(1/skip_size); % position difference per degree
        plotData(t_idx).diff = a;
        M(exp).binned.diff(t_idx) = a;
    end
    M(exp).binned.data = plotData;
    M(exp).binned.temp = tBins(2:end)+diff(tBins);
end

% PLOT FIGURE: 
for exp = 1:num.exp  
    n = ceil(sqrt(M(exp).nTemps+1));
    [r,c] = subplot_numbers(M(exp).nTemps+1,n);
    c_color = Color('dodgerblue');
    h_color = Color('red');

    x = M(exp).binned.temp;
    fig = getfig(grouped(exp).name,1);
    ylims = [];
    for tt = 0:M(exp).nTemps
        t_idx = tt+1;
        subplot(r,c,t_idx); hold on
        % plot cooling
        y = M(exp).binned.data(t_idx).y_cooling;
        plot_error_fills(true, x, y(:,1), y(:,2), c_color,  fig_type, 0.4);
        plot(x,y(:,1),'color',c_color,'linewidth',LW+1)
        % plot warming
        y = M(exp).binned.data(t_idx).y_heating;
        plot_error_fills(true, x, y(:,1), y(:,2), h_color,  fig_type, 0.4);
        plot(x,y(:,1),'color',h_color,'linewidth',LW+1)
        % labels
        title([num2str(tt) ' | ' num2str(M(exp).binned.data(t_idx).diff)])
        ylims = [ylims,ylim];
    end

    % add labels
    leftedge = 1:c:r*c;
    bottomrow = (r*c)-c+1:r*c;
    for rr = leftedge
        subplot(r,c,rr)
        ylabel('dist (mm)')
    end
    for rr = bottomrow
        subplot(r,c,rr)
        xlabel('temp (\circC)')
    end
    for rr = 1:r*c
        subplot(r,c,rr)
        if ~any(rr==leftedge)
            set(gca,'yticklabel',[])
        end
        ylim([min(ylims),max(ylims)])
    end
    % save the figure
    save_figure(fig,[figDir grouped(exp).name ' model tuning curves'],'-png');
end

% manual guesses at the best fits: 
delay = [0.5,2, 6, 22, 31];
delay = flip(delay);
fig = getfig('',1);
for i = 1:num.exp
    subplot(2,3,i)
    exp = expOrder(i);
    x = [0,M(exp).delta_T];
    y = M(exp).binned.diff;
    scatter(x,y,50,Color('black'),'filled')
    ylim([0, 50])
    v_line(delay(i))
    xlabel('future prediction (min)')
    ylabel('\Delta mm/\circC')
    
    curr_rate = abs(data(exp).G(1).TR.rates(1));
    future_temp(i) = delay(i)*curr_rate;
    disp([grouped(exp).name '  ' num2str(future_temp(i))])

    title({grouped(exp).name; [' ' num2str(future_temp(i)) '\circC ''predicted''']})
end
formatFig(fig,false,[2,3])
save_figure(fig,[figDir 'Simple Model Predictions'],'-png');


%% ANALYSIS & FIGURE: compare the linearity of the data: 
clearvars('-except',initial_vars{:})
% 
buff = 0;
buff = buff*fps;
% R = struct;

for exp = 1:num.exp
    tP = getTempTurnPoints(data(exp).temp_protocol);
    dur = round(mean([diff(tP.down,1,2); diff(tP.up,1,2)])); % duration of ramp up or down ONLY WORKS FOR SINGLE DIRECTION RAMPS
    coolingROI = [tP.down(:,2)-dur-buff, tP.down(:,2)];
    warmingROI = [tP.down(:,2), tP.down(:,2)+dur+buff];
    
    % region that will be considered for the fit correlation
    roi1 = []; roi2 = [];
    for ramp = 1:tP.nDown
        roi1 = [roi1, coolingROI(ramp,1)-dur:coolingROI(ramp,2)];
        roi2 = [roi2, warmingROI(ramp,1):warmingROI(ramp,2)+dur];
    end
    roi = [roi1,roi2];
    roi(roi>length(grouped(exp).time)) = []; % remove any points beyond the length of the experiment
    roi(roi<=0) = [];

    % R(exp).corr = nan([M(exp).nRates,M(exp).nTemps,num.trial(exp)]);
    for trial = 1:num.trial(exp)
        for rr = 1 : 1 % M(exp).nRates
            for tt = 1 : M(exp).nTemps
                % M(exp).delta_R(rr)
                % M(exp).delta_T(tt)
                x = M(exp).models(tt,rr).predicted_temp(roi);
                y = grouped(exp).dist.all(roi,trial);
                loc = isnan(x) | isnan(y);
                x(loc) = []; y(loc) = [];
                
                % Find the R2 value of the fit
                p = polyfit(x, y, 1); % Fit a linear model
                y_fit = polyval(p, x);
                SStot = sum((y - mean(y)).^2);  % Total sum of squares
                SSres = sum((y - y_fit).^2);    % Residual sum of squares
                R2 = 1 - SSres / SStot;        % Coefficient of determination
                R(exp).R2(rr,tt,trial) = R2;
                R(exp).MAE(rr,tt,trial) = mean(abs(y - y_fit)); % mean absolute error
                R(exp).RMSE(rr,tt,trial) = sqrt(mean((y - y_fit).^2)); % root mean square error

                % correlation between temp and distance:
                % R(exp).corr(rr,tt,trial) = corr(x(~loc),y(~loc));                
            end
        end
        x = grouped(exp).temp(roi);
        y = grouped(exp).dist.all(roi,trial);
        R(exp).g_corr(trial) = corr(x(~loc),y(~loc));
        fprintf(['\n ' num2str(trial)])
    end
    disp(['Done ' grouped(exp).name])
end

% Find the RMSE and MAE for the ground truth: 
for exp = 1:num.exp
    for trial = 1:num.trial(exp)
        x = grouped(exp).temp(roi);
        y = grouped(exp).dist.all(roi,trial);
        loc = isnan(x) | isnan(y);
        x(loc) = []; y(loc) = [];
        
        % Find the R2 value of the fit
        p = polyfit(x, y, 1); % Fit a linear model
        y_fit = polyval(p, x);
        SStot = sum((y - mean(y)).^2);  % Total sum of squares
        SSres = sum((y - y_fit).^2);    % Residual sum of squares
        R2 = 1 - SSres / SStot;        % Coefficient of determination
        R(exp).g_R2(trial) = R2;
        R(exp).g_MAE(trial) = mean(abs(y - y_fit)); % mean absolute error
        R(exp).g_RMSE(trial) = sqrt(mean((y - y_fit).^2)); % root mean square error
    end
end


%% 

% for rr = 1:1
% for rr = 1:M(exp).nRates
% 
% for ii = 1:num.exp
%     exp = expOrder(ii);
%     disp(grouped(exp).name)
% end

rr = 1;
% plot the correlations for each of the delays
sz = 10;
avg_sz = sz + 10;
kolor = Color('grey', 'blue', M(exp).nTemps);
r = 2; 
c = 3;
fig = getfig('',1);
for ii = 1:num.exp
    exp = expOrder(ii);
    disp(grouped(exp).name)
    subplot(r,c,ii)
    hold on
    for trial = 1:num.trial(exp)
        C = Color('grey');
        x = M(exp).delta_T;
        % y = squeeze(R(exp).RMSE(rr,:,trial));
        y = squeeze(R(exp).R2(rr,:,trial));
        plot(x,y,'color', C, 'linewidth',0.5)
    end
    % avg across trials
    % y_avg = mean(R(exp).RMSE(rr,:,:),3);
    y_avg = mean(R(exp).R2(rr,:,:),3);

    plot(x,y_avg,'color', 'k', 'linewidth',1)
    % % plot the original data
    x = zeros([1,num.trial(exp)]);
    % y = R(exp).g_RMSE;
    y = R(exp).g_R2;
    scatter(x,y,sz,'r')
    scatter(x(1), mean(y),avg_sz, 'k','filled')
    % h_line(mean(R(exp).g_RMSE),'r','-')
    h_line(mean(R(exp).g_R2),'r','-')

    % [~,idx] = min(y_avg);
    [~,idx] = max(y_avg);
    v_line(M(exp).delta_T(idx),'green','-')
    title([grouped(exp).name],'color', 'k')
end
% add labels
leftedge = 1:c:r*c;
bottomrow = (r*c)-c+1:r*c;
for rr = leftedge
    subplot(r,c,rr)
    % ylabel('RMSE')
    ylabel('R^2')
end
for rr = bottomrow
    subplot(r,c,rr)
    xlabel('future prediction (min)')
end
for rr = 1:r*c
    subplot(r,c,rr)
    if ~any(rr==bottomrow)
        set(gca,'xcolor','none')
    end
    if ~any(rr==leftedge)
        set(gca, 'ycolor', 'none')
    end
    % ylim([-0.9, -0.3])
    xlim([-2,45])
end
formatFig(fig, false, [r,c]);

save_figure(fig,[figDir 'Model R2 rate window ' num2str(M(exp).delta_R(rr))],'-png');


%% Plot the 'best' fit of the model vs the actual data


[low, high] = deal([]);
for exp = 1:num.exp
    for tt = 1:M(exp).nTemps
        for rr = 1:M(exp).nRates
            a = M(exp).models(tt,rr).predicted_temp;
            low = [low, min(a)];
            high = [high, max(a)];
        end
    end
end
bins = floor(min(low)):0.5:ceil(max(high));

rr = 1;
plotData = [];
for exp = 1:num.exp

    % TODO: plot this by temp & color for original heating vs cooling...
    cROI = []; hROI = [];
    for ramp = 1:tP.nDown
        % cooling
        roi = coolingROI(ramp,1):coolingROI(ramp,2);
        cROI = [cROI; roi'];
        % warming
        roi = warmingROI(ramp,1):warmingROI(ramp,2);
        hROI = [hROI; roi'];
    end

    % cooling
    x = grouped(exp).temp(cROI);
    y = grouped(exp).dist.all(cROI,:);
    a = discretize(x,bins);
    MT = nan(length(bins),2);
    for bin = 1:length(bins)
        loc = find(a==bin);
        if ~isempty(loc)
            b = mean(y(loc,:),2,'omitnan');
            MT(bin,1) = mean(b,'omitnan');
            MT(bin,2) = std(b,'omitnan');
        end
    end
    plotData(exp).control.cooling = MT;

    % warming
    x = grouped(exp).temp(hROI);
    y = grouped(exp).dist.all(hROI,:);
    a = discretize(x,bins);
    MT = nan(length(bins),2);
    for bin = 1:length(bins)
        loc = find(a==bin);
        if ~isempty(loc)
            b = mean(y(loc,:),2,'omitnan');
            MT(bin,1) = mean(b,'omitnan');
            MT(bin,2) = std(b,'omitnan');
        end
    end
    plotData(exp).control.warming = MT;

    % model data
    
    %find the best deltaT to use for the model
    [R2,idx] = max(mean(R(exp).R2(rr,:,:),3));
    plotData(exp).model.R2 = R2;
    x = M(exp).models(rr,idx).predicted_temp;
    y = grouped(exp).dist.all;
    cY = y(cROI,:);
    hY = y(hROI,:);
    cX = x(cROI);
    hX = x(hROI);
    %cooling
    a = discretize(cX,bins); % find fictive temp bins
    MT = nan(length(bins),2);
    for bin = 1:length(bins)
        loc = find(a==bin);
        if ~isempty(loc)
            b = mean(cY(loc,:),2,'omitnan');
            MT(bin,1) = mean(b,'omitnan');
            MT(bin,2) = std(b,'omitnan');
        end
    end
    plotData(exp).model.cooling = MT;
    %warming
    a = discretize(hX,bins); % find fictive temp bins
    MT = nan(length(bins),2);
    for bin = 1:length(bins)
        loc = find(a==bin);
        if ~isempty(loc)
            b = mean(hY(loc,:),2,'omitnan');
            MT(bin,1) = mean(b,'omitnan');
            MT(bin,2) = std(b,'omitnan');
        end
    end
    plotData(exp).model.warming = MT;

    sz = 35;
    fig = getfig('',1);
    subplot(1,2,1)
    hold on
    x = bins';
    y = plotData(exp).control.cooling(:,1);
    scatter(x,y,sz,Color('dodgerblue'),'filled')
    y = plotData(exp).control.warming(:,1);
    scatter(x,y,sz,Color('red'),'filled')
    xlabel('temp (\circC)')
    ylabel('distance to food (mm)')
    
    subplot(1,2,2)
    hold on
    x = bins';
    y = plotData(exp).model.cooling(:,1);
    scatter(x,y,sz,Color('dodgerblue'),'filled')
    y = plotData(exp).model.warming(:,1);
    scatter(x,y,sz,Color('red'),'filled')
    xlabel('predicted temp (\circC)')
    ylabel('distance to food (mm)')

    formatFig(fig,false,[1,2]);

    save_figure(fig,[figDir grouped(exp).name 'Model R2 vs real ' num2str(M(exp).delta_R(rr))],'-png');
end


% 
% save(['E:\model_predictions.mat'],'R', 'M')



%% WORKING BELOW HERE


%%

exp = 1;

% plot the correlations for each of the delays
sz = 10;
avg_sz = sz + 10;
kolor = Color('grey', 'blue', M(exp).nTemps);
r = 3; 
c = 4;
fig = getfig('',1);
for rr = 1:M(exp).nRates
    subplot(r,c,rr)
    hold on
    for trial = 1:num.trial(exp)
        C = Color('grey');
        x = M(exp).delta_T;
        y = squeeze(R(exp).RMSE(rr,:,trial));
        plot(x,y,'color', C, 'linewidth',0.5)
    end
    % avg across trials
    y_avg = mean(R(exp).RMSE(rr,:,:),3);
    plot(x,y_avg,'color', 'k', 'linewidth',1)
    % % plot the original data
    % x = zeros([1,num.trial(exp)]);
    % y = R(exp).g_corr;
    % scatter(x,y,sz,'r')
    % scatter(x(1), mean(y),avg_sz, 'k','filled')
    % h_line(mean(R(exp).g_corr),'r','-')
    [~,idx] = min(y_avg);
    v_line(M(exp).delta_T(idx),'green','-')
    title([num2str(M(exp).delta_R(rr)) 'min | min ' num2str(M(exp).delta_T(idx))],'color', 'k')
end
% add labels
leftedge = 1:c:r*c;
bottomrow = (r*c)-c+1:r*c;
for rr = leftedge
    subplot(r,c,rr)
    ylabel('correlation')
end
for rr = bottomrow
    subplot(r,c,rr)
    xlabel('future prediction (min)')
end
for rr = 1:r*c
    subplot(r,c,rr)
    if ~any(rr==bottomrow)
        set(gca,'xcolor','none')
    end
    if ~any(rr==leftedge)
        set(gca, 'ycolor', 'none')
    end
    % ylim([-0.9, -0.3])
    xlim([-2,30])
end
formatFig(fig, false, [r,c]);

save_figure(fig,[figDir grouped(exp).name ' model correlations for all times'],'-png');




% plot the best model fit for the data: 
for exp = 1:num.exp
    y = mean(R(exp).corr,3,'omitnan');
    minC = [];
    for rr = 1:M(exp).nRates
        [minC(rr), idx] = min(y(rr,:));
        R(exp).best_dt(rr) = idx; 
    end
    R(exp).optimal_dt = min(minC);
end

% for now, just plot for an instantaneous rate: 
fig = getfig('',1); hold on
% TODO: plot this by temp & color for original heating vs cooling...
for ramp = 1:tP.nDown
    % cooling
    roi = coolingROI(ramp,1):coolingROI(ramp,2);
    plot(x(roi),y(roi),'color',Color('dodgerBlue'))
    % warming
    roi = warmingROI(ramp,1):warmingROI(ramp,2);
    plot(x(roi),y(roi),'color',Color('red'))
end


% [~,idx] = min(min(y,1));
% 
%     x = M(exp).delta_T;
% 
%     plot(x,y,'color', C, 'linewidth',0.5)
% end
% % avg across trials
% y_avg = mean(R(exp).corr(rr,:,:),3);
% 
% 
% fig = getfig('',1); hold on
%     for ramp = 1:tP.nDown
%         % cooling
%         roi = coolingROI(ramp,1):coolingROI(ramp,2);
%         plot(x(roi),y(roi),'color',Color('dodgerBlue'))
%         % warming
%         roi = warmingROI(ramp,1):warmingROI(ramp,2);
%         plot(x(roi),y(roi),'color',Color('red'))
%     end

%% NEW: find and plot the MSE of each delta-T to determine the model with the best fit:





















%% ANALYSIS: correlation for measure of linearity

% Screen out the recovery period: 
% pre / post buffer: 16 mins
buff = 16;
buff = buff*fps;
tP = getTempTurnPoints(data(exp).temp_protocol);
roi = tP.down(:,1)-buff : tP.up(:,2)-buff; % selected region around temp ramps

% CAUTION: SLOW
correlations = nan([nRates,nTemps,num.trial(exp)]);
for rr = 1:nRates
    for tt = 1:nTemps
        A = models(tt,rr).predicted_temp(roi);
        for trial = 1:num.trial(exp)
            B = [A, grouped(exp).dist.all(roi,trial)];
            R = corrcoef(B,'Rows','pairwise');
            correlations(rr,tt,trial) = R(1,2);
        end
    end
    disp([num2str(rr) '/' num2str(nRates)])
end
% save([saveDir 'correlations'],'correlations','-v7.3')

% ground truth correlation: 
g_corr = nan([num.trial(exp),1]);
 for trial = 1:num.trial(exp)
    B = [grouped(exp).temp, grouped(exp).dist.all(:,trial)];
    R = corrcoef(B,'Rows','pairwise');
    g_corr(trial) = R(1,2);
end


% plot the correlations for each of the delays
sz = 10;
avg_sz = sz + 10;
kolor = Color('grey', 'blue', nTemps);
r = 3; 
c = 4;
fig = getfig('',1);
for rr = 1:nRates
    subplot(r,c,rr)
    hold on
    for trial = 1:num.trial(exp)
        C = Color('grey');
        x = delta_T;
        y = squeeze(correlations(rr,:,trial));
        plot(x,y,'color', C, 'linewidth',0.5)
    end
    % avg across trials
    y_avg = mean(correlations(rr,:,:),3);
    plot(x,y_avg,'color', 'k', 'linewidth',1)
    % plot the original data
    x = zeros([1,num.trial(exp)]);
    y = g_corr;
    scatter(x,y,sz,'r')
    scatter(x(1), mean(y),avg_sz, 'k','filled')
    h_line(mean(g_corr),'r','-')
    [~,idx] = min(y_avg);
    v_line(delta_T(idx),'green','-')
    title([num2str(delta_R(rr)) 's | min ' num2str(delta_T(idx))],'color', 'k')
end
% add labels
leftedge = 1:c:r*c;
bottomrow = (r*c)-c+1:r*c;
for rr = leftedge
    subplot(r,c,rr)
    ylabel('corr coeff')
end
for rr = bottomrow
    subplot(r,c,rr)
    xlabel('prediction delay (sec)')
end
for rr = 1:r*c
    subplot(r,c,rr)
    if ~any(rr==bottomrow)
        set(gca,'xcolor','none')
    end
    if ~any(rr==leftedge)
        set(gca, 'ycolor', 'none')
    end
    ylim([-1, 0])
    % xlim([-50,3600])
end
save_figure(fig,[figDir grouped(exp).name ' model correlations for all times'],'-png');


%  Plot the strongest correlation for each type/trial...
kolor = Color('blue','yellow',nRates);
min_corr = [];
past_time = [];
for rr = 1:nRates
     y_avg = mean(correlations(rr,:,:),3);
     [~,idx] = min(y_avg);
     min_corr(rr) = y_avg(idx);
     past_time(rr) = delta_T(idx);
end
corr_diff =  mean(g_corr) - min_corr;
[~, max_idx] = max(corr_diff); % 'best' model prediction time
delay_t = past_time(max_idx);
integration_t = delta_R(max_idx);

fig = getfig('',1);
hold on
scatter(past_time,corr_diff,50,kolor,"filled")
h_line(0,'k',':')
% plot the best model
scatter(past_time(max_idx),corr_diff(max_idx),100,'r')
title(['Highest correlation: ' num2str(delay_t) 's prior with ' num2str(integration_t) 's temp rate integration window'])
xlabel('time in past for current temp prediction (s)')
ylabel('Difference in correlation from baseline')
formatFig(fig);

save_figure(fig,[figDir grouped(exp).name ' model performance'],'-png')

% plot the behavior of the best model vs the real data...
tempRange = 10:0.1:30;
nBins = length(tempRange);

% find best model numbers: 
T_idx = find(delay_t==delta_T);
R_idx = max_idx;
x = models(T_idx,R_idx).predicted_temp; % predicted temp
y = grouped(exp).dist.all; %distance to food
tempBins = discretize(x,tempRange); % assign each temp point to a 0.5degC temp bin
%predicted version
[y_all,  y_err] = deal(nan(size(tempRange)));
for i = 1:nBins
    idx  = find(tempBins==i);
    if ~isempty(idx)
        y_all(i) = mean(mean(y(idx,:),'omitnan'),'omitnan');
        y_err(i) = std(mean(y(idx,:),'omitnan'),0,'omitnan');
    end
end
%real version
tempBins = discretize(grouped(exp).temp,tempRange); % assign each temp point to a 0.5degC temp bin
[yR_all,  yR_err] = deal(nan(size(tempRange)));
for i = 1:nBins
    idx  = find(tempBins==i);
    if ~isempty(idx)
        yR_all(i) = mean(mean(y(idx,:),'omitnan'),'omitnan');
        yR_err(i) = std(mean(y(idx,:),'omitnan'),0,'omitnan');
    end
end

fig = getfig('',1);
    hold on
    scatter(tempRange,y_all,50,'r','filled')
    scatter(tempRange,yR_all,50,'k','filled')
    xlabel('temp (\circC)')
    ylabel('distance to food (mm)')
    formatFig(fig); 
save_figure(fig,[figDir grouped(exp).name ' model vs real temp-dist'],'-png')









%%  Bin behavior by temperature...

bintime = 1; % seconds for averaging temp and temperature rate of change
tempRange = 10:0.5:30;
nBins = length(tempRange);

% calculate the instant temp rate: 
tic
for tt = 1:nTemps
    for rr = 1:nRates
        y = models(tt,rr).predicted_temp;
        tempBins = discretize(y,tempRange); % assign each temp point to a 0.5degC temp bin
        tempIdx = [];
        for n = 1:nBins
            tempIdx(n).idx = find(tempBins==n);
        end
        model(tt,rr).tempIdx = tempIdx;
    end
end
toc

% find the avg distance and error for each temp bin now... 
distances = struct;
y = grouped(exp).dist.all;
for rr = 1:nRates
    [distances(rr).avg, distances(rr).err] = deal(nan(nBins,nTemps));
    for tt = 1:nTemps  
          for n = 1:nBins
              idx = model(tt,rr).tempIdx(n).idx;
              if ~isempty(idx)
                  y_trials = mean(y(idx,:),1,'omitnan');
                  distances(rr).avg(n,tt) = mean(y_trials,'omitnan');
                  distances(rr).err(n,tt) = std(y_trials,0,'omitnan');
              end
          end
    end
end




kolor = Color('white', 'blue', nTemps);
x = tempRange;
r = 5; 
c = 6;
fig = getfig('',1);
for rr = 1:nRates
    subplot(r,c,rr)
    hold on
    for tt = nTemps:-10:1
        plot(x, distances(rr).avg(:,tt),'Color', kolor(tt,:))
    end
    plot(grouped(exp).dist.distavgbytemp(:,1), grouped(exp).dist.distavgbytemp(:,2),'color', 'r','LineWidth',1)
    % xlabel('time (min)')
    title([num2str(delta_R(rr)) ' s'],'color', 'k')
end
% add labels
leftedge = 1:r+1:r*c;
bottomrow = (r*c)-c+1:r*c;
for rr = leftedge
    subplot(r,c,rr)
    ylabel('dist (mm)')
end
for rr = bottomrow
    subplot(r,c,rr)
    xlabel('temp (\circC)')
end
for rr = 1:r*c
    subplot(r,c,rr)
    if ~any(rr==bottomrow)
        set(gca,'xcolor','none')
    end
    if ~any(rr==leftedge)
        set(gca, 'ycolor', 'none')
    end
    ylim([10,28])
end


























































































