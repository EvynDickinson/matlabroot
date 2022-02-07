

% Run this script after initial data processing with:
% QuadStep1.m
% QuadStep2.m
% GroupDataGUI
clear 

%% Load data file
% Load data for now from the first step of QuadStep3.m

clearvars('-except',initial_vars{:})

%% Trial by trial figures: 
disp_figs = true;
auto_save_figs = true;
ind_fig_loc = [figDir 'Trial by trial\'];
if ~isfolder(ind_fig_loc); mkdir(ind_fig_loc); end
vars = [initial_vars(:)', 'trial', 'threshHigh', 'threshLow', 'binSpace',...
        'G', 'disp_figs','auto_save_figs','ind_fig_loc'];

[threshHigh, threshLow] = getTempThresholds;
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 

G = struct;

% Get the temp rate, temp, and distance from food
for trial = 12:ntrials
    T_rates = [];
    % temperature
    temp = data(trial).occupancy.temp;
    plotData(:,1) = temp(1:end-1); 
    % distance from food
    for well = 1:4
        [kolor,num] = pullFoodColor(data(trial).wellLabels{well});
        if num==1
            y = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm;
            plotData(:,3) = y(1:end-1);
            break
        end
    end

    % rate of temperature change
    dT = diff(smooth(temp,180)); %180 = 1 minute smoothing kernal
    dt = diff(data(trial).occupancy.time); 
    plotData(:,2) = dT./dt'; 

    % find the mean temp rate during each ramp period:
    tPoints = getTempTurnPoints(T.TempProtocol{trial}); %accomodates multiple temp protocols within the data group
    for ii = 1:length(tPoints.transitions)-1
        roi = tPoints.transitions(ii):tPoints.transitions(ii+1);
        distD = median(plotData(roi,2));
        plotData(roi,4) = round(distD*ones(range(roi)+1,1),2);
    end

    % === Demo plot for single trial =====:
    title_str = [T.Date{trial} ' Arena ' T.Arena{trial}];
    if disp_figs 
        fig_file = [ind_fig_loc ExpGroup ' temp temp_rate dist ' title_str];
        if isfile([fig_file '.png']) %skip already generated figures
            continue
        end
        time = data(trial).occupancy.time(1:end-1);
        fig = figure; set(fig, 'pos', [689 48 1628 960])
        % temp
        subplot(3,1,1)
        plot(time, plotData(:,1),'color','w')
        v_line(time(tPoints.transitions),'y',':')
        ylabel('Temp (\circC)')
        title(title_str)
        % temp rate
        subplot(3,1,2); hold on
        plot(time, plotData(:,2),'color','w')
        plot(time, plotData(:,4),'color','r','linewidth', 2)
        v_line(time(tPoints.transitions),'y',':')
        ylabel('dT/dt (\circC/min)')
        % distance from food
        subplot(3,1,3)
        plot(time, plotData(:,3),'color',kolor)
        v_line(time(tPoints.transitions),'y',':')
        xlabel('time (min)')
        ylabel('distance (mm)')
        yyaxis right
        plot(time, smooth(data(trial).occupancy.occ(1:end-1,well),180),'color', 'w')
        ylabel('occupancy (prob)')
        fig = formatFig(fig, true, [3,1]);
        subplot(3,1,3)
        yyaxis left 
        set(gca,'YColor', 'w')
        save_figure(fig, fig_file, '-png',auto_save_figs);
    end
    % Temp-rate identification and sorting: 
    buffSize = 0.05;
    for ii = 1:tPoints.nRates
        edges(ii,:) = [tPoints.rates(ii)-buffSize, tPoints.rates(ii)+buffSize];
    end  
    rateData = plotData(:,4);
    nRates = tPoints.nRates;
    rateIdx = discretize(rateData,nRates);
    for ii = 1:nRates
        TD = rateData(rateIdx==ii);
        T_rates(ii) = round(mean(TD),2);
        % do these match the assumed rates?
        idx = find(T_rates(ii) > edges(:,1) & T_rates(ii) < edges(:,2));
        if isempty(idx)
            warndlg('Temp rate not within expected range')
            return
        end
        % set the rate value to the uniform cross-group value:
        rateData(rateIdx==ii) = tPoints.rates(ii);
        T_rates(ii) = tPoints.rates(ii);
%        figure; hold on
%        plot(rateData)
%        plot(plotData(:,4))

    end
    plotData(:,4) = rateData;
    plotData(:,5) = rateIdx; %rate bin
    G(trial).rateIdx = rateIdx;
    G(trial).rates = T_rates;
    G(trial).data = plotData;
    
    % Temperature range and bin size formatting
    t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
    if t_roi(end)<ceil(threshHigh)
        t_roi(end+1) = ceil(threshHigh) + binSpace;
    end
    nTemps = length(t_roi);
    tempIdx = discretize(plotData(:,1), t_roi);
    G(trial).nTemps = nTemps;
    G(trial).tempIdx = tempIdx;
    G(trial).temps = t_roi;

    % turn coordinates of heatmap into a vector to fill in 2D
    heatMapData = [];
    for col = 1:nTemps
        for row = 1:nRates
        % Pull the data that matches this category
        loc = (rateIdx==row) & (tempIdx==col);
        heatMapData(row,col) = mean(plotData(loc,3));
        end
    end
    G(trial).heatmap = heatMapData;

    % PLOT DISTANCE AS FUNCTION OF TEMP AND RATE OF TEMP
    if disp_figs
        fig_file = [ind_fig_loc ExpGroup ' temp temp_rate dist heatmap ' title_str];
        if isfile([fig_file '.png']) %skip already generated figures
            continue
        end
%         colormap default
        fig = figure; set(fig, 'pos', [560 127 983 417]);
        hold on
        imAlpha=ones(size(heatMapData));
        imAlpha(isnan(heatMapData))=0;
        imagesc(heatMapData,'AlphaData',imAlpha);
        set(gca,'color',0*[1 1 1]);
        axis tight
        title(title_str)
        % Axes formatting
        ax = gca;
        fig = formatFig(fig, true);
        XtickNum = ax.XTick;
        ax.XTickLabel = t_roi(XtickNum);
        YtickNum = ax.YTick;
        set(gca, 'ytick', 1:nRates,'YTickLabel',T_rates)
        ylabel('\DeltaT/dt (\circC/min)')
        xlabel('Temp (\circC)')
        % Colorbar formatting
        cbh = colorbar(); 
        cbh.Label.String = 'Distance from food (mm)';
        cbh.Color = Color('white');
        % flip colormap around to make yellow closer to food
        cmap = colormap;
        set(gca, 'colormap', flip(cmap))
        save_figure(fig, fig_file, '-png', auto_save_figs);
    end
    
clearvars('-except',vars{:}) 
end

% Save Data structure:
if questdlg('Save loaded data?')
    save([figDir ExpGroup ' temp rate'])
end



%% Grouped hysteresis heatmap output 
% TODO: reformat this to take advantage of the previously processed data
% Fly by fly average:
clearvars('-except',vars{:}) 


% Find the total number and id of temp rates:
allRates=[];
for trial = 1:ntrials
    allRates = [allRates,G(trial).rates];
end
tRates = sort(unique(allRates));
nRates = length(tRates);
nTemps = G(1).nTemps; %these should all be the same since they're held constant above

% Group the data for each temp rate
tempData = nan(nRates,nTemps,ntrials);
for trial = 1:ntrials
    for rr = 1:nRates
        idx = find(G(trial).rates==tRates(rr));
        if isempty(idx)
            continue
        end
        tempData(rr,:,trial) = G(trial).heatmap(idx,:);
    end
end

% Find the mean and err of each temp bin:
plotData.avg = mean(tempData,3,'omitnan');
plotData.err = std(tempData,0,3,'omitnan')./sqrt(ntrials);


% ========= HeatMap of dT/dt vs T =============
heatMapData = plotData.avg;
t_roi = G(trial).temps;
title_str = [ExpGroup ' (n = ' num2str(ntrials) ')'];

fig = figure; set(fig, 'pos', [560 127 983 417]);
    hold on
    imAlpha=ones(size(heatMapData));
    imAlpha(isnan(heatMapData))=0;
    imagesc(heatMapData,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    axis tight
    title(title_str)
    % Axes formatting
    ax = gca;
    fig = formatFig(fig, true);
    XtickNum = ax.XTick;
    ax.XTickLabel = t_roi(XtickNum);
    YtickNum = ax.YTick;
    ax.YTickLabel = tRates(YtickNum);
    ylabel('\DeltaT/dt (\circC/min)')
    xlabel('Temp (\circC)')
    % Colorbar formatting
    cbh = colorbar(); 
    cbh.Label.String = 'Distance from food (mm)';
    cbh.Color = Color('white');
    % flip colormap around to make yellow closer to food
    cmap = colormap;
    set(gca, 'colormap', flip(cmap))

save_figure(fig, [figDir ExpGroup ' temp temp_rate dist heatmap ' ExpGroup], '-png');

% ========== Line plots of each rate comparison ==============
cList = {'Orange', 'DarkViolet', 'Turquoise','white', 'Turquoise', 'DarkViolet', 'Orange'};
LS = {'--','-.','-'}; %cooling|stationary|heating

fig = figure;
hold on
for rr = 1:nRates
    if tRates(rr)>0
        lstyle = LS{3};
    elseif tRates(rr)<0
        lstyle = LS{1};
    elseif tRates(rr)==0
        continue
        lstyle = LS{2};
    end
    x = t_roi;
    y = plotData.avg(rr,:);
    plot(x,y,'color', Color(cList{rr}), 'linewidth', 2, 'linestyle', lstyle)
end
title(ExpGroup)
ylabel('Distance to food (mm)')
xlabel('Temperature (\circC)')

fig = formatFig(fig, true);
save_figure(fig, [figDir ExpGroup ' temp temp_rate dist all rates ' ExpGroup], '-png');



% ========== Line plots of separated by rate MANUAL ADJUST FOR MORE RATES ==============
cList = {'Orange', 'DarkViolet', 'Turquoise','white', 'Turquoise', 'DarkViolet', 'Orange'};
LS = {'--','-.','-'}; %cooling|stationary|heating
legStr = {'SEM','Cooling','SEM','Heating'};
ncol = 2;
nrow = 2;

fig = figure; set(fig, 'pos',[296 35 1224 956])
% All trials
    subplot(nrow, ncol, 1); hold on
    for rr = 1:nRates
        if tRates(rr)>0
            lstyle = LS{3};
        elseif tRates(rr)<0
            lstyle = LS{1};
        elseif tRates(rr)==0
            lstyle = LS{2};
            continue
        end
        x = t_roi;
        y = plotData.avg(rr,:);
        plot(x,y,'color', Color(cList{rr}), 'linewidth', 2, 'linestyle', lstyle)
    end
    title([ExpGroup ' (n = ' num2str(ntrials) ')'])
    ylabel('Distance to food (mm)')
    xlabel('Temperature (\circC)')
% Fast
    subplot(nrow, ncol, 2); hold on
    for rr = [1,7]
        if tRates(rr)>0
            lstyle = LS{3};
        elseif tRates(rr)<0
            lstyle = LS{1};
        elseif tRates(rr)==0
            lstyle = LS{2};
            continue
        end
        x = t_roi(1:end-1);
        y = plotData.avg(rr,1:end-1);
        kolor = Color(cList{rr});
        err = plotData.err(rr,1:end-1);
        fill_data = error_fill(x, y, err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
    end
    title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
    ylabel('Distance to food (mm)')
    xlabel('Temperature (\circC)')
    l = legend(legStr);
    set(l,'textcolor', 'w','box', 'off')
% Medium
    subplot(nrow, ncol, 3); hold on
    for rr = [2,6]
        if tRates(rr)>0
            lstyle = LS{3};
        elseif tRates(rr)<0
            lstyle = LS{1};
        elseif tRates(rr)==0
            lstyle = LS{2};
            continue
        end
        x = t_roi(1:end-1);
        y = plotData.avg(rr,1:end-1);
        kolor = Color(cList{rr});
        err = plotData.err(rr,1:end-1);
        fill_data = error_fill(x, y, err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
    end
    title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
    ylabel('Distance to food (mm)')
    xlabel('Temperature (\circC)')
    l = legend(legStr);
    set(l,'textcolor', 'w','box', 'off')
% Slow
    subplot(nrow, ncol, 4); hold on
    for rr = [3,5]
        if tRates(rr)>0
            lstyle = LS{3};
        elseif tRates(rr)<0
            lstyle = LS{1};
        elseif tRates(rr)==0
            lstyle = LS{2};
            continue
        end
        x = t_roi(1:end-1);
        y = plotData.avg(rr,1:end-1);
        kolor = Color(cList{rr});
        err = plotData.err(rr,1:end-1);
        fill_data = error_fill(x, y, err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
    end
    title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
    ylabel('Distance to food (mm)')
    xlabel('Temperature (\circC)')
    l = legend(legStr);
    set(l,'textcolor', 'w','box', 'off')
    
fig = formatFig(fig, true,[nrow,ncol]);
save_figure(fig, [figDir ExpGroup ' temp rate food preference dependence ' ExpGroup], '-png');

clearvars('-except',initial_vars{:})






%% CLOSER LOOK: confounding variables | circumstances
% TODO: set this to demo across individual trials to see how it correlates...

disp_figs = true;
auto_save_figs = true;
ind_fig_loc = [figDir 'Trial by trial\'];
if ~isfolder(ind_fig_loc); mkdir(ind_fig_loc); end
vars = [initial_vars(:)', 'trial', 'threshHigh', 'threshLow', 'binSpace',...
        'G', 'disp_figs','auto_save_figs','ind_fig_loc'];
tempRange = 1:30;

for trial = 1:ntrials
    title_str = [T.Date{trial} ' Arena ' T.Arena{trial}];
    fig_file = [ind_fig_loc ExpGroup ' fly count vs temp ' title_str];
    if isfile([fig_file '.png']) %skip already generated figures
        continue
    end
    fig = figure;
    idx = discretize(data(trial).occupancy.temp,tempRange);
    boxplot(data(trial).occupancy.flycount, idx,'Plotstyle', 'traditional','colors', 'w')
    h_line(T.NumFlies(trial), 'Teal', '-',3)
    xlabel('Temperature (\circC)')
    ylabel('Tracking Fly Count')
    fig = formatFig(fig, true);

    save_figure(fig, fig_file, '-png',true);
end
































