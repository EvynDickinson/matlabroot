

%% Jumping escape response analysis
% just run the time and tuning curves on the baseline 6.2 scripts

%% FIGURES: location of fly jumps
clearvars('-except',initial_var{:})

region_list = {'innerFoodQuad', 'OutterRing', 'innerEmptyQuad'};
region_names = {'food quad', 'escape ring', 'inner arena'};
color_list = {'vaporwavepurple', 'vaporwavegren', 'grey'};
sz = 50;
buff = 0.2;

fig = getfig('', 1, [530 718]); hold on
for ii = 1:length(region_list)
    % pull the number of jumps for each location type
    idx = data.jump;
    % absolute number of jumps per region
    loc = logical(replaceNaN(data.(region_list{ii}),false));
    jumps_in_roi = idx & loc;
    a = squeeze(sum(jumps_in_roi,1));
    nJumps = sum(a,1);
    % scatter jumps per trial
    x = ii*ones(size(nJumps));
    scatter(x, nJumps, sz, Color(color_list{ii}), "filled", ...
        MarkerFaceAlpha = 0.7,XJitter = "density",...
        XJitterWidth = 0.2)
    % average lines
    y_avg = mean(nJumps, 'omitnan');
    plot([ii-buff, ii+buff], [y_avg, y_avg], Color = foreColor,...
        LineWidth = 2)
end
% formatting:
formatFig(fig, blkbgd);
set(gca, XTick = 1:length(region_list), XTickLabel = region_names)
ylabel('Total number of jumps per region')

save_figure(fig, [figDir, 'Jumps per region scatter']);



% FIGURE: normalized within regions
% jump rate -- what percentage of time in the region are they jumping for
% escape? 


fig = getfig('', 1, [530 718]); hold on
for ii = 1:length(region_list)
    % pull the number of jumps for each location type
    idx = data.jump;
    % absolute number of jumps per region
    loc = logical(replaceNaN(data.(region_list{ii}),false));
    % total time in region: 
    frames_in_roi = sum(squeeze(sum(loc,1)),1);
    % total jumps 
    jumps_in_roi = idx & loc;
    a = squeeze(sum(jumps_in_roi,1));
    nJumps = sum(a,1);
    % jump rate: 
    jumpRate = (nJumps./frames_in_roi)*fps*60; %convert to jumps per minute
    % scatter jumps per trial
    x = ii*ones(size(nJumps));
    scatter(x, jumpRate, sz, Color(color_list{ii}), "filled", ...
        MarkerFaceAlpha = 0.7,XJitter = "density",...
        XJitterWidth = 0.2)
    % average lines
    y_avg = mean(jumpRate, 'omitnan');
    plot([ii-buff, ii+buff], [y_avg, y_avg], Color = foreColor,...
        LineWidth = 2)
end
% formatting:
formatFig(fig, blkbgd);
set(gca, XTick = 1:length(region_list), XTickLabel = region_names)
ylabel('Jumps Per Minute')
save_figure(fig, [figDir, 'Jump Rate scatter']);


% Total time spent in region across experiment: 
fig = getfig('', 1, [530 718]); hold on
for ii = 1:length(region_list)
    loc = logical(replaceNaN(data.(region_list{ii}),false));
    % total time in region: 
    frames_in_roi = sum(squeeze(sum(loc,1)),1);
    time_in_roi = frames_in_roi/(fps*60);
    % scatter jumps per trial
    x = ii*ones(size(nJumps));
    scatter(x, time_in_roi, sz, Color(color_list{ii}), "filled", ...
        MarkerFaceAlpha = 0.7,XJitter = "density",...
        XJitterWidth = 0.2)
    % average lines
    y_avg = mean(time_in_roi, 'omitnan');
    plot([ii-buff, ii+buff], [y_avg, y_avg], Color = foreColor,...
        LineWidth = 2)
end
% formatting:
formatFig(fig, blkbgd);
set(gca, XTick = 1:length(region_list), XTickLabel = region_names)
ylabel('Time in Region (min)')
save_figure(fig, [figDir, 'Time in region scatter']);

%% Is there a higher instance of jump response of females near males?
clearvars('-except',initial_var{:})
sz = 50; % scatter point size
FA = 0.7; % face alpha level
LW = 1.5; % line width

% compare the distribution of distances between flies to 
% the distribution of IFD during starts 

% extract the interfly distances for jumps: 
[IFD, f_IFD, m_IFD] = deal(data.IFD);
F_jump = squeeze(data.jump(:,F,:));
f_IFD(~F_jump) = nan;
M_jump = squeeze(data.jump(:,M,:));
m_IFD(~M_jump) = nan;
% averages
IFD_avg = mean(IFD,1, 'omitnan');
f_IFD_avg = mean(f_IFD, 1, 'omitnan');
m_IFD_avg = mean(m_IFD, 1, 'omitnan');

% compare the distributions with MCC
alpha = 0.5/num.trials; % multiple comparisons corrections
[~,mP] = ttest2(IFD,m_IFD);
m_sig = mP<=alpha; % signifant male jump IFD distributions
[~,fP] = ttest2(IFD,f_IFD);
f_sig = fP<=alpha; % signifant female jump IFD distributions

% plot type depends on significance: 
[m_type, f_type] = deal(repmat({'-'},[num.trials,1]));
m_type(~m_sig) = {'--'}; % non significant trials are dashed lines
f_type(~f_sig) = {'--'};


% quick plot of the different values: 
fig = getfig('', 1, [530 718]); hold on
    x = repmat([1,2],[num.trials,1]);
    fY = [IFD_avg; f_IFD_avg];
    mY = [IFD_avg; m_IFD_avg];
    scatter(1, IFD_avg, sz, foreColor, 'filled', MarkerFaceAlpha=FA)
    scatter(2, f_IFD_avg, sz, data.color(F,:),...
        'filled', MarkerFaceAlpha=FA)
    scatter(2, m_IFD_avg, sz, data.color(M,:),...
        'filled', MarkerFaceAlpha=FA)
    ff = plot(x', fY, Color=data.color(F,:), LineWidth=LW);
    set(ff, {'LineStyle'}, f_type);
    mm = plot(x', mY, Color=data.color(M,:), LineWidth=LW);
    set(mm, {'LineStyle'}, m_type);

    xlim([0.5,2.5])
    formatFig(fig, blkbgd);
    set(gca, XTick=[1,2],XTickLabel={'all', 'jumps only'})
    ylabel('distance between flies (mm)')
    save_figure(fig, [figDir, 'IFD during jumps scatter']);


% is there significance across all the trials?
[~, p] = ttest2(IFD_avg, m_IFD_avg);
fprintf('\n Male jump distance p-value: %2.4g',p)

[~, p] = ttest2(IFD_avg, f_IFD_avg);
fprintf('\n Female jump distance p-value: %2.4g',p)



%% Rate of jumping in food and starved flies

