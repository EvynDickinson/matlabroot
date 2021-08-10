clear
close all
clc  
 
file_root = 'C:\matlabroot';
fig_dir = [file_root '/Interneuron Lines/Figures/'];

% select the test data
STIM = loadStructure;
disp(STIM.structure_name)
% select the appropriate control data
CONTROL = loadStructure;
disp(CONTROL.structure_name)

initial_vars = {'initial_vars','CONTROL','file_root', 'fig_dir', 'STIM', 'Fig_8'};
clearvars('-except',initial_vars{:})


%% Build general parameters
condLength = [0, 0.03, 0.06, 0.09, 0.180, 0.320, 0.720];
type = 'walking'; % stationary
fps = 30;
kolor = [];

for dd = 1:2
   switch dd
       case 1 %SH control
           fly = CONTROL.fly;
           group = CONTROL.group;
           num = CONTROL.num;
       case 2 %IN stim
           fly = STIM.fly;
           group = STIM.group;
           num = STIM.num;
   end
    for cond = 1:7 
      for ifly = 1:num.fly  
        % select the data and generate a filter for the initial behavior state
        filter = [group(ifly).(type)(cond,:), group(ifly).(type)(cond+7,:)];  
        cw = [fly(ifly).Control.speed(cond).data(1:end-1,:); fly(ifly).Stim.speed(cond).data];
        ccw = [fly(ifly).Control.speed(cond+7).data(1:end-1,:); fly(ifly).Stim.speed(cond+7).data];
        input = [cw, ccw];
        data(cond).ALL(ifly).data = input;
        input(:,~filter) = nan; %nix the trials that didn't fit the movement type
        % if there aren't 2 trials min, all go NaN
        if sum(filter)<2 %2 ADJUSTING THIS NUMBER FOR A TRIAL
            input(:,filter) = nan;
        end
        %stats on the fly's data
        data(cond).(type).raw(ifly).data = input;
        data(cond).(type).raw(ifly).avg = nanmean(input,2);
        data(cond).(type).raw(ifly).err = sem(input,2);
        % add averages to the group data
        data(cond).(type).avg(:,ifly) = data(cond).(type).raw(ifly).avg;
        data(cond).(type).err(:,ifly) = data(cond).(type).raw(ifly).err;
      end
      % group avg
      data(cond).(type).AVG = nanmean(data(cond).(type).avg,2);
      data(cond).(type).ERR = sem(data(cond).(type).avg,2);
    end
    switch dd
        case 1
            Cdata = data;
        case 2
            Sdata = data;
    end
end

% create a selection tool here for which data sets are being imported...

% comp = [STIM.structure_name(1:3)];

% -- gtacr1 -- 
% comp = '10B vs SH gtACR1';
% comp = '9A vs SH gtACR1';
% comp = '13B vs SH gtACR1';
% InColor = {'003400', '004D00', '006700', '178017', '319A31', '4AB34A', '64CD64'}; % gtacr1 greens

% -- chrimson --
% comp = 'SH gtACR vs SH Chrimson';
% InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'}; %10B oranges
% 

% comp = 'BDP vs SH Chrimson';
% InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'}; %10B oranges

% comp = 'BDP-gtACR1 vs SH Chrimson';
% InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'}; %10B oranges

% comp = 'SH-new vs SH-old Chrimson';
% InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'}; %10B oranges
% 
% comp = '10B vs SH Chrimson';
% InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'}; %10B oranges

comp = '9A vs SH Chrimson';
InColor = {'16005F', '490092', '6400c8', '7f00ff', '9a36ff', 'B66dfe', 'd1a3ff'}; % 9A purples

% comp = '9A-II vs 9A Chrimson';
% InColor = {'16005F', '490092', '6400c8', '7f00ff', '9a36ff', 'B66dfe', 'd1a3ff'}; % 9A purples
% comp = '9A-II vs SH Chrimson';
% InColor = {'16005F', '490092', '6400c8', '7f00ff', '9a36ff', 'B66dfe', 'd1a3ff'}; % 9A purples

% comp = '13B vs SH Chrimson';
% InColor = {'00005C', '003AA8', '006DDb', '198bff', '50a7ff', '8ac4ff', 'b6daff'}; % 13B blues


% InColor = {'', '', '', '', '', '', ''}; % 
for cond = 1:7
    kolor(2).K(cond,:) = hex2rgb(InColor{cond});
end
kolor(1).K = Color('black', 'lightgrey', 7);

%% Figure: Time course speed analysis comparing light lengths between SH & IN
ymax = 1.8;
type = 'walking'; % stationary
 x =  -2:1/30:2;
sSpan = 3;
for cond = [1,4,6,7] %1:7 %1:2:7
    fig = getfig; hold all
    for ii = 1:2 
        y = []; 
        switch ii
            case 1 % control
                 y.avg = smooth(Cdata(cond).(type).AVG,sSpan);
                 y.err = smooth(Cdata(cond).(type).ERR,sSpan);
                 idx = 1;
            case 2 % stim (10b, etc)
                y.avg = smooth(Sdata(cond).(type).AVG,sSpan);
                y.err = smooth(Sdata(cond).(type).ERR,sSpan);
                idx = 4;
        end
        colorList = kolor(ii).K;
%         fill_data = error_fill(x,  y.avg, y.err);
%         h = fill(fill_data.X, fill_data.Y, colorList(cond,:), 'EdgeColor','none');
%         set(h, 'facealpha', 0.2)
        
        plot(x, y.avg, 'color', colorList(idx,:), 'linewidth', 2)
        plot(x, y.avg+y.err, 'color', colorList(idx,:), 'linewidth', 0.5)
        plot(x, y.avg-y.err, 'color', colorList(idx,:), 'linewidth', 0.5)
    end
    %axes & labels
    
    offset = .03;
    ylim([0 ymax])
    vline([0, condLength(cond)], 'k')
    xlabel('time (sec)')
    ylabel('Speed (cm/s)')
    title([comp ' headless ' type ' cond: ' num2str(cond)])
    set(gca,'TickDir','out');
    xlim([-1,1])
    figname = [fig_dir, comp ' speed timecourse cond ' num2str(cond)];
    % savefig(fig, figname);
    save_figure(fig, figname); 
end


clear y idx ymax offset sSpan

%% Time course figures for each line (separately)
ymax = 1.8;

for dd = 1:2
    switch dd
       case 1 %SH control
           structure_name = CONTROL.structure_name;
           data = Cdata;
       case 2 %IN stim
           structure_name = STIM.structure_name;
           data = Sdata;
   end
   colorList = kolor(dd).K;     
    
    % PLOT THE DATA %
    y = [];
    condLength = [0, 0.03, 0.06, 0.09, 0.180, 0.320, 0.720];
    x =  -2:1/30:2;
    sSpan = 3;
    fig = getfig;
    hold all
    %CONTROL
    y.avg = smooth(data(1).(type).AVG,sSpan);
    y.err = smooth(data(1).(type).ERR,sSpan);
    % fill_data = error_fill(x,  y.avg, y.err);
    % h = fill(fill_data.X, fill_data.Y, Color('Black'), 'EdgeColor','none');
    % set(h, 'facealpha', 0.2)
    plot(x, y.avg, 'color', Color('Black'), 'linewidth', 1)
    plot(x, y.avg+y.err, 'color', Color('Black'), 'linewidth', 0.5)
    plot(x, y.avg-y.err, 'color', Color('Black'), 'linewidth', 0.5)
    %STIM % 60, 180, 720 = 3,5,7;
    for cond = [1,4,5,6] %[3,5,7]  %
        y.avg = smooth(data(cond).(type).AVG,sSpan);
        y.err = smooth(data(cond).(type).ERR,sSpan);
    %     fill_data = error_fill(x, y.avg, y.err);
    %     h = fill(fill_data.X, fill_data.Y, colorList(cond,:), 'EdgeColor','none');
    %     set(h, 'facealpha', 0.2)
        plot(x, y.avg, 'color', colorList(cond,:), 'linewidth', 1)
        plot(x, y.avg+y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
        plot(x, y.avg-y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
        
        offset = .03;
        ylim([0 ymax])
        LX = [0,condLength(cond)];
        LY = [ymax-(offset*cond), ymax-(offset*cond)];
        plot(LX, LY, 'Color', colorList(cond,:), 'linewidth', 3)

    end
    %axes & labels
   
    % vline([0, 0, 0.6, 0.18, 0.72], 'g')
    xlabel('time (sec)')
    ylabel('Speed (cm/s)')
    title([structure_name ' ' type])
    set(gca,'TickDir','out');
    % save the data
    figname = [fig_dir, structure_name ' ' type ' speed timecourse'];
    % savefig(fig, figname);
    save_figure(fig, figname); 
end

clear cw ccw ans colorList cond condFrames controlRange data dd fig figname
clear fill_data filter fly group h ifly ii inColor input LX LY num offset sSpan 
clear structure_name windowEnd y ymax stimRange stppoint

%% Time course -- multiple conditions overlaid from SH and other
ymax = 2;

fig = getfig;
hold all

for dd = 1:2
    switch dd
       case 1 %SH control
           structure_name = CONTROL.structure_name;
           data = Cdata;
       case 2 %IN stim
           structure_name = STIM.structure_name;
           data = Sdata;
   end
   colorList = kolor(dd).K;
    % PLOT THE DATA %
    y = [];
    condLength = [0, 0.03, 0.06, 0.09, 0.180, 0.320, 0.720];
    x =  -2:1/30:2;
    sSpan = 3;

    offset = .03;
    ylim([0 ymax])
    xlim([-2,2])
    %STIM % 60, 180, 720 = 3,5,7;
    for cond = 1:7  %[3,5,7]
        y.avg = smooth(data(cond).(type).AVG,sSpan);
        plot(x, y.avg, 'color', colorList(cond,:), 'linewidth', 1)
%         y.err = smooth(data(cond).(type).ERR,sSpan);
%         plot(x, y.avg+y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
%         plot(x, y.avg-y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
        LX = [0,condLength(cond)];
        if dd == 1
            LY = [ymax-(offset*cond), ymax-(offset*cond)];
        else
            LY = [ymax-((offset*cond))-(offset/4), ymax-(offset*cond)-(offset/4)];
        end
        plot(LX, LY, 'color', colorList(cond,:))
    end
end
    % vline([0, 0, 0.6, 0.18, 0.72], 'g')
    xlabel('time (sec)')
    ylabel('Speed (cm/s)')
    title([comp ' ' type])
    set(gca,'TickDir','out');
    % save the data
    xlim([-0.5,1]);% ylim([0,1.6])
    figname = [fig_dir, structure_name ' ' type ' overlaid speed timecourse'];
    % savefig(fig, figname);
    save_figure(fig, figname); 

clear cw ccw ans colorList cond condFrames controlRange data dd fig figname
clear fill_data filter fly group h ifly ii inColor input LX LY num offset sSpan 
clear structure_name windowEnd y ymax stimRange stppoint

%% Bootstrapping statistics for the timecourse data
% ---------------------------------------------------------------------
% 1) calculate the difference in speed from stimulus onset to 200ms post
% stimulus offset (diff in area under the curv) between SH and IN
% ---------------------------------------------------------------------
% compare between each light length (0ms, 90ms, 720ms) in the split half
% and the interneuron (eg 10B, 13B or 9A)


% ---------------------------------------------------------------------
% 2) Mix and randomly pull (10,000 times) samples of the trials (with
% replacement) and assign to control | stim and then find the diff in speed
% ---------------------------------------------------------------------
 
% gather the orginal running trials --not separated by fly
for cond = 1:7 % 90ms, 720ms conditions
    % pull the STIM fly averages
    input = [];
    for ifly = 1:STIM.num.fly
        input = [input, Sdata(cond).walking.raw(ifly).data];
    end
    filter = isnan(input(1,:)); % remove flies with no trials
    input(:,filter) = [];
    Sdata(cond).stats.raw = input;
    Sdata(cond).stats.num = size(input,2);
    % pull the CONTROL fly averages
    input = [];
    for ifly = 1:CONTROL.num.fly
        input = [input, Cdata(cond).walking.raw(ifly).data];
    end
    filter = isnan(input(1,:)); % remove flies with no trials
    input(:,filter) = [];
    Cdata(cond).stats.raw = input;
    Cdata(cond).stats.num = size(input,2);
end

% find the diff between the change in running speed in SH and IN
stats = [];
sROI = 66:72; % 200ms period 200ms post stim start
cROI = 55:61; % 200ms pre stimstart
for cond = 1:7
    stats(cond).param.stimROI = sROI;
    stats(cond).param.controlROI = cROI;
   % SPLIT HALF difference in speed:
    numflies = sum(~isnan(Cdata(cond).walking.avg(1,:)));
    a = Cdata(cond).stats.raw; % all trials for this condition
    ROI = a(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    a(:,filter) = []; %remove non walkers
    [diff, err]= statsChange(a,cROI,sROI);
    % 'save' the data into the stats struct
    stats(cond).OG.SH_diff = diff;
    stats(cond).OG.SH_err = err;
    
    % INTERNEURON difference in speed:
    numflies = sum(~isnan(Sdata(cond).walking.avg(1,:)));
    a = Sdata(cond).stats.raw; % all trials for this condition
    ROI = a(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    a(:,filter) = []; %remove non walkers
    % calc the change in speed for that condition:
    [diff, err]= statsChange(a,cROI,sROI);
    % 'save' the data into the stats struct
    stats(cond).OG.IN_diff = diff;
    stats(cond).OG.IN_err = err;
    % find the difference between the two changes in speed:
    stats(cond).OG.diff = stats(cond).OG.IN_diff-stats(cond).OG.SH_diff;
end

% plot the change in value for all conditions:
sz = 50; offset = 0.2;
fig = getfig; idx = 0; hold all
for cond = 1:7
    idx = idx+1;
    % pull up the diff + err values for both SH and IN:
    SH_err = stats(cond).OG.SH_err;
    IN_err = stats(cond).OG.IN_err;
    SH_diff = stats(cond).OG.SH_diff;
    IN_diff = stats(cond).OG.IN_diff;
%     % find the err:
%     err = sqrt(SH_err^2+IN_err^2)/sqrt(min([CONTROL.num.fly, STIM.num.fly]));
%   % Plot the SH:
    scatter(idx, SH_diff, sz, kolor(1).K(1,:))
    errorbar(idx, SH_diff, SH_err, 'color', kolor(1).K(1,:))
    % Plot the IN: 
    scatter(idx+offset, IN_diff, sz, kolor(2).K(4,:))
    errorbar(idx+offset, IN_diff, IN_err, 'color', kolor(2).K(4,:)) 
end
xlim([0, idx+1])
set(gca,'TickDir','out');
xlabel('Condition')
ylabel('Difference in speed between SH & IN (cm/s)')
figname = [comp ' diff in speed change  between SH + IN'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);


% bootstrap the trials and make a distribution of possibilities
tic
fig = getfig;
N = 10E4;
for cond = 1:7
    % SPLIT HALF difference in speed:
    SH_flies = sum(~isnan(Cdata(cond).walking.avg(1,:)));
    a = Cdata(cond).stats.raw; % all trials for this condition
    ROI = a(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    a(:,filter) = []; %remove non walkers
    
    IN_flies = sum(~isnan(Sdata(cond).walking.avg(1,:)));
    b = Sdata(cond).stats.raw; % all trials for this condition
    ROI = b(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    b(:,filter) = []; %remove non walkers
    
    mixed_data = [a,b];
     test = [];
  for n = 1:N

    % draw 'new' data:
    randLoc = randperm(SH_flies+IN_flies);
    SH_loc = randLoc(1:SH_flies);
    IN_loc = randLoc(SH_flies+1:end);
    
    % calc the change in speed for SH:
    a = [];
    a = mixed_data(:,SH_loc);
    diff = statsChange(a,cROI,sROI);
    test.SH(n) = diff;
    
    % calc the change in speed for IN:
    a = [];
    a = mixed_data(:,IN_loc);
    diff = statsChange(a,cROI,sROI);
    test.IN(n) = diff;
    % calc the diff between SH and IN:
    test.diff(n) = test.IN(n)-test.SH(n);
  end
  % 'save' the data into the test struct
    stats(cond).distrb = test;
    rdistrb = test.diff;
    subplot(2,4,cond); hold all
    histogram(rdistrb)
    vline(stats(cond).OG.diff, 'r-')
    p = sum(abs(rdistrb)>=abs(stats(cond).OG.diff))/length(rdistrb);
    disp(p);
    title({['cond: ' num2str(cond)]; ['p = ' num2str(p)]})
    stats(cond).p = p;
end
toc
subplot(2,4,8)
title(comp)
save_figure(fig, [fig_dir, comp, ' 10E4 bootstrap distrb change in speed']);
save([fig_dir, comp, ' 10E4 bootstrap distrb change in speed'], 'stats')



% multiple comparisons test:
conds = [1:7]; %[1,6,7];% conditions to look at
idx = 0;
for cond = conds
    idx = idx+1;
    p_val(idx) = stats(cond).p;
end
 p_err = p_val;
 [P_err,loc] = sort(p_err);
 fprintf('\n All:')
for idx = 1:length(p_val)
    q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
    r = (P_err((idx)) > (idx/length(P_err))*.05);
    fprintf(['\nInsignificant change: ' num2str(q) ' vs ' num2str(r) ' Cond: ' num2str(conds(loc(idx)))])
end
fprintf('\n Done\n')


% find the number of trials:
conds = [1:7];% conditions to look at
fprintf('\n N''s by fly')
for cond = conds
   % find num trials
   
    SH_flies = sum(~isnan(Cdata(cond).walking.avg(1,:)));
    IN_flies = sum(~isnan(Sdata(cond).walking.avg(1,:)));
    fprintf(['\n COND: ' num2str(cond) ' SH: ' num2str(SH_flies) ' IN: ' num2str(IN_flies)])
end
fprintf('\n N''s by trial')
for cond = conds
    a = []; b = [];
    for ifly = 1:CONTROL.num.fly
        a = [a, Cdata(cond).walking.raw(ifly).data(1,:)];
    end 
    for ifly = 1:STIM.num.fly
        b = [b, Sdata(cond).walking.raw(ifly).data(1,:)];
    end
    SH_flies = sum(~isnan(a));
    IN_flies = sum(~isnan(b));
    fprintf(['\n COND: ' num2str(cond) ' SH: ' num2str(SH_flies) ' IN: ' num2str(IN_flies)])
end
fprintf('\n P''s by trial')
for cond = conds
    a = []; b = [];
    for ifly = 1:CONTROL.num.fly
        a = [a, Cdata(cond).walking.raw(ifly).data(1,:)];
    end 
    for ifly = 1:STIM.num.fly
        b = [b, Sdata(cond).walking.raw(ifly).data(1,:)];
    end
    SH_flies = sum(~isnan(a));
    IN_flies = sum(~isnan(b));
    fprintf(['\n COND: ' num2str(cond) ' SH: ' num2str(SH_flies) ' IN: ' num2str(IN_flies)])
end

% % change in speed comparison (plot)

% fig = figure; hold all
% sz = 50;
% for cond = 1:7
%    scatter(cond, stats(cond).OG.IN_diff, sz, kolor(2).K(3,:), 'filled')
%    errorbar(cond, stats(cond).OG.IN_diff,stats(cond).OG.IN_err, 'color', kolor(2).K(3,:))
%    scatter(cond, stats(cond).OG.SH_diff, sz, 'k', 'filled')
%    errorbar(cond, stats(cond).OG.SH_diff,stats(cond).OG.IN_err, 'color','k')
% end
% xlim([0,8])
% set(gca,'TickDir','out');
% xlabel('Condition')
% ylabel('Difference in speed change (cm/s)')
% figname = [comp ' diff in speed change'];
% title({figname;' ';' '})
% 
% save_figure(fig, [fig_dir, figname]);



% 
% % close all
% for cond = 1:7
% sROI = 66:72;
% cROI = 55:61;
% figure; subplot(1,2,1); hold all
%     a = Cdata(cond).stats.raw;
%     ROI = a(30:60,:);
%     filter = median(ROI)<=0.3;
%     a(:,filter) = [];
%     plot(a, 'color', 'k')
%     strt = mean(mean(a(cROI,:)));
%     stp = mean(mean(a(sROI,:)));
%     hline(stp, 'r-')
%     hline(strt, 'k-')
%     title({'SH control'; ['diff = ' num2str(strt-stp)]})
%     subplot(1,2,2); hold all
%     a = Sdata(cond).stats.raw;
%     filter = median(a(30:60,:))<=0.3;
%     a(:,filter) = [];
%     plot(a, 'color', 'b')
%     strt = median(mean(a(cROI,:)));
%     stp = median(mean(a(sROI,:)));
%     hline(stp, 'r-')
%     hline(strt, 'k-')
%     title({['IN stimulus cond: ' num2str(cond)]; ['diff = ' num2str(strt-stp)]})
% end 


clear a alltrials ans b cond conds cROI diff duration err fig filter
clear idx ii ifly IN_flies IN_loc INavg INcntl INnum INraw INstim loc
clear mixed_data loc n N numflies p p_err P_err p_val post post_err pre
clear pre_err pval q r randLoc raw rdistrb ROI SH_flies SH_loc SHavg
clear SHcntl SHnum SHstim sROI stp strt sz tEnd tRange x Z


%% Cumulative probability density function
type = 'walking';
min_speed = 0.3;

% percent freezing
stppoint = 0.500;
windowEnd = (60+stppoint*fps)*ones(1,7);

sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for dd = 1:2
    switch dd
        case 1 %SH control
            data = Cdata;
            group = CONTROL.group;
            num = CONTROL.num;
        case 2 %INT stim
            data = Sdata;
            group = STIM.group;
            num = STIM.num;
    end
    for cond = 1:7
        for ifly = 1:num.fly
            for rep = 1:6
                vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                      'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
                raw = data(cond).(type).raw(ifly).data(:,rep); 
                % if the data fits the behavior type, analyze
                if ~isnan(raw(1))
                    temp = smooth(raw, sSpan, 'rlowess');
                    vStart = nanmean(raw(54:60)); %value prestim
    %                 vChange = NaN;
                    % only include trials with greater than min speed
                    if vStart >= min_speed
                        % Change in speed
                        vStim = nanmean(raw(76:82)); %value prestim
                        vChange = vStim-vStart;
                        vFreeze(rep).speedchange = vChange;
                        %find locations that the fly isn't moving post stim start
                        vLow = find(temp(61:121) <= min_speed); 
                        if ~isempty(vLow)
                            vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                            % calculate the specific bits that qualify for 'freezing'
                            vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                            vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                            % find the start and end points of the 'freeze'
                            try % gaps present
                                vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                                vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                                vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                                if vFreeze(rep).start>vFreeze(rep).end
                                    vFreeze(rep).start = vFreeze(rep).end;
                                end
                            catch % no gaps in the 'frozen' data points
                                vFreeze(rep).start = 1;
                                vFreeze(rep).end = length(vFreeze(rep).pos);
                            end
                            vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                            vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                            if vFreeze(rep).idx(1)>windowEnd(cond)
                                vFreeze(rep).idx = [];
                                vFreeze(rep).start = [];
                                vFreeze(rep).end = [];
                            end 
                        end
                    end
                end
            end
            tot = (sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]));
            switch dd
                case 1
                    FreezeData(cond).CONTROL(ifly).data = vFreeze;
                    FreezeData(cond).CONTROL(ifly).totwalking = tot;
                case 2
                    FreezeData(cond).STIM(ifly).data = vFreeze;
                    FreezeData(cond).STIM(ifly).totwalking = tot;
            end
        end
    end
end

% Calculate the probabilty of a fly stopping...
MT = max([CONTROL.num.fly, STIM.num.fly]);
MT = NaN(MT,1);

for cond = 1:7
    for dd = 1:2
        temp = []; output = [];
        switch dd
            case 1
                temp = FreezeData(cond).CONTROL;
                num = CONTROL.num;
                output.raw = zeros(6,num.fly);
            case 2
                temp = FreezeData(cond).STIM;
                num = STIM.num;
                output.raw = zeros(6,num.fly);
        end
        for ifly = 1:num.fly
            for rep = 1:6
                if ~isempty(temp(ifly).data(rep).duration)
                    % find the start place
                    output.raw(rep,ifly) = temp(ifly).data(rep).pos(1);
                end
            end
            output.totwalking(ifly) = temp(ifly).totwalking;
        end
        % pull together data across flies
        output.totalpeak = sum(output.totwalking);
        a = (output.raw>0); % find all instances of freezing
        output.totalfreeze =  sum(sum(a));
        output.maxfreeze = output.totalfreeze/output.totalpeak;
        output.inc = output.maxfreeze/output.totalfreeze;
        output.loc = output.raw(a);
        
        % find the probability
        a = sort(output.loc);
        a = a-60;
        x = []; y = [];
        for ii = 1:length(a)
           x(ii) = a(ii); 
           y(ii) = ii*output.inc;
        end
        output.x = [0, x, 60];
        output.y = [0, y, output.maxfreeze];
        switch dd
            case 1
                cntl = output;           
            case 2
                stim = output;
        end
    end  
    %     fig = getfig;
    %     hold all
    %     plot(cntl.x, cntl.y, 'color', 'k')
    %     plot(stim.x, stim.y, 'color', 'r')
    FreezeData(cond).C_output = cntl;
    FreezeData(cond).S_output = stim;
end
   
% color selection:
controlColor = kolor(1).K;
stimColor = kolor(2).K;

% Plot the data together:
LW = 2;
fig = getfig; hold all

for cond = [1, 4, 6, 7]
    plot(FreezeData(cond).C_output.x, FreezeData(cond).C_output.y, 'color', controlColor(cond,:), 'linewidth', LW);
    plot(FreezeData(cond).S_output.x, FreezeData(cond).S_output.y, 'color', stimColor(cond,:), 'linewidth', LW);
end
xlim([0, 60])
ylim([0, 1])
set(gca,'TickDir','out');
xlabel('Time (fps)')
ylabel('Percent of freezes (%)')
figname = [comp ' freezing probabilty curve'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);

clear a cntl controlColor issort ix loc LW MT ouput raw sloc temp tmax 
clear tmin tot vChange vFreeze vLow vStart vStim x1 xq y1 yq xq1 stimColor
clear ans ifly ii cond dd rep sSpan x y windowEnd isort stppoint ColorList InColor


%% Heatmap of walking speed for only WALKING trials : max speed normalized
clear avg avg_speed c Ccombo clims cw err f fig figname filter MT
clear framelength h idx input L MED MT MX nbins num output plotdata positiveIdx
clear Scombo sitm sz Tduration total Trange tt xi xx Y YY ALL All_sorted
xlimits = [45, 90];
timemarks = 60+round(condLength*fps);

Scombo = []; 
MT = NaN(10,121);
%Stim data
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:STIM.num.fly
       a = Sdata(cond).walking.raw(ifly).data;
       loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Scombo = [Scombo; MT;All_sorted]; %
end
ALL = [];
% Z SCORE for the max speed on each trial
Z = zscore(Scombo,0,2);
peakvalue = max(Z,[],2);
ALL = Z./peakvalue;
normalizedZscore = ALL;

% normalized graph
ALL = [];
peakvalue = max(Scombo,[],2);
ALL = Scombo./peakvalue;
normalizedSpeed = ALL;

% PLOT THE interneuron DATA
fig = getfig; 
subplot(2,1,1)
    imagesc(normalizedSpeed)%(:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Normalized Speed (cm/s)';
%     xlabel('time')
    ylabel('Trials by Cond')  
    title([STIM.structure_name ' normalized speed'])
    xlim(xlimits)
    vline(timemarks, 'g')
    set(gca,'xtick',[])
    
subplot(2,1,2)
    imagesc(normalizedZscore)% (:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
    xlabel('time')
    ylabel('Trials by Cond')  
    title([STIM.structure_name ' normalized Z-Score'])
    xlim(xlimits)
    vline(timemarks, 'g')
    set(gca,'xtick',[])
    
figname  = [fig_dir, STIM.structure_name ' WALKING normalized heatmaps'];
save_figure(fig, figname);

% control data
Ccombo = [];
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:CONTROL.num.fly
       a = Cdata(cond).walking.raw(ifly).data;
       loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Ccombo = [Ccombo;MT; All_sorted]; %
end
ALL = [];
% Z SCORE for the max speed on each trial
Z = zscore(Ccombo,0,2);
peakvalue = max(Z,[],2);
ALL = Z./peakvalue;
normalizedZscore = ALL;

% normalized graph
ALL = [];
peakvalue = max(Ccombo,[],2);
ALL = Ccombo./peakvalue;
normalizedSpeed = ALL;

% PLOT THE interneuron DATA
fig = getfig; 
subplot(2,1,1)
    imagesc(normalizedSpeed)%(:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Normalized Speed (cm/s)';
%     xlabel('time')
    ylabel('Trials by Cond')  
    title([CONTROL.structure_name ' normalized speed'])
    xlim(xlimits)
    vline(timemarks, 'g')
    set(gca,'xtick',[])
subplot(2,1,2)
    imagesc(normalizedZscore)% (:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
    xlabel('time')
    ylabel('Trials by Cond')  
    title([CONTROL.structure_name ' normalized Z-Score'])
    xlim(xlimits)
    vline(timemarks, 'g')
    set(gca,'xtick',[])
    
figname  = [fig_dir, CONTROL.structure_name ' WALKING normalized heatmaps'];
save_figure(fig, figname);



%% JUMP COUNT:

% % Convert behaviors into numbers and then a filter:
% for kk = 1:length(group)
%    state = strcmpi(group(kk).behavior, '-');
%    output(kk).jump = ~state;  
%    output(kk).totalJumps = sum(sum(~state));
% end
% save('SH-gal4xUAS-csChrimson jump count', 'output', 'group');
% save('10B-04751-gal4xUAS-csChrimson jump count', 'output', 'group');

[IN_jump] = load('10B-04751-gal4xUAS-csChrimson jump count');
[SH_jump] = load('SH-gal4xUAS-csChrimson jump count');

% plot each number of jumps for both lines:

fig = figure; set(fig, 'color', 'w');
subplot(1,2,1)
hold all
num = length(SH_jump.output);
sz = 50; offset = 1/20;
for ifly = 1:num
    SH(ifly) = SH_jump.output(ifly).totalJumps;
    scatter(1+(ifly*offset), SH(ifly), sz, 'k')
end
plot([1,2], [mean(SH), mean(SH)], 'color', 'k')
num = length(IN_jump.output);
for ifly = 1:num
    IN(ifly) = IN_jump.output(ifly).totalJumps;
    scatter(3+(ifly*offset), IN(ifly), sz, 'k')
end
plot([3,4], [mean(IN), mean(IN)], 'color', 'k')
xlim([0,4])

p = ranksum(SH, IN);
title({''; '10B vs SH x Chrimson'; ['p = ' num2str(p)]})

set(gca,'TickDir','out');
xlabel('SH  vs   10B')
ylabel('Jumps (#)')

subplot(1,2,2)
hold all
num = length(SH_jump.output);
sz = 50; offset = 1/20;
for ifly = 1:num
    SH(ifly) = SH_jump.output(ifly).totalJumps/84;
    scatter(1+(ifly*offset), SH(ifly), sz, 'k')
end
plot([1,2], [mean(SH), mean(SH)], 'color', 'k')
num = length(IN_jump.output);
for ifly = 1:num
    IN(ifly) = IN_jump.output(ifly).totalJumps/84;
    scatter(3+(ifly*offset), IN(ifly), sz, 'k')
end
plot([3,4], [mean(IN), mean(IN)], 'color', 'k')   
xlim([0,4])
  
set(gca,'TickDir','out');
xlabel('SH  vs   10B')
ylabel('Percent of trials with jumps')

p = ranksum(SH, IN);
title({''; '10B vs SH x Chrimson'; ['p = ' num2str(p)]})
 
file_root = 'C:\matlabroot';
fig_dir = [file_root '/Interneuron Lines/Figures/'];
save_figure(fig, [fig_dir, '10B vs SH Jump comparison']);

% STATS:
% wilcox rank sum test




%% Unused sections below



%% Probability density function
type = 'walking';
min_speed = 0.3;

% percent freezing
stppoint = 0.500;
windowEnd = (60+stppoint*fps)*ones(1,7);

sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for dd = 1:2
    switch dd
        case 1 %SH control
            data = Cdata;
            group = CONTROL.group;
            num = CONTROL.num;
        case 2 %INT stim
            data = Sdata;
            group = STIM.group;
            num = STIM.num;
    end
    for cond = 1:7
        for ifly = 1:num.fly
            for rep = 1:6
                vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                      'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
                raw = data(cond).(type).raw(ifly).data(:,rep); 
                % if the data fits the behavior type, analyze
                if ~isnan(raw(1))
                    temp = smooth(raw, sSpan, 'rlowess');
                    vStart = nanmean(raw(54:60)); %value prestim
    %                 vChange = NaN;
                    % only include trials with greater than min speed
                    if vStart >= min_speed
                        % Change in speed
                        vStim = nanmean(raw(76:82)); %value prestim
                        vChange = vStim-vStart;
                        vFreeze(rep).speedchange = vChange;
                        %find locations that the fly isn't moving post stim start
                        vLow = find(temp(61:121) <= min_speed); 
                        if ~isempty(vLow)
                            vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                            % calculate the specific bits that qualify for 'freezing'
                            vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                            vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                            % find the start and end points of the 'freeze'
                            try % gaps present
                                vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                                vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                                vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                                if vFreeze(rep).start>vFreeze(rep).end
                                    vFreeze(rep).start = vFreeze(rep).end;
                                end
                            catch % no gaps in the 'frozen' data points
                                vFreeze(rep).start = 1;
                                vFreeze(rep).end = length(vFreeze(rep).pos);
                            end
                            vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                            vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                            if vFreeze(rep).idx(1)>windowEnd(cond)
                                vFreeze(rep).idx = [];
                                vFreeze(rep).start = [];
                                vFreeze(rep).end = [];
                            end 
                        end
                    end
                end
            end
            switch dd
                case 1
                    FreezeData(cond).CONTROL(ifly).data = vFreeze;
                case 2
                    FreezeData(cond).STIM(ifly).data = vFreeze;
            end
        end
    end
end

% Calculate the probabilty of a fly stopping...
MT = max([CONTROL.num.fly, STIM.num.fly]);
MT = NaN(MT,1);
for cond = 1:7
    for dd = 1:2
        temp = []; output = [];
        switch dd
            case 1
                temp = FreezeData(cond).CONTROL;
                num = CONTROL.num;
                output.raw = zeros(6,num.fly);
            case 2
                temp = FreezeData(cond).STIM;
                num = STIM.num;
                output.raw = zeros(6,num.fly);
        end
        for ifly = 1:num.fly
            for rep = 1:6
                if ~isempty(temp(ifly).data(rep).duration)
                    % find the start place
                    output.raw(rep,ifly) = temp(ifly).data(rep).pos(1);
                end
            end
        end
        % pull together data across flies
        a = (output.raw>0); % find all instances of freezing
        output.total = sum(sum(a));
        output.loc = output.raw(a);
        output.inc = 1/output.total;
        % find the probability
        a = sort(output.loc);
        a = a-60;
        x = []; y = [];
        for ii = 1:length(a)
           x(ii) = a(ii); 
           y(ii) = ii*output.inc;
        end
        output.x = [0, x, 60];
        output.y = [0, y, 1];
        switch dd
            case 1
                cntl = output;           
            case 2
                stim = output;
        end
    end  
    %     fig = getfig;
    %     hold all
    %     plot(cntl.x, cntl.y, 'color', 'k')
    %     plot(stim.x, stim.y, 'color', 'r')
    FreezeData(cond).C_output = cntl;
    FreezeData(cond).S_output = stim;
end
  
% color selection:
controlColor = kolor(1).K;
stimColor = kolor(2).K;

% Plot the data together:
LW = 2;
fig = getfig; hold all

for cond = [1, 4, 6, 7]
    plot(FreezeData(cond).C_output.x, FreezeData(cond).C_output.y, 'color', controlColor(cond,:), 'linewidth', LW);
    plot(FreezeData(cond).S_output.x, FreezeData(cond).S_output.y, 'color', stimColor(cond,:), 'linewidth', LW);
end
xlim([0, 61])

set(gca,'TickDir','out');
xlabel('Time (fps)')
ylabel('Percent of freezes (%)')
figname = [comp ' freezing probabilty curve'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);

clear a cntl controlColor issort ix loc LW MT ouput raw sloc temp tmax 
clear tmin tot vChange vFreeze vLow vStart vStim x1 xq y1 yq xq1 stimColor
clear ans ifly ii cond dd rep sSpan x y windowEnd isort stppoint ColorList InColor


%% heatmap of ALL flies' speed -- not just the walking flies

Scombo = []; 
MT = NaN(3,121);
%Stim data
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:STIM.num.fly
       a = Sdata(cond).ALL(ifly).data;
%        loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Scombo = [Scombo; MT;All_sorted];% 
end

ALL = [];
% Z SCORE for the max speed on each trial
Z = zscore(Scombo,0,2);
peakvalue = max(Z,[],2);
ALL = Z./peakvalue;
normalizedZscore = ALL;

% normalized graph
ALL = [];
peakvalue = max(Scombo,[],2);
ALL = Scombo./peakvalue;
normalizedSpeed = ALL;

% PLOT THE interneuron DATA
fig = getfig; 
subplot(2,1,1)
    imagesc(normalizedSpeed)%(:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    plasma = plasma();
    colormap(plasma);
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Normalized Speed (cm/s)';
%     xlabel('time')
    ylabel('Trials by Cond')  
    title([STIM.structure_name ' normalized speed'])
    xlim(xlimits)
subplot(2,1,2)
    imagesc(normalizedZscore)% (:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    plasma = plasma();
    colormap(plasma);
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
    xlabel('time')
    ylabel('Trials by Cond')  
    title([STIM.structure_name ' normalized Z-Score'])
    xlim(xlimits)
    
figname  = [fig_dir, STIM.structure_name ' ALL TRIALS normalized heatmaps'];
save_figure(fig, figname);


% control data
Ccombo = [];
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:CONTROL.num.fly
       a = Cdata(cond).ALL(ifly).data;
%        loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Ccombo = [Ccombo; MT;All_sorted];% 
end

ALL = [];
% Z SCORE for the max speed on each trial
Z = zscore(Ccombo,0,2);
peakvalue = max(Z,[],2);
ALL = Z./peakvalue;
normalizedZscore = ALL;

% normalized graph
ALL = [];
peakvalue = max(Ccombo,[],2);
ALL = Ccombo./peakvalue;
normalizedSpeed = ALL;

% PLOT THE interneuron DATA
fig = getfig; 
subplot(2,1,1)
    imagesc(normalizedSpeed)%(:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    plasma = plasma();
    colormap(plasma);
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Normalized Speed (cm/s)';
%     xlabel('time')
    ylabel('Trials by Cond')  
    title([CONTROL.structure_name ' normalized speed'])
    xlim(xlimits)
subplot(2,1,2)
    imagesc(normalizedZscore)% (:,30:90)
    clims =  [0, 1];
    caxis(clims);
    hold all
    plasma = plasma();
    colormap(plasma);
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
    xlabel('time')
    ylabel('Trials by Cond')  
    title([CONTROL.structure_name ' normalized Z-Score'])
    xlim(xlimits)
    
figname  = [fig_dir, CONTROL.structure_name ' ALL TRIALS normalized heatmaps'];
save_figure(fig, figname);



%% Avg delay to freeze for all flies  (no time limitations)
% Uses the FreezeData struct created for cumulative distribution

plotdata = [];
Tduration = 300; % how long after stim onset to look
Trange = 61+round(Tduration/fps);
TYPE = {'CONTROL', 'STIM'};
num = [CONTROL.num.fly, STIM.num.fly];
for tt = 1:2
    output = [];
    for cond = 1:7 
        for ifly = 1:num(tt)
            total = FreezeData(cond).(TYPE{tt})(ifly).totwalking;
            if total > 0 %if there are walking trials
                for rep = 1:6
                    try temp(rep) = FreezeData(cond).(TYPE{tt})(ifly).data(rep).pos(1);
                    catch temp(rep) = NaN; end
                end
                % does the pause happen within the time frame?
                positiveIdx = sum(temp < Trange);
                output(cond).frozen(ifly) = positiveIdx/total;
                output(cond).time(:,ifly) = temp';
                output(cond).walking(ifly) = true;
                output(cond).flyavg(ifly) = nanmean(temp);
            else
                output(cond).frozen(ifly) = NaN;
                output(cond).time(:,ifly) = NaN(6,1);
                output(cond).walking(ifly) = false;
                output(cond).flyavg(ifly) = NaN;
            end
        end
        loc = ~isnan(output(cond).time);
        a = output(cond).time(loc);
        output(cond).start = a;
    end
    plotdata.(TYPE{tt}) = output;
end


% MEAN TIME TO STOP DELAY -- all trials avg
fig = getfig; hold all
for cond = [1,4,6,7]
   x = condLength(cond);
   % control   
   total =  sum(plotdata.CONTROL(cond).walking); 
   a = plotdata.CONTROL(cond).start; 
   a = (a-60)/fps; % convert frames to seconds
   avg = mean(a);
   err = std(a)/sqrt(total); 
   
   scatter(x, avg, sz, kolor(1).K(1,:), 'filled')
   errorbar(x,avg,err, 'color', kolor(1).K(1,:))
   
   % stim
   total =  sum(plotdata.STIM(cond).walking); 
   a = plotdata.STIM(cond).start;
   a = (a-60)/fps; % convert frames to seconds
   avg = mean(a);  
   err = std(a)/sqrt(total); 
   
   scatter(x+.01, avg, sz, kolor(2).K(3,:), 'filled')
   errorbar(x+.01,avg,err, 'color', kolor(2).K(3,:))
end
xlim([-.10 0.8]); ylim([0, 2])


% MEAN TIME TO STOP DELAY - fly avg
% fig = getfig; hold all
% for cond = [1,4,6,7]
%    x = condLength(cond);
%    % control   
%    total =  sum(plotdata.CONTROL(cond).walking); 
%    a = plotdata.CONTROL(cond).flyavg; 
%    a = (a-60)/fps; % convert frames to seconds
%    avg = nanmean(a);
%    err = nanstd(a)/sqrt(total); 
%    
%    scatter(x, avg, sz, kolor(1).K(1,:), 'filled')
%    errorbar(x,avg,err, 'color', kolor(1).K(1,:))
%    
%    % stim
%    total =  sum(plotdata.STIM(cond).walking); 
%    a = plotdata.STIM(cond).flyavg;
%    a = (a-60)/fps; % convert frames to seconds
%    avg = nanmean(a);  
%    err = nanstd(a)/sqrt(total); 
%    
%    scatter(x+.01, avg, sz, kolor(2).K(3,:), 'filled')
%    errorbar(x+.01,avg,err, 'color', kolor(2).K(3,:))
% end
% xlim([-.10 0.8]); ylim([0, 2])

set(gca,'TickDir','out');
xlabel('stimulus duration (sec)')
ylabel('avg time to first freeze (s)')
figname = [comp ' avg delay to freeze'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);
save([fig_dir figname], 'plotdata')

%% Time to 50% total freezing
% Build the interpolated data line
for dd = 1:2
    for cond = 1:7
        switch dd
            case 1
                x = [FreezeData(cond).C_output.x(2:end), 61]; % 
                y = [FreezeData(cond).C_output.y(2:end), FreezeData(cond).C_output.maxfreeze];% 
            case 2
                x = [FreezeData(cond).S_output.x(2:end), 61]; % 
                y = [FreezeData(cond).S_output.y(2:end), FreezeData(cond).S_output.maxfreeze];% 
        end
        % remove repeat values (eg 60, 60; 0.6, 0.6)
        [fX,idx] = sort(x);
        Xloc = (diff(fX)==0);
        fY = y(idx);
        Yloc = (diff(fY)==0);
        repeatval = all([Yloc; Xloc]);
        x(repeatval) = []; y(repeatval) = [];
        xq = 1:61; xq(x) = [];
        ix = find(diff(x)==0); % where there are repeat numbers
        % run each section...
        sloc = 1; output = [];
        raw = [x',y'];
        for ii = 1:(length(ix)+1)
            try
                loc = (ix(ii));
            catch
                loc = length(x); 
            end
            x1 = x(sloc:loc);
            y1 = y(sloc:loc);
            tmin = xq>=x1(1); tmax = (xq<=x1(end));
            xq1 = xq(all([tmin; tmax]));
            yq = interp1(x1,y1,xq1);
            output(ii).yq = yq;
            output(ii).xq1 = xq1;
            sloc = loc+1;
            raw = [raw; output(ii).xq1', output(ii).yq'];
        end
        % sort the order 
        [~, isort] = sort(raw(:,1));
        temp = [];
        temp = raw(isort,:);
        % fill in the gaps from t=0 to first data point
        x1 = [0,temp(1,1)]; y1 = [0, (temp(1,2))]; 
        xq1 = (x1(1)+1):(x(1)-1);
        yq = interp1(x1,y1,xq1);
        raw = [xq1',yq';temp];
        sSpan = 6;
        y = smooth(raw(:,2),sSpan);
        x = raw(:,1);
%         midpoint = y(end)*.8;
        midpoint = y(end)/2;
        [~, midloc] = min(abs(y-midpoint));
%         figure; plot(x, y); hline(midpoint, 'r'); vline(midloc, 'r'); %visual check
        % Save the data back into the structure
        FreezeData(cond).data(dd).x = x;
        FreezeData(cond).data(dd).y = y;
        FreezeData(cond).data(dd).midpoint = midpoint;
        FreezeData(cond).data(dd).midloc = midloc;
    end
end       
        
% Compare the half freezing time points:
fig = getfig; hold all; ylim([0,35]); xlim([0, 0.8])
sz = 60;
for dd = 1:2
    temp = [];
    clr = kolor(dd).K(3,:);
    for cond = 1:7
        x = condLength(cond);
        y = FreezeData(cond).data(dd).midloc; %time length
        scatter(x,y, sz, clr, 'filled')
        temp(cond) = y;
    end
    plot(condLength, temp, 'color', clr, 'linewidth', 1)
end
set(gca,'TickDir','out');
xlabel('stimulus duration')
ylabel('Time to 50% freezes (fps)')
figname = [comp ' 50 percent freezing probabilty curve'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);

     

% Plotting 50% time vs freeze value
fig = getfig; hold all; %ylim([0,61]); xlim([0, 0.8])
sz = 60;
for dd = 1:2
    temp = [];
    clr = kolor(dd).K(3,:);
    for cond = 1:7
        x = FreezeData(cond).data(dd).midloc; %time length
        y = FreezeData(cond).data(dd).midpoint; %time length
        scatter(x,y, sz, clr, 'filled')
        temp.x(cond) = x;
        temp.y(cond) = y;
    end
    plot(temp.x, temp.y, 'color', clr, 'linewidth', 1)
end
set(gca,'TickDir','out');
xlabel('time to 50% freeze')
ylabel('percent of flies 50% freeze')
figname = [comp ' 50 percent freezing time vs total'];
title({figname;' ';' '})
save_figure(fig, [fig_dir, figname]);

% save the freezing data for comp across lines
save([comp ' FreezeData'], 'FreezeData')

clear x y minDistance loc maxpoint x1 y1 Xloc Yloc repeatval temp
clear ans C clr cond fig fX fY ia idx ii iy midloc midpoint output raw
clear sloc sSpan sz tmax tmin xq xq1 yq

%%  not yet combined for multiple data sets
% BUILD SOMETHING TO TEST...TIME TO 80% FROZEN? RATHER THAN DELAY TO
% FREEZE?



hlimit = 0.5;
itype = {'CONTROL', 'STIM'};
sz = 60;
MT = NaN(6,1);    
RNGE = 54:61;
fig = getfig; hold all; idx = 0;
num.fly = size(FreezeData(1).CONTROL,2);
for dd = 1:2 %SH vs IN
    if dd == 1
       data = Cdata;
       klr = 'k';
    else
        data = Sdata;
        klr = 'r';
    end
    for cond = 1:7 %[1, 4:6]  
        idx = idx+1;
        subplot(2,7,idx); hold all; ylim([0 2]); xlim([0 2.5])
        % load the delay and duration data into 'output'
        output = struct('delay', MT, 'duration', MT, 'speed', MT);
        for ifly = 1:num.fly
            inputdata = FreezeData(cond).(itype{dd})(ifly).data; 
            output.duration = MT; output.delay = MT; output.speed = MT;
            for rep = 1:6
                if ~isempty(inputdata(rep).duration)
                    output.duration(rep) = inputdata(rep).duration/fps;  
                    output.delay(rep) = (inputdata(rep).pos(1)-60)/fps;
                    % avg speed pre-stim onset
                    a = nanmean(data(cond).walking.raw(ifly).data(RNGE,rep));
                    output.speed(rep) = a;
                end
            end
            if sum(~isnan(output.speed))>1
                output.allduration(:,ifly) = output.duration;
                output.alldelay(:,ifly) = output.delay;
                output.allspeed(:,ifly) = output.speed;
            else
                output.allduration(:,ifly) = MT;
                output.alldelay(:,ifly) = MT;
                output.allspeed(:,ifly) = MT;
            end
        end
        % group averages for the delay and duration

        loc = ~isnan(output.allduration);
        %duration
        output.dur_list = output.allduration(loc);
        output.dur_avg = mean(output.dur_list);
        output.dur_err = nanstd(output.dur_list);
        %delay
        output.delay_list = output.alldelay(loc);
        output.delay_avg = mean(output.delay_list);
        output.delay_err = nanstd(output.delay_list);
        %speed
        output.speed_list = output.allspeed(loc);
        output.speed_avg = mean(output.speed_list);
        output.speed_err = nanstd(output.speed_list);
        % how many fall into the right delay group?
        TOTAL = length(output.delay_list);
        withinLim = sum(output.delay_list<=hlimit);
        percentLim = (withinLim/TOTAL)*100;
        title({[num2str(percentLim) '% under ' num2str(hlimit) 's'];''})

        % all points
        x = output.speed_list;
        y = output.delay_list;
        
        scatter(x,y,sz,kolor(dd).K(cond,:), 'filled')
%       
        h = lsline;
        h.Color = klr;
      
        set(gca,'TickDir','out');
        hline(condLength(cond), 'g')
        hline(hlimit, 'm')
        ylim([0 2]); xlim([0 2.5])
    end 
end
subplot(2,7,1)
xlabel('Speed at stim onset (cm/s)')
ylabel('Delay to freeze onset (s)')
figname  = [fig_dir, comp ' speed vs delay percentages'];
% title({[comp ' speed vs delay'];'';''})

save_figure(fig, figname);



%% Percent of walking flies that freeze (by fly) -- COMBINE FOR MULTIPLE DATA SETS
type = 'walking';
min_speed = 0.5;

% percent freezing
stppoint = 0.500;
windowEnd = (60+stppoint*fps)*ones(1,7);

condFrames = round(condLength*num.fps);
sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for dd = 1:2
    switch dd
        case 1 %SH control
            fly = CONTROL.fly;
            group = CONTROL.group;
            num = CONTROL.num;
        case 2 %IN stim
            fly = STIM.fly;
            group = STIM.group;
            num = STIM.num;
    end

    for cond = 1:7
        for ifly = 1:num.fly
            for rep = 1:6
                vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                      'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
                raw = data(cond).(type).raw(ifly).data(:,rep); 
                % if the data fits the behavior type, analyze
                if ~isnan(raw(1))
                    temp = smooth(raw, 'rlowess');
                    vStart = nanmean(raw(54:60)); %value prestim
    %                 vChange = NaN;

                    % only include trials with greater than min speed
                    if vStart >= min_speed
                        % Change in speed
                        vStim = nanmean(raw(76:82)); %value prestim
                        vChange = vStim-vStart;
                        vFreeze(rep).speedchange = vChange;
                        %find locations that the fly isn't moving post stim start
                        vLow = find(temp(61:121) <= min_speed); 
                        if ~isempty(vLow)
                            vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                            % calculate the specific bits that qualify for 'freezing'
                            vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                            vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                            % find the start and end points of the 'freeze'
                            try % gaps present
                                vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                                vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                                vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                                if vFreeze(rep).start>vFreeze(rep).end
                                    vFreeze(rep).start = vFreeze(rep).end;
                                end
                            catch % no gaps in the 'frozen' data points
                                vFreeze(rep).start = 1;
                                vFreeze(rep).end = length(vFreeze(rep).pos);
                            end
                            vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                            vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                            if vFreeze(rep).idx(1)>windowEnd(cond)
                                vFreeze(rep).idx = [];
                                vFreeze(rep).start = [];
                                vFreeze(rep).end = [];
                            end 
                        end
                    end
                end
            end
            FreezeData(cond).all(ifly).data = vFreeze;
        end
    end


    MT = NaN(6,1);
    for cond = 1:7   
        %empty data here -- so we don't factor in zeros
        FreezeData(cond).avgchange = NaN(1,num.fly);
        FreezeData(cond).avgdelay = NaN(1,num.fly);
        for ifly = 1:num.fly
            inputdata = FreezeData(cond).all(ifly).data;   
            output = struct('freezeresponse', MT, 'delay', MT, 'duration', MT, 'change', MT);
            for rep = 1:6
               idx = inputdata(rep).idx;
               % did the fly freeze on that trial?
               if isempty(idx)
                 output.freezeresponse(rep) = false;
               else
                 output.freezeresponse(rep) = true;  
                 output.delay(rep) = (idx(1)-60)/num.fps;
                 output.duration(rep) = length(idx)/num.fps;
                 output.change(rep) = inputdata(rep).speedchange;
               end
            end
            tot = (sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]));
            if tot >= 2 %must have a possibility of 2 trials for each fly to include
                output.freezecount = sum(output.freezeresponse);
                output.totalcount = sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]);
                output.avgdelay = nanmean(output.delay);
                output.avgduration = nanmean(output.duration);
                output.avgchange = nanmean(output.change);
            else
                output.freezecount = NaN;
                output.totalcount = NaN;
                output.avgdelay = NaN;
                output.avgduration = NaN;
                output.avgchange = NaN;
            end    
                FreezeData(cond).data(ifly) = output;
                FreezeData(cond).totalcount(ifly) = output.totalcount;
                FreezeData(cond).freezecount(ifly) = output.freezecount;
                FreezeData(cond).percent(ifly) = FreezeData(cond).freezecount(ifly)/FreezeData(cond).totalcount(ifly);
                FreezeData(cond).avgdelay(ifly) = output.avgdelay;
                FreezeData(cond).avgduration(ifly) = output.avgduration;
                FreezeData(cond).avgchange(ifly) = output.avgchange;
        end
    end
end

% -- line plot the frequency of freezing -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(FreezeData(cond).percent);
    y.err(cond) = sem(FreezeData(cond).percent,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([0, 1]);
set(gca,'TickDir','out');
title([filename(1:end-4) ' ' type ' freeze frequency'])
xlabel('Stimulation length (sec)')
ylabel('Frequency of freezing (%)')

% save the data
figname  = [fig_dir, structure_name ' ' type ' freezing lineplot endpoint ' num2str(stppoint)];
save_figure(fig, figname);



%---scatter plot Percent of flies that freeze
sz = 70;
Kolor = colorList(5,:);
fig = getfig; hold all
xlength = 0.2; xlim([0,8])
for cond = 1:7
    y = FreezeData(cond).percent;
    y(isnan(y)) = [];
    idx = length(y);
    x = cond:xlength/idx:cond+(xlength-xlength/idx);
    scatter(x, y, sz, Kolor, 'filled')
    avg = nanmean(y);
    err = sem(y,2);
%     errorbar(cond+(xlength/2), avg, err, 'Color', 'r')
    plot([cond-(xlength/2), cond+(xlength*1.5)], [avg,avg], 'Color', 'k', 'linewidth', 3)
%     errorbar(cond, avg, err, 'color', 'k')
end
title({[structure_name ' Frequency of freezing behavior. Window end: ' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Percent of trials freezing (%)')
figname  = [fig_dir, structure_name ' ' type ' percent freezing endpoint ' num2str(stppoint)];
set(gca,'TickDir','out');

save_figure(fig, figname);

%% Delay to freezing in walking flies (by fly)  -- COMBINE FOR MULTIPLE DATA SETS

condFrames = round(condLength*num.fps);
% Delay data
stppoint = 1.5;
windowEnd = (60+stppoint*num.fps)*ones(1,7);
% stppoint = 0.5;
% windowEnd = condFrames + 60 + round(stppoint*num.fps);
% windowEnd(1) = windowEnd(end);

condFrames = round(condLength*num.fps);
sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for cond = 1:7
    for ifly = 1:num.fly
        for rep = 1:6
            vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                  'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
            raw = data(cond).(type).raw(ifly).data(:,rep); 
            % if the data fits the behavior type, analyze
            if ~isnan(raw(1))
                temp = smooth(raw, 'rlowess');
                vStart = nanmean(raw(54:60)); %value prestim
%                 vChange = NaN;
                
                % only include trials with greater than min speed
                if vStart >= min_speed
                    % Change in speed
                    vStim = nanmean(raw(76:82)); %value prestim
                    vChange = vStim-vStart;
                    vFreeze(rep).speedchange = vChange;
                    %find locations that the fly isn't moving post stim start
                    vLow = find(temp(61:121) <= min_speed); 
                    if ~isempty(vLow)
                        vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                        % calculate the specific bits that qualify for 'freezing'
                        vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                        vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                        % find the start and end points of the 'freeze'
                        try % gaps present
                            vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                            vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                            vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                            if vFreeze(rep).start>vFreeze(rep).end
                                vFreeze(rep).start = vFreeze(rep).end;
                            end
                        catch % no gaps in the 'frozen' data points
                            vFreeze(rep).start = 1;
                            vFreeze(rep).end = length(vFreeze(rep).pos);
                        end
                        vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                        vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                        if vFreeze(rep).idx(1)>windowEnd(cond)
                            vFreeze(rep).idx = [];
                            vFreeze(rep).start = [];
                            vFreeze(rep).end = [];
                        end 
                    end
                end
            end
        end
        FreezeData(cond).all(ifly).data = vFreeze;
    end
end


MT = NaN(6,1);
for cond = 1:7   
    %empty data here -- so we don't factor in zeros
    FreezeData(cond).avgchange = NaN(1,num.fly);
    FreezeData(cond).avgdelay = NaN(1,num.fly);
    for ifly = 1:num.fly
        inputdata = FreezeData(cond).all(ifly).data;   
        output = struct('freezeresponse', MT, 'delay', MT, 'duration', MT, 'change', MT);
        for rep = 1:6
           idx = inputdata(rep).idx;
           % did the fly freeze on that trial?
           if isempty(idx)
             output.freezeresponse(rep) = false;
           else
             output.freezeresponse(rep) = true;  
             output.delay(rep) = (idx(1)-60)/num.fps;
             output.duration(rep) = length(idx)/num.fps;
             output.change(rep) = inputdata(rep).speedchange;
           end
        end
        tot = (sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]));
        if tot >= 2 %must have a possibility of 2 trials for each fly to include
            output.freezecount = sum(output.freezeresponse);
            output.totalcount = sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]);
            output.avgdelay = nanmean(output.delay);
            output.avgduration = nanmean(output.duration);
            output.avgchange = nanmean(output.change);
        else
            output.freezecount = NaN;
            output.totalcount = NaN;
            output.avgdelay = NaN;
            output.avgduration = NaN;
            output.avgchange = NaN;
        end    
            FreezeData(cond).data(ifly) = output;
            FreezeData(cond).totalcount(ifly) = output.totalcount;
            FreezeData(cond).freezecount(ifly) = output.freezecount;
            FreezeData(cond).percent(ifly) = FreezeData(cond).freezecount(ifly)/FreezeData(cond).totalcount(ifly);
            FreezeData(cond).avgdelay(ifly) = output.avgdelay;
            FreezeData(cond).avgduration(ifly) = output.avgduration;
            FreezeData(cond).avgchange(ifly) = output.avgchange;
    end
end


% -- line plot delay to freeze -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(FreezeData(cond).avgdelay);
    y.err(cond) = sem(FreezeData(cond).avgdelay,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([0, 1]);
title({[structure_name ' Freezing Delay. End point: stimend+' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Time to freeze initiation (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, structure_name ' ' type ' freezing delay linegraph ' num2str(stppoint)];

% save the data
figname  = [fig_dir, structure_name ' ' type ' freezing lineplot endpoint ' num2str(stppoint)];
save_figure(fig, figname);





% Time delay for freezing
fig = getfig; hold all
xlength = 0.2; xlim([0,8])
for cond = 1:7
    y = FreezeData(cond).avgdelay;
    y(isnan(y)) = [];
    idx = length(y);
    x = cond:xlength/idx:cond+(xlength-xlength/idx);
    scatter(x, y, sz, Kolor, 'filled')
    avg = nanmean(y);
    err = sem(y,2);
%     errorbar(cond+(xlength/2), avg, err, 'Color', 'r')
    plot([cond-(xlength/2), cond+(xlength*1.5)], [avg,avg], 'Color', 'k', 'linewidth', 3)
%     errorbar(cond, avg, err, 'color', 'k')
end
title({[structure_name ' Freezing Delay. End point: stimend+' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Time to freeze initiation (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, structure_name ' ' type ' freezing delay endpoint ' num2str(stppoint)];
ylim([0, 1.5])

save_figure(fig, figname);

%% Time to half-speed drop in walking flies
min_speed = 0.3;
fps = 30;
sSpan = 1;
stppoint = 0.2;
condFrames = round(condLength*fps);
windowEnd = condFrames + 60 + round(stppoint*fps);
windowEnd(1) = windowEnd(end);

controlRange = 55:61; %200ms preceeding stim onset
% stimRange = 61:88;

% find the half and min locations for change in speed for each fly
for dd = 1:2
    switch dd
       case 1 %SH control
           data = Cdata;
           num = CONTROL.num;
       case 2 %IN stim
           data = Sdata;
           num = STIM.num;
    end
    for cond = 1:7
        stimRange = [61:windowEnd(cond)];
        for ifly = 1:num.fly
            if ~isnan(nanmean(data(cond).walking.raw(ifly).data(1,:)))
                % include only >0.3cm/s walking flies
                temp = nanmean(data(cond).walking.raw(ifly).data(controlRange,:),1);
                filter = temp>min_speed;
                input = data(cond).walking.raw(ifly).data(:,filter);
                % find the average across all trials for each fly 
                input = nanmean(input,2);
                Y = diff(input)/(1/fps);
                Y = diff(input);
                [~, loc] = min(Y(stimRange));
                MinLoc = loc+60;
                data(cond).walking.peakslopeloc(ifly) = MinLoc;
%             figure; hold all; plot(input); plot(Y); vline(MinLoc)
            else
                data(cond).walking.peakslopeloc(ifly) = NaN;
            end
            
%             cntl = nanmean(data(cond).walking.raw(ifly).data(controlRange,:),2);
%             stim = nanmean(data(cond).walking.raw(ifly).data(stimRange,:),2);
%             if ~isnan(cntl(1))
%                 % min value:
%                 [~, loc] = (min(stim));
%                 MinLoc = loc+60;
%                 HLoc = round(loc/2)+60;
%                 data(cond).walking.MinLoc(ifly) = MinLoc;
%                 data(cond).walking.HLoc(ifly) = HLoc;
%                 % find the value at the min 
%                 data(cond).walking.CntlVal(ifly) = mean(data(cond).walking.avg(controlRange,ifly));
%                 data(cond).walking.MinVal(ifly) = data(cond).walking.avg(MinLoc,ifly);
%                 
% %                 %PLOT preview
% %                 fig = getfig; hold all
% %                 plot(x, nanmean(data(cond).walking.raw(ifly).data,2))
% %                 vline(x(HLoc))
% %                 vline([x(stimRange(1)), x(stimRange(end))], 'b')
% %                 vline([x(controlRange(1)), x(controlRange(end))], 'm')
% %                 scatter(x(stimRange), nanmean(data(cond).walking.raw(ifly).data(stimRange,:),2))
% %                 scatter(x(MinLoc), yMin, 60, 'r', 'filled')
%             else
%                 data(cond).walking.MinLoc(ifly) = NaN;
%                 data(cond).walking.HLoc(ifly) = NaN;
%                 data(cond).walking.CntlVal(ifly) =    
%                 data(cond).walking.MinVal(ifly) = data(cond).walking.avg(MinLoc,ifly);
%             end
        end
    end
    switch dd
       case 1 %SH control
           Cdata = data;
       case 2 %IN stim
           Sdata = data;
    end
    clear data
end
clear windowEnd stppoint stimRange stim num MinLoc loc InColor ifly HLoc cntl
   
% plot the data by condition -- across both lines
Colorlist = [hex2rgb('000000'); hex2rgb('f15a24')];

X = -0.5:(1.5/45):1;
xRange = 46:91;
fig = getfig; hold all
    for cond = 2:7
    % control
    sSpan = 1;
    loc = (nanmean(Cdata(cond).walking.peakslopeloc));
    steep = [floor(loc)+1, ceil(loc)+1];
    input = smooth(Cdata(cond).walking.AVG, sSpan);%, 'rlowess'
    plot(x(steep), input(steep), 'color', Color('yellow'), 'LineWidth', 3)
    plot(X, input(xRange), 'color', Colorlist(1,:), 'LineWidth', 1)
    
    % stim
    loc = (nanmean(Sdata(cond).walking.peakslopeloc));
    steep = [floor(loc)+1, ceil(loc)+1];
    input = smooth(Sdata(cond).walking.AVG, sSpan);%, 'rlowess'
    plot(x(steep), input(steep), 'color', 'y', 'LineWidth', 3)
    plot(X, input(xRange), 'color',  Colorlist(2,:), 'LineWidth', 1)
    end
    
figname = [comp ' location of greatest slowing value'];   
title({figname;' ';' '})
xlabel('Time (sec)')
ylabel('speed (cm/s)')
set(gca,'TickDir','out');

save_figure(fig, [fig_dir, figname]);   
  

% 
% % Plot of time delay to greatest slowing points
% 
% a = max([CONTROL.num.fly, STIM.num.fly]);
% MT = NaN(a,21);
% idx = 0;
% for cond = 1:7
%     idx = idx+1;
%     MT(1:CONTROL.num.fly,idx) = Cdata(cond).walking.peakslopeloc';
%     MT(1:STIM.num.fly,idx+1) = Sdata(cond).walking.peakslopeloc'; 
%     idx = idx+2;
% end
% fig = getfig; hold all       
% boxplot(MT)
% % title([' ])
%         
%         
%     % control
%     sSpan = 1;
%     loc = (nanmean(Cdata(cond).walking.peakslopeloc));
%     steep = [floor(loc)+1, ceil(loc)+1];
%     input = smooth(Cdata(cond).walking.AVG, sSpan);%, 'rlowess'
%     plot(x(steep), input(steep), 'color', Color('yellow'), 'LineWidth', 3)
%     plot(X, input(xRange), 'color', Colorlist(1,:), 'LineWidth', 1)
%     
%     % stim
%     loc = (nanmean(Sdata(cond).walking.peakslopeloc));
%     steep = [floor(loc)+1, ceil(loc)+1];
%     input = smooth(Sdata(cond).walking.AVG, sSpan);%, 'rlowess'
%     plot(x(steep), input(steep), 'color', 'y', 'LineWidth', 3)
%     plot(X, input(xRange), 'color',  Colorlist(2,:), 'LineWidth', 1)
%     end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %PLOT condition average
% fig = getfig;
%     hold all
%     cond = 1;
%     plot(x, data(cond).walking.AVG, 'k')
%     loc = round(nanmean(data(cond).walking.HLoc));
%     scatter(x(loc), data(cond).walking.AVG(loc), 60, 'r', 'filled')
%     
%     for cond = 3:2:7
%     plot(x, data(cond).walking.AVG, 'b')
%     loc = round(nanmean(data(cond).walking.HLoc));
%     scatter(x(loc), data(cond).walking.AVG(loc), 60, 'b', 'filled')
%     end
% 
% 
% 
% 
% set(gca,'TickDir','out');
% 
% figname  = [fig_dir, structure_name ' change in speed after stimulus'];
% save_figure(fig, figname);
% 

