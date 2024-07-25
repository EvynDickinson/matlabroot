
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











