
dataPath = 'G:\My Drive\Jeanne Lab\NOAA data\2022 SubHourly Environmental Data.xlsx';
data = readmatrix(dataPath);
time = data(:,3);
temp = data(:,9);
day = data(:,2);

% find the month timepoints
dayStart = find(time==5);
fiveminrate = (diff(temp))/5;
outliers = fiveminrate>30 | fiveminrate<-30;
exOut = fiveminrate(outliers);
fiveminrate(outliers) = nan;
DaysPerMonth = [1 31 28 31 30 31 30 31 31 30 31 30];
MonthStart = dayStart(cumsum(DaysPerMonth));
monthNames = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
% Planned temp rate changes
slow_rate = [0.1,-0.1];
mid_rate = [0.25, -0.25];
high_rate = [0.5, -0.5];

% PLOT FIGS
nrows = 3;
ncols = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;

fig = figure; set(fig, 'pos',[285 382 955 596])
% Temp rate histogram
subplot(nrows, ncols, sb(2).idx)
h = histogram(fiveminrate);
h.EdgeColor = Color('grey');
set(gca, 'YScale', 'log')
ylabel('Count (five minute periods)')
xlabel('\DeltaT (\circC/min)')
v_line(slow_rate, 'Cyan','-',2)
v_line(mid_rate, 'DarkViolet','-',2)
v_line(high_rate, 'Orange','-',2)

% Temp rate by day
subplot(nrows, ncols, sb(1).idx)
plot(fiveminrate,'linewidth', 2, 'Color',Color('grey'))
set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
ylabel('\DeltaT (\circC/min)')
title('Owls Head Maine 2021')
xlabel('Month')
axis tight
h_line(slow_rate, 'Cyan','-',1)
h_line(mid_rate, 'DarkViolet','-',1)
h_line(high_rate, 'Orange','-',1)
fig = formatFig(fig, true,[nrows,ncols], sb);

save_figure(fig, 'G:\My Drive\Jeanne Lab\NOAA data\Owls Head 2022 annual temp rate','-png');



%% Target temp ramp protocols

order = [1,2,3]; %S,M,F
% Ramp durations
fullTimes = [180,72,36];
halfTimes = [90, 36, 18];
cList = {'Cyan','DarkViolet','Orange'};
% Target Temperatures
maxTemp = 25;
minTemp = 7;
startTemp = 16;
startTime = 15;

% Protocol steps:
%START HOLD:
step(1).time = [0,startTime];
step(1).temp = [startTemp, startTemp];
% INCREASE TO TEMP MAX at order 1 rate
tEnd = step(1).time(end);
step(2).time = [tEnd,tEnd+halfTimes(order(1))];
step(2).temp = [startTemp, maxTemp];
% DECREASE TO TEMP MIN at order 1 rate
tEnd = step(1).time(end);
step(2).time = [tEnd,tEnd+halfTimes(order(1))];
step(2).temp = [startTemp, startTemp];

% Graph...
fig = figure; hold on
plot(step(1).time,step(1).temp,'color', 'w')
plot(step(2).time,step(2).temp,'color', Color(cList{order(1)}))


%% 

tempLog = readmatrix("G:\My Drive\Jeanne Lab\DATA\02.01.2022\TempRatePlantFood_RampLog(1).csv");
x = tempLog(:,1); %time (in seconds)
x = x./(60*60);
y1 = tempLog(:,3); %set temp
y2 = tempLog(:,2); %actual temp

fig = figure; hold on
plot(x,y1, 'Color',Color('red'), 'LineWidth',1, 'LineStyle',':')
plot(x,y2, 'Color',Color('white'), 'LineWidth',2)
ylabel('Temp (\circC)')
xlabel('Time (hr)')
fig = formatFig(fig, true);
save_figure(fig, 'G:\My Drive\Jeanne Lab\DATA\Temp Control\Exp 1 temp readout','-png');











