
clear; clc; close all

% Load fast logging temperature data
folder = getCloudPath;
folder = [folder(1:end-5) 'Temp logging/'];
fileList = dir([folder '*.csv']);
fileList = {fileList(:).name};
idx = listdlg("ListString",fileList,'PromptString','Select the temp log to load', 'ListSize',[300,300]);
if isempty(idx)
    return
end
fileName = fileList{idx};
% fileName = 'D3_-_2814794_Aug_11_2024_6_03_49_PM.csv';
% fileName = 'D3_-_2814794_Aug_11_2024_3_53_14_PM.csv';
temp_log = [folder fileName];

% Load the data into a matrix
raw_data = readcell(temp_log);

% pull the important data columns
raw_time = raw_data(6:end,2);
raw_temp = raw_data(6:end,3);
timeStrings = cellfun(@(x) char(x), raw_time, 'UniformOutput', false);

% Find the time difference from point to point
timeDurations = duration(timeStrings, 'InputFormat', 'hh:mm:ss');
timeDifferences = diff(timeDurations);
timeDifferences(timeDifferences < 0) = timeDifferences(timeDifferences < 0) + hours(12);
% quick check plot: 
% figure; plot(timeDifferences)

% Use cellfun to split each string and extract the appropriate data
% columns = {'temperature', 'humidity', 'pressure', 'heat_index', 'dew_point'};
splitData = cellfun(@(x) strsplit(x, ','), raw_temp, 'UniformOutput', false);
temperature_F = cellfun(@(x) str2double(x{2}), splitData);
humidity = cellfun(@(x) str2double(x{3}), splitData);
pressure = cellfun(@(x) str2double(x{4}), splitData);
heat_index = cellfun(@(x) str2double(x{5}), splitData);
dew_point = cellfun(@(x) str2double(x{6}), splitData);
AMPM = cellfun(@(x) (x{1}), splitData,'UniformOutput', false);
% figure; plot(temperature_F)

% Convert to Celcius
temperature_C = FtoC(temperature_F); 
time_diff = [nan; minutes(timeDifferences)];
tempRate = [nan; diff(temperature_C)./minutes(timeDifferences)];

% save processed data:
T = table(timeStrings, AMPM, temperature_C, temperature_F, ...
    humidity, pressure,heat_index, dew_point, time_diff, tempRate);

% Summary figure
r = 1;
c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;
kolor = Color('dodgerblue');

fig = getfig('',1,[772 556]); 
% temp rate histogram
subplot(r,c,sb(2).idx)
    histogram(T.tempRate,'FaceColor',kolor,'FaceAlpha',1); 
    xlabel('temp rate (\circC/min)'); 
    ylabel('count'); 
% temperature timecourse
subplot(r,c,sb(1).idx)
    plot(cumsum(T.time_diff,'omitnan'),T.temperature_C,'color', kolor,'linewidth', 1.5)
    xlabel('time (min)')
    ylabel('temp (\circC)')
    tempRange = range(T.temperature_C);
    tempSTD = std(T.temperature_C);
formatFig(fig, false,[r,c],sb);
subplot(r,c,sb(1).idx)
title({['Range: ' num2str(tempRange) '\circC | STD: ' num2str(tempSTD) '\circC'];...
        [T.timeStrings{1} T.AMPM{1} ' to ' T.timeStrings{end} T.AMPM{end}]},'fontsize', 12)

save_figure(fig, [temp_log(1:end-4) ' timecourse'],'-pdf');


% Save the data
if strcmp(questdlg('save data set to folder?'),'Yes')
    save([temp_log(1:end-4) '.mat'],'T')
end



%% Smoothed vs fully sampled

x = cumsum(T.time_diff,'omitnan');
y = T.tempRate;

% find 5-min period temp rate change
per = mean(T.time_diff,'omitnan'); %sample period
sSpan = ((1/per)*5);
y_smooth = smooth(T.temperature_C,sSpan, 'moving');

tempRange = range(T.temperature_C);
tempSTD = std(T.temperature_C);
title_str = {['Range: ' num2str(tempRange) '\circC | STD: ' num2str(tempSTD) '\circC'];...
        [T.timeStrings{1} T.AMPM{1} ' to ' T.timeStrings{end} T.AMPM{end}]};

% Plot zoomed in timecourse
fig = getfig('',1,[772 556]); 

% temp rate histogram
subplot(r,c,sb(2).idx)
    histogram(y,'FaceColor',kolor,'FaceAlpha',1); 


    v_line(mean(y),'r')
    %labels
    xlabel('temp rate (\circC/min)'); 
    ylabel('count'); 
    
% temperature timecourse
subplot(r,c,sb(1).idx); hold on
    plot(x,T.temperature_C,'color', kolor,'linewidth', 1.5)
    plot(x,y_smooth,'color', Color('black'),'linewidth', 1)%5 min smoothed
    %labels
    xlabel('time (min)')
    ylabel('temp (\circC)')
    xlim(xlimits)
% formatting and title    
formatFig(fig, false,[r,c],sb); subplot(r,c,sb(1).idx)
title(title_str,'fontsize', 12)

save_figure(fig, [temp_log(1:end-4) ' timecourse'],'-pdf');


%% Compile data logs

xlimits = [1,600]; %[215, 232]; %
x = cumsum(T.time_diff,'omitnan'); % cumulative time for whole time course
ROI = find(x<=xlimits(2) & x>=xlimits(1)); % find index for desired range
x_temp = x(ROI);

% fully sampled temp rate:
yF_rate = T.tempRate(ROI); 
yF_temp = T.temperature_C(ROI);

% find 5-min period temp rate change
per = mean(T.time_diff,'omitnan'); %sample period
sSpan = ((1/per)*5);
yS_temp = smooth(T.temperature_C(ROI),sSpan, 'moving');
yS_rate = (diff(yS_temp)./T.time_diff(ROI(2:end)));

tempRange = range(T.temperature_C(ROI));
tempSTD = std(T.temperature_C(ROI));
title_str = {['Range: ' num2str(tempRange) '\circC | STD: ' num2str(tempSTD) '\circC'];...
        [T.timeStrings{ROI(1)} T.AMPM{ROI(1)} ' to ' T.timeStrings{ROI(end)} T.AMPM{ROI(end)}]};

% Plot zoomed in timecourse
fig = getfig('',1,[772 556]); 

% temp rate histogram
subplot(r,c,sb(2).idx); hold on
    h = histogram(yF_rate,'FaceColor',kolor,'FaceAlpha',0.5,'EdgeColor',kolor); 
    histogram(yS_rate,'FaceColor',Color('grey'),'FaceAlpha',0.3,'EdgeColor',Color('grey')); 
    % v_line(mean(y),'r')
    %labels
    xlabel('temp rate (\circC/min)'); 
    ylabel('count'); 
    xlim([-25,25])
    
% temperature timecourse
subplot(r,c,sb(1).idx); hold on
    plot(x_temp,yF_temp,'color', kolor,'linewidth', 1.5)
    plot(x_temp,yS_temp,'color', Color('black'),'linewidth', 1)%5 min smoothed
    %labels
    xlabel('time (min)')
    ylabel('temp (\circC)')
    % xlim(xlimits)
% formatting and title    
formatFig(fig, false,[r,c],sb); subplot(r,c,sb(1).idx)
title(title_str,'fontsize', 12)

save_figure(fig, [temp_log(1:end-4) ' timecourse ROI ' num2str(xlimits(1)) '-' num2str(xlimits(2))],'-png');








