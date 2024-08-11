
clear; clc; close all

% Load fast logging temperature data
folder = getCloudPath;
folder = [folder(1:end-5) 'Temp logging/'];
fileName = 'D3_-_2814794_Aug_11_2024_6_03_49_PM.csv';
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
% figure; plot(temperature_F)

% Convert to Celcius
temperature_C = FtoC(temperature_F); 
time_diff = [nan; minutes(timeDifferences)];
tempRate = [nan; diff(temperature_C)./minutes(timeDifferences)];

% save processed data:
T = table(timeStrings, temperature_C, temperature_F, ...
    humidity, pressure,heat_index, dew_point, time_diff, tempRate);

save([temp_log(1:end-4) '.mat'],'T')

%% Summary figure
r = 1;
c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;
kolor = Color('dodgerblue');

% temp rate histogram
fig = getfig('',1,[772 420]); 

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
title(['Range: ' num2str(tempRange) '\circC | STD: ' num2str(tempSTD) '\circC'])

formatFig(fig, false,[r,c],sb);

save_figure(fig, [temp_log(1:end-4) ' timecourse'],'-pdf');

% tempRange = range(T.temperature_C);





%%

% sampling intrinsic features: 
% 2 second sampling, 0.

%% try a 5 minute smoothing of the data


test = data;


