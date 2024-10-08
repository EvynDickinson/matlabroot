
function zeitgeber_time = getCircadianTime(expTime, incubator, expDate)
% zeitgeber_time = getCircadianTime(expTime, incubator, expDate)
% convert the start time of an experiment to the circadian time 
% by accounting for the night shifted incubator and daylight savings time

% daylight savings regions: 
DLS.on = {'03.14.2021', '03.13.2022', '03.12.2023', '03.10.2024', '03.9.2025'};
DLS.off = {'10.07.2021', '10.06.2022', '10.05.2023', '10.03.2024', '10.02.2025'};

expYear = str2double(expDate(7:end));
trial_date = datetime(expDate, 'InputFormat','MM.dd.yyyy');

% find year
idx = expYear-2020;
TIMESHIFT_ON = datetime(DLS.on{idx}, 'InputFormat','MM.dd.yyyy');
TIMESHIFT_OFF = datetime(DLS.off{idx}, 'InputFormat','MM.dd.yyyy');

if trial_date>=TIMESHIFT_ON && trial_date<=TIMESHIFT_OFF
    hr_offset = 0; % these experiments have light-on at 8am
else 
    hr_offset = 1; % these experiements have light-on at 7am
end

% get the start time of the experiment in duration...
if ischar(expTime)
   startdur = duration(expTime);
else % contingency plan for non-duration values
    newTime = datetime(expTime,"ConvertFrom", "excel");
    [h,m,s] = hms(newTime); % Get the hour, minute, and second values
    startdur = duration(h,m,s); % Get the output as a duration() array
end

% night shifted incubator adjustment
if strcmpi(incubator,'N') 
    hr_offset = hr_offset + 7; % night incubator is 7 hours off
end

% account for daylight savings time: 
startdur = startdur - duration(hr_offset,0,0); 

% convert to zeitgeber time
startdur = startdur - duration(8,0,0);  % hours since lights on

if startdur<0
    startdur = startdur + duration(24,0,0); % jump the to positive side of things (loop around 24hrs)
end

zeitgeber_time = round(hours(startdur),4);