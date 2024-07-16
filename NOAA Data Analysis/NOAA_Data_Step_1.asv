
%% NOAA data column headers

% #    Name                                        Units
% 1    WBANNO                                  XXXXX
% 2    UTC_DATE                                 YYYYMMDD
% 3    UTC_TIME                                   HHmm
% 4    LST_DATE                                   YYYYMMDD
% 5    LST_TIME                                     HHmm
% 6    CRX_VN                                     XXXXXX
% 7    LONGITUDE                               Decimal_degrees
% 8    LATITUDE                                    Decimal_degrees
% 9    AIR_TEMPERATURE                    Celsius
% 10   PRECIPITATION                         mm
% 11   SOLAR_RADIATION                  W/m^2
% 12   SR_FLAG                                    X
% 13   SURFACE_TEMPERATURE         Celsius
% 14   ST_TYPE                                      X
% 15   ST_FLAG                                    X
% 16   RELATIVE_HUMIDITY                 %
% 17   RH_FLAG                                   X
% 18   SOIL_MOISTURE_5                    m^3/m^3
% 19   SOIL_TEMPERATURE_5             Celsius
% 20   WETNESS                                   Ohms
% 21   WET_FLAG                                X
% 22   WIND_1_5                                m/s
% 23   WIND_FLAG                             X
% ------------------------------------------------------
% 24 LIGHT CORRECTED                    (ESD 2024)
% 25 SURFACE   CORRECTED            (ESD 2024)
% 26 AIR_TEMP TEMP CORRECTED   (ESD 2024)
% 27 SURFACE TEMP RATE                 (ESD 2024)
% 28 AIR TEMP RATE                           (ESD 2024)

% X = linspace(-5,5,300); % these were the samples we used to generate the distributions


T_labels = {'WBANNO','UTC_DATE', 'UTC_TIME', 'LST_DATE', 'LST_TIME', 'CRX_VN', 'LONGITUDE', ...
                    'LATITUDE', 'AIR_TEMPERATURE','PRECIPITATION','SOLAR_RADIATION', 'SR_FLAG',...
                    'SURFACE_TEMPERATURE', 'ST_TYPE', 'ST_FLAG', 'RELATIVE_HUMIDITY', 'RH_FLAG', ...
                    'SOIL_MOISTURE_5', 'SOIL_TEMPERATURE_5', 'WETNESS', 'WET_FLAG', 'WIND_1_5',...
                    'WIND_FLAG', 'LIGHT_CORRECTED', 'SURFACE_TEMP_CORRECTED', 'AIR_TEMP_CORRECTED',...
                    'SURFACE_TEMP_RATE', 'AIR_TEMP_RATE'};
TH = [];
for i = 1:length(T_labels)
    TH.(T_labels{i}) = i;
end
clear T_labels i

%% Load data structure

% Folder structures: 
baseFolder = 'H:\NOAA Data\Fully sampled\'; % path to raw data base folder
load([baseFolder 'Fully sampled NOAA sites.mat'])

% pull some preliminary parameters
sites = unique({data(:).location});
n = [];
n.Sites = length(sites);
n.Samples = length(baseFolder);
n.yearList = 2006:2023;
n.years = length(n.yearList);

switch questdlg('Select figure parameters','Fig Type','White PDF', 'Black PNG','Cancel','White PDF')
    case 'White PDF'
        fig_type = '-pdf';
        blkbgd = false;
    case 'Black PNG'
        fig_type = '-png';
        blkbgd = true;
    case {'Cancel',''}
        return
end

initial_vars = who;
initial_vars{end+1} = 'initial_vars';

%% Plot the site locations on the map
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% need to make a sort index of some kind that we can just add to the data matrix
% and use to separate by latitude in the future

% get unique sampling names: 
[latitude, longitude] = deal([]);
all_sites = {data(:).location};
for i = 1:n.Sites
    loc = find(strcmp(sites{i},all_sites));
    [data(loc).siteID] = deal(i);
    latitude(i) = data(loc(1)).latitude;
    longitude(i) = data(loc(1)).longitude;
end 
% 
[~, lat_idx] = sort(latitude);
sorted_latitude = latitude(lat_idx);
sorted_longitude = longitude(lat_idx);

% color points by latitude (N-S axis)
CList = Color('indigo','lavender', n.Sites);

% Plot locations
fig = getfig('',1,[618 447]); 
    geoscatter(sorted_latitude, sorted_longitude, 75,[0 0 0],'filled','^')
    hold on
    geoscatter(sorted_latitude, sorted_longitude, 40, CList,'filled','^')
    geobasemap colorterrain
    set(fig, 'color', backColor); 
    % formatting
    ax = gca;
    set(ax, 'fontsize', 15, 'grid', 'off')
    ax.LongitudeLabel.Color = foreColor;
    ax.LatitudeLabel.Color = foreColor;
    ax.LongitudeAxis.Color = foreColor;
    ax.LatitudeAxis.Color = foreColor;

% ax.LongitudeLimits = [-175.3770 -62.7130];
% ax.LatitudeLimits = [18.2455 72.6145];

save_figure(fig,[baseFolder 'Figures/Sample sites map'],fig_type);

%% Plot the yearly temperature for example site: 
clearvars('-except',initial_vars{:})

CList = Color('Red', 'Gold', n.years);
LW = 0.5;

city = 1;
city_loc = find([data(:).siteID]==city);

fig = getfig('', 1); hold on
    for i = 2:length(city_loc)
        city = city_loc(i);
        y = data(city).T(:,TH.AIR_TEMP_CORRECTED);
        plot(y, 'Color', CList(i,:),'linewidth', LW)
        % h = cdfplot(data(city).T(:,TH.AIR_TEMP_RATE));  %
        % set(h,'color',CList(i,:),'linewidth',LW)
        % clear h
    end
formatFig(fig, blkbgd);
xlabel('Time (1 year)')
ylabel('Temperature (\circC)')


%% Site-based temprate histogram (Combine temperature rates for each site over the year)
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% CList = Color('Red', 'Gold', n.years);
LW = 0.5;

for city = 1:n.Sites
    city_loc = find([data(:).siteID]==city);
    plottingData = [];
    temp = struct; 
    for i = city_loc
        temp.R = data(city_loc(i)).T(:,TH.AIR_TEMP_RATE);
        temp.T = data(city_loc(i)).T(:,TH.AIR_TEMP_CORRECTED);
        plottingData = [plottingData; temp.T, temp.R];
    end
        
fig = getfig('',1)
 
fig = getfig('', 1); hold on
    for i = 2:length(city_loc)
        city = city_loc(i);
        y = data(city).T(:,TH.AIR_TEMP_CORRECTED);






















