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

%% Load data from the downloaded files from:
% [https://www.ncei.noaa.gov/pub/data/uscrn/products/subhourly01/]

% Folder structures: 
baseFolder = 'H:\NOAA Data\'; % path to raw data base folder
year_list = 2006:2023; % years that have data
nYears = length(year_list);

% Data check / cleaning values
minT = -90; % lowest recorded temperature on record
maxT = 57; %highest recorded temp on record 
light_err = -99999;

% master_matrix = [];

% load each data file and save to a processed data folder in the year folder  locations:
for year = 1:nYears
    % pull the names of the city sampling locations
    year_folder = [baseFolder num2str(year_list(year)) '\'];
    fileList = dir([year_folder 'Raw Data\*.txt']);
    city = struct;
    nCities = length(fileList);
    for i = 1:nCities
        city(i).name = fileList(i).name(1:end-4);

        fullPath = [year_folder 'Raw Data\' fileList(i).name];
        city_data = readmatrix(fullPath);

        % pull meta data for the city
        tempStr = strsplit(city(i).name,'-');
        city(i).location = tempStr{4};
        city(i).longitude = city_data(1,7);
        city(i).latitude = city_data(1,8);
        city(i).year = year_list(year);

        % +++++++++++++++ clean data for errors +++++++++++++++

        % catch if there is inconsistent data trial size
        if ~(size(city_data,2)==23) %not the correct number of data columns
            disp(['Missing number of columns for ' num2str(year_list(year)) ' in ' city(i).location])
            continue
        end
        
        % Solar radiation power spectrum
        light = city_data(:,11); % solar radiation data
        %clean light data
        light(light == light_err) = nan;
        light_filled = fillmissing(light, 'linear');
        city_data(:,24) = light_filled; % resave the data to the og structure
        clear light light_filled
    
        % Clean up the SURFACE temperature data
        temperature = city_data(:,13);
        error_loc = (temperature>=maxT | temperature<=minT);
        temperature(error_loc) = nan;
        % add a check for the gap being too large
        temperature_filled = fillmissing(temperature, 'linear');
        city_data(:,25) = temperature_filled;
        clear temperature_filled temperature
    
        % Clean up the AIR temperature data
        temperature = city_data(:,9);
        error_loc = find(temperature>=maxT | temperature<=minT);
        temperature(error_loc) = nan;
        temperature_filled = fillmissing(temperature, 'linear');
        city_data(:,26) = temperature_filled;
        clear temperature_filled temperature

        % +++++++++++++++ calculate new variables +++++++++++++++
        
        % make an index list for the day starts:
        time =  city_data(:,5);
        daily_date = city_data(:,4); % column with date
        day_starts = find(time==0);
        nDays = length(day_starts); % how many days of data       
        
        dateList = daily_date(day_starts);
        city(i).nDays = nDays;
        city(i).dates = dateList;
        city(i).dayStart_idx = day_starts;

        % Create a datetime array for the base date
        datetimes = datetime(daily_date, 'ConvertFrom', 'yyyymmdd');
        hours_part = floor(time / 100);  % Extract hours
        minutes_part = mod(time, 100);   % Extract minutes
        fullDateTime = datetimes + hours(hours_part) + minutes(minutes_part);
        time_gap = minutes(diff(fullDateTime));

        % Check if all differences are 5 minutes
        expectedDiff = 5;
        all_five_minutes_apart = all(time_gap == expectedDiff);
        city(i).noskippedsamples = all_five_minutes_apart;

        % Determine the rates of temperature change
        surf_temp = city_data(:,25);
        surface_TR = [diff(surf_temp)]./time_gap;
        city_data(:,27) = [nan; surface_TR];

        air_temp = city_data(:,26);
        air_TR = diff(air_temp)./time_gap;
        city_data(:,28) = [nan; air_TR];

        % calculate the PDF and CDF of temperature rate for the year (both surface and air)
        X = linspace(-5,5,300); 
        
        % Surface temperature data fit
        mu = mean(surface_TR,'omitnan');
        sigma = var(surface_TR,'omitnan');
        gm = gmdistribution(mu,sigma);
        dataPDF = pdf(gm,X');
        dataCDF = cdf(gm,X');

        city(i).surfTR.mu = mu;
        city(i).surfTR.sigma = sigma;
        city(i).surfTR.dataPDF = dataPDF;
        city(i).surfTR.dataCDF = dataCDF;

        % Air temperature data fit
        mu = mean(air_TR,'omitnan');
        sigma = var(air_TR,'omitnan');
        gm = gmdistribution(mu,sigma);
        dataPDF = pdf(gm,X');
        dataCDF = cdf(gm,X');

        city(i).airTR.mu = mu;
        city(i).airTR.sigma = sigma;
        city(i).airTR.dataPDF = dataPDF;
        city(i).airTR.dataCDF = dataCDF;

        city(i).T = city_data;
                
        % Save the data as both a large structure, and smaller ones broken down by data type
        
        local_path = [year_folder 'Processed Data\'];
        if ~isfolder(local_path)
            mkdir(local_path)
        end
        
        % save a local copy of the data: 
        data = city(i);
        save([local_path city(i).name '.mat'], 'data', "T_labels", '-v7.3')
            
        disp(['Finished ' city(i).location ' from ' num2str(city(i).year) ' file ' num2str(i) '/' num2str(nCities)])

    end

    % save a copy of the full year 
    save([local_path num2str(city(i).year) ' processed NOAA data.mat'], 'city', "T_labels", '-v7.3')
        
    disp(['Finished '  num2str(city(i).year)])
end


%% Quickly make a master list of the cities and their data

clear 

% Folder structures: 
baseFolder = 'H:\NOAA Data\'; % path to raw data base folder
year_list = 2006:2023; % years that have data
nYears = length(year_list);

all_locations = [];

% originally, there were 3003 sites across the year span of 2006-2023

% load each data file and save to a processed data folder in the year folder  locations:
for year = 1:nYears
    % pull the names of the city sampling locations
    year_folder = [baseFolder num2str(year_list(year)) '\Processed Data\'];
    fileList = dir([year_folder '*.mat']);
    fileList(1) = [];
    all_locations = [all_locations; {fileList(:).name}'];
end

% Functions to extract the desired parts and remove the '.mat'
extractPart = @(name) strsplit(name, {'-'});
extractAndTrimPart = @(name) name(1:end-4);
% Use cellfun to apply the function to each element
splitNames = cellfun(extractPart, all_locations, 'UniformOutput', false);
% Extract the desired part from each split file name
desiredParts  = cellfun(@(parts) parts{4}, splitNames, 'UniformOutput', false);
trimmedNames = cellfun(@(name) extractAndTrimPart(name), desiredParts, 'UniformOutput', false);
siteList = unique(trimmedNames);

save([baseFolder 'Site List.mat'],'siteList')


    
%% reload the files and create a master file for them

master_matrix = cell(length(siteList)+1,nYears+1);
master_matrix(2:end,1) = siteList;
% field_names = fields(data);

% load each data file and save to a processed data folder in the year folder  locations:

for year = 1:nYears
    idx = 0;
    % pull the names of the city sampling locations
    year_folder = [baseFolder num2str(year_list(year)) '\Processed Data\'];
    fileList = dir([year_folder '*.mat']);
    city = struct;
    nCities = length(fileList); 
    master_matrix{1,year+1} = year_list(year);
    for i = 2:nCities %skip the first file, which is the failed group file
        temp = load([year_folder fileList(i).name],'data');
        % all the data points were sampled
        dayCheck = (temp.data.nDays==365 ||temp.data.nDays==366) ;
        sampleCheck = temp.data.noskippedsamples;

        loc = find(strcmp(temp.data.location, siteList));
        if  sampleCheck && dayCheck
            idx = idx+1;
            if idx == 1
                data = temp.data;
            else
                data(idx) = temp.data;
            end
            master_matrix{loc+1,year+1} = true; % write yes into the matrix
        else
            master_matrix{loc+1,year+1} = false; % write no into the matrix
        end
        disp([num2str(i) '/' num2str(nCities)])
    end
    save([baseFolder num2str(year_list(year)) ' NOAA Data.mat'],'data',  '-v7.3')
end

save([baseFolder 'NOAA Data MasterList.mat'],'master_matrix', 'siteList', '-v7.3')



%% Find an example of a site that has had consistent temperature data for all the years and use that for a power analysis?

% get some prelim stats on the meta data: 
nSites = size(master_matrix,1)-1;

masterSheet = false(nSites, nYears);
for i = 1:nSites
    for t = 1:nYears
        dum = master_matrix{i+1, t+1};
        if dum
            masterSheet(i,t) = true;
        end
    end
end
            
sites.isdata = masterSheet;
sites.master_sheet = master_matrix;
sites.names = siteList;
sites.years = year_list;

fullDataSite = find(sum(masterSheet,2)==nYears);
sites.fullySampled = siteList(fullDataSite);


%% Make a data matrix with the 82 sites that are fully sampled?
fullySampled = siteList(fullDataSite);

data = [];

idx = 0;
for city = 11: length(fullySampled)
    city_name = fullySampled{city};
    for year = 1:nYears
        idx = idx + 1;
        % pull the names of the city sampling locations
        year_folder = [baseFolder num2str(year_list(year)) '\Processed Data\'];
        fileList = dir([year_folder '*' city_name '.mat']);
        temp  = load([year_folder, fileList.name],'data');
   
        if idx == 1
            data = temp.data;
        else
            data(idx) = temp.data;
        end
    end
end

save([baseFolder 'Fully sampled NOAA sites.mat'],'data',  '-v7.3')






















