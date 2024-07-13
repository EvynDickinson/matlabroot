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




