clear; clc
%% Format and process annual temperature data from NOAA

cityNames = {'Spokane WA', 'Boulder CO', 'Versailles KY', 'Owls Head ME', 'Austin TX'};


sheet = cityNames{4};
dataPath = 'G:\My Drive\Jeanne Lab\NOAA data\2021 SubHourly Environmental Data.xlsx';
data = readmatrix(dataPath,'Sheet', sheet);
time = data(:,3);
temp = data(:,9);
day = data(:,2);

% Temp data pull:
loc = temp>100 | temp<-50;
temp(loc) = nan;


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
ncols = 4;
sb(1).idx = 1:3;
sb(2).idx = [5:7,9:11];
sb(3).idx = [4,8,12];

fig = figure; set(fig, 'pos',[285 382 955 596])
% Temp rate histogram
subplot(nrows, ncols, sb(2).idx)
h = histogram(fiveminrate);
h.EdgeColor = Color('grey');
set(gca, 'YScale', 'log')
ylabel('Count (five minute periods)')
xlabel('\DeltaT (\circC/min)')
% v_line(slow_rate, 'Cyan','-',2)
% v_line(mid_rate, 'DarkViolet','-',2)
% v_line(high_rate, 'Orange','-',2)

% Temp rate by day
subplot(nrows, ncols, sb(1).idx)
plot(fiveminrate,'linewidth', 2, 'Color',Color('grey'))
set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
ylabel('\DeltaT (\circC/min)')
title_str = [sheet ' 2021'];
title(title_str)
xlabel('Month')
axis tight
% h_line(slow_rate, 'Cyan','-',1)
% h_line(mid_rate, 'DarkViolet','-',1)
% h_line(high_rate, 'Orange','-',1)

% Yearly temp histogram
subplot(nrows, ncols, sb(3).idx)
yyaxis left
set(gca, 'YColor','k')
yyaxis right
h = histogram(temp,50);
h.EdgeColor = Color('teal');
h.FaceColor = Color('teal');
ylabel('Count (five minute periods)')
xlabel('Temp (\circC)')
fig = formatFig(fig, true,[nrows,ncols], sb);

save_figure(fig, ['G:\My Drive\Jeanne Lab\NOAA data\' title_str ' annual temp rate no lines'],'-png');

%% NOAA data column headers
% #    Name                           Units
% 1    WBANNO                         XXXXX
% 2    UTC_DATE                       YYYYMMDD
% 3    UTC_TIME                       HHmm
% 4    LST_DATE                       YYYYMMDD
% 5    LST_TIME                       HHmm
% 6    CRX_VN                         XXXXXX
% 7    LONGITUDE                      Decimal_degrees
% 8    LATITUDE                       Decimal_degrees
% 9    AIR_TEMPERATURE                Celsius
% 10   PRECIPITATION                  mm
% 11   SOLAR_RADIATION                W/m^2
% 12   SR_FLAG                        X
% 13   SURFACE_TEMPERATURE            Celsius
% 14   ST_TYPE                        X
% 15   ST_FLAG                        X
% 16   RELATIVE_HUMIDITY              %
% 17   RH_FLAG                        X
% 18   SOIL_MOISTURE_5                m^3/m^3
% 19   SOIL_TEMPERATURE_5             Celsius
% 20   WETNESS                        Ohms
% 21   WET_FLAG                       X
% 22   WIND_1_5                       m/s
% 23   WIND_FLAG                      X


%% Convert data from excel to .mat files
baseFolder = getCloudPath;  
rootdir = [baseFolder 'NOAA Data/'];
fileNames = {'AK', 'AL', 'AR', 'AZ', 'CA', 'CO', 'FL', 'GA', 'HI',...
             'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'ME', 'MI', 'MN', 'MO', 'MS', 'MT',...
             'NC', 'ND', 'NE','NH','NM','NV','NY', 'OR','ON', 'OK', 'OH',...
             'PA','RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VA', 'WA', 'WI', 'WV', 'WY'};

for ii = 1:length(fileNames)
    fullPath = [rootdir '2021 ' fileNames{ii} ' SubHourly NOAA Data.xlsx'];
    sheets = sheetnames(fullPath);
    disp(ii)
    for sheet = 1:size(sheets,1)
        % Pull the longitude and latitude of each recording location:
        location = sheets{sheet};
        tic
        data = readmatrix(fullPath,'Sheet', location);
        toc
        save([rootdir location '.mat'],'data','location');
        clear data location
    end
end
beep

%% Pull data and start building structure

clear
baseFolder = getCloudPath;  
rootdir = [baseFolder 'NOAA Data/'];

% get list of files
list_dirs = dir([rootdir, '/*.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
nCity = length(list_dirs);

% load data
[longitude, latitude, city_name] = deal([]);

for city = 1:nCity
    
    load([rootdir list_dirs{city}]);
    % log city name
    city_name{city} = location;
    
    % assign data
    longitude(city) = data(1,7);
    latitude(city) = data(1,8);
    cityData(city).T = data;
    
    clear location data
end
save([rootdir 'Processed Data/2021 Subhourly Data'],'-v7.3');

%% ANALYSIS find the temp rates for each sampling location
dirPath = getCloudPath;
rootdir = [dirPath, 'NOAA data\'];
load([rootdir 'Processed Data\2021 Subhourly Data.mat'])
%   save([rootdir 'Processed Data\2021 Subhourly Data.mat'],'-v7.3')

latitude = [];
longitude = [];
annualT = [];
X = linspace(-1.5,1.5,100);
[dataCDF, dataPDF] = deal([]);

for ii = 1 : nCity

    data = cityData(ii).T; % Select data for a particular city:
    if size(data,1)<105120
        dataPDF(:,ii) = nan(1,length(X));
        dataCDF(:,ii) = nan(1,length(X));
        [mu(ii),sigma(ii),annualT(ii)] = deal(nan);
        continue
    end
    % Extract data from structure
    day = data(:,4);
    time = data(:,5);
    temp = data(:,9);

    % Remove data outliers:
    loc = temp>100 | temp<-50;
    temp(loc) = nan;
    % save avg annual temp
    annualT(ii) = mean(temp,'omitnan');

    % find the month timepoints
    dayStart = find(time==5);
    fiveminrate = (diff(temp))/5;
    outliers = fiveminrate>30 | fiveminrate<-30;
    exOut = fiveminrate(outliers);
    fiveminrate(outliers) = nan;
    cityData(ii).tempRate = fiveminrate;

    % Fit data
    mu(ii) = mean(fiveminrate,'omitnan');
    sigma(ii) = var(fiveminrate,'omitnan');
    gm = gmdistribution(mu(ii),sigma(ii));
    dataPDF(:,ii) = pdf(gm,X');
    dataCDF(:,ii) = cdf(gm,X');
    
    % Pull long / lat data again
    longitude(ii) = data(1,7);
    latitude(ii) = data(1,8);

end

clear loc ii gm fiveminrate day time temp city exOut outliers

initial_vars = who;
initial_vars{end+1} = 'initial_vars';

%% ANALYSIS: clean and process temperature and light data
minT = -90; % lowest recorded temperature on record
maxT = 57; %highest recorded temp on record
light_err = -99999;

[lightData, surftempData, airtempData] = deal([]);
for ii = 1:nCity
    
    % Solar radiation power spectrum
    light = cityData(ii).T(:,11); % solar radiation data
    %clean light data
    light(light == light_err) = nan;
    light_filled = fillmissing(light, 'linear');
    lightData = autoCat(lightData, light_filled, false);
    cityData(ii).light = light_filled; % resave the data to the og structure

    % Clean up the surface temperature data
    temperature = cityData(ii).T(:,13);
    error_loc = find(temperature>=maxT | temperature<=minT);
    temperature(error_loc) = nan;
    temperature_filled = fillmissing(temperature, 'linear');
    cityData(ii).surfTemp = temperature_filled;
    surftempData = autoCat(surftempData, temperature_filled, false);

    % Clean up the surface temperature data
    temperature = cityData(ii).T(:,9);
    error_loc = find(temperature>=maxT | temperature<=minT);
    temperature(error_loc) = nan;
    temperature_filled = fillmissing(temperature, 'linear');
    cityData(ii).airTemp = temperature_filled;
    airtempData = autoCat(airtempData, temperature_filled, false);

end

initial_vars{end+1} = airtempData;

%   save([rootdir 'Processed Data\2021 Subhourly Data.mat'],'-v7.3')



%% START BY LOADING DATA HERE: 
dirPath = getCloudPath;
rootdir = [dirPath, 'NOAA data/'];
load([rootdir 'Processed Data/2021 Subhourly Data.mat'])


%% ANALYSIS: fourier transform of light / temp
% note: only use trials with FULL data set
T = 5 * 60; % Sampling period in seconds
Fs = 1/T;    % Sampling frequency        

% surface temperature
surf_T = surftempData;
remove_loc = (sum(isnan(surftempData),1)>0);
surf_T(:,remove_loc) = [];

% Fourier transform on surface temperature data
L = size(surf_T,1); % Length of signal
f = Fs*(0:(L/2))/L;
Y = fft(surf_T);
P2 = abs(Y./L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2.*P1(2:end-1,:);
surfT_f = f;
surfT_P1 = mean(P1,2);

% % FIGURE DEMONSTRATING THE ALIGNMENT ACROSS DIFFERENT REGIONS
% fig = getfig('',1); hold on
% for ii = 1:size(P1,2)
%     plot(f, P1(:,ii),'linewidth', 0.25)
% end
% formatFig(fig, true);
% grid on
% set(gca, 'TickDir','out')
% ylabel('Power')
% xlabel('Frequency (Hz)')
% title('Surface Temperature','color', 'w')
% save_figure(fig, [rootdir, 'Figures\surface temp full FFT'],'-png',false,false)
% xlim([0,1e-4])
% ylim([0,20])
% save_figure(fig, [rootdir, 'Figures\surface temp zoom 1 FFT'],'-png',false,false)
% xlim([0,0.3e-4])
% ylim([0,2])
% save_figure(fig, [rootdir, 'Figures\surface temp zoom 2 FFT'],'-png',false,false)
% xlim([0.65e-4,0.75e-4])
% ylim([0,0.5])
% save_figure(fig, [rootdir, 'Figures\surface temp zoom 3 FFT'],'-png',false,false)    

fig = getfig('',1); hold on 
plot(surfT_f, surfT_P1,'color', Color('cyan'),'linewidth', 1)
formatFig(fig, true);
grid on
set(gca, 'TickDir','out')
ylabel('Power')
xlabel('Frequency (Hz)')
title('Surface Temperature','color', 'w')

save_figure(fig, [rootdir, 'Figures\avg surface temp full FFT'],'-png',false,false)
xlim([0,1e-4])
ylim([0,10])
save_figure(fig, [rootdir, 'Figures\avg surface temp zoom 1 FFT'],'-png',false,false)
ylim([0,1])
save_figure(fig, [rootdir, 'Figures\avg surface temp zoom 2 FFT'],'-png',false,false)

    L = length(temperature_filled);  % Length of signal
    f = Fs*(0:(L/2))/L;
    Y = fft(temperature_filled);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f_temperature = f;
    P1_temperature = P1;



%% Plot geographical locations of all the cities in the data set
fig_type = '-png';
blackbackground = true;
if blackbackground
    fColor = 'k';
    bColor = 'w';
else
    fColor = 'w';
    bColor = 'k';
end

% make new data matrix with relavant information
data = [annualT',sigma',latitude',longitude',dataPDF'];
data1 = sortrows(data,3);

% remove not complete data
loc = isnan(data1(:,1));
data1(loc,:) = [];

% color points by latitude (N-S axis)
% CList = Color('indigo','lavender', nCity);
% CList(loc,:) = [];
CList = Color('plum','indigo', nCity);
CList(loc,:) = [];


% Plot locations
fig = figure;
geoscatter(data1(:,3), data1(:,4), 36, CList,'filled','^')
geobasemap colorterrain
set(fig, 'color', fColor); 

ax = gca;
set(ax, 'fontsize', 15, 'grid', 'off')

ax.LongitudeLabel.Color = bColor;
ax.LatitudeLabel.Color = bColor;
ax.LongitudeAxis.Color = bColor;
ax.LatitudeAxis.Color = bColor;

% ax.LongitudeLimits = [-175.3770 -62.7130];
% ax.LatitudeLimits = [18.2455 72.6145];

save_figure(fig,[rootdir 'Figures/Sample locations purple'],fig_type);

%% FIGURE: temp rate PDFs for all NOAA stations

% Plot the distribution functions overlaid


% fig = figure; hold on; %set(fig,'pos',[-655 564 442 662]); hold on
fig = getfig('', 1); hold on
    temp = data1(:,5:end)';
    for ii = size(temp,2):-1:1
        plot(X, temp(:,ii),'color', Color('grey'),'linewidth',0.5) %CList(ii,:) ; Color('grey')
    end
    xlabel('\DeltaT (\circC/min)')
    ylabel('PDF')
    formatFig(fig, true);
    set(gca, 'fontsize', 18,'tickdir','out')

%     v_line([-0.5,0.5],'lime','-',1)
%     v_line([-0.1,0.1],'lime','-',1)
%     v_line([-0.5,0.5],'orange',':',1)
%     v_line([-0.25,0.25],'darkviolet',':',1)
%     v_line([-0.1,0.1],'turquoise',':',1)

xlim([-1,1])

save_figure(fig,[rootdir 'Figures/All locations temp-rate zoomed in'],'-png'); % rate lines

save_figure(fig,[rootdir 'Figures/All locations temp-rate PDF'],'-pdf'); % rate lines

%% FIGURE:  plot histogram of all the change rates...

allRates = [];
for city = 1:nCity
    allRates = [allRates; cityData(city).tempRate];
end

fig = figure;
h = histogram(allRates,100);
xlim([-1,1])

%% FIGURE: correlation of temp and light 
% over the year for each city organized by longitude
all_corr = [];
for ii = 1 : nCity

    data = cityData(ii).T; 
    y = data(:,[9,11]);
    loc = y(:,1)>150 | y(:,1)<-150;
    y(loc,:) = []; % remove temp outliers
    loc = any(isnan(y),2);
    y(loc,:) = []; %remove nans
    rho = corr(y);
    all_corr(ii) = rho(1,2);
end

fig = figure;
scatter(latitude,all_corr,40,'w')%longitude
xlabel('Latitude')
ylabel('Correlation coefficient')
formatFig(fig,true);
xlim([10,80])
save_figure(fig,[rootdir 'Figures/Temp-light correlation across cities'],'-png');
  
%% POWER ANALYSIS: 
% Fourier Transform to look at the power of frequencies of light and temperature
T = 5 * 60; % Sampling period in seconds
Fs = 1/T;    % Sampling frequency        

DaysPerMonth = [1 31 28 31 30 31 30 31 31 30 31 30];
MonthStart = dayStart(cumsum(DaysPerMonth));
monthNames = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
minT = -90; % lowest recorded temperature on record
maxT = 57; %highest recorded temp on record

for ii = 1: nCity
    % Solar radiation power spectrum
    light = cityData(ii).T(:,11); % solar radiation data
    %clean light data
    light(light==-99999) = nan;
    light_filled = fillmissing(light, 'linear');
    L = length(light_filled); % Length of signal
    f = Fs*(0:(L/2))/L;
    Y = fft(light_filled);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f_light = f;
    P1_light = P1;

    % Clean up the temperature data
    temperature = cityData(ii).T(:,9);
    error_loc = find(temperature>=maxT | temperature<=minT);
    temperature(error_loc) = nan;
    temperature_filled = fillmissing(temperature, 'linear');
    
  
    L = length(temperature_filled);  % Length of signal
    f = Fs*(0:(L/2))/L;
    Y = fft(temperature_filled);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f_temperature = f;
    P1_temperature = P1;

    r = 2; c = 2;
    fig = getfig('',1);
    % light levels
    subplot(r,c,1)
    % title('Solar radiation')
    plot(light_filled, 'color', Color('gold'),'linewidth', 0.5)
    set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
    ylabel('Solar Radiation')   
    
    % temp levels
    subplot(r,c,3)
    % title('Temperature')
    plot(temperature_filled, 'color', Color('dodgerblue'),'linewidth', 0.5)
    set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
    ylabel('temperature (\circC)')  

    % light power spectrum
    subplot(r,c,2)
    plot(f_light,P1_light,'Color', Color('gold'));
    xlim([0,1e-4])
    xlabel("Frequency (Hz)")
    ylabel("Power")

    % temperature power spectrum
    subplot(r,c,4)
    plot(f_temperature,P1_temperature,'Color', Color('dodgerblue'));
    xlim([0,1e-4])
    xlabel("Frequency (Hz)")
    ylabel("Power")

    figure; 
    plot(f_light,P1_light,'Color', Color('gold'));
    hold on
    plot(f_temperature, P1_temperature)
    title("Single-Sided Amplitude Spectrum of solar radiation")
    xlabel("Frequency (Hz)")
    ylabel("Power")
    grid on;
    
end 

%% POWER ANALYSIS: 
% Fourier Transform to look at the power of frequencies of light and temperature
T = 5 * 60; % Sampling period in seconds
Fs = 1/T;    % Sampling frequency        

DaysPerMonth = [1 31 28 31 30 31 30 31 31 30 31 30];
MonthStart = dayStart(cumsum(DaysPerMonth));
monthNames = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
minT = -90; % lowest recorded temperature on record
maxT = 57; %highest recorded temp on record

% Pull all the light and temperature data into a giant matrix
[lightraw, tempraw] = deal([]);



for ii = 1: nCity
    % Solar radiation power spectrum
    light = cityData(ii).T(:,11); % solar radiation data
    %clean light data
    light(light==-99999) = nan;
    light_filled = fillmissing(light, 'linear');
    L = length(light_filled); % Length of signal
    f = Fs*(0:(L/2))/L;
    Y = fft(light_filled);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f_light = f;
    P1_light = P1;

    % Clean up the temperature data
    temperature = cityData(ii).T(:,9);
    error_loc = find(temperature>=maxT | temperature<=minT);
    temperature(error_loc) = nan;
    temperature_filled = fillmissing(temperature, 'linear');
    
  
    L = length(temperature_filled);  % Length of signal
    f = Fs*(0:(L/2))/L;
    Y = fft(temperature_filled);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f_temperature = f;
    P1_temperature = P1;

    r = 2; c = 2;
    fig = getfig('',1);
    % light levels
    subplot(r,c,1)
    % title('Solar radiation')
    plot(light_filled, 'color', Color('gold'),'linewidth', 0.5)
    set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
    ylabel('Solar Radiation')   
    
    % temp levels
    subplot(r,c,3)
    % title('Temperature')
    plot(temperature_filled, 'color', Color('dodgerblue'),'linewidth', 0.5)
    set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
    ylabel('temperature (\circC)')  

    % light power spectrum
    subplot(r,c,2)
    plot(f_light,P1_light,'Color', Color('gold'));
    xlim([0,1e-4])
    xlabel("Frequency (Hz)")
    ylabel("Power")

    % temperature power spectrum
    subplot(r,c,4)
    plot(f_temperature,P1_temperature,'Color', Color('dodgerblue'));
    xlim([0,1e-4])
    xlabel("Frequency (Hz)")
    ylabel("Power")

    figure; 
    plot(f_light,P1_light,'Color', Color('gold'));
    hold on
    plot(f_temperature, P1_temperature)
    title("Single-Sided Amplitude Spectrum of solar radiation")
    xlabel("Frequency (Hz)")
    ylabel("Power")
    grid on;
    
end 

%% FIGURE: temp over a year (single plot example)
DaysPerMonth = [1 31 28 31 30 31 30 31 31 30 31 30];
MonthStart = dayStart(cumsum(DaysPerMonth));
monthNames = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

minT = -90; % lowest recorded temperature on record
maxT = 57; %highest recorded temp on record

city = 105; %Old Town ME
data = cityData(city).T;
time = data(:,5);
surf_temp = data(:,13);
air_temp = data(:,9);

% clean data
loc = surf_temp<minT | surf_temp>maxT;
surf_temp(loc) = nan;
loc = air_temp<minT | air_temp>maxT;
air_temp(loc) = nan;

fig = figure; set(fig,'pos',[-1024 265 944 306])
plot(surf_temp,'linewidth', 2, 'Color',Color('teal'))
set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
ylabel('Temperature (\circC)')
title_str = [city_name{city} ' 2021'];
title(title_str)
xlabel('Month')
axis tight
formatFig(fig,true);
set(gca, 'fontsize', 15)
save_figure(fig,[rootdir 'Figures/Annual temp in ' city_name{city}],'-png');

%% FIGURE: temp rate over a for a year (single plot example)

fig = figure;set(fig,'pos',[-1024 265 944 306])
plot(cityData(city).tempRate,'linewidth', 2, 'Color',Color('grey'))%Color('orangered')
set(gca,'XTick',MonthStart,'XTickLabel',monthNames)
ylabel('\DeltaT (\circC/min)')
title_str = [city_name{city} ' 2021'];
title(title_str)
xlabel('Month')
axis tight

formatFig(fig,true);
set(gca, 'fontsize', 15)
ylim([-1.5,1.5])
save_figure(fig,[rootdir 'Figures/Annual temprate in ' city_name{city}],'-png');

%% FIGURE: temp rate histogram over a year (sinple plot example)
xMin = -0.5;
xMax = 0.5;
nbins = 400; 

fig = figure;
h = histogram(cityData(city).tempRate,nbins);
xlim([xMin,xMax])
h.FaceColor = 'w';
h.EdgeColor = 'w';
xlabel('\DeltaT (\circC/min)')
ylabel('Count')
title_str = [city_name{city} ' 2021'];
title(title_str)

formatFig(fig,true);
set(gca, 'fontsize', 15)
save_figure(fig,[rootdir 'Figures/Annual temprate histogram ' title_str],'-png');


% overlay the PDF fuction for this dataset:
X = linspace(xMin,xMax,100);
dist = gmdistribution(mu(city),sigma(city));
Y = pdf(dist,X');

fig = figure;
    h = histogram(cityData(city).tempRate,nbins);
    xlim([-0.5,0.5])
    h.FaceColor = 'w';
    h.EdgeColor = 'w';
    xlabel('\DeltaT (\circC/min)')
    ylabel('Count')
yyaxis right
    plot(X,Y,'color', Color('dodgerblue'),'linewidth', 2)
    ylabel('PDF')

formatFig(fig,true);
set(gca, 'fontsize', 15)
yyaxis left
set(gca, 'YColor', 'w')
yyaxis right
set(gca, 'YColor', Color('dodgerblue'))

save_figure(fig,[rootdir 'Figures/Annual temp-rate histogram with PDF ' title_str],'-png');

% POSTER CODE
fig = figure;
    h = histogram(cityData(city).tempRate,nbins);
    xlim([-0.25,0.25])
    h.FaceColor = 'k';
    h.EdgeColor = 'k';
    xlabel('\DeltaT (\circC/min)')
    ylabel('Count')
yyaxis right
    plot(X,Y,'color', Color('dodgerblue'),'linewidth', 2)
    ylabel('PDF')

formatFig(fig,false);
set(gca, 'fontsize', 18)
yyaxis left
set(gca, 'YColor', 'k')
yyaxis right
set(gca, 'YColor', Color('dodgerblue'))

save_figure(fig,[rootdir 'Figures/Annual temp-rate histogram with PDF white background ' title_str],'-png');

%% FIGURE: show city sample location on the map (single plot example)
fig = figure; 
geoscatter(cityData(city).T(1,8), cityData(city).T(1,7), 75, Color('blue'),'filled','^');
geobasemap colorterrain
set(fig, 'color', 'k')
ax = gca;
set(ax, 'fontsize', 13, 'grid', 'off')
ax.LongitudeLabel.Color = 'w';
ax.LatitudeLabel.Color = 'w';
ax.LongitudeAxis.Color = 'w'; 
ax.LatitudeAxis.Color = 'w'; 
title(city_name{city},'color', 'w')
% 
% ax.LongitudeLimits = [-175.3770 -62.7130];
% ax.LatitudeLimits = [18.2455 72.6145];

save_figure(fig,[rootdir 'Figures/Sampling location for ' city_name{city}],'-png');

%% FIGURE: comparison of annual temp and temp rate variability across the US by latitude
latitude = [];
longitude = [];
for ii = 1:nCity
    latitude(ii) = cityData(ii).T(1,8);
    longitude(ii) = cityData(ii).T(1,7);
end

data = [annualT',sigma',latitude'];
data1 = sortrows(data,3);

% color points by latitude (N-S axis)
% CList = Color('indigo','lavender', nCity);
CList = Color('plum','indigo', nCity);

fig = figure;
scatter(data1(:,1),data1(:,2),75,CList,'filled')
ylim([0,0.02])
xlabel('average annual temperature (\circC)')
ylabel('temp rate variance')
formatFig(fig,false);
set(gca, 'fontsize', 18)

save_figure(fig,[rootdir 'Figures/Annual temp vs temp rate variance white background.pdf'],'-pdf');


% Linear regression:

mdl = fitlm(data1(:,1),data1(:,2));

fig = figure;
scatter(data1(:,3), data1(:,1), 75, CList, 'filled')
xlabel('Latitude (\circ)')
ylabel('average annual temperature (\circC)')
formatFig(fig,true);
set(gca, 'fontsize', 18)

save_figure(fig,[rootdir 'Figures/Annual temp vs latitude'],'-png');

%% FIGURE: temp rate outliers and their locations
loc = find(sigma>0.02);
lon = longitude(loc);
lat = latitude(loc);

fig = figure;
    g = geoscatter(latitude, longitude, 75, Color('blue'),'filled','^');
    hold on
    g = geoscatter(lat, lon, 75, Color('red'),'filled','^');
    geobasemap colorterrain
    set(fig, 'color', 'k')
    ax = gca;
    set(ax, 'fontsize', 13, 'grid', 'off')
    ax.LongitudeLabel.Color = 'w';
    ax.LatitudeLabel.Color = 'w';
    ax.LongitudeAxis.Color = 'w'; 
    ax.LatitudeAxis.Color = 'w'; 
    title('Temp rate outlier locations','color', 'w')

save_figure(fig,[rootdir 'Figures/Sampling location for temprate outliers'],'-png');

% scatter plot overlays with the outliers color coded
fig = figure;
hold on
scatter(data1(:,3), data1(:,1), 75, CList, 'filled')
scatter(latitude(loc), annualT(loc), 75, 'r', 'filled')
xlabel('Latitude (\circ)')
ylabel('average annual temperature (\circC)')
formatFig(fig,true);
set(gca, 'fontsize', 18)

save_figure(fig,[rootdir 'Figures/Annual temp vs latitude temprate outliers'],'-png');

fig = figure; hold on
scatter(data1(:,1),data1(:,2),75,CList,'filled')
scatter(annualT(loc), sigma(loc), 75, 'r', 'filled')
xlabel('average annual temperature (\circC)')
ylabel('temp rate variance')
formatFig(fig,true);
set(gca, 'fontsize', 18)

save_figure(fig,[rootdir 'Figures/Annual temp vs temp rate variance temprate outliers'],'-png');

%% How does morning temperature compare to expected light levels? (e.g time vs temp)
T = [];
% find the 'morning' hours -- i.e. show am vs pm
time = cityData(city).T(:,5); % local solar time
temp = cityData(city).T(:,9); % air temp
tempRate = cityData(city).tempRate;
radiation = cityData(city).T(:,11);
loc = temp>maxT | temp<minT; % find false temperatures
temp(loc) = nan; % remove temp anomalies

amLoc = time<1200; % hours in before noon
pmLoc = time>=1200; % hours after noon

% combine time and temp:
data = [time(1:end-1),temp(1:end-1),tempRate,radiation(1:end-1)];
morningTemp = data(amLoc,:);
% eveningTemp = data(pmLoc,:);

% pull out each day's time course:
amStarts = find(time==0);
pmStarts = find(time==1200);

if amStarts(1)<pmStarts(1)
    for ii = 1:364
        T.time(:,ii) = data(amStarts(ii):pmStarts(ii),1);
        T.temp(:,ii) = data(amStarts(ii):pmStarts(ii),2);
        T.tempR(:,ii) = data(amStarts(ii):pmStarts(ii),3);
        T.radi(:,ii) = data(amStarts(ii):pmStarts(ii),4);
    end
end

% check individual patterns

% plotList = [160 161 166 167 170];
% ii = 172;
    
ii = ii+1;
x = T.time(:,ii);
x = find(rem(x,100)==0); % pull hour start times
row = 2; col = 1;
fig = figure; 
    subplot(row,col,1)
        plot(T.tempR(:,ii),'color', 'w','linewidth', 1.5)
        axis tight
        hline(0)
        ax = gca;
        set(ax, 'XTick',x,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('\DeltaT (\circC/min)')
    subplot(row,col,2); hold on
        plot(T.temp(:,ii),'color', 'w','linewidth', 1.5)
        yyaxis right
        plot(T.radi(:,ii),'color', Color('gold'),'linewidth', 1)
        axis tight
        yyaxis left
        ax = gca;
        set(ax, 'XTick',x,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('Temperature (\circC)')
formatFig(fig,true, [row,col]);
subplot(row,col,2);
yyaxis right
ax = gca;
set(ax, 'YColor',Color('gold'))
ylabel('Solar radiation (W/m^2)')


save_figure(fig,[rootdir 'Figures/Morning temperature ' city_name{city} ' day ' num2str(ii)],'-png');

%% 
T = [];
% find the 'morning' hours -- i.e. show am vs pm
time = cityData(city).T(:,5); % local solar time
temp = cityData(city).T(:,9); % air temp
tempRate = cityData(city).tempRate;
radiation = cityData(city).T(:,11);
loc = temp>maxT | temp<minT; % find false temperatures
temp(loc) = nan; % remove temp anomalies

amLoc = time<1200; % hours in before noon
pmLoc = time>=1200; % hours after noon

% pull out each day's time course:
amStarts = find(time==0);
% pmStarts = find(time==1200);

% combine time and temp:
data = [time(1:end-1),temp(1:end-1),tempRate,radiation(1:end-1)];
for ii = 1:364
    T.time(:,ii) = data(amStarts(ii):amStarts(ii+1)-1,1);
    T.temp(:,ii) = data(amStarts(ii):amStarts(ii+1)-1,2);
    T.tempR(:,ii) = data(amStarts(ii):amStarts(ii+1)-1,3);
    T.radi(:,ii) = data(amStarts(ii):amStarts(ii+1)-1,4);
end



% Radiation vs temperature slowly add new lines
ii_range = 190:199;

n = length(T.radi(:,ii));
colorList = [Color('Blue','White',144); Color('White','Red',144)];
for kk = 1:length(ii_range)
    fig = figure; set(fig, 'pos', [-1072 430 1033 712]); hold on
    for ii = ii_range(1):ii_range(kk)
%         scatter(T.radi(:,ii),T.temp(:,ii),50,'filled')
        xx = [T.radi(:,ii);T.radi(1:10,ii+1)];
        yy = [T.temp(:,ii);T.temp(1:10,ii+1)];
        plot(xx,yy,'linewidth',0.5,'marker','o')
    end
    axis tight
    xlabel('Solar radiation (W/m^2)')
    ylabel('Temperature (\circC)')
    formatFig(fig, true);
    xlim([0 902])
    ylim([9.9, 26.5])
    save_figure(fig,[rootdir 'Figures/time series/radiation v temperature day ' num2str(kk)],'-png',true);
end
        
    

    
    
% Radiation vs temperature 
n = length(T.radi(:,ii));
colorList = [Color('Blue','White',73); Color('White','Red',72)];
fig = figure; hold on
    for ii = 1:15
        xx = [T.radi(:,ii);T.radi(1:10,ii+1)];
        yy = [T.temp(:,ii);T.temp(1:10,ii+1)];
        plot(xx,yy,'color', 'w')
    end
    axis tight
    xlabel('Solar radiation (W/m^2)')
    ylabel('Temperature (\circC)')
    formatFig(fig, true);

figure;
    scatter(1:length(xx),xx)
    
    
    
ii = ii+1;
% Overlay multiple days in one plot

% plotList = [160 161 166 167 170];
% CList = Color('orange', 'teal',length(plotList));

plotList = [162 163 164 165 168 169];
CList = Color('pink', 'blue',length(plotList));

x = T.time(:,plotList(1));
x = find(rem(x,100)==0); % pull hour start times
row = 2; col = 1;
fig = figure; 
idx = 0;
for ii = plotList
    idx = idx+1;
    subplot(row,col,1)
    hold on
        plot(T.tempR(:,ii),'color', CList(idx,:),'linewidth', 1.5)
    subplot(row,col,2)
    hold on
        plot(T.temp(:,ii),'color', CList(idx,:),'linewidth', 1.5)
    lStr{idx} = num2str(ii);
end
% formatting:
 subplot(row,col,1)
        axis tight
        hline(0)
        ax = gca;
        set(ax, 'XTick',x,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('\DeltaT (\circC/min)')
    subplot(row,col,2)
        axis tight
        ax = gca;
        set(ax, 'XTick',x,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('Temperature (\circC)')
formatFig(fig,true, [row,col]);
legend(lStr,'textcolor', 'w','box', 'off')

save_figure(fig,[rootdir 'Figures/Morning temperature ' city_name{city} ' temp rise'],'-png');

%% FIGURE:  overlay the 'up' and 'down' morning temp behaviors

% find the 'morning' hours -- i.e. show am vs pm
data = [];
time = cityData(city).T(:,5); % local solar time
temp = cityData(city).T(:,9); % air temp
tempRate = cityData(city).tempRate;
loc = temp>maxT | temp<minT; % find false temperatures
temp(loc) = nan; % remove temp anomalies

amLoc = find(time<1200); % hours in before noon
pmLoc = find(time>=1200); % hours after noon

% combine time and temp:
data = [time(1:end-1),temp(1:end-1),tempRate];
morningTemp = data(amLoc,:);
% eveningTemp = data(pmLoc,:);

% pull out each day's time course:
amStarts = find(time==0);
pmStarts = find(time==1200);

if amStarts(1)<pmStarts(1)
    for ii = 1:364
        T.time(:,ii) = data(amStarts(ii):pmStarts(ii),1);
        T.temp(:,ii) = data(amStarts(ii):pmStarts(ii),2);
        T.tempR(:,ii) = data(amStarts(ii):pmStarts(ii),3);
        T.solar(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),11);
        T.wind(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),22);
        T.surfT(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),13);
        T.rain(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),10);
        T.precip(:,ii) = sum(cityData(city).T(amStarts(ii):amStarts(ii+1),10));
    end
end

% Overlay multiple days in one plot
x = T.time(:,plotList(1));
x = find(rem(x,100)==0); % pull hour start times
row = 2; col = 1;
fig = figure; 

for type = 1:2
    switch type
        case 1
            plotList = [160 161 166 167 170];
            kolor = Color('orange');
        case 2
            plotList = [162 163 164 165 168 169];
            kolor = Color('dodgerblue');
    end
    idx = 0;
    for ii = plotList
        idx = idx+1;
        subplot(row,col,1)
        hold on
            plot(T.tempR(:,ii),'color',kolor,'linewidth', 1.5)
        subplot(row,col,2)
        hold on
            plot(T.temp(:,ii),'color', kolor,'linewidth', 1.5)
        lStr{idx} = num2str(ii);
    end
end
% formatting:
 subplot(row,col,1)
        axis tight
        hline(0)
        ax = gca;
        set(ax, 'XTick',x ,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('\DeltaT (\circC/min)')
    subplot(row,col,2)
        axis tight
        ax = gca;
        set(ax, 'XTick',x ,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('Temperature (\circC)')
formatFig(fig,true, [row,col]);

save_figure(fig,[rootdir 'Figures/Morning temperature ' city_name{city} ' temp overlays'],'-png');

% What other data is diff / the same about these days' groupings??
% Overlay multiple days in one plot
x = T.time(:,plotList(1));
x = find(rem(x,100)==0); % pull hour start times
row = 2; col = 3;
fig = figure; 
for type = 1:2
    switch type
        case 1
            plotList = [160 161 166 167 170];
            kolor = Color('orange');
        case 2
            plotList = [162 163 164 165 168 169];
            kolor = Color('dodgerblue');
    end
    for ii = plotList
        % 1 temp rate
        subplot(row,col,1)
        hold on
            plot(T.tempR(:,ii),'color',kolor,'linewidth', 1.5)
            ylabel('\DeltaT (\circC/min)')
        % 4 air temp
        subplot(row,col,4)
        hold on
            plot(T.temp(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Temperature (\circC)')
        % 3 wind speed (22)
        subplot(row,col,3)
        hold on
            plot(T.wind(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('\Wind speed (m/s)')
        % 2 precipitation (10)
        subplot(row,col,2)
        hold on
        scatter((rand(1)*0.3)+type,T.precip(ii),50,kolor,'filled')
%             plot(T.rain(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Precipitation (mm)')
        % 6 SOLAR_RADIATION
        subplot(row,col,6)
        hold on
            plot(T.solar(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Solar radiation (W/m^2)')
        % 5 SURFACE_TEMPERATURE
        subplot(row,col,5)
        hold on
            plot(T.surfT(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Surface temp (\circC)')
    end
end
% formatting:
for ii = 1:6
 subplot(row,col,ii)
        axis tight
        ax = gca;
        set(ax, 'XTick',x ,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
end
subplot(row,col,2)
xlim([0.5,2.5])
xlabel('')
set(gca,'XTickLabel',[])
formatFig(fig,true, [row,col]);
save_figure(fig,[rootdir 'Figures/Morning temperature ' city_name{city} ' multifactor overlays'],'-png');

%% TRIAL AND ERROR: try looking at temp patterns by precipitation

% find the 'morning' hours -- i.e. show am vs pm
data = [];
time = cityData(city).T(:,5); % local solar time
temp = cityData(city).T(:,9); % air temp
tempRate = cityData(city).tempRate;
loc = temp>maxT | temp<minT; % find false temperatures
temp(loc) = nan; % remove temp anomalies

amLoc = find(time<1200); % hours in before noon
pmLoc = find(time>=1200); % hours after noon

% combine time and temp:
data = [time(1:end-1),temp(1:end-1),tempRate];
morningTemp = data(amLoc,:);
% eveningTemp = data(pmLoc,:);

% pull out each day's time course:
amStarts = find(time==0);
pmStarts = find(time==1200);

if amStarts(1)<pmStarts(1)
    for ii = 1:364
        T.time(:,ii) = data(amStarts(ii):pmStarts(ii),1);
        T.temp(:,ii) = data(amStarts(ii):pmStarts(ii),2);
        T.tempR(:,ii) = data(amStarts(ii):pmStarts(ii),3);
        T.solar(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),11);
        T.wind(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),22);
        T.surfT(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),13);
        T.rain(:,ii) = cityData(city).T(amStarts(ii):pmStarts(ii),10);
        T.precip(:,ii) = sum(cityData(city).T(amStarts(ii):amStarts(ii+1),10));
    end
end


% Overlay multiple days in one plot
x = T.time(:,plotList(1));
x = find(rem(x,100)==0); % pull hour start times
row = 2; col = 1;
fig = figure; 


for ii = 125:200
    if T.precip(ii)>1
        type = 2;
    else type = 1;
    end
    switch type
        case 1
            plotList = [160 161 166 167 170];
            kolor = Color('orange');
        case 2
            plotList = [162 163 164 165 168 169];
            kolor = Color('dodgerblue');
    end
    idx = 0;
        idx = idx+1;
        subplot(row,col,1)
        hold on
            plot(T.tempR(:,ii),'color',kolor,'linewidth', 1.5)
        subplot(row,col,2)
        hold on
            plot(T.temp(:,ii),'color', kolor,'linewidth', 1.5)
        lStr{idx} = num2str(ii);
end
% formatting:
 subplot(row,col,1)
        axis tight
        hline(0)
        ax = gca;
        set(ax, 'XTick',x ,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('\DeltaT (\circC/min)')
    subplot(row,col,2)
        axis tight
        ax = gca;
        set(ax, 'XTick',x ,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
        ylabel('Temperature (\circC)')
formatFig(fig,true, [row,col]);

save_figure(fig,[rootdir 'Figures/Morning temperature ' city_name{city} ' temp overlays'],'-png');

% What other data is diff / the same about these days' groupings??
% Overlay multiple days in one plot
x = T.time(:,plotList(1));
x = find(rem(x,100)==0); % pull hour start times
row = 2; col = 3;
fig = figure; 
for type = 1:2
    switch type
        case 1
            plotList = [160 161 166 167 170];
            kolor = Color('orange');
        case 2
            plotList = [162 163 164 165 168 169];
            kolor = Color('dodgerblue');
    end
    for ii = plotList
        % 1 temp rate
        subplot(row,col,1)
        hold on
            plot(T.tempR(:,ii),'color',kolor,'linewidth', 1.5)
            ylabel('\DeltaT (\circC/min)')
        % 4 air temp
        subplot(row,col,4)
        hold on
            plot(T.temp(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Temperature (\circC)')
        % 3 wind speed (22)
        subplot(row,col,3)
        hold on
            plot(T.wind(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('\Wind speed (m/s)')
        % 2 precipitation (10)
        subplot(row,col,2)
        hold on
        scatter((rand(1)*0.3)+type,T.precip(ii),50,kolor,'filled')
%             plot(T.rain(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Precipitation (mm)')
        % 6 SOLAR_RADIATION
        subplot(row,col,6)
        hold on
            plot(T.solar(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Solar radiation (W/m^2)')
        % 5 SURFACE_TEMPERATURE
        subplot(row,col,5)
        hold on
            plot(T.surfT(:,ii),'color', kolor,'linewidth', 1.5)
            ylabel('Surface temp (\circC)')
    end
end
% formatting:
for ii = 1:6
 subplot(row,col,ii)
        axis tight
        ax = gca;
        set(ax, 'XTick',x ,'XTickLabels',T.time(x,ii)/100)
        xlabel('Time (hr)')
end
subplot(row,col,2)
xlim([0.5,2.5])
xlabel('')
set(gca,'XTickLabel',[])
formatFig(fig,true, [row,col]);
save_figure(fig,[rootdir 'Figures/Morning temperature ' city_name{city} ' multifactor overlays'],'-png');

%% How much into the future could the flies predict from knowing current temp and rate of temp change?

%% Light levels vs temperature for a full day
city = 5;
T = [];
% find the 'morning' hours -- i.e. show am vs pm
time = cityData(city).T(:,5); % local solar time
temp = cityData(city).T(:,9); % air temp
tempRate = cityData(city).tempRate;
radiation = cityData(city).T(:,11);
loc = temp>maxT | temp<minT; % find false temperatures
temp(loc) = nan; % remove temp anomalies

amLoc = time<1200; % hours in before noon
pmLoc = time>=1200; % hours after noon

% combine time and temp:
data = [time(1:end-1),temp(1:end-1),tempRate,radiation(1:end-1)];
morningTemp = data(amLoc,:);
% eveningTemp = data(pmLoc,:);

% pull out each day's time course:
amStarts = find(time==0);
% pmStarts = find(time==1200);

for ii = 1:364
    T.time(:,ii) = data(amStarts(ii):amStarts(ii+1),1);
    T.temp(:,ii) = data(amStarts(ii):amStarts(ii+1),2);
    T.tempR(:,ii) = data(amStarts(ii):amStarts(ii+1),3);
    T.radi(:,ii) = data(amStarts(ii):amStarts(ii+1),4);
    % change in solar radiation calculation
    changeT = (diff(T.time(:,ii)));
    changeT(changeT==45) = 5; %replace end of hour time gaps (ie 55->0 loop)
    changeT(end) = 5; % set last time gap to same rate
    if any(changeT>5 | changeT<5)
        warn('missing time data points')
        break
    end
    T.radiR(:,ii) = (diff(T.radi(:,ii)))./changeT;
end

% check individual patterns below

%% 
% plotList = [169 177 190 195 202];
ii = 202;
gap = 3;

x = T.time(:,ii);
x = find(rem(x,100)==0); % pull hour start times
row = 2; col = 1;
fig = figure; set(fig, 'pos', [-767 427 469 733])
% Rate of change
    subplot(row,col,1); hold on
        yyaxis left
        plot(T.tempR(:,ii),'color', Color('dodgerblue'),'linewidth', 1.5)
        ylabel('\DeltaT/t (\circC/min)')
        yyaxis right
        plot(T.radiR(:,ii),'color', Color('orange'),'linewidth', 1.5)
        ylabel('\Delta in radiation/minute')
        
        axis tight
        ax = gca;
        set(ax, 'XTick',x(1:gap:end),'XTickLabels',T.time(x(1:gap:end),ii)/100)
        xlabel('Time (hr)')
        
% Absolute values
    subplot(row,col,2); hold on
        plot(T.temp(:,ii),'color', Color('dodgerblue'),'linewidth', 1.5)
        ylabel('Temperature (\circC)')
        yyaxis right
        plot(T.radi(:,ii),'color', Color('orange'),'linewidth', 1)
        ylabel('Solar radiation (W/m^2)')
        axis tight
        ax = gca;
        set(ax, 'XTick',x(1:gap:end),'XTickLabels',T.time(x(1:gap:end),ii)/100)
        xlabel('Time (hr)')
% FORMATTING
formatFig(fig,false, [row,col]);
subplot(row,col,1);
yyaxis left
set(gca, 'YColor',Color('dodgerblue'))
yyaxis right
set(gca, 'YColor',Color('orange'))
subplot(row,col,2);
yyaxis left
set(gca, 'YColor',Color('dodgerblue'))
yyaxis right
set(gca, 'YColor',Color('orange'))


save_figure(fig,[rootdir 'Figures/Morning temperature ' city_name{city} ' day ' num2str(ii)],'-pdf');

%% FIGURE: 
% plot the temp vs solar radiation over time -- averaged across each month

monthLength = [31,28,31,30,31,30,31,31,30,31,30,29]; % shortened december by one day

startnum = 1;

fig = figure; set(fig,'pos',[-807 411 434 633])
hold on
for ii = 1:12

    roi = startnum : (startnum + monthLength(ii));
    startnum = startnum + monthLength(ii);

    xx = T.radi(:,roi);
    yy = T.temp(:,roi);
%     yy(yy<0) = nan; % screen outliers (clear errors)
%     xx(xx<0) = nan; % screen outliers (clear errors)
    % for ii = 1:length(roi)
    %     plot(xx(:,ii),yy(:,ii),'color',Color('grey'),'linewidth',0 5)
    % end
    plot(mean(xx,2),mean(yy,2),'linewidth',1.5)

end

% Formatting
xlabel('solar radiation (W/m^2)')
ylabel('temperature (\circC)')
formatFig(fig,false);
legend(monthNames,'box','off')

save_figure(fig,[rootdir 'Figures/solar vs temp month average for ' city_name{city}],'-pdf');

%% Long-range temp vs light overlay

dayROI = 198:204;

[radiData,tempData] = deal([]);
for ii = dayROI
  
  tempData = [tempData; T.time(:,ii), T.temp(:,ii),T.tempR(:,ii)];
  radiData = [radiData; T.time(:,ii), T.radi(:,ii),[T.radiR(:,ii);nan]];

end

tempData(:,3) = fillmissing(tempData(:,3),'movmedian',30000);  
radiData(:,3) = fillmissing(radiData(:,3),'movmedian',30000);  


x = 0:289:size(tempData,1);

% plotList = [169 177 190 195 202];
ii = 202;
gap = 3;

row = 2; col = 1;
fig = figure; set(fig, 'pos', [-951 661 891 499])
% Rate of change
    subplot(row,col,1); hold on
        plot(normalize(tempData(:,3)),'color', Color('dodgerblue'),'linewidth', 1.5)
        plot(normalize(radiData(:,3)),'color', Color('orange'),'linewidth', 1.5)  
        ylabel('rate of change z-score')
        axis tight
        set(gca, 'XTick',x,'XTickLabels',1:length(x))
        xlabel('Time (hr)')
       
% Absolute values
    subplot(row,col,2); hold on
        plot(tempData(:,2),'color', Color('dodgerblue'),'linewidth', 1.5)
        ylabel('Temperature (\circC)')
        yyaxis right
        plot(radiData(:,2),'color', Color('orange'),'linewidth', 1)
        ylabel('Solar radiation (W/m^2)')
        axis tight
        set(gca, 'XTick',x,'XTickLabels',1:length(x))
        xlabel('Time (hr)')
        
% FORMATTING
formatFig(fig,false, [row,col]);
subplot(row,col,1);
% yyaxis left
% set(gca, 'YColor',Color('dodgerblue'))
% yyaxis right
% set(gca, 'YColor',Color('orange'))
subplot(row,col,2);
yyaxis left
set(gca, 'YColor',Color('dodgerblue'))
yyaxis right
set(gca, 'YColor',Color('orange'))
save_figure(fig,[rootdir 'Figures/week long temp vs light ' city_name{city} ' days '...
    num2str(dayROI(1)) '-' num2str(dayROI(end))],'-pdf');





% FIGURE WITH CORRELATION COEFFICIENT
% Correlation coefficients...
sSpan = 80; % two hour rolling window
[corrRate,corrAbs] = deal(nan([size(tempData,1),1]));
for ii = sSpan/2:size(tempData,1)-sSpan
    roi = ii:ii+sSpan;
    R = corrcoef(tempData(roi,3),radiData(roi,3));
    corrRate(ii) = R(1,2);
    R = corrcoef(tempData(roi,2),radiData(roi,2));
    corrAbs(ii) = R(1,2);
end


x = 0:289:size(tempData,1);

% plotList = [169 177 190 195 202];
ii = 202;
gap = 3;

row = 3; col = 1;
fig = figure; set(fig, 'pos', [-1077 456 1050 704])
% Absolute values
    subplot(row,col,1); hold on
        plot(tempData(:,2),'color', Color('dodgerblue'),'linewidth', 1.5)
        ylabel('Temperature (\circC)')
        yyaxis right
        plot(radiData(:,2),'color', Color('orange'),'linewidth', 1)
        ylabel('Solar radiation (W/m^2)')
        axis tight
        xlimits = xlim;
        set(gca, 'XTick',x,'XTickLabels',1:length(x))
        xlabel('Time (hr)')
% Rate of change
    subplot(row,col,2); hold on
        plot(normalize(tempData(:,3)),'color', Color('dodgerblue'),'linewidth', 1.5)
        plot(normalize(radiData(:,3)),'color', Color('orange'),'linewidth', 1.5)  
        ylabel('rate of change z-score')
        axis tight
        xlim(xlimits)
        set(gca, 'XTick',x,'XTickLabels',1:length(x))
        xlabel('Time (hr)')
% Correlation
    subplot(row,col,3); hold on
        plot(corrAbs,'color', Color('black'),'linewidth', 1.5)
        plot(corrRate,'color', Color('grey'),'linewidth', 1.5)  
        ylabel('Correlation Coefficient')
        set(gca, 'XTick',x,'XTickLabels',1:length(x))
        xlabel('Time (hr)')
        xlim(xlimits)
% FORMATTING
formatFig(fig,false, [row,col]);
% subplot(row,col,1);
% yyaxis left
% set(gca, 'YColor',Color('dodgerblue'))
% yyaxis right
% set(gca, 'YColor',Color('orange'))
subplot(row,col,1);
yyaxis left
set(gca, 'YColor',Color('dodgerblue'))
yyaxis right
set(gca, 'YColor',Color('orange'))
save_figure(fig,[rootdir 'Figures/absolute rate and coefficient data for ' city_name{city} ' days '...
    num2str(dayROI(1)) '-' num2str(dayROI(end))],'-pdf');


%% Power analysis

workingData = fillmissing(surf_temp,'movmedian',30000);  



nSamples = length(workingData);
Fs = 0.083; %five minute bin sampling
f = 0:Fs/nSamples:Fs/2;     % frequency axis
a = fft(workingData);         % Fourier transform
P = a.*conj(a)./nSamples;   % power = square of the amplitude, normalized
S = P(1:floor(nSamples/2+1));      

% find the frequency with the most power:
transData = [f',S];
transData = sortrows(transData,2,'descend'); % first two values here are quite large
max_spec = transData(1:2,:); %pull out the large values 
disp('The frequencies the most power are: ') 
disp(max_spec(:,1))
disp('The power at these frequencies are: ') 
disp(max_spec(:,2))


% PLOT SPECTRAL FREQUENCY DATA
row = 2;
col = 1;
LW = 2;

fig = figure;
% Time domain signal (raw)
subplot(row, col, 1); hold on
plot(workingData,'color', Color('red'),'linewidth',LW)
plot(surf_temp,'color', Color('grey'),'linewidth',LW)
xlabel('time')
ylabel('temperature (\circC)')
axis tight
% Frequency domain signal
subplot(row,col,2)
semilogy(f,S,'color', Color('teal'),'linewidth',LW) %log scale
xlabel('Frequency (Hz)')
ylabel('Spectral power density')
% v_line(max_spec(:,1),'red','-.',0.5) 
axis tight

formatFig(fig, false, [row,col]);

%% Range of temperature change in a day

[B,I] = sort(latitude);

tRange = nan(366,nCity);
city_list = 1:nCity;
city_list = city_list(I);

for city = city_list
    if city==87
        continue
    end
    temp = cityData(city).T(:,9); % air temp
    time = cityData(city).T(:,5); % local solar time
    %remove temp outliers
    loc = temp(:,1)>57 | temp(:,1)<-50;
    temp(loc) = nan;   
    %find temp range per day
    dayStarts = find(time==0);
    for ii = 1:length(dayStarts)-1
        ROI = dayStarts(ii):dayStarts(ii+1);
        tRange(ii,city) = range(temp(ROI));
    end
end

CList = Color('plum','indigo', nCity);

% FIGURE
SZ = 10;
MT = ones(1,366);
fig = getfig('',true,[396 680]);  
hold on; %set(fig,'pos',[-1056 637 1010 496]);
for city = 1:nCity
    x = city*MT;
    y = tRange(:,city);
    scatter(x,y,SZ,CList(city,:))% 
end
for city = 1:nCity
    x = city*MT;
    y = tRange(:,city);
    scatter(city,mean(y,'omitnan'),SZ+10,'w','filled')
end

formatFig(fig,true);
ylabel('Daily temperature range (\circC)')
set(gca,'xtick',[],'xcolor','k','tickdir','out')
xlabel('2021 US cities','color','white')
ylim([0,40])
% optional plot nationwide average
national_avg = median(median(tRange,1,'omitnan'),'omitnan');
h_line(national_avg,'w','-',2)


save_figure(fig,[rootdir 'Figures/avg temperature range all cities'],'-pdf');


[~,city] = (max(mean(tRange,1,'omitnan')));
disp(city_name{city})

city = 49;
figure; plot(temp)

city_name{87}



