clear; clc

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

%% Format and process annual temperature data from NOAA

yearNames = {'2019', '2020', '2021', '2022'};
rootdir = 'G:\My Drive\Jeanne Lab\DATA\NOAA data\Excel Sheets AL\';
for i = 1:length(yearNames)
    dataPath = [rootdir yearNames{i} '.xlsx'];
    data = readmatrix(dataPath);
    Year = yearNames{i};
    save([rootdir yearNames{i} '.mat'],'data','Year');
    clear data Year
end
    
% get list of files
list_dirs = dir([rootdir, '/*.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
nYears = length(list_dirs);

% load data
yearData = [];
for yy = 1:nYears  
    load([rootdir list_dirs{yy}]);
    % assign data
    yearData(yy).T = data;
end
save([rootdir 'Subhourly Data'],'-v7.3');

%% 


%% FIGURE: temp over a year (single plot example)

minT = -90; % lowest recorded temperature on record
maxT = 57; %highest recorded temp on record

temp = [];
for yy = 1:nYears
    temp = [temp; yearData(yy).T(:,9)];
end
loc = temp>100 | temp<-50;
temp(loc) = nan;


buff = 1;
x_ticks = linspace(buff,length(temp)-buff,4);

fig = figure; set(fig,'pos',[-1058 666 999 481])
plot(smooth(temp,360,'moving','omitnan'),'linewidth',2,'color','w')
xlabel('Years')
set(gca,'xtick',x_ticks,'xticklabel',yearNames)
ylabel('Temperature (\circC)')
formatFig(fig,true);

save_figure(fig,[rootdir 'Annual temps'],'-png');


%% 

% minT = -90; % lowest recorded temperature on record
% maxT = 57; %highest recorded temp on record
% 
% temp = [];
% for yy = 1:nYears
%     temp = [temp; yearData(yy).T(:,9)];
% end
% loc = temp>maxT | temp<minT;
% temp(loc) = nan;
% 
% 
% buff = 1;
% x_ticks = linspace(buff,length(temp)-buff,4);
% 
% fig = figure; set(fig,'pos',[-1058 666 999 481])
% plot(smooth(temp,360,'moving','omitnan'),'linewidth',2,'color','w')
% xlabel('Years')
% set(gca,'xtick',x_ticks,'xticklabel',yearNames)
% ylabel('Temperature (\circC)')
% formatFig(fig,true);
% 
% save_figure(fig,[rootdir 'hourly temps'],'-png');



%%
latitude = [];
longitude = [];
annualT = [];
X = linspace(-1.5,1.5,100);
[dataCDF, dataPDF] = deal([]);

for ii = 1 : nCity

    data = yearData(ii).T; % Select data for a particular year:
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

    % Remove outliers:
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
    yearData(ii).tempRate = fiveminrate;

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


data = yearData(2).T(:,9); % 5 minute sampling bins

fig = figure; set(fig,'pos',[-858 675 802 487])
plot(data(1100:1120),'color','w','linewidth',2)
xlabel('time (min)')
set(gca,'xtick',[0,10,20],'xticklabels',{'0','30','60'})
xlim([0,22])
ylabel('Temperature (\circC)')
formatFig(fig,true);

save_figure(fig,[rootdir 'hour temps'],'-png');



