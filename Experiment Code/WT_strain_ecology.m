


% Berlin.min_temp = ;
% Berlin.max_temp = ;
% Berlin.temp_variance = ;
% Berlin.elevation = 43;



expList = {'Berlin WT','CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'};

% latitide organized:
expList = {'Swedish', 'Berlin WT', 'OregonR','CantonS','Malawi', 'Zimbabwe'};
latitudes = [60.1282,  52.5200,    43.8041,   40.7989, -13.2543, -19.0154];

Swedish.lat = 60.1282;
Berlin.lat = 52.5200;
OregonR.lat = 43.8041;
CantonS.lat = 40.7989;
Malawi.lat = -13.2543;
Zimbabwe.lat = -19.0154;




%% Saved geographical locations and names
clear
rootdir = 'G:\My Drive\Jeanne Lab\DATA\Ecological Data\';

stocks = [];

stocks(1).name = 'Berlin';
stocks(1).coor = [52.52, 13.4];
stocks(1).color = 'Purple';

stocks(2).name = 'CantonS';
stocks(2).coor = [40.79, -81.38];
stocks(2).color = 'Blue';

stocks(3).name = 'OregonR';
stocks(3).coor = [43.8,-120.55];
stocks(3).color = 'Green';

stocks(4).name = 'SwedishC';
stocks(4).coor = [60.13, 18.64];
stocks(4).color = 'Red';

stocks(5).name = 'Malawi';
stocks(5).coor = [-13.25, 34.3];
stocks(5).color = 'Gold';

stocks(6).name = 'Zimbabwe';
stocks(6).coor = [-19.02, 29.15];
stocks(6).color = 'White';
curr_stocks = length(stocks);

%%  potential new stocks
idx = curr_stocks+1;
for i = 1:idx-1
    stocks(i).color = 'White';
end

names = {'hikone','harwich','crimea','bogota','Zambia','lousiana','spain','smarkand', 'portugal','nyc','athens','urbana'};
coords = [35.27,136.26;...   hikone
          41.67,-70.06;...  harwich
          44.95, 34.10;...  crimea
          4.71, -74.07;...  bogota
          -16.92, 28;...    Zambia
          30.98, -91.96;... lousiana
          42.39, 1;...      spain
          39.63, 66.97;...  smarkand
          32.76, -16.96;... portugal
          40.71, -74;...    nyc
          37.98, 23.73;...  athens
          40.11, -88.21]; % urbana

for i = 1:length(names)
    stocks(idx).name = names{i};
    stocks(idx).coor = coords(i,:);
    stocks(idx).color = 'magenta';
    idx = idx+1;
end

%% Where evyn has lived...

evyn = [31.4887, -4.4015;...%JORF
        33.8730, -5.5407;...%MEKNES
        44.096, -69.97;...%BOWDOIN
        47.61, -122.33;...%SEATTLE
        41.3839, -72.9026;...%HAMDEN
        43.2965, 5.3698;...%MARSEILLES
        44.8026, -0.5881;...%TALENCE
        48.8566, 2.3522;...%PARIS
        -25.0225, 46.9854;...%FORT DAUPHIN
        -16.1706, 49.7652]; %MANANARA NORD


%% Geographical locations

blackbackground = true;
if blackbackground
    fColor = 'k';
    bColor = 'w';
else
    fColor = 'w';
    bColor = 'k';
end

% Plot locations
fig = figure; set(fig, 'pos',[534 166 956 663])
for i = 1:length(stocks)
    geoscatter(stocks(i).coor(1), stocks(i).coor(2), 36, Color(stocks(i).color),'filled','MarkerEdgeColor','black')
    hold on
end
geobasemap colorterrain
set(fig, 'color', fColor); 
% for fun
for i = 1:length(evyn)
    geoscatter(evyn(i,1), evyn(i,2), 36, Color('gold'),'filled','MarkerEdgeColor','black')
    hold on
end

ax = gca;
set(ax, 'fontsize', 15, 'grid', 'on')

ax.LongitudeLabel.Color = bColor;
ax.LatitudeLabel.Color = bColor;
ax.LongitudeAxis.Color = bColor;
ax.LatitudeAxis.Color = bColor;

% ax.LongitudeLimits = [-175.3770 -62.7130];
% ax.LatitudeLimits = [18.2455 72.6145];

save_figure(fig,[rootdir 'Figures/Genotype locations with evyn locations'],'-png');
save_figure(fig,[rootdir 'Figures/Genotype locations'],'-png');

%% 


fig = figure; set(fig, 'pos',[-843 481 292 584]); hold on
% x = shuffle_data(linspace(1,2,length(stocks)));
for i = 1:length(stocks)
    if i<=curr_stocks
        x = 1;
    else
        x = 2;
    end
    scatter(x, stocks(i).coor(1),50,stocks(i).color,'filled')
end
xlim([0.5,2.5])

% for fun
scatter(3*ones([1,length(evyn)]),evyn(:,1),50,Color('gold'),'filled')
xlim([0.5,3.5])

% serious to consider
idx = curr_stocks+[1,4,5,6,11];
for i = idx
    scatter(2, stocks(i).coor(1), 75, 'w')
end

ylabel('Latitude (\circ)')
formatFig(fig,true);
set(gca,'xcolor', 'k')
save_figure(fig,[rootdir 'Figures/Genotype latitudes'],'-png');














