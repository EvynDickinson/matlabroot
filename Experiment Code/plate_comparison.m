
%% Plate comparison 

%% Load data

% clear; warning off

%select data structure from folder names out of GroupDataGUI:
% Need both the Data Folder path (pulling data from here) and the
% structures folder (saving data to here)
baseFolder = getDataPath(5,0);
pathNames = getPathNames;


ExpGroup = selectFolder([baseFolder pathNames.grouped_trials],false,'Select data structure');
ExpGroup = ExpGroup{:};
% data structures folder:
figDir = [baseFolder pathNames.grouped_trials ExpGroup '/'];

rawDataFolder = [baseFolder, pathNames.single_trial];

% Does the data grouping already exist?
raw_file = [figDir ExpGroup ' raw.mat'];
if exist(raw_file,'file') == 2 %file exists
    oldList =  load(raw_file, 'T');
    newList = load([figDir 'fileList.mat'],'T');
end

% load the directory list:
load([figDir 'fileList.mat'])

% extract information
ntrials = size(T,1);
dates = T.Date;
arenas = T.Arena;
expID = T.ExperimentID;

% load data
n = cell(1,ntrials);
data = struct('x', n, 'y', n, 'ecent', n,'nflies',n,'wellcenters', n,'foodwell',n);  
var2load = fieldnames(data);

tic
for trial = 1:ntrials

    trial_ID = [dates{trial} '_' expID{trial} '_' arenas{trial}];
    filePath = [rawDataFolder trial_ID '/'];
    % filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/analysis/'];
    if ~isfolder(filePath)
        disp(['File not found ' trial_ID])
        continue
    end
    
    temp = load([filePath expID{trial} arenas{trial} ' timecourse data.mat'], 'data');
    data(trial).x = temp.data.x_loc;
    data(trial).y = temp.data.y_loc;
    data(trial).wellcenters = temp.data.wellcenters;
    data(trial).foodwell = temp.data.foodwell;

    disp([expID{trial} arenas{trial} ' ' num2str(trial) '/' num2str(ntrials)])
end
toc

% Save the structure: 

%Pull and reorganize data within the structure:
initial_vars = {'ExpGroup','baseFolder', 'T', 'data', 'figDir', 'filePath',...
                'initial_vars', 'folder', 'ntrials', 'pix2mm', 'FPS','grouped'};
clearvars('-except',initial_vars{:})
if questdlg('Save loaded data?')
    save([figDir ExpGroup ' raw.mat'],'-v7.3','data','ExpGroup', 'ntrials','T','initial_vars')
end
fprintf('Data loaded\n')
disp('next section')



%% Create measure of eccentricity 

exp = 4;
grouped(exp).name = ExpGroup;
grouped(exp).baseFolder = baseFolder;
grouped(exp).figDir = figDir;
grouped(exp).ntrials = ntrials;
grouped(exp).T = T;
grouped(exp).data = data;

data = [];

%% Save data structure: 

baseFolder = getDataPath(4,0);
pathNames = getPathNames;

figDir = createFolder([baseFolder 'plate comparisons/']);

save([figDir 'Full plate comparison raw.mat'],'-v7.3') 

initial_vars{end+1} = 'nexp';


%% Extract ecentricity for all the trials [takes a little time and requires large memory]
pix2mm = 12.8; %conversion from pixels to mm for these videos
tic
nexp = 4;
% extract the eccentricity for all the trials: 
for exp = 1:nexp
    temp = [grouped(exp).data(:).foodwell];
    grouped(exp).ntrials = sum(~isnan(temp));
    for trial = 1:grouped(exp).ntrials
        x = grouped(exp).data(trial).x; % x locations for the entire experiment
        y = grouped(exp).data(trial).y; % y locations for the entire experiment
        centre = grouped(exp).data(trial).wellcenters(:,5); %distance to center of arena 
        D = sqrt(((x-centre(1)).^2 + (y-centre(2)).^2)); %distance from center of arena
        D = D./pix2mm;
        grouped(exp).data(trial).ecent = D;
        disp([num2str(trial) '/' num2str(grouped(exp).ntrials)])
    end
end
toc


% extract all the ecentricity data into the grouped structure
for exp = 1:nexp
    temp = [20000*16*grouped(exp).ntrials,1];
    i = 1;
    for trial = 1:grouped(exp).ntrials
        y = grouped(exp).data(trial).ecent(:);
        i2 = i+length(y)-1;
        temp(i:i2) = y;
        i = i2+1;
    end
    grouped(exp).ecent = temp(~isnan(temp));
    disp(grouped(exp).name)
end


%% Figure:  trial demographics

% caviar:  plate 1 vs plate 2
kolors = {'magenta', 'dodgerblue', 'magenta', 'dodgerblue'};  
edgekolors = {'white', 'white', 'grey', 'grey'};

fig = figure; hold on
   h = bar([grouped(:).ntrials]);
   for i = 1:nexp
        h = bar(i, grouped(i).ntrials, 0.9);  % 0.9 bar width
        h.FaceColor = Color(kolors{i});
        h.EdgeColor = Color(edgekolors{i});
        h.LineWidth = 1.5;
   end
    hold off;
formatFig(fig, true);
ylabel('trial count')
set(gca, 'xcolor', 'none')

save_figure(fig, [figDir 'trial count by type'],'-png');

 
%% Figure:  histograms of ecentricity for each plate / condition
clearvars('-except',initial_vars{:})

kolors = {'magenta', 'dodgerblue', 'magenta', 'dodgerblue'};  
edgekolors = {'white', 'white', 'grey', 'grey'};

bins = 0:0.5:35;

fig = getfig('',1);
for exp = 1:nexp
    subplot(2,2,exp)
    hold on
    histogram(grouped(exp).ecent,bins,'EdgeColor','none','FaceColor',Color(kolors{exp}),'FaceAlpha',0.8)
    title(strrep(grouped(exp).name,'_',' '))
end

formatFig(fig, true,[2,2])
for exp = 1:4
    xlabel('dist to center (mm)')
    xlim([0,35])
end

save_figure(fig, [figDir 'eccentricity histogram each plate by condition'],'-png');


%% Figure: overlay of eccentricity histograms between plates
clearvars('-except',initial_vars{:})

% plate information: 
% R1 = 30; %mm for plate 1
% innerR1 = R1*sqrt(3/4); % radius of the inner 25% occupancy space R*sqrt(1/2)
% 
kolors = {'magenta', 'dodgerblue', 'magenta', 'dodgerblue'};  
edgekolors = {'white', 'white', 'grey', 'grey'};

bins = 0:0.5:35;

fig = getfig('',1);
subplot(1,2,1); hold on
for exp = 1:2
    histogram(grouped(exp).ecent,bins,'EdgeColor','none','FaceColor',Color(kolors{exp}),'FaceAlpha',0.6)
    yyaxis right
end
subplot(1,2,2); hold on
for exp = 3:4
    histogram(grouped(exp).ecent,bins,'EdgeColor','none','FaceColor',Color(kolors{exp}),'FaceAlpha',0.6)
    yyaxis right
end

formatFig(fig, true,[1,2]);
subplot(1,2,1);
    xlabel('dist to center (mm)')
    title('with food','color', 'w')
    yyaxis left
    set(gca, 'ycolor', 'none', 'box', 'off')
    yyaxis right
    set(gca, 'ycolor', 'none', 'box', 'off')
    v_line([15.6,21.6],'w',':',1) % edges of food well
    v_line(36.1,'magenta','-',1) % edges of OG arena
    v_line(33.6,'dodgerblue','-',1) % edges of new arena

subplot(1,2,2);
    xlabel('dist to center (mm)')
    title('empty','color', 'w')
    yyaxis left
    set(gca, 'ycolor', 'none', 'box', 'off')
    yyaxis right
    set(gca, 'ycolor', 'none', 'box', 'off')
    v_line([15.6,21.6],'w',':',1) % edges of food well
    v_line(36.1,'magenta','-',1) % edges of OG arena
    v_line(33.6,'dodgerblue','-',1) % edges of new arena



save_figure(fig, [figDir 'eccentricity histogram comparison for plates'],'-png');

%% Figure: how does food change the spatial distribution for a given plate?

clearvars('-except',initial_vars{:})
% 
% % plate information: 
% R1 = 30; %mm for plate 1
% innerR1 = R1*sqrt(3/4); % radius of the inner 25% occupancy space R*sqrt(1/2)

kolors = {'gold', 'gold', 'white', 'white'};  

bins = 0:0.5:35;

fig = getfig('',1);
subplot(1,2,1); hold on
for exp = 1:2:4
    histogram(grouped(exp).ecent,bins,'EdgeColor','none','FaceColor',Color(kolors{exp}),'FaceAlpha',0.6)
    yyaxis right
end
subplot(1,2,2); hold on
for exp = 2:2:4
    histogram(grouped(exp).ecent,bins,'EdgeColor','none','FaceColor',Color(kolors{exp}),'FaceAlpha',0.6)
    yyaxis right
end

formatFig(fig, true,[1,2]);
subplot(1,2,1);
    xlabel('dist to center (mm)')
    title('plate 1','color', 'w')
    yyaxis left
    set(gca, 'ycolor', 'none', 'box', 'off')
    yyaxis right
    set(gca, 'ycolor', 'none', 'box', 'off')
    v_line([15.6,21.6],'w',':',1) % edges of food well
    v_line(36.1,'magenta','-',1) % edges of OG arena
    % v_line(33.6,'dodgerblue','-',1) % edges of new arena

subplot(1,2,2);
    xlabel('dist to center (mm)')
    title('plate 2','color', 'w')
    yyaxis left
    set(gca, 'ycolor', 'none', 'box', 'off')
    yyaxis right
    set(gca, 'ycolor', 'none', 'box', 'off')
    v_line([15.6,21.6],'w',':',1) % edges of food well
    % v_line(36.1,'magenta','-',1) % edges of OG arena
    v_line(33.6,'dodgerblue','-',1) % edges of new arena


save_figure(fig, [figDir 'eccentricity histogram food vs no food'],'-png');


%% TODO: test the distributions of the trials that have the same distribution of temp protocol trials 


%% Determine the pixel to mm conversion from the entire population of trials: 

D1 = [];
for exp = [1,3] %plate 1
    for trial = 1:grouped(exp).ntrials
        well_loc = grouped(exp).data(trial).wellcenters;

        % distance between well 1 and 3
        d = [];
        well_1 = well_loc(:,1);
        well_3 = well_loc(:,3);
        d = [d;sum((well_1-well_3).^2).^0.5];
        % distance between well 2 and 4
        well_1 = well_loc(:,2);
        well_3 = well_loc(:,4);
        d = [d;sum((well_1-well_3).^2).^0.5];
        
        pixelsbetweenwells = mean(d); %pixels
        actualdistance =  36.2; %mm distance between the center of two wells criss-crossing the arena
        pix2mm = pixelsbetweenwells/actualdistance; % divisor
        D1 = [D1, pix2mm];
    end
end

D2 = [];
for exp = [2 4] %plate 1
    for trial = 1:grouped(exp).ntrials
        well_loc = grouped(exp).data(trial).wellcenters;

        % distance between well 1 and 3
        d = [];
        well_1 = well_loc(:,1);
        well_3 = well_loc(:,3);
        d = [d;sum((well_1-well_3).^2).^0.5];
        % distance between well 2 and 4
        well_1 = well_loc(:,2);
        well_3 = well_loc(:,4);
        d = [d;sum((well_1-well_3).^2).^0.5];
        
        pixelsbetweenwells = mean(d); %pixels
        actualdistance =  36.2; %mm distance between the center of two wells criss-crossing the arena
        pix2mm = pixelsbetweenwells/actualdistance; % divisor
        D2 = [D2, pix2mm];
    end
end

% pull the date strings for the data
numeric_dates1 = [datenum(grouped(1).T.Date, 'mm.dd.yyyy'); datenum(grouped(3).T.Date, 'mm.dd.yyyy')];
numeric_dates2 = [datenum(grouped(2).T.Date, 'mm.dd.yyyy'); datenum(grouped(4).T.Date, 'mm.dd.yyyy')];

cut_off_one = datenum('11.10.2023', 'mm.dd.yyyy');

% Figures: plot the pixel to mm ratio over time and see if there is a time
% trend or it has changed significantly at any point in time...

fig = getfig('',1);
    hold on
    scatter(numeric_dates1,D1,30, Color('magenta'),'filled')
    scatter(numeric_dates2,D2,30, Color('dodgerblue'),'filled')
    
    datetick('x', 'mm.dd.yyyy', 'keepticks');
    xlabel('Date');
    ylabel('pixel to mm conversion');
    formatFig(fig,true);
    v_line(cut_off_one,'w')
    save_figure(fig, [figDir 'pix2mm conversion distribution timeline'],'-png');

% filter and find the average value for each time period: 
% before 11.10.23 (pre move to 100 college st)

idx = numeric_dates1<=cut_off_one;

conversion(1).name = 'Cedar st plate 1';
conversion(1).cutoff = cut_off_one;
conversion(1).values = D1(idx);
conversion(1).pix2mm = median(conversion(1).values);
conversion(1).color = 'pink';
disp([conversion(1).name ' pixel to mm converion: ' num2str(conversion(1).pix2mm)])

conversion(2).name = 'College st plate 1';
conversion(2).cutoff = cut_off_one;
conversion(2).values = D1(~idx);
conversion(2).pix2mm = median(conversion(2).values);
conversion(2).color = 'magenta';
disp([conversion(2).name ' pixel to mm converion: ' num2str(conversion(2).pix2mm)])

conversion(3).name = 'Cedar st plate 2';
conversion(3).cutoff = cut_off_one;
conversion(3).values = D2;
conversion(3).pix2mm = median(conversion(3).values);
conversion(3).color = 'dodgerblue';
disp([conversion(3).name ' pixel to mm converion: ' num2str(conversion(3).pix2mm)])

fig = getfig('',1,[940 900]); hold on 
    for i = 1:3
        histogram(conversion(i).values,'FaceColor',Color(conversion(i).color))
        v_line(conversion(i).pix2mm,conversion(i).color,'-',2)
    end
    xlabel('pixel to mm conversion number')
    formatFig(fig, true);
    set(gca,'ycolor', 'none')
    save_figure(fig, [figDir 'pix2mm conversion distribution histogram'],'-png');

% Create a pixel to mm conversion that can be used for the forseable future
% with low resolution data trials
save([figDir 'pixel conversion.mat'],'conversion');


%% Visual plate comparison of the new test numbers: 

p1_path = 'S:\Evyn\DATA\Raw Data\05.07.2025\C2_hold_35_caviar_1.avi';
p2_path = 'S:\Evyn\DATA\Raw Data\05.07.2025\C1_hold_35_caviar_1.avi';

% PLATE 1: 
movieInfo = VideoReader(p1_path); %read in video
demoImg = (read(movieInfo,1));

pix2mm = conversion(2).pix2mm;
fig = figure;
    imshow(demoImg); hold on
    dia = nan(1,4);
    cen = nan(2,4);
    for i = 1:4
        roi = drawcircle; % manually add in the circle over the food well
        newR = roi.Radius/pix2mm;
        disp(['Plate 1: ' num2str(newR)])
        scatter(roi.Center(1),roi.Center(2),15,'y', 'filled')
        cen(:,i) = roi.Center;
        dia(i) = newR;
    end

    disp(['arena radius = ' num2str(mean(dia))])
    % test the playback of the arena size fit: 
    fig = figure;
    imshow(demoImg); hold on
    for i = 1:4
        R = mean(dia);
        viscircles(cen(:,i)',R*pix2mm,'color', 'w','LineWidth',0.25);
    end
plate1_R = 35.2; % full reach (from manual determined value) 

% PLATE 2: 
movieInfo = VideoReader(p2_path); %read in video
demoImg = (read(movieInfo,1));

pix2mm = conversion(3).pix2mm;
fig = figure;
    imshow(demoImg); hold on
    dia = nan(1,4);
    cen = nan(2,4);
    for i = 1:4
        roi = drawcircle; % manually add in the circle over the food well
        newR = roi.Radius/pix2mm;
        disp(['Plate 2: ' num2str(newR)])
        scatter(roi.Center(1),roi.Center(2),15,'y', 'filled')
        cen(:,i) = roi.Center;
        dia(i) = newR;
    end

    disp(['arena radius = ' num2str(mean(dia))])
    % test the playback of the arena size fit: 
    fig = figure;
    imshow(demoImg); hold on
    for i = 1:4
        R = mean(dia);
        viscircles(cen(:,i)',R*pix2mm,'color', 'w','LineWidth',0.25);
    end

plate2_R = 33.7; % full reach (from manual determined value) 
    



% % distance between well 1 and 3
% d = [];
% well_1 = well_loc(:,1);
% well_3 = well_loc(:,3);
% d = [d;sum((well_1-well_3).^2).^0.5];
% % distance between well 2 and 4
% well_1 = well_loc(:,2);
% well_3 = well_loc(:,4);
% d = [d;sum((well_1-well_3).^2).^0.5];

% pixelsbetweenwells = mean(d); %pixels
% actualdistance =  36.2; %mm distance between the center of two wells criss-crossing the arena
% pix2mm = actualdistance/pixelsbetweenwells; % multiplier

% for PLATE 1: 
% well_loc = readPoints(demoImg,4); % click on center positions of the wells in the arena
% WC = well_loc;
% x1 = WC(1,1:2:4);
% y1 = WC(2,1:2:4);
% x2 = WC(1,2:2:4);
% y2 = WC(2,2:2:4);
% [xi,yi] = polyxpoly(x1,y1,x2,y2);
% centre = [xi,yi];

% R = 35.9; % radius of plate 1:
pix2mm = conversion(3).pix2mm;
fig = figure;
    imshow(demoImg); hold on
    scatter(centre(1),centre(2),15,'y', 'filled')
    viscircles(centre,R*pix2mm,'color', 'w','LineWidth',0.25);
    % find actual plate size: 
    % roi = drawcircle; % manually add in the circle over the food well
    % newR = roi.Radius/pix2mm;
    newR = 33.99;
    viscircles(centre,newR/pix2mm,'color', 'w','LineWidth',0.25);

%% Determine the max accessible space for each of the plates: 

% since the data is already generated with potentially incorrect conversion
% ratios, we need to recalculate the distance from the center, which we can
% then use to extract the appropriate extreme reachable distances

% 1) calculate the eccentricity of the data trials: 

numeric_dates1 = [datenum(grouped(1).T.Date, 'mm.dd.yyyy'); datenum(grouped(3).T.Date, 'mm.dd.yyyy')];
numeric_dates2 = [datenum(grouped(2).T.Date, 'mm.dd.yyyy'); datenum(grouped(4).T.Date, 'mm.dd.yyyy')];

% organize data for forward compatability
con = getConversion;
cut_off_one = con(1).cutoff;

temp = nan([20000*16*300,1]);

[eccentricity(1).all,eccentricity(2).all,eccentricity(3).all] = deal(temp);
[eccentricity(1).N,eccentricity(2).N,eccentricity(3).N] = deal(0);
[eccentricity(1).trial,eccentricity(2).trial,eccentricity(3).trial] = deal([]);

% ECCENTRICITY
for exp = 1:nexp
    tic
    idx = [datenum(grouped(exp).T.Date, 'mm.dd.yyyy')]<=cut_off_one;
    for trial = 1:grouped(exp).ntrials
        % get the appropriate pixel to mm conversion for these trials: 
        switch exp
            case {1,3} % plate 1 food | no food
                if idx(trial) % old lab space (cedar st) 
                    pix2mm = con(1).pix2mm;
                    conType = 1;
                else % new lab space (college st) 
                    pix2mm = con(2).pix2mm;
                    conType = 2;
                end
            case {2,4}  % plate 2 food | no food
                pix2mm = con(3).pix2mm;
                conType = 3;
        end

        eccentricity(conType).trial = [eccentricity(conType).trial; exp, trial];

        xi = grouped(exp).data(trial).wellcenters(1,5);
        yi = grouped(exp).data(trial).wellcenters(2,5);
        D = sqrt((grouped(exp).data(trial).x-xi).^2 + (grouped(exp).data(trial).y-yi).^2)./pix2mm; %distance from center of arena

        i_start = eccentricity(conType).N+1;
        D = D(:);
        D(isnan(D)) = [];
        i_end = i_start+length(D)-1;
        eccentricity(conType).all(i_start:i_end) = D; 
        eccentricity(conType).N = i_end;

    end
    toc
end

for i = 1:3
    eccentricity(i).all(isnan(eccentricity(i).all)) = [];
end

% Now look for max eccentricities across the different arenas for each time
% period / functional setup 

% key point = the conType 1 and 2 should be the same after the appropriate
% pix2mm conversion pulled in

% Figure:  histograms of ecentricity for each plate / condition

% kolors = {'magenta', 'dodgerblue'};  
kolors = {'pink', 'magenta', 'dodgerblue'};  

bins = 0:0.5:35;

fig = getfig('',1); hold on
for i = 1:3
    histogram(eccentricity(i).all,bins,'EdgeColor','none','FaceColor',Color(kolors{i}),'FaceAlpha',0.8)
end
for i = 1:3
    v_line(con(i).circle75,kolors{i}, '-', 2)
end
formatFig(fig, true)
% xlim([25,35])
xlabel('dist to center (mm)')
ylabel('count')

set(gca, 'YScale', 'log')


save_figure(fig, [figDir 'eccentricity with outer ring cutoffs'],'-png');


% FIGURE: show each of the groups with their outer and inner 25% lines
fig = getfig('',1);
for i = 1:3
    subplot(1,3,i); hold on
    histogram(eccentricity(i).all,bins,'EdgeColor','none','FaceColor',Color(kolors{i}),'FaceAlpha',0.8)
    v_line(con(i).circle75,'w', '-', 2)
    v_line(con(i).R,'w', '-', 2)
    xlabel('dist to center (mm)')
end
formatFig(fig, true, [1,3]);
subplot(1,3,1)
ylabel('count')
save_figure(fig, [figDir 'eccentricity with outer ring cutoffs separated'],'-png');

%%  NEW: determine where the outlier points from the outer ring occupancy come from in the new plate food trials
% con = getConversion;
% 
% % find the eccentricities for all of these trials
% logNames = table;
% wellD = [];
% % allD = [];
% exp = 3; % new plate with food (since this has the outlier distances)
% for trial = 1:grouped(exp).ntrials
%     pix2mm = con(3).pix2mm;
%     wellcenters = grouped(exp).data(trial).wellcenters;
%     xi = wellcenters(1,5);
%     yi = wellcenters(2,5);
%     D = sqrt((grouped(exp).data(trial).x-xi).^2 + (grouped(exp).data(trial).y-yi).^2)./pix2mm; %distance from center of arena
%     % allD = [allD; D(:)];
%     % log trial if the D is greater than 30mm (or the achievable outer ring
%     % distance, obstensibly 
%     if any(any(D>con(3).R))
%         logNames(end+1,:) = grouped(exp).T(trial,1:4); % add file info to the log
% 
%         % determine distance between the wells for this trial
%         d1 = sqrt((wellcenters(1,1) - wellcenters(1,3)).^2 + (wellcenters(2,1) - wellcenters(2,3)).^2)./pix2mm; %distance from center of arena
%         d2 = sqrt((wellcenters(1,2) - wellcenters(1,4)).^2 + (wellcenters(2,2) - wellcenters(2,4)).^2)./pix2mm; %distance from center of arena
%         wellD = [wellD; d1, d2];
%     end
%     disp(trial)
% end
% 
% % histogram of the eccentricities to see why it doesn't match the earlier
% % error...
% bins = 0:0.5:35;
% figure; 
% histogram(allD,bins)
% tempLog = find(allD>30);
% 
% 
% 
% 
% % 
% % [eccentricity(1).all,eccentricity(2).all,eccentricity(3).all] = deal(temp);
% % [eccentricity(1).N,eccentricity(2).N,eccentricity(3).N] = deal(0);
% % [eccentricity(1).trial,eccentricity(2).trial,eccentricity(3).trial] = deal([]);
% % 
% 
% 






%% now if we go back and then compare these lines and distributions for the
% food vs no food trials, what does that look like?

temp = nan([20000*16*300,1]);
for i = 1:6
    eccentricity(i).all = deal(temp);
    eccentricity(i).N = deal(0);
    eccentricity(i).trial = deal([]);
end

% conTypes: 
% 1 = plate 1, cedar st, food
% 2 = plate 1, cedar st, no food
% 3 = plate 1, college st, food
% 4 = plate 1, college st, no food
% 5 = plate 2, college st, food
% 6 = plate 2, college st, no food

cut_off_one = con(1).cutoff;
% ECCENTRICITY
for exp = 1:nexp
    idx = [datenum(grouped(exp).T.Date, 'mm.dd.yyyy')]<=cut_off_one;
    for trial = 1:grouped(exp).ntrials   
        % get the appropriate pixel to mm conversion for these trials: 
        switch exp
            case 1 % 1 old plate, food
                if idx(trial) % cedar st 
                    pix2mm = con(1).pix2mm;
                    conType = 1;
                else % college st
                    pix2mm = con(2).pix2mm;
                    conType = 3;
                end
            case 3 % old plate, no food
                if idx(trial) % cedar st
                    pix2mm = con(1).pix2mm;
                    conType = 2;
                else % college st
                    pix2mm = con(2).pix2mm;
                    conType = 4;
                end
            case 2 % new plate, food
                pix2mm = con(3).pix2mm;
                conType = 5;
            case 4  % new plate, no food
                pix2mm = con(3).pix2mm;
                conType = 6;
        end
        eccentricity(conType).trial = [eccentricity(conType).trial; exp, trial];
        xi = grouped(exp).data(trial).wellcenters(1,5);
        yi = grouped(exp).data(trial).wellcenters(2,5);
        D = sqrt((grouped(exp).data(trial).x-xi).^2 + (grouped(exp).data(trial).y-yi).^2)./pix2mm; %distance from center of arena

        i_start = eccentricity(conType).N+1;
        D = D(:);
        D(isnan(D)) = [];
        i_end = i_start+length(D)-1;
        eccentricity(conType).all(i_start:i_end) = D; 
        eccentricity(conType).N = i_end;
    end
end

% FIGURE: plot out the data to see how it matches up:

bins = 0:0.5:35;
kolors = {'pink','pink','magenta', 'magenta', 'dodgerblue', 'dodgerblue'};  
cc = [1 1 2 2 3 3];

fig = getfig('',1);
for i = 1:6
    subplot(3,2,i)
    histogram(eccentricity(i).all, bins,'EdgeColor','none','FaceColor',Color(kolors{i}),'FaceAlpha',0.8)
    v_line(con(cc(i)).circle75,'w', '-', 2)
    v_line(con(cc(i)).R,'w', '-', 2)
end
formatFig(fig, true, [3,2]);
for i = 1:2:6
    subplot(3,2,i)
    ylabel('count')
end
for i = 5:6
    subplot(3,2,i)
    xlabel('mm to center')
end

% figDir = '/Users/evyndickinson/Documents/Jeanne Lab/DATA/Grouped Data Structures/plate comparisons/';

save_figure(fig, [figDir 'eccentricity with outer ring cutoffs separated by type'],'-png');

%% QUICK FIG FOR COMP WITH HIGH RES EXPERIMENTS
% figure of plate two high res experiments distribution to make sure it all
% aligns 

fig = getfig('',1);
i = 5; % new plate, with food 
histogram(eccentricity(i).all, bins,'EdgeColor','k','FaceColor',Color('dodgerblue'),'FaceAlpha',0.8)
v_line(con(cc(i)).circle75,'w', '-', 2)
v_line(con(cc(i)).R,'w', '-', 2)
formatFig(fig, true);
xlabel('distance to center (mm)')
ylabel('count')

save_figure(fig, [figDir, 'new plate with food eccentricity histogram'],'-png',0,0,'-r80');


%% TEST HIGH RESOLUTION PLATE VIDEO CONVERSION NUMBERS
% need to load data from a grouped structure in the high resolutoin courtship video
% structures
D = [];
for i = 1:num.trials
        well_loc = fly(i).well.center';

        % distance between well 1 and 3
        d = [];
        well_1 = well_loc(:,1);
        well_3 = well_loc(:,3);
        d = [d;sum((well_1-well_3).^2).^0.5];
        % distance between well 2 and 4
        well_1 = well_loc(:,2);
        well_3 = well_loc(:,4);
        d = [d;sum((well_1-well_3).^2).^0.5];
        
        pixelsbetweenwells = mean(d); %pixels
        actualdistance =  36.2; %mm distance between the center of two wells criss-crossing the arena
        pix2mm = pixelsbetweenwells/actualdistance; % divisor
        D = [D, pix2mm];
end

fig = getfig('',1);
histogram(D,'FaceColor',Color('dodgerblue'),'FaceAlpha',0.8)
v_line(median(D),'white','-',2)
formatFig(fig,true);
ylabel('count')
xlabel('pixel to mm conversion')
title(['High resolution plate 2 : ' num2str(median(D))], 'color', 'w')
save_figure(fig, [figDir 'pix2mm conversion distribution histogram'],'-png');

pix2mm = 30.0093;

% test to see how the new alignment data works with the described plate measurements
% from the above analyses and getConversion

% use the head position of the flies to plot out the histogram of their eccentricies
% and then overlay the cut-offs from the previous data % TODO :  working here  

% % temporary save: 
% load('K:\DATA\Courtship Videos\grouped\Berlin F LRR 25-17 caviar MF\working data.mat');

% determine eccentricity of the flies in the arenas with the new pix2mm conversion

eccentricity = [];
for trial = 1:num.trials
    x = [fly(trial).data(1).rawX(:,body.head); fly(trial).data(2).rawX(:,body.head)];
    y = [fly(trial).data(1).rawY(:,body.head); fly(trial).data(2).rawY(:,body.head)];

    xi = fly(trial).well.center(5,1);
    yi = fly(trial).well.center(5,2);
    D = sqrt((x-xi).^2 + (y-yi).^2)./pix2mm; %distance from center of arena
    D(isnan(D)) = [];
    eccentricity = [eccentricity; D];
    disp(trial)
end

bins = 0:0.5:35;

fig = getfig('',1);
    histogram(eccentricity,bins,'FaceColor', Color('dodgerblue'), 'FaceAlpha', 0.8)
    v_line(conversion(3).circle75,'w', '-', 2)
    v_line(conversion(3).R,'w', '-', 2)
    formatFig(fig, true);
    ylabel('count')
    xlabel('mm to center')
    save_figure(fig, [figDir, 'eccentricity all time points and bins'],fig_type);


%% Visual arena regions

% Low resolution trial example: 
vidpath = "S:\Evyn\DATA\Raw Data\05.01.2025\C2_hold_33_caviar_1.avi";
movieInfo = VideoReader(vidpath); %read in video
demoImg = (read(movieInfo,1));
well_loc = readPoints(demoImg,4); % click on center positions of the wells in the arena
conversion = getConversion;
pix2mm = conversion(2).pix2mm; % multiplier
R = conversion(2).R;
x1 = well_loc(1,1:2:4);
y1 = well_loc(2,1:2:4);
x2 = well_loc(1,2:2:4);
y2 = well_loc(2,2:2:4);
[xi,yi] = polyxpoly(x1,y1,x2,y2);
centre = [xi,yi];
disp('select center of the food well')
food_loc = readPoints(demoImg,1); % click on center positions of the wells in the arena

% crop image into the appropriate area around it
crop_dist = 37 * pix2mm; %mm
fig = getfig('',1,[900,900]);
    % a = imadjust(demoImg,[53/255,149/255]); % plate 1
    a = imadjust(demoImg,[53/255,175/255]); % plate 2
    imshow(a); hold on
    xlim([centre(1)-crop_dist, centre(1)+crop_dist])
    ylim([centre(2)-crop_dist, centre(2)+crop_dist])
    viscircles(centre,conversion(2).R*pix2mm,'color', 'w','LineWidth',0.25);
    viscircles(centre,conversion(2).circle75*pix2mm,'color', 'w','LineWidth',0.25);
    viscircles(food_loc',conversion(2).circle10*pix2mm,'color', Color('gold'),'LineWidth',0.25);
    viscircles(food_loc',conversion(2).circle7*pix2mm,'color', Color('red'),'LineWidth',0.25);
    viscircles(food_loc',conversion(2).circle5*pix2mm,'color', Color('purple'),'LineWidth',0.25);





fig = figure;
    imshow(demoImg); hold on
    viscircles(centre,R*pix2mm,'color', 'w','LineWidth',0.25);
    % find actual plate size: 
    % roi = drawcircle; % manually add in the circle over the food well
    % newR = roi.Radius/pix2mm;
    newR = 33.6;
    viscircles(centre,newR*pix2mm,'color', 'w','LineWidth',0.25);


vidpath = "S:\Evyn\DATA\Courtship Videos\02.05.2025\Berlin_courtship_F_LRR_caviar_ramp1\compiled_video_1.avi";
movieInfo = VideoReader(vidpath); %read in video
demoImg = (read(movieInfo,1));

fig = getfig('',1,[900,900]);
    imshow(demoImg); hold on
    viscircles(well.center(5,:),conversion(4).R*pix2mm,'color', 'w','LineWidth',0.25);
    viscircles(well.center(5,:),conversion(4).circle75*pix2mm,'color', 'w','LineWidth',0.25);
    viscircles(well.food,conversion(4).circle10*pix2mm,'color', Color('gold'),'LineWidth',0.25);
    viscircles(well.food,conversion(4).circle7*pix2mm,'color', Color('red'),'LineWidth',0.25);
    viscircles(well.food,conversion(4).circle5*pix2mm,'color', Color('purple'),'LineWidth',0.25);


% vidpath = "S:\Evyn\DATA\Raw Data\05.01.2025\C1_hold_33_caviar_1.avi";
































