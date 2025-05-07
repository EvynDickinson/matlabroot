
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









