
clear

%% Load in multiple trials that are grouped into a structure
% baseFolder = getCloudPath;
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
baseFolder = getCloudPath;

% Select structure to load:
[~,~,structInfo] = getExcelStructureNames(true);
ExpGroup = structInfo.StructName;
ntrials = structInfo.numTrials;

% Make a data folder for structure images
figDir = [baseFolder 'Data structures\' ExpGroup '\'];
if ~exist(figDir, 'dir'); mkdir(figDir); end

% Load data from each trial in the structure
data = [];
fprintf('\nLoading trials: \n')
for trial = 1:ntrials
    % print the experiments as they are loaded
    trialExpID = excelfile{structInfo.rowNum(trial), Excel.expID};
    trialDate = excelfile{structInfo.rowNum(trial), Excel.date};
    trialArena = excelfile{structInfo.rowNum(trial), Excel.arena};
    trialName = [trialExpID ' Arena ' trialArena ' ' trialDate];
    disp(trialName)

    % build the path for the trial data
    dirc = [baseFolder, trialDate, '\Arena ' trialArena '\analysis\' trialExpID trialArena  ' timecourse data.mat'];
    
    % load data
    todel = load(dirc);
    variList = fieldnames(todel);
%     todel = load(dirc, varList{:});
    data(trial).trialName = trialName;
    for ii = 1:length(variList)
        data(trial).(variList{ii}) = todel.(variList{ii});
    end
end; clear varList
initial_vars = {'baseFolder', 'data', 'ExpGroup', 'ntrials', 'initial_vars', 'figDir'};
clearvars('-except',initial_vars{:})
fprintf('Data loaded\n')

%% visual check of temperature alignment across the two experiments:
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.time;
    Y = data(trial).occupancy.temp;
    plot(X, Y, 'linewidth', 1)
end
xlabel('Time (min)')
ylabel('Temp (\circ)')
title({'temperature alignment across experiments';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true)

save_figure(fig, [figDir ExpGroup ' temperature alignment'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Distance to food vs temperature (no time component --  quick and dirty)
sSpan = 240;
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.temp;
    for well = 1:4
        Y = (data(trial).dist2wells(well).N(:,1));
        % reorder the data by temperature:
        [plotX, idx] = sort(X);
        plot(plotX, smooth(Y(idx),sSpan), 'linewidth', 1,'color', pullFoodColor(data(trial).wellLabels{well}))
    end
end
xlabel('Time (min)')
ylabel('Temp (\circ)')
title({'Location from food sources by temperature';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true)

save_figure(fig, [figDir ExpGroup ' all temp vs distance'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Average across trials:

strfind(data(1).wellLabels{4},'Plant')
sSpan = 240;
[plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        % plant food well:
        if strfind(data(trial).wellLabels{well},'Plant')
            x = data(trial).occupancy.temp';
            y = data(trial).dist2wells(well).N(:,1);
            plant = [plant; x, y];
        end
        % yeast food well:
        if strfind(data(trial).wellLabels{well},'Yeast')
            x = data(trial).occupancy.temp';
            y = data(trial).dist2wells(well).N(:,1);
            yeast = [yeast; x, y];
        end
        % empty food well:
        if strfind(data(trial).wellLabels{well},'Empty')
            x = data(trial).occupancy.temp';
            y = data(trial).dist2wells(well).N(:,1);
            empty = [empty; x, y];
        end 
    end
end 

threshHigh = 19.88;
threshLow = 8.02;
LW = 2;
fig = figure; hold on
% reorder the PLANT data by temperature:
pData = sort(plant,1);
pData(pData(:,1)<threshLow,:) = [];
pData(pData(:,1)>threshHigh,:) = [];
plot(pData(:,1), smooth(pData(:,2),sSpan), 'linewidth', LW,'color', pullFoodColor('Plant'))
% reorder the YEAST data by temperature:
pData = sort(yeast,1);
pData(pData(:,1)<threshLow,:) = [];
pData(pData(:,1)>threshHigh,:) = [];
plot(pData(:,1), smooth(pData(:,2),sSpan), 'linewidth', LW,'color', pullFoodColor('Yeast'))
% reorder the EMPTY data by temperature:
pData = sort(empty,1);
pData(pData(:,1)<threshLow,:) = [];
pData(pData(:,1)>threshHigh,:) = [];
plot(pData(:,1), smooth(pData(:,2),sSpan), 'linewidth', LW,'color', pullFoodColor('Empty'))
% Labels:
xlabel('temperature (\circC)')
ylabel('distance to well (a.u.)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true)

save_figure(fig, [figDir ExpGroup ' temp vs distance'], '-png');
clearvars('-except',initial_vars{:})

%% WILL NEED A SECTION ON WELL IDENTITY MATCHING FOR ONCE WELLS ARE SWAPPED
% AROUND

%% Determine a well-identity matching list:
well_contents = unique(data(1).wellLabels);

for trial = 1:ntrials
    for well = 1:length(well_contents)
        wellID(well).loc(:,trial) = find(strcmpi(data(trial).wellLabels,well_contents{well}));
        wellID(well).name = well_contents{well};
    end
end

initial_vars = [initial_vars {'well_contents', 'wellID'}];
clearvars('-except',initial_vars{:})

%% Group occupancy across the trials

% for now: average across the experiments (assuming they are TEMP locked)
[A,B,C,D] = deal([]);
for trial = 1:ntrials
    A = [A,data(trial).occupancy.occ(:,1)]; %(:,well) % well occupancy probability for each well
    B = [B, data(trial).occupancy.occ(:,2)];
    C = [C, data(trial).occupancy.occ(:,3)];
    D = [D, data(trial).occupancy.occ(:,4)];
end
% average across the trials TODO update the timing according to temp alignment
for trial = 1:ntrials
    temperature(trial).log = data(trial).tempLog(:,2);
end
% TEMP ALIGNED TIME (ONLY WORKS FOR TEMP HOLD OR SAME GROUP TRIALS
TAT = data(1).time();
temperature = data(1).occupancy.temp;
wellLabels = data(1).wellLabels;

nrows = 4;
ncols = 1;
subs(1).idx = 1;
subs(2).idx = 2:4;

fig = getfig;
subplot(nrows, ncols, subs(1).idx)
    plot(TAT,temperature, 'linewidth', 2, 'color', 'w')
    ylim([5,30])
    ylabel('Temp (\circC)')
subplot(nrows, ncols, subs(2).idx)
    sSpan = 240;
    hold on
    Y = [smooth(mean(A,2),sSpan), smooth(mean(B,2),sSpan),...
         smooth(mean(C,2),sSpan), smooth(mean(D,2),sSpan)];
    h = area(TAT,Y);
    for well = 1:4
        h(well).FaceColor = pullFoodColor(data(1).wellLabels{well});
    end
    set(gca, 'tickdir', 'out')
    l = legend(strrep(wellLabels,'_','-'));
    set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k',...
    'position', [0.7457 0.6520 0.0883 0.0772]);% [0.8780 0.8119 0.0963 0.1126])%

ylabel('Occupancy')
xlabel('Time (min)')
formatFig(fig, true, [nrows, ncols], subs)
subplot(nrows, ncols, subs(1).idx)
set(gca, 'XColor', 'k')
title(strrep(ExpGroup,'_','-'),'color', 'w')

save_figure(fig, [figDir 'Quadrant Occupation'], '-png');

clearvars('-except',initial_vars{:})

%% Temperature and occupancy relationship: [takes a few seconds to run]
wellLabels = data(1).wellLabels;
% organize the occupation by fly count in each temperature
% --> temperature for each fly at each timepoint
tic
[sumA,sumB,sumC,sumD,tempAll] = deal(zeros(30,1));
for trial = 1:ntrials
    [A,B,C,D] = deal([]);
    temperature = round(data(trial).occupancy.temp);
    raw = data(trial).occupancy.count;
    for ii = 1:size(raw,1) %for each data point
        A = [A; temperature(ii)*ones(raw(ii,1),1)];
        B = [B; temperature(ii)*ones(raw(ii,2),1)];
        C = [C; temperature(ii)*ones(raw(ii,3),1)];
        D = [D; temperature(ii)*ones(raw(ii,4),1)];
    end
    for temp = 5:30
        sumA(temp) = sumA(temp) + sum(A==temp);
        sumB(temp) = sumB(temp) + sum(B==temp);
        sumC(temp) = sumC(temp) + sum(C==temp);
        sumD(temp) = sumD(temp) + sum(D==temp);
        tempAll(temp) = tempAll(temp) + sum(temperature==temp);
    end
end
toc

% normalize the number of flies for each temperature but the total time
% spent at a given temperature:
tempFrac = tempAll./(sum(tempAll));
ratio = 1; %each point will represent 20 flies normalized for total time spent at each temperature
arenaA = round(sumA.*tempFrac./ratio);
arenaB = round(sumB.*tempFrac./ratio);
arenaC = round(sumC.*tempFrac./ratio);
arenaD = round(sumD.*tempFrac./ratio);


% increase points for each temperature

[yA,yB,yC,yD] = deal([]);
for temp = 1:length(arenaA)
    yA = [yA; temp*ones(arenaA(temp),1)];
    yB = [yB; temp*ones(arenaB(temp),1)];
    yC = [yC; temp*ones(arenaC(temp),1)];
    yD = [yD; temp*ones(arenaD(temp),1)];
end
xA = 1*ones(length(yA),1);
xB = 2*ones(length(yB),1);
xC = 3*ones(length(yC),1);
xD = 4*ones(length(yD),1);

% Plot the occupation by temperature plots
sz = 30;
fig = getfig;
    hold on
    swarmchart(xA,yA,sz, pullFoodColor(wellLabels{1}), 'filled');
    swarmchart(xB,yB,sz, pullFoodColor(wellLabels{2}), 'filled');
    swarmchart(xC,yC,sz, pullFoodColor(wellLabels{3}), 'filled');
    swarmchart(xD,yD,sz, pullFoodColor(wellLabels{4}), 'filled');

% swarmchart(x,y,sz,c)














%% Collapse Well ROIs across all trials for each temperature to create 'occupancyFrames'
gridSize = [50,50];
nvids = length(data(1).temp);

% loop through all trials
for trial = 1:ntrials
  temp = data(trial).temp; 
  for ii = 1:4 % well#
    % crop each of the frames from the video into well rois
    frames = data(trial).frame;
    mask = data(trial).quadMask(ii).mask;
    props = regionprops(~mask, 'BoundingBox');
    mask = repmat(mask,[1,1,size(frames,3)]); %mask for all frames
    frames(mask) = 0; %mask out other data outside ROI
    for vid = 1:nvids
        img = imcrop(frames(:,:,vid), props.BoundingBox); %crop to the circle
        maskedImage(:,:,vid) = img;
        roiImage(:,:,vid) = imresize(img, gridSize); %bin across pixels
    end
    well(ii).img = roiImage;
    well(ii).maskedImg = maskedImage;
    well(ii).temp = temp;
  end
  data(trial).well = well;
end

% align across videos roughly (TODO: refine this later)
Wells = [];
for roi = 1:4
    Wells(roi).name = data(1).trialData(1).params.(['well_' num2str(roi)]);
    Wells(roi).img = zeros(size(data(1).well(roi).img));
    for trial = 1:ntrials
        Wells(roi).img = Wells(roi).img + data(trial).well(roi).img;
        Wells(roi).occ(:,trial) = data(trial).y(:,roi); % determine this...
    end
    % find density range for the rois:
    M(roi) = max(Wells(roi).img,[],'all');
end

% Normalize the well density figure:
gridS = 10;
n = gridS+2;
well = [];
for roi = 1:4
    well(roi).img = Wells(roi).img / max(M);
    dummy = nan(n,n,nvids);
    img = imresize(well(roi).img, [gridS,gridS]);
    dummy(2:gridS+1,2:gridS+1,:) = img;
    well(roi).bin = dummy;
end
% grouped image:
occupancyFrames = [well(1).bin, well(2).bin;...
                   well(4).bin, well(3).bin];
               
occupancyNames = {Wells(1).name, Wells(2).name;...
                   Wells(4).name, Wells(3).name};              

% quick vid preview of any structure over time
fig = figure; set(fig, 'color', 'k');
disp(occupancyNames)
for vid = 1:nvids
    img = occupancyFrames(:,:,vid);
    imagesc(img)
    pause(0.1)
end
     






%% IDEAS for later:

% swarmchart(x,y,sz,c)



























