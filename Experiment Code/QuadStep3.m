

%% Load in multiple trials that are grouped into a structure
% baseFolder = getCloudPath;
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
baseFolder = getCloudPath;

% Select structure to load:
[~,~,structInfo] = getExcelStructureNames(true);
ExpGroup = structInfo.StructName;
ntrials = structInfo.numTrials;

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
    dirc = [baseFolder, trialDate, '\analysis\' trialExpID trialArena ' timecourse data.mat'];
    
    % load data
    todel = load(dirc);
    variList = fieldnames(todel);
%     todel = load(dirc, varList{:});
    data(trial).trialName = trialName;
    for ii = 1:length(variList)
        data(trial).(variList{ii}) = todel.(variList{ii});
    end
end; clear varList
initial_vars = {'baseFolder', 'data', 'ExpGroup', 'ntrials', 'initial_vars'};
clearvars('-except',initial_vars{:})
fprintf('Data loaded\n')

% WILL NEED A SECTION ON TEMPERATURE ALIGNMENT ACROSS EXPERIMENTS

% WILL NEED A SECTION ON WELL IDENTITY MATCHING FOR ONCE WELLS ARE SWAPPED
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
ylabel('Occupancy')
xlabel('Time (min)')
formatFig(fig, true, [nrows, ncols], subs)
subplot(nrows, ncols, subs(1).idx)
set(gca, 'XColor', 'k')

% titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
% title(titleName,'color', 'w')


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



























