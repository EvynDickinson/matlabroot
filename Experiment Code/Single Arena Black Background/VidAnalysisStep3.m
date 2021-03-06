% 
% % optional update cross numbers in excel sheet
% Export_Cross_Ns;

clear; close all

%% Load in multiple flies that are grouped together:
% baseFolder = getCloudPath;
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
baseFolder = getCloudPath;
varList = {'folder', 'quadMask','trialData', 'tframe','frame', 'temp','y','nflies'}; % variables to load from each trial

% Select structure to load:
[~,~,structInfo] = getExcelStructureNames(true);
ExpGroup = structInfo.StructName;
ntrials = structInfo.numTrials;

% Load data from each trial
fprintf('\nLoading trials: \n')
for trial = 1:ntrials
    % print the experiments as they are loaded
    trialName = excelfile{structInfo.rowNum(trial), Excel.expID};
    trialDate = excelfile{structInfo.rowNum(trial), Excel.date};
    disp([trialDate,'  ' trialName])

    % build the path for the trial data
    dirc = [baseFolder, trialDate, '\analysis\' trialName '_analysis'];
    
    % load data
    todel = load(dirc, varList{:});
    for ii = 1:length(varList)
        data(trial).(varList{ii}) = todel.(varList{ii});
    end
end; clear temp
fprintf('Data loaded\n')

%% Plot ROI occupancy across all trials




% % extract occupancy probability per region
% for vid = 1:nvids
%     frame(:,:,vid) = trialData(vid).videoData.occ_Prob/...
%                      trialData(vid).params.num.framesPerVid; % occupancy prob for this frame
%     temp(vid) = mean(trialData(vid).videoData.tempLog(:,2)); %avg temp during video
%     
%     for roi = 1:4
%         img = frame(:,:,vid); %blank image
%         img(quadMask(roi).mask) = 0;
%         y(vid,roi) = sum(sum(img));
%     end
% end
% 
% % PLOT THE CRAP OUTTA THAT OCCUPANCY!
% LW = 2;
% fig = getfig('');
% hold on
% for roi = 1:4
%     lngd{roi} = trialData(1).params.(['well_' num2str(roi)]);
%     switch lngd{roi}
%         case 'OTC'
%             kolor = Color('orange');
%         case 'Plant'
%             kolor = Color('aqua');
%         case 'Empty'
%             kolor = Color('grey');
%         case 'Yeast'
%             kolor = Color('fuchsia');
%     end
%     plot(flip(temp), flip(y(:,roi)), 'color', kolor, 'linewidth', LW)
% end
% xlim([8,22])
% 
% % labels
% title(strrep([expName ' occupation per region'], '_', '-'))
% xlabel('Temperature (\circC)')
% ylabel('Occupation Probability')
% l = legend(lngd); set(l, 'color', 'k','TextColor', 'w','EdgeColor', 'k')
% 
% % set the label sizes and color
% fig = formatFig(fig, true);
% labelHandles = findall(gca, 'type', 'text','handlevisibility', 'off');
% set(labelHandles,'FontSize', 16)
% 
% % save figure
% save_figure(fig, [expRoot, 'ROI occupancy'], '-png');


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
     


%% FIGURE: temp &  roi occupancy
genotype = data(1).trialData(1).params.genotype;
nflies = 0;
for trial = 1:ntrials
    temps(:,trial) = data(trial).temp;
    nflies = data(trial).nflies + nflies;
end
x = mean(temps,2); %mean temperature

fig = getfig(''); hold on
for roi = 1:4
    lngd{roi} = Wells(roi).name;
    switch lngd{roi}
        case 'OTC'
            kolor = Color('orange');
        case 'Plant'
            kolor = Color('DeepSkyBlue');
        case 'Empty'
            kolor = Color('Snow');
        case 'Yeast'
            kolor = Color('DeepPink');
    end
    y = mean(Wells(roi).occ,2);
    plot(x,y,'linewidth', 2, 'color', kolor) 
end

% labels
title({strrep([ExpGroup ' occupation per region'], '_', '-');...
      ['Fly genotype: ' genotype];...
      ['N trials: ' num2str(ntrials) ' | N flies: ' num2str(nflies)]})
xlabel('Temperature (\circC)')
ylabel('Occupation Probability')
l = legend(lngd); set(l, 'color', 'k','TextColor', 'w','EdgeColor', 'k')

% set the label sizes and color
fig = formatFig(fig, true);
labelHandles = findall(gca, 'type', 'text','handlevisibility', 'off');
set(labelHandles,'FontSize', 16)

% save figure
figName = [baseFolder, 'Analysis\' ExpGroup ' ' genotype ' ROI occupancy'];
save_figure(fig, strrep(figName, '>','-'), '-png');






        
               
               
               
               
               
               
%%  working space...              
fig = figure; set(fig, 'color', 'k')
s = pcolor(occupancyFrames(:,:,1));
s.FaceColor = 'interp';
shading flat;


x = repmat(1:n,[n,1]);
y = repmat((1:n)',[1,n]);
z = well(roi).bin(:,:,1);


figure;
for vid = 1:nvids
    z = well(roi).bin(:,:,vid);
    surf(x,y,z)
    pause(0.1)
end

figure
surf(x,y,z)



% quick vid preview of any structure over time
figure;
for vid = 1:nvids
    img = occupancyFrames(:,:,vid);
    imagesc(img)
    pause(0.1)
end
    



Vq = interpn(well(roi).bin,3);
fig = getfig('',1);
hold on
for ii = 1:size(Vq,3)
    imagesc(Vq(:,:,ii))
    pause(0.05)
end

    
%     for vid = 1:nvids
%         data(trial).frame(:,:,vid);
%     wells.temp
%     % pull in the 


% Align video frames by temperature...
figure; hold on
for trial = 1:ntrials
    scatter(1:length(data(trial).temp), data(trial).temp, 30)
end




















