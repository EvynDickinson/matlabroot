%% 
vidPath = 'E:\My Drive\Jeanne Lab\DATA\08.12.2021\DummyVid_20C_1.avi';
movieInfo = VideoReader(vidPath); %read in video
Img1 = read(movieInfo,1);
vidPath = 'E:\My Drive\Jeanne Lab\DATA\08.24.2021\PlantYeastChoice_N1_1.avi';
movieInfo = VideoReader(vidPath); %read in video
Img2 = read(movieInfo,1);


f = figure; set(f, 'color', 'k')
imshowpair(Img1, Img2,'montage'); title('Paint weathering')

save_figure(f, 'E:\My Drive\Jeanne Lab\DATA\Analysis\Paint weathering example', '-png');



%% get base folder pathway
baseFolder = getCloudPath;

params.well_1 = 'Plant';
params.well_2 = 'Empty';
params.well_3 = 'Empty';
params.well_4 = 'Empty';

params.genotype = select_cross; % choose the fly genotypes
params.protocol = select_protocol;  % choose experiment protocol

% move params into structure:
params.date = dirName;
params.expID = videoNames;
params.num = num;


writeExptoExcel(params) % save the data to the experiement list

% resave with same name: 
save([baseFolder dirName '\' videoNames 'dataMat'])


%% 
baseFolder = getCloudPath;

for i = 1:num.vids
    filename = [baseFolder dirName '\analysis\' videoNames '_' num2str(i) ' data'];
    load(filename);
    
    params.well_1 = 'Plant';
    
    save(filename, 'videoData', 'params')
end


%% Fix the temp log information being added...
clear
% manually load the exp data file
% Fix the MAIN file
[~,baseFolder] = getCloudPath(1);

tempLogStart = [(1:num.vids)', tempLogStart];
save([baseFolder, '\', videoNames, 'dataMat.mat'])

% load each vid file & adjust the temp log:
tempLog = readmatrix([baseFolder '\' videoNames '_RampLog']);
roi = [tempLogStart(:,3), tempLogEnd(:,3)];
filepath = [baseFolder, '\analysis\' videoNames '_'];

for i = 1:num.vids
    load([filepath num2str(i) ' data.mat']);
    
    % update file
    loc = tempLog(:,1)>=roi(i,1) & tempLog(:,1)<=roi(i,2);
    videoData.tempLog = tempLog(loc,:);
    
    % save file
    save([filepath num2str(i) ' data.mat'], 'params', 'videoData')
end

%% Video testing


Vq = interpn(tframe,3);
fig = getfig('',1);
hold on
for ii = 1:size(Vq,3)
    imagesc(Vq(:,:,ii))
    pause(0.05)
end


%% Crop out the regions around each food well to compare across trials...

figure;
mask = ~data(1).quadMask(1).mask;
imshow(mask);


for vid = 1:nvids
    Img = data(1).frame(:,:,vid); % Image
    
    
    imagesc(Img)
    
%     temp(vid) = mean(trialData(vid).videoData.tempLog(:,2)); %avg temp during video
    
    for roi = 1:4
        
        
        imgAdj = Img; % start with full image
        mask = data(1).quadMask(roi).mask;
        imgAdj(mask) = 0; % crop out frame
        %show masked image
        imagesc(imgAdj)
        
        props = regionprops(~mask, 'BoundingBox');
        maskedImage = imcrop(imgAdj, props.BoundingBox);
        % show cropped image
        imagesc(maskedImage);
        
        
        y(vid,roi) = sum(sum(imgAdj));
    end
end


props = regionprops(mask, 'BoundingBox');
maskedImage = imcrop(Img, props.BoundingBox);
imshow(maskedImage, []);


% Crop the image to the bounding box.
props = regionprops(mask, 'BoundingBox');
maskedImage = imcrop(maskedImage, props.BoundingBox);
% Display it in the lower right plot.
subplot(2, 2, 4);
imshow(maskedImage, []);


%%


%% Simple visualization: relationship between temp & well occupation
nrow = 4;
ncol = 1;
subplotInd(1).idx = 1;
subplotInd(2).idx = 2:4;
% group data across videos:
[plotX, plotY] = deal([]);
for vid = 1:nvids
    plotX = [plotX, data(vid).tempLog];
    plotY = [plotY; data(vid).well_counts];
end
% adj plant label bias:
plotY(:,3) = plotY(:,3)-3;

plotY = plotY./nflies;
sSpan = 180;
LW = 1;
time = linspace(1,(length(plotX)/3)/60,length(plotX));

fig = getfig('',1); 
    subplot(nrow,ncol,subplotInd(1).idx)
    plot(time, smooth(plotX,sSpan), 'linewidth', LW, 'color', 'w')
    ylabel('temp (\circC)')
    ylim([5,25])
    subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    for well = 1:4
        switch expData.params.(['well_' num2str(well)'])
            case 'Yeast'
                kolor = Color('gold');
            case 'Plant'
                kolor = Color('green');
            case 'Empty'
                kolor = Color('grey');
        end

        plot(time, smooth(plotY(:,well),sSpan), 'linewidth', LW, 'color', kolor);
    end
    xlabel('time (min)'); ylabel('occupation probability')
    
formatFig(fig, true, [nrow, ncol], subplotInd)
l = legend(wellLabels);
set(l, 'color', 'k', 'textcolor', 'w','position', [0.7947 0.6462 0.0963 0.1126])
subplot(nrow,ncol,subplotInd(1).idx)
set(gca, 'XColor', 'k')
titleName = strrep([folder ' ' expName], '_',' ');
title(titleName,'color', 'w')

% Save image
save_figure(fig, [analysisDir expName ' adjusted well occupation timcourse'], '-png');



