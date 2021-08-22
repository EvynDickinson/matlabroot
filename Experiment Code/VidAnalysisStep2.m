
% VidAnalysisStep2 --> group flies and plot graphs of temp & occupancy
% follows VidAnalysisStep1
clear

%% Select experiment data
% set directories
[baseFolder, folder] = getCloudPath(2);
dataFolder = [baseFolder, folder, '\analysis\'];

% select trial:
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
loc = find(strcmpi(excelfile(:,Excel.date), folder)); %rows with sel date
expList = excelfile(loc,Excel.expID); % list of trials on selected date
sel = listdlg('ListString', expList, 'ListSize', [300, 400]);
XLrow = loc(sel);   %row location in excel file for the selected trial
expName = expList{sel}; % name of selected trial
nvids = length(dir([dataFolder, expName '*data.mat'])); % number of vids in trial
expRoot = [dataFolder, expName, '_'];

%% Load selected data

if isfile([expRoot 'analysis.mat'])
    switch questdlg('Former file found, load that?')
        case 'Yes'
            load([expRoot 'analysis'])
            return
    end
else
    h = waitbar(0, 'loading videos');
    for vid = 1:nvids
        temp = load([expRoot num2str(vid) ' data.mat']);
        trialData(vid).videoData = temp.videoData;
        trialData(vid).params = temp.params;
        waitbar(vid/nvids,h)
    end
    close(h); clear h
end

%% Fly count / pixel normalization

% Get number of flies in this trial
XLflynum = excelfile{XLrow, Excel.numflies};
% pull up a still image from the video and manually count flies
movieInfo = VideoReader([baseFolder folder '\' expName '_1.avi']); %read in video
demoImg = rgb2gray(read(movieInfo,1));
demoImg = imsharpen(imadjust(demoImg),'Radius',2,'Amount',1);% increase contrast   
if isnan(XLflynum) %condition: no known fly number  
    fprintf('\n Count the number of flies in the picture by clicking them\n then hit ENTER\n')
    nflies = size(readPoints(demoImg),2);
    % write number of flies into the excel sheet
    xlswrite(xlFile, {num2str(nflies)}, 'Exp List', [Alphabet(Excel.numflies) num2str(XLrow)]);
elseif XLflynum>0
    nflies = XLflynum;
else, warndlg('Manually count fly numbers')
end

% find fly:pixel normalization
for vid = 1:nvids
    raw = (trialData(vid).videoData.currOcc);
    raw(raw>0) = 1;
    pixCount(vid) = sum(sum(raw));
end
pixRatio = (sum(pixCount)/nvids)/nflies; % rough number of pixels/fly

%% Partition pixels into regions for postional preference
radii = 200; % code at BOTTOM for how to recalc. this at a future date if cam specs change
m = size(demoImg,1);
n = size(demoImg,2);
badFit = true;

while badFit
    % Select the center of the four food wells in order 
    fprintf('\n Mark center of each well -- in order 1-4\n')
    pts = readPoints(demoImg,4);

    % draw demo circles around the wells and confirm placement
    Img = demoImg;
    fullmask = zeros(m,n);
    for ii = 1:4
        
        % Define coordinates and radius
        x1 = pts(1,ii);
        y1 = pts(2,ii);
        % Generate grid with binary mask representing the circle. 
        [xx,yy] = ndgrid((1:m)-y1,(1:n)-x1);
        Mask = (xx.^2 + yy.^2<radii^2);
        quadMask(ii).mask = ~Mask;
        fullmask(Mask) = 1;
        temp = demoImg;
        temp(~Mask) = 0;
        quadMask(ii).img = temp;
        % add to mask update:
        Img(Mask) = 0;
    end
    % offer mask fit approval:
    temp = demoImg;
    temp(~fullmask) = 0;
    f = figure;
    imshowpair(Img, temp, 'montage')
    h = questdlg('Good fit?');
    close(f)
    if strcmp(h, 'Yes')
       badFit = false;
    elseif strcmp(h, 'Cancel'); return
    end
end
clear temp

%% Calculate & plot occupancy per quadrant
% TODO: working here --> why aren't the probs summing to 1 per frame???
% extract occupancy probability per region
for vid = 1:nvids
    frame(:,:,vid) = trialData(vid).videoData.occ_Prob/...
                     trialData(vid).params.num.framesPerVid; % occupancy prob for this frame
    temp(vid) = mean(trialData(vid).videoData.tempLog(:,2)); %avg temp during video
    
    for roi = 1:4
        img = frame(:,:,vid); %blank image
        img(quadMask(roi).mask) = 0;
        y(vid,roi) = sum(sum(img));
    end
end

% PLOT THE CRAP OUTTA THAT OCCUPANCY!
LW = 2;
fig = getfig('',1);
hold on
for roi = 1:4
    lngd{roi} = trialData(1).params.(['well_' num2str(roi)]);
    switch lngd{roi}
        case 'OTC'
            kolor = Color('orange');
        case 'Plant'
            kolor = Color('aqua');
        case 'Empty'
            kolor = Color('grey');
        case 'Yeast'
            kolor = Color('fuchsia');
    end
    plot(flip(temp), flip(y(:,roi)), 'color', kolor, 'linewidth', LW)
end
xlim([8,22])

% labels
title(strrep([expName ' occupation per region'], '_', '-'))
xlabel('Temperature (\circC)')
ylabel('Occupation Probability')
l = legend(lngd); set(l, 'color', 'k','TextColor', 'w','EdgeColor', 'k')

% set the label sizes and color
fig = formatFig(fig, true);
labelHandles = findall(gca, 'type', 'text','handlevisibility', 'off');
set(labelHandles,'FontSize', 16)

% save figure
save_figure(fig, [expRoot, 'ROI occupancy'], '-png');


%% In progress ... 

%% Optional: Write a video of the occupational probability over the experiment

% set(0,'DefaultFigureVisible','off');
% Set up recording video parameters
vid_name = [expRoot(1:end-1) ' occupancy probability vid'];
FrameRate = 3; %trialData(1).params.num.fps*2; %playback rate
nframes = nvids;
for ii = 1:nvids
    tframe(:,:,ii) = imresize(frame(:,:,ii),[20,20]);
end
cRange = [min(min(min(tframe))), max(max(max(tframe)))]; % occupancy prob range
dummyFrame = cRange(1) * ones(size(tframe(:,:,1),1,2));

% set up plot spacing
subplotInd(1).idx = 1:9; % occupation plot
subplotInd(2).idx = 10:12; % temperature bar
nrow = 4;
ncol = 3;


v = VideoWriter(vid_name, 'Motion JPEG AVI');
v.Quality = 90;
v.FrameRate = FrameRate;
open(v);
fig = getfig('',1); set(fig, 'pos', [-1044 261 788 724]); %update this position for other computers


tic
h = waitbar(0,'reading & writing frames');
for ii = 1:nvids

    subplot(nrow, ncol, subplotInd(2).idx); hold on
        plot(time, temp, 'linewidth', 2, 'color', 'w')
        xlabel('Time (min)')
        ylabel('Temp \circC')
        ylim([5,25])
        xlim(xlimit)
        vline((ii*num.vidLength)-num.vidLength,'y')
        formatFig(fig, true);

    subplot(nrow, ncol, subplotInd(1).idx); hold on
        imagesc(tframe(:,:,ii),cRange); 
        colormap hot;
        c = colorbar;
        currTemp = floor(temp(ii));
        title(['Occupation probability at ' num2str(currTemp) '\circC']);
        formatFig(fig, true);
        set(gca,'xtick',[],'ytick',[],'xcolor', 'k', 'ycolor', 'k')
        set(c,'Color', 'w', 'FontName', 'Arial');

    % collect frame info:
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset') 
    waitbar(ii/nvids,h)
end
close(h)
close(v)
toc 

fprintf('\n Video Saved! \n')  
    
    %%
% 
%     % make a video with a temp bar at the bottom....
% num = trialData(1).params.num;
% time = 0:num.vidLength:num.sessionLength-num.vidLength;
% xlimit = [time(1)-num.vidLength, time(end)+num.vidLength];
% 

% 
% % 
% 
% fig = getfig('',1); set(fig, 'pos', [-1044 261 788 724]);
% subplot(nrow, ncol, subplotInd(2).idx)
% hold on
% plot(time, temp, 'linewidth', 2, 'color', 'w')
% xlabel('Time (min)')
% ylabel('Temp \circC')
% ylim([5,25])
% xlim(xlimit)
% vline(ii*num.vidLength,'y')
% formatFig(fig, true);
% 
% subplot(nrow, ncol, subplotInd(1).idx)
% h = imagesc(frame(:,:,ii),cRange); 
%     colormap hot;
%     c = colorbar;
%     currTemp = floor(temp(ii));
%     title(['Occupation probability at ' num2str(currTemp) '\circC']);
%     formatFig(fig, true);
%     set(gca,'xtick',[],'ytick',[],'xcolor', 'k', 'ycolor', 'k')
%     set(c,'Color', 'w', 'FontName', 'Arial');



    
    


fig = getfig; 
    hold on
    v = VideoWriter(vid_name, 'Motion JPEG AVI');
    v.Quality = 90;
    v.FrameRate = FrameRate;
    open(v);
    set(fig, 'Color', 'k')   
    % add scale bar
    imshow(zeros(txt.frameSize))
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset') 
 

%load all the data for the video & early process images:
n_tot = nframes;
tot_occ = (zeros(height, width));

tic
h = waitbar(0,'reading & writing frames');
for i = 1:n_tot 
%     % Read in current frame from raw video
%     frame = (rgb2gray(read(movieInfo,i)));
%     BWframe = imbinarize(im2double(frame),...
%               'adaptive','ForegroundPolarity','bright','Sensitivity',0.3);
%     [currProb, tot_occ] = (getOccProb(tot_occ,BWframe,i)); % spatial probabilty
%     currProb(currProb==1) = 0; %filter out light particles
    

    % load previous frame
    data = load([analysisDir vidName '_' num2str(i) ' data.mat']);
    
    
    % Build an image frame
    frame = videoData.occ_Prob/params.num.framesPerVid;
    
  
    % --- log spatial distribution ---
    fig = getfig('',1); 
  
        h = heatmap(frame); 
        colormap hot;
        set(h,'gridvisible', 'off')
        ax = gca;
        ax.XDisplayLabels = nan(size(ax.XDisplayData));
        ax.YDisplayLabels = nan(size(ax.YDisplayData));
%         ax.ColorScaling = 'log';
        set(ax,'FontColor', 'w', 'FontName', 'Arial');
        title(['Occupation probability ' num2str(2) '\circC']);
        
        formatFig(fig, true)
    subplot(1,2,2)
        h = heatmap(tot_probabilty); 
        colormap hot;
        set(h,'gridvisible', 'off')
        ax = gca;
        ax.XDisplayLabels = nan(size(ax.XDisplayData));
        ax.YDisplayLabels = nan(size(ax.YDisplayData));
        set(ax,'FontColor', 'w', 'FontName', 'Arial');
        title('Spatial occupation probability');


    
    
    
    
    
    
    
    newFrame = OccBuildFrame(im2double(frame), BWframe, currProb, txt);
    % Collect image info for movie file
    set(fig, 'Color', 'k')
    imshow(newFrame)
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')
    
    waitbar(i/n_tot,h)
end
close(h)
close(v)
toc 

fprintf('\n Video Saved! \n')  



%% Save and update fly list for finished processing

save([expRoot 'analysis'])


%%

% % find RADII distance between 2 wells & then calc. distance and divide by two <--
% % gives the radius of desired circles
% figure;
% imshow(demoImg)
% roi = drawline;
% pos = roi.Position;
% radii = sqrt(diff(pos(:,1))^2 + diff(pos(:,2))^2)/2;

















