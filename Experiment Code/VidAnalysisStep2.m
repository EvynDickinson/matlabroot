
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
nvids = length(dir([dataFolder, expName '*.mat'])); % number of vids in trial
expRoot = [dataFolder, expName, '_'];

%% Load selected data
h = waitbar(0, ['loading data ' num2str(vid) '/' num2str(nvids)]);
for vid = 1:nvids
    temp = load([expRoot num2str(vid) ' data.mat']);
    trialData(vid).videoData = temp.videoData;
    trialData(vid).params = temp.params;
    waitbar(vid/nvids,h)
end
close(h); clear h

%% Fly count / pixel normalization

% Get number of flies in this trial
XLflynum = excelfile{XLrow, Excel.numflies};
if isnan(XLflynum) %condition: no known fly number 
   
    % pull up a still image from the video and manually count flies
    movieInfo = VideoReader([baseFolder folder '\' expName '_1.avi']); %read in video
    demoImg = rgb2gray(read(movieInfo,1));
    demoImg = imsharpen(imadjust(demoImg),'Radius',2,'Amount',1);% increase contrast   
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

%% Calculate flies per quadrant
% TODO: working here --> why aren't the probs summing to 1 per frame???
% extract occupancy probability per region
for vid = 1:nvids
    frame(:,:,vid) = trialData(vid).videoData.occ_Prob; % occupancy prob for this frame
    temp(vid) = mean(trialData(vid).videoData.tempLog(:,2)); %avg temp during video
    
    for roi = 1:4
        img = frame(:,:,vid); %blank image
        img(quadMask(roi).mask) = 0;
        y(vid,roi) = sum(sum(img));
    end
end


% PLOT THE CRAP OUTTA THAT!
colorList = {'green', 'white', 'yellow', 'grey'};
LW = 2;
fig = getfig;
hold on
for roi = 1:4
    plot(temp, y(:,roi), 'color', Color(colorList{roi}), 'linewidth', LW)
end


fig = formatFig(fig, 1);




sum(sum(trialData(vid).videoData.occ_Prob))













%% Save and update fly list for finished processing









%%

% % find RADII distance between 2 wells & then calc. distance and divide by two <--
% % gives the radius of desired circles
% figure;
% imshow(demoImg)
% roi = drawline;
% pos = roi.Position;
% radii = sqrt(diff(pos(:,1))^2 + diff(pos(:,2))^2)/2;



