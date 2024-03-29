
clear; close all; clc
autoSave = true;
essentialfigs = true;
%% ---------------------Select Date & Experiment to Process ------------------------------

%load excel file:
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
switch questdlg('Use Excel sheet to select experiment?')
    case 'Yes'
        loc1 = cellfun(@isnan,excelfile(2:end,Excel.numflies));
        loc2 = cellfun(@ischar,excelfile(2:end,Excel.tracked));
        loc = loc1 & loc2;
        rownums = find(loc)+1;
        if isempty(rownums)
            warndlg('No available experiments')
            return
        end
        eligible_files = excelfile([false;loc],[Excel.date, Excel.arena, Excel.expID]);
        FileNames = join(eligible_files);
        fileIdx = listdlg('ListString', FileNames,'ListSize',[250,450]);   
        if isempty(fileIdx); return; end
        % get file info:
        baseFolder = getCloudPath;
        folder = eligible_files{fileIdx,1};
        arenaSel = eligible_files{fileIdx,2};
        expName = eligible_files{fileIdx,3};
        clear loc1 loc2 loc eligible_files FileNames rownums fileIdx
    case 'No'
        %get base folder pathway
        [baseFolder, folder] = getCloudPath(2);    
        %select arena to work with:
        arenaList = {'A', 'B', 'C', 'D'};
        arenaSel = arenaList{listdlg('ListString', arenaList)};
        % Select the complete experiments to process
        list_dirs = dir([baseFolder folder, '\*dataMat.mat']); %only matlab files
        list_dirs = {list_dirs(:).name};
        expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
        expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
        expName = expName(1:end-1);clear expNames
    case 'Cancel'
        return
end

% Saving and Loading Directories:
vidFolder = [folder '\Arena ' arenaSel];
analysisDir = [baseFolder vidFolder '\analysis\'];
if ~isfolder(analysisDir); mkdir(analysisDir); end
expPDF = [analysisDir folder expName arenaSel ' summary.pdf'];
XLrow = find(strcmpi(excelfile(:,Excel.date), folder) & ...
             strcmpi(excelfile(:,Excel.expID), expName) & ...
             strcmpi(excelfile(:,Excel.arena), arenaSel(end:end))); %rows with sel date with ARENA ID

% Load relevant data files (.mat, .csv, .h5)
warning off 
expData = load([baseFolder folder '\' expName ' dataMat.mat']);
tempLog = readmatrix([baseFolder folder '\' expName '_RampLog']);
nvids = expData.parameters.numVids;

% load tracking predictions     
data = struct;
for vid = 1:nvids
    filePath = [baseFolder vidFolder '\' expName '_' num2str(vid) '.h5'];
    data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
    data(vid).tracks = h5read(filePath,'/tracks');
end; clear filePath

initial_vars = who; initial_vars{end+1} = 'initial_vars';
fprintf('\nNext\n')

%% --------------------- Parameter extraction ----------------------------------------
% Number of flies:
nflies = excelfile{XLrow,Excel.numflies};
movieInfo = VideoReader([baseFolder vidFolder '\' expName   '_1.avi']); %read in video
ii = randi(size(data(1).occupancy_matrix,2),[3,1]); %random selection of frames to count fly number
demoImg = rgb2gray(read(movieInfo,1));
if isnan(nflies)
    % manual count of flies
    fprintf('\nCount the number of flies in the picture by clicking them\n then hit ENTER\n')
    T = true;
    while T % get the number of flies in the arena
        for jj = 1:3
            demoImg = rgb2gray(read(movieInfo,ii(jj)));
            nflies(jj) = size(readPoints(demoImg),2);
        end
        fprintf(['\nNumber of flies counted: ' num2str(nflies) '\n'])
        if sum(~diff(nflies)==0)==0
            nflies = nflies(1); 
            T = false;
        else
            switch questdlg(['Nonmatching fly counts, manually select number?: ' num2str(nflies)])
                case 'Yes'
                    nflies = str2double(inputdlg('Input number of flies'));
                    T = false;
                case 'No'
                    T = true;
                case 'Cancel'
                    return
            end
        end
    end
    % write number of flies into the excel sheet
    try
        xlswrite(xlFile, {num2str(nflies)}, 'Exp List', [Alphabet(Excel.numflies) num2str(XLrow)]);
    catch
        h = warndlg('Close Experiment Summary excel file and then close this warning box');
        uiwait(h)
        xlswrite(xlFile, {num2str(nflies)}, 'Exp List', [Alphabet(Excel.numflies) num2str(XLrow)]);
    end
end
fprintf(['\nNumber of flies: ' num2str(nflies) '\n'])

% Find food wells and their identities %TODO
h = warndlg('Select the well locations starting at "9 o''clock" and proceeding clock wise'); uiwait(h)
% label the wells:
wellLabels = {expData.parameters.(['Arena' arenaSel]).well_1;...
              expData.parameters.(['Arena' arenaSel]).well_2;...
              expData.parameters.(['Arena' arenaSel]).well_3;...
              expData.parameters.(['Arena' arenaSel]).well_4};   
disp(wellLabels)
wellcenters = readPoints(demoImg,4); % get locations of the wells from the image

initial_vars = [initial_vars; 'nflies'; 'wellLabels'; 'wellcenters'; 'demoImg'];
clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% ----------------------------- Data organization by video ---------------------------
% find number of flies per frame
currFrame = 0;
for vid = 1:nvids
    occupancy.frameROI(vid,1) = currFrame+1;
    nframes = size(data(vid).occupancy_matrix,2);
    currFrame = nframes + currFrame;
    occupancy.frameROI(vid,2) = currFrame;
end

%% --------------------------- Mask the food wells & arena ---------------------------
maskPath = [analysisDir expName arenaSel ' ArenaMask.mat'];
if ~isfile(maskPath)
    idx = find(~strcmpi('Empty', wellLabels));
    if ~isempty(idx)
        % draw out the masks:
        for ii = 1:length(idx)
            wellName = strrep(wellLabels{idx(ii)}, '_', '-');
            % label the ROI
            f = figure;
            imshow(demoImg)
            title(['Outline the ' wellName ' well'])
            roi = drawpolygon;
            mask(ii).n = roi.Position;
            uiwait(f)
        end
        % save the food masks:
        save(maskPath, 'mask')
    end
else
    load(maskPath)
end

% Check the video size
if all(size(demoImg)==946)
    r = 400; %radius of the arena for the quadbowl
else
    warndlg('Not expected video pixel size...manually set radius')
    return
    
    % manually select an arena radius?
    
    % % pic preview of mask circle:
    % figure; hold on
    % imshow(demoImg)
    % viscircles(centre', r);
end

% Screen any points outside the arena
x1 = wellcenters(1,1:2:4);
y1 = wellcenters(2,1:2:4);
x2 = wellcenters(1,2:2:4);
y2 = wellcenters(2,2:2:4);
[xi,yi] = polyxpoly(x1,y1,x2,y2);
centre = [xi;yi];
% 
% % pic preview of mask circle:
% figure; hold on
% imshow(demoImg)
% viscircles(centre', r);

% Tracking matrix locations: [frame, node, xy, fly]
scaleFactor = 1.31; % WHY IS THIS NECESSARY??? TODO

% MASK FOOD AND ARENA:
[raw_flyCount, mid_flyCount] = deal([]);
for vid = 1:nvids
    % number of flies labeled on each frame:
    raw_flyCount = [raw_flyCount; sum(data(vid).occupancy_matrix)'];
    
    % all data for head tracked location
    headData = squeeze(data(vid).tracks(:,2,:,:));
%     headData = squeeze(data(vid).tracks(:,1,:,:)); % TODO REVERT TO HEADPOINT
    %save the unadjusted data in the same format
    data(vid).x_loc_raw = squeeze(headData(:,1,:)).*scaleFactor;
    data(vid).y_loc_raw = squeeze(headData(:,2,:)).*scaleFactor;
    nframes = size(headData,1); 
    
    % x-y coordinates of flies for each frame
    x_loc = squeeze(headData(:,1,:)).*scaleFactor;
    y_loc = squeeze(headData(:,2,:)).*scaleFactor;

    % REMOVE FOOD TRACKED POINTS HERE 
    Xdim = size(x_loc);
    Ydim = size(y_loc);
    %resize the data:
    X = reshape(x_loc,[numel(x_loc),1]);
    Y = reshape(y_loc,[numel(y_loc),1]);
    % Find points within the masked region and turn to NaN
    if exist('mask','var')
        for ii = 1:length(mask)
            [in,on] = inpolygon(X,Y, mask(ii).n(:,1),mask(ii).n(:,2));   % Logical Matrix
            inon = in | on;                                    % Combine ‘in’ And ‘on’
            X(inon) = NaN;
            Y(inon) = NaN;
        end
    end
    % Remove any labeled points outside the arena:
    loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)>=r;
    X(loc) = NaN; Y(loc) = NaN;

    % Resize the data to OG structure:
    X = reshape(X,Xdim);
    Y = reshape(Y,Ydim);

    data(vid).x_loc_mid = X; % save for later convenience
    data(vid).y_loc_mid = Y;
    
    % New fly count measure:
    mid_flyCount = [mid_flyCount; sum(~isnan(X),2)];
end

%% --------------------------- Mask the overtracked points ---------------------------
removalRadius = 10;
% Find the video with the greatest avg over-count of fly points:
for vid = 1:nvids
    ROI = occupancy.frameROI(vid,:);
    numberCount(vid) = median(mid_flyCount(ROI(1):ROI(2)));
end
% find the top 4 frames and display them
ptsPath = [analysisDir expName arenaSel ' FalsePoints.mat'];
if ~isfile(ptsPath)
    pullFrames = 5;
    frameList = [];
    vid_idx = find(numberCount == max(numberCount)); vid_idx = vid_idx(1);
    ROI = occupancy.frameROI(vid_idx,:);
    numlist = mid_flyCount(ROI(1):ROI(2));
    [~,Idx] = (sort(numlist));
    frameList = Idx(end-pullFrames+1:end); %pull the four highest overcounts

    % pull up the images in order and click on points that ARE NOT FLIES:
    movieInfo = VideoReader([baseFolder vidFolder '\' expName '_' num2str(vid_idx) '.avi']); %read in video
    for ii = 1:length(frameList)
        frame = frameList(ii);
        img = read(movieInfo,frame);
        fig = figure;
        imshow(img);
            hold on
            x = data(vid_idx).x_loc_mid(frame,:); x(isnan(x)) = [];
            y = data(vid_idx).y_loc_mid(frame,:); y(isnan(y)) = [];
        scatter(x,y, 10, 'y')
        title(['select all points that are NOT flies ' num2str(ii) '/' num2str(pullFrames)])
        pointLabels(ii).coord = labelWrongPoints(fig);
    end

    % Now determine how they are related and if there is consistent overlap
    % that can be targeted for deletion...
    fig = figure; set(fig, 'color', 'k')
    imshow(img)
    hold on
    for ii = 1:length(frameList)
        scatter(pointLabels(ii).coord(1,:), pointLabels(ii).coord(2,:),20, 'filled')
    end
    title('Select points with overlap for masking')
    % save_figure(fig, [analysisDir expName arenaSel ' overtracking point IDs'], '-png');
    pointLabels(1).finalRound = labelWrongPoints(fig);
    nmaskpoints = length(frameList);

    % visualize the selected ROIs 
    fig = figure; set(fig, 'color', 'k');
    imshow(img)
    hold on
    for ii = 1:nmaskpoints
        scatter(pointLabels(ii).coord(1,:), pointLabels(ii).coord(2,:),20, 'filled')
    end
    R = removalRadius*ones(size(pointLabels(1).finalRound,2),1);
    viscircles(pointLabels(1).finalRound', R ,'Color','r');
    save_figure(fig,[analysisDir,expName,arenaSel,' overtracking points selected for deletion'], '-png');
    
    % save the food masks:
    save(ptsPath, 'frameList','vid_idx', 'pointLabels','nmaskpoints')
else
    load(ptsPath)
end

% how many flies are in those circles?? TODO integrate the last OG fly count...
flyCount = [];
for vid = 1:nvids
%resize the data:
    x_loc = data(vid).x_loc_mid; 
    y_loc = data(vid).y_loc_mid;
    
    % how many flies OG
    Xdim = size(x_loc); 
    Ydim = size(y_loc);
    X = reshape(x_loc,[numel(x_loc),1]); 
    Y = reshape(y_loc,[numel(y_loc),1]);

    % Find points within the masked region and turn to NaN
    if ~isempty(pointLabels(1).finalRound)
        for ii = 1:size(pointLabels(1).finalRound,2)
            loc = (((X-pointLabels(1).finalRound(1,ii)).^2 + (Y-pointLabels(1).finalRound(2,ii)).^2).^0.5)<=removalRadius;
            X(loc) = NaN; Y(loc) = NaN;
        end
    end
    % Resize the data to OG structure:
    X = reshape(X,Xdim);
    Y = reshape(Y,Ydim);
    flyCount = [flyCount; sum(~isnan(X),2)];
    
    data(vid).x_loc = X; % save for later convience
    data(vid).y_loc = Y;
end

% save long-term fly count data:
occupancy.time = linspace(1,(length(mid_flyCount)/3)/60,length(mid_flyCount));
occupancy.raw_flyCount = raw_flyCount;
occupancy.mid_flyCount = mid_flyCount;
occupancy.flycount = flyCount;

%% Temp log pull:

% Find number of flies that are near each well for each frame
full_temp = [];
for vid = 1:nvids
    radii = 165; %110; %distance must be less than this number to count for a well ROI 
    dims = size(data(vid).occupancy_matrix);
    % loop for all wells 
    for well = 1:4
        b = wellcenters(:,well); 
        % reshape data points from whole video for optimized maths
        x_loc = reshape(data(vid).x_loc,numel(data(vid).y_loc),1); 
        y_loc = reshape(data(vid).y_loc,numel(data(vid).y_loc),1); 
        % within well range?
        euDist = sqrt((b(1)-x_loc).^2 + (b(2)-y_loc).^2); %euclidian distance from well center
        well_loc = reshape((euDist<=radii),[dims(2),dims(1)]);
        welldata(well).loc = well_loc;
        % how many within the circle?
        welldata(well).count = sum(well_loc,2);
        data(vid).well_counts(:,well) = welldata(well).count;
    end
    % temperature alignment
    logROI(1) = find(tempLog(:,1)==expData.tempLogStart(vid,3));
    logROI(2) = find(tempLog(:,1)==expData.tempLogEnd(vid,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    x = round(linspace(1, dims(2), length(tempCourse)));
    % upsample the temperature log:
    fullTempList = interp1(x,tempCourse,1:dims(2),'spline');   
    data(vid).tempLog = fullTempList;
    full_temp = [full_temp, data(vid).tempLog]; % temperature
end

%%  ------------------------------ Visualization --------------------------------------
%read in video
movieInfo = VideoReader([baseFolder,vidFolder,'\',expName,'_',num2str(vid_idx),'.avi']); 
frame = frameList(1);
img = read(movieInfo,frame);
vid = vid_idx;
nrows = 4;
ncols = 3;
sb(1).idx = [1,2,4,5]; %imshow image
sb(2).idx = [7,8,10,11]; %time course frame count
sb(3).idx = 3:3:12; %histogram
    
fig = figure; set(fig, 'pos', [300 100 1171 826]); % X-off, Y-off, width, height
% ARENA IMAGE
subplot(nrows,ncols,sb(1).idx)
    imshow(img)
    axis tight square
    hold on
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        scatter(wellcenters(1,well),wellcenters(2,well), 75,...
            'MarkerFaceColor', kolor, 'MarkerEdgeColor', 'w') 
    end 
    % raw points:
    x = data(vid).x_loc_raw(frame,:);
    y = data(vid).y_loc_raw(frame,:);
    x(isnan(x)) = []; % remove empty tracks
    y(isnan(y)) = [];
    scatter(x,y, 20, 'k')
    scatter(x,y, 15, Color('grey'),'filled')
    
    % intermediate points:
    x = data(vid).x_loc_mid(frame,:);
    y = data(vid).y_loc_mid(frame,:);
    x(isnan(x)) = []; % remove empty tracks
    y(isnan(y)) = [];
    scatter(x,y, 15, Color('orangered'),'filled')
    
    % intermediate points:
    x = data(vid).x_loc(frame,:);
    y = data(vid).y_loc(frame,:);
    x(isnan(x)) = []; % remove empty tracks
    y(isnan(y)) = [];
    scatter(x,y, 15, Color('teal'),'filled')

% OVERTRACKING OVER TIME
subplot(nrows,ncols,sb(2).idx)
    hold on
    plot(occupancy.time, raw_flyCount, 'color', Color('grey'))
    plot(occupancy.time, mid_flyCount,'color', Color('orangered'))
    plot(occupancy.time, flyCount, 'color', Color('teal'))
    hline(nflies, 'w')
    ylabel('Fly count'); xlabel('time (min)')
    axis tight
    yyaxis right
    plot(occupancy.time,full_temp,'color', 'y')
    ylabel('temp (\circC)')
    xlim([0,occupancy.time(end)+10])

% FLY COUNT HISTOGRAM
subplot(nrows,ncols,sb(3).idx)
    yyaxis right
    hold on
    h = histogram(raw_flyCount);
    h.FaceColor = Color('grey');
    h = histogram(mid_flyCount);
    h.FaceColor = Color('orangered'); 
    h = histogram(flyCount);
    h.FaceColor = Color('teal');
    vline(nflies, 'w')
    xlabel('Number of flies')
    ylabel('Frame count')

% LABELS AND FORMATTING
fig = formatFig(fig, true, [nrows, ncols], sb);
l = legend({['SLEAP ' num2str(mean(raw_flyCount)) ' avg'],...
            ['wells & arena masks ' num2str(mean(mid_flyCount)) ' avg'],...
            ['wells, arena & manual ' num2str(mean(flyCount)) ' avg']});
set(l,'color','k','textcolor','w','edgecolor','k','Position', [0.0169 0.8981 0.3122 0.0877]);
subplot(nrows,ncols,sb(2).idx)
yyaxis left
set(gca, 'ycolor', 'w')
subplot(nrows,ncols,sb(3).idx)
yyaxis left
set(gca, 'ycolor', 'k')

% save and export figure:
if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
end  
save_figure(fig, [analysisDir expName arenaSel ' quality control'], '-png');


% ----------------------------- zoom-in-take on the tracking correction -------------
% % visual confirmation that the selected points are near the well:
% AllPoints = [];
% for vid = 1:nvids
%    X = reshape(data(vid).x_loc,numel(data(vid).x_loc),1); 
%    Y = reshape(data(vid).y_loc,numel(data(vid).y_loc),1); 
%    X(isnan(X)) = [];
%    Y(isnan(Y)) = []; 4
%    AllPoints = [AllPoints; X , Y];
% end
% fig = getfig; set(fig, 'color', 'k');
% hist2d(AllPoints(:,1),AllPoints(:,2), 'probability', 'tile')
% axis tight; axis square
% set(gca, 'visible', 'off')
% c = colorbar;
% c.Color = [1,1,1];
% hold on 
% % c = drawcircle('Center', wellcenters(:,5)', 'Radius', r);
% % hold on
% % viscircles(wellcenters(:,1:4)',[radii,radii,radii,radii])
% % b = wellcenters(:,well); % well 1 points
% save_figure(fig, [analysisDir expName arenaSel ' cropped out food locations quality control'], '-png');

%% Temp vs fly count

fig = figure;
idx = discretize(full_temp,1:30);
boxplot(flyCount, idx,'Plotstyle', 'traditional','colors', 'w')
h_line(nflies, 'Teal', '-',3)
xlabel('Temperature (\circC)')
ylabel('Tracking Fly Count')
fig = formatFig(fig, true);

save_figure(fig, [analysisDir expName arenaSel ' fly count v temperature'], '-png',true);

%% ------------------- Save preformatted data for QuadStep2 ------------------------
initial_vars = [initial_vars; 'radii'; 'welldata'; 'flyCount'; 'occupancy'];
clearvars('-except',initial_vars{:})
save([analysisDir expName arenaSel ' preformed data'])
disp('Formatted data saved')
disp('Done')

%% Older version (auto-video divide (pre VirtualDub))
% %% Select data file to make new folders for
% %get base folder pathway
% [baseFolder, folder] = getCloudPath(2); 
% 
% % Make a new video folder:
% videosDirA = [baseFolder folder '\Arena A\'];
% if ~isfolder(videosDirA); mkdir(videosDirA); end
% 
% videosDirB = [baseFolder folder '\Arena B\'];
% if ~isfolder(videosDirB); mkdir(videosDirB); end
% 
% videosDirC = [baseFolder folder '\Arena C\'];
% if ~isfolder(videosDirC); mkdir(videosDirC); end
% 
% videosDirD = [baseFolder folder '\Arena D\'];
% if ~isfolder(videosDirD); mkdir(videosDirD); end
% 
% 
% %% Visual comparison of video cropping:
% 
% fig = cropCheck;
% 
% 
% %%
% 
% 
% % LOAD VIDEOS FOR AN EXPERIMENT AND CROP THEM INTO THE FOUR INDEPENDENT
% % QUADRANTS / ARENAS FOR TRACKING AND POST PROCESSING
% clear; clc  
%  
% %% Select data to split:
% %get base folder pathway
% [baseFolder, folder] = getCloudPath(2); 
% 
% 
% % Select the complete experiments to process
% list_dirs = dir([baseFolder folder, '\*.mat']); %only matlab files
% list_dirs = {list_dirs(:).name};
% expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
% expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
% expName = expName(1:end-1);
% clear expNames
% % Pull fly summary sheet information on selected experiment
% [excelfile, Excel, xlFile] = load_QuadBowlExperiments;
% 
% % Make a new video folder:
% videosDir = [baseFolder folder '\Split Videos\'];
% if ~isfolder(videosDir); mkdir(videosDir); end
% 
% %% Assign quadrants to split
% %Load vid template and manually align vid split
% params = load([baseFolder folder '\' expName ' parameters.mat']);
% nvids = params.parameters.numVids;
% movieInfo = VideoReader([baseFolder folder '\' expName '_1.avi']); %read in video
% demoImg = rgb2gray(read(movieInfo,1));
% % figure; imshow(demoImg)
% % [J,rect] = imcrop(demoImg);
% 
% % use fixed crop size for all images and videos
% cropWidth = 946;
% width = 1920;
% % RGB = insertShape(demoImg,'rectangle',[0,0,cropWidth, cropWidth]);
% % imshow(RGB)
% % cropping indexes
% roiA = [width-cropWidth, width, 1, cropWidth];
% roiB = [1, cropWidth, 1, cropWidth];
% roiC = [width-cropWidth, width, width-cropWidth, width];
% roiD = [1, cropWidth, width-cropWidth, width];
% 
% % crop the image (will be used for later processing)
% imgA = demoImg(roiA(1):roiA(2), roiA(3):roiA(4));
% imgB = demoImg(roiB(1):roiB(2), roiB(3):roiB(4));
% imgC = demoImg(roiC(1):roiC(2), roiC(3):roiC(4));
% imgD = demoImg(roiD(1):roiD(2), roiD(3):roiD(4));
% 
% % Preview selecion:
% nrow = 2; ncol = 4;
% fig = figure; set(fig, 'pos', [267 233 959 420], 'color', 'k');
% hold on;
% subplot(nrow, ncol, 1) %Arena B
% imshow(imgB)
% subplot(nrow, ncol, 2) %Arena D
% imshow(imgD)
% subplot(nrow, ncol, 5) %Arena A
% imshow(imgA)
% subplot(nrow, ncol, 6) %Arena C
% imshow(imgC)
% subplot(nrow, ncol, [3:4,7:8]) %Arena B
% imshow(demoImg)
% save_figure(fig, [baseFolder folder '\' expName ' Video Cropping'], '-pdf')
% 
% switch questdlg('Acceptable quadrants?')
%     case 'Yes'
%     case 'No'
%         return
%     case 'Cancel'
%         return
% end
% 
% %% Split the videos
% switch questdlg('Group AB or Group CD?','', 'AB', 'CD','Cancel') %
%     case 'AB'
%         tic
%         for vid = 1:nvids
%             movieInfo = VideoReader([baseFolder folder '\' expName '_' num2str(vid) '.avi']); %read in video
%             nframes = movieInfo.NumFrames;
%             h = waitbar(0,['Cropping video ' num2str(vid) '/' num2str(nvids)]);
%             % initiate arena A writer
%             videoA = [baseFolder folder '\Split Videos\' expName 'A_' num2str(vid) '.avi'];
%             vidA = VideoWriter(videoA, 'Motion JPEG AVI');
%             vidA.Quality = 75;
%             vidA.FrameRate = 6;
%             open(vidA);
%             % initiate arena B writer
%             videoB = [baseFolder folder '\Split Videos\' expName 'B_' num2str(vid) '.avi'];
%             vidB = VideoWriter(videoB, 'Motion JPEG AVI');
%             vidB.Quality = 75;
%             vidB.FrameRate = 6;
%             open(vidB);
%             fig = figure; set(fig, 'color', 'k','visible', 'off');
% 
%             for frame = 1:nframes
%                 %Extract and crop frames
%                 demoImg = rgb2gray(read(movieInfo,frame));
%                 imgA = demoImg(roiA(1):roiA(2), roiA(3):roiA(4));
%                 imgB = demoImg(roiB(1):roiB(2), roiB(3):roiB(4));
% 
%                 %Save A frames 
%                 imshow(imgA)
%                 f = getframe(fig);
%                 writeVideo(vidA, f)  
% 
%                 %Save B frames
%                 imshow(imgB)
%                 f = getframe(fig);
%                 writeVideo(vidB, f)  
% 
%                 % update the waitbar every 10 frames
%                 if rem(frame,10) == 0
%                     waitbar(frame/nframes,h)
%                 end   
%             end
%             close(vidA); close(vidB)
%              close(h); close(fig)   
%         end
%         toc
%     case 'CD'
%         tic
%         for vid = 1:nvids
%             movieInfo = VideoReader([baseFolder folder '\' expName '_' num2str(vid) '.avi']); %read in video
%             nframes = movieInfo.duration;
%             h = waitbar(0,['Cropping video ' num2str(vid) '/' num2str(nvids)]);
% 
%             % initiate arena C writer
%             videoC = [baseFolder folder '\Split Videos\' expName 'C_' num2str(vid) '.avi'];
%             vidC = VideoWriter(videoC, 'Motion JPEG AVI');
%             vidC.Quality = 75;
%             vidC.FrameRate = 6;
%             open(vidC);
%             % initiate arena D writer
%             videoD = [baseFolder folder '\Split Videos\' expName 'D_' num2str(vid) '.avi'];
%             vidD = VideoWriter(videoD, 'Motion JPEG AVI');
%             vidD.Quality = 75;
%             vidD.FrameRate = 6;
%             open(vidD);
% 
%             fig = figure; set(fig, 'color', 'k','visible', 'off');
%             for frame = 1:nframes
%                 %Extract and crop frames
%                 demoImg = rgb2gray(read(movieInfo,frame));
%                 imgC = demoImg(roiC(1):roiC(2), roiC(3):roiC(4));
%                 imgD = demoImg(roiD(1):roiD(2), roiD(3):roiD(4));
%   
%                 %Save C frames
%                 imshow(imgC)
%                 f = getframe(fig);
%                 writeVideo(vidC, f)        
%     
%                 %Save D frames
%                 imshow(imgD)
%                 f = getframe(fig);
%                 writeVideo(vidD, f)  
% 
%                 % update the waitbar every 10 frames
%                 if rem(frame,10) == 0
%                     waitbar(frame/nframes,h)
%                 end   
%             end
%             close(vidC); close(vidD)
%              close(h); close(fig)   
%         end
%         toc
% end
