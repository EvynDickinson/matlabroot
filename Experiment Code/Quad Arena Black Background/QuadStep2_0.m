
% TODO: automate so it can loop through 'autosaving' after the initial manual work

%% -------- Find the files that haven't been analyized yet and run them -------------
clear; close all; clc
autoSave = true;
essentialfigs = true; 
excelWrite = true;

%load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;
loc = cellfun(@isnan,excelfile(2:end,Excel.numflies));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.arena, Excel.expID, Excel.processed]);
loc1 = cellfun(@isnan,eligible_files(:,4));
c = cellfun(@string,eligible_files);
c(loc1,4) = ' ';
FileNames = join(c);
fileIdx = listdlg('ListString', FileNames,'ListSize',[250,450]);
%pull the list of dates and arenas to be 
List.date = eligible_files(fileIdx,1);
List.arena = eligible_files(fileIdx,2);
List.expID = eligible_files(fileIdx,3); 

%get base folder pathway
baseFolder = getCloudPath;
for ii = 1:length(fileIdx)
    inputPath = [baseFolder List.date{ii} '/Arena ' List.arena{ii} '/analysis/'...
                 List.expID{ii} List.arena{ii} ' preformed data.mat'];
%     results = runQuadStep2_excludePoorFrames(inputPath,autoSave,essentialfigs);
    results = runQuadStep2(inputPath,autoSave,essentialfigs); % Run the basic figures
    if excelWrite == true
        if strcmpi(results, 'Saved Data')
            XLrow = rownums(fileIdx(ii));
            % write number of flies into the excel sheet
            try
                xlswrite(XL, {'Y'}, 'Exp List', [Alphabet(Excel.processed) num2str(XLrow)]);
            catch
                h = warndlg('Close Experiment Summary excel file and then close this warning box');
                uiwait(h)
                xlswrite(XL, {'Y'}, 'Exp List', [Alphabet(Excel.processed) num2str(XLrow)]);
            end
        end
    end
disp(['Finished ' FileNames(fileIdx(ii))])         
end

return % prevent auto running the whole script

%% Select Date & Experiment to Process
clear; close all; clc
autoSave = true;
essentialfigs = true;
excelWrite = true;

%get base folder pathway
[baseFolder, folder] = getCloudPath(2);    
%select arena to work with:
arenaList = {'A', 'B', 'C', 'D'};
arenaSel = arenaList{listdlg('ListString', arenaList)};
% vidFolder = [folder '\Split Videos']; % old  
vidFolder = [folder '\Arena ' arenaSel];

% Select the complete experiments to process
list_dirs = dir([baseFolder folder, '\*dataMat.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
expName = expName(1:end-1);clear expNames
% Pull fly summary sheet information on selected experiment
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;

% Make new analyzed file directory
analysisDir = [baseFolder vidFolder '\analysis\'];
if ~isfolder(analysisDir); mkdir(analysisDir); end

% Load relevant data files (.mat, .csv, .h5)
warning off
% load matlab file for experiment 
expData = load([baseFolder folder '\' expName ' dataMat.mat']);
% load the temperature log for the experiment
tempLog = readmatrix([baseFolder folder '\' expName '_RampLog']);
% Number of videos in the experiment (same for all arenas)
nvids = expData.parameters.numVids;

% Select which Arena to process
% arenaSel = Alphabet(listdlg('ListString', {'A', 'B', 'C', 'D'}, 'SelectionMode', 'Single'));
expPDF = [analysisDir folder expName arenaSel ' summary.pdf'];
XLrow = find(strcmpi(excelfile(:,Excel.date), folder) & ...
             strcmpi(excelfile(:,Excel.expID), expName) & ...
             strcmpi(excelfile(:,Excel.arena), arenaSel(end:end))); %rows with sel date with ARENA ID

% load tracking predictions         
for vid = 1:nvids
    filePath = [baseFolder vidFolder '\' expName '_' num2str(vid) '.h5'];
    data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
    data(vid).tracks = h5read(filePath,'/tracks');
end; clear filePath

initial_vars = who; initial_vars{end+1} = 'initial_vars';
fprintf('\nNext\n')

% --------------------- Parameter extraction ----------------------------------------
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

%----------------------------- Data organization by video ---------------------------
% find number of 
currFrame = 0;
for vid = 1:nvids
    occupancy.frameROI(vid,1) = currFrame+1;
    nframes = size(data(vid).occupancy_matrix,2);
    currFrame = nframes + currFrame;
    occupancy.frameROI(vid,2) = currFrame;
end

% --------------------------- Mask the food wells & arena ---------------------------
maskPath = [analysisDir expName arenaSel ' ArenaMask.mat'];
if ~isfile(maskPath)
    idx = find(~strcmpi('Empty', wellLabels));
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
else
    load(maskPath)
end

% Screen any points outside the arena
r = 400; %pixel arena radius
x1 = wellcenters(1,1:2:4);
y1 = wellcenters(2,1:2:4);
x2 = wellcenters(1,2:2:4);
y2 = wellcenters(2,2:2:4);
[xi,yi] = polyxpoly(x1,y1,x2,y2);
centre = [xi;yi];

% Tracking matrix locations: [frame, node, xy, fly]
scaleFactor = 1.31; % WHY IS THIS NECESSARY??? TODO

% MASK FOOD AND ARENA:
[raw_flyCount, mid_flyCount] = deal([]);
for vid = 1:nvids
    % number of flies labeled on each frame:
    raw_flyCount = [raw_flyCount; sum(data(vid).occupancy_matrix)'];
    
    % all data for head tracked location
    headData = squeeze(data(vid).tracks(:,1,:,:));
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
    for ii = 1:length(mask)
        [in,on] = inpolygon(X,Y, mask(ii).n(:,1),mask(ii).n(:,2));   % Logical Matrix
        inon = in | on;                                    % Combine ‘in’ And ‘on’
        X(inon) = NaN;
        Y(inon) = NaN;
    end
    % Remove any labeled points outside the arena:
    loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)>=r;
    X(loc) = NaN; Y(loc) = NaN;

    % Resize the data to OG structure:
    X = reshape(X,Xdim);
    Y = reshape(Y,Ydim);

    data(vid).x_loc_mid = X; % save for later convience
    data(vid).y_loc_mid = Y;
    
    % New fly count measure:
    mid_flyCount = [mid_flyCount; sum(~isnan(X),2)];
end

% --------------------------- Mask the overtracked points ---------------------------
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


% ------------------------------ Visualization --------------------------------------
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
    
fig = figure; set(fig, 'pos', [37 518 1171 826]);
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

% FLY COUNT HISTOGRAM
subplot(nrows,ncols,sb(3).idx)
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
%    Y(isnan(Y)) = []; 
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


initial_vars = [initial_vars; 'radii'; 'welldata'; 'flyCount'; 'occupancy'];
clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

%% Automated here on down to just accept the pics and move on...
%% --------Find well count + distance from each well + temperature-------------------
% Find number of flies that are near each well for each frame
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
end

clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

% ------------------------------ Summary figure -------------------------------------
nrow = 5; ncol = 4;
subplotInd(1).idx = 5:7; % temperature
subplotInd(2).idx = [9:11,13:15,17:19]; % occupation
subplotInd(3).idx = 1:3; % fly count
subplotInd(4).idx = 4:4:20; % histogram

% group data across videos:
[plotX, plotY] = deal([]);
for vid = 1:nvids
    plotX = [plotX, data(vid).tempLog]; % temperature
    plotY = [plotY; data(vid).well_counts]; % flies per well
end
plotZ = occupancy.flycount;
plotY = plotY./nflies;
sSpan = 180;
LW = 2;
time = occupancy.time;

fig = getfig(''); 
 % tracking accuracy proxy (# flies)
 subplot(nrow,ncol,subplotInd(3).idx)
    y = smooth(plotZ,sSpan);
    roi = 2:length(y)-1;
    plot(time(roi), y(roi), 'linewidth', LW, 'color', Color('grey'))
    hline(nflies, 'w--')
    ylabel('fly count')
 
 % temperature over time
 subplot(nrow,ncol,subplotInd(1).idx)
    y = smooth(plotX,sSpan);
    roi = 2:length(y)-1;
    plot(time(roi), y(roi), 'linewidth', LW, 'color', 'w')
    ylabel('temp (\circC)')
    ylim([5,26])
    
 % occupation probability
 subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        y = smooth(plotY(:,well),sSpan);
        roi = 2:length(y)-1;
        plot(time(roi), y(roi), 'linewidth', LW, 'color', kolor);
    end
    xlabel('time (min)'); ylabel('occupation probability')
    
 % fly count histogram
 subplot(nrow,ncol,subplotInd(4).idx)
    yyaxis left
    set(gca,'YColor', 'k')
    yyaxis right
    h = histogram(flyCount); vline(nflies, 'w--'); 
    set(h, 'facecolor', Color('grey'))
    xlabel('Number tracked flies'); ylabel('Frame count')
%     labelHandles = findall(gca, 'type', 'text', 'handlevisibility', 'off');
%     set(labelHandles,'FontSize', 18);
%     set(gca,'fontsize',15,'FontWeight','normal');
    
formatFig(fig, true, [nrow, ncol], subplotInd);
subplot(nrow,ncol,subplotInd(2).idx)
l = legend(strrep(wellLabels,'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','position', [0.6463 0.5633 0.0963 0.1126]);% [0.8780 0.8119 0.0963 0.1126])%
subplot(nrow,ncol,subplotInd(1).idx)
set(gca, 'XColor', 'k')
subplot(nrow,ncol,subplotInd(3).idx)
set(gca, 'XColor', 'k')
titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
title(titleName,'color', 'w')

% Save image
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
    end
end  
save_figure(fig, [analysisDir expName arenaSel ' summary figure'], '-png', autoSave);


clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

% -------Simple visualization: relationship between temp & well occupation-----------
if essentialfigs == false
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
    plotY = plotY./nflies;
    sSpan = 180;
    LW = 1;
    time = occupancy.time;

    fig = getfig(''); 
        subplot(nrow,ncol,subplotInd(1).idx)
        y = smooth(plotX,sSpan);
        plot(time(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'w')
        ylabel('temp (\circC)')
        ylim([5,27])

        subplot(nrow,ncol,subplotInd(2).idx)
        hold on
        % error fills
        for well = 1:4
            kolor = pullFoodColor(wellLabels{well}); % plotting color for food
            y_avg(:,well) = smooth(plotY(:,well),sSpan);
            y_err = movstd(plotY(:,well),sSpan);
            fill_data = error_fill(time, y_avg(:,well), y_err);
            h(well) = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h(well), 'facealpha', 0.2)
        end
        % average line
        for well = 1:4
            kolor = pullFoodColor(wellLabels{well});
            plot(time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
        end
        xlabel('time (min)'); ylabel('occupation probability')

    formatFig(fig, true, [nrow, ncol], subplotInd);
    l = legend([{'';'';'';''};strrep(wellLabels,'_','-')]);
    set(l, 'color', 'k', 'textcolor', 'w','position', [0.7947 0.6462 0.0963 0.1126])
    for well = 1:4
        set(get(get(h(well),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    subplot(nrow,ncol,subplotInd(1).idx)
    set(gca, 'XColor', 'k')
    titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
    title(titleName,'color', 'w')


    % Save image
    save_figure(fig, [analysisDir expName arenaSel ' well occupation timcourse'], '-png', autoSave);

    initial_vars = [initial_vars; 'time'];
    clearvars('-except',initial_vars{:})
    fprintf('\nNext\n')
end

% ----------------------Organize info from experiment--------------------------------
% group data across videos:
[plotX, plotY] = deal([]);
for vid = 1:nvids
    plotX = [plotX, data(vid).tempLog];
    plotY = [plotY; data(vid).well_counts];
end
occupancy.temp = plotX;
occupancy.occ = plotY./nflies;
occupancy.count = plotY;
% add total well occupancy to occupancy structure
occupancy.allwellOcc = sum(occupancy.occ,2);


% Measure of eccentricity
%how far are the flies from the center of the arena, on average?
% auto generate the arena center from the middlepoint of the 4 wells?
x1 = wellcenters(1,1:2:4);
y1 = wellcenters(2,1:2:4);
x2 = wellcenters(1,2:2:4);
y2 = wellcenters(2,2:2:4);
[xi,yi] = polyxpoly(x1,y1,x2,y2);
wellcenters(:,5) = [xi;yi];
% 
% %visualize the arena center
% fig = getfig; set(fig, 'color', 'k')
% imshow(demoImg); axis tight square
% hold on
% scatter(xi,yi, 45, 'y', 'filled')

% find distance from center for each fly:
EDist = [];
for vid = 1:nvids
    N = [];
    for frame = 1:length(data(vid).tempLog)
        Y = data(vid).y_loc(frame,:);
        X = data(vid).x_loc(frame,:);
        centre = wellcenters(:,5);
        D = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5); %distance from center of arena
        D(isnan(D)) = [];
        N(frame,:) = [mean(D), std(D)];
    end
    data(vid).eccentricity = N;
    EDist = [EDist;N];
end
occupancy.eccentricity = EDist;

% Plot the results of the eccentricity ** could be cool to do this overlaid
% onto the arena image like a polar plot but with temp being color
% coordinated ** 
% 
% % % Try polar plot vectorization:
% function [r,theta]=cart2polar(x,y)
%     r=sqrt(x.^2+y.^2);
%     theta=atan(y./x);
% end

clearvars('-except',initial_vars{:})
fprintf('Next\n')

% -----------------------Calculate degree of clustering------------------------------
% pull info for each video
LDist = []; LD_err = [];
for vid = 1:nvids
    flyDistance = []; D_err = [];
    for frame = 1:length(data(vid).tempLog)
        test = [data(vid).x_loc(frame,:)',data(vid).y_loc(frame,:)'];
        D = pdist(test);
        D(isnan(D)) = [];
        flyDistance(frame) = mean(D);
        D_err(frame) = std(D);
    end
    data(vid).flyDistance = flyDistance;
    LDist = [LDist,flyDistance]; 
    LD_err = [LD_err,D_err];
end

%add inter-fly-distance to the occupancy time course structure
occupancy.IFD = LDist;
occupancy.IFD_err = LD_err;

% demo the clustering proxy
if nvids>8
    divisor = round(nvids/6);
    vidList = 1:divisor:nvids;
else
    vidList = 1:nvids;
end
nrow = 2; ncol = length(vidList); ii = 0; 
[~,minidx] = deal([]);

% VISUALIZE a demo of the clustering accuracy
if essentialfigs == false
    fig = getfig(''); set(fig, 'pos',[120 331 1244 368], 'color', 'k');
    for vid = vidList
        ii = ii+1;
        movieInfo = VideoReader([baseFolder vidFolder '\' expName '_' num2str(vid) '.avi']); %read in video
        headData = squeeze(data(vid).tracks(:,1,:,:));

        % Most clustered:
        [M,minidx(ii)] = min(data(vid).flyDistance);
        img = read(movieInfo,minidx(ii));
        subplot(nrow, ncol, ii)
        imshow(img); hold on
        axis tight square
        set(gca, 'visible', 'off')
        x = data(vid).x_loc(minidx(ii),:); x(isnan(x)) = [];
        y = data(vid).y_loc(minidx(ii),:); y(isnan(y)) = [];
    %     x = squeeze(headData(minidx(ii), 1, :)); x(isnan(x)) = [];
    %     y = squeeze(headData(minidx(ii), 2, :)); y(isnan(y)) = [];
        scatter(x,y, 10, 'y', 'filled')
        title(num2str(M))
        % overlay 'size bar' for min dist:
        plot([100,100+M], [20,20], 'linewidth', 0.5, 'color', 'r')

        % Least clustered:
        [M,maxidx(ii)] = max(data(vid).flyDistance);
        img = read(movieInfo,maxidx(ii));
        subplot(nrow, ncol, ii+length(vidList))
        imshow(img); hold on
        axis tight square
        set(gca, 'visible', 'off')
        x = data(vid).x_loc(maxidx(ii),:); x(isnan(x)) = [];
        y = data(vid).y_loc(maxidx(ii),:); y(isnan(y)) = [];
    %     x = squeeze(headData(maxidx(ii), 1, :)); x(isnan(x)) = [];
    %     y = squeeze(headData(maxidx(ii), 2, :)); y(isnan(y)) = [];
        scatter(x,y, 10, 'y', 'filled')
        title(num2str(M))
        % overlay 'size bar' for min dist:
        plot([100,100+M], [20,20], 'linewidth', 0.5, 'color', 'r')
    end
    labelHandles = findall(gcf, 'type', 'text', 'handlevisibility', 'off');
    set(labelHandles,'FontSize', 15, 'color', 'w')

    save_figure(fig, [analysisDir expName arenaSel ' linear clustering demo'], '-pdf', autoSave);
end

clearvars('-except',initial_vars{:})
fprintf('Next\n')

% ------------------------------Movement analysis------------------------------------
%bin frame and check bin occupation changes across frames as proxy for
%movement --> won't give an accurate 'speed' etc. but it might give
%movement.
tic
nbins = 100;
spaceChange = [];

for vid = 1:nvids   
    N = [];
    for frame = 1:length(data(vid).tempLog)
        X = data(vid).x_loc(frame,:); X(isnan(X)) = [];
        Y = data(vid).y_loc(frame,:); Y(isnan(Y)) = [];
        N(:,:,frame) = histcounts2(X,Y,nbins);
    end
    if vid>1 % add last frame of previous vid to prevent frame loss 
        N = cat(3, lastFrame, N);
    end
    lastFrame =  N(:,:,end);
    binDiff = diff(N,1,3);
    binDiff(binDiff<0) = 0;
    BinChange = squeeze(sum(sum(binDiff,1),2));
    % save pertinent information:
    data(vid).BinChange = BinChange;
    spaceChange = [spaceChange; BinChange];
end
occupancy.movement = spaceChange;
toc

% preview the movement average
if essentialfigs==false
    fig = getfig;  hold on
    plot(occupancy.time(2:end), occupancy.movement,'color', Color('teal'))
    plot(occupancy.time(2:end), smooth(occupancy.movement,180),...
        'linewidth', 2, 'color', 'w')
    xlabel('time (min)'), ylabel('movement (a.u.)')
    formatFig(fig, true);
    save_figure(fig, [analysisDir expName arenaSel ' movement over time'], '-pdf', autoSave);
end

clearvars('-except',initial_vars{:})
fprintf('Next\n')

% -----------------------Time course feature comparison------------------------------
% plot parameters:
nrow = 8; ncol = 1;
sSpan = 180; %1 minute filter length
sbpts(1).idx = 1;
sbpts(2).idx = 2:4;
sbpts(3).idx = 5:6;
sbpts(4).idx = 7:8;
LW = 1.5;

% FIGURE:
fig = getfig; set(fig, 'pos', [157 86 1232 878])
% TEMPERATURE
subplot(nrow,ncol,sbpts(1).idx)
plot(occupancy.time,smooth(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
ylabel('(\circ)')
title('temperature')
ax = gca;
x_lim = ax.XLim;

% OCCUPANCY
subplot(nrow,ncol,sbpts(2).idx)
hold on
y = [];
for well = 1:4
   y = [y, smooth(occupancy.occ(:,well),sSpan)];
end
h = area(occupancy.time,y);
for well = 1:4
    h(well).FaceColor = pullFoodColor(wellLabels{well});
end
xlim(x_lim)
set(gca, 'tickdir', 'out')
l = legend(strrep(wellLabels,'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k',...
'position', [0.7972 0.7194 0.1039 0.0792]);% [0.8780 0.8119 0.0963 0.1126])%
ylabel('occupancy probability')
title('individual well occupancy')

% 
% hold on
% for well = 1:4
%     kolor = pullFoodColor(wellLabels{well});
%     plot(time,smooth(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
% end
% plot(time,smooth(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
% ylabel('occupancy probability')
% title('individual well occupancy')
% legend(wellLabels)
% l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
% set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
%     'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%

% CLUSTERING
subplot(nrow,ncol,sbpts(3).idx); hold on
y_avg = smooth(occupancy.IFD,sSpan);
y_err = smooth(occupancy.IFD_err,sSpan);
x = occupancy.time;
fill_data = error_fill(x, y_avg, y_err);
h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(occupancy.time,smooth(occupancy.IFD,sSpan),'linewidth', LW, 'color', 'w')
ylabel('pixels')
title('inter-fly-distance')

% MOVEMENT
subplot(nrow,ncol,sbpts(4).idx); hold on
y_avg = smooth(occupancy.movement,sSpan);
y_err = movstd(occupancy.movement,sSpan);
x = occupancy.time(2:end);
fill_data = error_fill(x, y_avg, y_err);
h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(x,y_avg, 'linewidth', LW, 'color', 'w')  
ylabel('(a.u.)')
title('movement')
xlabel('time (min)')

formatFig(fig,true,[nrow, ncol], sbpts);
for ii = 1:3
    subplot(nrow,ncol,sbpts(ii).idx)
    set(gca, 'XColor', 'k')
end

% save and export figure
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
    end
end  
save_figure(fig, [analysisDir expName arenaSel ' timecourse features'], '-png',autoSave);

clearvars('-except',initial_vars{:})
fprintf('Next\n')

% ---------------------Position & Movement Summary Figure----------------------------
% **TODO remove total well occupancy
% figure parameters:
nrow = 10; ncol = 1;
sSpan = 180; %1 minute filter length
sbpts(1).idx = 1;    % temp
sbpts(2).idx = 2:5;  % occupancy
sbpts(3).idx = 6:8;  % clustering & eccentricity
sbpts(4).idx = 9:10; % movement
LW = 1.5;

% FIGURE:
fig = getfig; set(fig, 'pos', [157 86 1232 878])
subplot(nrow,ncol,sbpts(1).idx)
plot(occupancy.time,smooth(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
ylabel('Temp(\circ)')
% title('temperature')

subplot(nrow,ncol,sbpts(2).idx)
hold on
for well = 1:4
    kolor = pullFoodColor(wellLabels{well});
    plot(occupancy.time,smooth(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
end
% plot(occupancy.time,smooth(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
ylabel('occupancy probability')
% title('individual well occupancy')
legend(wellLabels)
% l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
l = legend(strrep(wellLabels,'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
    'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%

% CLUSTERING
subplot(nrow,ncol,sbpts(3).idx); hold on
y_avg = smooth(occupancy.IFD,sSpan);
y_err = smooth(occupancy.IFD_err,sSpan);
x = occupancy.time;
fill_data = error_fill(x, y_avg, y_err);
h(1) = fill(fill_data.X, fill_data.Y, get_color('teal'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(occupancy.time,smooth(occupancy.IFD,sSpan),'linewidth', LW, 'color', Color('teal'))
% ECCENTRICITY
y_avg = smooth(occupancy.eccentricity(:,1),sSpan);
y_err = smooth(occupancy.eccentricity(:,2),sSpan);
x = occupancy.time;
fill_data = error_fill(x, y_avg, y_err);
h(2) = fill(fill_data.X, fill_data.Y, get_color('orange'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(x,y_avg, 'linewidth', LW, 'color', Color('orange'))  
ylabel('pixels')
l2 = legend({'','inter fly distance', '', 'eccentricity'});
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(l2, 'color', 'k', 'textcolor', 'w','FontSize', 8,'edgecolor', 'k',...
    'position', [0.1350 0.4638 0.1039 0.0370]);% 

% MOVEMENT
subplot(nrow,ncol,sbpts(4).idx); hold on
y_avg = smooth(occupancy.movement,sSpan);
y_err = movstd(occupancy.movement,sSpan);
x = occupancy.time(2:end);
fill_data = error_fill(x, y_avg, y_err);
h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(x,y_avg, 'linewidth', LW, 'color', 'w')  
ylabel('movement (a.u.)')
% title('movement')
xlabel('time (min)')

formatFig(fig,true,[nrow, ncol], sbpts);
for ii = 1:3
    subplot(nrow,ncol,sbpts(ii).idx)
    set(gca, 'XColor', 'k')
end

% save and export figure
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
    end
end
save_figure(fig, [analysisDir expName arenaSel ' timecourse features'], '-png',autoSave);

clearvars('-except',initial_vars{:})
fprintf('Next\n')

% ----------------Average distance between each fly and each food source-------------
%(takes a full 2 mins to run)

% find distance from center for each fly:
tic
FDist = [];
idx = 0;
for vid = 1:nvids
    for frame = 1:length(data(vid).tempLog)
        idx = idx+1;
        for well = 1:4
            test = [wellcenters(:,well)'; data(vid).x_loc(frame,:)',data(vid).y_loc(frame,:)'];
            D = squareform(pdist(test));
            D = D(2:end,1);
            D(isnan(D)) = [];
            FDist(well).N(idx,:) = [median(D), std(D)];
        end
    end
end
toc

occupancy.dist2wells = FDist;

% occupancy.eccentricity = EDist;
nrow = 4;
ncol = 1;
subplotInd(1).idx = 1;
subplotInd(2).idx = 2:4;
% group data across videos:
plotX = occupancy.time;
sSpan = 180;
LW = 1;

fig = getfig(''); 
    subplot(nrow,ncol,subplotInd(1).idx)
    y = smooth(occupancy.temp,sSpan);
    plot(plotX(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'w')
    ylabel('temp (\circC)')
    ylim([5,27])
    
    subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    % error fills
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well}); % plotting color for food
        y_avg(:,well) = smooth(FDist(well).N(:,1),sSpan);
        y_err(:,well) = smooth(FDist(well).N(:,2),sSpan);
        fill_data = error_fill(plotX, y_avg(:,well), y_err(:,well));
        h(well) = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
        set(h(well), 'facealpha', 0.2)
    end
    % average line
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        plot(occupancy.time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
    end
    xlabel('time (min)'); ylabel('avg distance between fly and food (a.u.)')
    
formatFig(fig, true, [nrow, ncol], subplotInd);
l = legend([{'';'';'';''};strrep(wellLabels,'_','-')]);
set(l, 'color', 'k', 'textcolor', 'w','position', [0.7947 0.6462 0.0963 0.1126])
for well = 1:4
    set(get(get(h(well),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
subplot(nrow,ncol,subplotInd(1).idx)
set(gca, 'XColor', 'k')
titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
title(titleName,'color', 'w')
 
% save and export figure
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
    end
end
save_figure(fig, [analysisDir expName arenaSel ' avg distance to well'], '-png', autoSave);

clearvars('-except',initial_vars{:})
fprintf('Next\n')

% ----------------------- Save experiment data thus far -----------------------------
if autoSave == true
    clearvars('-except',initial_vars{:})
    initial_vars = unique(initial_vars);
    save([analysisDir expName arenaSel ' timecourse data'])
    copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')
    fprintf('Experiment data saved\n')
else
    switch questdlg('Save experiment analysis?')
        case 'Yes'
            clearvars('-except',initial_vars{:})
            initial_vars = unique(initial_vars);
            save([analysisDir expName arenaSel ' timecourse data'])
            copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')
            fprintf('Experiment data saved\n')

        case 'No'
            return
        case 'Cancel'
            return
    end
end


%% Open former data:


%% Visual comparison of tracked frames to look for outliers / patterns
% 
% switch questdlg('Make tracking example videos?','', 'Yes', 'No', 'No')
%     case 'No'
%         return
%     case 'Cancel'
%         return
%     case 'Yes'  
%         switch questdlg('Demo first 300 frames?')
%             case 'Yes'
%                 nframes = 300;
%             case 'No'
%                 nframes = inf;
%             case 'Cancel'
%                 return
%         end
% end
%     
% if nvids>5
%     divisor = round(nvids/10);
%     vidList = 1:divisor:nvids;
% else
%     vidList = 1:nvids;
% end
% ii = 0; 
% 
% 
% vidpath = [baseFolder folder '\labeled vid examples\'];
% if ~isfolder(vidpath)
%     mkdir(vidpath)
% end
% % find frames with largest fly count offset
% for vid = vidList
%     headdata = squeeze(data(vid).tracks(:,1,:,:));
% 
%     tempVidName = [baseFolder folder '\' expName '_' num2str(vid) '.avi'];
%     movieInfo = VideoReader(tempVidName);
% 
%     fig = getfig('',1); set(fig, 'color', 'k','pos',[-851 329 707 656]);
%     set(gca, 'visible', 'off')
% 
%     vid_name = [vidpath expName '_' num2str(vid) ' predicted frames'];
% 
%     if isinf(nframes)
%         n = movieInfo.NumFrames;
%         FrameRate = 30; %playback rate
%         vid_name = [vid_name ' all frames'];
%     else
%         n = nframes;
%         FrameRate = 6; %playback rate
%     end
% 
%     % initiate video writer
%     v = VideoWriter(vid_name, 'Motion JPEG AVI');
%     v.Quality = 75;
%     v.FrameRate = FrameRate;
%     open(v);
% 
% 
%     for frame = 1:n
%         img = read(movieInfo,frame);
%         imagesc(img)
%         set(fig, 'Color', 'k') 
%         hold on
%         x = squeeze(headdata(frame, 1, :));
%         y = squeeze(headdata(frame, 2, :));
%         x(isnan(x)) = []; % remove empty tracks
%         y(isnan(y)) = [];
%         scatter(x,y, 30, 'y')
%         axis tight; axis square
%         set(gca, 'visible', 'off')
%         f = getframe(fig);
%         writeVideo(v, f)  
%         clf('reset') 
%     end
%     close(v)
%     close(fig)
% end

%% Demo random selection of frames and their tracking points (QC)
% vid = 6;
% headdata = squeeze(data(vid).tracks(:,1,:,:));
% 
% tempVidName = [baseFolder vidFolder '\' expName '_' num2str(vid) '.avi'];
% movieInfo = VideoReader(tempVidName);
% 
% fig = getfig('',1); set(fig, 'color', 'k','visible', 'on');
% set(gca, 'visible', 'off')
% 
% % vid_name = [baseFolder vidFolder '\' expName arenaSel '_' num2str(vid) ' predicted frames'];
% n = floor(movieInfo.NumFrames);
% 
% for frame = 1:15:200
%     img = rgb2gray(read(movieInfo,frame));
%     imshow(img)
%     hold on
%     x = squeeze(headdata(frame, 1, :));
%     y = squeeze(headdata(frame, 2, :));
%     x(isnan(x)) = []; % remove empty tracks
%     y(isnan(y)) = [];
% %     scatter(x,y, 30, 'y')
%     scatter(x.*1.31, y.*1.31, 30, 'r')
%     pause(0.1)
% end

%%
% 
% [offset, indx] = sort(plotZ-nflies);   
% flyCounts = sum(data(vid).occupancy_matrix); %num tracked flies
% [offset, indx] = sort(flyCounts-nflies); %find max incorrect number
% % load image
% flyCounts(indx(end))

    
    
    % add scale bar
    
    
    
%     
% 
% figure;
% plot(indx)
% 
% 
% % Basic sorting of the data...
% % demo only -- load in the camera image for the frame and overlay...
% demovid = 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\NonlinearCooling_1.avi';
% movieInfo = VideoReader(demovid);
% % mat locations: 
% ML.frame = 1;   % frames are 1st dim
% ML.node = 2;    % 
% ML.xy = 3;      % xy coordinates
% ML.fly = 4;     % tracks for 'a' fly
% 
% % all data for head tracked location
% headdata = squeeze(tracks_matrix(:,1,:,:));
% 
% frame = 5;
% 
% % tracked points
% x = squeeze(headdata(frame, 1, :));
% y = squeeze(headdata(frame, 2, :));
% x(isnan(x)) = []; % remove empty tracks
% y(isnan(y)) = [];
% % image
% img = read(movieInfo,frame);
% 
% % overlay:
% fig = getfig; set(fig, 'color', 'k')
% hold on
% imagesc(img)
% scatter(x,y, 30, 'y')
% axis tight; axis square
% set(gca, 'visible', 'off')

%%  
% 
% 
% %% vvvv OLD STUFF BELOW vvvv
% %% 
% 
% fig = getfig('',1);
% hold on
% 
% scatter(plotX, plotY(:,1))
% 
% 
% 
% %% Basic sorting of the data...
% % demo only -- load in the camera image for the frame and overlay...
% demovid = 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\NonlinearCooling_1.avi';
% movieInfo = VideoReader(demovid);
% % mat locations: 
% ML.frame = 1;   % frames are 1st dim
% ML.node = 2;    % 
% ML.xy = 3;      % xy coordinates
% ML.fly = 4;     % tracks for 'a' fly
% 
% % all data for head tracked location
% headdata = squeeze(tracks_matrix(:,1,:,:));
% 
% frame = 5;
% 
% % tracked points
% x = squeeze(headdata(frame, 1, :));
% y = squeeze(headdata(frame, 2, :));
% x(isnan(x)) = []; % remove empty tracks
% y(isnan(y)) = [];
% % image
% img = read(movieInfo,frame);
% 
% % overlay:
% fig = getfig; set(fig, 'color', 'k')
% hold on
% imagesc(img)
% scatter(x,y, 30, 'y')
% axis tight; axis square
% set(gca, 'visible', 'off')
% 
% % save_figure(fig, 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\demo tracking image','-png');
% 
% 
% %% 
% pts = readPoints(img,4); % get locations of the wells from the image
% 
% x = squeeze(headdata(:, 1, :)); % frame, instance
% y = squeeze(headdata(:, 2, :));
% 
% X = reshape(x,numel(x),1);
% Y = reshape(y,numel(y),1);
% loc = isnan(X);
% X(loc) = [];
% Y(loc) = [];
% % squash all points together over the video:
% 
% 
% % what if we do a min radius from the coordinates of the food bowl?  
% fig = getfig; set(fig, 'color', 'k');
% hist2d(X,Y, 'probability', 'tile')
% axis tight; axis square
% set(gca, 'visible', 'off')
% c = colorbar;
% c.Color = [1,1,1];
% hold on
% scatter(pts(1,:),pts(2,:), 75, 'r', 'filled') % could automate the color on this
% 
% % save_figure(fig, 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\demo histogram vid 1','-png');
% 
% 
% %% Find flies within ROI of each well
% 
% % how many points within a specific radius of the well?
% radii = 200; %distance must be less than this number to count for a well ROI
% nInst = length(X);
% % loop for all wells
% for well = 1:4
%     b = pts(:,well); % well 1 points
%     euDist = sqrt((b(1)-X).^2 + (b(2)-Y).^2);
%     data(well).well_loc = (euDist<=radii);
%     data(well).well_count = sum((euDist<=radii));
%     
%     % occupation probability
%     data(well).occ_prob = data(well).well_count/nInst;
%     data(well).nInst = nInst;
% end
% 
% 
% 
% 
% 
% 
% 
% % % visual confirmation that the selected points are near the well:
% % fig = getfig; set(fig, 'color', 'k');
% % hist2d(X,Y, 'probability', 'tile')
% % axis tight; axis square
% % set(gca, 'visible', 'off')
% % c = colorbar;
% % c.Color = [1,1,1];
% % hold on
% % scatter(X(well_loc),Y(well_loc), 50, 'y') % could automate the color on this
% % scatter(pts(1,:),pts(2,:), 75, 'r', 'filled') % could automate the color on this
% % viscircles(b',radii)
% % 
% % save_figure(fig, 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\demo histogram point sel','-png');
% % 
% 
% %% 
% 
% 






   








