

%% Select Date & Experiment to Process
clear; close all; clc
%get base folder pathway
[baseFolder, folder] = getCloudPath(2); 

% Select the complete experiments to process
list_dirs = dir([baseFolder folder, '\*.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
clear expNames
% Pull fly summary sheet information on selected experiment
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
XLrow = find(strcmpi(excelfile(:,Excel.date), folder) & ...
           strcmpi(excelfile(:,Excel.expID), expName)); %rows with sel date

% Make new analyzed file directory
analysisDir = [baseFolder folder '\analysis\'];
if ~isfolder(analysisDir); mkdir(analysisDir); end
expPDF = [analysisDir expName ' summary.pdf'];

% Load relevant data files (.mat, .csv, .h5)
warning off
% load matlab file for experiment 
expData = load([baseFolder folder '\' expName 'dataMat.mat']);

% load the temperature log for the experiment
tempLog = readmatrix([baseFolder folder '\' expName '_RampLog']);

nvids = expData.num.vids;
for vid = 1:nvids
    % load tracking predictions
    filePath = [baseFolder folder '\' expName '_' num2str(vid) '.h5'];
    data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
    data(vid).tracks = h5read(filePath,'/tracks');
end

initial_vars = who; initial_vars{end+1} = 'initial_vars';
fprintf('\nNext\n')

%% Parameter extraction
% Number of flies:
nflies = excelfile{XLrow,Excel.numflies};
movieInfo = VideoReader([baseFolder folder '\' expName '_1.avi']); %read in video
ii = randi(size(data(1).occupancy_matrix,2),[3,1]); %random selection of frames to count fly number
demoImg = rgb2gray(read(movieInfo,1));
if isnan(nflies)
    % manual count of flies
    fprintf('\nCount the number of flies in the picture by clicking them\n then hit ENTER\n')
    for jj = 1:3
        demoImg = rgb2gray(read(movieInfo,ii(jj)));
        nflies(jj) = size(readPoints(demoImg),2);
    end
    if sum(~diff(nflies)==0)==0
        nflies = nflies(1);
    else
        warndlg('Check fly counts again')
        fprintf(['\nNumber of flies counted: ' num2str(nflies) '\n'])
        return
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
h = warndlg('Select the well identities in their printed order'); uiwait(h)
% label the wells:
wellLabels = {expData.params.well_1; expData.params.well_2;...
          expData.params.well_3; expData.params.well_4};   
disp(wellLabels)
wellcenters = readPoints(demoImg,4); % get locations of the wells from the image

initial_vars = [initial_vars; 'nflies'; 'wellLabels'; 'wellcenters'; 'demoImg'];
clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Data organization by video
% Tracking matrix locations: [frame, node, xy, fly]

% pull info for each video
flyCount = []; 
for vid = 1:nvids
    % number of flies labeled on each frame:
    flyCount = [flyCount; sum(data(vid).occupancy_matrix)'];
    
    % all data for head tracked location
    headData = squeeze(data(vid).tracks(:,1,:,:));
    nframes = size(headData,1); 
    
    % x-y coordinates of flies for each frame
    x_loc = squeeze(headData(:,1,:));
    y_loc = squeeze(headData(:,2,:));
    data(vid).x_loc = x_loc; % save for later convience
    data(vid).y_loc = y_loc;
    
    % temperature alignment
    
    logROI(1) = find(tempLog(:,1)==expData.tempLogStart(vid,3));
    logROI(2) = find(tempLog(:,1)==expData.tempLogEnd(vid,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    x = round(linspace(1, nframes, length(tempCourse)));
    % upsample the temperature log:
    fullTempList = interp1(x,tempCourse,1:nframes,'spline');   

    % Find number of flies that are near each well for each frame
    radii = 200; %distance must be less than this number to count for a well ROI
    
    % loop for all wells
    for well = 1:4
        b = wellcenters(:,well); % well 1 points
        % reshape data points from whole video for optimized maths
        ntracks = size(x_loc,2);
        X = reshape(x_loc,numel(x_loc),1); 
        Y = reshape(y_loc,numel(y_loc),1); 
        % within well range?
        euDist = sqrt((b(1)-X).^2 + (b(2)-Y).^2); %euclidian distance from well center
        well_loc = reshape((euDist<=radii),[nframes,ntracks]);
        welldata(well).loc = well_loc;
        % how many within the circle?
        welldata(well).count = sum(well_loc,2);
        data(vid).well_counts(:,well) = welldata(well).count;
    end
    data(vid).tempLog = fullTempList;
end

% % visual confirmation that the selected points are near the well:
% AllPoints = [];
% for vid = 1:nvids
%    headData = squeeze(data(vid).tracks(:,1,:,:));
%    x_loc = squeeze(headData(:,1,:));
%    y_loc = squeeze(headData(:,2,:));
%    X = reshape(x_loc,numel(x_loc),1); 
%    Y = reshape(y_loc,numel(y_loc),1); 
%    X(isnan(X)) = [];
%    Y(isnan(Y)) = []; 
%    AllPoints = [AllPoints; X , Y];
% end
% 
% fig = getfig; set(fig, 'color', 'k');
% hist2d(AllPoints(:,1),AllPoints(:,2), 'probability', 'tile')
% axis tight; axis square
% set(gca, 'visible', 'off')
% c = colorbar;
% c.Color = [1,1,1];
% hold on
% % scatter(X(well_loc),Y(well_loc), 50, 'y') % could automate the color on this
% % scatter(pts(1,:),pts(2,:), 75, 'r', 'filled') % could automate the color on this
% viscircles(pts',[radii,radii,radii,radii])
% b = pts(:,well); % well 1 points

% visual check that the number of tracked flies is close to the actual number
nrow = 1; ncol = 2;
fig = getfig; set(fig, 'pos', [195 157 1413 646]);%[2160 185 1413 646]) <--evynpc
subplot(nrow, ncol, 2)
    histogram(flyCount); vline(nflies, 'w-'); 
    xlabel('Number tracked flies'); ylabel('Frame count')
subplot(nrow, ncol, 1)
    % Take and save a snapshop pic of the arena w/ labels & well roi
    outerColors = {'red', 'yellow', 'blue', 'green'};
    imshow(demoImg)
    axis tight square
    hold on
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        scatter(wellcenters(1,well),wellcenters(2,well), 75,...
            'MarkerFaceColor', kolor, 'MarkerEdgeColor', Color(outerColors{well})) 
    end
    
    l = legend(strrep(wellLabels,'_','-')); 
    set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.0959 0.7409 0.1271 0.1428]);
    viscircles(wellcenters',ones(1,4)*radii);
    titleName = strrep([folder ' ' expName], '_',' ');
    title(titleName,'color', 'w')
    formatFig(fig, true,[nrow, ncol]);
    
% save and export figure:
if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
end  
save_figure(fig, [analysisDir expName ' quality control'], '-png');

initial_vars = [initial_vars; 'radii'; 'welldata'; 'flyCount'];
clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

%% Summary figure
nrow = 5; ncol = 4;
subplotInd(1).idx = 5:7; % temperature
subplotInd(2).idx = [9:11,13:15,17:19]; % occupation
subplotInd(3).idx = 1:3; % fly count
subplotInd(4).idx = 4:4:20; % histogram

% group data across videos:
[plotX, plotY, plotZ] = deal([]);
for vid = 1:nvids
    plotX = [plotX, data(vid).tempLog]; % temperature
    plotY = [plotY; data(vid).well_counts]; % flies per well
    plotZ = [plotZ, sum(data(vid).occupancy_matrix)]; %number of flies tracked
end
plotY = plotY./nflies;
sSpan = 180;
LW = 1;
time = linspace(1,(length(plotX)/3)/60,length(plotX));

fig = getfig(''); 
 % tracking accuracy proxy (# flies)
 subplot(nrow,ncol,subplotInd(3).idx)
    y = moving_average(plotZ,sSpan);
    roi = 2:length(y)-1;
    plot(time(roi), y(roi), 'linewidth', LW, 'color', Color('grey'))
    hline(nflies, 'w--')
    ylabel('fly count')
 
 % temperature over time
 subplot(nrow,ncol,subplotInd(1).idx)
    y = moving_average(plotX,sSpan);
    roi = 2:length(y)-1;
    plot(time(roi), y(roi), 'linewidth', LW, 'color', 'w')
    ylabel('temp (\circC)')
    ylim([5,26])
    
 % occupation probability
 subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        y = moving_average(plotY(:,well),sSpan);
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
titleName = strrep([folder ' ' expName], '_',' ');
title(titleName,'color', 'w')

% Save image
if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
end  
save_figure(fig, [analysisDir expName ' summary figure'], '-png');

clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

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
plotY = plotY./nflies;
sSpan = 180;
LW = 1;
time = linspace(1,(length(plotX)/3)/60,length(plotX));

fig = getfig(''); 
    subplot(nrow,ncol,subplotInd(1).idx)
    y = moving_average(plotX,sSpan);
    plot(time(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'w')
    ylabel('temp (\circC)')
    ylim([5,27])
    
    subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    % error fills
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well}); % plotting color for food
        y_avg(:,well) = moving_average(plotY(:,well),sSpan);
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
titleName = strrep([folder ' ' expName], '_',' ');
title(titleName,'color', 'w')


% Save image
save_figure(fig, [analysisDir expName ' well occupation timcourse'], '-png');

initial_vars = [initial_vars; 'time'];
clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

%% Organize info from experiment
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

initial_vars = [initial_vars; 'occupancy'];
clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Measure of eccentricity
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
        test = [wellcenters(:,5)'; data(vid).x_loc(frame,:)',data(vid).y_loc(frame,:)'];
        D = squareform(pdist(test));
        D = D(2:end,1);
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

%% Calculate degree of clustering
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
% time = linspace(1,(length(LDist)/3)/60,length(LDist));
occupancy.time = time;
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
[minidx,minidx] = deal([]);

% VISUALIZE a demo of the clustering accuracy
fig = getfig(''); set(fig, 'pos',[120 331 1244 368], 'color', 'k');
for vid = vidList
    ii = ii+1;
    movieInfo = VideoReader([baseFolder folder '\' expName '_' num2str(vid) '.avi']); %read in video
    headData = squeeze(data(vid).tracks(:,1,:,:));
    
    % Most clustered:
    [M,minidx(ii)] = min(data(vid).flyDistance);
    img = read(movieInfo,minidx(ii));
    subplot(nrow, ncol, ii)
    imshow(img); hold on
    axis tight square
    set(gca, 'visible', 'off')
    x = squeeze(headData(minidx(ii), 1, :)); x(isnan(x)) = [];
    y = squeeze(headData(minidx(ii), 2, :)); y(isnan(y)) = [];
    scatter(x,y, 10, 'y', 'filled')
    title(num2str(M))
    % overlay 'size bar' for min dist:
    plot([10,10+M], [20,20], 'linewidth', 0.5, 'color', 'r')
    
    % Least clustered:
    [M,maxidx(ii)] = max(data(vid).flyDistance);
    img = read(movieInfo,maxidx(ii));
    subplot(nrow, ncol, ii+length(vidList))
    imshow(img); hold on
    axis tight square
    set(gca, 'visible', 'off')
    x = squeeze(headData(maxidx(ii), 1, :)); x(isnan(x)) = [];
    y = squeeze(headData(maxidx(ii), 2, :)); y(isnan(y)) = [];
    scatter(x,y, 10, 'y', 'filled')
    title(num2str(M))
    % overlay 'size bar' for min dist:
    plot([10,10+M], [20,20], 'linewidth', 0.5, 'color', 'r')
end
labelHandles = findall(gcf, 'type', 'text', 'handlevisibility', 'off');
set(labelHandles,'FontSize', 15, 'color', 'w')

save_figure(fig, [analysisDir expName ' linear clustering demo'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Movement analysis

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
fig = getfig;  hold on
plot(occupancy.time(2:end), occupancy.movement,'color', Color('teal'))
plot(occupancy.time(2:end), moving_average(occupancy.movement,180),...
    'linewidth', 2, 'color', 'w')
xlabel('time (min)'), ylabel('movement (a.u.)')
formatFig(fig, true);
save_figure(fig, [analysisDir expName ' movement over time'], '-pdf');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Time course feature comparison
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
subplot(nrow,ncol,sbpts(1).idx)
plot(time,moving_average(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
ylabel('(\circ)')
title('temperature')

subplot(nrow,ncol,sbpts(2).idx)
hold on
for well = 1:4
    kolor = pullFoodColor(wellLabels{well});
    plot(time,moving_average(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
end
plot(time,moving_average(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
ylabel('occupancy probability')
title('individual well occupancy')
legend(wellLabels)
l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
    'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%

% CLUSTERING
subplot(nrow,ncol,sbpts(3).idx); hold on
y_avg = moving_average(occupancy.IFD,sSpan);
y_err = moving_average(occupancy.IFD_err,sSpan);
x = time;
fill_data = error_fill(x, y_avg, y_err);
h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(time,moving_average(occupancy.IFD,sSpan),'linewidth', LW, 'color', 'w')
ylabel('pixels')
title('inter-fly-distance')

% MOVEMENT
subplot(nrow,ncol,sbpts(4).idx); hold on
y_avg = moving_average(occupancy.movement,sSpan);
y_err = movstd(occupancy.movement,sSpan);
x = time(2:end);
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
if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
end  
save_figure(fig, [analysisDir expName ' timecourse features'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Position & Movement Summary Figure 
%**TODO add in the eccentricity, remove total well occupancy
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
plot(time,moving_average(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
ylabel('Temp(\circ)')
% title('temperature')

subplot(nrow,ncol,sbpts(2).idx)
hold on
for well = 1:4
    kolor = pullFoodColor(wellLabels{well});
    plot(time,moving_average(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
end
% plot(time,moving_average(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
ylabel('occupancy probability')
% title('individual well occupancy')
legend(wellLabels)
% l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
l = legend(strrep(wellLabels,'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
    'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%

% CLUSTERING
subplot(nrow,ncol,sbpts(3).idx); hold on
y_avg = moving_average(occupancy.IFD,sSpan);
y_err = moving_average(occupancy.IFD_err,sSpan);
x = time;
fill_data = error_fill(x, y_avg, y_err);
h(1) = fill(fill_data.X, fill_data.Y, get_color('teal'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(time,moving_average(occupancy.IFD,sSpan),'linewidth', LW, 'color', Color('teal'))
% ECCENTRICITY
y_avg = moving_average(occupancy.eccentricity(:,1),sSpan);
y_err = moving_average(occupancy.eccentricity(:,2),sSpan);
x = time;
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
y_avg = moving_average(occupancy.movement,sSpan);
y_err = movstd(occupancy.movement,sSpan);
x = time(2:end);
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
if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
end  
save_figure(fig, [analysisDir expName ' timecourse features'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Save experiment data thus far:

switch questdlg('Save experiment analysis?')
    case 'Yes'
        clearvars('-except',initial_vars{:})
        initial_vars = unique(initial_vars);
        save([analysisDir expName ' timecourse data'])
        fprintf('Experiment data saved\n')
    case 'No'
        return
    case 'Cancel'
        return
end



%% Visual comparison of tracked frames to look for outliers / patterns

switch questdlg('Make tracking example videos?')
    case 'No'
        return
    case 'Cancel'
        return
    case 'Yes'  
        switch questdlg('Demo first 300 frames?')
            case 'Yes'
                nframes = 300;
            case 'No'
                nframes = inf;
            case 'Cancel'
                return
        end
end
    
if nvids>5
    divisor = round(nvids/10);
    vidList = 1:divisor:nvids;
else
    vidList = 1:nvids;
end
ii = 0; 


vidpath = [baseFolder folder '\labeled vid examples\'];
if ~isfolder(vidpath)
    mkdir(vidpath)
end
% find frames with largest fly count offset
for vid = vidList
    headdata = squeeze(data(vid).tracks(:,1,:,:));

    tempVidName = [baseFolder folder '\' expName '_' num2str(vid) '.avi'];
    movieInfo = VideoReader(tempVidName);

    fig = getfig('',1); set(fig, 'color', 'k','pos',[-851 329 707 656]);
    set(gca, 'visible', 'off')

    vid_name = [vidpath expName '_' num2str(vid) ' predicted frames'];

    if isinf(nframes)
        n = movieInfo.NumFrames;
        FrameRate = 30; %playback rate
        vid_name = [vid_name ' all frames'];
    else
        n = nframes;
        FrameRate = 6; %playback rate
    end

    % initiate video writer
    v = VideoWriter(vid_name, 'Motion JPEG AVI');
    v.Quality = 75;
    v.FrameRate = FrameRate;
    open(v);


    for frame = 1:n
        img = read(movieInfo,frame);
        imagesc(img)
        set(fig, 'Color', 'k') 
        hold on
        x = squeeze(headdata(frame, 1, :));
        y = squeeze(headdata(frame, 2, :));
        x(isnan(x)) = []; % remove empty tracks
        y(isnan(y)) = [];
        scatter(x,y, 30, 'y')
        axis tight; axis square
        set(gca, 'visible', 'off')
        f = getframe(fig);
        writeVideo(v, f)  
        clf('reset') 
    end
    close(v)
    close(fig)
end

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

%% Demo random selection of frames and their tracking points (QC)

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






   








