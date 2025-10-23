
%% Update name in high res experiments

parameters.expID = 'Berlin_courtship_F_LRR_no_food_ramp1';
parameters.videoName = 'Berlin_courtship_F_LRR_no_food_ramp1';

expName = 'Berlin_courtship_F_LRR_no_food_ramp1';


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

%% change food type
parameters.ArenaA.well_1 = 'Caviar';
parameters.ArenaB.well_1 = 'Caviar';
parameters.ArenaC.well_1 = 'Caviar';
parameters.ArenaD.well_1 = 'Caviar';

%% change genotype
parameters.ArenaA.genotype = 'Empty_dbd-Gal4>UAS-Kir2.1_A1';
parameters.ArenaB.genotype = 'Empty_dbd-Gal4>UAS-Kir2.1_A1';
parameters.ArenaC.genotype = 'Empty_dbd-Gal4>UAS-Kir2.1_A1';
parameters.ArenaD.genotype = 'Empty_dbd-Gal4>UAS-Kir2.1_A1';

%% change sex
parameters.ArenaA.sex = 'Mixed';
parameters.ArenaB.sex = 'Mixed';
parameters.ArenaC.sex = 'Mixed';
parameters.ArenaD.sex = 'Mixed';

%% change temp protocol
parameters.protocol = 'survival_hold_with_recovery_35-25';

%% get base folder pathway
% baseFolder = getCloudPath;
let = {'A' 'B' 'C' 'D'};
for i = 1:4
    for well = 1:4
        parameters.(['Arena' let{i}]).(['well_' num2str(well)]) = 'Empty';
    end
end
clear let i well

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


%%  Compare the temp passive warming ramps

% load the temperature log for the experiment

fileList = {'G:\My Drive\Jeanne Lab\DATA\09.22.2021\PlantYeastChoice_N1_RampLog(1).csv',...
            'G:\My Drive\Jeanne Lab\DATA\09.21.2021\PlantChoice_N1_RampLog.csv',...
            'G:\My Drive\Jeanne Lab\DATA\09.08.2021\PlantYeastChoice_N1_RampLog.csv',...
            'G:\My Drive\Jeanne Lab\DATA\09.08.2021\PlantYeastChoice_N2_RampLog.csv',...
            'G:\My Drive\Jeanne Lab\DATA\09.03.2021\PlantYeastChoice_N1_RampLog.csv'};

nlogs = length(fileList);
for ii = 1:nlogs
    tempLog = readmatrix(fileList{ii});
    temp(ii).log = tempLog(:,2);
end



% 
X = [8 12 16 18 20 22 23 23.7 24 25];
Y = [572 576 582 586 592 602 614 630 650 763];

X = [8 20 22 23.5 24 25]
Y = [191 197 201 209 215 255];

time = round((Y*4)/12);
test = time - (time(1)-15);


figure;
hold on
for ii = 1:nlogs
    y = temp(ii).log;
    plot(y)
end
vline(time)

figure;
hold on
for ii = 1:nlogs
    x = linspace(0,length(temp(ii).log)/12,length(temp(ii).log));
    y = temp(ii).log;
    plot(x,y)
end
plot(Y,X, 'color', 'k', 'linewidth', 1)
xlabel('time (min)')



%% Supervised tracking cleaning:

% tracks(time, body part, x|y, animalID)

size(data(vid).tracks)



% try to find frames that don't move for Threshold time during video
x_loc = squeeze(data(vid).tracks(:,1,1,:));
y_loc = squeeze(data(vid).tracks(:,1,2,:));




%%
newVars = {'dayList', 'root', 'day', 'ii', 'newVars'};
dayList = {'11.11.2021', '11.10.2021', '11.09.2021', '11.08.2021'};
root = 'G:\My Drive\Jeanne Lab\DATA\';

for day = 1:4 %exeriment days
  for ii = 1:4 %wells
    clearvars('-except',newVars{:})

    load([root dayList{day} '/Arena ' Alphabet(ii) '/analysis/PlantYeastChoice' Alphabet(ii) ' timecourse data.mat'])

    % find distance from center for each fly:
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
            plot(time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
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
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');

    export_fig(fig, [analysisDir expName arenaSel ' avg distance to well.png'], '-png', '-nocrop', '-r300' , '-painters', '-rgb');
    close(fig)

      
    dist2wells = FDist;
    save([analysisDir expName arenaSel ' timecourse data.mat'],'dist2wells','-append')

    disp(['Done ' dayList{day} ' arena ' arenaSel])
  end
end

%% Find the pixel to mm conversion:

vidInfo = VideoReader("G:\My Drive\Jeanne Lab\DATA\11.11.2021\Arena A\PlantYeastChoice_1.avi");

demoImg = rgb2gray(read(vidInfo,1));
datapoints = readPoints(demoImg);

% find distance between the two sets of points: 
inner = abs(datapoints(2,1)-datapoints(2,2)); % should = 41.2mm 
outer = abs(datapoints(2,3)-datapoints(2,4)); % should  =  72 mm

pixelperMM_1 = inner/42;
pixelperMM_2 = outer/70;

disp(pixelperMM_1)
disp(pixelperMM_2)

% best of both = the pix to mm conversion

%% try 'cleaning' the videos and see how the tracking compares
% make some sort of analytic comparison of over/under tracking


nflies = excelfile{XLrow,Excel.numflies};
movieInfo = VideoReader([baseFolder vidFolder '\' expName   '_1.avi']); %read in video
demoImg = rgb2gray(read(movieInfo,1));


f = figure;
imshow(demoImg)
title('Outline the food well')
roi = drawpolygon;
mask1 = roi.Position;

% determine if any tracking points fall within the bounded region??

% all data for head tracked location
headData = squeeze(data(1).tracks(:,1,:,:));
nframes = size(headData,1); 

% x-y coordinates of flies for each frame
scaleFactor = 1.31;
x_loc = squeeze(headData(:,1,:)).*scaleFactor;
y_loc = squeeze(headData(:,2,:)).*scaleFactor;
X = reshape(x_loc,[numel(x_loc),1]);
Y = reshape(y_loc,[numel(y_loc),1]);

[in,on] = inpolygon(X,Y, mask1(:,1),mask1(:,2));   % Logical Matrix
inon = in | on;                                    % Combine ‘in’ And ‘on’
X(inon) = NaN;
Y(inon) = NaN;
                                        

%%

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

% Tracking matrix locations: [frame, node, xy, fly]
scaleFactor = 1.31; % WHY IS THIS NECESSARY??? TODO
% pull info for each video
[OG_flyCount, flyCount] = deal([]);
for vid = 1:nvids
    % number of flies labeled on each frame:
    OG_flyCount = [OG_flyCount; sum(data(vid).occupancy_matrix)'];
    
    % all data for head tracked location
    headData = squeeze(data(vid).tracks(:,1,:,:));
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
    for ii = 1:length(idx)
        [in,on] = inpolygon(X,Y, mask(ii).n(:,1),mask(ii).n(:,2));   % Logical Matrix
        inon = in | on;                                    % Combine ‘in’ And ‘on’
        X(inon) = NaN;
        Y(inon) = NaN;
    end

    % Resize the data to OG structure:
    X = reshape(X,Xdim);
    Y = reshape(Y,Ydim);

    data(vid).x_loc = X; % save for later convience
    data(vid).y_loc = Y;
    
    % New fly count measure:
    flyCount = [flyCount; sum(~isnan(X),2)];

    % temperature alignment
    logROI(1) = find(tempLog(:,1)==expData.tempLogStart(vid,3));
    logROI(2) = find(tempLog(:,1)==expData.tempLogEnd(vid,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    x = round(linspace(1, nframes, length(tempCourse)));
    % upsample the temperature log:
    fullTempList = interp1(x,tempCourse,1:nframes,'spline');   

    % Find number of flies that are near each well for each frame
    radii = 165; %110; %distance must be less than this number to count for a well ROI 
    % loop for all wells
    for well = 1:4
        b = wellcenters(:,well); 
        % reshape data points from whole video for optimized maths
        ntracks = size(data(vid).x_loc,2);
        x_loc = reshape(data(vid).x_loc,numel(data(vid).y_loc),1); 
        y_loc = reshape(data(vid).y_loc,numel(data(vid).y_loc),1); 
        % within well range?
        euDist = sqrt((b(1)-x_loc).^2 + (b(2)-y_loc).^2); %euclidian distance from well center
        well_loc = reshape((euDist<=radii),[nframes,ntracks]);
        welldata(well).loc = well_loc;
        % how many within the circle?
        welldata(well).count = sum(well_loc,2);
        data(vid).well_counts(:,well) = welldata(well).count;
    end
    data(vid).tempLog = fullTempList;
end

overtrackedFlies = sum(OG_flyCount-flyCount);
 
vid = 1; frame = 1;
movieInfo = VideoReader([baseFolder vidFolder '\' expName  '_' num2str(vid) '.avi']); %read in video
demoImg = rgb2gray(read(movieInfo,frame));

% visual check that the number of tracked flies is close to the actual number
nrow = 1; ncol = 2;
fig = getfig; set(fig, 'pos', [195 157 1413 646]);%[2160 185 1413 646]) <--evynpc
subplot(nrow, ncol, 2)
    h = histogram(OG_flyCount);
    h.FaceColor = Color('grey'); hold on
    histogram(flyCount); vline(nflies, 'w-'); 
    xlabel('Number tracked flies'); ylabel('Frame count')
    title(['Removed ' num2str(overtrackedFlies) ' overtracking points'])
subplot(nrow, ncol, 1)
% Take and save a snapshop pic of the arena w/ labels & well roi
%     outerColors = {'red', 'yellow', 'blue', 'green'};
    imshow(demoImg)
    axis tight square
    hold on
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        scatter(wellcenters(1,well),wellcenters(2,well), 75,...
            'MarkerFaceColor', kolor, 'MarkerEdgeColor', 'w') 
    end

    x = data(vid).x_loc(frame,:);
    y = data(vid).y_loc(frame,:);
    x(isnan(x)) = []; % remove empty tracks
    y(isnan(y)) = [];
    
    scatter(x,y, 15, 'b')

    l = legend(strrep(wellLabels,'_','-')); 
    set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.0959 0.7409 0.1271 0.1428]);
    viscircles(wellcenters',ones(1,4)*radii);
    titleName = strrep([folder ' ' expName arenaSel], '_',' ');
    title(titleName,'color', 'w')
    formatFig(fig, true,[nrow, ncol]);
   
% save and export figure:
if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
end  
save_figure(fig, [analysisDir expName arenaSel ' quality control'], '-png');

initial_vars = [initial_vars; 'radii'; 'welldata'; 'flyCount'];
clearvars('-except',initial_vars{:})
fprintf('\nNext\n')







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


%% Manual tracking cleanup option
currFrame = 0;
for vid = 1:nvids
    occupancy.frameROI(vid,1) = currFrame+1;
    nframes = size(data(vid).occupancy_matrix,2);
    currFrame = nframes + currFrame;
    occupancy.frameROI(vid,2) = currFrame;
end

% Find the video with the greatest avg over-count of fly points:
for vid = 1:nvids
    ROI = occupancy.frameROI(vid,:);
    numberCount(vid) = median(occupancy.flycount(ROI(1):ROI(2)));
end
% find the top 4 frames and display them
pullFrames = 5;
frameList = [];
vid = find(numberCount == max(numberCount));
ROI = occupancy.frameROI(vid,:);
numlist = occupancy.flycount(ROI(1):ROI(2));
[~,Idx] = (sort(numlist));
frameList = Idx(end-pullFrames+1:end); %pull the four highest overcounts

% pull up the images in order and click on points that ARE NOT FLIES:
movieInfo = VideoReader([baseFolder vidFolder '\' expName '_' num2str(vid) '.avi']); %read in video
for ii = 1:length(frameList)
    frame = frameList(ii);
    img = read(movieInfo,frame);
    fig = figure;
    imshow(img);
        hold on
        x = data(vid).x_loc(frame,:); x(isnan(x)) = [];
        y = data(vid).y_loc(frame,:); y(isnan(y)) = [];
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

% draw circles around the selected ROIs 
removalRadius = 10;
fig = figure; set(fig, 'color', 'k');
imshow(img)
hold on
for ii = 1:nmaskpoints
    scatter(pointLabels(ii).coord(1,:), pointLabels(ii).coord(2,:),20, 'filled')
end
R = removalRadius*ones(size(pointLabels(1).finalRound,2),1);
viscircles(pointLabels(1).finalRound', R ,'Color','r');
save_figure(fig, [analysisDir expName arenaSel ' overtracking points selected for deletion'], '-png');

% how many flies are in those circles?? TODO integrate the last OG fly count...
[og_flyCount, new_flyCount, raw_flyCount] = deal([]);
for vid = 1:nvids
%resize the data:
    x_loc = data(vid).x_loc; 
    y_loc = data(vid).y_loc;
    
    % how many flies OG
    og_flyCount = [og_flyCount; sum(~isnan(x_loc),2)];
    
    Xdim = size(x_loc); 
    Ydim = size(y_loc);
    X = reshape(x_loc,[numel(x_loc),1]); 
    Y = reshape(y_loc,[numel(y_loc),1]);

    % Find points within the masked region and turn to NaN
    for ii = 1:nmaskpoints
        loc = (((X-pointLabels(1).finalRound(1,ii)).^2 + (Y-pointLabels(1).finalRound(2,ii)).^2).^0.5)<=removalRadius;
        X(loc) = NaN; Y(loc) = NaN;
    end
    % Resize the data to OG structure:
    X = reshape(X,Xdim);
    Y = reshape(Y,Ydim);
    new_flyCount = [new_flyCount; sum(~isnan(X),2)];
    % raw fly count:
    raw_flyCount = [raw_flyCount; sum(data(vid).occupancy_matrix)'];
end
% number of tracked points:
overtracked(1) = sum(raw_flyCount)-sum(og_flyCount);
overtracked(2) = sum(og_flyCount)-sum(new_flyCount);
overtracked(3) = sum(raw_flyCount)-sum(new_flyCount);




% compare number of points same/different
nrow = 1;
ncol = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;
% PLOT TRACKING COUNT OVER TIME
fig = figure; set(fig, 'pos', [754 631 1936 707])
subplot(nrow, ncol, sb(1).idx)
hold on
plot(occupancy.time, raw_flyCount, 'color', Color('grey'))
plot(occupancy.time, og_flyCount,'color', Color('orangered'))
plot(occupancy.time, new_flyCount, 'color', Color('teal'))
hline(nflies, 'w')
ylabel('Fly count'); xlabel('time (min)')
subplot(nrow, ncol, sb(2).idx)
hold on
h = histogram(raw_flyCount);
h.FaceColor = Color('grey');
h = histogram(og_flyCount);
h.FaceColor = Color('orangered'); 
h = histogram(new_flyCount);
h.FaceColor = Color('teal');
vline(nflies, 'w')
xlabel('Number of flies')
ylabel('Frame count')
l = legend({['SLEAP ' num2str(mean(raw_flyCount)) ' avg'],...
            ['wells & arena masks ' num2str(mean(og_flyCount)) ' avg'],...
            ['wells, arena & manual ' num2str(mean(new_flyCount)) ' avg']});
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.7941 0.8152 0.1524 0.0877]); %

fig = formatFig(fig, true, [nrow, ncol], sb);

save_figure(fig, [analysisDir expName arenaSel ' overtracking points overview'], '-png');



vid = 1; frame = 1;
movieInfo = VideoReader([baseFolder vidFolder '\' expName  '_' num2str(vid) '.avi']); %read in video
demoImg = rgb2gray(read(movieInfo,frame));

% visual check that the number of tracked flies is close to the actual number
nrow = 1; ncol = 2;
fig = getfig; set(fig, 'pos', [195 157 1413 646]);%[2160 185 1413 646]) <--evynpc
subplot(nrow, ncol, 2)
    h = histogram(raw_flyCount);
    h.FaceColor = Color('grey'); hold on
    histogram(flyCount); vline(nflies, 'w-'); 
    xlabel('Number tracked flies'); ylabel('Frame count')
    title(['Removed ' num2str(overtrackedFlies) ' overtracking points'])
subplot(nrow, ncol, 1)
% Take and save a snapshop pic of the arena w/ labels & well roi
%     outerColors = {'red', 'yellow', 'blue', 'green'};
    imshow(demoImg)
    axis tight square
    hold on
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        scatter(wellcenters(1,well),wellcenters(2,well), 75,...
            'MarkerFaceColor', kolor, 'MarkerEdgeColor', 'w') 
    end

    x = data(vid).x_loc(frame,:);
    y = data(vid).y_loc(frame,:);
    x(isnan(x)) = []; % remove empty tracks
    y(isnan(y)) = [];
    
    scatter(x,y, 15, 'b')

    l = legend(strrep(wellLabels,'_','-')); 
    set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.0959 0.7409 0.1271 0.1428]);
    viscircles(wellcenters',ones(1,4)*radii);
    titleName = strrep([folder ' ' expName arenaSel], '_',' ');
    title(titleName,'color', 'w')
    formatFig(fig, true,[nrow, ncol]);

%%   
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
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.0169 0.8981 0.3122 0.0877]); %


%% Finding specific flies within the jumble

% Exemplary Trials
% ===================================================================================
loc = true(ntrials,1);

% take only positive plant relationship
loc(food(1).slope<0)=false;
% take only negative yeast slopes
loc(food(2).slope>0)=false;

flies = T(loc,:);

flyData = [food(1).slope(loc),food(2).slope(loc)];

fig = figure; hold on
for ii = 1:sum(loc)
    plot([1,2], flyData(ii,:),'Color', Color('powderblue'), 'linewidth', 1)
end
hline(0,'r--')
xlim([0.8, 2.2])
ylabel('Temp-distance relation')
xlabel('Food Type')
set(gca, 'xtick', [1,2], 'xticklabels', {'Plant', 'Yeast'})
title('Exemplary trials')
fig = formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' exemplar temp-distance trials'], '-png');


% Non-Exemplary Trials
% ===================================================================================
loc = true(ntrials,1);

% take only positive plant relationship
loc(food(1).slope>=0)=false;
% take only negative yeast slopes
loc(food(2).slope<=0)=false;

flies = T(loc,:);

flyData = [food(1).slope(loc),food(2).slope(loc)];

fig = figure; hold on
for ii = 1:sum(loc)
    plot([1,2], flyData(ii,:),'Color', Color('mistyrose'), 'linewidth', 1)
end
hline(0,'r--')
xlim([0.8, 2.2])
ylabel('Temp-distance relation')
xlabel('Food Type')
set(gca, 'xtick', [1,2], 'xticklabels', {'Plant', 'Yeast'})
title('Non-Exemplary trials')
fig = formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' non-exemplar temp-distance trials'], '-png');


% Positive attraction to food at warm temp Trials
% ===================================================================================
loc = true(ntrials,1);

% take only negative plant relationship
loc(food(1).slope>=0)=false;
% take only negative yeast slopes
loc(food(2).slope>=0)=false;

flies = T(loc,:);

flyData = [food(1).slope(loc),food(2).slope(loc)];

fig = figure; hold on
for ii = 1:sum(loc)
    plot([1,2], flyData(ii,:),'Color', Color('grey'), 'linewidth', 1)
end
hline(0,'r--')
xlim([0.8, 2.2])
ylabel('Temp-distance relation')
xlabel('Food Type')
set(gca, 'xtick', [1,2], 'xticklabels', {'Plant', 'Yeast'})
title('Stronger attraction to all food at warm temp trials')
fig = formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' positive warm relationship temp-distance trials'], '-png');



% Positive attraction to food at cold temp Trials
% ===================================================================================
loc = true(ntrials,1);

% take only negative plant relationship
loc(food(1).slope<=0)=false;
% take only negative yeast slopes
loc(food(2).slope<=0)=false;
sum(loc)
flies = T(loc,:);

flyData = [food(1).slope(loc),food(2).slope(loc)];

fig = figure; hold on
for ii = 1:sum(loc)
    plot([1,2], flyData(ii,:),'Color', Color('grey'), 'linewidth', 1)
end
hline(0,'r--')
xlim([0.8, 2.2])
ylabel('Temp-distance relation')
xlabel('Food Type')
set(gca, 'xtick', [1,2], 'xticklabels', {'Plant', 'Yeast'})
title('Stronger attraction to all food at COLD temps trials')
fig = formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' positive cold relationship temp-distance trials'], '-png');

%% Fly clustering by temp bin:

ROI = 1:40490; % take the first 225 minutes only | [] =  all data

if ntrials > 15
    warndlg('There are too many data points for this analysis'); return
end
sSpan = 360;
timeLen = zeros(1,ntrials);
for trial = 1:ntrials
    timeLen(trial) = length(data(trial).occupancy.flycount);
end
len = max(timeLen); %max number of frames (or time) within the selected experiments
uni_time = linspace(1,(len/3)/60,len); % universal time component...

sSpan = 360;
[clust, temp] = deal(nan([len,ntrials]));
all = [];
for trial = 1:ntrials
    
    x = data(trial).occupancy.temp';
    y = data(trial).occupancy.IFD./pix2mm;  % interfly distance in mm
    
    if ~isempty(ROI)
        x = x(ROI);
        y = y(ROI);
    end

%     %normalize to the minimum distance (1=most clustered)...this helps deal with the
%     %fly density variability
%     minY = min(y);
%     maxY = max(y);
%     z = abs(1-((y-minY)./(maxY-minY))); %1=most clustered 0 = no clustering
    % TODO: group the reponses by temp degree bins
    % clustering data
    clust(1:length(z),trial) = z;
    % temperature data
    temp(1:length(x),trial) = x;
    % both:
    all = [all; x,z'];
end
% sort the 'all' data by temperature dependence:
[~, idx] = sort(all(:,1));
sortedData = all(idx,:);
% avg and err of all the trials
Y_avg = nanmean(clust,2);
Y_err = std(clust,0,2);
nanloc = isnan(Y_err);
Y_err(nanloc) = [];
Y_avg(nanloc) = [];

xlimits = [8,20]; %TODO: change this if the temp range changes drastically!
ylimits = [0,1];
nrows = 1;
ncols = 4;
sb(1).idx = 1:2;
sb(2).idx = 3:4;
kolor = Color('slateblue');


%% Movement quantification & update....
% curr_file = 'G:\My Drive\Jeanne Lab\DATA\01.24.2022\Arena A\analysis\PlantOnly_1A timecourse data.mat';
clear 

pix2mm = 12.8;

%load excel file:
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
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
List.expID = eligible_files(fileIdx,3); 
List.arena = eligible_files(fileIdx,2); 

%get base folder pathway
baseFolder = getCloudPath;

initial_vars = who; initial_vars{end+1} = 'initial_vars'; initial_vars{end+1} = 'ii';

for ii = 1:length(fileIdx)
    inputPath = [baseFolder List.date{ii} '/Arena ' List.arena{ii} '/' List.expID{ii} List.arena{ii} ' timecourse data.mat'];
    
    % load data:
    tic
    load(inputPath)
    
    % --------------------- DATA UPDATE ---------------------
    % get movement size spacing
    xmin = data.centre(1)-data.r;
    xmax = data.centre(1)+data.r;
    ymin = data.centre(2)-data.r;
    ymax = data.centre(2)+data.r;
    nbins = 100;
    xedge = linspace(xmin,xmax,nbins+1);
    yedge = linspace(ymin,ymax,nbins+1);
    N = []; 
    for frame = 1:data.T.frame(end)
        X = data.x_loc(frame,:); X(isnan(X)) = [];
        Y = data.y_loc(frame,:); Y(isnan(Y)) = [];
        N(:,:,frame) = histcounts2(X,Y,xedge,yedge);
    end
    binDiff = diff(N,1,3);
    binDiff(binDiff<0) = 0;
    movement = squeeze(sum(sum(binDiff,1),2));
    binSizeinMM = ((2*data.r)/nbins)/pix2mm;
    movement = movement./binSizeinMM;
    movement = movement*3; % convert to per second WTIH 3fps videos
    % normalize for the number of flies actually tracked on the page:
    avg_movement = movement./data.T.flyCount(1:end-1);
    
    % plot change in movement 
    sSpan = 180;
    fig = figure;
    hold on
    plot(smooth(avg_movement,sSpan),'color', 'r')
    plot(smooth(data.occupancy.movement,sSpan),'color', 'k')
    legend({'new', 'old'})
    save_figure(fig,[baseFolder List.date{ii} '/Arena ' List.arena{ii} '/Movement Update'],'-png',true);
    
    % update data structures
    data.occupancy.movement = avg_movement;
    data.occupancy.baseMovement = movement;
    data.T.movement(1:end-1) = avg_movement;
    
    save(inputPath,'data','expData', 'expName', 'tempLog');
    
    toc
    disp(['Finished ' FileNames(fileIdx(ii))])  
    disp([num2str(ii) '/' num2str(length(fileIdx))])
    clearvars('-except',initial_vars{:})
end


%% Select experiments and update something about them

clear; close all; clc

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
fileIdx = listdlg('ListString', FileNames,'ListSize',[300,450]);
%pull the list of dates and arenas to be 
List.date = eligible_files(fileIdx,1);
List.expID = eligible_files(fileIdx,3); 


%get base folder pathway
finishedFiles = [];
baseFolder = getCloudPath;
for ii = 1:length(fileIdx)
    results = [];
    clear speed trackROI speedTracks SPEED  
    inputPath = [baseFolder List.date{ii} '/Analysis/' List.expID{ii} ' speed data.mat'];
    if ~any(strcmp(finishedFiles,[List.date{ii} ' ' List.expID{ii}]))
        load(inputPath)
        

        % ------ CHANGE THIS ----------
        %trackroi not saved, load from individual prior data set
        if ~exist('trackROI','var')
            loadInd = true;
        end

        SPEED = speed;
        for arena = 1:4
            dataPath = [baseFolder List.date{ii} '/Arena ' Alphabet(arena) '/' List.expID{ii} ' speed data.mat'];
            speed = SPEED(arena);
            
            if loadInd
                speedTracks = load(dataPath,'speedTracks');
            else
                speedTracks = trackROI(:,arena);
            end
            save(dataPath,'speed', 'speedTracks');
        end
        results = 'Saved Data';
        % -----------------------------
    end
    finishedFiles{ii} = [List.date{ii} ' ' List.expID{ii}];
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

disp('Done with full set')



%% 

theta = deg2rad(35);
h = (150)/25.4;

D = 2*tan(theta)*h;
disp(D)


%% 

% PRE STEP 1
parameters.ArenaA.well_1 = 'Caviar';
parameters.ArenaB.well_1 = 'Caviar';
parameters.ArenaC.well_1 = 'Caviar';
parameters.ArenaD.well_1 = 'Caviar';

% PRE STEP 2
arenaData(1).wellLabels{1} = 'Caviar';
arenaData(2).wellLabels{1} = 'Caviar';
arenaData(3).wellLabels{1} = 'Caviar';
arenaData(4).wellLabels{1} = 'Caviar';

expData.parameters.ArenaA.well_1 = 'Caviar';
expData.parameters.ArenaB.well_1 = 'Caviar';
expData.parameters.ArenaC.well_1 = 'Caviar';
expData.parameters.ArenaD.well_1 = 'Caviar';


% DONT REMEMBER
data.wellLabels{4} = 'Empty';
    
    

    
%% Get formatFig to work with boxplots etc 

% colorsforboxplot=colorsforboxplot.';
fig = figure;
boxplot([(1:25)',(1:25)',(1:25)'],[2,4,6])
formatFig(fig)

odorCombos = {'dummy1','dummy2','dummy3'};
colorsforboxplot = Color('green','red',3);
allmeanMax = [(1:25)',(1:25)',(1:25)'];

fig4= figure(4);
z = 1:3;
boxplot(allmeanMax,z)
h = findobj(gca,'Tag','Box');
for zz=1:length(h)
    patch(get(h(zz),'XData'),get(h(zz),'YData'),colorsforboxplot(zz,:),'FaceAlpha',.5);
end
hold on

ylabel('normalized max df/f')
xlabel('odor')
set(gca,'XTickLabel',odorCombos)

formatFig(fig4,false,[1,1]);




%%
vidPath = "F:\Evyn\DATA\10.03.2023\C2ExpName_1.avi";
movieInfo = VideoReader(vidPath); %read in video
disp(['Dimensions: ' num2str(movieInfo.Width) ' x ' num2str(movieInfo.Height)])



%% UPDATE THE NAMES IN FOOD WELLS

%Blank all of them
for i = 1:4
    for ii = 1:4
        parameters.(['Arena' Alphabet(i)]).(['well_' num2str(ii)]) = 'Empty';
    end
end

%Readd the appropriate food labels:
parameters.ArenaA.well_3 = '0.5M_glucose_15%ACV';
parameters.ArenaB.well_3 = '0.5M_glucose_15%ACV';
parameters.ArenaC.well_1 = '0.5M_glucose_5%ACV';
parameters.ArenaD.well_1 = '0.5M_glucose_5%ACV';

% Clear extra variables:
clear i ii ans

%% UPDATE THE GENOTYPE

%Blank all of them
for i = 1:4
        parameters.(['Arena' Alphabet(i)]).genotype = 'UAS-Kir2.1_B_F8';
end

% %Custom (not all) updates to genotype the appropriate food labels:
% parameters.ArenaA.genotype = '0.5M_glucose_15%ACV';
% parameters.ArenaB.genotype = '0.5M_glucose_15%ACV';
% parameters.ArenaC.genotype = '0.5M_glucose_5%ACV';
% parameters.ArenaD.genotype = '0.5M_glucose_5%ACV';

% Clear extra variables:
clear i ii ans

%% Double templog 
clear 

basefolder = 'G:\My Drive\Jeanne Lab\DATA\01.21.2024';

% 1) get the list of experiments in the date folder
videoList = dir([basefolder '\*_1.avi']);
[exp_names, start_times] = deal([]);
for i = 1:length(videoList)
    exp_names{i} =  videoList(i).name(1:end-6);
    exp_start_time{i} = videoList(i).date(end-7:end-3);
end

% 2) get the list of temperature log files
tempLogList = dir([basefolder '\*_RampLog.csv']);
% get experiment names that have a temp log file
tempLog_names = [];
for i = 1:length(tempLogList)
    tempLog_names{i} =  tempLogList(i).name(1:end-12);
end

% 3) Find the experiments that do not have a temperature log
idx = 0;
for i = 1:length(exp_names)
    if ~any(strcmp(exp_names{i}, tempLog_names)) %no matches
       %note the unmatched experiment name
       idx = idx + 1;
       no_log_list{idx} = exp_names{i};
    end
end

% 4) Find any experiments that match the start time of the missing log experiments:
if idx > 0
    for i = 1:idx %for each unmatched experiment
        sel_start_time =  exp_start_time{strcmp(no_log_list{i},exp_names)}; % start time of the selected experiment
        search_loc = find(~strcmp(no_log_list{i},exp_names)); % index of other experiments
        for ii = search_loc
              if str2double(strrep(sel_start_time,':','')) == str2double(strrep(exp_start_time{ii},':','')) % start time match
                    source_file = [basefolder '\' exp_names{ii} '_RampLog.csv'];
                    destination_file = [basefolder '\' no_log_list{i} '_RampLog.csv'];
                    copyStatus = copyfile(source_file, destination_file);
                    if ~copyStatus
                        h = warndlg([exp_names{ii} ' not copied']);
                        uiwait(h)
                    end
                    return
              end
        end
    end
end
                 


                
                
%% Static temperature hold tuning curve

% get the avg distance over the whole experiment
roi = [1000 159935];
[dist_mean, dist_err] = deal([]);
for i = 1:num.exp
    dist_all = mean(grouped(i).dist.all(roi(1):roi(2),:),1,'omitnan');
    dist_err(i) = std(dist_all);
    dist_mean(i) = mean(dist_all);
end

e = errorbar([17,20,23,25],dist_mean,dist_err);
e.Marker  = 'o';
e.Color = 'r';
e.MarkerSize = 10;



%% 

% % Create a figure and axes
% f = figure; 
% imshow(read(movieInfo,ii(jj)));
% 
% % Display an image or plot something
% % For demonstration, let's just use a blank axis
% axis(ax, [0 10 0 10]);
% 
% % Enable interactive addition of points
% % The 'PositionConstraintFcn' can be used to restrict movement.
% h = drawpoint('Color','r', 'Selectable', true);
% 
% % Example to keep the point within bounds
% % Uncomment the following line to see how it works
% % h.PositionConstraintFcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
% 
% % ------------------------------------------------------------- 
% % Initial point coordinates
% x = 5;
% y = 5;
% 
% % Plot the point
% hPoint = plot(x, y, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% 
% % Enable dragging functionality
% set(hPoint, 'ButtonDownFcn', @(src, event) startDragFcn(src, ax));
% set(gcf, 'WindowButtonUpFcn', @stopDragFcn);
% 
% % Dragging function
% function startDragFcn(src, ax)
%     set(gcf, 'WindowButtonMotionFcn', @(fig, event) draggingFcn(fig, src, ax));
% end
% 
% % Function to execute while dragging
% function draggingFcn(fig, src, ax)
%     pt = ax.CurrentPoint;
%     src.XData = pt(1,1);
%     src.YData = pt(1,2);
% end
% 
% % Stop dragging
% function stopDragFcn(fig, event)
%     set(fig, 'WindowButtonMotionFcn', '');
% end


%% Mac testing
% fig = open('/Users/evyndickinson/Desktop/untitled folder/demofig.fig');
% 
% save_loc = '/Users/evyndickinson/Desktop/untitled folder/demofig';
% export_fig(fig, save_loc, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');


%% Write to excel replacement functions

% 
% xlswrite(XL, {videoStartTime}, 'Exp List', [Alphabet(Excel.starttime) num2str(XLrow)]);
% 
% writematrix(M,'M.xls','Sheet',2,'Range','A3:E8')
% 
% testFile = '/Users/evyndickinson/Desktop/Quad Bowl Experiments.xlsx';
% 
% writematrix(nflies(arena),testFile,'Sheet','Exp List','Range',[Alphabet(Excel.numflies) num2str(XLrow(arena))])
% 
% 
% % videoStartTime = 
% 
% writecell({videoStartTime},XL,'Sheet','Exp List','Range',[Alphabet(Excel.starttime) num2str(XLrow)])


%% Select data groups to compare
% add matlabroot folder to the directory path
% addpath(genpath('C:\matlabroot'));

% % clear; close all; clc
% baseFolder = getCloudPath;
% 
% switch questdlg('Load existing data?','Quad Step 4 data processing','Yes','No','Cancel','Yes')
%     case 'Cancel'
%         return
%     case 'Yes' % have user select which preexisting structure to load
%         list_dirs = dir([baseFolder 'Grouped Data Structures\']);
%         list_dirs = {list_dirs(:).name};
%         list_dirs(1:2) = [];
%         dirIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[350,650]);
% 
%         try expGroup = list_dirs{dirIdx}; %load experiment groups selected
%         catch
%             disp('No group selected')
%             return
%         end
%         saveDir = [baseFolder 'Grouped Data Structures\' expGroup '\'];
%         load([saveDir expGroup ' data.mat']);
%         disp([expGroup ' loaded'])
%         %check for temperature protocol assignments (make back-compatible)
%         for i = 1:num.exp
%             if ~isfield(data(i),'temp_protocol')
%                 data(i).temp_protocol = data(i).T.TempProtocol{1};
%             end
%             if isempty(data(i).temp_protocol)
%                 data(i).temp_protocol = data(i).T.TempProtocol{1};
%             end
%         end
%     case 'No'
%         % Select processed data structures to compare:
%         structFolder = [baseFolder 'Data structures\'];
%         list_dirs = dir(structFolder);
%         list_dirs = {list_dirs(:).name};
%         list_dirs(1:2) = [];
%         expIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'multiple','ListSize',[300,450]);
%         expNames = list_dirs(expIdx); %name of experiment groups selected
%         num.exp = length(expIdx);  %number of groups selected
% 
%         % Load selected experiment data groups
%         for i = 1:num.exp
%             % get field list for loading data:
%             dummy = load([structFolder expNames{i} '\' expNames{i} ' post 3.1 data.mat']);
%              if ~isfield(dummy,'hold_exp') % 1/24/24 updates for hold temperature experiments
%                    dummy.hold_exp = false; % account for new data structures
%                    dummy.temp_protocol = dummy.T.TempProtocol{1};
%              end
%             data(i) = dummy;
%             % try data(i) = dummy;
%             % catch % TODO -- better automate search for missing fields and options
%             % % data(i) = load([structFolder expNames{i} '\' expNames{i} ' post 3.1 data.mat']);
%             %
%             %     data(i) = dummy;
%             % end
%         end
% 
%         clear list_dirs expIdx dirIdx
%         % Set up base variables
%         initial_vars = who;
%         initial_vars = [initial_vars(:); 'initial_vars'; 'grouped'; 'expGroup'; 'saveDir'; 'mat';'expOrder'];
%         initial_vars = unique(initial_vars);
% 
%         % Save data / make new grouped data folder
%         switch questdlg('Select data saving format:','','new structure','existing structure', 'cancel','new structure')
%             case 'new structure'
%                 expGroup = char(inputdlg('Structure name:'));
%                 saveDir = [baseFolder 'Grouped Data Structures\' expGroup '\'];
%                 if ~exist(saveDir,'dir')
%                     mkdir(saveDir);
%                 end
%                 save([saveDir expGroup ' data.mat'],'-v7.3');
%                 disp([expGroup ' saved'])
%             case 'existing structure'
%                 list_dirs = dir([baseFolder 'Grouped Data Structures\']);
%                 list_dirs = {list_dirs(:).name};
%                 list_dirs(1:2) = [];
%                 dirIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[300,450]);
%                 expGroup = list_dirs{dirIdx}; %name of experiment groups selected
%                 saveDir = [baseFolder 'Grouped Data Structures\' expGroup '\'];
%                 save([saveDir expGroup ' data.mat'],'-v7.3');
%                 disp([expGroup ' saved'])
%             case 'cancel'
%                 return
%         end
% end
% disp(expNames')


%% Can we pair speed and position?
% use the trackROI data to determine the location or distance from center to
% calculate speed for each location?
% 
% temp = load([filePath expID{trial} arenas{trial} ' timecourse data.mat'], var2load{:});
% 
% % load speed data
%             temp = load([filePath expID{trial} ' speed data.mat']);
%             data(trial).speed = temp.speed;


%% 
% 
% for i = 1:num.exp
%     exp = expOrder(i)
%     disp(grouped(exp).name)
% end
% 

%% Extract timestamps and templogs from HighRes data sets

% look for a raw folder
% or look for regular files in the main folder (*file#.mat)
% print out a list of files that still have the raw files that need to be
% deleted -- add this as a column to the excel file
% clear; close all; clc

[excelfile, Excel, xlFile] = load_HighResExperiments;

rootDir = getDataPath(5, 2, 'Select location for data');
paths = getPathNames;

for i = 2:length(excelfile)
    date = excelfile{i,Excel.date};
    folder = excelfile{i,Excel.expID};
    rootFolder = [rootDir, paths.courtship, date '/' folder '/'];

    % Look for an existing timestamps data file:
    if exist([rootFolder 'raw timestamps.mat'],'file')
        disp([date ' ' folder ' Y'])
        if ~isExcelFileOpen(xlFile,true)
            try writecell({'Y'},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.raw) num2str(i)]);
            catch 
                disp('couldn''t write to excel:  manually update')
            end
        else
            disp('couldn''t write to excel:  manually update')
        end
        continue
    end

    % Look for raw data: 
    if isfolder([rootFolder 'raw'])
        basePath = [rootFolder 'raw/'];
    else
        basePath = rootFolder;
    end
    a = dir([basePath '*file*.mat']);
    if ~isempty(a)
        tic
        timestamps = struct;
        timestamps.time = [];
        timestamps.tempLog = [];
        for ii = 1:length(a)
            dummy = load([basePath 'file' num2str(ii) '.mat'],'timestamp', 'tempLog');
            timestamps.time{ii,1} = dummy.timestamp;
            timestamps.tempLog(ii,:) = dummy.tempLog;
        end
        % save the data 
        save([rootFolder 'raw timestamps.mat'],'timestamps');
        statusStr = 'Y'; % as in the data is now saved
        disp([date ' ' folder ' ' statusStr])
        toc
    else
        statusStr = 'N'; % as in this data does not exist any longer
        disp([date ' ' folder ' ' statusStr])
    end

    % write that the data has been processed
    if ~isExcelFileOpen(xlFile,true)
        try writecell({statusStr},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.raw) num2str(i)]);
        catch 
            disp('couldn''t write to excel:  manually update')
        end
    else
        disp('couldn''t write to excel:  manually update')
    end
end
    

%% 

%% test the difference in start times for an experiment
startDir = uigetdir;
% startDir = ['I:\Data\02.28.2025\Berlin_courtship_F_LRR_caviar_ramp4\'];

tic
time = NaT([1,768]);
for i = 1:768
    dummy = load([startDir, 'file' num2str(i)],'timestamp');
    time(i) = dummy.timestamp; 
    disp(i)
end
toc

a = diff(time);

figure; plot(a)

%% Quick visual check of the arena plate and named assignment 4/22/2025
clear

baseFolder = getDataPath(2,0);

[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
processed = [];

loc1 = find(cellfun(@isnan,excelfile(2:end,Excel.plate)))+1; % initial rows with no plate information
for i = 1:length(loc1)

    % find all trials from that same plate: 
    loc = loc1(i);

    % check that the information hasn't been updated already! 
    if ~any(processed==loc)

        % find all sister trials
        date_str = excelfile{loc, Excel.date};
        expName = excelfile{loc, Excel.expID};
        if any(isnan(expName))
            processed = [processed; loc];
            continue
        end
        loc2 = strcmp(excelfile(:,Excel.date),date_str); % locations with same date
        loc3 = strcmp(excelfile(:,Excel.expID),expName); % locations with same exp name
        loc_all = find(loc3 & loc2); % location of all trials for this plate 
        
        % display image of the arena & choose the plate name
        movieInfo = VideoReader([baseFolder date_str '\' expName '_1.avi']); %#ok<*TNMLP> (error suppression message)
        img = rgb2gray(read(movieInfo,1));
        fig = figure;
        imshow(img)
        title(strrep(expName,'_',' '))
        switch questdlg('What arena is this?', expName, 'Plate 1 (old)', 'Plate 2 (new)','Skip', 'Plate 1 (old)')
            case 'Plate 1 (old)'
                    plate_str = 1;
            case 'Plate 2 (new)'
                    plate_str = 2;
            case {'Skip', ''}
                return
        end
        
        close(fig)
        
        isExcelFileOpen(xlFile);
        ntrials = length(loc_all);
        xlRange = [Alphabet(Excel.plate) num2str(loc_all(1)) ':' Alphabet(Excel.plate) num2str(loc_all(end))];
        plateList = repmat({plate_str}, ntrials, 1);
        writecell(plateList, xlFile, 'Sheet', 'Exp List', 'Range', xlRange); % write information to excel file

        disp([date_str ' | ' xlRange ' | plate: ' num2str(plate_str) ' i: ' num2str(i)])

        % add to processed list
        processed = [processed; loc_all];
    end
end

%% Look for videos with poor tracking (low res setup)

fig = getfig;
for i = 1:4 % for each arena
    kolor = arenaData(i).color;
    subplot(2,2,i)
    hold on
    x = T.vidNums;
    y = T.flyCount(:,i);
    scatter(x,y,50,kolor,'filled')
    xlabel('video number')
    ylabel('number of tracked flies')
    title(['Arena ' Alphabet(i)],'color', 'w')
end
formatFig(fig,true,[2 2]);

% determine the avg number of tracked flies per video: 
fig = getfig;
for vid = 1:nvids
    loc = T.vidNums==vid;
    a = flyCount(loc,:);
    b = mean(a,1,'omitnan');
    c = std(a,0,1,'omitnan');
    for i = 1:4
        kolor = arenaData(i).color;
        subplot(2,2,i); hold on
        errorbar(vid,b(i),c(i),'Color',kolor)
        scatter(vid, b(i),50, kolor,'filled')
    end
end
formatFig(fig, true,[2 2]);
for i = 1:4
    subplot(2,2,i)
    xlabel('video number')
    ylabel('avg tracked flies')
    title(['Arena ' Alphabet(i)],'color', 'w')
end


%% 
frame = frames(rr).idx;
x = grouped(exp).(pos_type).trial(trial).x;


%% Determine how many trials we may need to scrap if they are Long and have food present on plate 1...

[excelfile, Excel, xlFile] = load_QuadBowlExperiments;

% look for trials that meet these conditions: [LTS or temp hold], [caviar], [plate 2]
protocols = excelfile(2:end,Excel.protocol);
proto_idx = strcmp(protocols, 'Large_temp_sweep_15_35');
food_idx = false(size(proto_idx)); %initialize as an empty false vector
for well = 1:4
    food = excelfile(2:end,Excel.(['well_' num2str(well)]));
    food_idx = food_idx | strcmpi(food, 'Caviar');
end
plate = [excelfile{2:end,Excel.plate}]';
plate_idx = plate==2;

% trials to toss: 
bad_idx = proto_idx & food_idx & plate_idx;
sum(bad_idx)

% what are the trials that we would need to scrap: 
badTrials = [excelfile([false; bad_idx],Excel.trialID), excelfile([false; bad_idx],Excel.genotype)];





















