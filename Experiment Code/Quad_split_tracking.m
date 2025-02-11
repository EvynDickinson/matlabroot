
clear
% Split the tracking by arena region and assign it to the appropriate arena

%% Load data
dataDir = 'G:\My Drive\Jeanne Lab\DATA\02.25.2022\Tracking\';
arenaData = struct;

% Find the number of videos:
vidList = dir([dataDir '\*.avi']);
nvids = length(vidList);

% Load tracking data
data = struct;
for vid = 1:nvids
    filePath = [dataDir vidList(vid).name '.predictions.analysis.h5']; % [baseFolder vidFolder '\' expName '_' num2str(vid) '.h5'];
    data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
    data(vid).tracks = h5read(filePath,'/tracks');
end; clear filePath


% Pull image from first video to use as arena template
movieInfo = VideoReader([dataDir vidList(1).name]); %read in video
demoImg = rgb2gray(read(movieInfo,1));

% Select the well centers from top right corner and go clockwise:
wellcenters = readPoints(demoImg,16); % get locations of the wells from the image
r = 435; % radius of the arena
arenaIdx = {'B', 'D', 'A', 'C'}; % order of arenas on video going clockwise from upper left corner

fig = figure; 
imshow(demoImg); set(fig, 'color', 'k')
hold on
% Screen any points outside the arena
for arena = 1:4
    arenaData(arena).name = ['Arena ' arenaIdx{arena}];
    % Divide well centers by arena
    roi = (arena-1)*4 + 1 : (arena-1)*4 + 4;
    arenaData(arena).wellcenters = wellcenters(:,roi);
    
    % Find the center of each arena
    x1 = arenaData(arena).wellcenters(1,1:2:4);
    y1 = arenaData(arena).wellcenters(2,1:2:4);
    x2 = arenaData(arena).wellcenters(1,2:2:4);
    y2 = arenaData(arena).wellcenters(2,2:2:4);
    [xi,yi] = polyxpoly(x1,y1,x2,y2);
    arenaData(arena).centre = [xi;yi];

    % Overlay the arena screen:
    viscircles(arenaData(arena).centre', r, 'color', Color('white'));
end

% TODO: save this image to the analysis folder

%% Get fly counts:

% Number of flies:
nflies = []; %TODO -- update this to actually read the information from excel
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
            PT(jj).frame = readPoints(demoImg);
            for arena = 1:4
                centre = arenaData(arena).centre;
                X = PT(jj).frame(1,:);
                Y = PT(jj).frame(2,:);
                % find points within each arena sphere
                loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)<=r;
                nflies(arena,jj) = sum(loc);
            end
        end
        disp('Number of flies counted')
        A = nflies(4,:)'; B = nflies(1,:)'; C = nflies(3,:)'; D = nflies(2,:)';
        numFlies = [A,B,C,D];
        disp(table(A,B,C,D))
        count_match = (diff(numFlies,2,1)==0);
        if all(count_match)
            nflies = nflies(1); 
            T = false;
        else
            switch questdlg(['Nonmatching fly counts, manually select number?: ' num2str(nflies)])
                case 'Yes'
                    for arena = 1:4
                        arenaNums{arena} = num2str(median(numFlies(:,arena)));
                        if ~count_match(arena)
                            arenaNums{arena} = 'NaN';
                        end
                    end
                    prompt = {'Arena A','Arena B', 'Arena C', 'Arena D'}; dlgtitle = 'Input';
                    answer = inputdlg(prompt,dlgtitle,[1 25],arenaNums);
                    for arena = 1:4
                        arenaData(arena).nflies = str2double(answer{arena});
                    end
                    T = false;
                case 'No'
                    T = true;
                case 'Cancel'
                    return
            end
        end
    end
    % write number of flies into the excel sheet %TODO : update this to
    % work with multiple arenas
    try
        xlswrite(xlFile, {num2str(nflies)}, 'Exp List', [Alphabet(Excel.numflies) num2str(XLrow)]);
    catch
        h = warndlg('Close Experiment Summary excel file and then close this warning box');
        uiwait(h)
        xlswrite(xlFile, {num2str(nflies)}, 'Exp List', [Alphabet(Excel.numflies) num2str(XLrow)]);
    end
end
fprintf(['\nNumber of flies: ' num2str(nflies) '\n'])

%% Check tracking quality: %TODO this whole section...

% Find food wells and their identities %TODO
h = warndlg('Select the well locations starting at "9 o''clock" and proceeding clock wise'); uiwait(h)
% label the wells:
for arena = 1:4
    wellLabels = {expData.parameters.(['Arena' arenaSel]).well_1;...
                  expData.parameters.(['Arena' arenaSel]).well_2;...
                  expData.parameters.(['Arena' arenaSel]).well_3;...
                  expData.parameters.(['Arena' arenaSel]).well_4};   
    disp(wellLabels)
end
arenaData(arena).wellLabels = wellLabels;
% wellcenters = readPoints(demoImg,4); % get locations of the wells from the image

initial_vars = [initial_vars; 'arenaData'; 'wellLabels'; 'wellcenters'; 'demoImg'];
clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% find number of frames per video
currFrame = 0;
for vid = 1:nvids
    occupancy.frameROI(vid,1) = currFrame+1;
    nframes = size(data(vid).occupancy_matrix,2);
    currFrame = nframes + currFrame;
    occupancy.frameROI(vid,2) = currFrame;
end

%% Show tracking points on frames:
scaleFactor = 1;
vid = 1;
frame = 1;
movieInfo = VideoReader([dataDir vidList(vid).name]); %read in video
demoImg = rgb2gray(read(movieInfo,frame));



for vid = 1:nvids
    [raw_flyCount, mid_flyCount] = deal([]);
    
    % raw_flyCount = [raw_flyCount; sum(data(vid).occupancy_matrix)'];
    
    % all data for head tracked location
    headData = squeeze(data(vid).tracks(:,1,:,:));
    
    % x-y coordinates of flies for each frame
    x_loc = squeeze(headData(:,1,:)).*scaleFactor;
    y_loc = squeeze(headData(:,2,:)).*scaleFactor;
    data(vid).x_loc_raw = x_loc;
    data(vid).y_loc_raw = y_loc;
    nframes = size(headData,1); 
    
    % --------- Find data points that lie within each arena ---------
    Xdim = size(x_loc);
    Ydim = size(y_loc);
    %resize the data:

    for arena = 1:4
        X = reshape(x_loc,[numel(x_loc),1]);
        Y = reshape(y_loc,[numel(y_loc),1]);
        centre = arenaData(arena).centre;
        % find points within each arena sphere
        loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)>=r;
        X(loc) = NaN;
        Y(loc) = NaN;
        %resize the data:
        X = reshape(X,Xdim);
        Y = reshape(Y,Ydim);
    
        % save tracking points into structure:
        data(vid).(arenaIdx{arena}).x_loc = X;
        data(vid).(arenaIdx{arena}).y_loc = Y;
    
        % number of flies within this arena:
        raw_flyCount = [raw_flyCount; sum(~isnan(X),2)];
        
    end
end
    
arenaData(arena).flyCount = raw_flyCount;


%% Check tracking locations:

fig = figure;
imshow(demoImg)
hold on
scatter(data(1).x_loc_raw(frame,:),data(1).y_loc_raw(frame,:),10,'yellow', 'filled')

% TODO:  show a random smattering of these frames and the alignment of the
% tracked points

%%
arena = 1;

% maskPath = [analysisDir expName arenaSel ' ArenaMask.mat'];
% if ~isfile(maskPath)
%     idx = find(~strcmpi('Empty', wellLabels));
%     if ~isempty(idx)
%         % draw out the masks:
%         for ii = 1:length(idx)
%             wellName = strrep(wellLabels{idx(ii)}, '_', '-');
%             % label the ROI
%             f = figure;
%             imshow(demoImg)
%             title(['Outline the ' wellName ' well'])
%             roi = drawpolygon;
%             mask(ii).n = roi.Position;
%             uiwait(f)
%         end
%         % save the food masks:
%         save(maskPath, 'mask')
%     end
% else
%     load(maskPath)
% end

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
scaleFactor = 1.31; % offset from the video encoding compression

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



%% 

% Load tracking data'
vidPath = 'G:\My Drive\Jeanne Lab\DATA\02.24.2022\PlantFood_LinearRamps_3.avi';
filePath = [vidPath '.predictions.analysis.h5']; % [baseFolder vidFolder '\' expName '_' num2str(vid) '.h5'];
occupancy_matrix = h5read(filePath,'/track_occupancy');
tracks = h5read(filePath,'/tracks');

















