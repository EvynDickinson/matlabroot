

%% TODOs 
% load data
% load the parameter files
% check the number of tracks
% check that the same fly is aligned across vidoes
% add position data to the temperature matrix
% calculate new parameters from the data
% plot data
% create predictions of courtship periods

%% PARAMETERS TO SHOW: 

% 2) each fly speed
% 3) body angle between flies
% 4) male fly wing angle

%% Load tracking points for courtship
clear; clc;
baseFolder = getDataPath(6,0);
dateDir = (selectFolder(baseFolder));
trialDir = selectFolder([baseFolder, dateDir{:}]);
baseDir = [baseFolder, dateDir{:} '\', trialDir{:} '\'];

fileList = dir([baseDir '*alignment table.mat']);
load([baseDir, fileList(1).name]) % load the parameters and temp table
nvids = parameters.nVids; % number of videos
nBP = 5; % num of body parts


% load the video tracks
fileList = dir([baseDir '*.h5']);
data = [];
% TODO: quick check to make sure this matches the expected number of videos...
for vid = 1:nvids
    filePath = [baseDir, 'compiled_video_' num2str(vid) '.avi.predictions.slp.h5'];
    data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
    data(vid).tracks = h5read(filePath,'/tracks');
    data(vid).node_names = h5read(filePath, '/node_names');
    data(vid).track_names = h5read(filePath, '/track_names');
    dataIn = true; %log successful data load
end

node_names = {'head', 'center', 'abdomen', 'right_wing', 'left_wing'}; % currently assuming these are stable across all videos

% for vid = 1:nvids
%    vidBase = [baseFolder folder '/' expName '_' num2str(vid)]; 
%    fp_endings = {'.h5', '.avi.predictions.slp.h5', '.avi.predictions.analysis.h5'};
%    dataIn = false; % switch for logging data
%    for suf = 1:length(fp_endings)
%        filePath = [vidBase fp_endings{suf}];
%        if exist(filePath,'file')
%            data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
%            data(vid).tracks = h5read(filePath,'/tracks');
%            dataIn = true; %log successful data load
%             continue
%        end
%    end
%    if ~dataIn
%         disp(vidBase)
%         h = warndlg('Warning: file not found');
%         uiwait(h)
%         return
%         % If file is corrupt or permanently missing...fill with blank NaN
%         % data
% %         data(vid).occupancy_matrix = [];
% %         data(vid).tracks = [];
%    end
% end; clear filePath fp_endings suf vidBase

%% Quality control

ntracks = [];
for vid = 1:parameters.nVids
    ntracks(vid) = size(data(vid).occupancy_matrix,1);
end


%% 
% tracks (frame, body points, XY, fly)

%% Find the x-y coordinates for each fly

% Set up empty tracks for each of the body parts that are being labeled?
w = 0;
for vid = 1:nvids
    if isempty(data(vid).occupancy_matrix)
        continue
    end
    w = max([size(data(vid).occupancy_matrix,1),w]);
end

% load fly locations into the X & Y matrices
positions = [];
for i = 1:nBP
    [positions.(node_names{i}).X, positions.(node_names{i}).Y]  = deal([]);
end
for vid = 1:nvids
    if isempty(data(vid).occupancy_matrix)
        continue
    end
    % ---- fly tracked locations -------
    for bp = 1:length(data(vid).node_names)
        raw = squeeze(data(vid).tracks(:,bp,:,:));
        % x-y coordinates of flies for each frame
        x_loc = squeeze(raw(:,1,:));
        y_loc = squeeze(raw(:,2,:));
        data(vid).(node_names{bp}).x =  x_loc;
        data(vid).(node_names{bp}).y =  y_loc;
        positions.(node_names{bp}).X = autoCat(positions.(node_names{bp}).X, x_loc,true);
        positions.(node_names{bp}).Y = autoCat(positions.(node_names{bp}).Y, y_loc,true);
    end
end


%% Align the tracks from video to video: 
% display image of fly postion at the end of the first video and then the
% next chunk for the subsequent video to see if they align
buff = 4;
vid_align = false([1,nvids]);
vid_align(1) = true;
for vid = 1:nvids-1
    frame_end = size(data(vid).occupancy_matrix,2);
    eROI = frame_end-buff:frame_end;
    sROI = 1:1+buff;
    
    fig = figure;
    hold on
    for i = eROI
        plotFlySkeleton(fig, data(vid).tracks(i,:,1,1),data(vid).tracks(i,:,2,1),Color('dodgerblue'),true);
        plotFlySkeleton(fig, data(vid).tracks(i,:,1,2),data(vid).tracks(i,:,2,2),Color('deeppink'),true); 
    end
    for i = sROI
        plotFlySkeleton(fig, data(vid+1).tracks(i,:,1,1),data(vid+1).tracks(i,:,2,1),Color('blue'),false);
        plotFlySkeleton(fig, data(vid+1).tracks(i,:,1,2),data(vid+1).tracks(i,:,2,2),Color('pink'),false); 
    end
    axis square
    switch questdlg('Are the tracks aligned?')
        case 'Yes'
        vid_align(vid+1) = true;
        case 'Cancel'
            return
        case ''
            return
    end
    close(fig)
    disp([num2str(vid) '/' num2str(nvids-1)])
end

state_switch = false;
state = true;
for i = 2:nvids
    a = vid_align(i);
    b = state(i-1);
    if a && b
        state_switch(i) = false;
        state(i) = true;
    elseif ~a && ~b
        state_switch(i) = false;
        state(i) = true;
    elseif ~a && b
        state_switch(i) = true;
        state(i) = false;
    elseif a && ~b
        state_switch(i) = true;
         state(i) = false;
    end
    % disp(state);
end
disp(state)
switch_state = ~state;




%% Targeted look at heating and cooling
pix2mm = 0.0305; %conversion from pixels to mm for these videos
IFD = [];
speed = [];
for vid = 1:nvids
% TODO determine which track is likely the male and which is likely the female :
% do this by comparing the avg size of the flies
    % inter-fly-distance from the fly's center point
    x1 = data(vid).center.x(:,1);
    y1 = data(vid).center.y(:,1);
    x2 = data(vid).center.x(:,2);
    y2 = data(vid).center.y(:,2);
    IFD = [IFD; (sqrt((x1-x2).^2 + (y1-y2).^2)).*pix2mm];

    % speeds
    D1 = (sqrt((x1(1:end-1)-x1(2:end)).^2 + (y1(1:end-1)-y1(2:end)).^2)).*pix2mm; 
    D2 = (sqrt((x2(1:end-1)-x2(2:end)).^2 + (y2(1:end-1)-y2(2:end)).^2)).*pix2mm;
    speed = [speed; (D1./(1/parameters.FPS)),(D2./(1/parameters.FPS))];
end


fig = getfig('',1); 

plot(T.time,IFD,'color', Color('teal'),'linewidth', 1)
h_line(5,'r','--',1)


%%






%% Distance, speed, and speed correcation between flies over full video course
% use the center body point to determine between-fly distance
fps = parameters.FPS;
time = T.time;
sSpan = 90; % 3 seconds
[foreColor,backColor] = formattingColors(false); %get background colors

r = 3; c = 1;
xlimit = [0,120];

fig = getfig('',1,[1032 1042]);
% INTERFLY DISTANCE
subplot(r,c,1) 
    plot(T.time,IFD,'color', foreColor,'LineWidth', 1.5)
    % xlabel('time (s)')
    ylabel('inter-fly distance (mm)')
    % xlim(xlimit)
% FLY SPEED
subplot(r,c,2); hold on 
    sSpan = fps; %single second smoothing
    x = time(1:end-1);
    y1 = smooth(speed(:,1),sSpan,'moving');
    y2 = smooth(speed(:,2),sSpan,'moving');
    plot(x, y1, 'color', Color('dodgerblue'),'LineWidth', 2)
    plot(x, y2, 'color', Color('deeppink'),'LineWidth', 2)
    ylabel('fly speed (mm/s)')
    xlim(xlimit)

% FLY SPEED CORRELATION
subplot(r,c,3); hold on 
    tSpan = 3; % time window in seconds
    pSpan = tSpan * fps; % frames in the sliding window
    y = runningCorrelation(speed, pSpan);
    y = smooth(y,pSpan,'moving');
    offset = ceil(pSpan/2);
    x = time(offset:offset+length(y)-1)';
    plot(x, y, 'color', foreColor,'LineWidth', 2)
    xlabel('time (s)')
    ylabel('M & F speed correlation')
    h_line(0,'r','--',1)
    xlim(xlimit)

formatFig(fig, blkbnd,[r,c]);
subplot(r,c,1)
set(gca, 'xcolor', backColor)
subplot(r,c,2)
set(gca, 'xcolor', backColor)

save_figure(fig,[baseFolder 'Figures/speed and distance timecourse'], fig_type);




%%


dataFile = 'labels.v001.analysis.h5';
filePath = [baseFolder, dataFile];

occupancy_matrix = h5read(filePath,'/track_occupancy');
tracks = h5read(filePath,'/tracks');
node_names = h5read(filePath, '/node_names');
% edge_names = h5read(filePath, '/edge_names');
track_names = h5read(filePath, '/track_names');

data = [];
blkbnd = false;
fig_type = '-pdf';
[foreColor,backColor] = formattingColors(blkbnd); %get background colors

%% Create data labels and structure organization 
frames = 1:length(occupancy_matrix);
nframes = frames(end);

% tracks (frame, body points, XY, fly)
m = squeeze(tracks(:,:,:,2)); %male fly
f = squeeze(tracks(:,:,:,3)); %female fly

fps = 80;
pix2mm = 0.0305;

% time
vid_length = (nframes/fps); %time in seconds
time = linspace(0,vid_length,nframes);

% inter-fly-distance from the fly's center point
x1 = m(:,2,1);
y1 = m(:,2,2);
x2 = f(:,2,1);
y2 = f(:,2,2);
IFD = (sqrt((x1-x2).^2 + (y1-y2).^2)).*pix2mm;

% fly speed
speed = [];
D = (sqrt((x1(1:end-1)-x1(2:end)).^2 + (y1(1:end-1)-y1(2:end)).^2)).*pix2mm; % male speed
speed(:,1) = (D./(1/fps));
D = (sqrt((x2(1:end-1)-x2(2:end)).^2 + (y2(1:end-1)-y2(2:end)).^2)).*pix2mm; % male speed
speed(:,2) = (D./(1/fps));
clear x1 x2 y1 y2 D

initial_var = who;
initial_var{end+1} = 'initial_var';
initial_var{end+1} = 'position';
initial_var{end+1} = 'wing';
clearvars('-except',initial_var{:})

%% Determine the conversion between pixels to mm

% vidpath = "S:\Evyn\Courtship Videos\9.12.2024\Courtship 0001.avi";
% movieInfo = VideoReader(vidpath); %read in video
% demoImg = (read(movieInfo,1));
% well_loc = readPoints(demoImg,4); % click on center positions of the wells in the arena
% 
% % distance between well 1 and 3
% d = [];
% well_1 = well_loc(:,1);
% well_3 = well_loc(:,3);
% d = [d;sum((well_1-well_3).^2).^0.5];
% % distance between well 2 and 4
% well_1 = well_loc(:,2);
% well_3 = well_loc(:,4);
% d = [d;sum((well_1-well_3).^2).^0.5];
% 
% pixelsbetweenwells = mean(d); %pixels
% actualdistance =  31.2; %mm
% pix2mm = actualdistance/pixelsbetweenwells; % multiplier


%% Plot fly positions for a given frame
frame = 1;

% plot fly positions: 
fig = getfig('',1,[560 420]);
    hold on
    plotFlySkeleton(fig, m(frame, :,1),m(frame, :,2),Color('dodgerblue'),true);
    plotFlySkeleton(fig, f(frame, :,1),f(frame, :,2),Color('deeppink'),true);
    axis square equal
set(gca, 'xcolor', backColor,'ycolor',backColor)


%% Distance, speed, and speed correcation between flies over full video course
% use the center body point to determine between-fly distance

r = 3; c = 1;
xlimit = [0,120];

fig = getfig('',1,[1032 1042]);
% INTERFLY DISTANCE
subplot(r,c,1) 
    plot(time,IFD,'color', foreColor,'LineWidth', 2)
    % xlabel('time (s)')
    ylabel('inter-fly distance (mm)')
    xlim(xlimit)
% FLY SPEED
subplot(r,c,2); hold on 
    sSpan = fps; %single second smoothing
    x = time(1:end-1);
    y1 = smooth(speed(:,1),sSpan,'moving');
    y2 = smooth(speed(:,2),sSpan,'moving');
    plot(x, y1, 'color', Color('dodgerblue'),'LineWidth', 2)
    plot(x, y2, 'color', Color('deeppink'),'LineWidth', 2)
    ylabel('fly speed (mm/s)')
    xlim(xlimit)

% FLY SPEED CORRELATION
subplot(r,c,3); hold on 
    tSpan = 3; % time window in seconds
    pSpan = tSpan * fps; % frames in the sliding window
    y = runningCorrelation(speed, pSpan);
    y = smooth(y,pSpan,'moving');
    offset = ceil(pSpan/2);
    x = time(offset:offset+length(y)-1)';
    plot(x, y, 'color', foreColor,'LineWidth', 2)
    xlabel('time (s)')
    ylabel('M & F speed correlation')
    h_line(0,'r','--',1)
    xlim(xlimit)

formatFig(fig, blkbnd,[r,c]);
subplot(r,c,1)
set(gca, 'xcolor', backColor)
subplot(r,c,2)
set(gca, 'xcolor', backColor)

save_figure(fig,[baseFolder 'Figures/speed and distance timecourse'], fig_type);

%% Screen for frames with funky wing positions
clearvars('-except',initial_var{:})

% postions: 1 -head, 2-center, 3-abdo, 4-left wing, 5-right wing
position = [];
for sex = 1:2
    switch sex
        case 1
            sex_type = 'male';
            x = m(:,:,1);
            y = m(:,:,2);
        case 2
            sex_type = 'female';
            x = f(:,:,1);
            y = f(:,:,2);
    end
    % 1) Shift each frame to the origin
    x_offset = x(:,2);
    y_offset = y(:,2);
    X = x-x_offset;
    Y = y-y_offset;
    
    % 2) Rotate each set of points to head and center on the Y axis
    newX = nan([length(X),5]);
    newY = nan([length(Y),5]);
    for frame = 1:length(X)
        points = [X(frame,:)',Y(frame,:)'];
        rotatedPoints = rotateToVertical(points, 2,1,false);
        newX(frame,:) = rotatedPoints(:,1);
        newY(frame,:) = rotatedPoints(:,2);
    end
    
    % plot fly positions for all time: 
    SZ = 5; r = 1; c = 2;
    fig = getfig(sex_type,1,[1058 781]);
    subplot(r,c,1)
        hold on
        CList = {'teal', 'red', 'grey', 'blue', 'gold'};
        for i = 1:5
            scatter(newX(:,i),newY(:,i),SZ, Color(CList{i}),'filled')
        end
        axis square equal
        set(gca, 'xcolor', backColor, 'ycolor',backColor)
        title('Before Correction')
    
    % 3) eliminate points that are absolutely incorrect
    failstate = newY(:,1)<0; % head is below the center
    failstate = [failstate, newY(:,3:5)>0]; %abdomen or wings are above the center point (on head side)
    loc = any(failstate,2);
    disp([sex_type ' fly elimination points: ' num2str(sum(loc))])
    newX(loc,:) = nan;
    newY(loc,:) = nan;
    
    % 4) flip R and L wing for mirror image ones
    flipstate = newX(:,5)<newX(:,3) & newX(:,5)<newX(:,4); % right wing is to the left of the abdomen and L wing
    flipstate = [flipstate, newX(:,4)>newX(:,3) & newX(:,4)>newX(:,5)]; % left wing is to the right of the abdomen and R wing
    loc = any(flipstate,2);
    disp([sex_type ' fly flip wing points: ' num2str(sum(loc))])
    wing_L = newX(loc,5);
    wing_R = newX(loc,4);
    newX(loc,4:5) = [wing_L,wing_R];
    
    % 5) remove frames with both wings on one side of the abdoment
    flipstate = newX(:,5)<newX(:,3) & newX(:,4)<newX(:,3); % wings the left of the abdomen
    flipstate = [flipstate, newX(:,5)>newX(:,3) & newX(:,4)>newX(:,3)]; %wings to the right of the abdomen 
    loc = any(flipstate,2);
    disp([sex_type ' fly wing-abdo elimination points: ' num2str(sum(loc))])
    newX(loc,:) = nan;
    newY(loc,:) = nan;
    
    % plot fly positions after adjustments: 
    subplot(r,c,2)
        hold on
        CList = {'teal', 'red', 'grey', 'blue', 'gold'};
        for i = 1:5
            scatter(newX(:,i),newY(:,i),SZ, Color(CList{i}),'filled')
        end
        axis square equal
        set(gca, 'xcolor',backColor, 'ycolor',backColor)
        title('After Correction')
    % save the data
    position(sex).x = newX;
    position(sex).y = newY;
end

% plot male and female on the same scale
 fig = getfig(sex_type,1,[1058 781]);
 for sex = 1:2
    subplot(r,c,sex)
    hold on
        CList = {'teal', 'red', 'grey', 'blue', 'gold'};
        for i = 1:5
            scatter(position(sex).x(:,i),position(sex).y(:,i),SZ, Color(CList{i}),'filled')
        end
        axis square equal
        set(gca, 'xcolor',backColor, 'ycolor',backColor)
 end
fig = matchAxis(fig,true);

%% Calculate the angle between the fly wings
clearvars('-except',initial_var{:})
wing = [];
for sex = 1:2
    for w = 1:2
        switch w 
            case 1
                P1 = [position(sex).x(:,4),position(sex).y(:,4)]; % left
            case 2
                P1 = [position(sex).x(:,5),position(sex).y(:,5)]; % right
        end
        
        P2 = [position(sex).x(:,2),position(sex).y(:,2)]; % center
        P3 = [position(sex).x(:,3),position(sex).y(:,3)]; % abdomen
        
        % 1: Calculate vectors from P2 to P1 and from P2 to P3
        v1 = P1 - P2;  % Nx2 matrix, vector from P2 to P1 for each time step
        v2 = P3 - P2;  % Nx2 matrix, vector from P2 to P3 for each time step
        
        % 2: Calculate the dot product of v1 and v2 for each time step
        dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);
        
        % 3: Calculate the magnitudes of v1 and v2 for each time step
        mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2);
        mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2);
        
        % 4: Calculate the cosine of the angle for each time step
        cosTheta = dotProduct ./ (mag_v1 .* mag_v2);
        
        % 5: Compute the angle in radians, and then convert to degrees
        anglesRadians = acos(cosTheta);  % Nx1 vector of angles in radians
        anglesDegrees = rad2deg(anglesRadians);  % Convert to degrees
        
        wing(sex).angle(:,w) = anglesDegrees;
    end
end

% compare male L and R wing angles
fig = getfig('',true, [1032 300]); 
    hold on
    plot(time, wing(1).angle(:,1),'color', Color('blue'),'linewidth', 1) % left wing
    plot(time, wing(1).angle(:,2),'color', Color('gold'),'linewidth', 1) % right wing
    xlabel('time (s)')
    ylabel('wing angle (\circ)')
formatFig(fig, blkbnd);
set(gca, 'xcolor',foreColor)
save_figure(fig,[baseFolder 'Figures/male wing angle'],fig_type);

% compare male and female wing angles
fig = getfig('',true, [1032 300]); 
hold on
    plot(time, wing(1).angle(:,1),'color', Color('dodgerblue'),'linewidth', 1) % left wing
    plot(time, wing(1).angle(:,2),'color', Color('dodgerblue'),'linewidth', 1,'LineStyle','--') % right wing
    plot(time, wing(2).angle(:,1),'color', Color('deeppink'),'linewidth', 1) % left wing
    plot(time, wing(2).angle(:,2),'color', Color('deeppink'),'linewidth', 1,'LineStyle','--') % right wing
xlabel('time (s)')
ylabel('wing angle (\circ)')
formatFig(fig, blkbnd);
save_figure(fig,[baseFolder 'Figures/M and F wing angles'],fig_type);

%% Angle between the flies
clearvars('-except',initial_var{:})

% 1) Create a vector defined by the head and center of the male fly
P1 =[m(:,1,1), m(:,1,2)]; % male head x-y vector
P2 =[m(:,2,1), m(:,2,2)]; % male center x-y vector
v1 = P2 - P1;  % Vector for male fly

% 2) Create a vector defined by the head and center of the male fly
P3 =[f(:,1,1), f(:,1,2)]; % female head x-y vector
P4 =[f(:,2,1), f(:,2,2)]; % female center x-y vector
v2 = P3 - P1;  % Vector for female fly

% 3) Dot product of the two vectors
dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);    

% 4) Compute the magnitudes of the vectors
mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 

% 5) Calculate the cosine of the angle
cosTheta = dotProduct ./ (mag_v1 .* mag_v2);

% 6) Compute the angle in radians and convert to degrees
angleRadians = acos(cosTheta);  % Angle in radians
angleDegrees = rad2deg(angleRadians);  % Convert to degrees
data.btwnflyangle = angleDegrees;

% 7) Plot the angle between the two flies over time
fig = getfig('',true, [1032 300]); 
    hold on
    plot(time, angleDegrees,'color', foreColor,'linewidth', 1) % left wing
    xlabel('time (s)') 
    ylabel('angle between flies (\circ)')
formatFig(fig, blkbnd);
set(gca, 'xcolor', foreColor)
save_figure(fig,[baseFolder 'Figures/angle between flies'],'-png');





%% Distance, speed, and speed correlation between flies over full video course
% use the center body point to determine between-fly distance
clearvars('-except',initial_var{:})

fig_type = '-pdf';
blkbnd = false;
[foreColor,backColor] = formattingColors(blkbnd); %get background colors

frame = 5100;
sz = 50;
frame_skip = 5;
windowsize = 4; % seconds
roi = windowsize*80;
ROI = frame-roi:frame_skip:frame;

r = 5; c = 1;
% xlimit = [50,70];
xlimit = [time(ROI(1)),time(ROI(end))];

fig = getfig('',1,[1032 1042]);
% INTERFLY DISTANCE
subplot(r,c,1) 
    plot(time,IFD,'color', foreColor,'LineWidth', 2)
    % xlabel('time (s)')
    ylabel('inter-fly distance (mm)')
    xlim(xlimit)

% FLY SPEED
subplot(r,c,2); hold on 
    sSpan = fps; %single second smoothing
    x = time(1:end-1);
    y1 = smooth(speed(:,1),sSpan,'moving');
    y2 = smooth(speed(:,2),sSpan,'moving');
    plot(x, y1, 'color', Color('dodgerblue'),'LineWidth', 2)
    plot(x, y2, 'color', Color('deeppink'),'LineWidth', 2)
    ylabel('fly speed (mm/s)')
    xlim(xlimit)

% FLY SPEED CORRELATION
subplot(r,c,3); hold on 
    tSpan = 3; % time window in seconds
    pSpan = tSpan * fps; % frames in the sliding window
    y = runningCorrelation(speed, pSpan);
    y = smooth(y,pSpan,'moving');
    offset = ceil(pSpan/2);
    x = time(offset:offset+length(y)-1)';
    plot(x, y, 'color', foreColor,'LineWidth', 2)
    % xlabel('time (s)')
    ylabel('M & F speed correlation')
    h_line(0,'r','--',1)
    xlim(xlimit)

% FLY WING ANGLE
subplot(r,c,4); hold on 
    plot(time, wing(1).angle(:,1),'color', Color('dodgerblue'),'linewidth', 1) % left wing
    plot(time, wing(1).angle(:,2),'color', Color('dodgerblue'),'linewidth', 1,'LineStyle','--') % right wing
    plot(time, wing(2).angle(:,1),'color', Color('deeppink'),'linewidth', 1) % left wing
    plot(time, wing(2).angle(:,2),'color', Color('deeppink'),'linewidth', 1,'LineStyle','--') % right wing
    ylabel('wing angle (\circ)')
    xlim(xlimit)

% FLY WING ANGLE
subplot(r,c,5); hold on 
    plot(time, data.btwnflyangle,'color', foreColor,'linewidth', 1) % left wing
    xlabel('time (s)') 
    ylabel('angle between flies (\circ)')
     xlim(xlimit)

formatFig(fig, blkbnd,[r,c]);
for i = 1:4
    subplot(r,c,i)
    set(gca, 'xcolor', backColor)
end

save_figure(fig,[baseFolder 'Figures/full timecourse zoom in'], fig_type);


%%  Plot a frame with tracked points overlaid AND maybe a zoom in of the the behavior at that point?
clearvars('-except',initial_var{:})

vidpath = "S:\Evyn\DATA\Courtship Videos\Jaime Grant Figure\9.12.2024\Courtship 0001.avi";
movieInfo = VideoReader(vidpath); %read in video

% plot image with selected number of previously tracked points -- have a
% zoom in on the fly skeleton?

frame = 5100;
sz = 50;
frame_skip = 5;
demoImg = (read(movieInfo,frame));
img = imadjust(demoImg,[72/255, 180/255]);

[xlimits, ylimits] = deal([]);
fig = getfig; 
    % plot the image
    imshow(img)
    % plot the center points of the flies from the past 10 seconds
    windowsize = 4; % seconds
    roi = windowsize*80;
    ROI = frame-roi:frame_skip:frame;
    %MALE
    x1 = m(ROI, 1,1);
    y1 = m(ROI, 1,2);
    hold on
    scatter(x1,y1,sz,Color('dodgerblue'), "filled")
    %FEMALE
    x2 = f(ROI, 1,1);
    y2 = f(ROI, 1,2);
    scatter(x2,y2,sz,Color('deeppink'), "filled")

    lineROI = drawline(gca);

save_figure(fig,[baseFolder 'Figures/full frame image with flies ' num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type,false, false);

    % Zoom in on the flies
    xlimits(1) = min([x1;x2]);
    xlimits(2) = max([x1;x2]);
    ylimits(1) = min([y1;y2]);
    ylimits(2) = max([y1;y2]); 
    buff = 50;
    xlim([xlimits(1)-buff, xlimits(2)+buff])
    ylim([ylimits(1)-buff, ylimits(2)+buff])

    % overlay the current body position of the male fly
    x = m(frame, :,1);
    y = m(frame, :,2);
    scatter(x,y,sz,Color('black'), "filled")
    skeleton = [1,2; 2,3; 2,4; 2,5];
    for i = 1:size(skeleton,1)
        plot(x(skeleton(i,:)),y(skeleton(i,:)), 'color', Color('black'),'LineWidth', 1.5)
    end

    % overlay the current body position of the female fly
    x = f(frame, :,1);
    y = f(frame, :,2);
    scatter(x,y,sz,Color('black'), "filled")
    skeleton = [1,2; 2,3; 2,4; 2,5];
    for i = 1:size(skeleton,1)
        plot(x(skeleton(i,:)),y(skeleton(i,:)), 'color', Color('black'),'LineWidth', 1.5)
    end

save_figure(fig,[baseFolder 'Figures/zoom frame image with flies ' num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type);












































































