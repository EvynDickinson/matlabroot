%% 
% tracks (frame, body points, XY, fly)
% body points (head, center, abdomen, right wing, left wing)



%% PARAMETERS TO SHOW: 

% 2) each fly speed
% 3) body angle between flies
% 4) male fly wing angle

%% Load tracking points
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
proofedFlag = [];
% TODO: quick check to make sure this matches the expected number of videos...
for vid = 1:nvids
    filePath = [baseDir, 'compiled_video_' num2str(vid) '.avi.predictions.slp.h5'];
    proofedFilePath = [baseDir, 'compiled_video_' num2str(vid) '.avi.predictions.analysis.h5'];
    try data(vid).occupancy_matrix = h5read(proofedFilePath,'/track_occupancy');
        data(vid).tracks = h5read(proofedFilePath,'/tracks');
        data(vid).node_names = h5read(proofedFilePath, '/node_names');
        data(vid).track_names = h5read(proofedFilePath, '/track_names');
        proofedFlag(end+1) = vid;
    catch data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
        data(vid).tracks = h5read(filePath,'/tracks');
        data(vid).node_names = h5read(filePath, '/node_names');
        data(vid).track_names = h5read(filePath, '/track_names');
    end

    dataIn = true; %log successful data load
end

node_names = {'head', 'center', 'abdomen', 'right_wing', 'left_wing'}; % currently assuming these are stable across all videos
disp('data loaded')

%% Post proofing cleaning

for vid = 1:length(proofedFlag)
    i = proofedFlag(vid);
    x = squeeze(data(i).tracks(:,2,1,:));
    y = squeeze(data(i).tracks(:,2,2,:));
    buffersize = 3;
    
    % Quality control 1: screening for locations with more than 2 track points
    nanloc = isnan(x);
    a = sum(nanloc,2)>2;
    if sum(a)>0
        response = questdlg({[num2str(sum(a)) ' rows with 2+ points found.']; 'Replace with nans?'});
        if strcmpi(response,'Yes')
            x(a,:) = nan;
            y(a,:) = nan;
            data(i).tracks(a,:,:,:) = nan;
        else 
            return
        end  
    end
    
    % Quality control 2: condense multiple tracks into 2
    % Find where columns 1 and 2 have nans, but row sum = 2
    approval = [];
    nanloc = isnan(x);
    
    for flies = 1:2
        loc = find(isnan(x(:,flies)) & sum(~isnan(x),2)==2);
    
        if sum(loc)>0
            response = questdlg([num2str(length(loc)) ' instances found. What would you like to do?'],...
                '','Visually approve','Nan replace','Cancel','Visually approve');
            if strcmpi(response,'Visually approve') % option 1: visually approve
                for ii = loc
                    clist = {'blue','magenta'};
                    % Get x and y values for both flies 3 frames before and after frame of interest
                    X1 = squeeze(data(i).tracks(ii-buffersize:ii+buffersize,:,1,1));
                    Y1 = squeeze(data(i).tracks(ii-buffersize:ii+buffersize,:,2,1));
                    X2 = squeeze(data(i).tracks(ii-buffersize:ii+buffersize,:,1,2));
                    Y2 = squeeze(data(i).tracks(ii-buffersize:ii+buffersize,:,2,2));
                    
                    % Find track that has x or y value outside tracks 1 or 2
                    mysterypoint = find(~isnan(x(ii,3:end)))+2;
                    X3 = squeeze(data(i).tracks(ii,:,1,mysterypoint));
                    Y3 = squeeze(data(i).tracks(ii,:,2,mysterypoint));
        
                    fig = getfig(" ",1);
                    hold on
                    for frame = 1:size(X1,1)
                        plotFlySkeleton(fig, X1(frame,:),Y1(frame,:),Color('dodgerblue'),0); % plot Fly 1 positions before and after frame of interest
                        plotFlySkeleton(fig, X2(frame,:),Y2(frame,:),Color('pink'),0);  % plot Fly 2 positions before and after frame of interest
                    end
                    plotFlySkeleton(fig, X3,Y3,Color(clist{flies}),1); % plot mystery fly 
                    axis square
        
                    % Is it the correct fly identity? If yes, save which location,
                    % which track it came from, and which track it went to
                    approval(end+1,:) = [strcmp(questdlg('Correct fly identity?'),'Yes'),ii,mysterypoint,flies];
                    close(fig);
                end
            % elseif strcmpi(response,'Nan replace') % option 2: replace with nans/skip
            %     x(loc,:) = nan;
            % else
                return % option 3: stop and go to sleap
            end
        end
    end
    
    for ii = 1:size(approval,1)
        if approval(ii,1)==1
            disp([i,approval(ii,2)])
            % fill fly track with visually identified fly
            data(i).tracks(approval(ii,2),:,:,approval(ii,4)) = data(i).tracks(approval(ii,2),:,:,approval(ii,3));
            % replace extra fly track with nans
            data(i).tracks(approval(ii,2),:,:,approval(ii,3)) = nan;
        else
            % replace all fly tracks at that location with nans
            data(i).tracks(approval(ii,2),:,:,:) = nan;
        end
    end
end


% Confirm only 2 data tracks for each video
for vid = 1:nvids
    occmat = squeeze(~isnan(data(vid).tracks(:,2,1,:)));
    occsum = sum(occmat,1);
    overcount = sum(occsum(3:end));
    if overcount>0
        response = questdlg([num2str(overcount) ' remaining points outside of primary tracks. Ignore and delete data?']);
        if ~strcmp('Yes',response)
            return
        end
    end
    
    data(vid).occupancy_matrix = occmat(:,1:2);
end

ntracks = [];
for vid = 1:parameters.nVids
    ntracks(vid) = size(data(vid).occupancy_matrix,2);
end


%% Align the tracks from video to video: 
% display image of fly postion at the end of the first video and then the
% next chunk for the subsequent video to see if they align
buff = 4;
vid_align = false([1,nvids]);
vid_align(1) = true;
for vid = 1:nvids-1
    frame_end = size(data(vid).occupancy_matrix,1);
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
parameters.switch_state = switch_state;

% Align fly tracks
[d1,d2] = deal(struct);

for vid = 1:nvids
    if switch_state(vid)
        d1(vid).tracks = data(vid).tracks(:,:,:,2);
        d2(vid).tracks = data(vid).tracks(:,:,:,1);
    else
       d1(vid).tracks = data(vid).tracks(:,:,:,1);
       d2(vid).tracks = data(vid).tracks(:,:,:,2); 
    end
end

[D1.pos,D2.pos] = deal([]);

for vid = 1:nvids
    D1.pos = [D1.pos; d1(vid).tracks];
    D2.pos = [D2.pos; d2(vid).tracks];
end


%% Confirming fly identities

% Make male and female matrices

wing = [];
for sex = 1:2
    switch sex 
        case 1
            x = D1.pos(:,:,1);
            y = D1.pos(:,:,2);
        case 2
            x = D2.pos(:,:,1);
            y = D2.pos(:,:,2);
    end
    for w = 1:2
        switch w 
            case 1
                P1 = [x(:,4),y(:,4)]; % left
            case 2
                P1 = [x(:,5),y(:,5)]; % right
        end
        
        P2 = [x(:,2),y(:,2)]; % center
        P3 = [x(:,3),y(:,3)]; % abdomen
        
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

% Compare male and female wing angles
time = T.time;
fig = getfig('',true, [1032 300]); 
hold on
    plot(time, wing(1).angle(:,1),'color', Color('dodgerblue'),'linewidth', 1) % left wing
    plot(time, wing(1).angle(:,2),'color', Color('dodgerblue'),'linewidth', 1,'LineStyle','--') % right wing
    plot(time, wing(2).angle(:,1),'color', Color('deeppink'),'linewidth', 1) % left wing
    plot(time, wing(2).angle(:,2),'color', Color('deeppink'),'linewidth', 1,'LineStyle','--') % right wing
xlabel('time (s)')
ylabel('wing angle (\circ)')
formatFig(fig, true);
% save_figure(fig,[baseFolder 'Figures/M and F wing angles'],fig_type);

fig = getfig;
hold on
histogram(wing(1).angle(:))
histogram(wing(2).angle(:),'FaceColor',Color('deeppink'))

% Compare fly speed
speed = [];
for sex = 1:2
    switch sex 
        case 1
            x = D1.pos(:,2,1);
            y = D1.pos(:,2,2);
        case 2
            x = D2.pos(:,2,1);
            y = D2.pos(:,2,2);
    end
S = (sqrt((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2)).*pix2mm; % male speed
speed(sex).speed = [0;(S./(1/fps))];
end

time = T.time;
fig = getfig('',true, [1032 300]); 
hold on
    plot(time, speed(1).speed(:,1),'color', Color('dodgerblue'),'linewidth', 1)
    plot(time, speed(2).speed(:,1),'color', Color('deeppink'),'linewidth', 1)
xlabel('time (s)')
ylabel('speed (mm/s)')
formatFig(fig, true);

% View fly image
vidname = [baseDir, 'compiled_video_1.avi'];
vidh = VideoReader(vidname);
frame = 3500;

img = read(vidh,frame);

fig = getfig;
imshow(img)
hold on
plotFlySkeleton(fig, D1.pos(frame,:,1),D1.pos(frame,:,2),Color('dodgerblue'),1); 
plotFlySkeleton(fig, D2.pos(frame,:,1),D2.pos(frame,:,2),Color('pink'),1); 


%% Screen for frames with funky wing positions
clearvars('-except',initial_var{:})

% postions: 1-head, 2-center, 3-abdomen, 4-left wing, 5-right wing

% Create variable to hold X and Y data for both male and female
data = [];
data(M).rawX = m.pos(:,:,1);
data(M).rawY = m.pos(:,:,2);
data(F).rawX = f.pos(:,:,1);
data(F).rawY = f.pos(:,:,2);
initial_var{end+1} = 'data';

% For male and female, shift and rotate frame
for sex = 1:2
    % for each loop/sex, switch x and y to appropriate values
    switch sex
        case 1
            sex_type = 'male';
            x = m.pos(:,:,1);
            y = m.pos(:,:,2);
        case 2
            sex_type = 'female';
            x = f.pos(:,:,1);
            y = f.pos(:,:,2);
    end

    % 1) Shift each frame so that the center is at the origin (0,0)
    x_offset = x(:,2); % x values for center
    y_offset = y(:,2);
    % subtract the offset from the x values for each body point
    X = x-x_offset; % new shifted x data set
    Y = y-y_offset;
    
    % 2) Rotate each frame so that the head and center are aligned on the Y axis
    newX = nan([length(X),5]); % establish variable to hold new rotated x data set
    newY = nan([length(Y),5]);
    for frame = 1:length(X)
        % identify coordinates for each body point in that frame
        points = [X(frame,:)',Y(frame,:)'];
        % rotate points so head and center align with y axis
        rotatedPoints = rotateToVertical(points, 2,1,false);
        % load new rotated points into variable
        newX(frame,:) = rotatedPoints(:,1);
        newY(frame,:) = rotatedPoints(:,2);
    end
    
    % Plot body point positions for all frames 
    SZ = 5; r = 1; c = 2;
    if disp_fig
        fig = getfig(sex_type,1,[1058 781]);
        % before corrections plot
        subplot(r,c,1)
            hold on
            CList = {'teal', 'red', 'grey', 'blue', 'gold'};
            % plot each body point using new shifted and rotated data
            for i = 1:5
                scatter(newX(:,i),newY(:,i),SZ, Color(CList{i}),'filled')
            end
            % format figure
            axis square equal
            set(gca, 'xcolor', 'none', 'ycolor','none')
            title('Before Correction')
    end
    % 3) Eliminate frames that are absolutely incorrect
    failstate = newY(:,1)<0; % state where head is below the center
    failstate = [failstate, newY(:,3:5)>0]; % add state where abdomen or wings are above the center point (on head side)
    % identify which frames have a failstate (where there are nonzeros)
    loc = any(failstate,2); 
    % display number of frames to be eliminated for that fly
    disp([sex_type ' fly elimination points: ' num2str(sum(loc))])
    % replace elimination frames with nans
    newX(loc,:) = nan;
    newY(loc,:) = nan;
    
    % 4) Flip R and L wing for mirror image ones
    flipstate = newX(:,5)<newX(:,3) & newX(:,5)<newX(:,4); % state where right wing is to the left of the abdomen and L wing
    flipstate = [flipstate, newX(:,4)>newX(:,3) & newX(:,4)>newX(:,5)]; % add state where left wing is to the right of the abdomen and R wing
    % identify which frames have a failstate (where there are nonzeros)
    loc = any(flipstate,2);
    % display number of wing points to be flipped for that fly
    disp([sex_type ' fly flip wing points: ' num2str(sum(loc))])
    % flip x values for wing points needing flipped
    wing_L = newX(loc,5);
    wing_R = newX(loc,4);
    newX(loc,4:5) = [wing_L,wing_R];
    
    % 5) Eliminate frames with both wings on one side of the abdomen
    flipstate = newX(:,5)<newX(:,3) & newX(:,4)<newX(:,3); % state where wings are to the left of the abdomen
    flipstate = [flipstate, newX(:,5)>newX(:,3) & newX(:,4)>newX(:,3)]; % add state where wings are to the right of the abdomen 
    % identify which frames have a failstate (where there are nonzeros)
    loc = any(flipstate,2);
    % display number of frames to be eliminated for that fly
    disp([sex_type ' fly wing-abdo elimination points: ' num2str(sum(loc))])
    % replace elimination frames with nans
    newX(loc,:) = nan;
    newY(loc,:) = nan;
    % Save the cleaned, adjusted data
    data(sex).x = newX;
    data(sex).y = newY;
    
    if disp_fig
        % Plot fly positions after adjustments: 
        % after corrections plot
        subplot(r,c,2)
            hold on
            CList = {'teal', 'red', 'grey', 'blue', 'gold'};
            % plot each body point using cleaned data
            for i = 1:5
                scatter(newX(:,i),newY(:,i),SZ, Color(CList{i}),'filled')
            end
            % format figure
            axis square equal
            set(gca, 'xcolor','none', 'ycolor','none')
            title('After Correction')

        % Save figures
        save_figure(fig,[figDir  sex_type ' wing position correction'], fig_type);
    end
end


if disp_fig
    % Plot male and female body points on the same figure
     fig = getfig('Male vs Female',1,[1058 781]);
     for sex = 1:2
         switch sex
            case 1
                sex_type = 'Male';
            case 2
                sex_type = 'Female';
        end
        subplot(r,c,sex)
        hold on
            CList = {'teal', 'red', 'grey', 'blue', 'gold'};
            % plot each body point
            for i = 1:5
                scatter(data(sex).x(:,i),data(sex).y(:,i),SZ, Color(CList{i}),'filled')
            end
            % format figure
            axis square equal
            set(gca, 'xcolor','none', 'ycolor','none')
            title(sex_type)
     end
    fig = matchAxis(fig,true);
    
    % Save figure
    save_figure(fig,[figDir  'M vs F wing position scatter'], fig_type);
end

%% Save male and female identities + set base parameters

response = questdlg('What is pink?','Fly Sex','Female','Male','Cancel','Female');
switch response
    case 'Female'
        f = D2;
        f.wing = wing(2);
        m = D1;
        m.wing = wing(1);
        parameters.OGFlyOrder = 'MF';
    case 'Male'
        f = D1;
        f.wing = wing(1);
        m  = D2;
        m.wing = wing(2);
        parameters.OGFlyOrder = 'FM';
    case 'Cancel'
        return
end
parameters.node_names = node_names;

% --------- write parameters to excel sheet --------- 
[excelfile, Excel, xlFile] = load_CourtshipExperiments;
trialID = [dateDir{1},'_',trialDir{1}];
parameters.trialID = trialID;
% DATE LOCATION:
date_loc = (strcmp(excelfile(:,Excel.date),dateDir{1})); % find the rows in the excel sheet that match the current exp date
% RAMP LOCATION:
% find the ramp number using the diff between the folder and exp name
rampnum = trialDir{1}(end);
try rampnum = str2double(rampnum);
    expName_short = parameters.expID(1:end-1);
catch; rampnum = 1;
    expName_short = parameters.expID;
end
ramp_col = excelfile(2:end,Excel.ramp);
ramp_col = [ramp_col{:}];
ramploc = [false, (ramp_col==rampnum)];
% EXP ID LOCATION
exploc = (strcmp(excelfile(:,Excel.expID),expName_short));
% find location where these all align: 
xlRow = find(all([date_loc, exploc,  ramploc'],2));

% write the trial ID into the excel sheet
writecell({trialID},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.trialID) num2str(xlRow)]);

% save parameters
c = [baseFolder,'Trial Data/',trialID,'/'];
if ~exist(c,'dir')
    mkdir(c)
end
save([c,'basic data.mat'],'T','m','f','parameters')






