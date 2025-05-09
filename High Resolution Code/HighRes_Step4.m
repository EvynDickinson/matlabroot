%% 
% tracks (frame, body points, XY, fly)
% body points (head, center, abdomen, right wing, left wing)


%% Load tracking points
clear; clc;
baseFolder = getDataPath(6,0);

% Find files that can be run
[excelfile, Excel, xlFile] = load_HighResExperiments;

% Find trials that have nothing in their Proofed Column & have been tracked: 
loc = cellfun(@isnan,excelfile(2:end,Excel.tracked));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.expID, Excel.ramp, Excel.proofed]);
FileNames = format_eligible_files(eligible_files);

fileIdx = listdlg('ListString', FileNames,'ListSize',[350,450],'promptstring', 'Select trial to process');
if isempty(fileIdx)
    disp('No trials selected')
    return
end
% pull the list of dates and arenas to be 
dateDir = eligible_files(fileIdx,1);
trialDir = eligible_files(fileIdx,2); 
baseDir = [baseFolder, dateDir{:} '\', trialDir{:} '\'];

% initial_var = who;
% initial_var{end+1} = 'initial_var';

fileList = dir([baseDir '*alignment table.mat']);
load([baseDir, fileList(1).name]) % load the parameters and temp table
nvids = parameters.nVids; % number of videos
nBP = 5; % num of body parts

% load the video tracks
fileList = dir([baseDir '*.h5']);
data = [];
proofedFlag = [];
for vid = 1:nvids
    proofedFileCheck = dir([baseDir, 'compiled_video_', num2str(vid), '.*analysis.h5']);
    if ~isempty(proofedFileCheck) % file was edited 
        filePath = [baseDir proofedFileCheck(1).name];
        proofedFlag(end+1) = vid;
    else % no proofed file found
        filePath = [baseDir, 'compiled_video_' num2str(vid) '.avi.predictions.slp.h5'];
    end
    % load data from each video file: 
    data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
    data(vid).tracks = h5read(filePath,'/tracks');
    data(vid).node_names = h5read(filePath, '/node_names');
    data(vid).track_names = h5read(filePath, '/track_names');
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
        queststr = {['Vid ' num2str(vid)]; [num2str(overcount) ' points outside of primary tracks.']; 'Ignore and delete data?'};
        response = questdlg(queststr);
        if ~strcmp('Yes',response)
            disp(['Return to manual check of compiled video #' num2str(vid)])
            return
        end
    end
    data(vid).occupancy_matrix = occmat(:,1:2);
end

ntracks = [];
for vid = 1:parameters.nVids
    ntracks(vid) = size(data(vid).occupancy_matrix,2);
end

disp('cleaned data proofed')

%% Align the tracks from video to video: 
% display image of fly postion at the end of the first video and then the
% next chunk for the subsequent video to see if they align
buff = 4;
vid_align = false([1,nvids]);
vid_align(1) = true;
for vid = 1:nvids-1
    frame_end = length(data(vid).occupancy_matrix);
    eROI = frame_end-buff:frame_end;
    sROI = 1:1+buff;
    
    fig = getfig(['Vid ' num2str(vid)],0,[882 694]);
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

% Option to change order if a track was answered incorrectly: 
% TODO: (optional) incorporate a video verification of the alignment for
% the incorrectly assigned videos
switch questdlg('Were all the tracks correctly answered?')
    case 'Yes'
    case 'No'
        list_str = inputdlg('Write in the incorrect video numbers, separated by a space (e.g., ''22 45'')');
        incorrect_vids = str2double(strsplit(list_str{:},' '));
        for i = incorrect_vids
            vid_align(i) = ~vid_align(i);
        end
        disp('Video alignments updated')
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

%% Confirm and save fly identities
if ~exist([baseDir '/Figures'], 'dir')
    mkdir([baseDir '/Figures'])
end
pix2mm = 0.0289; % calculated on the new back setup 11/7/24
fps = parameters.FPS;

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


% Figure
r = 4;
c = 7;
sb(1).idx = 1:3; % wing angles over time
sb(2).idx = 8:10; % speed over time
sb(3).idx = [15:17, 22:24]; % wing angle histogram
sb(4).idx = [4:7,11:14, 18:21, 25:28]; % arena image

unclear = true;
time = T.time;

while unclear == true
fig = getfig('Fly identification',true);
    % Wing angles over time
    subplot(r, c, sb(1).idx)
        hold on
        plot(time, wing(1).angle(:,1),'color', Color('dodgerblue'),'linewidth', 1) % left wing
        plot(time, wing(1).angle(:,2),'color', Color('dodgerblue'),'linewidth', 1,'LineStyle','--') % right wing
        plot(time, wing(2).angle(:,1),'color', Color('deeppink'),'linewidth', 1) % left wing
        plot(time, wing(2).angle(:,2),'color', Color('deeppink'),'linewidth', 1,'LineStyle','--') % right wing
        ylabel('wing angle')
        xlabel('time (min)')
    % Speed over time
    subplot(r, c, sb(2).idx)
        hold on
        plot(time, speed(1).speed(:,1),'color', Color('dodgerblue'),'linewidth', 1)
        plot(time, speed(2).speed(:,1),'color', Color('deeppink'),'linewidth', 1)
        ylabel('speed (mm/s)')
        xlabel('time (min)')
    % Wing angle histogram
    subplot(r,c,sb(3).idx)
        hold on
        histogram(wing(1).angle(:)) %,'FaceColor', Color('dodgerblue'))
        histogram(wing(2).angle(:),'FaceColor',Color('deeppink'))
        xlim([0,90])
        ylabel('frequency')
        xlabel('wing angle')
    % Arena image
    subplot(r,c,sb(4).idx)
        vid = 1;
        vidname = [baseDir, 'compiled_video_', num2str(vid), '.avi'];
        vidh = VideoReader(vidname);
        frame = randi(vidh.NumFrames);
        img = read(vidh,frame);
        a = T.vidNums == vid;
        b = T.vidFrame == frame;
        exp_frame = T.frame(find(a & b));
        imshow(img)
        plotFlySkeleton(fig, D1.pos(exp_frame,:,1),D1.pos(exp_frame,:,2),Color('dodgerblue'),1); 
        plotFlySkeleton(fig, D2.pos(exp_frame,:,1),D2.pos(exp_frame,:,2),Color('pink'),1); 
        % Zoom in on pink fly
        xlimits = D2.pos(exp_frame,2,1);
        ylimits = D2.pos(exp_frame,2,2);
        buff = 250;
        xlim([xlimits-buff, xlimits+buff])
        ylim([ylimits-buff, ylimits+buff])
uiwait(fig)


% Save fly identities
response = questdlg('Which is pink?','Fly Sex','Female','Male','I dont know','Female');
switch response
    case 'Female'
        unclear = false;
        f = D2;
        f.wing = wing(2);
        m = D1;
        m.wing = wing(1);
        parameters.OGFlyOrder = 'MF';
        parameters.node_names = node_names;
        close all
    case 'Male'
        unclear = false;
        f = D1;
        f.wing = wing(1);
        m = D2;
        m.wing = wing(2);
        parameters.OGFlyOrder = 'FM';
        parameters.node_names = node_names;
        close all
    case 'I dont know'
        unclear = true;
end
end


%% Randomized track switch check point

response = questdlg('Double check for track switches?', '','Yes', 'No', 'Yes');
switch response
    case 'Yes'
    row = 2;
    col = 4;
    nfigs = ceil(parameters.nVids/col);
    
    frame = [];
    vid = 0;
    
    for figs = 1:nfigs
        fig = getfig;
        for pics = 1:col
            vid = vid + 1;
            if vid > parameters.nVids
                return
            end
            vidname = [baseDir, 'compiled_video_', num2str(vid), '.avi'];
            vidh = VideoReader(vidname);
            frame = randi(vidh.NumFrames); 
            % mnotrack = isnan(m.pos(exp_frame,2,1));
            % fnotrack = isnan(f.pos(exp_frame,2,1)); 
            % while mnotrack == true | fnotrack == true
            % frame = randi(vidh.NumFrames); 
            % mnotrack = isnan(m.pos(exp_frame,2,1));
            % fnotrack = isnan(f.pos(exp_frame,2,1)); 
            % end
            % if mnotrack == true | fnotrack == true
            %     frame = randi(vidh.NumFrames);
            % end
            img = read(vidh,frame);
            a = T.vidNums == vid;
            b = T.vidFrame == frame;
            exp_frame = T.frame(find(a & b));
            
            % Male fly
            subplot(row, col, pics)
                imshow(img)
                plotFlySkeleton(fig, m.pos(exp_frame,:,1),m.pos(exp_frame,:,2),Color('dodgerblue'),0); 
                % Zoom in on fly
                xlimits = m.pos(exp_frame,2,1);
                ylimits = m.pos(exp_frame,2,2);
                buff = 75;
                xlim([xlimits-buff, xlimits+buff])
                ylim([ylimits-buff, ylimits+buff])
                title(['video ', num2str(vid), ', frame ', num2str(exp_frame)])
            % Female fly
            subplot(row, col, (pics + col))
                imshow(img)
                plotFlySkeleton(fig, f.pos(exp_frame,:,1),f.pos(exp_frame,:,2),Color('pink'),0);
                % Zoom in on pink fly
                xlimits = f.pos(exp_frame,2,1);
                ylimits = f.pos(exp_frame,2,2);
                buff = 75;
                xlim([xlimits-buff, xlimits+buff])
                ylim([ylimits-buff, ylimits+buff])
        end
    end
    case 'No'
        return
end

%% Screen for frames with funky wing positions
disp_fig = true;
% postions: 1-head, 2-center, 3-abdomen, 4-left wing, 5-right wing
% Create variable to hold X and Y data for both male and female
DATA = [];
M = 1; % male fly index number
F = 2; % female fly index number
DATA(M).rawX = m.pos(:,:,1);
DATA(M).rawY = m.pos(:,:,2);
DATA(F).rawX = f.pos(:,:,1);
DATA(F).rawY = f.pos(:,:,2);

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
    DATA(sex).x = newX;
    DATA(sex).y = newY;
    
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
        save_figure(fig,[baseDir  'Figures/ ' sex_type ' wing position correction'], '-pdf');
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
                scatter(DATA(sex).x(:,i),DATA(sex).y(:,i),SZ, Color(CList{i}),'filled')
            end
            % format figure
            axis square equal
            set(gca, 'xcolor','none', 'ycolor','none')
            title(sex_type)
     end
    fig = matchAxis(fig,true);
    
    % Save figure
    save_figure(fig,[baseDir  'Figures/M vs F wing position scatter'], '-pdf');
end


%% Write to excel and Save base parameters
old_data = data;

data = DATA;

%  --------- save parameters ---------
trialID = [dateDir{1},'_',trialDir{1}];
c = [baseFolder,'Trial Data/',trialID,'/'];
if ~exist(c,'dir')
    mkdir(c)
end
save([c,'basic data.mat'],'T','m','f','parameters','data')

% --------- write parameters to excel sheet --------- 
[excelfile, Excel, xlFile] = load_HighResExperiments;

parameters.trialID = trialID;
% DATE LOCATION:
date_loc = (strcmp(excelfile(:,Excel.date),dateDir{1})); % find the rows in the excel sheet that match the current exp date
name_loc = (strcmp(excelfile(:,Excel.expID),expName)); % find name match
loc = find(date_loc & name_loc);
if isempty(loc)
    warndlg('Can''t find the row location for the file to add the trial ID information. Input manually.')
    disp(trialID)
else
    % write the trial ID into the excel sheet
    isExcelFileOpen(xlFile);
    writecell({trialID},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.trialID) num2str(loc)]);
    writecell({'Y'},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.basicfigs) num2str(loc)]);
end

disp(['Finished processing ' trialID])








