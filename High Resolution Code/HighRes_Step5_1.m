
% Load single trial data and run final processing and figures 

%% LOAD data
clear; clc;
warning off

% updates for updating the fly data with no swaps:
[excelfile, Excel, xlFile] = load_HighResExperiments;
% find trials that can be updated: 
done_loc = (strcmpi('Y', excelfile(:,Excel.swapcorrected)) | strcmpi('NA', excelfile(:,Excel.swapcorrected)));
loc_ready = strcmpi('Y', excelfile(:,Excel.groupready)) & ~done_loc;
loc = find(loc_ready);
% select trial to update: 
trial_options = excelfile(loc,Excel.trialID);
idx = listdlg('PromptString','Select experiments to update:','ListString',trial_options,...
                    'SelectionMode','single','ListSize',[350,300]);
excel_loc = loc(idx); % location in excel file of trials to process
path = getDataPath(6,2);
baseFolder = [path,'Trial Data/'];
trialDir = trial_options{idx};
baseDir = [baseFolder, trialDir '/']; % full folder directory for that trial
figDir = createFolder([baseDir,'Figures/']);

% % Manual selection (commented out for auto-revamp of all data on swaps)
% path = getDataPath(6,0);
% baseFolder = [path,'Trial Data/'];
% trialDir = selectFolder(baseFolder); 
% baseDir = [baseFolder, trialDir{:} '/']; % full folder directory for that trial


processed_path = [baseDir 'post-5.1.2 data.mat'];
if isfile(processed_path) % DATA ALREADY EXISTS
    if strcmp('Yes',questdlg('Processed data file found, load that?'))
        curr_baseFolder = baseFolder;
        curr_baseDir = baseDir;
        load(processed_path)
        baseFolder = curr_baseFolder;
        baseDir = curr_baseDir;
        conversion = getConversion; 
        pix2mm = conversion(4).pix2mm; 
        initial_var = add_var(initial_var, 'pix2mm');
        disp('Data loaded!')
        clearvars('-except',initial_var{:})
    else % cancel if data not currently needed
        return
    end
else % DATA DOES NOT YET EXIST
    if strcmp('Yes', questdlg('Run basic analysis now?'))
        load([baseDir, 'basic data.mat']) % load the parameters and temp table
        disp('data loaded')
        tic 
        % Experiment parameters
        nvids = parameters.nVids; % number of videos
        fps = parameters.FPS;
        time = T.time;
        blkbnd = true;
        fig_type = '-png';
        [foreColor,backColor] = formattingColors(blkbnd); % get background colors
        conversion = getConversion; 
        pix2mm = conversion(4).pix2mm; %updates 5.12.25
        
        % Create variable to store body points
        body = [];
        for i = 1:length(parameters.node_names)
            body.(parameters.node_names{i}) = i;
        end
        
        M = 1; % male fly index number
        F = 2; % female fly index number
        
        % Initial variables
        initial_var = who; % who = all variables created so far
        initial_var{end+1} = 'initial_var';
        initial_var{end+1} = 'well';
        initial_var{end+1} = 'path';
        
        disp_fig = false; % display baseline figures?
        initial_var{end+1} = 'disp_fig';
        toc
        disp('Data loaded, continuing evaluation')
    end
end

%% ANALYSIS: Extract calculated variables
clearvars('-except',initial_var{:})

% Inter-fly-distance from the fly's center point
x1 = m.pos(:,body.center,1); % x location for male center
y1 = m.pos(:,body.center,2);
x2 = f.pos(:,body.center,1); % x location for female center
y2 = f.pos(:,body.center,2);
% calculate interfly distance : distance between fly body centers
IFD = hypot((x2-x1),(y2-y1))./pix2mm;
% add interfly distance to the data table
T.IFD = IFD; 

% Fly speed (position one ahead-current position/time)
% Male Speed:
D = hypot(diff(x1), diff(y1)); % find the distance in pixels one frame to the next
D = D./pix2mm; % convert from pixels/frame to mm/frame
D = D.*fps; % convert to mm/sec
m.speed = [0; D];
% Female speed: 
D = hypot(diff(x2), diff(y2)); % find the distance in pixels one frame to the next
D = D./pix2mm; % convert from pixels/frame to mm/frame
D = D.*fps; % convert to mm/sec
f.speed = [0; D];

%% Find high-speed short-term swap frames and reverse them: 
clearvars('-except',initial_var{:})

speed_thresh = 35; %mm/s speed threshold
skip_threshold = 3; % how many frames for a confident swap pair in time
nNeighbors = 5; % how many preceeding frames to look at for close distance?

% Find possible swap locations: 
m_loc = m.speed>=speed_thresh;
f_loc = f.speed>=speed_thresh;

% speed above threshold for both male and female
double_speed_loc = m_loc & f_loc;
swap_ConfidentFrames = find(double_speed_loc);

% speed above threshold and nan for other fly
a = f_loc & isnan(m.speed);
b = m_loc & isnan(f.speed);
possible_swapFrames = find(a | b);

% combine all possible swap frames: 
allSwaps = sort([swap_ConfidentFrames; possible_swapFrames]);

% PAIR LIKELY SWAP FRAMES FOR EACH OF THE TRIALS
% Find differences between consecutive frames
frame_diffs = diff(allSwaps);
% Find where frames are close together (potential pairs)
is_close = frame_diffs <= skip_threshold;
% Extract pairs
pairs = [];
ii = 1; % index location within frame list
while ii <= length(is_close)
    if is_close(ii)
        % Found a pair: frames a(i) and a(i+1)
        pairs = [pairs; allSwaps(ii), allSwaps(ii+1)];
        
        % Skip the next position to avoid overlapping pairs
        % (e.g., if frames 1,2,3 are all close, we want pairs [1,2] and not [2,3])
        ii = ii + 2;
    else
        ii = ii + 1;
    end
end

tot_frames = size(T,1);
m_framesP = mean((sum(m_loc)/tot_frames)*100);
f_framesP = mean((sum(f_loc)/tot_frames)*100);
fprintf('\n%2.3g  percent of total frames are high speed\n', (m_framesP + f_framesP))

if ~isempty(pairs)
    % Second filter: body size and position relative to past trajectory
    mPos = squeeze(m.pos(:,body.center,:));
    fPos = squeeze(f.pos(:,body.center,:));
    
    % for each of the likely speed frame pairs
    % likely_pairs = pairs;
    rois = pairs(:,1) + ((-nNeighbors-1): -1);
    
    % check location alignment for each swap pair
    likelyswitch = false([size(pairs, 1), 1]);
    for ii = 1:size(pairs,1)
        
        % frame locations for pre and during swaps
        pre_roi = rois(ii,:);
        dur_roi = pairs(ii,1);
        
        % distance of each during swap roi to pre-swap roi relative
        x = mPos(dur_roi,1) - mPos(pre_roi,1); 
        y = mPos(dur_roi,2) - mPos(pre_roi,2);
        dM = mean(hypot(x,y),'omitnan');  % male track distance to assigned male (distance to male)
        x = mPos(dur_roi,1) - fPos(pre_roi,1); 
        y = mPos(dur_roi,2) - fPos(pre_roi,2);
        dF = mean(hypot(x,y),'omitnan');  % male track distance to assigned female (distance to female)
    
        % determine if likely switch
        likelyswitch(ii) = dM>=dF; % female would be closer than correct male
    end
    
    % extract the likely swap locations
    if any(likelyswitch)
        swap_pairs = pairs(likelyswitch,:); % pull only the likely pairs
    end

    % relative percentages
    swap_percentage =  mean((numel(swap_pairs) ./ tot_frames)*100,'omitnan');
    fprintf('\n%2.5g  percent of total frames are confident swaps \n',swap_percentage);
else 
    swap_pairs = [];
    disp('UPDATE: NO SWAPPED FRAMES FOUND')

    % determine if processing a full data set or just updating the data
    % with the new name...
    processed_file = [baseDir 'post-5.1.1 data.mat'];
    if exist(processed_file,"file")==2 % old processed data already exists -- change the name now: 
        newpath =  [baseDir 'post-5.1.2 data.mat'];

        status = movefile(processed_file, newpath);
        if status == 1
            disp('File renamed successfully');
        else
            disp(['Error: ' msg]);
        end
        % update that it was processed in the excel sheet:     
        isExcelFileOpen(xlFile);
        writecell({'NA'},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.swapcorrected) num2str(excel_loc)]);  
        % NA for no swaps, not applicible 
        disp('Data updated and name changed.')
        return

    end

end

%% FIX FLY ID SWAPS: (TODO: 2/26)
if ~isempty(swap_pairs)

nSwaps = size(swap_pairs,1);
frameBuffer = 3; % how many frames before/after swap to show
mColor = Color('vaporwaveblue');
fColor = Color('vaporwavepink');
[r, c] = subplot_numbers(nSwaps);
sz = 15; % node size for fly skeleton


% Preview the swaps: 
fig = getfig('', 1); 
for ii = 1:nSwaps
   subplot(r, c, ii); hold on
   % plot the swapped frames
   pre_roi = swap_pairs(ii,1)-frameBuffer:swap_pairs(ii,1)-1;
   post_roi = swap_pairs(ii,2) : swap_pairs(ii,2)+frameBuffer;
   OG_roi = [pre_roi, post_roi]; % combined before and after
   swap_roi = swap_pairs(ii,1):swap_pairs(ii,2)-1; % swap frames to color opposite

   % --- pull body locations for plotting ---

   % Male tracks: 
   mX_swap = data(F).rawX(swap_roi,:);
   mY_swap = data(F).rawY(swap_roi,:);
   mX = data(M).rawX(OG_roi,:);
   mY = data(M).rawY(OG_roi,:);
   plotFlySkeleton(fig, [mX; mX_swap], [mY; mY_swap], mColor, true, sz);
   % highlight corrected frames
   scatter(mX_swap(:), mY_swap(:), sz+10, foreColor, 'filled')

   % Female tracks: 
   fX = data(F).rawX(OG_roi,:);
   fY = data(F).rawY(OG_roi,:);
   fX_swap = data(M).rawX(swap_roi,:);
   fY_swap = data(M).rawY(swap_roi,:);
   plotFlySkeleton(fig, [fX; fX_swap], [fY; fY_swap], fColor, true, sz);
   % highlight corrected frames
   scatter(fX_swap(:), fY_swap(:), sz+10, foreColor, 'filled')
   
   % formatting and labels for mismatched alignments
   title(sprintf('%s | %i-%i', Alphabet(ii),swap_pairs(ii,1),swap_pairs(ii,2)-1))
end
formatFig(fig, blkbnd,[r,c]);
for ii = 1:nSwaps
     subplot(r, c, ii)
     set(gca, 'xcolor', 'none', 'ycolor', 'none')
     axis equal
end

% Manually approve the swaps: 
switch questdlg('Are all tracks in alignment?')
    case 'Yes'
        disp('swapping data information')
        % raw alignment...
        swap_roi = [];
        for ii = 1:nSwaps
            swap_roi = [swap_roi, swap_pairs(ii,1) : swap_pairs(ii,2)-1];
        end
        
        % Swap 'data' sections: 
        data_fields = {'rawX', 'rawY', 'x', 'y'};
        for ii = 1:length(data_fields)
            correct_M = data(F).(data_fields{ii})(swap_roi,:);
            correct_F = data(M).(data_fields{ii})(swap_roi,:);
            % swap data: 
            data(M).(data_fields{ii})(swap_roi,:) = correct_M;
            data(F).(data_fields{ii})(swap_roi,:) = correct_F;
        end
        
        % Swap 'f' and 'm' sections
        correct_M = f.pos(swap_roi,:,:);
        correct_F = m.pos(swap_roi,:,:);
        m.pos(swap_roi,:,:) = correct_M;
        f.pos(swap_roi,:,:) = correct_F;
        
        correct_M = f.wing.angle(swap_roi,:);
        correct_F = m.wing.angle(swap_roi,:);
        m.wing.angle(swap_roi,:) = correct_M;
        f.wing.angle(swap_roi,:) = correct_F;

        % recalculate speed: 
        % Inter-fly-distance from the fly's center point
        x1 = m.pos(:,body.center,1); % x location for male center
        y1 = m.pos(:,body.center,2);
        x2 = f.pos(:,body.center,1); % x location for female center
        y2 = f.pos(:,body.center,2);
        % Male Speed:
        D = hypot(diff(x1), diff(y1)); % find the distance in pixels one frame to the next
        D = D./pix2mm; % convert from pixels/frame to mm/frame
        D = D.*fps; % convert to mm/sec
        m.speed = [0; D];
        % Female speed: 
        D = hypot(diff(x2), diff(y2)); % find the distance in pixels one frame to the next
        D = D./pix2mm; % convert from pixels/frame to mm/frame
        D = D.*fps; % convert to mm/sec
        f.speed = [0; D];

    case 'No'
        % % give list of options that do not need to be swapped? Then show
        % % the original version of those frames...
        % idxList = Alphabet(1:nSwaps);
        % fix_idx = listdlg('PromptString', 'Select the incorrectly handled swap trials',...
        %                 'ListString',idxList,'SelectionMode','multiple');
        % 
        % % look at the original vs the corrected swap versions for the trials that were suspect: 
        % r = 1; c = 2;
        % for ii = 1:length(fix_idx)
        %     fig = getfig('', 1); 
        %            subplot(r, c, 1); 
        %            hold on
        %            % plot OG frames: 
        %            roi = swap_pairs(ii,1)-frameBuffer:swap_pairs(ii,2)+frameBuffer;
        %            mX = data(M).rawX(roi,:);
        %            mY = data(M).rawY(roi,:);
        %            fX = data(F).rawX(roi,:);
        %            fY = data(F).rawY(roi,:);
        %            plotFlySkeleton(fig, mX, mY, mColor, true, sz);
        %            plotFlySkeleton(fig, fX, fY, fColor, true, sz);
        % 
        %            % plot 'fixed' frames:
        %            % plot the swapped frames
        %            pre_roi = swap_pairs(ii,1)-frameBuffer:swap_pairs(ii,1)-1;
        %            post_roi = swap_pairs(ii,2) : swap_pairs(ii,2)+frameBuffer;
        %            OG_roi = [pre_roi, post_roi]; % combined before and after
        %            swap_roi = swap_pairs(ii,1):swap_pairs(ii,2)-1; % swap frames to color opposite
        % 
        %            % --- pull body locations for plotting ---
        % 
        %            % Male tracks: 
        %            mX = data(M).rawX(OG_roi,:);
        %            mY = data(M).rawY(OG_roi,:);
        %            % add middle frames: 
        %            mX = [mX; data(F).rawX(swap_roi,:)];
        %            mY = [mY; data(F).rawY(swap_roi,:)];
        % 
        %            % Female tracks: 
        %            fX = data(F).rawX(OG_roi,:);
        %            fY = data(F).rawY(OG_roi,:);
        %            % add middle frames: 
        %            fX = [fX; data(M).rawX(swap_roi,:)];
        %            fY = [fY; data(M).rawY(swap_roi,:)];
        % 
        %            % plot both sets of tracks: 
        %            plotFlySkeleton(fig, mX, mY, mColor, true, sz);
        %            plotFlySkeleton(fig, fX, fY, fColor, true, sz);
        % 
        %            % formatting and labels for mismatched alignments
        %            title([Alphabet(ii)])

        
        % TODO: tell Evyn if one of these pops up so we can make the code
        % to account for it 

        % (give number pair to the swaps, show all the swaps on a single image)
        % (question approve the swaps) 
        % (option to UNCORRECT a swap if not actually a swap)

        warndlg('Manually check and fix the misaligned tracks...')
        disp('Manual alignment of the tracks is still in progress....')
        return
    case {'Cancel', ''}
        warndlg('Manually check and fix the misaligned tracks...')
        disp('Manual alignment of the tracks is still in progress....')
end

end

%% ANALYSIS: Calculate M and F wing angles
clearvars('-except',initial_var{:})

% Calculate wing angles for male and female
wing = [];
for sex = 1:2
    % for each loop/wing, switch x and y to appropriate values
    for w = 1:2
        switch w 
            case 1
                P1 = [data(sex).x(:,body.right_wing),data(sex).y(:,body.right_wing)]; % 
            case 2
                P1 = [data(sex).x(:,body.left_wing),data(sex).y(:,body.left_wing)]; % right
        end
        
        P2 = [data(sex).x(:,body.center),data(sex).y(:,body.center)]; % center
        P3 = [data(sex).x(:,body.abdomen),data(sex).y(:,body.abdomen)]; % abdomen
        
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
        anglesDegrees = rad2deg(anglesRadians);  % convert to degrees
        
        % Save wing angle values as variable
        data(sex).wingangle(:,w) = anglesDegrees;
    end

     % identify where there are nans for wing angle for either wing
     loc = any(isnan(data(sex).wingangle),2);
     % create wing spread variable by summing wingangles
     data(sex).wingspread = sum(data(sex).wingangle,2);
     % replace wing spread values with nans where wing angles have nans
     data(sex).wingspread(loc,:) = nan;
end

%% ANALYSIS: Calculate the body angle between the flies
clearvars('-except',initial_var{:})

x = 1;
y = 2;

% positions of M and F head and center
P1 = [m.pos(:,body.center,x),m.pos(:,body.center,y)]; % male center
P2 = [m.pos(:,body.head,x),m.pos(:,body.head,y)]; % male head
P3 = [f.pos(:,body.center,x),f.pos(:,body.center,y)]; % female center
P4 = [f.pos(:,body.head,x),f.pos(:,body.head,y)]; % female head

% 1: Calculate body vectors
v1 = P1 - P2;  % Nx2 matrix, vector for male body vector
v2 = P3 - P4;  % Nx2 matrix, vector for female body vector

% 2: Calculate the dot product of v1 and v2 for each time step
dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);

% 3) Compute the magnitudes of the vectors
mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 

% 4) Calculate the cosine of the angle
cosTheta = dotProduct ./ (mag_v1 .* mag_v2);

% 5) Compute the angle in radians and convert to degrees
angleRadians = acos(cosTheta);  % angle in radians
angleDegrees = rad2deg(angleRadians);  % convert to degrees
data(1).mfbodyangle = angleDegrees;
data(2).mfbodyangle = angleDegrees;

%% ANALYSIS: Determine male position relative to female

clearvars('-except',initial_var{:})
% center all points to female fly
% align all points to female fly heading (center = 0)

response = questdlg('Visualize male positions too?','','Yes','No','No');

data(M).color = Color('dodgerblue');
data(F).color = Color('deeppink');

% ------------------------- 1) Shift each frame to the origin for female fly ------------------------- 
x_offset = data(F).rawX(:,body.center); % x values for female center
y_offset = data(F).rawY(:,body.center);
mx = data(M).rawX-x_offset; % subtract the offset from the x values for each body point
my = data(M).rawY-y_offset;
fx = data(F).rawX-x_offset;
fy = data(F).rawY-y_offset;

% ------------------------- 2) Rotate each set of points to head and center on the Y axis ------------------------- 
[mX,mY,fX,fY] = deal(nan(size(mx)));
for frame = 1:length(mX)
    % identify coordinates for each body point in that frame - female only
    fpoints = [fx(frame,:)',fy(frame,:)']; 
    % rotate points so head and center align with y axis
    [rotatedPoints, R] = rotateToVertical(fpoints, body.center,body.head,false);
    % load new rotated points into variable
    fX(frame,:) = rotatedPoints(:,1);
    fY(frame,:) = rotatedPoints(:,2);
    % identify coordinates for each body point in that frame - male only
    mpoints = [mx(frame,:)',my(frame,:)'];
    % rotate points in relation to female (using female R)
    mpoints = (R * mpoints')';
    % load new rotated points into variable
    mX(frame,:) = mpoints(:,1);
    mY(frame,:) = mpoints(:,2);
end
initial_var{end+1} = 'mX';
initial_var{end+1} = 'mY';
initial_var{end+1} = 'fX';
initial_var{end+1} = 'fY';

% Calulate the angle between fly bodies (both - and +)
theta = data(M).mfbodyangle; 
test = theta<90; % when is male less than 90 deg from female

% ------------------------- 3) Establish likely and unlikely courtship positions ------------------------- 

% Determine fly body length (~2.5mm)
BL = 2.5*pix2mm; % in pixels

% Body position location rules
r.a = mX(:,body.center) >= 0; % everything to right of y axis
r.b = mX(:,body.center) <= 0; % everything to left of y axis
r.c = mY(:,body.center) >= 0; % everything above x axis
r.d = mY(:,body.center) <= 0; % everything below x axis
q1 = r.b & r.c; % quadrant one occupancy
q2 = r.a & r.c; % quadrant two occupancy
q3 = r.d & r.b; % quadrant three occupancy
q4 = r.a & r.d; % quadrant four occupancy

% Heading direction rules
r.e = mY(:,body.head) >= mY(:,body.center); % heading direction facing north
r.f = mY(:,body.head) <= mY(:,body.center); % heading direction facing south
r.g = mX(:,body.head) >= mX(:,body.center); % heading direction facing east
r.h= mX(:,body.head) <= mX(:,body.center); % heading direction facing west
ne = r.e & r.g; % north east direction
nw = r.e & r.h; % north west direction
se = r.f & r.g; % south east direction
sw = r.f & r.h; % south west direction


% Likely rules
BLidx = abs(mX(:,body.center)) <= (5*BL) & abs(mY(:,body.center)) <= (5*BL); % limits likely roi's to 5 body lengths
% 1-4)
roi1 = q1 & se & BLidx; % quadrant 1, southeast
roi2 = q2 & sw & BLidx; % quadrant 2, southwest
roi3 = q3 & ne & BLidx; % quadrant 3, northeast
roi4 = q4 & nw & BLidx; % quadrant 4, northwest

% Gray/maybe rules
r.i = theta < 90;
r.j = theta > 90;
deg = [45, 20, 10];

    % 5) quadrant 1, southwest
    roi5 = [];
    for i = 1:3 % each angle
        % body center is within appropriate BL, body angle is between [deg] and 180, and inside quadrant 1
        dummy = abs(mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q1 & sw;
        roi5(:,i) = dummy;
    end
    roi5 = any(roi5,2);
    
    % 6) quadrant 2, southeast
    roi6 = [];
    for i = 1:3
        dummy = abs(mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q2 & se;
        roi6(:,i) = dummy;
    end
    roi6 = any(roi6,2);
    
    % 7) quadrant 3, northwest
    roi7 = [];
    for i = 1:3
        dummy = abs(mX(:,body.center)) <= i*BL & theta < deg(i) & q3 & nw;
        roi7(:,i) = dummy;
    end
    roi7 = any(roi7,2);
    
    % 8) quadrant 4, northeast
    roi8 = [];
    for i = 1:3
        dummy = abs(mX(:,body.center)) <= i*BL & theta < deg(i) & q4 & ne;
        roi8(:,i) = dummy;
    end
    roi8 = any(roi8,2);
    
    % 9) quadrant 1, northeast
    roi9 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q1 & ne;
        roi9(:,i) = dummy;
    end
    roi9 = any(roi9,2);
    
    % 10) quadrant 2, northwest
    roi10 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q2 & nw;
        roi10(:,i) = dummy;
    end
    roi10 = any(roi10,2);
    
    % 11) quadrant 3, southeast
    roi11 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q3 & se;
        roi11(:,i) = dummy;
    end
    roi11 = any(roi11,2);
    
    % 12) quadrant 4, southwest
    roi12 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q4 & sw;
        roi12(:,i) = dummy;
    end
    roi12 = any(roi12,2);

% ROI's that are likely courtship positions
likely = roi1 | roi2 | roi3 | roi4;
% ROI's that are maybe courtship positions
gray = roi5 | roi6 | roi7 | roi8; % exceptions moving towards X axis
gray2 = roi9 | roi10 | roi11 | roi12; %exceptions moving towards Y axis

% Saving ROIs
position = [];
position.L1 = roi1;
position.L2 = roi2;
position.L3 = roi3;
position.L4 = roi4;
position.GX1 = roi5;
position.GX2 = roi6;
position.GX3 = roi7;
position.GX4 = roi8;
position.GY1 = roi9;
position.GY2 = roi10;
position.GY3 = roi11;
position.GY4 = roi12;


% Create variable for likely and maybe courtship positions
loc = likely | gray | gray2;
T.courtposition = loc;  % yes = green, no = red on graph vv
likelyidx = find(loc);
unlikelyidx = find(~loc);
position.likely = likely;
position.GX = gray;
position.GY = gray2;
position.all_likely = loc;
position.unlikely = ~loc;
initial_var{end+1} = 'position';

% ------------------------------------------ 4) Visualizations ------------------------------------------ 

switch response
    case 'Yes'

    % 4.1) Demo selected angles in male positions relative to female fly
        zoom = [-250,250];
        % adjust skip size based on total number of points
        displayNum = 500; 
        if sum(test)<displayNum
            displayNum = sum(test);
        end
        plotLoc = randi(sum(test), [displayNum,1]);
        
        fig = getfig('Random selection of all male positions relative to female',1,[1075 871]);
        hold on
        
        % plot male coordinates for head and body within test constraints
        x = mX(test,[body.head,body.center]);
        y = mY(test,[body.head,body.center]);
        % only plot every [skip value] points to reduce volume
        x = x(plotLoc,:);
        y = y(plotLoc,:);
        plot(x',y','color',data(M).color)
        scatter(x(:,1),y(:,1),15,data(M).color,"filled","^") % arrow head on male
        
        % plot female coordiates for head, body, and abdomen
        x = fX(plotLoc,[body.head,body.center,body.abdomen]);
        y = fY(plotLoc,[body.head,body.center,body.abdomen]);
        plot(x',y','color',data(F).color, 'LineWidth', 2)
        
        % format figure
        axis  equal square
        h_line(0,'gray',':',2)
        v_line(0,'grey',':',2)
        xlim(zoom)
        ylim(zoom)
        formatFig(fig,blkbnd);
        set(gca,'XColor','none','YColor','none')
        
        save_figure(fig,[figDir ' all male fly body positions relative to female'],'-png');    

        %---------------------------------------------------------------------------------------------------------

        % 4.2) Visualize likely and unlikely courtship positions 
        zoom = [-250,250];
        skip = 20;

        
        likelyidx = likelyidx(1:skip:end);
        unlikelyidx = unlikelyidx(1:skip:end);
        allidx = [likelyidx; unlikelyidx];
        
        fig = getfig('Likely and unlikely courtship positions',1,[1075 871]);
        hold on
        
        % % plot all fly positions unlikely courtship male fly body positions
        % kolor = Color('grey');
        % x = mX(allidx,[body.head,body.center]);
        % y = mY(allidx,[body.head,body.center]);
        % plot(x',y','color',kolor)
        % scatter(x(:,1),y(:,1),15,kolor,"filled","^")
        % % axis equal square
        % formatFig(fig,blkbnd)
        % set(gca,'XColor','none','YColor','none')
        
        % plot the unlikely courtship male fly body positions
        kolor = Color('red');
        x = mX(unlikelyidx,[body.head,body.center]);
        y = mY(unlikelyidx,[body.head,body.center]);
        plot(x',y','color',kolor)
        scatter(x(:,1),y(:,1),15,kolor,"filled","^")
        
        % plot the likely courtship male fly body positions
        kolor = Color('limegreen');
        x = mX(likelyidx,[body.head,body.center]);
        y = mY(likelyidx,[body.head,body.center]);
        plot(x',y','color',kolor)
        scatter(x(:,1),y(:,1),15,kolor,"filled","^")
        
        % plot female head, center, and abdoment
        x = fX(1:skip:end,[body.head,body.center,body.abdomen]);
        y = fY(1:skip:end,[body.head,body.center,body.abdomen]);
        plot(x',y','color',foreColor, 'LineWidth', 2)
        % xlim([-2000,1500]); ylim([-2000,2000])
        
        % format figure
        axis  equal square
        h_line(0,'gray',':',2)
        v_line(0,'grey',':',2)
        xlim(zoom)
        ylim(zoom)
        formatFig(fig,blkbnd);
        set(gca,'XColor','none','YColor','none')
        
        formatFig(fig,blkbnd);
        rectangle('Position',[zoom(1),zoom(1) sum(abs(zoom)) sum(abs(zoom))],'edgecolor',foreColor,'linewidth', 1)
        
        save_figure(fig,[figDir 'likely and unlikely male body pos. relative to female'],'-png',0,1,'-r100');

    % case 'No'
    %     return
end

%% ANALYSIS: Identify food well and calulate distance to food

% determine if the well outlines already exist
well_file = [baseDir 'well locations.mat'];
if ~exist(well_file, 'file')    
    % pull up picture to find food well
    vidpath = [path, parameters.date, '/', parameters.videoName, '/compiled_video_1.avi'];
    % vidpath = '/Volumes/OnTheGoData/Courtship Videos/09.26.2024/Berlin_courtship_F_LRR_caviar_ramp/compiled_video_1.avi';
    movieInfo = VideoReader(vidpath); %read in video
    demoImg = (read(movieInfo,T.vidFrame(1)));
    img = imadjust(demoImg,[72/255, 215/255]); % adjust the contrast of the image
    
    % save food well location
    txt = {'12', '3', '6', '9'};
    well = struct; % initialize the new well structure
    fig = getfig('');
        imshow(img)
        for i = 1:4 
            h = warndlg(['Outline the ' txt{i} ' oclock well']);
            uiwait(h)
            roi = drawcircle; % manually add in the circle over the food well
            well.radius(i) = roi.Radius;
            well.center(i,:) = roi.Center;
        end
    h = warndlg('Close this when finished with the final well location');
    uiwait(h)    
    save_figure(fig,[figDir 'well outlines'],'-pdf',0,1, '-r100');

    % select the well with the food
    fig = getfig('');
        imshow(img); hold on
        for i = 1:4
            viscircles(well.center(i,:), well.radius(i),'color',Color('white'));
            % drawcircle('Center',[well.center(i,:)],'Radius',well.radius(i),'StripeColor','red','Color',Color('grey'));
        end
        title('Click the food well','FontSize',18)
        [xi, yi] = crosshairs(1,{'black','black','yellow','yellow'});      % get a point
        % find the well with the selected data point
        % distance between each well center and the selected point
        [~, index] = min((sqrt((xi-well.center(:,1)).^2 + (yi-well.center(:,2)).^2))./pix2mm);
        
        % plot the selected point
        scatter(xi,yi,60,"white",'filled')
        scatter(xi,yi,30,"red",'filled')
        % change the color of the well outline to highlight selection
        viscircles(well.center(index,:), well.radius(index),'color',Color('red'));
    if strcmp(questdlg('okay food well selection?'),'Yes')
        well.food_idx = index;
        well.food = well.center(index,:);
        close(fig)
    else
        warndlg('Rerun the food well identification')
        return
    end
    % SAVE DATA?
    if strcmp(questdlg('save well locations?'),'Yes')
        save(well_file,'well')
    else 
    end
else
    load(well_file,'well')
    disp('Loaded prior well locations')
    % load('/Users/evyndickinson/Documents/Jeanne Lab/Courtship Videos/Trial Data/01.14.2025_Berlin_courtship_F_LRR_caviar_ramp1/well locations.mat','well');
end

well.R = 2.5; % this is a fixed value based on the arena size 

% Center of the arena
WC = well.center';
N = [];
x1 = WC(1,1:2:4);
y1 = WC(2,1:2:4);
x2 = WC(1,2:2:4);
y2 = WC(2,2:2:4);
[xi,yi] = polyxpoly(x1,y1,x2,y2);
well.center(5,:) = [xi,yi];

% calculate distance to food (from fly head)
x1 = m.pos(:,body.head,1); % x location for male center
y1 = m.pos(:,body.head,2);
x2 = f.pos(:,body.head,1);
y2 = f.pos(:,body.head,2);
c1 = well.center(1);
c2 = well.center(2);

m.dist2food = (sqrt((x1-c1).^2 + (y1-c2).^2))./pix2mm;
f.dist2food = (sqrt((x2-c1).^2 + (y2-c2).^2))./pix2mm;
T.dist2food = [m.dist2food, f.dist2food];

%% ANALYSIS: Determine fly occupancy in ring, quadrant, food circle  
clearvars('-except',initial_var{:})

% Max Distance from Center of Arena
R = conversion(4).R; % outer distance limits
innerR = conversion(4).circle75; % radius of a circle occupying 75% of the arena
maxR = conversion(4).circle10; % radius of a circle occupying 10% of the arena

% find distance from center for each fly center point: 
xi = well.center(5,1);
yi = well.center(5,2);

% eccentricity, outter ring occupancy of the flies
for sex = 1:2
    x_loc = data(sex).rawX(:,body.center); % x position of the fly
    y_loc = data(sex).rawY(:,body.center); % y position of the fly
    
    %distance from center of arena
    D = sqrt(((x_loc-xi).^2 + (y_loc-yi).^2))./pix2mm; 
    data(sex).eccentricity = D;

    % outter ring occupancy
    data(sex).OutterRing = D<=R & D>=innerR; % find the locations that are between edge and inner R

    % quadrant occupancy 
    center = well.center(5,:);
    % r = data(exp).data(trial).data.r;
    foodWell = well.food; %data(exp).T.foodLoc(trial);

    % Adjust the X and Y coordinates relative to the new center
    adjustedX = x_loc - center(1);
    adjustedY = y_loc - center(2);
    
    % Initialize matrix to hold quadrant classification (same size as input matrices)
    quadrantMatrix = zeros(size(x_loc));
    
    % Define quadrant masks based on the new center
    Q = [];
    Q(1).Mask = (adjustedY > adjustedX) & (adjustedY <= -adjustedX);  % Top
    Q(2).Mask = (adjustedY <= adjustedX) & (adjustedY <= -adjustedX); % Bottom
    Q(3).Mask = (adjustedY <= adjustedX) & (adjustedY > -adjustedX);  % Left
    Q(4).Mask = (adjustedY > adjustedX) & (adjustedY > -adjustedX);   % Right
    
    % Determine the well locations (which determine quadrant assignment):
    adjusted_wx = foodWell(1) - center(1);
    adjusted_wy = foodWell(2) - center(2);
    
    idx_loc = false(1,4);
    % Find the food quadrant (find location with the food well coordinates included)
    idx_loc(1) = (adjusted_wy > adjusted_wx) & (adjusted_wy <= -adjusted_wx);  % top
    idx_loc(2) = (adjusted_wy <= adjusted_wx) & (adjusted_wy <= -adjusted_wx); % right
    idx_loc(3) = (adjusted_wy <= adjusted_wx) & (adjusted_wy > -adjusted_wx);  % bottom
    idx_loc(4) = (adjusted_wy > adjusted_wx) & (adjusted_wy > -adjusted_wx);   % left
    quad_loc = find(idx_loc);

    fly_loc = ~isnan(x_loc); %gives logical for all fly positions in the position matrix
    foodQuad = Q(quad_loc).Mask & fly_loc; % flies in food quad
    
    data(sex).foodQuad = foodQuad;
    
    % 10% Food Circle
    D = sqrt(((x_loc-foodWell(1)).^2 + (y_loc-foodWell(2)).^2)); % distance from center of food well
    D = D./pix2mm;
    data(sex).foodcircle = D<=maxR ; % flies within 10% circle around the food

end

% Fly on the food: 
T.FlyOnFood = T.dist2food<=well.R; % fly head must be within the food circle
disp('Flies on food: M & F')
sum(T.FlyOnFood)

%% ANALYSIS: Determine temperature bins and directions
% TODO (2/26) update this to work for temp holds and trials that don't have
% a temperature protocol name:
clearvars('-except',initial_var{:})

% Check if there is a paratameter protocol: 
try disp(parameters.protocol)
catch
    h = warndlg('no temperature protocol found -- please manually select an option:');
    uiwait(h);
    [excelfile, Excel, xlFile] = load_HighResExperiments;
    tempProtocolList = unique(excelfile(2:end,Excel.protocol));
    startIdx = find(strcmp(tempProtocolList,'courtship_F_LRR_25-17'));
    idx = listdlg('PromptString','Select temp protocol for this trial', 'SelectionMode', 'single', 'InitialValue',startIdx,...
        'ListString',tempProtocolList,'ListSize',[200 200]);
    if isempty(idx)
        return
    end
    parameters.protocol = tempProtocolList{idx};
end
tp = getTempParamsHighRes(parameters.protocol);

initial_var{end+1} = 'tRate';
temp_file = [baseDir 'temp regions.mat'];
if ~exist(temp_file, 'file') 
    % manually select the different time points
    x = []; 
    fig = getfig; 
    plot(time, T.temperature,'color', 'k')
    title(tp.labelstr)
    if tp.ntrans>0 % for all dynamic trials
        for i = 1:tp.ntrans
            [xi, ~] = crosshairs(1,{'black','black','red','red'});  % get a point
            [~,xidx] = min(abs(time-xi));
            x(i) = xidx(1);
            v_line(time(x(i)),'r','--',1)
        end
        % find the x-time value for each time period
        tRate = tp.tRate;
        tRate(1).idx = [1, x(1)-1];
        for i = 2:tp.ntrans
            tRate(i).idx = [x(i-1), x(i)-1];
        end
        tRate(end).idx = [ x(i),T.frame(end)];
    else % static trials
        tRate = tp.tRate;
        tRate(1).idx = [1, T.frame(end)];
        y_avg = mean(T.temperature);
        ylim([y_avg-2, y_avg+2]) 
    end
    
    % plot the regions onto the graph
    hold on
    ylims = ylim;
    for i = 1:length(tRate)
        roi = time(tRate(i).idx);
        rectangle('Position', [roi(1), ylims(1), diff(roi),diff(ylims)], 'FaceColor', tRate(i).color);
    end
    plot(time, smooth(T.temperature,fps*10,'moving'),'color', 'k','linewidth', 5)
    xlabel('time (min)')
    ylabel('temperature (\circC)')
    formatFig(fig,false);
    title('temperature conditions')
    
     if strcmp(questdlg('save temperature regions?'),'Yes')
         save_figure(fig, [figDir, 'temperature regions'],'-pdf',1);
         save(temp_file,'tRate')
         disp('saved temperature regions')
     else 
         close(fig)
    end
else
    load(temp_file,'tRate')
    % load('/Users/evyndickinson/Documents/Jeanne Lab/Courtship Videos/Trial Data/01.14.2025_Berlin_courtship_F_LRR_caviar_ramp1/temp regions.mat','tRate');
    disp('loaded temperature regions')
end

% Set up logicals for each of the temperature regions
MT = false(size(time));
%cooling
idx = find(strcmp('cooling',{tRate(:).name}));
T.cooling = MT;
if ~isempty(idx)
    for i = 1:length(idx)
        T.cooling(tRate(idx(i)).idx(1):tRate(idx(i)).idx(2)) = true;
    end
end
%warming
idx = find(strcmp('warming',{tRate(:).name}));
T.warming = MT;
if ~isempty(idx)
    for i = 1:length(idx)
        T.warming(tRate(idx(i)).idx(1):tRate(idx(i)).idx(2)) = true;
    end
end
%holds
idx = [find(strcmp('start hold',{tRate(:).name})) , find(strcmp('end hold',{tRate(:).name})), find(strcmp('hold',{tRate(:).name}))];
T.hold = MT;
if ~isempty(idx)
    for i = 1:length(idx)
        T.hold(tRate(idx(i)).idx(1):tRate(idx(i)).idx(2)) = true;
    end
end

%% ANALYSIS: Calculate male wing extension
clearvars('-except',initial_var{:})

Lwing = [];
Rwing = [];
% Determine which positions require which wing to be extended 
L_items = {'L1', 'L4', 'GX1', 'GX4', 'GY2', 'GY3'};
R_items = {'L2', 'L3', 'GX2', 'GX3', 'GY1', 'GY4'};
% Identify if male is in an appropriate position for each wing direction across each item
for i = 1:length(L_items)
    Lwing = [Lwing, position.(L_items{i})];
    Rwing = [Rwing, position.(R_items{i})];
end
% Condense to identify if male is in any of the appropriate positions
Lwing = any(Lwing,2);
Rwing = any(Rwing,2);

% Pull wing angles equal or greater than extension minimum for L and R
wa_cutoff = 50; % minimum wing extension angle for courtship
wing_ext = (Lwing & (m.wing.angle(:,1) >= wa_cutoff)) | (Rwing & (m.wing.angle(:,2) >= wa_cutoff)); % wing must be facing the female fly
% Each value subtracted by the value before it (1 = ext starts, -1 = ext stops, 0 = no state change)
a = diff(wing_ext); 
% Add the first extension value to the list to account for the starting condition
b = [wing_ext(1); a]; 
% Locations in wing_ext where extension period starts/end
ext_start = find(b == 1); 
ext_stop = find(b == -1);
% If wing ext doesn't stop by end, add stop location at end of ext_stop (loc = length of experiment value)
if wing_ext(end)
    ext_stop(end + 1) = length(time);
end
% Calculate the length of each wing ext bout
ext_dur = ext_stop - ext_start;
% Find where wing ext lasts longer than 1sec
dur_loc = find(ext_dur > fps);

% Create new courtship matrix with only true wing ext for bouts longer than 1sec
mt = false(size(time));
for i = 1:length(dur_loc)
    ii = dur_loc(i);
    mt(ext_start(ii):ext_stop(ii)) = true;
end
T.wing_ext = mt;
T.wing_ext_all = wing_ext;

% % TODO WORKING HERE 
% wing_ext_filled = imclose(wing_ext, ones(max_gap + 1, 1));


%% ANALYSIS: Chase identification
% < 120 deg area behind female x
% facing female x
% 7mm between M center and F center x
% female speed > 0 x
% if all are 1 then courtship
% diff to see if chasing lasts > 2sec x
clearvars('-except',initial_var{:})

x = 1;
y = 2;

% positions of M head and F head and center
P1 = [m.pos(:,body.head,x),m.pos(:,body.head,y)]; % male head
P2 = [f.pos(:,body.center,x),f.pos(:,body.center,y)]; % female center
P3 = [f.pos(:,body.head,x),f.pos(:,body.head,y)]; % female head

% 1) Calculate body vectors
v1 = P3 - P1;  % Nx2 matrix, vector for female head to male head
v2 = P3 - P2;  % Nx2 matrix, vector for female head to female center

% 2) Calculate the dot product of v1 and v2 for each time step
dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);

% 3) Compute the magnitudes of the vectors
mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 

% 4) Calculate the cosine of the angle
cosTheta = dotProduct ./ (mag_v1 .* mag_v2);

% 5) Compute the angle in radians and convert to degrees
angleRadians = acos(cosTheta);  % angle in radians
angleDegrees = rad2deg(angleRadians);  % convert to degrees
mfpos_angle = angleDegrees;

% Identify when male position angle is less than 60 degrees from female
pos_angle = abs(mfpos_angle) <= 60;

% Identify when male is facing female
facing = [];
m_items = {'L3', 'L4', 'GX3', 'GX4', 'GY3', 'GY4'};
% Identify if male is in an appropriate position for each wing direction across each item
for i = 1:length(m_items)
    facing = [facing, position.(m_items{i})];
end
facing = any(facing,2);

% Identify when male is behind female AND facing her
mbehindf = (facing & pos_angle);

% Identify when male is within 7mm of female
close_dist = T.IFD <= 7; % mm

% Identify when female is moving
fmoving = f.speed >= 0.1; % min speed up for debate

% Identify when male is moving
mmoving = m.speed >= 0.01; % min speed diff than F in order to include true chase bouts

chase = (mbehindf & close_dist & fmoving & mmoving); % requirements for 'chase all'
a = diff(chase); 
% Add the first chase value to the list to account for the starting condition
b = [chase(1); a]; 
% Locations in chase where chasing period starts/end
ch_start = find(b == 1); 
ch_stop = find(b == -1);
% If chasing doesn't stop by end, add stop location at end of ch_stop (loc = length of experiment value)
if chase(end)
    ch_stop(end + 1) = length(time);
end
% Calculate the length of each chasing bout
ch_dur = ch_stop - ch_start;
% Find where chasing lasts longer than 2sec
dur_loc = find(ch_dur > (2*fps));

% Create new courtship matrix with only true chasing bouts longer than 2sec
mt = false(size(time));
if isempty(dur_loc)
    m.chaseroi = [];
else
    for i = 1:length(dur_loc)
        ii = dur_loc(i);
        mt(ch_start(ii):ch_stop(ii)) = true;
        m.chaseroi(i,:) = [ch_start(ii), ch_stop(ii)];
    end
end
T.court_chase = mt; % time restriction 2 seconds
T.chase_all = chase; % NO time limit

%% ANALYSIS: Circling behavior
% M head within 3mm? [head_dist]
% M facing female [position.likely]
% M velocity constant within 1sec [const_var]
% F not moving
% if all are 1 then courtship
% diff to see if circling lasts > 1sec

clearvars('-except',initial_var{:})

% Distance between male head and female
x1 = m.pos(:,body.head,1); % x location for male head
y1 = m.pos(:,body.head,2);
x2 = f.pos(:,body.center,1); % x location for female center
y2 = f.pos(:,body.center,2);
% Calculate male head - female center distance
d = ((sqrt((x1-x2).^2 + (y1-y2).^2))./pix2mm);
head_dist = d <= 3; % mm
y = nan(size(head_dist));
y(head_dist) = 1;

% (for fun) how much time did the male fly spend really close to the female fly?
percent_time_within_3mm = (sum(head_dist)/length(time))*100;
disp(['The male fly spent ' num2str(percent_time_within_3mm) '% of the experiment close to the female'])

% position.all_likely

% female fly is not moving
f_speed_cut = (f.speed<=0.2);

% Determine if semiconstant male velocity
% 0.5sec smooth, then 2 sec std
binwidth = fps; %  1 second bins to look for avg speed
n = length(m.speed) - binwidth; % number of iterations
offset = ceil(binwidth * 0.5);
smoothspeed = smooth(m.speed,floor(0.5*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.5*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.3*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.2*fps),'moving');
sp_var = nan(size(time)); % initialize the speed variability variable
for i = 1:n % loop through each bin
    roi = i:i+binwidth; % frames from i to one second later
    sp_var(i+offset) = std(smoothspeed(roi));
end

testlim = 0.7;
const_var = sp_var>testlim;
ok_var = m.speed;
ok_var(const_var) = nan;% all data points with acceptable variance

% -- View the threshold for consistent speed --
% r = 2;
% c = 1;
% 
% fig = figure('Position',[1 10 1478 845]);
% subplot(r,c,1); hold on
%     plot(time,m.speed,'color', Color('red'))
%     plot(time, ok_var, 'color', foreColor)
% subplot(r,c,2)
%     plot(time,sp_var,'color', foreColor)
% % formatting
% xlimits = [1,2];
% subplot(r,c,1)
%     xlim(xlimits)
%     ylabel('male fly speed (mm/s)')
% subplot(r,c,2)
%     xlim(xlimits)
%     h_line(testlim, 'r','-',1)
%     ylabel('speed variance')
% formatFig(fig,blkbnd,[r,c])

% Full selection criteria: 
V = position.likely & head_dist & const_var & f_speed_cut;

%  --- Find periods longer than 1 second ---
a = diff(V); % when does the speed switch between stability and instability 
% Add the first chase value to the list to account for the starting condition
b = [V(1); a]; 
% Find when stability periods starts/end
v_start = find(b == 1); 
v_stop = find(b == -1);
% If speed stability doesn't stop by end, add stop location at end of v_stop (loc = length of experiment value)
if V(end)
    v_stop(end + 1) = length(time);
end
% Calculate the length of each speed bout
v_dur = v_stop - v_start;
% Find where speed stability lasts longer than 1sec
dur_loc = find(v_dur > fps);
% save the constant speed regions into a new metric
constant_velocity = false(size(time));
for i = 1:length(dur_loc)
    idx = dur_loc(i);
    constant_velocity(v_start(idx): v_stop(idx)) = true;
end

T.circling_all = V; % when the male fly is circling the female (no time restriction)
T.circling_1sec = constant_velocity; % when circling is longer than 1 second

% COURTSHIP INDEX
CI = any([T.court_chase,T.wing_ext,T.circling_1sec],2); % courtship index
T.CI = CI;

%% ANALYSIS: Fly turning 

clearvars('-except',initial_var{:})
for sex = 1:2
    % extract the x and y head and center positions to calculcate the slope of the body line
    x = data(sex).rawX(:,body.head:body.center);
    y = data(sex).rawY(:,body.head:body.center);
    
    % zero the flys center (actually unneccessary) 
    X = x - x(:,2);
    Y = y - y(:,2);
    % slope for each point in time
    slope = (Y(:,2)-Y(:,1))./(X(:,2)-X(:,1));
    m1 = slope(1:end-1);   % 'past' heading
    m2 = slope(2:end);      % 'current' heading
    
    % calculate the angle between the two slopes 
    theta = atan((m1-m2)./(1+(m1.*m2)));
    theta = rad2deg(theta);
    data(sex).turning = [nan; theta].*(fps);
end

%% ANALYSIS: Extract sleep data
clearvars('-except',initial_var{:})

bout = 5*60*parameters.FPS; %(5-min period minimum for sleep threshold)
dummy = [];

% Extract sleep bouts from position data
for sex = 1:2
    switch sex
        case 1
            x = m.pos(:,body.center,1); % male center
        case 2
            x = f.pos(:,body.center,1); % female center
    end
    % Calculate difference between all x values
    x_diff = diff(x); 
    % Identify when position is not changing
    u = abs(x_diff)<= 1;
    % Each value subtracted by the value before it (1 = ext starts, -1 = ext stops, 0 = no state change)
    a = diff(u);
    % Add the first position value to the list to account for the starting condition
    b = [u(1); a]; 
    % Frames where 'position-no-change' period starts/end
    slp_start = find(b == 1); 
    slp_stop = find(b == -1);
    % If sleep doesn't stop by end, add stop location at end of slp_stop
    if u(end)
        slp_stop(end + 1) = length(time);
    end
    % Calculate the length of each 'position-no-change' bout
    slp_dur = slp_stop - slp_start;
    % Find where bout lasts longer than 5min (when is sleep)
    slp_loc = find(slp_dur > bout);

    % Create dummy matrix with only true sleep bouts
    mt = false(size(time));
    for i = 1:length(slp_loc)
        ii = slp_loc(i);
        mt(slp_start(ii):slp_stop(ii)) = true;
    end
    dummy(sex).sleep = mt;
end

% Save sleep data
m.sleep = dummy(1).sleep;
f.sleep = dummy(2).sleep;

g = find(m.sleep);
h = find(f.sleep);
if isempty(g) && isempty(h)
    disp('no male or female sleep found')
end

%% ANALYSIS: Behavior probability map
clearvars('-except',initial_var{:})

% 1 = on food
% 2 = sleeping
% 3 = edge-occupancy
% 4 = courtship
b_list = {'food', 'sleeping', 'edge','m_courtship'};
nstates = length(b_list);
ntime = length(time);
time_buff = 10; % number of frames that can be 'skipped' before the state is considered switched

for sex = 1:2 % male and female
   
    tic
    % Initialize behavior states
    behavior = nan(ntime, 4);
    behavior(T.FlyOnFood(:, sex), 1) = 1;
    
    if sex == 1
        behavior(m.sleep, 2) = 2;
    else
        behavior(f.sleep, 2) = 2;
    end
    
    behavior(data(sex).OutterRing, 3) = 3;
    behavior(T.CI, 4) = 4;

    % Sleep takes priority when overlapping with other states
    sleep_override = ~isnan(behavior(:, 2)) & any(~isnan(behavior(:, [1, 3, 4])), 2);
    behavior(sleep_override, [1, 3, 4]) = nan;

    % Assign a single state per time point
    beh = nan(size(time));
    for i = [4, 3, 1, 2] % Sleep is assigned last to take priority
        beh(~isnan(behavior(:, i))) = i;
    end

    % Find state start indices
    state_starts = find(~isnan(beh));

    % Preallocate a table for state transitions
    max_transitions = numel(state_starts); % Overestimate, we'll trim later
    ST = table(nan(max_transitions, 1), nan(max_transitions, 1), nan(max_transitions, 1), ...
        nan(max_transitions, 1), nan(max_transitions, 1), nan(max_transitions, 1), nan(max_transitions, 1), ...
        'VariableNames', {'state1', 'state1_start', 'state1_end', 'state2', 'state2_start', 'state2_end', 'transition'});

    % Identify state transitions
    idx = 1;
    i = 1;
    while i < numel(state_starts)
        curr_state = beh(state_starts(i));
        state1_start = state_starts(i);
        state1_end = state1_start;

        % Move forward until we find a new state or exceed time_buff
        for j = i + 1:numel(state_starts)
            t = state_starts(j);
            next_state = beh(t);
            gap = t - state1_end;

            if curr_state == next_state && gap <= time_buff
                state1_end = t; % Extend the current state
            else
                % Store state transition
                ST.state1(idx) = curr_state;
                ST.state1_start(idx) = state1_start;
                ST.state1_end(idx) = state1_end;
                ST.state2(idx) = next_state;
                ST.state2_start(idx) = t;
                ST.transition(idx) = str2double([num2str(curr_state), num2str(next_state)]);

                % Update previous state's end time
                if idx > 1
                    ST.state2_end(idx - 1) = state1_end;
                end

                % Prepare for next transition
                idx = idx + 1;
                break;
            end
        end

        i = j; % Move to next segment
    end

    % Trim empty rows
    ST(idx:end, :) = [];
   
    % Make the probability map as a heatmap with the transition probability as a shaded square with a number
    % then it will also be easy to make a 'difference' map between each of the temperature periods
    
    transition_list = nan([nstates^2,1]);
    idx = 1;
    for i = 1:nstates
        for j = 1:nstates
            transition_list(idx) = str2double([num2str(i) num2str(j)]); 
            idx = idx + 1;
        end
    end
    
    transitions = [];
    for i = 1:length(transition_list)
        transitions(i) = sum(ST.transition==transition_list(i));
    end
    transitions = reshape(transitions,[nstates, nstates])';

    % Save some of the data to the data structure for future use:
    data(sex).states.ST = ST;
    data(sex).states.b_list = b_list;
    data(sex).states.nstates = nstates;
    data(sex).states.beh = beh;
    data(sex).states.behavior = behavior;
    data(sex).states.transitions = transitions;
    data(sex).states.transition_list = transition_list;
    toc
end

% % How long to each behavior state after each temp change and what order?
% for sex = 1:2 % for each of the two flies 
%     freq = [];
%     for r = 1:size(tRate,2) % for each of the different temperature regimes
%         roi = tRate(r).idx; % get the temperature regime frame indeces 
%         for i = 1:data(sex).states.nstates % for each type of behavior | state
%             freq(r,i).name = data(sex).states.b_list{i};
%             total_frames = diff(roi);
% 
%             % -- make a list of all the times for this behavior (so we can do a frequency graph etc) 
%             loc = data(sex).states.ST.state1==i;
%             a = find(loc);
%             if isempty(a)
%                 [freq(r,i).n_states, freq(r,i).n_instances, freq(r,i).perc_time_in_state] = deal(0);
%                 [freq(r,i).time_to_first_instance, freq(r,i).temp_at_first_instance] = deal(nan);
%                 continue
%             end
%             first_frame = data(sex).states.ST.state1_start(a(1));
%             freq(r,i).n_states = sum(loc); % number of times the behavior occured (not the duration of them)
% 
%             % -- total time in the state and the percent of the temp regime in that state
%             freq(r,i).n_instances = sum(data(sex).states.beh==i); % number of time points with the observed behavior
%             freq(r,i).perc_time_in_state = (freq(r,i).n_instances/total_frames)*100;
% 
%             % -- pull time of the first one
%             freq(r,i).time_to_first_instance = (first_frame - roi(1))/parameters.FPS; % in seconds
% 
%             % -- pull the temp for the first one
%             freq(r,i).temp_at_first_instance = mean(T.temperature(first_frame:first_frame+3));
%         end
%     end
%     data(sex).states.freq = freq;
% end

% How long to each behavior state after each temp change and what order?
for sex = 1:2  % for each of the two flies
    nstates = data(sex).states.nstates;
    nregimes = numel(tRate);

    % Preallocate freq struct
    freq(nregimes, nstates) = struct(...
        'name', '', ...
        'n_states', 0, ...
        'n_instances', 0, ...
        'perc_time_in_state', 0, ...
        'time_to_first_instance', nan, ...
        'temp_at_first_instance', nan);

    for r = 1:nregimes  % for each temperature regime
        roi = tRate(r).idx;  % get the temperature regime frame indices
        total_frames = diff(roi);

        for i = 1:nstates  % for each type of behavior/state
            freq(r, i).name = data(sex).states.b_list{i};
            loc = (data(sex).states.ST.state1 == i);
            
            if ~any(loc)
                continue  % Skip if behavior never occurred
            end
            
            % Find first occurrence
            first_idx = find(loc, 1, 'first');
            first_frame = data(sex).states.ST.state1_start(first_idx);

            % Number of occurrences
            freq(r, i).n_states = sum(loc);

            % Total time in state
            beh_match = (data(sex).states.beh == i);
            freq(r, i).n_instances = sum(beh_match);
            freq(r, i).perc_time_in_state = (freq(r, i).n_instances / total_frames) * 100;

            % Time to first instance in seconds
            freq(r, i).time_to_first_instance = (first_frame - roi(1)) / parameters.FPS;

            % Temperature at first occurrence (mean over 3 frames)
            freq(r, i).temp_at_first_instance = mean(T.temperature(first_frame:first_frame+3));
        end
    end
    
    % Store results in data structure
    data(sex).states.freq = freq;
end


% FIGURES: number and distribution of behaviors that were observed
% quick figure on the total set of observed behaviors
clearvars('-except',initial_var{:})

sexes = {'male', 'female'};
nstates = data(1).states.nstates;
b_list = strrep(data(1).states.b_list,'_',' ');
fig = getfig('',1,[ 937 406]);
for sex = 1:2
    subplot(1,2,sex)
    bar(1:nstates, sum(~isnan(data(sex).states.behavior)),'FaceColor',data(sex).color)
    set(gca, "XTickLabel", b_list)
    xlabel('behavior state')
    ylabel('full trial count')
end
formatFig(fig,false,[1,2]);
save_figure(fig,[figDir 'Behavior state bar graph full ramp  M and F'],fig_type);

% time course of the behavior states
fig = getfig('',1,[ 937 406]);
for sex = 1:2
    subplot(1,2,sex);
    hold on
    % h = rectangle('Position', [roi(1), ylims(1), diff(roi),diff(ylims)], 'FaceColor', tRate(i).color);

    for i = 1:4
        y = data(sex).states.behavior(:,i);
        scatter(time, y, 35, data(sex).color, 'filled')
    end
    xlabel('time (min)')
    ylabel('behavior state')
    ylim([0.5,nstates+0.5])
    set(gca, "YTick",1:nstates, 'YTickLabel',  b_list, 'ydir', 'reverse')
end
formatFig(fig,false, [1,2]);
save_figure(fig,[figDir 'Behavior state timecourse M and F'],fig_type);

% FIGURE: transition matrix heatmap of the states
fig = getfig('',1);
for sex = 1:2
    subplot(1,2,sex);
    heatmap(data(sex).states.transitions)
    xlabel('State 2')
    ylabel('State 1')
    set(gca, 'XDisplayLabels',b_list,'YDisplayLabels',b_list)
    title(sexes{sex})
end
save_figure(fig,[figDir 'Behavior state transitions heatmap M and F'],fig_type);

% Figures: Histogram of the different state transitions & total count of the behaviors (~time)
fig = getfig('',1,[1064 473]);
for sex = 1:2
    subplot(1,2,sex)
    bin_edges = min(data(sex).states.ST.transition)-1:max(data(sex).states.ST.transition)+1;
    histogram(data(sex).states.ST.transition,bin_edges,'FaceColor',data(sex).color)
    xlabel('transition')
    ylabel('count (#)')
end
fig = formatFig(fig,false,[1,2]);
save_figure(fig,[figDir 'Behavior state transitions histogram M and F'],fig_type);

% Figure:  pie chart of different behaviors across the full experiment
fig = getfig('',1);
for sex = 1:2
    subplot(1,2,sex)
    y = data(sex).states.beh;
    x = sum(isnan(y));
    for i = 1:data(sex).states.nstates
        x = [x, sum(y==i)];
    end
    xplode = true(size(x));
    xplode(1) = false; % do not explode out the unknown states/positions
    pie(x,xplode)
    title({sexes{sex}; '   '})
end
save_figure(fig,[figDir 'Behavior state transitions pie chart M and F'],fig_type);

% when (if) does each state arise in each temp regime? At what temp is the
% behavior observed?

%%  Save the 'analyzed' data package:
clearvars('-except',initial_var{:})

switch questdlg('Save processed data?')
    case 'Yes'
        save([baseDir 'post-5.1.2 data.mat'],'-v7.3')
        % update that it was processed in the excel sheet:     
        isExcelFileOpen(xlFile);
        writecell({'Y'},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.swapcorrected) num2str(excel_loc)]);  
        disp('Saved data file')
    case 'No'
        return
    case 'Cancel'
        return
end


% parameters.trialID = trialID;
% % DATE LOCATION:
% date_loc = (strcmp(excelfile(:,Excel.date),dateDir{1})); % find the rows in the excel sheet that match the current exp date
% name_loc = (strcmp(excelfile(:,Excel.expID),expName)); % find name match
% loc = find(date_loc & name_loc);
% if isempty(loc)
%     warndlg('Can''t find the row location for the file to add the trial ID information. Input manually.')
%     disp(trialID)
% else
%     % write the trial ID into the excel sheet
%     isExcelFileOpen(xlFile);
%     writecell({trialID},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.trialID) num2str(loc)]);
%     writecell({'Y'},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.basicfigs) num2str(loc)]);
% end



% NOTE: 2.27.26
% this code has been update to include a swaped-frames catch for single and
% double frame identity swap situations. Moving forward, data will be saved
% as 'post-5.1.2.mat' to indicate the updated version of the data

% NOTE: 5.12.25 
% this code has been updated to include the new distance measures for the
% plate size and thus moving forward, data will be saved as
% 'post-5.1.1.mat' rather than just 'post-5.1.mat' data. Do not use data
% from the original set, as it has an incorrect distance scaling factor. 



%% (DONE -- but save for any arena changes) Determine the conversion between pixels to mm for any new video configuration
% vidpath = "S:\Evyn\DATA\Courtship Videos\02.05.2025\Berlin_courtship_F_LRR_caviar_ramp1\compiled_video_1.avi";
% % vidpath = "S:\Evyn\DATA\Courtship Videos\09.26.2024\Berlin_courtship_F_LRR_caviar_ramp\compiled_video_2.avi";
% movieInfo = VideoReader(vidpath); %read in video
% demoImg = (read(movieInfo,1));
% well_loc = readPoints(demoImg,4); % click on center positions of the wells in the arena

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
% actualdistance =  36.2; %mm
% pix2mm = actualdistance/pixelsbetweenwells; % multiplier
% 
% R = 36.1;
% WC = well_loc;
% x1 = WC(1,1:2:4);
% y1 = WC(2,1:2:4);
% x2 = WC(1,2:2:4);
% y2 = WC(2,2:2:4);
% [xi,yi] = polyxpoly(x1,y1,x2,y2);
% centre = [xi,yi];
% 
% fig = figure;
%     imshow(demoImg); hold on
%     viscircles(centre,R*pix2mm,'color', 'w','LineWidth',0.25);
%     % find actual plate size: 
%     % roi = drawcircle; % manually add in the circle over the food well
%     % newR = roi.Radius/pix2mm;
%     newR = 33.6;
%     viscircles(centre,newR*pix2mm,'color', 'w','LineWidth',0.25);

%% Visual check of the area sizes for the arena

% vidpath = "S:\Evyn\DATA\Courtship Videos\02.05.2025\Berlin_courtship_F_LRR_caviar_ramp1\compiled_video_1.avi";
% movieInfo = VideoReader(vidpath); %read in video
% demoImg = (read(movieInfo,1));
% 
% fig = getfig('',1,[900,900]);
%     imshow(demoImg); hold on
%     viscircles(well.center(5,:),conversion(4).R*pix2mm,'color', 'w','LineWidth',0.25);
%     viscircles(well.center(5,:),conversion(4).circle75*pix2mm,'color', 'w','LineWidth',0.25);
%     viscircles(well.food,conversion(4).circle10*pix2mm,'color', Color('gold'),'LineWidth',0.25);
%     viscircles(well.food,conversion(4).circle7*pix2mm,'color', Color('red'),'LineWidth',0.25);
%     viscircles(well.food,conversion(4).circle5*pix2mm,'color', Color('purple'),'LineWidth',0.25);
% 
% 
% % vidpath = "S:\Evyn\DATA\Raw Data\05.01.2025\C1_hold_33_caviar_1.avi";


%% RUN THIS TO CHECK THE FIT OF THE OUTER RING 
% TODO (3/11) run this on all the seperate trials between high res and low res to be
% sure that they are aligned and captureing the appropriate space limits

% % Current trial movie path: 
% loc_str = ['Experiment: ' parameters.date ' ' parameters.expID];
% disp(loc_str)
% dummypath = uigetdir;
% movieInfo = VideoReader([dummypath '/compiled_video_1.avi']); %read in video
% demoImg = (read(movieInfo,1));
% 
% % ArenaArea = 2827.43;
% R =  25.6; % (full arena edge: 29.5)   %functional distance they can reach in the circle
% innerR = R*sqrt(3/4); % radius of the inner 50% occupancy space R*sqrt(1/2)
% dist_from_edge = (R - innerR); % distance acceptable for start of outer 25%
% maxR = R*sqrt(0.1); % radius of a circle occupying 10% of the arena
% 
% centre = well.center(5,:);
% foodWell = well.center(well.food_idx,:);
% figure; imshow(demoImg);
%     viscircles(centre,R*pix2mm,'color', foreColor,'LineWidth',0.25);
%     viscircles(centre,innerR/pix2mm,'color', foreColor,'LineWidth',0.25);
%     viscircles(foodWell,maxR/pix2mm,'color', foreColor,'LineWidth',0.25);
%     hold on
% % plot all the position points within this experiment: 
% 
% x = data(M).rawX(:,body.head);
% y = data(M).rawY(:,body.head);
% scatter(x,y,1,'b','filled')
% x = data(F).rawX(:,body.head);
% y = data(F).rawY(:,body.head);
% scatter(x,y,1,'m','filled')








