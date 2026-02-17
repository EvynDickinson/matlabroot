
% Code to find and sort the high resolutions skip frame errors in tracking

%% FIGURES: look at speed distributions

fig = figure; 
for ii = 1:num.trials
    subplot(7,2,ii);
    hold on
    histogram(fly(ii).m.speed,'FaceColor','b', 'EdgeColor','none')
    histogram(fly(ii).f.speed, 'FaceColor', 'm', 'EdgeColor', 'none')
    set(gca, 'xscale', 'log', 'yscale', 'log')
end
formatFig(fig, blkbgd, [7,2]);
matchAxis(fig);
xlabel('speed (mm/s)')


fig = figure; 
for ii = 1:num.trials
    subplot(7,2,ii);
    hold on
    histogram(fly(ii).m.speed,'FaceColor','b', 'EdgeColor','none')
    histogram(fly(ii).f.speed, 'FaceColor', 'm', 'EdgeColor','none')
    set(gca, 'xscale', 'log', 'yscale', 'log')
    % xlim([0,50])
end
formatFig(fig, blkbgd, [7,2]);
% matchAxis(fig);
for ii = 1:num.trials
    subplot(7,2,ii);
    % set(gca, 'yscale', 'log')
    set(gca, 'xscale', 'linear', 'yscale', 'log')
    xlim([-5, 200])
end
xlabel('speed (mm/s)')


% Combined group histogram: 
all_speed = [data.speed(:)];

fig = figure; 
histogram(all_speed,0:10:1000,'FaceColor',Color('vaporwavegren'), 'FaceAlpha', 1)
set(gca, 'xscale', 'linear', 'yscale', 'linear')
formatFig(fig, blkbgd);
xlabel('speed (mm/s)')
ylabel('all flys (frame count)')

thresh = 50;
tot = sum(all_speed>=thresh);
percent_frames_above_threshold = (tot/length(all_speed)*100);

% Where in time over the experiment are the jumps the fastest?

all_speed = [squeeze(data.speed(:,M,:)), squeeze(data.speed(:,F,:))];
tot = sum(all_speed>=thresh,2);
full_total = sum(tot);

frames = [];
for ii = 1:num.trials*2
    frames = [frames; find(all_speed(:,ii)>=thresh)];
end
time = data.time(frames);


r = 5;
c = 1;
sb(1).i = 1;
sb(2).i = 2:5;
xlims = [];

fig = figure; 
subplot(r,c,sb(1).i)
    plot(data.time, data.temp, 'color', Color('vaporwavegren'), 'linewidth', 2)
    ylabel('temp (\circC)')
    xlims = [xlims, xlim];
    title(sprintf('Speed over %d mm/s', thresh),'color', foreColor)
subplot(r,c,sb(2).i)
    histogram(time,'FaceColor',Color('vaporwavegren'), 'FaceAlpha', 1)
    xlabel('time in experiment (min)')
    ylabel('speed over threshold count')
    xlims = [xlims, xlim];

formatFig(fig, blkbgd,[r,c],sb);
subplot(r,c,sb(1).i)
xlim([min(xlims), max(xlims)])
set(gca, 'xcolor', 'none')
subplot(r,c,sb(2).i)
xlim([min(xlims), max(xlims)])

%% FIGURE: targeted look at the nth fastest speed for a selected video trial
clearvars('-except',initial_var{:})

saveDir = createFolder([figDir, 'extreme speed troubleshooting/']);

% Inter-fly-distance from the fly's center point
trial = 9;
maxtrial = 1; % what peak speed do you want to look at? max = 1, second fastest = 2 etc....

speed_thresh = 100;
x1 = fly(trial).m.pos(:,body.center,1); % x location for male center
y1 = fly(trial).m.pos(:,body.center,2);
x2 = fly(trial).f.pos(:,body.center,1); % x location for female center
y2 = fly(trial).f.pos(:,body.center,2);

D = hypot(diff(x1), diff(y1)); 
test = smooth(D,6,'moving');
a = find(D>=speed_thresh); % where does speed exceed 50mm/s
[maxSpeed, b] = maxk(D,maxtrial); % custom selected range that has two high speed frames
if isempty(b)
    disp('no examples of excess speed')
    return
end
B = b(maxtrial)-3:b(maxtrial)+3;
printstring = '\n\t trial had %i frames above %d mm/s\n\t max speed: %i mm/s\n';
titleString = sprintf(printstring, length(a), speed_thresh, round(maxSpeed(maxtrial)));
% fprintf(printstring, length(a), speed_thresh, maxSpeed)
disp(titleString)

% fig = figure; hold on
%     plot(D,'color', Color('vaporwavepurple'))
%     plot(test, 'color', Color('vaporwaveyellow'))
%     h_line(50,'vaporwaveblue', '-',2)
%     v_line([B(1),B(end)], 'metrored','-',2)
%     xlabel('frame number')
%     ylabel('fly speed (mm/s)')
%     formatFig(fig, blkbgd);
%     title(titleString,'Color',foreColor)
% 
%     save_figure(fig,[saveDir  ' ' fly(trial).name  ' speed over time trial ' num2str(maxtrial)]);


% plot out the fly locations for these frames
fig = figure; hold on
for ii = 1:length(B)
    mColor = Color('vaporwaveblue');
    fColor = (Color('vaporwavepink'));
    if B(ii)==b(maxtrial) || B(ii)==b(maxtrial)+1
        mColor = Color('floralblue');
        fColor = Color('floraldarkpink');
    end
    sex = 'm';
    x =  fly(trial).(sex).pos(B(ii),:,1);
    y =  fly(trial).(sex).pos(B(ii),:,2);
    plotFlySkeleton(fig, x,y,mColor,true);
    sex = 'f';
    x =  fly(trial).(sex).pos(B(ii),:,1);
    y =  fly(trial).(sex).pos(B(ii),:,2);
    plotFlySkeleton(fig, x,y,fColor,true);
end
% format
title(titleString,'color', foreColor)
axis equal
formatFig(fig, blkbgd);
set(gca, 'xcolor', 'none', 'ycolor', 'none')
disp('Speeds')
disp(D(B))
disp('Frame Numbers')
disp(B')

save_figure(fig,[saveDir  ' ' fly(trial).name  ' speed-max trial ' num2str(maxtrial)]);

% print location information: 
fullFrame = b(maxtrial);
fprintf('\n\t skip frame: %i\n\t video number: %i\n\t video frame: %i\n',...
            fullFrame,fly(trial).T.vidNums(fullFrame), fly(trial).T.vidFrame(fullFrame))

%% FIGURE: Plot out the 8 fastest speed frames to see if they look like skips and their frame numbers for each trial
clearvars('-except',initial_var{:})

saveDir = createFolder([figDir, 'extreme speed troubleshooting/']);

% Inter-fly-distance from the fly's center point
nFrames = 8; % 2x4 top 8 male and top 8 female 
r = 2;
c = 4;
sz = 10;

for trial = 1:num.trials
    x1 = fly(trial).m.pos(:,body.center,1); % x location for male center
    y1 = fly(trial).m.pos(:,body.center,2);
    x2 = fly(trial).f.pos(:,body.center,1); % x location for female center
    y2 = fly(trial).f.pos(:,body.center,2);
    sexList = {'M','F'};
    
    for sex = 1:1
        titleString = {};
         if sex==M
                % male speed: 
                D = hypot(diff(x1), diff(y1)); 
            else
                D = hypot(diff(x2), diff(y2)); 
         end
        [maxSpeed, b] = maxk(D,nFrames); % custom selected range that has two high speed frames
    
         % PLOT DATA
        fig = getfig('',1);   
        for ff = 1:nFrames
            subplot(r,c,ff); hold on
            roi = b(ff)-3:b(ff)+3; % frames to plot
            for ii = 1:length(roi)
                mColor = Color('vaporwaveblue');
                fColor = Color('vaporwavepink');
                if roi(ii)==b(ff) || roi(ii)==b(ff)+1
                    mColor = Color('floralblue');
                    fColor = Color('floraldarkpink');
                end
                x =  fly(trial).m.pos(roi(ii),:,1);
                y =  fly(trial).m.pos(roi(ii),:,2);
                plotFlySkeleton(fig, x,y,mColor,true,sz);
                x =  fly(trial).f.pos(roi(ii),:,1);
                y =  fly(trial).f.pos(roi(ii),:,2);
                plotFlySkeleton(fig, x,y,fColor,true,sz);
            end
            titleString{ff} = {sprintf('%s | %i mm/s', sexList{sex}, round(maxSpeed(ff)));...
                                      sprintf('frame %d', b(ff))};
        end
        % format figure: 
        formatFig(fig, blkbgd,[r,c]);
         for ff = 1:nFrames
            subplot(r,c,ff)
            axis tight square 
            title(titleString{ff}, 'color', foreColor,'fontsize', 10)
            set(gca, 'xcolor', 'none', 'ycolor', 'none')
         end
    
         % Save the figure: 
        save_figure(fig,[saveDir  fly(trial).name  ' max speed trials ' sexList{sex}],fig_type, true,true, '-r100');
    end
end

%% TODO: pull up image frames of the lower speed, non-swap jumps to see if they are indeed a ballistic jump

%% FIGURES: double speed locations
clearvars('-except',initial_var{:})
saveDir = createFolder([figDir, 'extreme speed troubleshooting/']);

trial = 2;
sz = 10;
edgeFrames = 5;
speed_thresh = 50; %mm/s speed threshold

mBaseColor = Color('vaporwaveblue');
fBaseColor = Color('vaporwavepink');
mColor_jump = Color('floralblue');
fColor_jump = Color('floraldarkpink');

% locations based on both flies having speed jumps
double_speed = squeeze(data.speed(:,:,trial)); % locations for speeding in both sexes
double_speed_loc = double_speed>=speed_thresh;
confident_frames = find(sum(double_speed_loc,2)==2);
% diff(confident_frames)
fprintf('\nThere are %d likely swap frames\n', length(confident_frames))

maxExmp = 8; % how many different examples to show
r = 2; c = 4; % rows and columns for the figure
% how many figures with subplots will it take show all the different examples
nFigs = ceil(length(confident_frames)/(r*c)); 

% code efficient preallocation
mPos = fly(trial).m.pos;
fPos = fly(trial).f.pos;
allROI = confident_frames + (-edgeFrames:edgeFrames);
% create matrix of the frames for each figure:
frameList = nan([1,nFigs * r * c]);
frameList(1:length(confident_frames)) = confident_frames;
frameList = (reshape(frameList, [r*c, nFigs]))';


ax = [];
idx = 0;
for gg = 1:nFigs

    % PLOT DATA
    fig = getfig('',1);   
    % set(fig, 'Visible', 'off') % hold off on drawing the figure until it is all plotted
    for ff = 1:maxExmp
        subplot(r,c,ff); hold on
        idx = idx + 1;
        start_frame = frameList(gg,ff);
        if isnan(start_frame)
            continue
        end

        % Initialize fly skeletons for plotting later
        ax = gca;
        % hM = initFlySkeleton(ax, mBaseColor, true, sz);
        % hF = initFlySkeleton(ax, fBaseColor, true, sz);

        roi = allROI(idx,:);% frames to plot
        for ii = 1 : length(roi)
            % select the appropriate plotting color
            if roi(ii)==start_frame %|| roi(ii)==stop_frame
                mColor = mColor_jump;
                fColor = fColor_jump;
            else
                mColor = mBaseColor;
                fColor = fBaseColor;
            end
            x =  mPos(roi(ii),:,1);
            y =  mPos(roi(ii),:,2);
            plotFlySkeleton(fig, x, y, mColor, true, sz);
            % updateFlySkeleton(hM, x, y, mColor);
            x =  fPos(roi(ii),:,1);
            y =  fPos(roi(ii),:,2);
            % updateFlySkeleton(hF, x, y, fColor);
            plotFlySkeleton(fig, x, y, fColor, true, sz);
        end
    end

    % format figure: 
    formatFig(fig, blkbgd,[r,c]);
     for ff = 1:maxExmp
        subplot(r,c,ff);
        axis equal tight
        set(gca, 'xcolor', 'none', 'ycolor', 'none')
        start_frame = frameList(gg,ff);
        if isnan(start_frame)
            continue
        end
        title(num2str(frameList(gg,ff)), 'color', foreColor,'fontsize', 10)
        % find the size and build a rectagle 'frame'
        ax = axis; % [xmin xmax ymin ymax]
        rectangle('Position',[ax(1) ax(3) ax(2)-ax(1) ax(4)-ax(3)],...
              'EdgeColor',foreColor,'LineWidth',0.5, 'Curvature',[0.1 0.1],...
              'linestyle', ':')
     end
     set(fig, 'Visible', 'on');
     % drawnow
    % save_figure(fig, [saveDir fly(trial).name  ' speed error trials ' num2str(gg)],fig_type, true,false, '-r100');
end


%% ANALYSIS & FIGURES: identify likely frame swap locations
if ~any(strcmp('keyFrames', initial_var))
    initial_var{end+1} = 'keyFrames';
end
clearvars('-except',initial_var{:})
saveDir = createFolder([figDir, 'extreme speed troubleshooting/']);

speed_thresh = 35; %mm/s speed threshold
skip_threshold = 3; % how many frames for a confident swap pair in time

% initialize empty variables
% Total_High_Speed = nan([num.trials, 2]);
% [Confident_Swaps, Likely_Swaps, All_Speed, Confident_Speed, Likely_Speed] = ...
%     deal(nan([num.trials, 1]));

keyFrames = []; % data structure holding the frame data 

for trial = 1:num.trials
    % MALE speed data for current trial 
    mspeed = squeeze(data.speed(:,M,trial));
    mspeed_loc = mspeed>=speed_thresh; % frames with above threshold speed
    mspeed_frames = find(mspeed_loc); % pull frames above the limit speed
    keyFrames(trial).mspeed_frames = mspeed_frames;
    
    % FEMALE speed data for current trial 
    fspeed = squeeze(data.speed(:,F,trial));
    fspeed_loc = fspeed>=speed_thresh; % frames with above threshold speed
    fspeed_frames = find(fspeed_loc); % pull frames above the limit speed
    keyFrames(trial).fspeed_frames = fspeed_frames;
    
    % Possible swap locations: (based on high speed for both flies)
    double_speed = squeeze(data.speed(:,:,trial)); % speed in both sexes
    double_speed_loc = double_speed>=speed_thresh;
    swap_ConfidentFrames = find(sum(double_speed_loc,2)==2); % where are both M&F speeding
    
    % Possible swap locations if one of the sexes does not have a labeled skeleton
    a = fspeed_loc & isnan(mspeed);
    b = mspeed_loc & isnan(fspeed);
    possible_swapFrames = find(a | b);
    if ~isempty(possible_swapFrames)
        allSwaps = sort([swap_ConfidentFrames; possible_swapFrames]); % all frames with double speed MF or single+nan
        keyFrames(trial).swap_LikelyFrames = allSwaps; % all frames with double speed M F
    else % no changes due to new pairings, so the same frames as confident
        keyFrames(trial).swap_LikelyFrames = swap_ConfidentFrames;
    end
end

% PAIR LIKELY SWAP FRAMES FOR EACH OF THE TRIALS
for trial = 1:num.trials
    a = keyFrames(trial).swap_LikelyFrames;% current list of possible frames for pairs
    % Find differences between consecutive frames
    frame_diffs = diff(a);
    % Find where frames are close together (potential pairs)
    is_close = frame_diffs <= skip_threshold;
    % Extract pairs
    pairs = [];
    ii = 1; % index location within frame list
    while ii <= length(is_close)
        if is_close(ii)
            % Found a pair: frames a(i) and a(i+1)
            pairs = [pairs; a(ii), a(ii+1)];
            
            % Skip the next position to avoid overlapping pairs
            % (e.g., if frames 1,2,3 are all close, we want pairs [1,2] and not [2,3])
            ii = ii + 2;
        else
            ii = ii + 1;
        end
    end
    % Store the pairs
    keyFrames(trial).frame_pairs = pairs;
end

% Identify nearest neighbors to see if the swap is likely for the
% identified swap locations
% look at the avg distance to the fly center for each of the
% flies in the preceeding frames 
nNeighbors = 5; % how many preceeding frames to look at for close distance?
for trial = 1:num.trials 
    keyFrames(trial).frame_pairs_swap = [];
    mPos = squeeze(fly(trial).m.pos(:,body.center,:));
    fPos = squeeze(fly(trial).f.pos(:,body.center,:));

    % for each of the likely speed frame pairs
    likely_pairs = keyFrames(trial).frame_pairs;
    rois = likely_pairs(:,1) + ((-nNeighbors-1): -1);
   
    for ii = 1:size(likely_pairs,1)
        
        % frame locations for pre and during swaps
        pre_roi = rois(ii,:);
        dur_roi = likely_pairs(ii,1);
        
        % distance of each during swap roi to pre-swap roi relative
        x = mPos(dur_roi,1) - mPos(pre_roi,1); 
        y = mPos(dur_roi,2) - mPos(pre_roi,2);
        dM = mean(hypot(x,y),'omitnan');  % male track distance to assigned male (distance to male)
        x = mPos(dur_roi,1) - fPos(pre_roi,1); 
        y = mPos(dur_roi,2) - fPos(pre_roi,2);
        dF = mean(hypot(x,y),'omitnan');  % male track distance to assigned female (distance to female)

        % determine if likely switch
        likelyswitch = dM>=dF; % female would be closer than correct male
        keyFrames(trial).frame_pairs_swap(ii) = likelyswitch;
    end
end


% ======== FIGURE =========
% Quick look at the jump vs swap stats
[m_frames, f_frames,double_frames, maybe_pair_frames, final_pair_frames] = deal(nan(num.trials, 1));
for trial = 1:num.trials
    m_frames(trial) = size(keyFrames(trial).mspeed_frames,1);
    f_frames(trial) = size(keyFrames(trial).fspeed_frames,1);
    double_frames(trial) = size(keyFrames(trial).swap_LikelyFrames,1);
    maybe_pair_frames(trial) = numel(keyFrames(trial).frame_pairs);
    final_pair_frames(trial) = sum(keyFrames(trial).frame_pairs_swap)*2;
end
% relative percentages
tot_frames = size(fly(1).T,1);
m_framesP = mean((m_frames/tot_frames)*100);
f_framesP = mean((f_frames/tot_frames)*100);
swap_percentage =  mean((final_pair_frames ./ tot_frames)*100,'omitnan');
fprintf('\n%2.3g  percent of total frames are high speed \n', (m_framesP + f_framesP))
fprintf('\n%2.5g  percent of total frames are confident swaps',swap_percentage);


% plotting parameters
mt = ones([num.trials,1]);
r = 1; c = 2;
sz = 50; FA = 0.7;
mColor = Color('vaporwaveblue');
fColor = Color('vaporwavepink');

fig = getfig('',1);
    % total high speed frames for male and female flies
    subplot(r,c,1); hold on
        plot([mt, 2*mt]', [m_frames, f_frames]', 'color', foreColor, 'LineWidth', 1.5)
        scatter(mt, m_frames, sz, mColor, "filled", 'MarkerFaceAlpha',FA)
        scatter(2*mt, f_frames, sz, fColor, "filled", 'MarkerFaceAlpha',FA)
        xlim([0.7, 2.2])
        set(gca, 'xtick', 1:2, 'xticklabel', {'male', 'female'})
        ylabel('high speed frame counts')
        title('all high speed events')
    
    % likely swap pairs
    subplot(r,c,2); hold on
        plot([mt, 2*mt, 3*mt]', [double_frames, maybe_pair_frames, final_pair_frames]', 'color', foreColor, 'LineWidth', 1.5)
        scatter([mt; 2*mt; 3*mt], [double_frames; maybe_pair_frames; final_pair_frames]', sz, Color('vaporwavepurple'), "filled", 'MarkerFaceAlpha',FA)
        xlim([0.7, 3.2])
        set(gca, 'xtick', 1:3, 'xticklabel', {'all', 'likely', 'confident'})
        title('Swap Frame Counts')
        ylabel('high speed frame counts')
    formatFig(fig, blkbgd, [r,c]);
save_figure(fig, [saveDir  'Jump vs swap frame number comparisons']);





% FIGURES: plot out all the fly skeletons for the confident trials based
% on the double fly speed criteria
maxExmp = 50; % how many different examples to show
r = 5; c = 10; % rows and columns for the figure
% how many figures with subplots will it take show all the different examples
for trial = 1:num.trials
    frames = keyFrames(trial).swap_LikelyFrames;
    frame_buff = 3; % how many frames on either side of the target 
    nFigs = ceil(length(frames)/(r*c)); 
    
    % code efficient preallocation
    mPos = fly(trial).m.pos;
    fPos = fly(trial).f.pos;
    roi = frames + (-frame_buff:+frame_buff);% frames to plot
    % create matrix of the frames for each figure:
    frameList = nan([1, nFigs * r * c]);
    frameList(1:length(frames)) = frames;
    frameList = (reshape(frameList, [r*c, nFigs]))';
    % color choices
    mBase = Color('vaporwaveblue');
    fBase = Color('vaporwavepink');
    mHighlight = Color('floralblue');
    fHighlight = Color('floraldarkpink');
    
    frame_pairs = keyFrames(trial).frame_pairs;
    
    % -------------- PLOT DATA ----------------
    for gg = 1:nFigs
        fig_str = sprintf('Trial %d | fig %d', trial, gg);
        fig = getfig(fig_str, 0, [3306 1319]);   
        % subplots for each specific instance of high speed
        for ff = 1:maxExmp 
            subplot(r,c,ff); hold on
            loc = frameList(gg,ff); % identified target frame
            if isnan(loc) % skip this plot if there isnt data
                continue
            end
            % plot data for the frames around the key frame
            idx = ((gg-1)*(r*c)) + ff; % what frame in the long list is this
            % fprintf('\nCurrent Frame Count %i', idx)
            for ii = 1:size(roi,2) % for all frames around the target
                curr_frame = roi(idx,ii); % frame to plot
                % select the color based on the frame number
                if loc==frames(idx)
                    mColor = mHighlight;
                    fColor = fHighlight;
                else 
                    mColor = mBase;
                    fColor = fBase;
                end
                % extract and plot the location data for the flies
                x =  mPos(curr_frame,:,1);
                y =  mPos(curr_frame,:,2);
                plotFlySkeleton(fig, x, y, mColor, false);
                x =  fPos(curr_frame,:,1);
                y =  fPos(curr_frame,:,2);
                plotFlySkeleton(fig, x, y, fColor, false);
            end
        end
    
        % format figure: 
        formatFig(fig, blkbgd,[r,c]);
         for ff = 1:maxExmp
            subplot(r,c,ff)
            frame_num = frameList(gg,ff);
            set(gca, 'xcolor', 'none', 'ycolor', 'none')
            if isnan(frame_num)
                continue
            end
            axis equal tight
            title(frame_num, 'color', foreColor,'fontsize', 10)
            % give frame pair numbers if it exists
            frame_pair_loc =  find(frame_num==frame_pairs(:,1) | frame_num==frame_pairs(:,2)) ;
            if ~isempty(frame_pair_loc)
                % now check if it made it past second layer distance filter: 
                if keyFrames(trial).frame_pairs_swap(frame_pair_loc)
                    kolor = Color('vaporwavegren');
                else
                    kolor = Color('vaporwaveyellow');
                end
                    title_str = [num2str(frame_num) ' | ' Alphabet(frame_pair_loc)];
                    title(title_str, 'color', foreColor,'fontsize', 10)
            else
                kolor = foreColor;
            end
    
            % find the size and build a rectagle 'frame'
            ax = axis; % [xmin xmax ymin ymax]
            rectangle('Position',[ax(1) ax(3) ax(2)-ax(1) ax(4)-ax(3)],...
                  'EdgeColor',kolor,'LineWidth',0.5, 'Curvature',[0.1 0.1],...
                  'linestyle', ':')
         end
        % save_figure(fig, [saveDir fly(trial).name  ' speed error trials ' num2str(gg)],fig_type, true,false, '-r100');
    end

end


%% PULL UP IMAGES OF HIGH SPEED FRAMES FROM LOCAL DRIVE DATA
clearvars('-except',initial_var{:})

% load saved frame identity information or generate it: 
% get computer name: 
if ismac 
    computerName = getenv('USER');
else 
    computerName = getenv('COMPUTERNAME');
end
file_name = [groupDir, computerName ' speed labeling.mat'];
% load in preexisting data or make a new blank template
if ~(exist(file_name, 'file')==2)
    speedTest = struct;
    nFrames = 50; % how many frames to test for this fly for female and male
    for trial = 1:num.trials
        % select random frames to label
        possible_frames = keyFrames(trial).fspeed_frames;
        idx_M = randi(length(possible_frames),[nFrames, 1]);
        frames_M = possible_frames(idx_M);
        possible_frames = keyFrames(trial).mspeed_frames;
        idx_F = randi(length(possible_frames),[nFrames, 1]);
        frames_F = possible_frames(idx_F);
        speedTest(trial).idx = sort([frames_M; frames_F]);
        % create empty label structure
        speedTest(trial).label = cell([nFrames*2, 1]);
        speedTest(trial).labeled = false([nFrames*2, 1]);
    end
    save(file_name,'speedTest');
else 
    switch questdlg('Load speed labeling data?')
        case 'Yes'
            load(file_name);
            disp('loaded data file')
        case {'No', 'Cancel', ''}
            return
    end
end


% SELECT A TRIAL TO LABEL:
% which trials are on a drive that is present at this computer?
expList = {fly(:).name};
ExpOnDrive = cell(size(expList)); % cell with location addresses for the videos
ExpsFound = false(size(expList)); % logical for local copy of the videos
% find storage drives
drvStr = 'Data Storage';
possible_drives =  findDriveByName(drvStr); % find drive letter associated with data storage, if possible
% find list of experiments on the local drives
if ~isempty(possible_drives)
    for ii  = 1:length(possible_drives)
        rootPath = [possible_drives{ii}, 'Courtship Videos/'];
        allFolders = dir(rootPath);
        expDates = {allFolders(:).name}; % experiments on the drive
        % compare between current experiments and the list of names
        for jj = 1:length(expList)
            datestring = expList{jj}(1:10); % current experiment to look for
            if any(contains(expDates, datestring))
                % extract the full file path to the videos
                ExpOnDrive{jj} = [rootPath datestring '/' expList{jj}(12:end) '/'];
                ExpsFound(jj) = true;
            end
        end
    end
end

% give option to pull up instances from one of the possible trials
videoList = {fly(ExpsFound).name};
idx = listdlg("PromptString",'Select which experiment to demo images from',...
            'SelectionMode','single',...
            'ListString',videoList,...
            'ListSize',[350,400]);
trial = find(strcmp(videoList{idx},{fly(:).name})); % trial number selected
% location of the video files for this trial on the current computer: 
video_root = [ExpOnDrive{trial}  'AVI_Files/compiled_video_'];

% LABEL FRAMES FOR THE SELECTED FILE

% figure parameters
pixbuff = 60; % pixel border buffer size
nPre = 2; % before event frames to show
nPost = 3;  % after event frames to show
shiftedframeList = speedTest(trial).idx;
% adjust frameList for possible shift in timing across videos
frameList = data.frame(shiftedframeList,trial);
c = nPre + nPost +1;
r = 2; % female and male

% extact data for this trial
mXPos = squeeze(fly(trial).m.pos(:,:,1));
mYPos = squeeze(fly(trial).m.pos(:,:,2));
fXPos = squeeze(fly(trial).f.pos(:,:,1));
fYPos = squeeze(fly(trial).f.pos(:,:,2));
mColor = Color('vaporwaveblue');
fColor = Color('Vaporwavepink');

for ii = 1:5% length(frameList)
    % only run for unlabeled frames
    if ~speedTest(trial).labeled(ii)
        frame = frameList(ii);
        trial_frames = frame + (-nPre:nPost);
        vidNum = fly(trial).T.vidNums(frame);
        vidFrames = fly(trial).T.vidFrame(frame) + (-nPre:nPost);

        % check if frames could actually exist - if not, skip this one
        if any(vidFrames<=0 | vidFrames>5400)
            continue
        end

        % load video frames
        vidPath = [video_root num2str(vidNum) '.avi'];
        movieInfo = VideoReader(vidPath); %read in video
        img = read(movieInfo, [vidFrames(1), vidFrames(end)]);

        % find the position data of the male fly
        mX = mXPos(trial_frames,:);
        mY = mYPos(trial_frames,:);
        mZoomX = [min(min(mX)) - pixbuff, max(max(mX)) + pixbuff];
        mZoomY = [min(min(mY)) - pixbuff, max(max(mY)) + pixbuff];

        % find the position data of the female fly
        fX = fXPos(trial_frames,:);
        fY = fYPos(trial_frames,:);
        fZoomX = [min(min(fX)) - pixbuff, max(max(fX)) + pixbuff];
        fZoomY = [min(min(fY)) - pixbuff, max(max(fY)) + pixbuff];
        

        fig = getfig('',1); set(fig, 'color', 'k');
        for ff = 1 : c
            % male fly
            subplot(r, c, ff)
            imshow(img(:,:,:,ff))
            xlim(mZoomX)
            ylim(mZoomY)
            % plot skeleton on top?
            hold on
            plotFlySkeleton(fig, mX(ff,:)', mY(ff,:)', mColor, false);

            % female fly
            subplot(r, c, ff+c)
            imshow(img(:,:,:,ff))
            xlim(fZoomX)
            ylim(fZoomY)
            hold on
            plotFlySkeleton(fig, fX(ff,:)', fY(ff,:)', fColor, false);
        end

        % label behavior: 
        switch questdlg('Behavior?', '', 'Jump', 'Swap', 'Other', 'Jump')
            case 'Jump'
                speedTest(trial).label = 'Jump';
            case 'Swap'
                speedTest(trial).label = 'Swap';
            case 'Other'
                speedTest(trial).label = 'Other';
            case '' % canceled 
                disp('Labeling Canceled')
                return
        end
        
        % behavior was labeled
        speedTest(trial).labeled(ii) = true;
        close(fig)
    end

end


% SAVE DATA STRUCTURE WHEN DONE
save(file_name,'speedTest');
 




%% HIGH SPEED TYPE SORTING

% which of the high speed frames are jumps vs which are swaps? 
% pull up pic evidence for the highest n number of trials
% pull up pic evidence for a random selection of the lower trials

% quantify the number of jumps vs swaps in the random sample

% Base Rules: 
% jumps: single fly above speed threshold and IS data for other fly body 
% swaps: both flies above speed threshold and PAIRED frames OR missing data
% for one of the other flies but still close by frame with above speed
% threshold for the non-nan fly


%% MANY FIGURES: Check close frame high speeds for possible switch patterns: 
clearvars('-except',initial_var{:})
saveDir = createFolder([figDir, 'extreme speed troubleshooting/']);

trial = 2;
sz = 10;
sex = 1;
speed_thresh = 30; %mm/s speed threshold

% speed data for current trial 
speed = squeeze(data.speed(:,sex,trial));
speed_loc = speed>=speed_thresh;

% pull frames above the limit speed
speed_frames = find(speed_loc); 

% locations based on both flies having speed jumps
double_speed = squeeze(data.speed(:,:,trial)); % locations for speeding in both sexes
double_speed_loc = double_speed>=speed_thresh;
confident_frames = find(sum(double_speed_loc,2)==2);
fprintf('\nThere are %d likely swap frames\n', length(confident_frames))

% what portion of total over-speed instances is this
proportion_of_frames = (length(confident_frames)/sum(speed_loc))*100;
% how similar is the speed between the male and female flies when they swap
% locations?
double_speed(confident_frames,:)

% how many of these are likely swaps?
a = find(diff(speed_frames)<=3); % possible start of switch location in speed frames
b = a+1;
possible_switches_locs = speed_frames(a); % where are there speed jumps nearly back to back
reverse_switch_locs = speed_frames(b); 

maxExmp = 8; % how many different examples to show
r = 2; c = 4; % rows and columns for the figure
% how many figures with subplots will it take show all the different examples
nFigs = ceil(length(possible_switches_locs)/(r*c)); 
titleString = {};
% check how these look 
curr_frame = 0; % counter for switch location
for gg = 1:nFigs

    % PLOT DATA
    fig = getfig('',1);   
    for ff = 1:maxExmp
        subplot(r,c,ff)
        curr_frame = curr_frame + 1;
        start_frame = possible_switches_locs(curr_frame);
        stop_frame = reverse_switch_locs(curr_frame);
        roi = start_frame-3:stop_frame+3;% frames to plot

        for ii = 1:length(roi)
            mColor = Color('vaporwaveblue');
            fColor = Color('vaporwavepink');
            if roi(ii)==start_frame || roi(ii)==stop_frame
                mColor = Color('floralblue');
                fColor = Color('floraldarkpink');
            end
            x =  fly(trial).m.pos(roi(ii),:,1);
            y =  fly(trial).m.pos(roi(ii),:,2);
            plotFlySkeleton(fig, x,y,mColor,true,sz);
            x =  fly(trial).f.pos(roi(ii),:,1);
            y =  fly(trial).f.pos(roi(ii),:,2);
            plotFlySkeleton(fig, x,y,fColor,true,sz);
        end
        if any(start_frame==confident_frames) || any(stop_frame==confident_frames)
            titleString{ff} = 'Confident';
        else 
            titleString{ff} = '';
        end
    end

    % format figure: 
    formatFig(fig, blkbgd,[r,c]);
     for ff = 1:maxExmp
        subplot(r,c,ff)
        axis equal tight
        title(titleString{ff}, 'color', foreColor,'fontsize', 10)
        set(gca, 'xcolor', 'none', 'ycolor', 'none')
        % find the size and build a rectagle 'frame'
        ax = axis; % [xmin xmax ymin ymax]
        rectangle('Position',[ax(1) ax(3) ax(2)-ax(1) ax(4)-ax(3)],...
              'EdgeColor',foreColor,'LineWidth',0.5, 'Curvature',[0.1 0.1],...
              'linestyle', ':')
     end
    save_figure(fig, [saveDir fly(trial).name  ' speed error trials ' num2str(gg)],fig_type, true,false, '-r100');
end

%% ANALYSIS :  can we pair skip frames together?
% how long is the avg duration of a missed swapped identity?
% how many frames above 50mm/s?

clearvars('-except',initial_var{:})

trial = 2;

nNeighbors = 5; % how many past frames to compare for location
sz = 10;
buff = 0.3;
% pix2mm = conversion(4).pix2mm;

% find all frames above a speed threshold
sex = 1;
sexList = {'m', 'f'};
speed_thresh = 50; %mm/s speed threshold

% distance data for current trial 
x1 = fly(trial).m.pos(:,body.center,1); % x location for male center
y1 = fly(trial).m.pos(:,body.center,2);
x2 = fly(trial).f.pos(:,body.center,1); % x location for female center
y2 = fly(trial).f.pos(:,body.center,2);
dM = hypot(diff(x1), diff(y1));  % male distance
dF = hypot(diff(x2), diff(y2));  % female distance

speed = squeeze(data.speed(:,sex,trial));

speed_loc = speed>=speed_thresh;
speed_frames = find(speed_loc); % what frames have above limit speeds?

% how many of these are likely swaps?
a = find(diff(speed_frames)<=3); % possible start of switch location in speed frames
b = a+1;
possible_switches_locs = speed_frames(a); % where are there speed jumps nearly back to back
reverse_switch_locs = speed_frames(b); 

% check how these look ...


[possible_switches_locs, reverse_switch_locs]

% is there an average higher speed for the likely swaps based the step size between
% speed frames? 


% look at the avg distance to the fly center for each of the flies in the
% preceeding frames 
for ii = 1:length(speed_frames)
    % current speeding frame: 
    frame = speed_frames(ii);

    switch sex 
        case M
            cX = x1(frame); % male x body center
            cY = y1(frame); % male y body center
            c2X = x1(frame+1); % subsrequent frame male x body center
            c2Y = y1(frame+1); % subsrequent frame male y body center
        case F
            cX = x2(frame); % female x body center
            cY = y2(frame); % female y body center
            c2X = x2(frame+1); % subsrequent frame female x body center
            c2Y = y2(frame+1); % subsrequent frame female y body center
    end

    % find the distance within sex points vs. the across sex points
    past_roi = frame-nNeighbors:frame-1; % frames behind the fast speed frame
    dM = hypot(cX-x1(past_roi), cY-y1(past_roi));  % male distance to target fly point
    dMavg = mean(dM, 'omitnan');
    dF = hypot(cX-x2(past_roi), cY-y2(past_roi));  % female distance to target fly point
    dFavg = mean(dF, 'omitnan');
    switch sex
        case M
            likelyswitch = dMavg>=dFavg; % female would be closer than correct male
        case F 
            likelyswitch = dMavg<=dFavg; % male would be closer than correct female
    end

    % if frame is likely a switch, check if next frame is a switch back...
    if likelyswitch
        % find the distance within sex points vs the across sex points for
        % the subsequent frame to check if it switched back
        d2M = hypot(c2X-x1(past_roi), c2Y-y1(past_roi));  % male distance to target fly point (pixels)
        d2Mavg = mean(d2M, 'omitnan');
        d2F = hypot(c2X-x2(past_roi), c2Y-y2(past_roi));  % female distance to target fly point
        d2Favg = mean(d2F, 'omitnan');
        switch sex
            case M
                likelyreverseswitch = d2Mavg<=d2Favg; % female would be closer than correct male
            case F           
                likelyreverseswitch = d2Mavg>=d2Favg; % male would be closer than correct female
        end
    end
    
    % if the next point switches back, check to see if it is flagged in the speed
    % data list:
    if likelyreverseswitch
       loc = (frame+1==speed_frames);
       speed_frames(loc)
       
       diff(speed_frames)
    
    
    
    
    
    % quick visual test of the currrent frame figure: 
    r = 1; 
    c = 3;
    sb(1).i = 1:2;
    sb(2).i = 3;
    fig = getfig('',1,[688 353]);   
    % fly positions over time
    subplot(r,c,sb(1).i); hold on
        roi = frame-nNeighbors:frame+3; % frames to plot
        for ff = 1:length(roi) % frames before and after the high speed frame
            mColor = Color('vaporwaveblue');
            fColor = Color('vaporwavepink');
            if roi(ff)==frame || roi(ff)==frame+1
                mColor = Color('floralblue');
                fColor = Color('floraldarkpink');
            end
            x =  fly(trial).m.pos(roi(ff),:,1);
            y =  fly(trial).m.pos(roi(ff),:,2);
            plotFlySkeleton(fig, x,y,mColor,true,sz);
            x =  fly(trial).f.pos(roi(ff),:,1);
            y =  fly(trial).f.pos(roi(ff),:,2);
            plotFlySkeleton(fig, x,y,fColor,true,sz);
            % plot the fly with speed frame in forecolor:
            x =  fly(trial).(sexList{sex}).pos(frame,:,1);
            y =  fly(trial).(sexList{sex}).pos(frame,:,2);
            plotFlySkeleton(fig, x,y,Color('gold'),false);
            scatter(x(body.center), y(body.center), 45, Color('gold'), 'filled')
            if likelyswitch
                x =  fly(trial).(sexList{sex}).pos([frame-1,frame+1],body.center,1);
                y =  fly(trial).(sexList{sex}).pos([frame-1,frame+1],body.center,2);
                scatter(x,y,20,foreColor,"filled")
            end
        end
        axis equal
    % quick test of the within and across fly distance
    subplot(r,c,sb(2).i);
    hold on
        scatter(F*ones(size(dM)), dF, 55, fColor, 'filled', 'MarkerFaceAlpha', 0.7, 'xjitter', 'density', 'xjitterwidth', 0.3)
        scatter(M*ones(size(dM)), dM, 55, mColor, 'filled', 'MarkerFaceAlpha', 0.7, 'xjitter', 'density', 'xjitterwidth', 0.3)
        plot([F-buff, F+buff], [dFavg, dFavg], 'color',  foreColor, 'linewidth', 2)
        plot([M-buff, M+buff], [dMavg, dMavg], 'color', foreColor, 'linewidth', 2)
        xlim([0,3])
        ylabel('distance to fast fly (pixels)')
        set(gca, 'xtick', sex, 'xticklabel', 'Within')
        % set(gca, 'xcolor', foreColor)
    % formatting
    formatFig(fig, blkbgd,[r,c], sb);
    subplot(r,c,sb(1).i);
    set(gca,'xcolor', 'none', 'ycolor', 'none')
    if likelyswitch
        title_str = sprintf('%i mm/s | switch likely', round(speed(frame)));
    else
        title_str = sprintf('%i mm/s |switch unlikely', round(speed(frame)));
    end
    title(title_str,'color', foreColor,'FontSize',12)
    


    
    
    


end








































































