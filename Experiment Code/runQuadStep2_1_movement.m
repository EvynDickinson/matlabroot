
function results = runQuadStep2_1_movement(inputPath, autoSave, essentialfigs)

% inputPath = 'G:\My Drive\Jeanne Lab\DATA\06.28.2022\analysis\half processed data.mat';
% essentialfigs = true; % only run the essential figures automatically

%% Load data
warning off
load(inputPath);

disp(['Starting figures for ' folder ' ' expName])
baseFolder = getCloudPath;
analysisDir = [baseFolder folder '/analysis/'];
nArenas = 4;
initial_vars{end+1} = 'nArenas';
initial_vars{end+1} = 'figDir';
figDir = [baseFolder folder '/'];
for arena = 1:nArenas
    arenaData(arena).figDir = [figDir arenaData(arena).name '/'];
end

% 
% tPoints = getTempTurnPoints(expData.parameters.protocol);
% threshHigh = tPoints.threshHigh;
% threshLow = tPoints.threshLow;
% binSpace = 1; %temp degree bin widths

initial_vars = [initial_vars(:)', 'arena', 'threshHigh', 'threshLow', 'binSpace'];

%% ANALYSIS: Segment data into functional tracks 
clearvars('-except',initial_vars{:})
close all


trackROI = struct;

for vid = 1:nvids
    % pull and select data
    if isempty(data(vid).tracks) % add nans for the 'processed data'
        trackROI(vid,arena).start = [];
        trackROI(vid,arena).stop = [];
        vars = {'speed','speedAvg','speedErr'};

        vidLength = length(data(vid).tempLog);
        MT = nan([vidLength,1]);
        for arena = 1:4
            trackROI(vid,arena).tracks = [];
            trackROI(vid,arena).ntrackedlines = zeros([vidLength,1]);
            for i = 1:length(vars)
                trackROI(vid,arena).(vars{i}) = MT;
            end
            data(vid).x_loc_raw = MT;
            data(vid).y_loc_raw = MT;
        end
        continue
    end


    ntracks = size(data(vid).tracks,4);
    flynums = 1:ntracks;
    
    % CList = {'yellow', 'pink', 'purple', 'teal', 'navy', 'green'};
    
    % SEGMENT TRACKS INTO BETWEEN NAN SECTIONS
    % pull image
    movieInfo = VideoReader([figDir,'/',expName,'_',num2str(vid),'.avi']); 
%     img = read(movieInfo,1);
    nframes = movieInfo.NumFrames;
    
    % set an empty structure for tracks to load into:
    for arena = 1:4
        trackROI(vid,arena).tracks = [];
    end
    
    for fly = flynums    
        plotData = squeeze(data(vid).tracks(:,1,:,fly));    
        % APPLY 4 CONDITIONS FOR POINT BEING CONTINOUS WITH PREVIOUS FRAME
        % 0) tracked point
        nanLoc = isnan(plotData(:,1));
        % 1) arena identity?
        X = data(vid).tracks(:,1,1,fly);
        Y = data(vid).tracks(:,1,2,fly);
        % Find points within arena:
        arenaLoc = false(nframes,4);
        for arena = 1:4
            arenaLoc(:,arena) = sqrt((X-arenaData(arena).centre(1)).^2 + (Y-arenaData(arena).centre(2)).^2)<=r; 
        end
        % 2) previous frame is a tracked point
        lastPoint = [true; ~nanLoc(1:end-1)];
        
        % 2a) previous two frames is a tracked point (min 3 points to make
        % a track)
        secondPoint = [true; true; ~nanLoc(1:end-2)];

        % 3) previous frame in the same arena...
        lastArena = false(nframes,4);
        for arena = 1:4
            lastArena(:,arena) = [true; arenaLoc(1:end-1,arena)];
        end
        
        % Are all the conditions met?
        validLoc = false(nframes,4);
        for arena = 1:4 
            validLoc(:,arena) = ~nanLoc & arenaLoc(:,arena) & lastPoint & lastArena(:,arena) & secondPoint;
            % if first frame is good but second is bad, then skip the first
            if validLoc(1,arena)==true && validLoc(2,arena)==false
                validLoc(1,arena) = false;
            end

            % SEGMENT THE TRACKS INTO UNIQUE STRETCHES
            trackROI(vid,arena).start = find(diff(validLoc(:,arena))==1);
            trackROI(vid,arena).stop = find(diff(validLoc(:,arena))==-1);
            % add points for first and last frame if they are valid tracking
            if validLoc(1,arena)
                trackROI(vid,arena).start = [1;trackROI(vid,arena).start];
            end
            if validLoc(end,arena)
               trackROI(vid,arena).stop = [trackROI(vid,arena).stop;nframes];
            end 
            % make sure the number of starts and stops aligns 
            if ~(length(trackROI(vid,arena).start)==length(trackROI(vid,arena).stop))
                warndlg('mismatched start and stop of a tracked segment')
                return
            end
            % make sure that the start comes before each stop...
            if ~all(trackROI(vid,arena).start<trackROI(vid,arena).stop)
                warndlg('mismatched start and stop of a tracked segment')
                return
            end
    
            % sort into segments
            for tt = 1:length(trackROI(vid,arena).start)
                MT = nan(nframes,2);
                roi = trackROI(vid,arena).start(tt):trackROI(vid,arena).stop(tt);
                MT(roi,:) = [X(roi),Y(roi)];
                trackROI(vid,arena).tracks(:,:,end+1) = MT;
            end
        end
    end
    
    % % Remove the column of zeros 
    % for arena = 1:4
    %     X = squeeze(trackROI(vid,arena).tracks(:,1,:));
    %     loc = sum(X,1,'omitnan')==0;
    %     trackROI(vid,arena).tracks(:,1,:)
    % end
    
    % Visual confirmation that this roughly maps onto the number of flies that
    % are 'counted' 
    % fig = figure;
    % hold on
    for arena = 1:4
        ntrackedlines = sum(~isnan(squeeze(trackROI(vid,arena).tracks(:,1,:))),2);
    %     plot(ntrackedlines,'color', arenaData(arena).color)
        trackROI(vid,arena).ntrackedlines = ntrackedlines;
    end
        

    % Make a speed array
    for arena = 1:4
        X = squeeze(trackROI(vid,arena).tracks(:,1,:));
        Y = squeeze(trackROI(vid,arena).tracks(:,2,:));
        loc = sum(X,1,'omitnan')==0;
        X(:,loc) = []; %remove a column of all zeros
        Y(:,loc) = []; %remove a column of all zeros

        % calculate speed: 
        dX = diff(X,1,1);
        dY = diff(Y,1,1);
        speed = sqrt(dX.^2 + dY.^2);
        speed = (speed./pix2mm)*expData.parameters.FPS;
        trackROI(vid,arena).speed = speed;
        
        % extract the avg and err of speed for a given frame
        speedAvg = mean(speed,2,'omitnan');
        speedErr = std(speed,0,2,'omitnan');
    
        trackROI(vid,arena).speedAvg = [speedAvg; nan]; % nan for the final frame whose speed we can't calculate
        trackROI(vid,arena).speedErr = [speedErr; nan];

%         figure; histogram(reshape(speed,1,numel(speed)))
%         figure; plot(speedAvg)

    end
end

initial_vars{end+1} = 'trackROI';
clearvars('-except',initial_vars{:})
close all

%% FIGURE: plot the tracks for a given video 
if essentialfigs == false
    fig = figure;
    for arena = 1:4
        hold on
        vid = 1;
        z = trackROI(vid, arena).tracks;
        for ii = 1:size(z,3)
            plot(z(:,1,ii),z(:,2,ii),'linewidth', 0.5)
        end
    end
    axis square
    formatFig(fig,true);
    set(gca, 'XColor', 'k','YColor','k')
    
    save_figure(fig, [analysisDir expName ' track overlay video' num2str(vid)], '-png', autoSave,true,'-r80');
    clearvars('-except',initial_vars{:})
end

%% FIGURE: Avg length of a trace for each video
if essentialfigs == false
    fig = figure; hold on
    for arena = 1:4
        for vid = 1:nvids
            dummy = squeeze(trackROI(vid, arena).tracks(:,1,:));
            X = sum(~isnan(dummy),1);
            scatter(vid*ones(1,numel(X)),X,25,arenaData(arena).color,"filled")
            meanLength(vid,arena) = mean(X,2,'omitnan');
        end
    end
    for arena = 1:4
        plot(1:nvids,meanLength(:,arena),'color','w','linewidth', 2.5)
        plot(1:nvids,meanLength(:,arena),'color',arenaData(arena).color,'linewidth', 1.5)
    end
    
    xlabel('Video number')
    ylabel('track length (frames)')
    formatFig(fig,true);
    save_figure(fig, [analysisDir expName ' track length across videos'], '-png', autoSave,true,'-r80');
    
    clearvars('-except',initial_vars{:})
end

%% ANALYSIS: group the speed data over time

% find the maximum number of tracks for speed data
for arena = 1:4
    nMax = [];
    for vid = 1:nvids
        nMax = max([nMax, max(trackROI(vid, arena).ntrackedlines)]);
    end
    maxTrack(arena) = nMax;
end
disp('Max tracks for each arena:')
disp(maxTrack)

% combine all the speed points into a single structure
speed = struct;
for arena = 1:4
    speed(arena).raw = [];
    for vid = 1:nvids
        plotdata = (trackROI(vid, arena).speed);
        % move all the data points to an area that can be appended
        plotdata = sort(plotdata,2); 
        if size(plotdata,2)<maxTrack(arena)
            addData = nan(size(plotdata,1),maxTrack(arena));
            addData(:,1:size(plotdata,2)) = plotdata;
        else
            addData = plotdata(:,1:maxTrack(arena));
        end
        speed(arena).raw = [speed(arena).raw; nan(1,maxTrack(arena));addData];
    end
%     speed(arena).raw
end

% layout the average speed & error over time
for arena = 1:4
    speed(arena).avg = mean(speed(arena).raw,2,'omitnan');
    speed(arena).err = std(speed(arena).raw,0,2,'omitnan')./sqrt(nflies(arena));
end

% define walking as flies moving faster than 1.5mm/s
walkLim = 1.5;
for arena = 1:4
    loc = speed(arena).raw>=walkLim;
    speed(arena).walkNum = sum(loc,2);
    loc = speed(arena).raw<walkLim;
    speed(arena).restNum = sum(loc,2);
end

% Save the speed data into a separate structure for the experiments
initial_vars{end+1} = 'speed';
clearvars('-except',initial_vars{:})


%% FIGURE: speed summary figure
sSpan = 180;
row = 5;
col = 1;
sb(1).idx = 1;
sb(2).idx = 2:5;

fig = figure; set(fig, 'pos', [1941 145 998 601])
    subplot(row, col, sb(1).idx)
    plot(occupancy.time, occupancy.temp,'color', 'w', 'LineWidth',1)
    ylabel('\circC')
    subplot(row, col, sb(2).idx)
    hold on
    for arena = 1:4
       try
           plot(occupancy.time, smooth(speed(arena).avg(1:end-1),sSpan,'moving'),'color', arenaData(arena).color,'LineWidth',1)
       catch
           plot(occupancy.time, smooth(speed(arena).avg,sSpan,'moving'),'color', arenaData(arena).color,'LineWidth',1)
       end
    end
    xlabel('time (min)')
    ylabel('speed (mm/s)')
formatFig(fig, true,[row,col],sb);

save_figure(fig, [analysisDir expName ' avg speed over time all arenas'], '-png', autoSave,true,'-r80');
clearvars('-except',initial_vars{:})

%% FIGURE: speed histogram
if essentialfigs == false
    speedMax = 20;
    
    fig = figure; set(fig, 'pos', [87 258 1230 720])
    for arena = 1:4
        subplot(2,2,arena)
        X = speed(arena).raw;
        h(arena).data = histogram(X,'EdgeColor','w');
        set(gca, 'YScale','log')
        h_line(10,'gold')
        v_line(20,'red')
        xlabel('speed (mm/s)')
        ylabel('frame count')
        % find the portion of frames above the allowable instance:
        raw = reshape(X,numel(X),1);
        raw(isnan(raw)) = [];
        percentHigh(arena) = (sum(raw>speedMax)/length(raw))*100;
        title(['Over limit: ' num2str(percentHigh(arena)), '%'])
        
    end
    formatFig(fig,true,[2,2]);
    
    % Quality control
    if any(percentHigh>1)
        switch questdlg('More than 1% of frames are above speed limit, continue?')
            case 'Yes'
                disp('Continuing with analysis')
            case 'No'
                return
            case 'Cancel'
                close(fig)
                return
        end
    end
    
    save_figure(fig, [analysisDir expName ' speed histogram across arenas'], '-png', autoSave,true,'-r80');
    clearvars('-except',initial_vars{:})
end

%% FIGURE: percentage of flies walking
if essentialfigs == false
    sSpan = 1;
    row = 5;
    col = 1;
    sb(1).idx = 1;
    sb(2).idx = 2:5;
    
    for arena = 1:4
    
        fig = figure; set(fig, 'pos', [2035 443 908 600])
            subplot(row, col, sb(1).idx)
                plot(occupancy.time, occupancy.temp,'color', 'w','linewidth', 2)
                ylabel('\circC')
                xlabel('time')
                set(gca,'TickDir','out')
                axis tight
            subplot(row, col, sb(2).idx)
                plotdata = [smooth(speed(arena).walkNum,sSpan,'moving'),...
                            smooth(speed(arena).restNum,sSpan,'moving')];
                y = area(plotdata);
                y(2).EdgeColor = 'none';
                y(2).FaceColor = Color('darkslategrey'); % resting
                y(1).EdgeColor = 'none';
                y(1).FaceColor = Color('springgreen'); % walking
            %labels
            axis tight
            ylabel('flies (#)')
            formatFig(fig, true,[row,col],sb);
            subplot(row, col, sb(2).idx)
            set(gca,'XColor','k','TickDir','out')
        
        save_figure(fig, [figDir 'Arena ' Alphabet(arena) '/' expName ' walking vs resting fly numbers' ], '-png',autoSave,true,'-r80');
    
    end
    
    clearvars('-except',initial_vars{:})
end

%% SAVE: data stored in each subfolder for arenas

% Save group data into combo folder: 
save([figDir 'analysis/' expName ' speed data.mat'],'speed', 'trackROI');

% Save into each group folder:
SPEED = speed;
for arena = 1:nArenas
    speed = SPEED(arena);
    speedTracks = trackROI(:,arena);
    save([figDir 'Arena ' Alphabet(arena) '/' expName ' speed data.mat'],'speed', 'speedTracks');
end

results = 'Saved Data';
display(['Processed ' folder ' ' expName])

