
%% WATCH VIDEO: watch clip of fly tracks
clearvars('-except',initial_vars{:})
close all
%  Data Structure: data(video).tracks(frame number, body nodes, X Y coordinates, indivual fly tracks)

% pull and select data
vid = 1;
roi = 1:1797;
flynums = [3];
% flynums = [2,10,16,22,33,41,50];
test = squeeze(data(vid).tracks(:,1,:,flynums));

% pull image
movieInfo = VideoReader([figDir,'\',expName,'_',num2str(vid),'.avi']); 
% PLOT IMAGE
fig = figure; set(fig, 'pos', [2040 499 814 733],'color', 'k'); % X-off, Y-off, width, height
nanLoc = [nan,nan];
for ii = roi
    %pull current data
    img = read(movieInfo,ii);
    plotdata = squeeze(data(vid).tracks(ii,1,:,:));
    
    imshow(img)
    axis tight square
    hold on
    
    scatter(plotdata(1,:),plotdata(2,:), 15,'w','filled')
    % plot the history of 5 random flies...
    for fly = 1:length(flynums)
        x = test(ii,1,fly);
        y = test(ii,2,fly);
        plot(test(roi(1):ii,1,fly),test(roi(1):ii,2,fly),'Color', Color('yellow'),'LineWidth',0.5)
        scatter(x,y,15,'r','filled')
        if isnan(x) 
            % save the nan locations so that all future slides can also include the points
            nanLoc = [nanLoc; test(ii-1,:,fly)];
        end
        scatter(nanLoc(:,1),nanLoc(:,2),20,Color('orange')) 
        
    end

%     % throw an orange dot for any nan locations
%     if any(plotdata())
%         
%         test(roi(1):ii,1,fly)
%     % TODO...
%     end
    pause(0.05)
end

%% FIGURE: plot the full trace for a given fly over the course of a movie...
clearvars('-except',initial_vars{:})
close all

% pull and select data
flynums = 2; %[4,5,6]; %[2,3];
CList = {'yellow', 'pink', 'purple', 'teal', 'navy', 'green'};

% pull image
movieInfo = VideoReader([figDir,'\',expName,'_',num2str(vid),'.avi']); 
img = read(movieInfo,1);

% PLOT IMAGE
fig = figure; set(fig, 'pos', [2040 499 814 733],'color', 'k'); % X-off, Y-off, width, height
imshow(img)
axis tight square
hold on
% PLOT MOVEMENT TRACE
for vid = 1:expData.parameters.numVids
    lastPoint = [nan,nan];
    for ii = 1:length(flynums)
        fly = flynums(ii);
        plotData = squeeze(data(vid).tracks(:,1,:,fly));
        plot(plotData(:,1),plotData(:,2),'Color', Color(CList{ii}),'LineWidth',0.5) %
        % find the location of the last point before a NaN
        nanLoc = isnan(plotData(:,1));
        lastPoint = find(diff(nanLoc)==1);
        scatter(plotData(lastPoint,1),plotData(lastPoint,2),20,Color('red'),'filled') 
    end
end

save_figure(fig, [analysisDir 'individual fly track ' num2str(flynums)], '-png');


%% SEGMENT TRACKS INTO SINGLE TRUSTWORTHY TRACKS
clearvars('-except',initial_vars{:})
close all


trackROI = struct;

for vid = 1:nvids
    % pull and select data
    ntracks = size(data(vid).tracks,4);
    flynums = 1:ntracks;
    % CList = {'yellow', 'pink', 'purple', 'teal', 'navy', 'green'};
    
    % SEGMENT TRACKS INTO BETWEEN NAN SECTIONS
    % pull image
    movieInfo = VideoReader([figDir,'\',expName,'_',num2str(vid),'.avi']); 
    img = read(movieInfo,1);
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

%% FIGURE: number of tracks for each arena across the whole experiment

fig = figure;
    for arena = 1:4
        subplot(2,2,arena)
        hold on
        for vid = 1:nvids
            y = trackROI(vid,arena).ntrackedlines;
            scatter(vid*ones(1, length(y)), trackROI(vid,arena).ntrackedlines, 15,arenaData(arena).color,'filled')
        end
        % Labels
        xlabel('video number')
        ylabel('track count')
        h_line(nflies(arena),'w')
    end
formatFig(fig,true,[2,2]);
save_figure(fig, [analysisDir 'track count by video'], '-png');

clearvars('-except',initial_vars{:})

%% FIGURE: plot the tracks for a given video 

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

save_figure(fig, [analysisDir 'track overlay video' num2str(vid)], '-png');

clearvars('-except',initial_vars{:})

%% FIGURE: Avg length of a trace for each video
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
save_figure(fig, [analysisDir 'track length across videos'], '-png');

clearvars('-except',initial_vars{:})

%% ANALYSIS: group the speed data over time

% find the maximum number of tracks for speed data
for arena = 1:4
    nMax = [];
    for vid = 1:nvids
        nMax = max([nMax, max(trackROI(vid, arena).ntrackedlines)]);
    end
    maxTrack(arena) = nMax;
end
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


occupancy.speed = speed;

clearvars('-except',initial_vars{:})

%% FIGURE: show speed over the course of the experiment

arena = 1;
SZ = 15;

fig = figure; 
    hold on
    Y = occupancy.speed.raw;
    X = repmat(occupancy.time,1,size(Y,2));
    scatter(X,Y,SZ,'w')
    formatFig(fig,true);
xlabel('time (min)')
ylabel('speed (mm/s)')

save_figure(fig, [analysisDir 'speed over time arena ' Alphabet(arena)], '-png');

clearvars('-except',initial_vars{:})


%% FIGURE: average speed
% vidSpeed = [];
% for vid = 1:nvids
%     for arena = 1:4
%         vidSpeed(vid,arena) = mean(trackROI(vid,arena).speedAvg,'omitnan');
%     end
% end
% 
% figure;
% plot(vidSpeed)

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
        plot(occupancy.time, smooth(occupancy.speed(arena).avg(:),sSpan,'moving'),'color', arenaData(arena).color,'LineWidth',1)
    end
    xlabel('time (min)')
    ylabel('speed (mm/s)')
formatFig(fig, true,[row,col],sb);

save_figure(fig, [analysisDir 'avg speed over time all arenas'], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: speed histogram
speedMax = 20;

fig = figure; set(fig, 'pos', [87 258 1230 720])
for arena = 1:4
    subplot(2,2,arena)
    X = occupancy.speed(arena).raw;
    h(arena).data = histogram(X,'EdgeColor','w');
    set(gca, 'YScale','log')
    h_line(10,'gold')
    v_line(20,'red')
    xlabel('speed (mm/s)')
    ylabel('frame count')
    % find the portion of frames above the allowable instance:
    raw = reshape(X,numel(X),1);
    raw(isnan(raw)) = [];
    percentHigh = (sum(raw>speedMax)/length(raw))*100;
    title(['Over limit: ' num2str(percentHigh), '%'])
end
formatFig(fig,true,[2,2]);

save_figure(fig, [analysisDir 'speed histogram across arenas'], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: Average speed with speed cap
speedMax = 20;
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
        
        X = occupancy.speed(arena).raw;
        X(X>speedMax) = nan;
        Y = mean(X,2,'omitnan');

        plot(occupancy.time, smooth(Y,sSpan,'moving'),'color', arenaData(arena).color,'LineWidth',1)
    end
    xlabel('time (min)')
    ylabel('speed (mm/s)')
formatFig(fig, true,[row,col],sb);

save_figure(fig, [analysisDir 'avg speed with cap over time all arenas'], '-png');
clearvars('-except',initial_vars{:})


%% FIGURE: overlay of avg speed with and without speed cap
speedMax = 20;
sSpan = 180;


fig = figure; set(fig, 'pos', [1941 145 998 601])
    for arena = 1:4
        subplot(2,2,arena)
        X = occupancy.speed(arena).raw;
        plot(occupancy.time, smooth(mean(X,2,'omitnan'),sSpan,'moving'),'color', 'w','LineWidth',1)
        
        hold on
        X(X>speedMax) = nan;
        Y = mean(X,2,'omitnan');
        plot(occupancy.time, smooth(Y,sSpan,'moving'),'color', arenaData(arena).color,'LineWidth',1)
    end
    xlabel('time (min)')
    ylabel('speed (mm/s)')
formatFig(fig, true,[2,2]);

save_figure(fig, [analysisDir 'avg speed comparision between speed cap'], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: how does speed compare to movement? Quite nicely! 

speedMax = 20;
sSpan = 180;

row = 5;
col = 1;
sb(1).idx = 1;
sb(2).idx = 2:5;

for arena = 1:4

    fig = figure; set(fig, 'pos', [1941 145 998 601])
        subplot(row, col, sb(1).idx)
        plot(occupancy.time, occupancy.temp,'color', 'w', 'LineWidth',1)
        ylabel('\circC')
        subplot(row, col, sb(2).idx)
        hold on
        
        % plot speed
        X = occupancy.speed(arena).raw;
        X(X>speedMax) = nan;
        Y = mean(X,2,'omitnan');
        plot(occupancy.time, smooth(Y,sSpan,'moving'),'color', Color('gold'),'LineWidth',1)
        
        % plot movement
        Z = arenaData(arena).occupancy.movement;
        plot(occupancy.time(1:end-1), smooth(Z,sSpan,'moving'),'color', Color('royalblue'), 'LineWidth',1)
        
        % labels
        xlabel('time (min)')
        ylabel('speed (mm/s)')
        legend({'Speed', 'Movement'},'TextColor','w', 'box', 'off')

        % Correlation coefficient:
        subplot(row, col, sb(1).idx)        
        loc = isnan(Y(2:end)) | isnan(Z);
        Z(loc) = [];
        Y(1) = [];
        Y(loc) = [];
        R = corrcoef(Y,Z);
        title(['Correlation Coefficient: ' num2str(R(2,1))])

        
        
    formatFig(fig, true,[row,col],sb);
    
    save_figure(fig, [analysisDir 'speed vs movement arena ' Alphabet(arena)], '-png');

end

clearvars('-except',initial_vars{:})

%% Determine the percentage of flies that are walking: (percent of those tracked)

% define walking as flies moving faster than 1.5mm/s
walkLim = 1.5;

for arena = 1:4
    loc = occupancy.speed(arena).raw>=walkLim;
    occupancy.speed(arena).walkNum = sum(loc,2);
    loc = occupancy.speed(arena).raw<walkLim;
    occupancy.speed(arena).restNum = sum(loc,2);
end

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
            plotdata = [smooth(occupancy.speed(arena).restNum,sSpan,'moving'),...
                        smooth(occupancy.speed(arena).walkNum,sSpan,'moving')];
            y = area(plotdata);
            y(1).EdgeColor = 'none';
            y(1).FaceColor = Color('purple'); % resting
            y(2).EdgeColor = 'none';
            y(2).FaceColor = Color('orange'); % walking
        %labels
        axis tight
        ylabel('flies (#)')
        formatFig(fig, true,[row,col],sb);
        subplot(row, col, sb(2).idx)
        set(gca,'XColor','k','TickDir','out')
    
    save_figure(fig, [analysisDir 'walking vs resting ' Alphabet(arena)], '-png');

end

clearvars('-except',initial_vars{:})


%% 



































%% Scratch pad area

for arena = 1:4
    plot(1:nvids,meanLength(:,arena),'color','w','linewidth', 2.5)
    plot(1:nvids,meanLength(:,arena),'color',arenaData(arena).color,'linewidth', 1.5)
end

xlabel('Video number')
ylabel('track length (frames)')
formatFig(fig,true);
save_figure(fig, [analysisDir 'track length across videos'], '-png');


clearvars('-except',initial_vars{:})


% fig = figure;
%     plot(speedAvg,'color', arenaData(arena).color,'linewidth', 2)

arena = 4;
xMax = max(max(trackROI(vid,arena).tracks(:,1,:))) + 50;
xMin = min(min(trackROI(vid,arena).tracks(:,1,:))) - 50;
yMax = max(max(trackROI(vid,arena).tracks(:,2,:))) + 50;
yMin = min(min(trackROI(vid,arena).tracks(:,2,:))) - 50;

fig = figure; 
for tt = 1:1797
   
    X = squeeze(trackROI(vid,arena).tracks(tt,1,:));
    Y = squeeze(trackROI(vid,arena).tracks(tt,2,:));
    scatter(X,Y, 20,'black', 'filled')
    xlim([xMin,xMax])
    ylim([yMin,yMax])
    pause(0.01)
end
close(fig)

fig = figure; hold on
tt = 345;
for arena = 1:4
    X = squeeze(trackROI(vid,arena).tracks(tt,1,:));
    Y = squeeze(trackROI(vid,arena).tracks(tt,2,:));
    scatter(X,Y,30,arenaData(arena).color, 'filled')
    disp(sum(~isnan(X)))
end
formatFig(fig,true);
axis square
save_figure(fig, [analysisDir 'track fly positions frame ' num2str(tt)], '-png');



X = squeeze(trackROI(vid,2).tracks(:,1,:));

sum_X = sum(X,1,'omitnan');

% how many identical tracks are there??

ident_X = diff(X,1,2); % if a column is all zeros, that means the track to the right is identical
sum_ident_X = sum(ident_X,1,'omitnan'); % if the column total is zero, the track can be deleted as a duplicate

sum(sum_ident_X==0)





figure;
plot(squeeze(max(trackROI(vid,arena).tracks(:,2,:))))


save_figure(fig, [analysisDir 'corrected fly track numbers'], '-png');


save_figure(fig, [analysisDir 'track fly positions frame 345'], '-png');


% 
% 
% % PLOT IMAGE
% fig = figure; set(fig, 'pos', [2040 499 814 733],'color', 'k'); % X-off, Y-off, width, height
% imshow(img)
% axis tight square
% hold on
% % PLOT MOVEMENT TRACE
% for vid = 1:expData.parameters.numVids
%     lastPoint = [nan,nan];
%     for ii = 1:length(flynums)
%         fly = flynums(ii);
%         plotData = squeeze(data(vid).tracks(:,1,:,fly));
%         plot(plotData(:,1),plotData(:,2),'Color', Color(CList{ii}),'LineWidth',0.5) %
%         % find the location of the last point before a NaN
%         nanLoc = isnan(plotData(:,1));
%         lastPoint = find(diff(nanLoc)==1);
%         scatter(plotData(lastPoint,1),plotData(lastPoint,2),20,Color('red'),'filled') 
%     end
% end
% 
% save_figure(fig, [analysisDir 'individual fly track ' num2str(flynums)], '-png');







%% RECORD VIDEO OF SELECTED FLY TRACES THROUGH TIME

clearvars('-except',initial_vars{:})
close all
videoName = 'G:\My Drive\Jeanne Lab\DATA\06.22.2022\analysis\Track 3 vid 1';

% pull and select data
vid = 1;
roi = 1:1797;
flynums = [3];
test = squeeze(data(vid).tracks(:,1,:,flynums));

% pull video information
movieInfo = VideoReader([figDir,'\',expName,'_',num2str(vid),'.avi']); 
buffer_frame = read(movieInfo,1);
buffer_frame(:,:,:) = 1;

% Set up video parameters
fig = figure; set(fig, 'pos', [2040 499 814 733],'color', 'k'); % X-off, Y-off, width, height
nanLoc = [nan,nan];

    hold on
    v = VideoWriter(videoName, 'Uncompressed AVI');
    v.FrameRate = 30;
    open(v);
    imshow(buffer_frame); axis tight square
    currAxes.Visible = 'off';
    % Collect image info for movie file
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')

    % Draw video frames...
    for ii = roi
        set(fig, 'color','k'); 
        %pull current data
        img = read(movieInfo,ii);
        plotdata = squeeze(data(vid).tracks(ii,1,:,:));
        imshow(img)
        axis tight square
        hold on
        scatter(plotdata(1,:),plotdata(2,:), 15,'w','filled')
        % plot the history of 5 random flies...
        for fly = 1:length(flynums)
            x = test(ii,1,fly);
            y = test(ii,2,fly);
            plot(test(roi(1):ii,1,fly),test(roi(1):ii,2,fly),'Color', Color('yellow'),'LineWidth',0.5)
            scatter(x,y,15,'r','filled')
            if isnan(x) 
                % save the nan locations so that all future slides can also include the points
                nanLoc = [nanLoc; test(ii-1,:,fly)];
            end
            scatter(nanLoc(:,1),nanLoc(:,2),20,Color('orange'))     
        end
    


        currAxes.Visible = 'off';
        % Collect image info for movie file
        f = getframe(fig);
        writeVideo(v, f)  
        clf('reset')
    end

    
close(v)
close all
fprintf('\n Video Saved! \n')  



%% RECORD VIDEO OF ALL FLY TRACES THROUGH TIME for a given time history and for a single video

clearvars('-except',initial_vars{:})
close all
vid = 56;
videoName = [figDir 'analysis\All tracks vid ' num2str(vid)];

% pull video information
movieInfo = VideoReader([figDir,'\',expName,'_',num2str(vid),'.avi']); 
buffer_frame = read(movieInfo,1);
buffer_frame(:,:,:) = 1;
vidLength = movieInfo.NumFrames;

% pull and select data
traceLength = 30; % how many previous points to demo for each fly track
roi = 1:vidLength-2;
% roi = 1:1795;
% rawdata = trackROI(vid,arena).tracks;
walkLim = 1.5;

% Set up video parameters
fig = figure; set(fig, 'pos', [2040 499 814 733],'color', 'k'); % X-off, Y-off, width, height
nanLoc = [nan,nan];

    hold on
    v = VideoWriter(videoName, 'Uncompressed AVI');
    v.FrameRate = 30;
    open(v);
    imshow(buffer_frame); axis tight square
    currAxes.Visible = 'off';
    % Collect image info for movie file
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')

    % Draw video frames...
    for ii = roi
        set(fig, 'color','k'); 
        %pull current data
        img = read(movieInfo,ii);
%         plotdata = squeeze(data(vid).tracks(ii,1,:,:));
        imshow(img)
        axis tight square
        hold on

        % find the history plot region of interest
        if ii<=traceLength
            pROI = 1:ii;
        else
            pROI = ii-traceLength:ii;
        end

        for arena = 1:4
            rawdata = trackROI(vid,arena).tracks;
            % plot the current position of the flies
            x = squeeze(rawdata(ii,1,:));
            y = squeeze(rawdata(ii,2,:));
            % plot all fly positions
            scatter(x,y,30,arenaData(arena).color,'filled')
            % Check the walking speed of the flies and color code if they
            % are 'walking'
            loc = [false,trackROI(vid,arena).speed(ii,:)>=walkLim];
            scatter(x(loc'),y(loc'),60,'red','filled')
            scatter(x(loc'),y(loc'),30,'yellow','filled')
           
            % plot the past positions of the flies
            if length(pROI)>1
                plot(squeeze(rawdata(pROI,1,:)),squeeze(rawdata(pROI,2,:)),'color', 'w','LineWidth',0.5)

            end
        end
        
%         pause(0.01)
%     end
       
        currAxes.Visible = 'off';
        % Collect image info for movie file
        f = getframe(fig);
        writeVideo(v, f)  
        clf('reset')
    end

    
close(v)
close all
fprintf('\n Video Saved! \n')  


    



  
   
        










        
















