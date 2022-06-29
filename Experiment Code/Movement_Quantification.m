
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

vid = 1;
trackROI = struct;

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
    % 3) previous frame in the same arena...
    lastArena = false(nframes,4);
    for arena = 1:4
        lastArena(:,arena) = [true; arenaLoc(1:end-1,arena)];
    end
    
    % Are all the conditions met?
    validLoc = false(nframes,4);
    for arena = 1:4 
        validLoc(:,arena) = ~nanLoc & arenaLoc(:,arena) & lastPoint & lastArena(:,arena);

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
            break
        end
        % make sure that the start comes before each stop...
        if ~all(trackROI(vid,arena).start<trackROI(vid,arena).stop)
            warndlg('mismatched start and stop of a tracked segment')
            break
        end

        % sort into segments
        for tt = 1:length(trackROI(vid,arena).start)
            MT = nan(nframes,2);
            roi = trackROI(vid,arena).start:trackROI(vid,arena).stop;
            MT(roi,:) = [X(roi),Y(roi)];
            trackROI(vid,arena).tracks(:,:,end+1) = MT;
        end
    end

end


% Visual confirmation that this roughly maps onto the number of flies that
% are 'counted' 

% CList = {'yellow', 'pink', 'purple', 'teal', 'navy', 'green'};

fig = figure;
hold on
for arena = 1:4
    ntrackedlines = sum(~isnan(squeeze(trackROI(vid,arena).tracks(:,1,:))),2);
    plot(ntrackedlines,'color', arenaData(arena).color)
end
    

    










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





    



  
   
        










        
















