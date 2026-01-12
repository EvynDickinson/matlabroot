
% TransferProcessedDataToServer

%% 
% pos (frame, body points, XY)
% body points (head, center, abdomen, right wing, left wing)

%% LOAD DATA

% TODOs
% identify periods of temperature cooling, warming, binned temp.
% create predictions of courtship periods
% how long are they in courtship
% distance to food
% body angle between flies
% outer ring occupancy
% sleep
% make a plot of all behaviors that meet everything but the time minimum &
% plot them with the 'official' behavior periods (raster plot style)
% circling behavior
% courtship index  
% courtship index vs temperature

% Prep data
if ~strcmp(questdlg('Is data from 5.1.1 already loaded?'),'Yes')
    clear; clc;
    startloc = getDataPath(6,0);
    if isempty(startloc)
        return
    end
    baseFolder_1 = [startloc,'Trial Data/'];
    trialDir_1 = selectFolder(baseFolder_1); 
    baseDir_1 = [baseFolder_1, trialDir_1{:} '/']; % full folder directory for that trial
    
    load([baseDir_1, 'post-5.1.1 data.mat']) % load the parameters and temp table
    baseDir = baseDir_1;
    baseFolder = baseFolder_1;
    trialDir = trialDir_1;
    clear("baseDir_1","trialDir_1","baseFolder_1")
    figDir = [baseDir,'Figures/']; 
    if ~exist(figDir, 'dir')
        mkdir(figDir)
    end
    disp('data loaded')
end

% Experiment parameters
blkbnd = true;
fig_type = '-png';

% Initial variables
initial_var = who; % who = all variables created so far
initial_var{end+1} = 'initial_var';
initial_var{end+1} = 'well';
initial_var{end+1} = 'path';
            
disp_fig = false; % display baseline figures?
initial_var{end+1} = 'disp_fig';
            
%% FIGURE: Compare wing angles within M and between M and F
clearvars('-except',initial_var{:})
% Compare male L and R wing angles
fig = getfig('comparing  male L and R wing angles',true, [1032 300]); 
    hold on
    plot(time, data(M).wingangle(:,1),'color', Color('blue'),'linewidth', 1) % left wing
    plot(time, data(M).wingangle(:,2),'color', Color('gold'),'linewidth', 1) % right wing
% format figure
xlabel('time (s)')
ylabel('wing angle (\circ)')
formatFig(fig, blkbnd);
set(gca, 'xcolor',foreColor)
save_figure(fig,[figDir 'male wing angle'],fig_type);

% Compare male and female wing angles
fig = getfig('Compare M and F wing angles',true, [1032 300]); 
hold on
    plot(time, data(M).wingspread,'color', Color('dodgerblue'),'linewidth', 1) % male wing spread
    plot(time, data(F).wingspread,'color', Color('deeppink'),'linewidth', 1) % female wing spread
% format figure
xlabel('time (s)')
ylabel('wing angle (\circ)')
formatFig(fig, blkbnd);
save_figure(fig,[figDir 'M and F wing angles'],fig_type);

%% FIGURE: Distance to food histogram
clearvars('-except',initial_var{:})
fig = figure; 
    histogram(T.dist2food(:,M),'FaceColor',data(M).color,'FaceAlpha',0.8)
    hold on
    histogram(T.dist2food(:,F),'FaceColor',data(F).color,'FaceAlpha',0.8)
    xlabel('distance to food (mm)')
    formatFig(fig, blkbnd);
    set(gca,'ycolor', 'none')
    save_figure(fig,[figDir 'distance to food histogram M vs F'],fig_type);

% flies ON food
% T.FlyOnFood = T.dist2food<=well.R; % fly head must be within the food circle
% disp('Flies on food: M & F')
% sum(T.FlyOnFood)

%% FIGURE: Fly distance to food and flies on food timecourse
r = 4;
c = 1;
sb(1).idx = 1; % temperature
sb(2).idx = 2:4; % distance to food
lw = 2;
sSpan = 5*fps;
spike_H = 2;    %height of each raster
trial_space = 1;    %gap between trial lines
spike_W = 0.5;    % raster line width

fig =  getfig('',1); 
    subplot(r,c,sb(1).idx); hold on
        plot(time, T.temperature, 'color', foreColor,'linewidth', lw)
        ylabel('temp (\circC)')
    subplot(r,c,sb(2).idx); hold on
        for sex = 1:2
            % plot the smoothed average distance
            plot(time, smooth(T.dist2food(:,sex),sSpan, 'moving'), 'color', data(sex).color,'linewidth', lw)
        end
        y_base = rangeLine(fig, 5);
        for sex = 1:2 
            % plot raster of when the fly is on the food
            x = time(T.FlyOnFood(:,sex));
            X = [x';x'];
            Y = repmat([y_base;y_base+spike_H],[1,size(X,2)]);           
            plot(X,Y,'color',data(sex).color,'linewidth',spike_W)
            
            % Update y-offset
            y_base = y_base+spike_H+trial_space;

        end
        ylabel('distance to food (mm)')
        xlabel('time (min)')

    formatFig(fig, blkbnd, [r,c],sb);
    save_figure(fig,[figDir 'Distance to food with points on food M vs F'],fig_type);

%% FIGURE: M body position during wing extension
% wing angle > 60 deg
% L wing angle in appropriate quadrants
% R wing angle in appropriate quandrants
% matrix for yes courtship wing extension, no not courtship
% diff to see if extension lasts > 1sec

clearvars('-except',initial_var{:})

tic
frames = find(T.court_chase==1);
targetn = 10;
skip = floor(length(frames)/targetn);
if skip <1
    skip = 1;
end

% demo images of wing extension:
frames = find(T.wing_ext==1);
if isempty(frames)
    disp('no bouts of wing extension found')
    return
else
    endidx = (find(diff(frames)>1)); % where the does first 'extension' bout end?
    try 
        % find frame where extension bout ends
        roi = frames(1):frames(endidx(1));
    catch
        % if extension bout doesn't end by the end, use last frame
        roi = frames(1):frames(end);
    end
    
    fig = getfig('',1); hold on
    % plot the female fly
    for ff = 1:length(roi)
        frame = roi(ff);
        for sex = 1:2
            x = data(sex).rawX(frame,:);
            y = data(sex).rawY(frame,:);
            kolor = data(sex).color;
            plotFlySkeleton(fig, x,y,kolor,false);
            scatter(x,y, 15, Color('grey'),'filled')
            scatter(x(body.head),y(body.head), 35, Color('yellow'),'filled', 'Marker','^')
            % scatter(x(body.center),y(body.center), 15, Color('grey'),'filled')
            
            if sex==1
                scatter(x(body.left_wing),y(body.left_wing),35,foreColor,'filled')
            end
        end
    end
    formatFig(fig,blkbnd);
    set(gca, 'xcolor', 'none', 'ycolor', 'none')
    save_figure(fig, [figDir, 'Wing extension example 1'],'-png');
end

% Pull point locations that will be plotted
if strcmp(parameters.protocol,'high_res_LTS_35-15')
    skip = 100;
else
    skip = 20;
end
zoom = [-250,250];

% screening = close_dist;

fig = getfig('Male positions during wing extension',1,[1075 871]);
hold on
% Plot all male body positions
kolor = Color('grey');
x = mX(1:skip:end,[body.head,body.center]);
y = mY(1:skip:end,[body.head,body.center]);
plot(x',y','color',kolor)
scatter(x(:,1),y(:,1),15,kolor,"filled","^")
% axis equal square
formatFig(fig,blkbnd);
set(gca,'XColor','none','YColor','none')

% Plot the all male body positions under chase instances
kolor = Color('lime');
x = mX(T.wing_ext,[body.head,body.center]);
y = mY(T.wing_ext,[body.head,body.center]);
plot(x',y','color',kolor)
scatter(x(:,1),y(:,1),15,kolor,"filled","^")

% % Screening
% kolor = Color('green');
% x = mX(screening,[body.head,body.center]);
% y = mY(screening,[body.head,body.center]);
% plot(x',y','color',kolor)
% scatter(x(:,1),y(:,1),15,kolor,"filled","^")

% Plot female body
x = fX(1:skip:end,[body.head,body.center,body.abdomen]);
y = fY(1:skip:end,[body.head,body.center,body.abdomen]);
plot(x',y','color',foreColor, 'LineWidth', 2)
% xlim([-2000,1500]); ylim([-2000,2000])

% Format figure
axis  equal square
h_line(0,'gray',':',2)
v_line(0,'grey',':',2)
xlim(zoom)
ylim(zoom)
formatFig(fig,blkbnd);
set(gca,'XColor','none','YColor','none')

formatFig(fig,blkbnd);
% rectangle('Position',[zoom(1),zoom(1) sum(abs(zoom)) sum(abs(zoom))],'edgecolor',foreColor,'linewidth', 1)
% save_figure(fig, 'G:\My Drive\Jeanne Lab\Presentations\Data Presentation 12.6.2024\All Chase Positions',fig_type);

% Save figure
save_figure(fig,[figDir 'wing extension positions M fly'],fig_type);



%% FIGURE: M body positions during chase
% Pull point locations that will be plotted
tic
frames = find(T.court_chase==1);
targetn = 10;
skip = floor(length(frames)/targetn);
if skip <1
    skip = 1;
end

% if strcmp(parameters.protocol,'high_res_LTS_35-15')
%     skip = 400;
% else
%     skip = 20;
% end
zoom = [-250,250];

% screening = close_dist;

if isempty(frames)
    disp('no bouts of chase found')
    return
else
    fig = getfig('Male positions during chase',1,[1075 871]);
    hold on
    % Plot all male body positions
    kolor = Color('grey');
    x = mX(1:skip:end,[body.head,body.center]);
    y = mY(1:skip:end,[body.head,body.center]);
    plot(x',y','color',kolor)
    scatter(x(:,1),y(:,1),15,kolor,"filled","^")
    axis equal square
    formatFig(fig,blkbnd);
    set(gca,'XColor','none','YColor','none')
    
    % Plot the all male body positions under chase instances
    kolor = Color('lime');
    x = mX(T.court_chase,[body.head,body.center]);
    y = mY(T.court_chase,[body.head,body.center]);
    plot(x',y','color',kolor)
    scatter(x(:,1),y(:,1),15,kolor,"filled","^")
    
    % % Screening
    % kolor = Color('green');
    % x = mX(screening,[body.head,body.center]);
    % y = mY(screening,[body.head,body.center]);
    % plot(x',y','color',kolor)
    % scatter(x(:,1),y(:,1),15,kolor,"filled","^")
    
    % Plot female body
    x = fX(1:5,[body.head,body.center,body.abdomen]);
    y = fY(1:5,[body.head,body.center,body.abdomen]);
    plot(x',y','color',foreColor, 'LineWidth', 2)
    % xlim([-2000,1500]); ylim([-2000,2000])
    
    % Format figure
    axis  equal square
    h_line(0,'gray',':',2)
    v_line(0,'grey',':',2)
    xlim(zoom)
    ylim(zoom)
    formatFig(fig,blkbnd);
    set(gca,'XColor','none','YColor','none')
    
    formatFig(fig,blkbnd);
    % rectangle('Position',[zoom(1),zoom(1) sum(abs(zoom)) sum(abs(zoom))],'edgecolor',foreColor,'linewidth', 1)
    % save_figure(fig, 'G:\My Drive\Jeanne Lab\Presentations\Data Presentation 12.6.2024\All Chase Positions',fig_type);
    
    % Save figure
    save_figure(fig,[figDir 'chase positions M fly'],fig_type);
end
toc



%% FIGURE: Chase overlaid on arena with temp timecourse
%  OLD CODE--- NOW INCORPORATED INTO THE NEXT SECTION
% switch questdlg(['Show all ' num2str(size(m.chaseroi,1)) ' instances of chase?'])
%     case 'Yes'
%     case 'No'
%         return
%     case 'Cancel'
%         return
% end
% % Subplots
% r = 5;
% c = 1;
% sb(1).idx = 1;
% sb(2).idx = 2:5;
% 
% % For each bout of chase, plot M and F movement on arena image
% for i = 1:size(m.chaseroi,1)
%     % Frame numbers at start and end of each chase bout
%     plotroi = m.chaseroi(i,1):m.chaseroi(i,2);
%     % Identify which frame (last in each bout) and video to pull 
%     frame = m.chaseroi(i,2);
%     vidnum = T.vidNums(frame);
% 
%     % Pull and read in video
%     vidpath = [path, parameters.date, '\', parameters.videoName, '\compiled_video_', num2str(vidnum), '.avi'];
%     movieInfo = VideoReader(vidpath);
%     demoImg = (read(movieInfo,T.vidFrame(frame)));
%     img = imadjust(demoImg,[72/255, 215/255]);
% 
%     % Scatter point size and linewidth
%     sz = 10;
%     lw = 1;
%     % Plot body position at every __ frame
%     frame_skip = 1;
% 
%     % Create x and y limit variable for zoomed figures
%     [xlimits, ylimits] = deal([]);
% 
%     % Plot body positions over chase bout
%     fig = getfig(' ', true, [759 900]); 
%     % 1) Temperature
%     subplot(r, c, sb(1).idx)
%         hold on
%         % Plot temp timecourse
%         x = time;
%         y = T.temperature;
%         plot(x,y,'color', foreColor,'LineWidth', lw)
%         % Plot vertical lines at each chase bout (orange = current bout shown)
%         v_line(time(m.chaseroi(:)),'teal',':',2)
%         v_line(time(m.chaseroi(i,:)),'orange',':',2)
%         % Axes labels and limits
%         xlabel('time (s)')
%         ylabel('\circC')
%         xlim([0,time(end)])
%     % 2) Body positions overlaid on arena image 
%     subplot(r, c, sb(2).idx)
%         % Plot arena image
%         imshow(img)
%         % Plot fly centers over the course of the chase bout
%         ROI = plotroi(1:frame_skip:end);
%         % Male positions
%         x1 = m.pos(ROI, 1,1);
%         y1 = m.pos(ROI, 1,2);
%         hold on
%         scatter(x1,y1,sz,Color('dodgerblue'), "filled")
%         % Female positions
%         x2 = f.pos(ROI, 1,1);
%         y2 = f.pos(ROI, 1,2);
%         scatter(x2,y2,sz,Color('deeppink'), "filled")
% 
%     % Format figure
%     formatFig(fig, blkbnd, [r,c], sb);
%     % lineROI = drawline(gca);
% 
%     % Save figure
%     save_figure(fig,[figDir 'Chase_', num2str(i), ' ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type,false, false);
% 
%     % ---------------------------------------- Zoom in on arena and display skeletons -----------------------------------
%         % Zoom in on the flies
%         xlimits(1) = min([x1; x2; m.pos(frame, :,1)'; f.pos(frame, :,1)']);
%         xlimits(2) = max([x1; x2; m.pos(frame, :,1)'; f.pos(frame, :,1)']);
%         ylimits(1) = min([y1;y2; m.pos(frame, :,2)'; f.pos(frame, :,2)']);
%         ylimits(2) = max([y1;y2; m.pos(frame, :,2)'; f.pos(frame, :,2)']);
%         buff = 50;
%         xlim([xlimits(1)-buff, xlimits(2)+buff])
%         ylim([ylimits(1)-buff, ylimits(2)+buff])
% 
%         % Overlay the current body position of the male fly
%         x = m.pos(frame, :,1);
%         y = m.pos(frame, :,2);
%         scatter(x,y,sz,Color('black'), "filled")
%         skeleton = [1,2; 2,3; 2,4; 2,5];
%         for ii = 1:size(skeleton,1)
%             plot(x(skeleton(ii,:)),y(skeleton(ii,:)), 'color', Color('black'),'LineWidth', 1.5)
%         end
% 
%         % Overlay the current body position of the female fly
%         x = f.pos(frame, :,1);
%         y = f.pos(frame, :,2);
%         scatter(x,y,sz,Color('black'), "filled")
%         skeleton = [1,2; 2,3; 2,4; 2,5];
%         for ii = 1:size(skeleton,1)
%             plot(x(skeleton(ii,:)),y(skeleton(ii,:)), 'color', Color('black'),'LineWidth', 1.5)
%         end
% 
%     % Save figure
%     save_figure(fig,[figDir 'Chase_zoom_', num2str(i), ' ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type);
% end

%% FIGURES: Visualize chasing overlaid on arena image, with side zoom on other parameters 
clearvars('-except',initial_var{:})
chaseDir = createFolder([figDir, 'Chase Figures\']);
fullArena_R = 29.5;
% TODO (3/11) update this to have an option to just display the tracks but not the
% images (espcially useful for trials with large numbers of chase bouts) -- this
% would also allow this section to be run without access to the raw videos

nChaseBouts = size(m.chaseroi,1);
disp(['There are ' num2str(nChaseBouts) ' chase bouts in this trial'])
switch questdlg('How do you want to visualize the chase bouts?', '', 'Image overlay', 'Schematic Summary', 'Both','Schematic Summary')
    case 'Image overlay'
        dType_1 = true;
        dType_2 = false;
    case 'Schematic Summary'
        dType_1 = false;
        dType_2 = true;
    case 'Both'
        dType_1 = true;
        dType_2 = true;
end

timebuff = 3; % time buffer before and after chase (in seconds)
timebuff = timebuff/60; % set to minutes

if dType_1
    % Subplots
    r = 7;
    c = 9;
    sb(1).idx = 2:4; %  temperature
    sb(2).idx = [19:23, 28:32, 37:41, 46:50, 55:59]; % arena image
    sb(3).idx = 6:9; % zoom in temp 
    sb(4).idx = [15:18, 24:27]; % distance between flies
    sb(5).idx = [33:36, 42:45]; % speed correlation
    sb(6).idx = [51:54, 60:63]; % male wing angles
    
    
    
    % Plot body movements and timecourses for each chase bout
    for i = 1:nChaseBouts
        % Establish x limits for non-temp timecourses
        xlimit = [time(m.chaseroi(i,1))-timebuff,time(m.chaseroi(i,2))+timebuff];
    
        % Frame numbers at start and end of each chase bout
        plotroi = m.chaseroi(i,1):m.chaseroi(i,2);
        % Identify which frame (last in each bout) and video to pull
        frame = m.chaseroi(i,2);
        vidnum = T.vidNums(frame);
        
        vidpath = [path, parameters.date, '\', parameters.videoName, '\compiled_video_', num2str(vidnum), '.avi'];
        movieInfo = VideoReader(vidpath); %read in video
        demoImg = (read(movieInfo,T.vidFrame(frame)));
        img = imadjust(demoImg,[72/255, 215/255]); % adjust the contrast of the image
    
        % Scatter point size and linewidth
        sz = 10;
        lw = 2;
        % Plot body position at every __ frame
        frame_skip = 1;
      
        % Plot body positions and timecourses
        fig = getfig(' ', false, [2100 1065]); 
    
        % 1) Temperature 
        subplot(r, c, sb(1).idx)
            hold on
            % Plot temp timecourse
            x = time;
            y = T.temperature;
            plot(x,y,'color', foreColor,'LineWidth', lw)
            % Plot vertical lines at each chase bout (orange = current bout shown)
            v_line(time(m.chaseroi(:)),'teal',':',2)
            v_line(time(m.chaseroi(i,:)),'orange',':',2)
            % Axes labels and limits
            xlabel('time (min)')
            ylabel('(\circC)')
            xlim([0,time(end)])
    
        % 2) Image of flies in the arena
        subplot(r,c,sb(2).idx)
            imshow(img)
            % Plot fly centers over the course of the chase bout
            ROI = plotroi(1:frame_skip:end);
            % Male positions
            x1 = m.pos(ROI, 1,1);
            y1 = m.pos(ROI, 1,2);
            hold on
            scatter(x1,y1,sz,Color('dodgerblue'), "filled")
            % Female positions
            x2 = f.pos(ROI, 1,1);
            y2 = f.pos(ROI, 1,2);
            scatter(x2,y2,sz,Color('deeppink'), "filled")
            set(gca, 'xcolor', 'none', 'ycolor', 'none')
    
        % 3) Mini temp zoom in 
        subplot(r,c,sb(3).idx)
            x = time;
            y = T.temperature;
            plot(x,y,'color', foreColor,'LineWidth', lw)
            ylabel('\circC')
    
        % 4) Zoom of IFD
         subplot(r,c,sb(4).idx)
            plot(time,T.IFD,'color', foreColor,'LineWidth', lw)
            % Axes labels and limits
            ylabel('IFD (mm)')
    
        % 5) Fly speed correlation
        subplot(r,c,sb(5).idx); 
            tSpan = 3; % time window (sec)
            pSpan = tSpan * fps; % frames in the sliding window
            % Plot speed corr timecourse
            y = runningCorrelation([m.speed,f.speed], pSpan);
            y = smooth(y,pSpan,'moving');
            offset = ceil(pSpan/2);
            x = time(offset:offset+length(y)-1)';
            plot(x, y, 'color', foreColor,'LineWidth', lw)
            ylabel('Speed Corr')
            h_line(0,'grey','--',1)
    
        % 6) Male wingspread 
        subplot(r,c,sb(6).idx); hold on
            plot(time,data(M).wingangle(:,1),'color', Color('dodgerblue'),'LineWidth', lw)
            plot(time,data(M).wingangle(:,2),'color', Color('gold'),'LineWidth', lw)
            h_line(50,'grey', ':',1)
            % Axes labels and limits
            ylabel('M wing angle (\circ)')
            xlabel('time (min)')
    
        % Format figure    
        formatFig(fig, blkbnd, [r,c], sb);
        for subby = 3:5
            subplot(r,c,sb(subby).idx)
            set(gca, 'xcolor', 'none')
        end
        % set x limits and the start and stop time of the 'official' chase period
        for subby = 3:6
          subplot(r,c,sb(subby).idx)
            v_line(time(m.chaseroi(i,:)),'orange',':',2)
            xlim(xlimit)
        end
        
        % Save figure
        save_figure(fig,[chaseDir 'Chase Bout_', num2str(i), ' from ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type,false, false);
    
        % ---------------------------------------- Zoom in on arena and display skeletons -----------------------------------
    
            % Zoom in on the flies
            subplot(r,c,sb(2).idx); hold on
            [xlimits, ylimits] = deal([]);% Create x and y limit variable for zoomed figures
            xlimits(1) = min([x1;x2]);
            xlimits(2) = max([x1;x2]);
            ylimits(1) = min([y1;y2]);
            ylimits(2) = max([y1;y2]); 
            buff = 50;
            xlim([xlimits(1)-buff, xlimits(2)+buff])
            ylim([ylimits(1)-buff, ylimits(2)+buff])
    
            % Overlay the current body position of the male fly
            x = m.pos(frame, :,1);
            y = m.pos(frame, :,2);
            scatter(x,y,sz,Color('black'), "filled")
            skeleton = [1,2; 2,3; 2,4; 2,5];
            for ii = 1:size(skeleton,1)
                plot(x(skeleton(ii,:)),y(skeleton(ii,:)), 'color', Color('black'),'LineWidth', 1.5)
            end
    
            % Overlay the current body position of the female fly
            x = f.pos(frame, :,1);
            y = f.pos(frame, :,2);
            scatter(x,y,sz,Color('black'), "filled")
            skeleton = [1,2; 2,3; 2,4; 2,5];
            for ii = 1:size(skeleton,1)
                plot(x(skeleton(ii,:)),y(skeleton(ii,:)), 'color', Color('black'),'LineWidth', 1.5)
            end
    
        % Save figure
        save_figure(fig,[chaseDir 'Chase Bout Zoom_', num2str(i), ' from ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type);
    end
end

if dType_2
    r = floor(sqrt(nChaseBouts));
    c = ceil(nChaseBouts/r);

     % Scatter point size and linewidth
    sz = 10;
    lw = 2;
    % Plot body position at every __ frame
    frame_skip = 1;

    fig = getfig('',1); 
    % Plot body movements and timecourses for each chase bout
    for i = 1:nChaseBouts
        subplot(r,c,i); hold on
        % Establish x limits for non-temp timecourses
        xlimit = [time(m.chaseroi(i,1))-timebuff,time(m.chaseroi(i,2))+timebuff];
        % Frame numbers at start and end of each chase bout
        plotroi = m.chaseroi(i,1):m.chaseroi(i,2);
        % Plot fly centers over the course of the chase bout
        ROI = plotroi(1:frame_skip:end);
        % Male positions
        x1 = m.pos(ROI, 1,1);
        y1 = m.pos(ROI, 1,2);
        hold on
        scatter(x1,y1,sz,Color('dodgerblue'), "filled")
        % Female positions
        x2 = f.pos(ROI, 1,1);
        y2 = f.pos(ROI, 1,2);
        scatter(x2,y2,sz,Color('deeppink'), "filled")
        % plot the arena circle: 
        foodWell = well.center(well.food_idx,:);
        centre = well.center(5,:);
        scatter(foodWell(1), foodWell(2),30,foreColor,'filled')
        viscircles(centre, fullArena_R/pix2mm,'color', foreColor,'linewidth', 0.5)
    end
    HCH = {'C', 'H', 'S'};
    % figure formatting: 
    formatFig(fig, blkbnd,[r,c]);
    for i = 1:nChaseBouts
        subplot(r,c,i); hold on
        set(gca, 'xcolor', 'none', 'ycolor', 'none')
        axis square
        idx = m.chaseroi(i,1);
        temptype = HCH([T.cooling(idx), T.warming(idx), T.hold(idx)]);
        title([temptype{:} ': ' num2str(T.temperature(idx))],'color', foreColor,'FontSize',10)
    end
    % SAVE FIGURE: 
    save_figure(fig,[chaseDir 'All chase bouts full arena schematic'], fig_type);

    % Zoomed in schematic summary version: 
    fig = getfig('',1); 
    % Plot body movements and timecourses for each chase bout
    for i = 1:nChaseBouts
        subplot(r,c,i); hold on
        % Establish x limits for non-temp timecourses
        xlimit = [time(m.chaseroi(i,1))-timebuff,time(m.chaseroi(i,2))+timebuff];
        % Frame numbers at start and end of each chase bout
        plotroi = m.chaseroi(i,1):m.chaseroi(i,2);  
        % Plot fly centers over the course of the chase bout
        ROI = plotroi(1:frame_skip:end);
        % Male positions
        x1 = m.pos(ROI, 1,1);
        y1 = m.pos(ROI, 1,2);
        hold on
        scatter(x1,y1,sz,Color('dodgerblue'), "filled")
        scatter(x1(end),y1(end),sz+5,foreColor, "filled")
        % Female positions
        x2 = f.pos(ROI, 1,1);
        y2 = f.pos(ROI, 1,2);
        scatter(x2,y2,sz,Color('deeppink'), "filled")
        scatter(x2(end),y2(end),sz+5,foreColor, "filled")
    end
    
    HCH = {'C', 'H', 'S'};
    % figure formatting: 
    formatFig(fig, blkbnd,[r,c]);
    for i = 1:nChaseBouts
        subplot(r,c,i); hold on
        set(gca, 'xcolor', 'none', 'ycolor', 'none')
        axis equal
        idx = m.chaseroi(i,1);
        temptype = HCH([T.cooling(idx), T.warming(idx), T.hold(idx)]);
        title([temptype{:} ': ' num2str(T.temperature(idx))],'color', foreColor,'FontSize',10)
    end
    % SAVE FIGURE: 
    save_figure(fig,[chaseDir 'All chase bouts zoomed in schematic'], fig_type);
end

%% FIGURE: (TODO) M body positions during circling
% Pull point locations that will be plotted
skip = 20;
zoom = [-250,250];

% find the longest circling behavior example: 
[~,idx] = max(smooth(T.circling_all,fps,'moving'));
roi = idx-fps:idx+fps;
% figure; plot(roi,T.circling_all(roi))

% roi = 81961:82022+fps; % uncomment for manual frame targeted display

% 
% frames = find(T.wing_ext==1);
% endidx = (find(diff(frames)>1)); % where the does first 'extension' bout end?
% roi = frames(1):frames(endidx(1));

fig = getfig('Male positions during circling',1); hold on
% plot the female fly
for ff = 1:length(roi)
    frame = roi(ff);
    for sex = 1:2
        x = data(sex).rawX(frame,:);
        y = data(sex).rawY(frame,:);
        kolor = data(sex).color;
        if T.circling_all(frame) && sex==M
            plotFlySkeleton(fig, x,y,foreColor,false);
            %plotFlySkeleton(fig, x,y,kolor,false); %plot all frames same color regardless of all vs 1sec category
        else
           plotFlySkeleton(fig, x,y,kolor,false);
        end
        scatter(x,y, 15, Color('grey'),'filled')
        scatter(x(body.head),y(body.head), 35, Color('yellow'),'filled', 'Marker','^')

    end
end
formatFig(fig,blkbnd);
set(gca, 'xcolor', 'none', 'ycolor', 'none')

save_figure(fig, [figDir, 'circling example 1'],'-png'); % fix name of saved figure - possible override of wg ext figs
% 
% xlim(xlims)
% ylim(ylims)
% 



% screening = close_dist;
frames = find(T.circling_1sec==1);
if isempty(frames)
    disp('no bouts of circling found')
    return
else
fig = getfig('',1,[1075 871]);
hold on
% Plot all male body positions
kolor = Color('grey');
x = mX(1:skip:end,[body.head,body.center]);
y = mY(1:skip:end,[body.head,body.center]);
plot(x',y','color',kolor)
scatter(x(:,1),y(:,1),15,kolor,"filled","^")
% axis equal square
formatFig(fig,blkbnd);
set(gca,'XColor','none','YColor','none')

% Plot the all male body positions under chase instances
kolor = Color('green');
x = mX(T.circling_1sec,[body.head,body.center]);
y = mY(T.circling_1sec,[body.head,body.center]);
plot(x',y','color',kolor)
scatter(x(:,1),y(:,1),15,kolor,"filled","^")

% % Screening
% kolor = Color('green');
% x = mX(screening,[body.head,body.center]);
% y = mY(screening,[body.head,body.center]);
% plot(x',y','color',kolor)
% scatter(x(:,1),y(:,1),15,kolor,"filled","^")

% Plot female body
x = fX(1:skip:end,[body.head,body.center,body.abdomen]);
y = fY(1:skip:end,[body.head,body.center,body.abdomen]);
plot(x',y','color',foreColor, 'LineWidth', 2)
% xlim([-2000,1500]); ylim([-2000,2000])

% Format figure
axis  equal square
h_line(0,'gray',':',2)
v_line(0,'grey',':',2)
xlim(zoom)
ylim(zoom)
formatFig(fig,blkbnd);
set(gca,'XColor','none','YColor','none')

formatFig(fig,blkbnd);
% rectangle('Position',[zoom(1),zoom(1) sum(abs(zoom)) sum(abs(zoom))],'edgecolor',foreColor,'linewidth', 1)
% save_figure(fig, 'G:\My Drive\Jeanne Lab\Presentations\Data Presentation 12.6.2024\All Chase Positions',fig_type);

% Save figure
save_figure(fig,[figDir 'Circling positions M fly'],fig_type);
end

%% FIGURE: CI & behavior components
clearvars('-except',initial_var{:})
% TODO: add descriptions to the raster plot so we know what each one
% corresponds to! 3/31
tickH = 1; % tick height
LS = 0.5; % vertical space between ticks
LW = 1; % tick line width
CI = any([T.court_chase,T.wing_ext,T.circling_1sec],2); % courtship index
T.CI = CI;

idx = [T.circling_all, T.circling_1sec,T.wing_ext_all,...
       T.wing_ext,T.chase_all,T.court_chase,CI];
kolor = {'grey', 'red'};
kolor = repmat(kolor, 1,3);
% kolor = {'grey','gold','grey','gold','grey','gold'}; % chase, wing ext, circling
y1 = 1;

r = 4;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2:4;

fig = getfig('',1);
% Temperature plot
subplot(r,c,sb(1).idx)
    plot(time, T.temperature, 'color', foreColor, 'linewidth', LW)
% Raster plot of the different features that contribute to 
subplot(r,c,sb(2).idx)
hold on
for i = 1:7
    if i<7
        C = Color(kolor{i});
    else
        C = foreColor;
    end
    x = time(idx(:,i));
    y = ones(size(x));
    X = [x,x];
    Y = [y1.*y,(y1+tickH).*y];
    plot(X',Y','color',C,'linewidth',LW) % all instances--not time restricted
    % increase the line height for the 
    if rem((i+2),2)==0 
       y1 = y1+tickH+LS;
    end
end

formatFig(fig, blkbnd,[r,c],sb);
subplot(r,c,sb(1).idx)
set(gca, 'xcolor', 'none','ycolor','none')
ylabel('\circC','color', foreColor)
subplot(r,c,sb(2).idx)
xlabel('time (min)','color', foreColor)
set(gca, 'ycolor', 'none')
ylabel('Courtship Metrics','color', foreColor)
ylim([0,y1+1+tickH])
save_figure(fig,[figDir 'Courtship Index timecourse'],fig_type);

%% FIGURE: Fly turning over time 
clearvars('-except',initial_var{:})

r = 4;
c = 1;
sb(1).idx = 1; %  temperature
sb(2).idx = 2:4; % circing

fig = getfig('',1);
subplot(r, c, sb(1).idx)
    hold on
    plot(time,T.temperature,'color', backColor,'LineWidth', 1)
    ylabel('Temp (\circC)')
subplot(r, c, sb(2).idx)
    hold on
    for sex = 1:2
        plot(time,smooth(data(sex).turning,fps*30,'moving'),'color', data(sex).color)
    end
    h_line(0,'grey','--')
    xlabel('Time (min)')
    ylabel('Turning (\circ/s)')

save_figure(fig,[figDir 'Fly turning over time'],fig_type);

%% FIGURE: Sleep
% update the figure to plot both male and female sleep over time
clearvars('-except',initial_var{:})

frames = find(m.sleep==1 & f.sleep==1);

if isempty(frames)
    disp('no bouts of sleep found')
    return
else
    bout = 5*60*parameters.FPS;
    dummy = [];
    
    % Extract sleep bouts from position data
    for sex = 1:2
        switch sex
            case 1
                x = m.pos(:,2,1); % male center
            case 2
                x = f.pos(:,2,1); % female center
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
    
    % Figure
    lw = 2;
    
    fig = getfig('',1);
    hold on
        sSpan = fps; % single second smoothing
        y1 = smooth(m.sleep,sSpan,'moving'); % male
        y2 = smooth(f.sleep,sSpan,'moving'); % female
        % plot(time, y1, 'color', Color('dodgerblue'),'LineWidth', lw)
        plot(time, y2, 'color', Color('deeppink'),'LineWidth', lw)
        plot(time, y1, 'color', Color('dodgerblue'),'LineWidth', lw)
        xlabel('time (min)')
    formatFig(fig);
    ylim([-0.2,1.2])
    set(gca, 'ytick', [0,1],'YTickLabel', {'Awake', 'Sleep'})
    
    save_figure(fig,[figDir 'Fly sleep over time'],fig_type);
    
    g = find(m.sleep);
    h = find(f.sleep);
    if isempty(g) && isempty(h)
        disp('no male or female sleep found')
    else
        return
    end
end

%% FIGURE: (WORKING 1/15) simple summary of fly positions across the trial...
% 
% % ba
% 
% % plot the occupancy of the different regions over time
% fig = getfig('', 1);
% 
% hold on




%% FIGURE: Distance, speed, and speed correlation between flies over full video course
% use the center body point to determine between-fly distance
clearvars('-except',initial_var{:})

% frame = 5100;
sz = 50;
lw = 1;
% frame_skip = 5;
% windowsize = 4; % seconds
% roi = windowsize*80;
% ROI = frame-roi:frame_skip:frame;

r = 5; c = 1;
% xlimit = [50,70];
% xlimit = [time(ROI(1)),time(ROI(end))];

fig = getfig('',1,[1032 1042]);
% TEMPERATURE
subplot(r,c,1)
    plot(time,T.temperature,'color', foreColor,'LineWidth', lw)
    % xlabel('time (s)')
    ylabel('Temp (\circC)')

% INTERFLY DISTANCE
subplot(r,c,2)
hold on
    plot(time,T.IFD,'color', foreColor,'LineWidth', lw)
    y = rangeLine(fig, 20, false);
    y1 = double(T.courtposition);
    zeroloc = T.courtposition == 0;
    y1(zeroloc) = nan;
    plot(T.time, y*y1, 'color', Color('gray'), 'LineWidth', 3)
    % xlabel('time (min)')
    ylabel('inter-fly distance (mm)')

% FLY SPEED
subplot(r,c,3); hold on 
    sSpan = fps; %single second smoothing
    y1 = smooth(m.speed,sSpan,'moving'); % male
    y2 = smooth(f.speed,sSpan,'moving'); % female
    plot(time, y1, 'color', Color('dodgerblue'),'LineWidth', lw)
    plot(time, y2, 'color', Color('deeppink'),'LineWidth', lw)
    ylabel('speed (mm/s)')

    % FLY SLEEP (in speed fig)
    % Plot instances of male sleep
    y = rangeLine(fig, 20, false);
    y1 = double(m.sleep);
    mslp = m.sleep == 0;
    y1(mslp) = nan;
    plot(T.time, y*y1, 'color', Color('dodgerblue'), 'LineWidth', 2)

    % Plot instances of female sleep
    y1 = double(f.sleep);
    fslp = f.sleep == 0;
    y1(fslp) = nan;
    plot(T.time, y*y1, 'color', Color('deeppink'), 'LineWidth', 2)

% FLY SPEED CORRELATION
subplot(r,c,4); hold on 
    tSpan = 3; % time window in seconds
    pSpan = tSpan * fps; % frames in the sliding window
    y = runningCorrelation([m.speed,f.speed], pSpan);
    y = smooth(y,pSpan,'moving');
    offset = ceil(pSpan/2);
    x = time(offset:offset+length(y)-1)';
    plot(x, y, 'color', foreColor,'LineWidth', lw)
    ylabel('M & F speed corr')
    h_line(0,'r','--',1)

% FLY WING ANGLE
subplot(r,c,5); hold on 
    plot(time, data(1).wingspread,'color', Color('dodgerblue'),'linewidth', lw) % male wingspread
    plot(time, data(2).wingspread,'color', Color('deeppink'),'linewidth', lw) % female wingspread
    ylabel('wing angle (\circ)')

% % FLY WING ANGLE
% subplot(r,c,6); hold on 
%     plot(time, data.btwnflyangle,'color', foreColor,'linewidth', 1) % left wing
%     xlabel('time (s)') 
%     ylabel('angle between flies (\circ)')

oglim = xlim;
formatFig(fig, blkbnd,[r,c]);

%xlimit = [0 64];
for i = 1:4
    subplot(r,c,i)
    set(gca, 'xcolor', 'none')
    %xlim(xlimit)
end
subplot(r,c,5)
%xlim(xlimit)

save_figure(fig,[figDir, 'Full timecourse zoom in'], fig_type);

%% FIGURE: Plot fly positions for given frames
clearvars('-except',initial_var{:})
time_roi = [58.5,59.5]; % time range (not frames) to seek
frames = find(time>=time_roi(1) & time<=time_roi(2));

% plot fly positions: 
fig = getfig('',1,[560 420]); hold on

for frame = frames
        plotFlySkeleton(fig, data(1).rawX(frame, :),data(1).rawY(frame, :),Color('dodgerblue'),true);

        plotFlySkeleton(fig, data(2).rawX(frame, :),data(2).rawY(frame, :),Color('deeppink'),true);
    axis square equal
set(gca, 'xcolor', backColor,'ycolor',backColor)
% title([num2str(data(1).mfbodyangle(frame))])
end

%     close all
%% FIGURE: Distance, speed, and speed correlation between flies over full video course
% use the center body point to determine between-fly distance
clearvars('-except',initial_var{:})

sSpan = 90; % 3 seconds
% [foreColor,backColor] = formattingColors(false); %get background colors

r = 3; c = 1;
xlimit = [0,120];

fig = getfig('',1,[1032 1042]);
% INTERFLY DISTANCE
subplot(r,c,1) 
    plot(T.time,T.IFD,'color', foreColor,'LineWidth', 1.5)
    % xlabel('time (s)')
    ylabel('inter-fly distance (mm)')
    % xlim(xlimit)
% FLY SPEED
subplot(r,c,2); hold on 
    sSpan = fps; %single second smoothing
    x = time;
    y1 = smooth(m.speed,sSpan,'moving');
    y2 = smooth(f.speed,sSpan,'moving');
    plot(x, y1, 'color', Color('dodgerblue'),'LineWidth', 2)
    plot(x, y2, 'color', Color('deeppink'),'LineWidth', 2)
    ylabel('fly speed (mm/s)')
    xlim(xlimit)

% FLY SPEED CORRELATION
subplot(r,c,3); hold on 
    tSpan = 3; % time window in seconds
    pSpan = tSpan * fps; % frames in the sliding window
    y = runningCorrelation([m.speed, f.speed], pSpan);
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

save_figure(fig,[figDir 'speed and interfly distance timecourse'], fig_type);

%% FIGURE: Plot fly positions for a given frame
frame = 1;

% plot fly positions: 
fig = getfig('',1,[560 420]);
    hold on
    plotFlySkeleton(fig, m.pos(frame, :,1),m.pos(frame, :,2),Color('dodgerblue'),true);
    plotFlySkeleton(fig, f.pos(frame, :,1),f.pos(frame, :,2),Color('deeppink'),true);
    axis square equal
set(gca, 'xcolor', backColor,'ycolor',backColor)


%% FIGURE: Distance, speed, and speed correlation between flies over full video course
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





%% FIGURE: Angle between the flies
clearvars('-except',initial_var{:})

% 1) Create a vector defined by the head and center of the male fly
P1 =[m(:,1,1), m(:,1,2)]; % male head x-y vector
P2 =[m(:,2,1), m(:,2,2)]; % male center x-y vector
v1 = P2 - P1;  % Vector for male fly

% 2) Create a vector defined by the head and center of the male fly
P3 =[f(:,1,1), f(:,1,2)]; % female head x-y vector
P4 =[f(:,2,1), f(:,2,2)]; % female center x-y vector
v2 = P3 - P4;  % Vector for female fly

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








%% FIGURE: Plot a frame with tracked points overlaid AND maybe a zoom in of the the behavior at that point?
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

%% FIGURE: Visualize body angles over time
clearvars('-except',initial_var{:})

% Compare male and female wing angles
sSpan = 1;%0*fps;
fig = getfig('',true, [1032 300]); 
hold on
  scatter(time, smooth(data(1).mfbodyangle,sSpan,'moving'),1, foreColor)
    % plot(time, smooth(angleDegrees,sSpan,'moving'), 'color', Color('deeppink'),'linewidth', 1)
xlabel('time (s)')
ylabel('body angle (\circ)')
formatFig(fig, blkbnd);
save_figure(fig,[figDir 'Body angle between M and F'],fig_type);

% figure;
% polarhistogram(angleDegrees)
%%
