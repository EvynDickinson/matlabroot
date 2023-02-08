

% QuadStep4_Videos

% Start by using QuadStep4_1_TempRate_TuningCurve.m
% Manually load or autoload data then run first step of processing with
% first 'Analysis' module in temprate tuningcurve code

%% ANALYSIS FOR VIDEO: position binned by temp and hysteresis

clearvars('-except',initial_vars{:})

types = {'hold', 'up', 'down'};

for i = 1:num.exp

  % get temp color distribution (1/4 degree incrememts)
  tp = getTempTurnPoints(data(i).T.TempProtocol{1});
  tempBins = tp.threshLow:0.25:tp.threshHigh;
  ntemps = length(tempBins)-1;

  all_data = struct;  
  for trial = 1:num.trial(i)
    % --pull full experiment data--
    % positions for all 'flies' over time
    x = data(i).data(trial).data.x_loc; 
    y = data(i).data(trial).data.y_loc;
    %distance to food
    food_dist = data(i).G(trial).TR.data(:,3); 
    % temperature for each frame
    temperature = data(i).data(trial).occupancy.temp;
    temp_rate = data(i).G(trial).TR.data(:,4); %temp rate info

    % get arena information
    well_loc = data(i).T.foodLoc(trial);
    C = data(i).data(trial).data.centre;
    r = data(i).data(trial).data.r;
    well_C = data(i).data(trial).data.wellcenters(:,well_loc);
    arena_C = data(i).data(trial).data.wellcenters(:,5);   
    wells = data(i).data(trial).data.wellcenters(:,1:4);

    % Make food well the origin
    x_offset = wells(1,well_loc);
    y_offset = wells(2,well_loc);
    wells_x = wells(1,:)-x_offset;
    wells_y = wells(2,:)-y_offset;
    X = x-x_offset;
    Y = y-y_offset;

    % pull the data points for the 'correct' conditions
    for tt = 1:3 %hold, up, down, etc.
        for ii = 1:ntemps-1
            all_data(tt,ii).G(trial).pos = [nan,nan];
            
            loc_1 = temperature>=tempBins(ii) & temperature<tempBins(ii+1);
            switch types{tt}
                case 'hold'
                    loc_2 = temp_rate==0;
                case 'up'
                    loc_2 = temp_rate>0;
                case 'down'
                    loc_2 = temp_rate<0;
            end
            loc_2 = [loc_2; false];
            loc = loc_1 & loc_2;
            % skip if  there are no frames that match
            if ~any(loc)
               continue
            end

            % pull all points that fit this temp bin and temp rate (for
            % this trial)
            dist = food_dist(loc);
            newX = X(loc,:);
            newY = Y(loc,:);
            newX = reshape(newX,[numel(newX),1]); %switch matrix to vector
            newY = reshape(newY,[numel(newY),1]); 
            [plotData, WELLS] = deal([]);

            % Rotate to correct orientation
            switch well_loc
                case 1
                    plotData(:,1) = newY;
                    plotData(:,2) = -newX;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 2 
                    plotData(:,1) = newX;
                    plotData(:,2) = -newY;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 3
                    plotData(:,1) = -newY;
                    plotData(:,2) = newX;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
                case 4 
                    plotData(:,1) = newX;
                    plotData(:,2) = newY;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
            end
            all_data(tt,ii).G(trial).pos = plotData;
            all_data(tt,ii).G(trial).wells = WELLS;
            all_data(tt,ii).G(trial).dist = dist;
        end
    end
  end
  
  % Plot avg fly positions...
  X_Edges = [-400 -354 -308 -262 -216 -170 -124 -78 -32 14 60 106 152 198 244 290 336 382 428];
  Y_Edges = [-680 -633 -586 -539 -492 -445 -398 -351 -304 -257 -210 -163 -116 -69 -22 25 72 119 166];

  % Group all fly positions across trials
  for tt = 1:3
    for ii = 1:ntemps-1
        [wells, temp, dist] = deal([]);
        for trial = 1:num.trial(i)
            x = all_data(tt,ii).G(trial).pos;
            if all(isnan(x))
                dist(trial) = nan;
                continue
            end
            temp = autoCat(temp,x);
            well_pos = all_data(tt,ii).G(trial).wells;
            wells = autoCat(wells,well_pos);
            dist(trial) = mean(all_data(tt,ii).G(trial).dist,'omitnan');
        end
        all_data(tt,ii).pos = temp;
        all_data(tt,ii).wells = wells;
        all_data(tt,ii).dist = mean(dist,'omitnan');
        all_data(tt,ii).dist_err = std(dist,0,'omitnan');
        all_data(tt,ii).temp = tempBins(ii);
        if ~isempty(temp)
            all_data(tt,ii).N = histcounts2(all_data(tt,ii).pos(:,1),all_data(tt,ii).pos(:,2),X_Edges,Y_Edges);
        end
    end
  end

  grouped(i).video.all_data = all_data;
  grouped(i).video.temperatures = tempBins;
end  

%% MAKE VIDEO OF POSITIONS OVER TEMP CHANGES
clearvars('-except',initial_vars{:})

exp = 1;

all_data = grouped(exp).video.all_data;
ntemps = length(grouped(exp).video.temperatures)-1;

% Find max count across all frames & pull distance data & organize images
% into movie order
idx = 1;
[dist, n, imageStack] = deal([]);
% decreasing temp data first
tt = 3; 
for ii = ntemps-1:-1:2
    tot = sum(all_data(tt,ii).N,'all');
    % Pull the distance data into a single structure
    dist  = autoCat(dist,[all_data(tt,ii).dist,all_data(tt,ii).dist_err,all_data(tt,ii).temp,-1]);
    % collate image data
    img = all_data(tt,ii).N;
    imageStack(:,:,idx) = rot90((img./tot)*100);
    n(idx) = max(imageStack(:,:,idx),[],'all');
    idx = idx+1;
end
% increasing temp data next
tt = 2; 
for ii = 2:ntemps-1
    tot = sum(all_data(tt,ii).N,'all');
    % Pull the distance data into a single structure
    dist  = autoCat(dist,[all_data(tt,ii).dist,all_data(tt,ii).dist_err,all_data(tt,ii).temp,1]);
    % collate image data
    img = all_data(tt,ii).N;
    imageStack(:,:,idx) = rot90((img./tot)*100);
    n(idx) = max(imageStack(:,:,idx),[],'all');
    idx = idx+1;
end; clear idx tt ii

% % Find max count across all frames & pull distance data & organize images
% % into movie order
% idx = 1;
% [dist, n, imageStack] = deal([]);
% % decreasing temp data first
% tt = 3; 
% for ii = ntemps-1:-1:2
%     n(idx) = max(all_data(tt,ii).N,[],'all');
%     % Pull the distance data into a single structure
%     dist  = autoCat(dist,[all_data(tt,ii).dist,all_data(tt,ii).dist_err,all_data(tt,ii).temp,-1]);
%     % collate image data
%     imageStack(:,:,idx) = rot90(all_data(tt,ii).N);
%     idx = idx+1;
% end
% % increasing temp data next
% tt = 2; 
% for ii = 2:ntemps-1
%     n(idx) = max(all_data(tt,ii).N,[],'all');
%     % Pull the distance data into a single structure
%     dist  = autoCat(dist,[all_data(tt,ii).dist,all_data(tt,ii).dist_err,all_data(tt,ii).temp,1]);
%     % collate image data
%     imageStack(:,:,idx) = rot90(all_data(tt,ii).N);
%     idx = idx+1;
% end; clear idx tt ii


% Determine plot limits
colorLim = [0,ceil(max(n))];
% colorLim = [0,5];
nFrames = size(imageStack,3);
tp = getTempTurnPoints(data(exp).T.TempProtocol{1});
tempLim = [tp.threshLow,tp.threshHigh];
distLim = [floor(min(dist(:,1)-dist(:,2))), ceil(max(dist(:,1)+dist(:,2)))];
rateShift = find(diff(dist(:,4))>1);

timeStep = median(abs(diff(dist(:,3))))/abs(tp.rates(1)); %temp step/exp temp rate
time = [0;cumsum(abs(diff(dist(:,3)))./abs(tp.rates(1)))];


% Video formatting

video_name = [saveDir expNames{exp} ' summary video.avi'];   % name of the video
video_fps = 5;        

% Plot formatting
LW = 1.5;
SZ = 50;
r = 4;
c = 3;
sb(1).idx = 1:c;   % current temperature
sb(2).idx = [(4:c:r*c),(5:c:r*c)]; % image stack
sb(3).idx = (6:c:r*c); % distance temp relationship

% Open video object to which you will write new frames
v = VideoWriter(video_name, 'Uncompressed AVI');
v.FrameRate = video_fps;
open(v);

fig = figure; set(fig,'color','k','pos',[1939 638 1003 641])

for i = 1:nFrames

% current temperature
subplot(r,c,sb(1).idx)
    cla
    set(gca,'box','off','color','k','xcolor','w','ycolor','w','Linewidth',LW); hold on
    ROI = 1:i;
%     x = (ROI-1)*2;
    x = time(ROI);
    y = dist(ROI,3); % temperature data
    plot(x,y,'linewidth',LW,'color',Color('dodgerblue'))
    if i > rateShift
        ROI = rateShift+1:i;
        plot(x(ROI),y(ROI),'linewidth',LW,'color',Color('red'))
    end
    scatter(x(end),y(end),SZ,'w','filled')
    xlim([-2,nFrames*2+2])
    ylim(tempLim)
    xlabel('time (min)','color','w','fontsize',12)
    ylabel('\circC','color','w','fontsize',12)

% image stack
subplot(r,c,sb(2).idx)
    img = imageStack(:,:,i);
    imagesc(img)
    axis square
    set(gca, 'color','k','xtick',[],'ytick',[], 'box','off')
    clim(colorLim)
    cb = colorbar;
    set(cb,'Location', 'westoutside','color','w');
    set(cb,'Ticks', colorLim) 
    ylabel(cb,'Occupation (%)','FontSize',12)
    hold on
    scatter(9.5,4.5,20,'r','filled')
    hold off
    
% distance temp relationship
subplot(r,c,sb(3).idx)
    [h_ROI, c_ROI] = deal([]);
    cla
    hold on
    if i>1
        if i<=rateShift && i>1 %cooling data
            c_ROI = 1:i;
        elseif i>rateShift
            c_ROI = 1:rateShift;
        end
        % plot cooling data
        x = dist(c_ROI,3);
        y = dist(c_ROI,1);
        y_err = dist(c_ROI,2);
        fill_data = error_fill(x, y, y_err);
        h = fill(fill_data.X, fill_data.Y, Color('dodgerblue'), 'EdgeColor','none');
        set(h, 'facealpha', 0.3)
        plot(x,y,'linewidth',LW,'color',Color('dodgerblue'))
    
        if i>rateShift
            h_ROI = rateShift+1:i;
            % plot current heating data
            x = dist(h_ROI,3);
            y = dist(h_ROI,1);
            y_err = dist(h_ROI,2);
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, Color('red'), 'EdgeColor','none');
            set(h, 'facealpha', 0.3)
            plot(x,y,'linewidth',LW,'color',Color('red'))
        end
    end
    scatter(dist(i,3),dist(i,1),SZ,'w','filled')
   
    set(gca,'box','off','color','k','xcolor','w','ycolor','w','ydir','reverse','Linewidth',LW)
    xlim(tempLim)
    ylim(distLim)
    xlabel('temp (\circC)','color','w','fontsize',12)
    ylabel('proximity to food (mm)','color','w','fontsize',12)

    % save the current figure as new frame in the video
    f = getframe(fig);
    writeVideo(v, f)

end


close(v) %close the movie writer object
close(fig)

%%
%% FIGURE: surface map of distance-temp tuning curve 
clearvars('-except',initial_vars{:})

blkbgd = true;

if blkbgd
    foreColor = 'w';
    backColor = 'k';
    gridalpha = 0.3;
else %white background
    foreColor = 'k';
    backColor = 'w';
    gridalpha = 0.5;
%     gridColor = [0.15 0.15 0.15];
end

x_limits = [0,num.exp+1];
y_limits = [14, 23]; % TODO: update this auto
z_limits = [14,28]; % TODO: update this auto
colorLimits = [15.7617 26.2390];

% get strain names
strain_label = [];
for exp = 1:num.exp
    i = expOrder(exp);
    temp = strsplit(expNames{i},' ');
    strain_label{exp} = temp{1};
end

buff = 0.3;

fig = figure; set(fig,'color',backColor,"Position",[1934 222 1061 836]);%[1934 468 1061 590])

for exp = 1:num.exp

    i = expOrder(exp);

%     temp = strsplit(expNames{i},' ');
%     strain_label{exp} = temp{1};


    z = grouped(i).dist.distavgbytemp(:,2); % distance
    loc = isnan(z); %no-data areas (temp bins outside perview of this exp)
    y = grouped(i).dist.distavgbytemp(:,1); % temperature
    % remove nans
    z(loc) = [];
    y(loc) = [];
    
    % Adjust temp and distance vectors for a surface plot
    Y = repmat(y,[1,2]);
    Z = repmat(z,[1,2]);
    
    % Build x-location map for this experiment group
    x = [exp-buff,exp+buff];
    X = repmat(x,[length(y),1]);

    surf(X,Y,Z,'edgecolor','none','facecolor','interp');
    
    % save the surf shape into a new folder...
    plotData = [];
    plotData.X = X; plotData.Y = Y; plotData.Z = Z;
    save([saveDir expNames{i} ' distance tuning curve surf map data.mat'],'plotData','-mat');

    % next experimental group
    hold on
end
% FORMATTING
colormap(flipud(parula))
set(gca,'zdir','reverse')        
% xlabel('Strain')
ylabel('temperature (\circC)')
zlabel('proximity to food (mm)')
set(gca,'color',backColor,'xcolor',foreColor,'ycolor',foreColor,'zcolor',foreColor)
set(gca,'gridAlpha',gridalpha,'CLim',colorLimits)
set(gca,'linewidth',2,'fontsize',12,'fontname','arial')

xlim(x_limits)
ylim(y_limits)
zlim(z_limits)
set(gca,'XTick',1:num.exp,'xticklabel',strain_label)

% ROTATE THE OBJECT...
% for i = 1:36
%    camorbit(10,0,'data',[0 1 0])
%    drawnow
%    pause(0.1)
% end

% % save figure
% save_figure(fig,[saveDir expGroup ' distance tuning curve surf map_2'],'-png',true,false);


%%
% data visualization techniques

% Single angle view
CameraPosition = [63.8311 18.4141 9.1428];
CameraPositionMode = 'manual';
CameraTarget = [3.5000 18.5000 21];
CameraTargetMode = 'manual';
CameraUpVector = [0 0 -1];
CameraUpVectorMode = 'manual';
CameraViewAngle = 7.2224;
CameraViewAngleMode = 'manual';
DataAspectRatio = [1 1.2857 2];
DataAspectRatioMode = 'manual';

% Turn on manual camera views: 
mode_set = 'manual';
set(gca,...
     'CameraPositionMode',mode_set,...
     'CameraTargetMode', mode_set,...
     'CameraUpVectorMode', mode_set,...
     'CameraViewAngleMode', mode_set,...
     'DataAspectRatioMode', mode_set)
set(gca,...
    'CameraPosition',CameraPosition,...
     'CameraTarget', CameraTarget,...
     'CameraUpVector', CameraUpVector,...
     'CameraViewAngle', CameraViewAngle,...
     'DataAspectRatio', DataAspectRatio)


% set(gca,'xdir','reverse')


%% Automate sliding between positions:
%baseline

axis vis3d
% Set specific views
set(gca,'CameraPosition',[3.5000 18.5000 -100.2440],'CameraUpVector', [0  1  0])

% % % Get current view positions:
% ax = gca;
% fprintf(['\nCameraPosition = [' num2str(ax.CameraPosition) ']; \n'])
% fprintf(['CameraUpVector = [' num2str(ax.CameraUpVector) ']; \n'])


% for jj = 1:nViews
%     set(gca,'CameraPosition',positions(jj,:),'CameraUpVector',camVectors(jj,:))
%     set(gca,'ytick',y_limits(1):y_limits(2))
%     pause(0.3)
% end



% VIEW 1: Temperature and strain ONLY
positions = [3.5000      18.5000     -100.2440;...% VIEW 1:
%              -15.3356     -22.0952     -75.4012;...
             -29.915     -27.1312     -51.0786]; 

cam_vectors = [0  1  0;...
%                0  0 -1;...
               0  0 -1]; 

nViews = size(positions,1);
intervals = 20; % frames between viewpoints
time = 1:intervals:intervals*nViews;
qTime = 1:time(end);

[camPos, camVectors] = deal([]);
for jj = 1:3
    camPos(:,jj) = interp1(time,positions(:,jj),qTime,'linear');
    camVectors(:,jj) = interp1(time,cam_vectors(:,jj),qTime,'linear');
end


for jj = 1:length(qTime)
    set(gca,'CameraPosition',camPos(jj,:),'CameraUpVector',camVectors(jj,:))
    set(gca,'ytick',y_limits(1):y_limits(2))
    pause(0.05)
end




% 

%%

% VIEW 1: Temperature and strain ONLY
positions = [3.5000      18.5000     -100.2440;...% VIEW 1:
             -15.3356     -22.0952     -75.4012;...
             -29.915     -27.1312     -51.0786;...
             -56.3038      10.4542       5.5926]; 

cam_vectors = [0  1  0;...
               0  0 -1;...
               0  0 -1;...
               0  0 -1]; 

nViews = size(positions,1);
intervals = 20; % frames between viewpoints
time = 1:intervals:intervals*nViews;
qTime = 1:time(end);

[camPos, camVectors] = deal([]);
for jj = 1:3
    camPos(:,jj) = interp1(time,positions(:,jj),qTime,'spline');
    camVectors(:,jj) = interp1(time,cam_vectors(:,jj),qTime,'spline');
end


for jj = 1:length(qTime)
    set(gca,'CameraPosition',camPos(jj,:),'CameraUpVector',camVectors(jj,:))
    set(gca,'ytick',y_limits(1):y_limits(2))
    pause(0.05)
end








