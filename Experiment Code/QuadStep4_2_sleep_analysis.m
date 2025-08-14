
% Find 'sleeping' points in data

%% (DONT RUN) visualize grid spacing & single trial visualization
% nbins = 50;
% i = 1;
% trial = 1;
% vid = 1;
% frame = 1;
% 
% % pull info for the first trial:
% dataDate = data(i).T.Date{trial};
% vid_name = data(i).T.ExperimentID{trial};
% vidDir = [baseFolder dataDate '/' vid_name '_'];
% videoPath = [vidDir num2str(vid) '.avi'];
% movieInfo = VideoReader(videoPath); %read in video

% % Set axis limits for the selected arena
% x = data(i).data(trial).data.centre(1);
% y = data(i).data(trial).data.centre(2);
% r = data(i).data(trial).data.r;
% xlimit = [x-(r+50),x+(r+50)];
% ylimit = [y-(r+50),y+50+r];
% 
% % Plot image of video frame
% fig = figure; set(fig,'pos',[-1030 279 772 1009],'color','k');
% currentImg = rgb2gray(read(movieInfo,frame));
% imshow(currentImg)
% xlim(xlimit); ylim(ylimit);   
% 
% % find the 'auto bin' lines
% xedge = linspace(xlimit(1),xlimit(2),nbins+1);
% yedge = linspace(ylimit(1),ylimit(2),nbins+1);
% 
% % Plot the bin outline edges:
% h_line(yedge,'yellow','-',0.25) 
% v_line(xedge,'yellow','-',0.25)

% %% Single trial quant and figure

% i = 1;
% 
% sleepData = struct;
% for trial = 1:num.trial(i)
% 
%     % Set axis limits for the selected arena
%     x = data(i).data(trial).data.centre(1);
%     y = data(i).data(trial).data.centre(2);
%     r = data(i).data(trial).data.r;
%     xlimit = [x-(r+50),x+(r+50)];
%     ylimit = [y-(r+50),y+50+r];
% 
%     % find the 'auto bin' lines
%     xedge = linspace(xlimit(1),xlimit(2),nbins+1);
%     yedge = linspace(ylimit(1),ylimit(2),nbins+1);
%     
%     % pull the fly locations during the trial
%     x_loc = data(i).data(trial).data.x_loc;
%     y_loc = data(i).data(trial).data.y_loc;
%     trial_length = size(y_loc,1);
%     
%     N = [];
%     for frame = 1:trial_length
%         X = x_loc(frame,:); X(isnan(X)) = [];
%         Y = y_loc(frame,:); Y(isnan(Y)) = [];
%         N(:,:,frame) = histcounts2(X,Y,xedge,yedge);
%     end
%     
%     % find grid space that have continuous occupation for more than min_duration 
%     frameCount = N(:,:,1);
%     for frame = 2:trial_length
%         currFrame = N(:,:,frame); % current frame locations
%         resetLoc = currFrame==0; % locations that do not have flies and thus need a count reset
%         
%         tempCount = frameCount(:,:,frame-1)+currFrame; % add current frames to list
%         tempCount(resetLoc) = 0; % reset counts for spots with no flies
%         
%         frameCount(:,:,frame) = tempCount; % add current count into the saving structure
%     end
%     
%     % Find instances of 'sleep'?
%     [sleepingCount,sleepLoc] = deal([]);
%     for frame = 1:trial_length
%         sleepLoc(:,:,frame) = frameCount(:,:,frame) > min_duration;
%         sleepingCount(frame) = sum(sum(sleepLoc(:,:,frame)));
%     end
% 
%     % Save into group structures
%     sleepData(trial).frameCount = frameCount;
%     sleepData(trial).sleepLoc = sleepLoc;
%     sleepData(trial).sleepingCount = sleepingCount;
% %     fprintf(['\nDone exp ' num2str(i) ' trial ' num2str(trial)])
% end

% sb(1).idx = 1;
% sb(2).idx = 2:3;
% r = 3;
% c = 1;
% 
% fig = figure;
% subplot(r,c,sb(1).idx)
%     time = data(i).data(trial).occupancy.time;
%     temp = data(i).data(trial).occupancy.temp;
%     plot(time, temp, 'color','w', 'linewidth', 1)
%     xlabel('time (min)')
%     ylabel('temp (\circC)')
% subplot(r,c,sb(2).idx)
%     plot(time, sleepingCount, 'color','w', 'linewidth', 1)
%     ylabel('# flies sleeping')
%      xlabel('time (min)')
% 
% formatFig(fig,true,[r,c],sb);
%% (DONT RUN) SLEEP METHODS

% % ---- Vectorize the data (find the flies that are sleeping....) -----
% 
% % find food well location for distance capture later
% foodWellLoc = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));
% c1 = foodWellLoc(1);
% c2 = foodWellLoc(2);
% 
% % create empty matrixes for the x and y positions of the sleeping flies
% sleeping = struct;
% [sleeping.X, sleeping.Y, sleeping.all_distance] = deal(nan(trial_length,data(i).T.NumFlies(trial)));
% sleeping.sleepNum = zeros(trial_length,1);
% [sleeping.dist_avg, sleeping.dist_err] = deal(nan(trial_length,1));
% % assign data by frame
% for frame = 1:trial_length
%     frame_data = sleepLoc(:,:,frame);
%     binLoc = find(frame_data>0);
%     
%     % Find the coordinates of the sleeping flies bins from the discretized data
%     y_row = ceil(binLoc/nbins);
%     x_row = rem(binLoc-1,nbins)+1;
%     x_position = (xedge(x_row) + xedge(x_row+1))/2;
%     y_position = (yedge(y_row) + yedge(y_row+1))/2;
%     
%     % add position data to the matrix:
%     if ~isempty(binLoc)
%         % number of flies sleeping
%         sleepNum = length(x_position);
%         sleeping.sleepNum(frame) = sleepNum;
%         % location of sleeping flies
%         sleeping.X(frame,1:sleepNum) = x_position;
%         sleeping.Y(frame,1:sleepNum) = y_position;
%         % distance to food...
%         temp_dist = sqrt((x_position-c1).^2 + (y_position-c2).^2)./pix2mm;
%         sleeping.all_distance(frame,1:sleepNum) = temp_dist;
%         % average distance:
%         sleeping.dist_avg(frame) = mean(temp_dist);
%         sleeping.dist_err(frame) = std(temp_dist);
%     end
% end
% 
% % Vectorize the data (find the flies that are sleeping....)
% % find x & y of flies that are sleeping...
% for frame = 1:trial_length
%     frame_data = sleepLoc(:,:,frame);
%     sleeping = find(frame_data>0);
% 
%     % Find the coordinates of the sleeping flies bins from the discretized data
%     y_row = ceil(sleeping/nbins);
%     x_row = rem(sleeping-1,nbins)+1;
%     x_position = (xedge(x_row) + xedge(x_row+1))/2;
%     y_position = (yedge(y_row) + yedge(y_row+1))/2;
% 
% end
% 
% 
% % Visualize the 'sleeping flies'
% vid = data(i).data(trial).data.T.vidNums(frame);
% vidFrame = data(i).data(trial).data.T.vidFrame(frame);
% 
% % pull info for the first trial:
% dataDate = data(i).T.Date{trial};
% vid_name = data(i).T.ExperimentID{trial};
% vidDir = [baseFolder dataDate '/' vid_name '_'];
% videoPath = [vidDir num2str(vid) '.avi'];
% movieInfo = VideoReader(videoPath); %read in video
% 
% % Set axis limits for the selected arena
% x = data(i).data(trial).data.centre(1);
% y = data(i).data(trial).data.centre(2);
% r = data(i).data(trial).data.r;
% xlimit = [x-(r+50),x+(r+50)];
% ylimit = [y-(r+50),y+50+r];
% 
% % Plot image of video frame
% fig = figure; set(fig,'pos',[-1030 279 772 1009],'color','k');
% currentImg = rgb2gray(read(movieInfo,vidFrame));
% imshow(currentImg)
% xlim(xlimit); ylim(ylimit);  
% 
% 
% 
% % find the 'auto bin' lines
% xedge = linspace(xlimit(1),xlimit(2),nbins+1);
% yedge = linspace(ylimit(1),ylimit(2),nbins+1);
% 
% % Plot the bin outline edges:
% h_line(yedge,'yellow','-',0.25) 
% v_line(xedge,'yellow','-',0.25)
% 
% 
% x_points = [xedge(x_row); xedge(x_row+1);  xedge(x_row+1); xedge(x_row)];
% y_points = [yedge(y_row); yedge(y_row); yedge(y_row+1); yedge(y_row+1)];
% 
% patch(x_points,y_points,Color('gold'),'FaceAlpha',.5,'EdgeColor','none');
% 
% 
%% ANALYSIS: Generate and save the sleep quanitification
clearvars('-except',initial_vars{:})
paths = getPathNames;
nbins = 50;

% Create sleep data for unprocessed files (trial by trial)
for i = 1:num.exp
    
    TP = getTempTurnPoints(data(i).temp_protocol);
    fps = TP.fps;

    % How long does a fly need to be still to count as 'sleep'
    min_duration = 5*fps*60; % 5 mins * 3fps*60sec = data point number that must be met 'unmoving'
    
    T = data(i).T;
    for trial = 1:num.trial(i)
        trial_ID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
        sleep_file = [baseFolder paths.single_trial '/' trial_ID  '/' T.ExperimentID{trial} ' sleeping data.mat'];  
        
        if ~exist(sleep_file,"file")
            sleepData = struct;
            %preallocate for speed and space
            trial_length = length(data(i).data(trial).occupancy.time);
            [N,frameCount,sleepLoc] = deal(nan(nbins,nbins,trial_length));
            sleepingCount = zeros(trial_length,1);

            % Set axis limits for the selected arena for bins that will dictate single fly sizes
            x = data(i).data(trial).data.centre(1);
            y = data(i).data(trial).data.centre(2);
            r = data(i).data(trial).data.r;
            xlimit = [x-(r+50),x+(r+50)];
            ylimit = [y-(r+50),y+50+r];

            % find the 'auto bin' lines (aka the bin edge coordinates)
            xedge = linspace(xlimit(1),xlimit(2),nbins+1); 
            yedge = linspace(ylimit(1),ylimit(2),nbins+1);

            % pull the fly locations during the trial
            x_loc = data(i).data(trial).data.x_loc;
            y_loc = data(i).data(trial).data.y_loc;
            % find the bins in which flies were present for this frame and add to the giant 'N' structure
            for frame = 1:trial_length
                X = x_loc(frame,:); % x-coordinates for all flies on this camera frame
                X(isnan(X)) = []; % removes any nan locations from the list
                Y = y_loc(frame,:); 
                Y(isnan(Y)) = [];
                N(:,:,frame) = histcounts2(X,Y,xedge,yedge); 
            end
            % (at this point, N is a large matrix where we will look for flies that
            % stay in the same x-y space bin for longer than 5 minutes (the def of
            % sleep in flies))

            % find grid space that have continuous occupation for more than min_duration
            frameCount(:,:,1) = N(:,:,1); % frameCount is the running count for # of frames that have a fly
            for frame = 2:trial_length
                currFrame = N(:,:,frame); % current frame locations
                resetLoc = currFrame==0; % locations that do not have flies and thus need a count reset

                tempCount = frameCount(:,:,frame-1)+currFrame; % add 1 to the frame count from previous frame count
                tempCount(resetLoc) = 0; % reset counts for spots with no flies

                frameCount(:,:,frame) = tempCount; % add current count into the saving structure
            end

            % ---- Vectorize the data (find the flies that are sleeping....) -----

            % pull coordinates of the food well for distance capture later
            foodWellLoc = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));
            c1 = foodWellLoc(1); % x-coordinate for the center of the food well
            c2 = foodWellLoc(2); % y-coordinate for the center of the food well
            
            % create empty matrixes for the x and y positions of the sleeping flies
            sleeping = struct;
            [sleeping.X, sleeping.Y, sleeping.all_distance] = deal(nan(trial_length,data(i).T.NumFlies(trial))); 
            sleeping.sleepNum = zeros(trial_length,1);
            [sleeping.dist_avg, sleeping.dist_err] = deal(nan(trial_length,1));

            % assign data by frame
            for frame = 1:trial_length
                frame_data = frameCount(:,:,frame) > min_duration; % puts a 'T' for any bin location that has had a stationary fly for >5 mins
                binLoc = find(frame_data>0); % vectorize the location % TODO HERE: resume sweep
                
                % Find the coordinates of the sleeping flies bins from the discretized data
                y_row = ceil(binLoc/nbins);
                x_row = rem(binLoc-1,nbins)+1;
                x_position = (xedge(x_row) + xedge(x_row+1))/2; %x-location is the middle of the x-bins
                y_position = (yedge(y_row) + yedge(y_row+1))/2; %y-location is the middle of the y-bins
                
                % add position data to the matrix:
                if ~isempty(binLoc)
                    % number of flies sleeping
                    sleepNum = length(x_position);
                    sleeping.sleepNum(frame) = sleepNum;
                    % location of sleeping flies
                    sleeping.X(frame,1:sleepNum) = x_position;
                    sleeping.Y(frame,1:sleepNum) = y_position;
                    % distance to food...
                    temp_dist = sqrt((x_position-c1).^2 + (y_position-c2).^2)./pix2mm;
                    sleeping.all_distance(frame,1:sleepNum) = temp_dist;
                    % average distance:
                    sleeping.dist_avg(frame) = mean(temp_dist);
                    sleeping.dist_err(frame) = std(temp_dist);
                end
            end

            save(sleep_file,'sleeping','-v7.3'); 
        end   
        disp([num2str(i) ' | ' num2str(trial)])
        clear N preallocate frameCount sleepingCount sleepLoc resetLoc tempCount
    end
    disp(['Done exp ' expNames{i}])
end
clearvars('-except',initial_vars{:})

%% ANALYSIS: Load previously created sleep data files and process data:
paths = getPathNames;
sleep = struct;
for i = 1:num.exp
    T = data(i).T;
    for trial = 1:num.trial(i)
        trial_ID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
        sleep_file = [baseFolder paths.single_trial '/' trial_ID  '/' T.ExperimentID{trial} ' sleeping data.mat'];  

        if exist(sleep_file,"file")
            load(sleep_file,'sleeping');
            sleep(i).trial(trial) = sleeping;
            clear sleeping
        else
            h = warndlg(['Warning: missing sleep data for ' expNames{i}]);
            uiwait(h)
        end
    end
end
disp('Loaded all sleep data')

% Process and prep the data for further analysis
for i = 1:num.exp
    [sleep(i).num, sleep(i).fract_sleep] = deal(zeros(length(grouped(i).temp),num.trial(i)));
    sleep(i).distance = nan(length(grouped(i).temp),num.trial(i));
    for trial = 1:num.trial(i)
        %number of sleeping flies
        inputdata = sleep(i).trial(trial).sleepNum;
        sleep(i).num(1:length(inputdata),trial) = inputdata;
        sleep(i).fract_sleep(1:length(inputdata),trial) = inputdata/data(i).T.NumFlies(trial);
        %distance to food for sleeping flies
        inputdata = sleep(i).trial(trial).dist_avg;
        sleep(i).distance(1:length(inputdata),trial) = inputdata;
    end
    sleep(i).sleepfract_avg = mean(sleep(i).fract_sleep,2,'omitnan');
    sleep(i).sleepfract_err = std(sleep(i).fract_sleep,0,2,'omitnan');
end

% Cluster the sleeping flies by temperature
for i = 1:num.exp  
    temps = unique(data(i).G(1).TR.temps);
    rateIdx = data(i).G(1).TR.rateIdx;
    tempIdx = data(i).G(1).TR.tempIdx;
    % find rate index
    heatRate = find(data(i).G(1).TR.rates>0);
    coolRate = find(data(i).G(1).TR.rates<0);
    try 
        holdRate = find(data(i).G(1).TR.rates==0);
        ntypes = 3;
    catch
        ntypes = 2;
    end
    
    for temp = 1:length(temps)
        for type = 1:ntypes
            switch type
                case 1 %heating
                    g_name = 'increasing';
                    idxSex = heatRate;
                case 2 %cooling
                    g_name = 'decreasing';
                    idxSex = coolRate;
                case 3 %holding
                    g_name = 'holding';
                    idxSex = holdRate;
            end
            %fraction of flies sleeping
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            sleep(i).(g_name)(temp,1) = mean(mean(sleep(i).fract_sleep(loc,:),2,'omitnan'),'omitnan'); %avg 
            sleep(i).(g_name)(temp,2) = std(mean(sleep(i).fract_sleep(loc,:),1,'omitnan'),'omitnan');%./num.trial(i); %err
            %distance of sleeping flies
            sleep(i).([g_name '_dist'])(temp,1) = mean(mean(sleep(i).distance(loc,:),2,'omitnan'),'omitnan');
            sleep(i).([g_name '_dist'])(temp,2) = std(mean(sleep(i).distance(loc,:),1,'omitnan'),'omitnan');
        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        sleep(i).temp_all(temp,1) = mean(mean(sleep(i).fract_sleep(loc,:),2,'omitnan'),'omitnan'); %avg 
        sleep(i).temp_all(temp,2) = std(mean(sleep(i).fract_sleep(loc,:),1,'omitnan'),'omitnan')./num.trial(i);% %err
        %distance
        sleep(i).tempBinDist(temp,1) = mean(mean(sleep(i).distance(loc,:),2,'omitnan'),'omitnan');
        sleep(i).tempBinDist(temp,2) = std(mean(sleep(i).distance(loc,:),1,'omitnan'),'omitnan')./num.trial(i);
    end
    sleep(i).temps = temps;
end


% Thermal threat quantifcation


% Thermal threat quantification and avg sleep quantity
for i = 1:num.exp
    tempProto = data(i).temp_protocol;
    if data(i).hold_exp
        temp_protocol  = 'linear_ramp_F_25-17';
    else
        temp_protocol = data(i).temp_protocol;
    end
    tPoints = getTempTurnPoints(temp_protocol); 
    fps = tPoints.fps; 

    if strcmp(data(i).temp_protocol,'Large_temp_sweep_15_35') ||...
       strcmp(data(i).temp_protocol,'Large_temp_sweep_15_35_FPS6') 
       sleep(i).avg_quant = nan;
       sleep(i).thermalThreat = nan;
       disp(['Skipping thermal threat for ' grouped(i).name])
       continue
    end
    
    demoRamp_idx = tPoints.down(1,1):tPoints.up(1,2); % points for a full ramp down and up
    
    % Thermal threat
    temp_ramp = grouped(i).temp(demoRamp_idx);
    thermalThreat = (sum(25-temp_ramp)/fps)/1000; % quant of time not at 25 *only punishes cold though...
    
    % Sleep duration
    
    nFlies = zeros(1,num.trial(i));
    nRamps = size(tPoints.up,1);
    for ramp = 1:nRamps
        rampROI = tPoints.down(ramp,1):tPoints.up(ramp,2);
        nFlies = nFlies + sum(sleep(i).num(rampROI,:),1);
    end
    avg_sleepDuration = ((nFlies/nRamps)./(data(i).T.NumFlies'))./fps; % avg fly sleep duration (in sec) during a single temp ramp
    % save data to sleep structure
    sleep(i).avg_quant = avg_sleepDuration;
    sleep(i).thermalThreat = thermalThreat;
end

initial_vars{end+1} = 'sleep';
initial_vars{end+1} = 'fps';
clearvars('-except',initial_vars{:})

%ANALYSIS: Sleep duration & start and stop of sleep
% Avg sleep duration
for i = 1:num.exp
    timing = struct;
    [sleepON,sleepOFF,sleepLength] = deal([]);
    for trial = 1:num.trial(i)
        dist_avg = sleep(i).trial(trial).all_distance;
        loc = diff(dist_avg)==0; % find frames where flies are stationary (where the distance to food does not change)
        duration_matrix = zeros(size(loc));
        duration_matrix = [zeros(1,size(loc,2)); duration_matrix]; %account for the first step where there is no sleep
        old_row = loc(1,:);
        for frame = 2:size(loc,1) % add frame count for each next moment of sleep
            curr_row = loc(frame,:);
            addLoc = ~(curr_row==0);
            duration_matrix(frame,addLoc) = old_row(addLoc) + curr_row(addLoc);
            old_row = duration_matrix(frame,:);
        end
        
        sleepDurationMatrix = duration_matrix./(fps*60); %sleep duration in minutes
        timing(trial).sleepduration = sleepDurationMatrix;
        % what is the max sleep length for these segments? 
        sleepDurationMatrix(sleepDurationMatrix==0)=nan;
%         fig = getfig('',1);
%         hold on
%         time = data(i).data(trial).occupancy.time;
%         for trial = 1:num.trial(i)
%             scatter(time, sleepDurationMatrix(:,trial),10)
%         end
%         xlabel('time (min)')
%         ylabel('sleep duration (min)')
%         formatFig(fig,true);
%         title(expNames{i},'color','w')
%         xlim([0,1000])

        % Find the 'end' of sleep duration and stop location (time)
        startLoc = (diff(diff(duration_matrix)))<0; % Find the start of sleep (time)
        stopLoc = (diff(duration_matrix))<0;
        total_duration = duration_matrix(stopLoc)/(fps*60); %duration in minutes

        [sleepStoppedFrame, sleepStartFrame] = deal([]); %when in time (frame) did the end of sleep occur?
        for fly = 1:size(stopLoc,2)
              sleepStoppedFrame = autoCat(sleepStoppedFrame,find(stopLoc(:,fly)),false);
              sleepStartFrame = autoCat(sleepStartFrame,find(startLoc(:,fly)),false);
        end
        % save individual start / stop sleep timing (for later alignment
        % with distance and location within the arena)
        sleep(i).trial(trial).sleepON = sleepStartFrame;
        sleep(i).trial(trial).sleepOFF = sleepStoppedFrame;
        sleep(i).trial(trial).sleepLoc = startLoc;
        
        % save trial information into larger matrix
        temp = sleepStartFrame(:);
        temp(isnan(temp)) = [];
        sleepON = autoCat(sleepON,temp,false);
        temp = sleepStoppedFrame(:);
        temp(isnan(temp)) = [];
        sleepOFF = autoCat(sleepOFF,temp,false);
        sleepLength = autoCat(sleepLength,total_duration,false);
        
    end
    
    % save into sleep structure
    sleep(i).sleepON = sleepON;
    sleep(i).sleepOFF = sleepOFF;
    sleep(i).sleepLength = sleepLength;
    sleep(i).durationMatrix = timing;

end

if ~exist([saveDir 'Sleep'],'dir')
    mkdir([saveDir 'Sleep'])
end

clearvars('-except',initial_vars{:})
disp('All finished')

%% FIGURE: Sleeping over time
clearvars('-except',initial_vars{:})

plot_err = false;
LW = 1.5;
sb(1).idx = 1;
sb(2).idx = 2:3;
r = 3;
c = 1;
sSpan = 180; %1 minute smoothing

fig = getfig('',1);
for i = 1:num.exp

    subplot(r,c,sb(1).idx); hold on
        time = grouped(i).time;
        temp = grouped(i).temp;
        kolor = grouped(i).color;
        plot(time, temp, 'color',kolor, 'linewidth', LW)
        ylabel('temp (\circC)')
        % xlim([0,365])
    subplot(r,c,sb(2).idx); hold on
        y = smooth(sleep(i).sleepfract_avg,sSpan,'moving');
        y_err = smooth(sleep(i).sleepfract_err,sSpan,'moving');
        plot_error_fills(plot_err, time, y, y_err, kolor,  fig_type, 0.4);
        plot(time,y,'color',kolor,'linewidth',LW)
        ylabel('fraction of flies sleeping')
        xlabel('time (min)')
        % xlim([0,365])
end

formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca,'xcolor','none') 
subplot(r,c,sb(2).idx)
legend(expNames,'box','off')

save_figure(fig,[saveDir 'Sleep\' expGroup ' sleep timecourse'],fig_type);

%% FIGURES: Sleeping over time with sleep temp tuning curve
clearvars('-except',initial_vars{:})
plot_err = false;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %sleeping flies timecourse
sb(3).idx = 3:c:r*c; %binned distance alignment

LW = 0.75;
sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true); 
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor)
    
    %number of flies sleeping over time
    subplot(r,c,sb(2).idx); hold on
        y = smooth(sleep(i).sleepfract_avg,sSpan,'moving');
        y_err = smooth(sleep(i).sleepfract_err,sSpan,'moving');
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.4);
        plot(x,y,'color',kolor,'linewidth',LW)
        ylabel('fraction flies sleeping')
        xlabel('time (min)')
%         if ~autoLim
%             ylim(num_lim)
%         end

    %temp dependent distance
    subplot(r,c,sb(3).idx); hold on
        x = sleep(i).temps;
        y = sleep(i).temp_all(:,1);
        y_err = sleep(i).temp_all(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
 
        plot(x,y,'color',kolor,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.4);
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx) 
ylabel('\circC')
set(gca,"XColor",backColor)
% distance
subplot(r,c,sb(2).idx) 
ylabel('fraction of flies sleeping')
xlabel('time (min)')
set(gca,"XColor",foreColor)
% temp-distance relationship 
subplot(r,c,sb(3).idx) 
ylabel('fraction of flies sleeping')
xlabel('temp (\circC)')
% if ~autoLim
%     ylim(num_temp_lim)
% end
% 
% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

save_figure(fig,[saveDir 'Sleep\' expGroup ' sleeping flies summary'],fig_type);

%% FIGURE: Flies sleeping divided by heating / cooling
plot_err = true;
equalLim = true;
LW = 1.5;
r = 1;
c = 2;
nMax = num.exp;

% FIGURE:
fig = getfig('',true); 
% AVG
subplot(r,c,1)
hold on
for i = 1:nMax
    kolor = grouped(i).color;
    x = sleep(i).temps;
    y = sleep(i).temp_all(:,1)*100;
    y_err = sleep(i).temp_all(:,2)*100;
    plot(x,y,'color',kolor,'linewidth',LW+1)
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
end
ylimits = ylim;
ylim([0,ylimits(2)])
xlabel('Temperature (\circC)')
ylabel('flies sleeping (%)')
% SEP HEAT / COOL
subplot(r,c,2)
hold on
for i = 1:nMax
    kolor = grouped(i).color;   
    for type = 1:2 %increasing | decreasing 
        switch type
            case 1
                section_type = 'increasing';
                l_style = '-';
            case 2
                section_type = 'decreasing';
                l_style = '--';
        end
        x = sleep(i).temps;
        y = sleep(i).(section_type)(:,1)*100;
        y_err = sleep(i).(section_type)(:,2)*100;
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',l_style)
    end
end
ylimits = ylim;
ylim([0,ylimits(2)])
xlabel('Temperature (\circC)')
ylabel('flies sleeping (%)')
    
% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
if equalLim
    fig = matchAxis(fig,true);
end
% ylim(num_temp_lim)

% save figure
save_figure(fig,[saveDir 'Sleep\Flies sleeping during heating and cooling'],fig_type);

%% FIGURE: Distance and sleep overlaid tuning curves
clearvars('-except',initial_vars{:})
% Set params:
LW = 1;
r = 1; %rows
c = 2; %columns
plot_err = false;
[foreColor,~] = formattingColors(blkbgd);
equalLim = true;

% Build figure:
fig = getfig('',true); 
for i = 1:num.exp
kolor = grouped(i).color;
    for type = 1:2 %increasing | decreasing 
        switch type
            case 2
                section_type = 'increasing';
%                 l_style = '-';
            case 1
                section_type = 'decreasing';
%                 l_style = '--';
        end
    
        % COOLING
        yyaxis left %Distance Data
        subplot(r,c,type);   hold on
        x = grouped(i).(section_type).temps;
        y = grouped(i).(section_type).avg;
        y_err = grouped(i).(section_type).err;
        loc = isnan(y) | isnan(y_err);% remove nans 
        y(loc) = []; x(loc) = []; y_err(loc) = [];
    
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle','-','Marker','none')
        set(gca,'ydir','reverse')
        
  
         
        yyaxis right
        % Sleeping Data
        x = sleep(i).temps;
        y = sleep(i).(section_type)(:,1);
        y_err = sleep(i).(section_type)(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',':','Marker','.','MarkerSize',10)

    end
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
[yleft,yright] = deal([]);
for type = 1:2
    subplot(r,c,type); 
    yyaxis left
        set(gca,'xcolor',foreColor, 'YColor',foreColor)
        yleft(:,type) = ylim;
        ylabel('food proximity (mm)')
        
    yyaxis right
        set(gca,'YColor',foreColor)
        yright(:,type) = ylim;
        ylabel('fraction of sleeping flies')
    xlabel('temp (\circC)')
end

if equalLim 
    for type = 1:2
        subplot(r,c,type); 
        yyaxis left
        ylim([min(yleft(1,:)), max(yleft(2,:))])
            yyaxis right
            ylim([min(yright(1,:)), max(yright(2,:))])
    end
end
% flip axis on right plot to 'match time'
subplot(r,c,2)
set(gca,'XDir','reverse')

% save figure
save_figure(fig,[saveDir 'Sleep\sleep and distance during heating and cooling'],fig_type);

%% FIGURES: Experiment separated overlay of distance to food and sleeping flies
clearvars('-except',initial_vars{:})
sleepColor = Color('cyan');
distColor = Color('orange');
LW = 1;
r = 1; %rows
c = 2; %columns
plot_err = true;
[foreColor,backColor] = formattingColors(blkbgd);
equalLim = true;

% Organize data
sSpan = 180; %1 minute smoothing

% Build figures:
for i = 1:num.exp
    
    fig = getfig('',true); 
    % kolor = grouped(i).color;
        for type = 1:2 %increasing | decreasing 
            switch type
                case 2
                    section_type = 'increasing';
    %                 l_style = '-';
                case 1
                    section_type = 'decreasing';
    %                 l_style = '--';
            end
        
            yyaxis left
            subplot(r,c,type);   hold on
            %Distance Data 
            x = grouped(i).(section_type).temps;
            y = grouped(i).(section_type).avg;
            y_err = grouped(i).(section_type).err;
            loc = isnan(y) | isnan(y_err);% remove nans 
            y(loc) = []; x(loc) = []; y_err(loc) = [];
            plot_error_fills(plot_err, x, y, y_err, distColor,  fig_type, 0.2);     
            plot(x,y,'color',distColor,'linewidth',LW+1,'linestyle','-','Marker','none')
            set(gca,'ydir','reverse')
            
            yyaxis right
            % Sleeping Data
            x = sleep(i).temps;
            y = sleep(i).(section_type)(:,1);
            y_err = sleep(i).(section_type)(:,2);
            loc = isnan(y)|isnan(y_err);
            x(loc) = [];
            y(loc) = [];
            y_err(loc) = [];
            plot_error_fills(plot_err, x, y, y_err, sleepColor,  fig_type, 0.2);
            plot(x,y,'color',sleepColor,'linewidth',LW+1,'linestyle',':','Marker','.','MarkerSize',10)
        end
    
        % FORMATING AND LABELS
        formatFig(fig,blkbgd,[r,c]);
        [yleft,yright] = deal([]);
        % decreasing temp formatting
        subplot(r,c,1); 
        title('Cooling','Color',foreColor)
        set(gca,'XDir','reverse'); % flip temp direction to match temporal sequence
        yyaxis left
        set(gca,'xcolor',foreColor, 'YColor',distColor,'TickDir','out')
        yleft(:,1) = ylim;
        ylabel('proximity to food (mm)')
        
        xlabel('temp (\circC)')
        yyaxis right
        set(gca,'YColor',backColor)
        yright(:,1) = ylim;
        % increasing temp formatting
        subplot(r,c,2); 
        title('Warming','Color',foreColor)
        yyaxis left
        set(gca,'YColor',backColor,'xcolor',foreColor, 'TickDir','out')
        yleft(:,2) = ylim;
        yyaxis right
        set(gca,'YColor',sleepColor)
        ylabel('fraction of flies sleeping')
        
        xlabel('temp (\circC)')
        yright(:,2) = ylim;

        if equalLim 
            for type = 1:2
                subplot(r,c,type); 
                yyaxis left
%                 ylim([min(yleft(1,:)), max(yleft(2,:))])
                ylim([12 28])
                yyaxis right
%                 ylim([min(yright(1,:)), max(yright(2,:))])
                ylim([-.05 .4])
            end
        end
        
    % save figure
    save_figure(fig,[saveDir 'Sleep\' expNames{i} ' sleep and distance during heating and cooling'],fig_type);
%         save_figure(fig,[saveDir expNames{i} ' sleep and distance during heating and cooling NO SLEEP'],fig_type);

end

%% FIGURE: Sleep vs thermal threat
% Can we easily predict behavior based on a simple cold-exposure metric?
clearvars('-except',initial_vars{:})
% quantify the cold exposure during a single ramp as sum of difference in
% temp from preferred temp: (might need to make it exponential????)

% Pull data
[sleepDuration, thermalThreat] = deal([]);
for i = 1:num.exp
    thermalThreat(i) = sleep(i).thermalThreat;
    sleepDuration = autoCat(sleepDuration,sleep(i).avg_quant',false);
end


% FIGURE:
SZ = 60;
buff = 2;
fig = getfig('',1,[600 680]); hold on
% for i = 1:num.exp
%     kolor = grouped(i).color;
%     if i == 1
%         kolor = Color('cyan');
%     end
%     y = sleepDuration(:,i);
%     y(isnan(y)) = [];
%     y_avg = median(y);
%     x = shuffle_data(linspace(thermalThreat(i)-buff,thermalThreat(i)+buff,length(y)));
% %     x = ones(1,length(y))*thermalThreat(i));
%     scatter(x,y,SZ,kolor,"filled","o")
%     plot([thermalThreat(i)-(buff*1.5),thermalThreat(i)+(buff*1.5)],[y_avg,y_avg],'Color',kolor,'linewidth',2)
% end

xlabel('Thermal Stress')
ylabel('Sleep duration per fly (sec)')
xlim([0, 70])
ylim([0, 1200])
formatFig(fig,blkbgd);
ax = gca;
set(ax, 'tickDir', 'in')
set(ax,'xcolor', 'none')

% set(ax,'ytick',0:300:1200)
% set(ax,'yticklabel',{'0','5','10',  '15', '20'})
% ylabel('Sleep duration per fly (min)')
% set(ax,'xcolor', 'w')
% save_figure(fig,[saveDir 'Sleep\Sleep duration by thermal stress empty axes'],'-png');


save_figure(fig,[saveDir 'Sleep\Sleep duration by thermal stress'],'-pdf',true,false);
save_figure(fig,[saveDir 'Sleep\Sleep duration by thermal stress'],'-png');



%% FIGURE: single group sleep duration per fly
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);

% Pull data
[sleepDuration, thermalThreat] = deal([]);
for i = 1:num.exp
    thermalThreat(i) = sleep(i).thermalThreat;
    sleepDuration = autoCat(sleepDuration,sleep(i).avg_quant',false);
end

% Single figures of the sleep duration per fly
SZ = 60;
buff = 0.3;

for i = 1:num.exp
    fig = getfig('',1,[386 680]); hold on
        kolor = grouped(i).color;
        y = sleepDuration(:,i);
        y(isnan(y)) = [];
        y_avg = mean(y);
        x = shuffle_data(linspace(1-buff,1+buff,length(y)));
        scatter(x,y,SZ,kolor,"filled","o")
        plot([1-buff,1+buff],[y_avg,y_avg],'Color',kolor,'linewidth',2)
        % Labels and formatting
        xlim([0,2])
        ylim([0, 2500])
        ylabel('Sleep duration per fly (sec)')
        ylim([0, 2500])
        formatFig(fig,blkbgd);
        set(gca,'XColor',backColor)
        xlabel(expNames{i},'Color',foreColor,'FontSize',12)
        save_figure(fig,[saveDir 'Sleep\' expNames{i} ' Sleep duration per fly'],fig_type);
end

%% FIGURE: FIX THIS! Avg quantity of sleep per fly

% Can we easily predict behavior based on a simple cold-exposure metric?
clearvars('-except',initial_vars{:})
% quantify the cold exposure during a single ramp as sum of difference in
% temp from preferred temp: (might need to make it exponential????)


sleepDuration = [];
for i = 1:num.exp
    % Thermal threat
    tPoints = getTempTurnPoints(data(i).temp_protocol);
    demoRamp_idx = tPoints.down(1,1):tPoints.up(1,2);

    temp_ramp = grouped(i).temp(demoRamp_idx);
    thermalThreat(i) = sum(25-temp_ramp)/(fps*rampDur);
    
    % Sleep duration
    nFlies = zeros(1,num.trial(i));
    nRamps = size(tPoints.up,1);
    for ramp = 1:nRamps
        rampROI = tPoints.down(ramp,1):tPoints.up(ramp,2);
        nFlies = nFlies + sum(sleep(i).num(rampROI,:),1);
    end
    avg_sleepDuration = ((nFlies/nRamps)./(data(i).T.NumFlies'))./fps; % avg fly sleep duration (in sec) during a single temp ramp
    sleepDuration = autoCat(sleepDuration,avg_sleepDuration',false);

end


% FIGURE:
SZ = 60;
buff = 0.05;
fig = getfig('',1,[600 680]); hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    y = sleepDuration(:,i);
    y(isnan(y)) = [];
    y_avg = mean(y);
    x = shuffle_data(linspace(thermalThreat(i)-buff,thermalThreat(i)+buff,length(y)));
%     x = ones(1,length(y))*thermalThreat(i));
    scatter(x,y,SZ,kolor,"filled","o")
    plot([thermalThreat(i)-buff,thermalThreat(i)+buff],[y_avg,y_avg],'Color',kolor,'linewidth',2)
end

xlabel('Thermal Stress Rate')
ylabel('Sleep duration per fly (sec)')
xlim([0,3])
formatFig(fig,blkbgd);

save_figure(fig,[saveDir 'Sleep\Sleep duration by thermal stress rate'],fig_type);

%% FIGURE: Sleeping distance to food across during cooling and warming
clearvars('-except',initial_vars{:})
LW = 2;
r = 1; %rows
c = 2; %columns
plot_err = true;
[foreColor,~] = formattingColors(blkbgd);
equalLim = true;

% where are the flies choosing to sleep? Distance from food? --> bin by temp
yLimits = [];
fig = getfig('',1); 
for tt = 1:2 %increasing | decreasing 
    switch tt
        case 2
            section_type = 'increasing';
            axis_dir = 'normal';
        case 1
            section_type = 'decreasing';
            axis_dir = 'reverse';
    end
    subplot(r,c,tt); hold on
    for i = 1:num.exp
        kolor = grouped(i).color;
        x = sleep(i).temps;
        y = sleep(i).([section_type '_dist'])(:,1);
        y_err = sleep(i).([section_type '_dist'])(:,2);
        
        loc = isnan(y) | isnan(y_err); % remove nans 
        y(loc) = []; x(loc) = []; y_err(loc) = [];

        % plot data
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
        plot(x,y,'color',kolor,'linewidth',LW,'linestyle','-','Marker','none')
        
    end
    % formatting and labeling
    set(gca,'ydir','reverse','xdir',axis_dir)
    xlabel('temp (\circC)')
    yLimits(:,tt) = ylim;
    title(section_type)
%     axis tight
end
% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
subplot(r,c,1)
ylabel('distance to well (mm)')
if equalLim
    ylim([min(min(yLimits)), max(max(yLimits))])
end
subplot(r,c,2)
if equalLim
    ylim([min(min(yLimits)), max(max(yLimits))])
end
set(gca,'YColor',foreColor)


% save figure
save_figure(fig,[saveDir 'Sleep\sleeping distance to food'],fig_type);

%% FIGURE: Histogram of sleep proximity to food
clearvars('-except',initial_vars{:})
fig_dir = [saveDir 'Sleep\sleeping histograms overlay\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
binedges = 0:3:51;
% ylimits = [-50 500]; % [-50,1400] 15-23 WT
ylimits = [-50 800];

for idx = 1:num.exp
    fig = getfig('',1); hold on
    for ii = 1:idx
        i = expOrder(ii);
        plotData = [];
        for trial = 1:num.trial(i)
            loc = sleep(i).trial(trial).sleepLoc;
            test = sleep(i).trial(trial).all_distance(loc);
            plotData = [plotData; test];
        end
        h(i) = histogram(plotData,binedges,'FaceColor',grouped(i).color,'FaceAlpha',0.7);
    end
    formatFig(fig,blkbgd);
    xlabel('distance to well (mm)','FontSize',20)
    ylabel('sleeping flies (#)','FontSize',20)
    set(gca,'TickDir','out')
    xlim([-1.5000   51.5000])
    ylim(ylimits)
    
    % save figure
    save_figure(fig,[fig_dir 'sleeping distance to food histogram ' num2str(idx)],fig_type);
end
%% FIGURE: TODO Sleeping distance histogram normalized by amount of sleeping....
clearvars('-except',initial_vars{:})
fig_dir = [saveDir 'Sleep\normalized sleeping histograms overlay\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

% load([baseFolder 'Fundamentals\Distance_to_well_distribution']);
initial_vars{end+1} = 'null_dist';

binedges = 0:2:51;
LW = 3;
ylimits = [-50 800];

% null_distribution for normalization
[N_null,edges] = histcounts(null_dist.distance,binedges,'Normalization','PDF');
% bar(edges(1:end-1),N,'FaceColor','y','FaceAlpha',0.3)

% for idx = 1:num.exp
idx = num.exp;
    fig = getfig('',1); hold on
    plot(bin_center, N_null,'color', Color('grey'),'linewidth', LW)
    for ii = 1:idx
        i = expOrder(ii);
        plotData = [];
        for trial = 1:num.trial(i)
            loc = sleep(i).trial(trial).sleepLoc;
            test = sleep(i).trial(trial).all_distance(loc);
            plotData = [plotData; test];
        end

        % [N,edges] = histcounts(plotData,binedges,'Normalization','probability');
        [N,edges] = histcounts(plotData,binedges,'Normalization','PDF');

        bin_center = (binedges(1:end-1) + binedges(2:end))/2;
        plot(bin_center, N,'color', grouped(i).color,'linewidth', LW)

        % N_test = N./N_null;
        % bar(edges(2:end),N,'FaceColor',grouped(i).color,'FaceAlpha',0.3)
        % bar(edges(1:end-1),N_test,'FaceColor',grouped(i).color,'FaceAlpha',0.3)
       
    end
    

    formatFig(fig,blkbgd);
    xlabel('distance to well (mm)','FontSize',20)
    ylabel('sleep occupancy probability','FontSize',20)
    set(gca,'TickDir','out')
    % xlim([-1.5000   51.5000])
    % ylim(ylimits)
     save_figure(fig,[fig_dir 'sleeping occupancy vs distance to food'],fig_type);

    % % save figure
    % save_figure(fig,[fig_dir 'sleeping occupancy vs distance to food ' num2str(idx)],fig_type);
% end



% FIGURE: probability difference in sleeping position from null distribution

N_null = histcounts(null_dist.distance,binedges,'Normalization','PDF');
bin_center = (binedges(1:end-1) + binedges(2:end))/2;

fig = getfig('',1); hold on
    for i = 1:num.exp
        plotData = [];
        for trial = 1:num.trial(i)
            loc = sleep(i).trial(trial).sleepLoc;
            test = sleep(i).trial(trial).all_distance(loc);
            plotData = [plotData; test];
        end

        N = histcounts(plotData,binedges,'Normalization','PDF');
        N = N-N_null;
        plot(bin_center, N,'color', grouped(i).color,'linewidth', LW)
    end
    h_line(0,'grey','--',2)
    ylim([-.05,0.15])
    formatFig(fig,blkbgd);
    xlabel('distance to well (mm)','FontSize',20)
    ylabel('sleep occupancy diff in probability from null','FontSize',20)
    set(gca,'TickDir','out')
save_figure(fig,[fig_dir 'sleeping null diff in occupancy vs distance to food'],fig_type);


%% ANALYSIS: TODO fly sleeping locations in the arena (oriented to same heading)
clearvars('-except',initial_vars{:})
% HAVENT SCREENED FOR ONLY RAMPSSSSSSSS!!! 
% TODO : 
% change the orientation/rotation to match across trials
% color code the dot to indicate the temperature for each sleeping fly
% location


for i = 1:num.exp
    all_data = [];
    temperature = grouped(i).temp; % temperatures for each time index point
    for trial = 1:num.trial(i)
        [x_pos,y_pos,temp_pos,idx_pos,plotData] = deal([]);
        % fly locations (concatenate across all the 'fly' tracks)
        x = sleep(i).trial(trial).X;
        y = sleep(i).trial(trial).Y;
        for idx = 1:size(x,2) 
            loc = sleep(i).trial(trial).sleepON(:,idx);
            loc(isnan(loc)) = [];
            x_pos = [x_pos; x(loc,idx)];
            y_pos = [y_pos; y(loc,idx)];
            temp_pos = [temp_pos; temperature(loc)];
            idx_pos = [idx_pos; loc];
        end
        wells = data(i).data(trial).data.wellcenters;
        wellLoc = data(i).T.foodLoc(trial);
        % Make food well the origin
        x_offset = wells(1,wellLoc);
        y_offset = wells(2,wellLoc);
        wells_x = wells(1,:)-x_offset;
        wells_y = wells(2,:)-y_offset;
        X = x_pos-x_offset;
        Y = y_pos-y_offset;
        
        % Rotate to correct orientation
        switch wellLoc
            case 1
                plotData(:,1) = Y;
                plotData(:,2) = -X;
                WELLS(:,1) = wells_y;
                WELLS(:,2) = -wells_x;
            case 2 
                plotData(:,1) = X;
                plotData(:,2) = -Y;
                WELLS(:,1) = wells_x;
                WELLS(:,2) = -wells_y;
            case 3
                plotData(:,1) = -Y;
                plotData(:,2) = X;
                WELLS(:,1) = -wells_y;
                WELLS(:,2) = wells_x;
            case 4 
                plotData(:,1) = X;
                plotData(:,2) = Y;
                WELLS(:,1) = wells_x;
                WELLS(:,2) = wells_y;
        end
        % save the data into larger group structure
        all_data = [all_data; plotData, temp_pos,idx_pos];
    end
    sleep(i).location.all_data = all_data;
    sleep(i).location.foodWell = [0,0];
    sleep(i).location.arenaCenter = [WELLS(5,1),WELLS(5,2)];
    sleep(i).location.r = data(i).data(trial).data.r;
    sleep(i).location.wells = WELLS;
end

% NULL DISTRIBUTION (this should be identical for all trials/experiments in the same arena)
% Generate the null distribution of distances to the food well:
nbins = 50; %This is same as used for determining sleep positions
pixel_buffer = 50; % edge of arena pixel buffer
i = 1;
% Set axis limits for the selected arena
c1 = data(i).data(trial).data.centre(1);
c2 = data(i).data(trial).data.centre(2);
r = data(i).data(trial).data.r;
xlimit = [c1-(r+pixel_buffer),c1+(r+pixel_buffer)];
ylimit = [c2-(r+pixel_buffer),c2+pixel_buffer+r];

% find the 'auto bin' lines
xedge = linspace(xlimit(1),xlimit(2),nbins+1);
yedge = linspace(ylimit(1),ylimit(2),nbins+1);
X = []; Y = [];
idx = 1;
for y_loc = 1:nbins
    for x_loc = 1:nbins
        X(idx) = xedge(x_loc);
        Y(idx) = yedge(y_loc);
        idx = idx + 1;
    end
end
% screen out the units outside the arena circle
temp_dist = sqrt((X-c1).^2 + (Y-c2).^2);
loc = temp_dist>r;
X(loc) = [];
Y(loc) = [];

% find food well location for distance capture later 
%  the specific well does not matter for the null distribution because the arena is
%  symmetrical
foodWellLoc = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));
c1 = foodWellLoc(1);
c2 = foodWellLoc(2);
null_distance = sqrt((X-c1).^2 + (Y-c2).^2)./pix2mm;

initial_vars{end+1}  = 'null_distance';
clearvars('-except',initial_vars{:})

%% FIGURES: fly sleeping positions in the arena
clearvars('-except',initial_vars{:})
% HAVENT SCREENED FOR ONLY RAMPSSSSSSSS!!! 
% TODO HERE...
fig_dir = [saveDir 'Sleep\sleeping locations\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

[foreColor,backColor] = formattingColors(blkbgd);
SZ = 40;
% Plot scatter point figure:
for i = 1:num.exp
    all_data = sleep(i).location.all_data;
    kolor = grouped(i).color;
    
    fig = getfig('',1); hold on
  
    % PLOT
    scatter(sleep(i).location.wells(1:4,1),sleep(i).location.wells(1:4,2),SZ+20,Color('grey'),'filled')
    scatter(0,0,SZ+20,'green','filled')
    scatter(all_data(:,1),all_data(:,2),15,kolor,'filled')
   
    viscircles([sleep(i).location.arenaCenter(1),sleep(i).location.arenaCenter(2)],sleep(i).location.r,'Color',grouped(i).color);
    axis square; axis equal;  
    formatFig(fig,blkbgd);
    set(gca,'XColor',backColor,'YColor',backColor);
    save_figure(fig,[fig_dir grouped(i).name ' sleeping positions in arena'],fig_type);

end


% FIGURE: Sleeping position heatmap
nBins = 15;
for i = 1:num.exp
    c1 = sleep(i).location.arenaCenter(1);
    c2 = sleep(i).location.arenaCenter(2);
    r = sleep(i).location.r;

    xlimit = [c1-(r+50),c1+(r+50)];
    ylimit = [c2-(r+50),c2+(50+r)];

    % find the 'auto bin' lines
    xedge = linspace(xlimit(1),xlimit(2),nBins+1);
    yedge = linspace(ylimit(1),ylimit(2),nBins+1);

    % pull the fly locations during the trial
    x_loc = sleep(i).location.all_data(:,1);
    y_loc = sleep(i).location.all_data(:,2);
    N = histcounts2(x_loc,y_loc,xedge,yedge);
    newN = imrotate(N,90);

    bw = median(diff(yedge)); %bin-width

    fig = getfig('',1); hold on
        imagesc([xlimit(1)+bw, xlimit(end)-bw], [ylimit(end)-bw, ylimit(1)+bw], newN)
        hold on
    
        scatter(sleep(i).location.wells(1:4,1),sleep(i).location.wells(1:4,2),SZ+20,'k','filled')
        scatter(0,0,SZ+20,'red','filled')
    %     scatter(sleep(i).location.all_data(:,1),sleep(i).location.all_data(:,2),15,'w','filled')
        viscircles([sleep(i).location.arenaCenter(1),sleep(i).location.arenaCenter(2)],sleep(i).location.r,'Color',grouped(i).color);
        axis square; axis equal;  
        c = colorbar;
        formatFig(fig,blkbgd);
        c.Color = foreColor;
        set(gca,'xcolor',backColor,'ycolor',backColor)
        c.Label.String = '# Flies';
        c.Label.Color = foreColor;
        c.Label.FontSize = 20;
        save_figure(fig,[fig_dir grouped(i).name ' sleeping positions in arena heatmap'],fig_type);

end

%% FIGURE: histogram of sleeping position by temperature...
clearvars('-except',initial_vars{:})
LW = 2;
r = 1; %rows
c = 2; %columns
plot_err = true;
[foreColor,~] = formattingColors(blkbgd);
equalLim = true;

temp_limits = [nan nan];
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    temp_limits = [min([tp.threshLow, temp_limits(1)]), max([tp.threshHigh, temp_limits(2)])];
end
nColors = floor(temp_limits(1)):1:ceil(temp_limits(2));
nReds = floor(length(nColors)/2);
nBlues = length(nColors) - nReds + 1;

cBlues = Color('DodgerBlue', 'lightyellow',nBlues);
cReds = Color('lightyellow', 'coral', nReds);
cMap = [cBlues(1:nBlues-1,:);cReds];
cMap(cMap==1) = 1;

% Need to make a stacked histogram....
distance_bins = 0:3:51;
for i = 1:num.exp
    all_data = sleep(i).location.all_data;
    all_data(:,5) = sqrt((all_data(:,1)).^2 + (all_data(:,2).^2))./pix2mm; %distance to food for sleeping flies
    all_data(:,6) = discretize(all_data(:,3),nColors); % temp bins for each data point
    nanLoc = isnan(all_data(:,6));
    all_data(nanLoc,:) = []; %remove data outside the temp constraints
    all_data(:,7) = discretize(all_data(:,5),distance_bins); % distance bins
    
    % Sort instances into distance then temp bins
    plotData = nan(length(distance_bins),length(nColors));
    for d = 1:length(distance_bins)
        for t = 1:length(nColors)
            plotData(d,t) = sum(all_data(:,6)==t & all_data(:,7)==d);
        end
    end
    
    x = distance_bins(1:end-1) + median(diff(distance_bins)/2);

    fig = getfig('',1); hold on
    h = bar(distance_bins,plotData,'stacked');
    for k = 1:length(nColors)
        try
            h(k).FaceColor = cMap(k,:);%grouped(i).color; %
        catch 
            h(k).FaceColor = 'r';
        end
        h(k).FaceAlpha = 1;
        h(k).EdgeColor = 'none';
    end
    % formats and labels
    formatFig(fig,blkbgd);
    xlabel('distance to well (mm)','FontSize',20)
    ylabel('sleeping flies (#)','FontSize',20)
    set(gca,'TickDir','out')
    % colorbar formatting
    colormap(cMap);
    c = colorbar; 
    c.Color = foreColor;
    clim([nColors(1), nColors(end)])
    c.Label.String = 'Temperature (\circC)';

    % normalized:
    ylim([0,800])

% save figure
save_figure(fig,[saveDir 'Sleep\' expNames{i} ' sleeping distance to food histogram temp colored normalized'],fig_type);
% save_figure(fig,[saveDir expNames{i} ' sleeping distance to food histogram group colored normalized'],fig_type);

end

fig = figure;
colormap('default'); %return map to normal
close(fig)

%% FIGURE: sleeping distance (binned warming and cooling)
clearvars('-except',initial_vars{:})
LW = 2;
plot_err = true;

% where are the flies choosing to sleep? Distance from food? --> bin by temp
fig = getfig('',1); hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    x = sleep(i).temps;
    y = sleep(i).tempBinDist(:,1);
    y_err = sleep(i).tempBinDist(:,2);
    
    loc = isnan(y) | isnan(y_err); % remove nans 
    y(loc) = []; x(loc) = []; y_err(loc) = [];

    % plot data
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.2);
    plot(x,y,'color',kolor,'linewidth',LW,'linestyle','-','Marker','none')
end
% formatting and labeling
set(gca,'ydir','reverse','TickDir','out')
xlabel('temp (\circC)')
ylabel('distance to well (mm)')

% FORMATING AND LABELS
formatFig(fig,blkbgd);

% save figure
save_figure(fig,[saveDir 'Sleep\avg sleeping distance to food'],fig_type);

%% FIGURE: Sleep duration and statistics across groups
[foreColor,~] = formattingColors(blkbgd);

% How does the avg duration of sleep change across fly groups? 
[y, kolor,gLabel] = deal([]);
for idx = 1:num.exp
    i = expOrder(idx);
    kolor(idx,:) = grouped(i).color;
    temp = sleep(i).sleepLength(:);
    temp(isnan(temp)) = [];
    y = autoCat(y,sleep(i).sleepLength(:),false);
%     gLabel{idx} = expNames{i};
end

fig = getfig('',1);
[h,L,MX,MED,bw] = violin(y,'xlabel',{''},'facecolor',kolor,'edgecolor',foreColor,'alp',1,...
       'bw',0.35,'mc',foreColor,'medc','r--');

ylabel('sleep duration (min)','FontSize',18)
% set(gca,'xcolor',backColor)
set(L, 'TextColor',foreColor);
formatFig(fig, blkbgd);

ylim([-2,20])

save_figure(fig,[saveDir 'Sleep duration across groups violin plot'],fig_type);

% FIGURE: cumulative distribution function of sleep
% clearvars('-except',initial_vars{:})
% autoSave = true;
autoLim = false;
xlimit = [0,20];
LW = 2;

fig = getfig('',1); hold on
for i = 1:num.exp
    h = cdfplot(y(:,i));
    set(h,'color',grouped(i).color,'linewidth',LW)
end

xlabel('Sleep duration (min)')
ylabel('Empirical CDF')
title('')
formatFig(fig, blkbgd);   
set(gca,'xgrid','off','ygrid','off')
set(gca,'ytick',0:0.2:1)
if ~autoLim
    xlim(xlimit)
end
save_figure(fig,[saveDir 'Sleep\' expGroup ' sleep CDF'],fig_type);

% % Histogram version...
% fig = getfig('',1); hold on
% for i = 1:num.exp
%     kolor = grouped(i).color;
%     histogram(sleep(i).sleepLength(:),'FaceColor',kolor)
% end
% formatFig(fig,true);
% xlabel('sleep duration (min)')
% ylabel('count')
% 
% % Could make this a sleep PDF?
% 
% median(sleep(1).sleepLength(:),'omitnan')
% median(sleep(2).sleepLength(:),'omitnan')

%% FIGURE: Sleep vs thermal threat by genotype
% Can we easily predict behavior based on a simple cold-exposure metric?
clearvars('-except',initial_vars{:})
% quantify the cold exposure during a single ramp as sum of difference in
% temp from preferred temp: (might need to make it exponential????)

% expList = {'Berlin WT', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
% colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
% kolor = Color('White');

% Pull data
[sleepDuration, thermalThreat] = deal([]);
for i = 1:num.exp
    thermalThreat(i) = sleep(i).thermalThreat;
    sleepDuration = autoCat(sleepDuration,sleep(i).avg_quant',false);
end


% FIGURE:
SZ = 60;
buff = 2;
fig = getfig('',1,[600 680]); hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    y = sleepDuration(:,i);
    y(isnan(y)) = [];
    y_avg = mean(y);
    x = shuffle_data(linspace(thermalThreat(i)-buff,thermalThreat(i)+buff,length(y)));
%     x = ones(1,length(y))*thermalThreat(i));
    scatter(x,y,SZ,kolor,"filled","o")
    plot([thermalThreat(i)-buff,thermalThreat(i)+buff],[y_avg,y_avg],'Color',kolor,'linewidth',2)
end

xlabel('Thermal Stress')
ylabel('Sleep duration per fly (sec)')
xlim([0, 70])
ylim([0, 2500])
formatFig(fig,blkbgd);

save_figure(fig,[saveDir 'Sleep\Sleep duration by thermal stress norm axes'],fig_type);

% save_figure(fig,[saveDir 'Sleep duration by thermal stress norm axes short'],fig_type);

% 

%% FIGURE: portion of sleep during warming cooling and hold
clearvars('-except',initial_vars{:})
buffer = 0.15;
[foreColor,~] = formattingColors(blkbgd);

SZ = 50;

for i = 1:num.exp
    all_roi = [];
    fig = getfig('',1,[481 680]); hold on
    tPoints = getTempTurnPoints(data(i).temp_protocol);
    % TOTALSLEEP = sum(sleep(i).num);
    
     for type = 1:3
        switch type
            case 1 %decreasing
                roi = tPoints.DownROI;
            case 2 % increasing
                roi = tPoints.UpROI;
            case 3 % holding
                roi = tPoints.HoldROI;
        end
        all_roi = autoCat(all_roi, roi,false);
     end
    TOTALSLEEP = sum(sleep(i).num(all_roi,:));

    for type = 1:3
        switch type
            case 1 %decreasing
                roi = tPoints.DownROI;
                kolor = Color('dodgerblue');
            case 2 % increasing
                roi = tPoints.UpROI;
                kolor = Color('red');
            case 3 % holding
                roi = tPoints.HoldROI;
                kolor = Color('grey');
        end

        sleepCount = sleep(i).num(roi,:);
        
        % normalize sleep to TOTAL sleep
        frames_of_sleep = sum(sleepCount);
        sleepPortion = (frames_of_sleep./TOTALSLEEP)*100;
        sleepPortion(isnan(sleepPortion)) = [];

        % Plot data
        %average
        bar(type, mean(sleepPortion),'FaceColor',kolor,'FaceAlpha',1,'EdgeColor',foreColor)
        %scatterpoints
        x = shuffle_data(linspace(type-buffer,type+buffer,length(sleepPortion)));
        scatter(x,sleepPortion,SZ,foreColor,"filled")

        % hold the data
        plotData(:,type) = sleepPortion;
    end
    % formatting and labels
    set(gca,'XTick',1:3,'XTickLabel',{'cooling','heating','hold'})
    ylabel('Portion of total sleep (%)')
    xlim([0,4]) 
    ylim([0,100])
    formatFig(fig,blkbgd);
    
    % Save figure
    save_figure(fig,[saveDir 'Sleep\' expNames{i} ' sleeping by temp type'],fig_type);
    
end

%% FIGURE: probabiltiy of sleep in each temperature regime
% uses Bayes' Formula to calculate the probability of sleep events
% Formula: P(W|S) = P(S|W)P(W) / (P(S|W)P(W) + P(S|C)P(C) + P(S|H)P(H))
% aka: probability of it being 'warming' if a fly is sleeping

clearvars('-except',initial_vars{:})
buffer = 0.15;
[foreColor,~] = formattingColors(blkbgd);

SZ = 50;

for i = 1:num.exp
    plotData = [];
    all_roi = [];
    tPoints = getTempTurnPoints(data(i).temp_protocol);
    
    roiC = tPoints.DownROI;
    roiW = tPoints.UpROI;
    roiH = tPoints.HoldROI;
    all_roi =[roiC,roiW,roiH];
            
    % variables needed to know:
    pW = length(roiW)/length(all_roi); % prob of warming
    pC = length(roiC)/length(all_roi); %prob of cooling
    pH = length(roiH)/length(all_roi); %prob of holding
    
    pSgW = sum(sleep(i).num(roiW,:))/length(all_roi); %prob of sleep given warming condition
    pSgC = sum(sleep(i).num(roiC,:))/length(all_roi); %prob of sleep given cooling condition
    pSgH = sum(sleep(i).num(roiH,:))/length(all_roi); %prob of sleep given holding condition

    % Calculate the probabilities
    probSW = pSgW.*pW ./ ( pSgW.*pW + pSgC.*pC + pSgH.*pH);
    probSC = pSgC.*pC ./ ( pSgW.*pW + pSgC.*pC + pSgH.*pH);
    probSH = pSgH.*pH ./ ( pSgW.*pW + pSgC.*pC + pSgH.*pH);
    
    plotData = [probSC',probSW',probSH']; %cool, warm, hold ploting
    
    % plot the probabilities
    fig = getfig('',1,[481 680]); hold on
     for type = 1:3
        switch type
            case 1 %decreasing
                kolor = Color('dodgerblue');
            case 2 % increasing
                kolor = Color('red');
            case 3 % holding
                kolor = Color('grey');
        end
        %average
        bar(type, mean(plotData(:,type),'omitnan'),'FaceColor',kolor,'FaceAlpha',1,'EdgeColor',foreColor)
        %scatterpoints
        x = shuffle_data(linspace(type-buffer,type+buffer,num.trial(i)));
        scatter(x,plotData(:,type),SZ,foreColor)
     end

    % formatting and labels
    set(gca,'XTick',1:3,'XTickLabel',{'cooling','heating','hold'})
    ylabel('Prob of sleep for each temp condition')
    xlim([0,4])
    ylim([0,1])
    formatFig(fig,blkbgd);
    title(grouped(i).name)
    % Save figure
    save_figure(fig,[saveDir 'Sleep\' expNames{i} 'Bayes probabilty of sleeping by temp type'],fig_type);
end




%% FIGURE: Onset of sleep raster plot 
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);

% parameters
spike_H = 2;    %height of each raster
trial_space = 1;    %gap between trial lines
spike_W = 0.5;    % raster line width

fig_height = 50 + num.exp * 100;
fig = getfig('',1,[1064 fig_height]); hold on

y_base = 1;
for ii = 1:num.exp
    i = expOrder(ii);
    time = grouped(i).time;
    kolor = grouped(i).color;
    
    for trial = 1:num.trial(i)
        plotdata = sleep(i).sleepON(:,trial);
        plotdata(isnan(plotdata)) = []; %remove dummy spaces
        plotdata = sort(plotdata);
        % remove any simultaneous starts (shift one point left or right)
        idx = 0;
        while any(diff(plotdata)==0)
            loc = find(diff(plotdata)==0);
            plotdata(loc+1) = plotdata(loc)+1;
            idx = idx+1;
            if idx>(fps*12)
                warndlg('Points shifting more than 10 seconds')
                return 
            end
        end
        % get time for each sleep start
        plotdata(plotdata>length(time)) = []; %remove those outside the universal time
        x = time(plotdata);
        X = [x';x'];
        % get y bin slot locations
        Y = repmat([y_base;y_base+spike_H],[1,size(X,2)]);

        % PLOT SPIKES
        plot(X,Y,'color',kolor,'linewidth',spike_W)
        
        % Update y-offset
        y_base = y_base+spike_H+trial_space;
    end
    % plot ramp times...
    y_base = y_base+trial_space;
    tp = getTempTurnPoints(data(i).temp_protocol);
    y = [y_base+(spike_H/2);y_base+(spike_H/2)];
    
    % if i == num.exp
        plot(time(tp.down)',repmat(y,[1,tp.nDown]),'color',Color('dodgerblue'),'linewidth',spike_H*2) %decreasing
        plot(time(tp.up)',repmat(y,[1,tp.nUp]),'color',Color('red'),'linewidth',spike_H*2) %increasing
        plot(time(tp.hold)',repmat(y,[1,tp.nHold]),'color',Color('grey'),'linewidth',spike_H*2) %decreasing
    % end
    y_base = y(2)+(2+trial_space);
end

% Plot 
ylim([-10,y_base+3])
formatFig(fig,blkbgd);
set(gca,'TickDir','out')
xlabel('time (min)')
set(gca,'ycolor',backColor)

% Save figure
save_figure(fig,[saveDir 'Sleep\' 'sleeping onset raster'],fig_type);    

%% FIGURE: onset to sleeping after warming start
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);
LW = 1.5;
autoSave = true;
fig_dir = [saveDir 'Sleep\' 'Sleep Onset Frequency\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

for idx = 1:num.exp

    fig = getfig('',1,[1064, 305]); hold on
    for ii = 1:idx
        i = expOrder(ii);
        time = grouped(i).time;
        kolor = grouped(i).color;
        tp = getTempTurnPoints(data(i).T.TempProtocol{1});
        
        plotdata = sleep(i).sleepON(:);
        plotdata(isnan(plotdata)) = [];
        x = sort(time(plotdata));
        binedges = 0:ceil(time(end));
    
        N = discretize(x,binedges);
        freq = zeros(length(binedges),1);
        for t = 1:length(binedges)
            freq(t) = sum(N==t);
        end
        freq = freq./sum(data(i).T.NumFlies(:));
        plot(binedges,freq,'color',kolor,'linewidth',LW)
    end
    formatFig(fig,blkbgd);
    xlabel('time (min)')
    ylabel('sleep frequency (fly/min)') 
    set(gca,'TickDir','out')
    % legend({grouped(:).name},'color', foreColor)
    % xlim([0, 400])
    % ylim([0, 0.3])

    % Save figure
    save_figure(fig,[fig_dir 'sleeping onset frequency ' num2str(idx)],fig_type,autoSave);   

end


%% FIGURE: event-aligned sleep onset
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
fig_height = 249*num.exp;
timeBins = 0:3:150; % 3 minute bins

% Pull data together
y = struct;
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).temp_protocol);
    rampDur = min(tp.up(:,2)-tp.up(:,1)); %shortest ramp duration for overlay
    plotdata = [];
    sleepData = sleep(i).sleepON;
    % Heating onset focus
    for n = 1:tp.nUp
        loc = (sleepData>=tp.up(n,1)) & (sleepData<=tp.up(n,2));
        sleepStarts = (sleepData(loc)-tp.up(n,1)); %frame offset from heating event
        plotdata = [plotdata; sleepStarts];
    end
    plotdata = plotdata/(fps*60); %time in minutes
    y(i).data = plotdata;
    y(i).rampDur = rampDur/(fps*60);
end


% FIGURE 1 Plot histogram distribution
fig = getfig('',1,[1064 fig_height]); 
for ii = 1:num.exp
    subplot(num.exp,1,ii)
    i = expOrder(ii);
    histogram(y(i).data,timeBins, 'FaceColor',grouped(i).color,'FaceAlpha',1,'EdgeColor',backColor)
    v_line(5,foreColor)
end
formatFig(fig,blkbgd,[num.exp,1]);
for i = 1:num.exp-1
    subplot(num.exp,1,i)
    set(gca,'XColor',foreColor,'xtick',[])
    ylimits = ylim;
    ylim([-5,ylimits(2)])
    ylabel('sleep onset count')
end
subplot(num.exp,1,num.exp)
ylimits = ylim;
ylim([-5,ylimits(2)])
xlabel('time since heating onset (min)')
save_figure(fig,[saveDir expGroup ' sleep onset distribution'],fig_type);

% FIGURE 2 Plot cumulative distribution function
fig = getfig('',1); hold on
for i = 1:num.exp
    h = cdfplot(y(i).data);
    set(h,'color',grouped(i).color,'linewidth',2)
end

xlabel('Sleep duration (min)')
ylabel('Empirical CDF')
title('')
formatFig(fig, blkbgd);   
set(gca,'xgrid','off','ygrid','off')
set(gca,'ytick',0:0.2:1)

save_figure(fig,[saveDir 'Sleep\' expGroup ' sleep onset CDF'],fig_type);


%% FIGURE: avg sleep per hour per fly
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
buffer = 0.25;
SZ = 50;
LW = 2;
total_sleep = [];
fig_H = 220 + (100*num.exp);

fig = getfig('',1,[fig_H,680]); hold on
for ii = 1:num.exp
    i = expOrder(ii);
    tp = getTempTurnPoints(data(i).temp_protocol);
    fps = tp.fps;
    roi = [tp.DownROI, tp.UpROI, tp.HoldROI];
    roi_min = length(roi)/(fps*60);
    y = sum(sleep(i).num(roi,:));
    y = y ./ data(i).T.NumFlies'; % frames of sleep per fly
    y = y ./ fps; % seconds of sleep per fly
    y = y ./ 60; % minutes of sleep per fly
    y = (y ./ roi_min) * 60; % mins sleep per hour of the experiment

    y_avg = median(y,'omitnan');
    x = shuffle_data(linspace(ii-buffer,ii+buffer,num.trial(i)));
    % save data into structure for loading later
    total_sleep(ii).x = x; total_sleep(ii).y = y; total_sleep(ii).y_avg = y_avg;
    total_sleep(ii).protocol = data(i).temp_protocol; total_sleep(ii).buff = buffer;
    total_sleep(ii).name = data(i).ExpGroup; total_sleep(ii).color = grouped(i).color;
    
    % plot data
    k = grouped(i).color;
    scatter(x,y,SZ,grouped(i).color,'filled')
    % plot([ii-buffer,ii+buffer],[y_avg,y_avg],'color',k,'linewidth',LW)

    boxchart(ii*ones([length(y),1]), y',"BoxFaceColor",k,"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')

end

% formats and labels
ylabel('minutes of sleep per hour')
formatFig(fig,blkbgd);
set(gca,'xcolor','none','XTick',[],'XTickLabel',[])
set(gca,'fontsize', 20,'tickDir','in')

% save point data:
save([saveDir expGroup ' total sleep.mat'],'total_sleep');

save_figure(fig,[saveDir 'Sleep\' expGroup ' total sleep'],'-png',true,false);
save_figure(fig,[saveDir 'Sleep\' expGroup ' total sleep'],'-pdf');

%% FIGURE: Null distribution of fly distances to food:
% TODO : 
% 1) normalize the histograms to compare between the null and experiment
% 2) determine the median distance from the food (if flies were normally
% distributed across the arena


clearvars('-except',initial_vars{:})

% [foreColor,backColor] = formattingColors(blkbgd);
LW = 3;
if ~exist('null_dist','var')
    pathy = getDataPath(5);
    load([pathy(1:end-5) 'Fundamentals\Distance_to_well_distribution']);
    initial_vars{end+1} = 'null_dist';
end

% PLOT NULL DISTRIBUTION:
fig_dir = [saveDir 'Sleep\' 'sleeping histograms overlay\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
binedges = 0:3:51;
% ylimits = [-50 500]; % [-50,1400] 15-23 WT
% ylimits = [-50 800];

for i = 1:num.exp 
    fig = getfig('',1); hold on
    plotData = [];
    for trial = 1:num.trial(i)
        loc = sleep(i).trial(trial).sleepLoc;
        test = sleep(i).trial(trial).all_distance(loc);
        plotData = [plotData; test];
    end
    % Exp Data
    h(i) = histogram(plotData,binedges,'FaceColor',grouped(i).color,'FaceAlpha',0.7);

    % NULL DISTRIBUTION: 
    % null_h = histogram(null_dist,binedges,'FaceColor',Color('grey'),'FaceAlpha',0.7);

    % histogram(null_dist.distance,null_dist.binedges,'FaceColor',Color('grey'),'FaceAlpha',1)
    yyaxis right
    plot(null_dist.binedges,null_dist.prob,'color', Color('teal'),'linewidth',LW,'linestyle','-')

    formatFig(fig,blkbgd);
    xlabel('distance to well (mm)','FontSize',20)
    ylabel('sleeping flies (#)','FontSize',20)
    set(gca,'TickDir','out')
    xlim([-1.5000   51.5000])
%     ylim(ylimits)
    
    % save figure
    save_figure(fig,[fig_dir expNames{i} ' sleeping distance null histogram'],fig_type);
end

% save_figure(fig,[fig_dir 'Null sleeping distance distribution'],fig_type);

%% FIGURE: Cumulative distribution of sleeping distances to food overlaid with null distribution
clearvars('-except',initial_vars{:})
% [foreColor,backColor] = formattingColors(blkbgd);
LW = 2;

% PLOT NULL DISTRIBUTION:
fig_dir = [saveDir 'Sleep\' 'sleeping histograms overlay\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir) 
end 

% FIGURE:
fig = getfig('',1);  
hold on

% NULL DISTRIBUTION:
null_CDF = cdfplot(null_dist.distance);
null_CDF.Color = Color('red');
null_CDF.LineWidth = LW; 

for i = 1:num.exp 
    % Exp Data
    plotData = [];
    for trial = 1:num.trial(i)
        loc = sleep(i).trial(trial).sleepLoc;
        test = sleep(i).trial(trial).all_distance(loc);
        plotData = [plotData; test];
    end
    data_CDF = cdfplot(plotData);
    data_CDF.Color = grouped(i).color;
    data_CDF.LineWidth = LW; 
end
 
formatFig(fig,blkbgd);
xlabel('distance to well (mm)','FontSize',20)
ylabel('Empirical Cumulative Distribution','FontSize',20)
set(gca,'TickDir','out','XGrid', 'off','YGrid', 'off')


% save figure
save_figure(fig,[fig_dir  'Sleeping distance CDF with null distribution'],fig_type);


% clearvars('-except',initial_vars{:})
% LW = 2;
% [~,backColor] = formattingColors(blkbgd);
% 
% % FIGURE:
% fig = getfig('',1,[687 680]);  
% hold on
% % NULL DISTRIBUTION:
% null_CDF = cdfplot(null_dist.distance);
% null_CDF.Color = Color('teal');
% null_CDF.LineWidth = 3; 
% 
% formatFig(fig,blkbgd);
% xlabel('distance to well (mm)','FontSize',20)
% ylabel('Empirical Cumulative Distribution','FontSize',20)
% set(gca,'TickDir','out')
% set(gca,'GridColor',backColor)
% title('')
% 



%% Find probability of occupancy for each distance

clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);

% PLOT NULL DISTRIBUTION:
fig_dir = [saveDir 'Sleep\' 'sleeping histograms overlay\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end 

binSpacing = 1; %mm spacing for bins
binedges = 0:binSpacing:55; 
sSpan = 5;
LW = 2;

null_distance = null_dist.distance;

% find probability of occupancy for each bin from the null distribution
N = discretize(null_distance,binedges);
N_tot = length(null_distance);
prob = zeros(length(binedges),1);
for i = 1:length(binedges)
    prob(i) = (sum(N==i))/N_tot;
end

fig = getfig('',1); hold on
yyaxis left
histogram(null_distance,binedges,'FaceColor',Color('grey'),'FaceAlpha',1)
ylabel('count')
yyaxis right
plot(binedges,prob,'color', foreColor,'linewidth',LW,'linestyle','-')
plot(binedges,smooth(prob,sSpan,'moving'),'color', Color('cyan'),'linewidth',LW+2,'linestyle','-')
formatFig(fig,blkbgd);
xlabel('distance to well (mm)')
ylabel('occupancy probability')
yyaxis left
set(gca,'Ycolor',foreColor)

% Save null distribution probability
save_figure(fig,[fig_dir  'Distance probabilty null distribution'],fig_type);


 
%% FIGURE: Null distribution vs experimental fly distances to food:
% TODO : 
% 1) normalize the histograms to compare between the null and experiment
% 2) determine the median distance from the food (if flies were normally
% distributed across the arena


clearvars('-except',initial_vars{:})

% [foreColor,backColor] = formattingColors(blkbgd);
LW = 3;
if ~exist('null_dist','var')
    load([getCloudPath 'Fundamentals\Distance_to_well_distribution']);
    initial_vars{end+1} = 'null_dist';
end

% PLOT NULL DISTRIBUTION:
fig_dir = [saveDir 'Sleep\' 'sleeping histograms overlay\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

binlow = 0;
binhigh = 51;
binwidth = 3;
binedges = binlow:binwidth:binhigh;
bincenters = (binwidth-binlow)/2:binwidth:binhigh;

% ylimits = [-50 500]; % [-50,1400] 15-23 WT
% ylimits = [-50 800];


for i = 1:num.exp 
    fig = getfig('',1); hold on
    plotData = [];
    for trial = 1:num.trial(i)
        loc = sleep(i).trial(trial).sleepLoc;
        test = sleep(i).trial(trial).all_distance(loc);
        plotData = [plotData; test];
    end

    % NULL DISTRIBUTION: 
    tic
    N = 1000; % number of selections from the null distribution to make the display null distr
    null_h_temp = nan([N, length(binedges)-1]);
    for jj = 1:N
        idx = randi(length(null_dist.distance),[1,length(plotData)]);
        null_h_temp(jj,:) = histcounts(null_dist.distance(idx), binedges);
    end
    null_h = mean(null_h_temp);
    toc
    bar(bincenters, null_h,'BarWidth', 1,'FaceColor',Color('grey'),'FaceAlpha',0.7);

    % Exp Data
    histogram(plotData,binedges,'FaceColor',grouped(i).color,'FaceAlpha',0.7);

    % % histogram(null_dist.distance,null_dist.binedges,'FaceColor',Color('grey'),'FaceAlpha',1)
    % yyaxis right
    % plot(null_dist.binedges,null_dist.prob,'color', Color('teal'),'linewidth',LW,'linestyle','-')

    formatFig(fig,blkbgd);
    xlabel('distance to well (mm)','FontSize',20)
    ylabel('sleeping flies (#)','FontSize',20)
    set(gca,'TickDir','out')
    xlim([-1.5000   51.5000])
%     ylim(ylimits)
    
    % save figure
    save_figure(fig,[fig_dir expNames{i} ' sleeping distance vs. normalized null histogram'],fig_type);
end

% save_figure(fig,[fig_dir 'Null sleeping distance distribution'],fig_type);

%% FIGURE: Probability of sleeping for each distance to food well
% 1) normalize the histograms to compare between the null and experiment
% 2) determine the median distance from the food (if flies were normally
% distributed across the arena
clearvars('-except',initial_vars{:})

[foreColor,backColor] = formattingColors(blkbgd);
LW = 3;

if ~exist('null_dist','var')
    load([getCloudPath 'Fundamentals\Distance_to_well_distribution']);
    initial_vars{end+1} = 'null_dist';
end

binlow = 0;
binhigh = 51;
binwidth = 3;
binedges = binlow:binwidth:binhigh;
bincenters = (binwidth-binlow)/2:binwidth:binhigh;

% ylimits = [-50 500]; % [-50,1400] 15-23 WT
% ylimits = [-50 800];

stats = [];

idx = [4 1 2 3 6 5];
icount = 0;
fig = getfig('',1); hold on

for i = idx
    icount = icount+1;
    % Pull the data
    plotData = [];
    for trial = 1:num.trial(i)
        loc = sleep(i).trial(trial).sleepLoc;
        test = sleep(i).trial(trial).all_distance(loc);
        plotData = [plotData; test];
    end
    stats(icount).name = grouped(i).name;
    stats(icount).data = plotData;
    % sort number of sleeping flies by distance
    h = discretize(plotData, binedges);
    bCount = [];
    for bin = 1:length(binedges)-1
        bCount(bin) = sum(h==bin);
    end
    bProb = (bCount/length(plotData))*100;

    % plot the percent occupancy
   switch i 
        case 4 % intact & dynamic caviar
            kolor = Color('dodgerblue');
        case 1 % intact & dynamic caviar
            kolor = Color('black');
        case 2 % intact & dynamic caviar
            kolor = Color('darkslategray');
       case 3 % intact & dynamic caviar
            kolor = Color('dimgray');
        case 6 % intact & dynamic caviar
            kolor = Color('silver');
        case 5 % intact & dynamic caviar
            kolor = Color('gold');
    end
    plot(bincenters, bProb, 'color', kolor, 'linewidth', LW)
end

% NULL DISTRIBUTION (smooth) 
h = discretize(null_dist.distance, binedges);
nCount = [];
for bin = 1:length(binedges)-1
    nCount(bin) = sum(h==bin);
end
nProb = (nCount/length(null_dist.distance))*100;
plot(bincenters, nProb, 'color', foreColor, 'linewidth', LW,'linestyle',':')

formatFig(fig,blkbgd);
xlabel('distance to well (mm)','FontSize',20)
ylabel('sleeping flies (%)','FontSize',20)
set(gca,'TickDir','out')
xlim([-1.5000   51.5000])
    
% save figure
save_figure(fig,[saveDir  'Sleeping flies at each distance from food'],fig_type);

% STATS: Ask if the distribution of distances are different from each other?
alpha = 0.05; % Original significance level
adjustedAlpha = alpha / (length(stats)); % Adjusted significance level for num comparisons
p = []; h = [];
for i = 1:length(stats)
    idx = randi(length(null_dist.distance),[1,length(stats(i).data)]);
    [p(i), ~, stats(i).results] = ranksum(stats(i).data, null_dist.distance(idx));
    h(i) = p(i)<adjustedAlpha;  % Bonferroni correction
    % display results:
    fprintf('Wilcoxon rank-sum test (%s)): h = %d, p = %f\n', stats(i).name, h(i), p(i));
end



 
%% [Temp-Shift ONLY] How does 'post-prandial sleep' change compared to post-temperature sleep? 
% Compare the 'quantity' of sleep between trials with food and without food
clearvars('-except',initial_vars{:})
clc
% compare the sleep quantity during the period preceding the ramp & the
% down ramp vs ramp up vs post ramp up
[foreColor,backColor] = formattingColors(blkbgd);

timeBuff = 10; %minutes of inclusion before the ramp
plotData = [];

for exp = 1:num.exp
    % pull the time ROIs
    tp = getTempTurnPoints(data(exp).temp_protocol);
    buff = tp.fps*timeBuff*60; % how many sampling points within the buffer time region
    % [pre_ROI, post_ROI] = deal(nan(tp.nDown,2));
    [pre_ROI, post_ROI] = deal([]);
    for i = 1:tp.nDown % for each ramp pull the values
        pre_ROI = [pre_ROI, tp.down(i,1)-buff : tp.down(i,2)];
        post_ROI = [post_ROI, tp.up(i,1) : tp.up(i,2) + buff];
    end
    plotData(exp).preROI = pre_ROI;
    plotData(exp).postROI = post_ROI;
    % find the fraction of flies that are sleeping for each trial during
    % this period of time...
    plotData(exp).preAll = mean(sleep(exp).fract_sleep(pre_ROI,:),'omitnan')*100;
    plotData(exp).postAll = mean(sleep(exp).fract_sleep(post_ROI,:),'omitnan')*100;
    plotData(exp).preAvg = median(plotData(exp).preAll);
    plotData(exp).postAvg = median(plotData(exp).postAll);
end

% Plot Figure: 
r = 1;
c = 2;
sz = 100;
LW = 3;

spacer = 3;
x_loc = 1:num.exp*spacer;
x_loc(spacer:spacer:end) = [];
exp_sel = [2 6 3 4 5 1];
b = 0.5;

fig = getfig('',1);
for i = 1:num.exp
    exp = (exp_sel(i));
    kolor = grouped(exp).color;
    x = x_loc(i)*ones(1,num.trial(exp));

    % pre-sleep measure
    subplot(r,c,1); hold on
    y = plotData(exp).preAll;
    y_err = std(y);
    y_avg = plotData(exp).preAvg;
    scatter(x,y,floor(sz/4),kolor,'xjitter', 'density')
    % plot([x_loc(i)-b,x_loc(i)+b],[y_avg, y_avg],'color', foreColor,'linewidth',LW)
    errorbar(x_loc(i),y_avg, y_err, 'color', foreColor,'linewidth', LW)
    scatter(x_loc(i),y_avg,sz,foreColor, 'filled')

     % post-sleep measure
    subplot(r,c,2); hold on
    y = plotData(exp).postAll;
    y_err = std(y);
    y_avg = plotData(exp).postAvg;
    scatter(x,y,floor(sz/4),kolor,'xjitter', 'density')
    % plot([x_loc(i)-b,x_loc(i)+b],[plotData(exp).postAvg, plotData(exp).postAvg],'color', foreColor,'linewidth',LW)
    errorbar(x_loc(i),y_avg, y_err, 'color', foreColor,'linewidth', LW)
    scatter(x_loc(i),y_avg,sz,foreColor, 'filled')

end
% format figure
formatFig(fig,blkbgd,[r,c]); 
subplot(r,c,1)
set(gca, 'xcolor', 'none')
ylabel('Avg % flies sleeping pre cooling and during cooling')
ylim([0,50])
subplot(r,c,2)
set(gca, 'xcolor', 'none')
ylabel('Avg % flies sleeping during warming and post warming')
ylim([0,50])

% Run quick stats on the comparisons
[h_pre,p_pre,h_post,p_post] = deal([]);
idx = 1;
for i = 1:3
    food = exp_sel(idx);
    nofood = exp_sel(idx+1);
    idx = idx+2;

    % pre
    y = plotData(food).preAll;
    z = plotData(nofood).preAll;
    [h_pre(i),p_pre(i)] = ttest2(y,z);
    disp([grouped(food).name ' vs' grouped(nofood).name ' pre:'])
    disp([num2str(p_pre(i)) '  ' num2str(h_pre(i))])

    % post
    y = plotData(food).postAll;
    z = plotData(nofood).postAll;
    [h_post(i),p_post(i)] = ttest2(y,z);
    disp([grouped(food).name ' vs' grouped(nofood).name ' post:'])
    disp([num2str(p_post(i)) '  ' num2str(h_post(i))])
end


% Add stats to figure
figure(fig)
y = rangeLine(fig,4,true);
y2 = rangeLine(fig,1.5,true);
idx = 1;
for i = 1:3
   x1 = idx-b;
   x2 = idx+1+b;
   idx = idx+spacer;
   subplot(r,c,1) % pre
   if h_pre(i)
       plot([x1,x2],[y,y],'color', foreColor,'linewidth', LW)
       scatter(mean([x1,x2]),y2,80,foreColor,"filled","*","MarkerEdgeColor",foreColor)
   end
   subplot(r,c,2) % post
   if h_post(i)
       plot([x1,x2],[y,y],'color', foreColor,'linewidth', LW)
       scatter(mean([x1,x2]),y2,80,foreColor,"filled","*","MarkerEdgeColor",foreColor)
   end
end
 
% Save Figure: 
save_figure(fig,[saveDir 'Sleep\Post prandial sleep comparision'],fig_type);


% two way anova example

% y = [52.7 57.5 45.9 44.5 53.0 57.0 45.9 44.0]';
% g1 = [1 2 1 2 1 2 1 2]; 
% g2 = {'hi';'hi';'lo';'lo';'hi';'hi';'lo';'lo'}; 
% g3 = {'may';'may';'may';'may';'june';'june';'june';'june'};
% 
% 
% p = anovan(y,{g1 g2 g3},'model','interaction','varnames',{'g1','g2','g3'})
% 
% [~,~,stats] = anovan(y,{g1 g2 g3},"Model","interaction", ...
%     "Varnames",["g1","g2","g3"]);


% ANOVA Stats
% Run quick stats on the comparisons
[h_pre,p_pre,h_post,p_post] = deal([]);
idx = 1;
for i = 1:3
    food = exp_sel(idx);
    nofood = exp_sel(idx+1);
    idx = idx+2;

    % pre
    y = plotData(food).preAll;
    z = plotData(nofood).preAll;
    [h_pre(i),p_pre(i)] = ttest2(y,z);
    disp([grouped(food).name ' vs' grouped(nofood).name ' pre:'])
    disp([num2str(p_pre(i)) '  ' num2str(h_pre(i))])

    % post
    y = plotData(food).postAll;
    z = plotData(nofood).postAll;
    [h_post(i),p_post(i)] = ttest2(y,z);
    disp([grouped(food).name ' vs' grouped(nofood).name ' post:'])
    disp([num2str(p_post(i)) '  ' num2str(h_post(i))])
end

% post-warming data anova
temp_proto_list = [1,3,2,2,1,3]; 
food_list = [0,1,1,0,1,0];
[y,g1,g2] = deal([]);
for i = 1:num.exp
    if i == 3 || i == 4
        continue
    end
    y = [y, plotData(i).postAll];
    g1 = [g1, temp_proto_list(i)*ones([1,num.trial(i)])];
    g2 = [g2, food_list(i)*ones([1,num.trial(i)])];
end

p = anovan(y,{g1 g2},'model','interaction','varnames',{'temp','food'});

% [p,~,stats] = anovan(y,{g1 g2},"Model","interaction", "Varnames",['temp';'food']);

% post-warming data anova
[y,g1,g2] = deal([]);
for i = 1:num.exp
    y = [y, plotData(i).preAll];
    g1 = [g1, temp_proto_list(i)*ones([1,num.trial(i)])];
    g2 = [g2, food_list(i)*ones([1,num.trial(i)])];
end

p = anovan(y,{g1 g2},'model','interaction','varnames',{'temp','food'});





cList = {'dodgerblue', 'grey'};
cList = repmat(cList,[1,ceil(num.exp/2)]);
fig = getfig('',1);

for i = 1:num.exp
    exp = (exp_sel(i));
    kolor = grouped(exp).color;
    x = x_loc(i)*ones(1,num.trial(exp));

    % pre-sleep measure
    subplot(r,c,1); hold on
    y = plotData(exp).preAll;
    y_err = std(y);
    y_avg = plotData(exp).preAvg;
      scatter(x,y,floor(sz/4),Color(cList{i}),'filled','xjitter', 'density')
    boxchart(x, y,"BoxFaceColor",Color(cList{i}),"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor,'MarkerColor',Color(cList{i}),...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
  
   
     % post-sleep measure
    subplot(r,c,2); hold on
    y = plotData(exp).postAll;
    y_err = std(y);
    y_avg = plotData(exp).postAvg;
    scatter(x,y,floor(sz/4),Color(cList{i}),'filled','xjitter', 'density')
    boxchart(x, y,"BoxFaceColor",Color(cList{i}),"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor, 'MarkerColor',Color(cList{i}),...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    

end

x_max = 30;
% format figure
formatFig(fig,blkbgd,[r,c]); 
subplot(r,c,1)
title('Cooling')
set(gca, 'xcolor', 'none')
ylabel({'Sleeping flies (%)'; '   '})
ylim([0,x_max])
set(gca,'ytick', 0:10:x_max,'FontSize', 20,'TickDir','in')
subplot(r,c,2)
title('Warming')
set(gca, 'xcolor', 'none')
ylabel({'Sleeping flies (%)'; '   '})
ylim([0,x_max])
set(gca,'ytick', 0:10:x_max,'FontSize', 20,'TickDir','in')

% Save Figure: 
save_figure(fig,[saveDir 'Sleep\Post prandial sleep comparision box plots'],'-png');

%% FIGURE: percent of flies asleep in different spatial regions
clearvars('-except',initial_vars{:})
clc
[foreColor,backColor] = formattingColors(blkbgd);
sz = 50;
r = 1; c = 2;
buff = 0.2;
LW = 2;

% expList = 1:num.exp;
expList = expOrder;

fig = getfig('',1,[1064 837]);
groupNames = []; h = []; p = [];
for idx = 1:length(expList)
    i = expList(idx);
    k = grouped(i).color;
    tp = getTempTurnPoints(grouped(i).fictivetemp);
    trange = [tp.hold(:); tp.down(:); tp.up(:)];
    x_roi = min(trange):max(trange);
    subplot(r,c,1) % quadrant
    hold on
    y = mean(grouped(i).sleep.quad.all(x_roi,:),'omitnan');
    x = idx*ones(length(y),1); 
    boxchart(x, y',"BoxFaceColor",k,"BoxFaceAlpha",1,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    scatter(x,y,sz, k, 'filled', 'xjitter', 'density')

    % check if the distribution is significantly different from the null distribution
    [h(idx),p(idx)] = ttest(y,25);
    

    subplot(r,c,2) % ring
    hold on
    y = mean(grouped(i).sleep.ring.percent(x_roi,:),'omitnan');

     boxchart(x, y',"BoxFaceColor",k,"BoxFaceAlpha",1,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
     scatter(i*ones(length(y),1),y,sz, k, 'filled', 'xjitter', 'density')
    groupNames{end+1} = grouped(i).name;
end

formatFig(fig, blkbgd,[r,c]);
title_str = {'quadrant', 'ring'};
for ii = 1:2
    subplot(r,c,ii)
    set(gca, 'xtick', 1:num.exp,'xticklabel',groupNames,'XTickLabelRotation',45)
    ylim([0,100])
    set(gca, 'ytick',0:20:100)
    ylabel('flies sleeping in region (%)')
    title(title_str{ii},'color', foreColor)
    h_line(25,'grey', '--',1.5)
end

save_figure(fig,[saveDir 'sleeping flies in quad and ring'],'-png',true,false);
save_figure(fig,[saveDir 'sleeping flies in quad and ring'],'-pdf',true,true);

%%


% check significantly different from the null distribution. 



















    