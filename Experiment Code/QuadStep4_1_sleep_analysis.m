
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
% How long does a fly need to be still to count as 'sleep'
min_duration = 5*3*60; % 5 mins * 3fps*60sec = data point number that must be met 'unmoving'
nbins = 50;

% Create sleep data for unprocessed files (trial by trial)
for i = 1:num.exp
    T = data(i).T;
    for trial = 1:num.trial(i)
        sleep_file = [baseFolder T.Date{trial} '\Arena ' T.Arena{trial} '\' T.ExperimentID{trial} ' sleeping data.mat'];  
        if ~exist(sleep_file,"file")
            sleepData = struct;
            %preallocate for speed and space
            trial_length = length(data(i).data(trial).occupancy.time);
            [N,frameCount,sleepLoc] = deal(nan(nbins,nbins,trial_length));
            sleepingCount = zeros(trial_length,1);

            % Set axis limits for the selected arena
            x = data(i).data(trial).data.centre(1);
            y = data(i).data(trial).data.centre(2);
            r = data(i).data(trial).data.r;
            xlimit = [x-(r+50),x+(r+50)];
            ylimit = [y-(r+50),y+50+r];

            % find the 'auto bin' lines
            xedge = linspace(xlimit(1),xlimit(2),nbins+1);
            yedge = linspace(ylimit(1),ylimit(2),nbins+1);

            % pull the fly locations during the trial
            x_loc = data(i).data(trial).data.x_loc;
            y_loc = data(i).data(trial).data.y_loc;
            for frame = 1:trial_length
                X = x_loc(frame,:); X(isnan(X)) = [];
                Y = y_loc(frame,:); Y(isnan(Y)) = [];
                N(:,:,frame) = histcounts2(X,Y,xedge,yedge);
            end

            % find grid space that have continuous occupation for more than min_duration
            frameCount(:,:,1) = N(:,:,1);
            for frame = 2:trial_length
                currFrame = N(:,:,frame); % current frame locations
                resetLoc = currFrame==0; % locations that do not have flies and thus need a count reset

                tempCount = frameCount(:,:,frame-1)+currFrame; % add current frames to list
                tempCount(resetLoc) = 0; % reset counts for spots with no flies

                frameCount(:,:,frame) = tempCount; % add current count into the saving structure
            end

            % ---- Vectorize the data (find the flies that are sleeping....) -----

            % find food well location for distance capture later
            foodWellLoc = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));
            c1 = foodWellLoc(1);
            c2 = foodWellLoc(2);
            
            % create empty matrixes for the x and y positions of the sleeping flies
            sleeping = struct;
            [sleeping.X, sleeping.Y, sleeping.all_distance] = deal(nan(trial_length,data(i).T.NumFlies(trial)));
            sleeping.sleepNum = zeros(trial_length,1);
            [sleeping.dist_avg, sleeping.dist_err] = deal(nan(trial_length,1));
            % assign data by frame
            for frame = 1:trial_length
                frame_data = frameCount(:,:,frame) > min_duration;
                binLoc = find(frame_data>0);
                
                % Find the coordinates of the sleeping flies bins from the discretized data
                y_row = ceil(binLoc/nbins);
                x_row = rem(binLoc-1,nbins)+1;
                x_position = (xedge(x_row) + xedge(x_row+1))/2;
                y_position = (yedge(y_row) + yedge(y_row+1))/2;
                
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

sleep = struct;
for i = 1:num.exp
    T = data(i).T;
    for trial = 1:num.trial(i)
        sleep_file = [baseFolder T.Date{trial} '\Arena ' T.Arena{trial} '\' T.ExperimentID{trial} ' sleeping data.mat'];  
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
fps = 3;
% Thermal threat quantification and avg sleep quantity
for i = 1:num.exp
    tPoints = getTempTurnPoints(data(i).T.TempProtocol{1});
    demoRamp_idx = tPoints.down(1,1):tPoints.up(1,2);
    
    % Thermal threat
    temp_ramp = grouped(i).temp(demoRamp_idx);
    thermalThreat = (sum(25-temp_ramp)/fps)/1000;
    
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

clearvars('-except',initial_vars{:})
disp('All finished')

%% FIGURE: Sleeping over time
plot_err = false;
[~,backColor] = formattingColors(blkbgd); %get background colors
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
    subplot(r,c,sb(2).idx); hold on
        y = smooth(sleep(i).sleepfract_avg,sSpan,'moving');
        y_err = smooth(sleep(i).sleepfract_err,sSpan,'moving');
        if plot_err
            fill_data = error_fill(time, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h, 'facealpha', 0.4)
        end
        plot(time,y,'color',kolor,'linewidth',LW)
        ylabel('fraction of flies sleeping')
        xlabel('time (min)')
end

formatFig(fig,true,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca,'xcolor',backColor) 

save_figure(fig,[saveDir expGroup ' sleep timecourse'],fig_type);

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
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h, 'facealpha', 0.4)
        end
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
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.35)
        end
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

save_figure(fig,[saveDir expGroup ' sleeping flies summary'],fig_type);

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
    y = sleep(i).temp_all(:,1);
    y_err = sleep(i).temp_all(:,2);
    loc = isnan(y)|isnan(y_err);
    x(loc) = [];
    y(loc) = [];
    y_err(loc) = [];

    plot(x,y,'color',kolor,'linewidth',LW+1)
    if plot_err
        fill_data = error_fill(x, y, y_err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
        set(h, 'facealpha', 0.35)
    end
end
xlabel('Temperature (\circC)')
ylabel('fraction of flies sleeping')
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
        y = sleep(i).(section_type)(:,1);
        y_err = sleep(i).(section_type)(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
%         if plot_err
%             fill_data = error_fill(x, y, y_err);
%             h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
%             set(h, 'facealpha', 0.2)
%         end
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',l_style)
    end
end
xlabel('Temperature (\circC)')
ylabel('fraction of flies sleeping')
    
% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
if equalLim
    fig = matchAxis(fig,true);
end
% ylim(num_temp_lim)

% save figure
save_figure(fig,[saveDir 'Flies sleeping during heating and cooling'],fig_type);

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
    
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
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
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
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

% save figure
save_figure(fig,[saveDir 'sleep and distance during heating and cooling'],fig_type);

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
        
            % COOLING
            yyaxis left
            subplot(r,c,type);   hold on
            %Distance Data 
            x = grouped(i).(section_type).temps;
            y = grouped(i).(section_type).avg;
            y_err = grouped(i).(section_type).err;
            loc = isnan(y) | isnan(y_err);% remove nans 
            y(loc) = []; x(loc) = []; y_err(loc) = [];
        
            if plot_err
                fill_data = error_fill(x, y, y_err);
                h = fill(fill_data.X, fill_data.Y, distColor, 'EdgeColor','none','HandleVisibility','off');
                set(h, 'facealpha', 0.2)
            end
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
            if plot_err
                fill_data = error_fill(x, y, y_err);
                h = fill(fill_data.X, fill_data.Y, sleepColor, 'EdgeColor','none','HandleVisibility','off');
                set(h, 'facealpha', 0.2)
            end
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
                ylim([min(yleft(1,:)), max(yleft(2,:))])
                yyaxis right
                ylim([min(yright(1,:)), max(yright(2,:))])
            end
        end
        
    % save figure
    save_figure(fig,[saveDir expNames{i} ' sleep and distance during heating and cooling'],fig_type);

end


%% FIGURE: Sleep vs thermal threat
% Can we easily predict behavior based on a simple cold-exposure metric?
clearvars('-except',initial_vars{:})
% quantify the cold exposure during a single ramp as sum of difference in
% temp from preferred temp: (might need to make it exponential????)
fps = 3;

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
formatFig(fig,true);

save_figure(fig,[saveDir 'Sleep duration by thermal stress'],fig_type);

%% FIGURE: Avg quantity of sleep per fly

% Can we easily predict behavior based on a simple cold-exposure metric?
clearvars('-except',initial_vars{:})
% quantify the cold exposure during a single ramp as sum of difference in
% temp from preferred temp: (might need to make it exponential????)
fps = 3;

sleepDuration = [];
for i = 1:num.exp
    % Thermal threat
    tPoints = getTempTurnPoints(data(i).T.TempProtocol{1});
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
formatFig(fig,true);

save_figure(fig,[saveDir 'Sleep duration by thermal stress rate'],fig_type);



%% Sleep location? Where do flies choose to sleep?
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
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
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
formatFig(fig,true,[r,c]);
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
save_figure(fig,[saveDir 'sleeping distance to food'],fig_type);


%% FIGURE: Histogram of sleep proximity to food
clearvars('-except',initial_vars{:})
LW = 2;
r = 1; %rows
c = 2; %columns
plot_err = true;
[foreColor,~] = formattingColors(blkbgd);
equalLim = true;
% 
% TODO: make these each their own plots on the same figure? (no overlaid?)
% or some min num is overlaid only?
fig = getfig('',1); hold on
for i = 1:num.exp
    plotData = [];
    for trial = 1:num.trial(i)
        loc = sleep(i).trial(trial).sleepLoc;
        test = sleep(i).trial(trial).all_distance(loc);
        plotData = [plotData; test];
    end
    h(i) = histogram(plotData,'FaceColor',grouped(i).color,'FaceAlpha',0.7);
end
formatFig(fig,true);
xlabel('distance to well (mm)','FontSize',20)
ylabel('sleeping flies (#)','FontSize',20)
set(gca,'TickDir','out')

% save figure
save_figure(fig,[saveDir 'sleeping distance to food histogram'],fig_type);

%% FIGURE: fly sleeping positions in the arena
clearvars('-except',initial_vars{:})
% HAVENT SCREENED FOR ONLY RAMPSSSSSSSS!!! 
% TODO : 
% change the orientation/rotation to match across trials
% color code the dot to indicate the temperature for each sleeping fly
% location
[foreColor,backColor] = formattingColors(blkbgd);
SZ = 40;

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


% Determine the color choices
temp_limits = [nan nan];
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    temp_limits = [min([tp.threshLow, temp_limits(1)]), max([tp.threshHigh, temp_limits(2)])];
end
nColors = floor(temp_limits(1)):1:ceil(temp_limits(2));
cMap = turbo(length(nColors));

% TODO HERE...
    
% Plot scatter point figure:
for i = 1:num.exp
    all_data = sleep(i).location.all_data;
    kolor = grouped(i).color;
    
    % fancy color
    temp = sleep(i).location.all_data(:,3);
    color_loc = discretize(temp,nColors);
    
    all_data(:,5) = color_loc;

    find(isnan(color_loc))

    fig = getfig('',1); hold on
  
    % PLOT
    scatter(sleep(i).location.wells(1:4,1),sleep(i).location.wells(1:4,2),SZ+20,Color('grey'),'filled')
    scatter(0,0,SZ+20,'green','filled')
    scatter(sleep(i).location.all_data(:,1),sleep(i).location.all_data(:,2),15,kolor,'filled')
   
    viscircles([sleep(i).location.arenaCenter(1),sleep(i).location.arenaCenter(2)],sleep(i).location.r,'Color',grouped(i).color);
    axis square; axis equal;  
    formatFig(fig,blkbgd);
    set(gca,'XColor',backColor,'YColor',backColor);
    save_figure(fig,[saveDir grouped(i).name ' sleeping positions in arena'],fig_type);

end


% Plot heatmap figure?
nBins = 15;

for i = 1:num.exp

%     kolor = grouped(i).color;
    
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
        formatFig(fig,true);
        c.Color = foreColor;
        set(gca,'xcolor',backColor,'ycolor',backColor)
        c.Label.String = '# Flies';
        c.Label.Color = foreColor;
        c.Label.FontSize = 20;
        save_figure(fig,[saveDir grouped(i).name ' sleeping positions in arena heatmap'],fig_type);

end





%%

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
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
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
formatFig(fig,true,[r,c]);
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
save_figure(fig,[saveDir 'sleeping distance to food'],fig_type);







%% Sleep duration and statistics across groups
[foreColor,backColor] = formattingColors(blkbgd);

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
formatFig(fig, true);

ylim([-2,12])

save_figure(fig,[saveDir 'Sleep duration across groups violin plot'],fig_type);






 













%% FIGURE: Sleep vs thermal threat by genotype
% Can we easily predict behavior based on a simple cold-exposure metric?
clearvars('-except',initial_vars{:})
% quantify the cold exposure during a single ramp as sum of difference in
% temp from preferred temp: (might need to make it exponential????)

% expList = {'Berlin WT', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
% colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
kolor = Color('White');

fps = 3;

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
%     kolor = grouped(i).color;
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
formatFig(fig,true);

save_figure(fig,[saveDir 'Sleep duration by thermal stress norm axes'],fig_type);

% save_figure(fig,[saveDir 'Sleep duration by thermal stress norm axes short'],fig_type);

% 















