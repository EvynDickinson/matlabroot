
% Find 'sleep' points in data

clearvars('-except',initial_vars{:})
% How long does a fly need to be still to count as 'sleep'
min_duration = 5*3*60; % 5 mins * 3fps*60sec = data point number that must be met 'unmoving'
nbins = 50;


%% visualize grid spacing
nbins = 50;
i = 1;
trial = 1;
vid = 1;
frame = 1;

% pull info for the first trial:
dataDate = data(i).T.Date{trial};
vid_name = data(i).T.ExperimentID{trial};
vidDir = [baseFolder dataDate '/' vid_name '_'];
videoPath = [vidDir num2str(vid) '.avi'];
movieInfo = VideoReader(videoPath); %read in video

% Set axis limits for the selected arena
x = data(i).data(trial).data.centre(1);
y = data(i).data(trial).data.centre(2);
r = data(i).data(trial).data.r;
xlimit = [x-(r+50),x+(r+50)];
ylimit = [y-(r+50),y+50+r];

% Plot image of video frame
fig = figure; set(fig,'pos',[-1030 279 772 1009],'color','k');
currentImg = rgb2gray(read(movieInfo,frame));
imshow(currentImg)
xlim(xlimit); ylim(ylimit);   

% find the 'auto bin' lines
xedge = linspace(xlimit(1),xlimit(2),nbins+1);
yedge = linspace(ylimit(1),ylimit(2),nbins+1);

% Plot the bin outline edges:
h_line(yedge,'yellow','-',0.25) 
v_line(xedge,'yellow','-',0.25)




%% Single trial quant and figure

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


%% ANALYSIS: run and save the sleep quanitification

for i = 1:num.exp
    sleep_file = [baseFolder,'Data structures\',expNames{i},...
               '\',expNames{i},' sleeping.mat'];
    if exist(sleep_file,"file")
        load(sleep_file,'sleepData')
        grouped(i).sleep.all = sleepData;
    else
        sleepData = struct;
        for trial = 1:num.trial(i)
        
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
            trial_length = size(y_loc,1);
            
            N = [];
            for frame = 1:trial_length
                X = x_loc(frame,:); X(isnan(X)) = [];
                Y = y_loc(frame,:); Y(isnan(Y)) = [];
                N(:,:,frame) = histcounts2(X,Y,xedge,yedge);
            end
            
            % find grid space that have continuous occupation for more than min_duration 
            frameCount = N(:,:,1);
            for frame = 2:trial_length
                currFrame = N(:,:,frame); % current frame locations
                resetLoc = currFrame==0; % locations that do not have flies and thus need a count reset
                
                tempCount = frameCount(:,:,frame-1)+currFrame; % add current frames to list
                tempCount(resetLoc) = 0; % reset counts for spots with no flies
                
                frameCount(:,:,frame) = tempCount; % add current count into the saving structure
            end
            
            % Find instances of 'sleep'?
            [sleepingCount,sleepLoc] = deal([]);
            for frame = 1:trial_length
                sleepLoc(:,:,frame) = frameCount(:,:,frame) > min_duration;
                sleepingCount(frame) = sum(sum(sleepLoc(:,:,frame)));
            end
        
            % Save into group structures
            sleepData(trial).frameCount = frameCount;
            sleepData(trial).sleepLoc = sleepLoc;
            sleepData(trial).sleepingCount = sleepingCount;
            fprintf(['\nDone exp ' num2str(i) ' trial ' num2str(trial)])
        end

        % Save data into external dataset
        grouped(i).sleep.all = sleepData;
        save(sleep_file,'sleepData','-v7.3'); clear sleepData
    end
    disp(['Done exp ' num2str(i)])
end
disp('All finished')


%% ANALYSIS: find average sleep values

for i = 1:num.exp
    all_sleep = [];
    for trial = 1:num.trial(i)
        sleepData = grouped(i).sleep.all;
        sleepCount = sleepData(trial).sleepingCount./data(i).T.NumFlies(trial);
        all_sleep = autoCat(all_sleep,sleepCount,true);
    end
    grouped(i).sleep.sleep_avg = mean(all_sleep,1,'omitnan');
    grouped(i).sleep.sleep_err = std(all_sleep,0,1,'omitnan');
    grouped(i).sleep.all_sleep = all_sleep;
end


%% FIGURE: sleeping over time
plot_err = true;
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
        y = smooth(grouped(i).sleep.sleep_avg,sSpan,'moving');
        y_err = smooth(grouped(i).sleep.sleep_err,sSpan,'moving');
        if plot_err
            fill_data = error_fill(time, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h, 'facealpha', 0.4)
        end
        plot(time,y,'color',kolor,'linewidth',LW)
        ylabel('fraction flies sleeping')
        xlabel('time (min)')
end

formatFig(fig,true,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca,'xcolor',backColor) 

save_figure(fig,[saveDir expGroup ' sleep timecourse'],fig_type);


%% ANALYSIS & FIGURE: sleeping tuning curve...
clearvars('-except',initial_vars{:})
plot_err = true;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors
LW = 1.5;
sSpan = 180; %1 minute smoothing

plotData = struct;

% Cluster the flies on food by temperature?
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
            % increasing rates:
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            plotData(i).(g_name)(temp,1) = mean(mean(grouped(i).sleep.all_sleep(:,loc),2,'omitnan'),'omitnan'); %avg 
            plotData(i).(g_name)(temp,2) = std(mean(grouped(i).sleep.all_sleep(:,loc),2,'omitnan'),'omitnan');%./num.trial(i); %err
        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        plotData(i).temp_all(temp,1) = mean(mean(grouped(i).sleep.all_sleep(:,loc),2,'omitnan'),'omitnan'); %avg 
        plotData(i).temp_all(temp,2) = std(mean(grouped(i).sleep.all_sleep(:,loc),2,'omitnan'),'omitnan')./num.trial(i);% %err
    end
    plotData(i).temps = temps;
end
disp('All finished')


% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive? 
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
         y = smooth(grouped(i).sleep.sleep_avg,sSpan,'moving');
        y_err = smooth(grouped(i).sleep.sleep_err,sSpan,'moving');
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
        x = plotData(i).temps;
        y = plotData(i).temp_all(:,1);
        y_err = plotData(i).temp_all(:,2);
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


% FIGURE: Flies on food by heating / cooling
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
    x = plotData(i).temps;
    y = plotData(i).temp_all(:,1);
    y_err = plotData(i).temp_all(:,2);
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
        x = plotData(i).temps;
        y = plotData(i).(section_type)(:,1);
        y_err = plotData(i).(section_type)(:,2);
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
















