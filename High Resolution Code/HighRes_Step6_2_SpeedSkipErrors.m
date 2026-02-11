
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
            subplot(r,c,ff)
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

%% TODO: determine method and approach for dealing with the swap frames

%% FIGURE CHECK: swap frames identified by both flies over the speed limit
clearvars('-except',initial_var{:})
saveDir = createFolder([figDir, 'extreme speed troubleshooting/']);

trial = 2;
sz = 10;
speed_thresh = 30; %mm/s speed threshold

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
allROI = confident_frames + (-3:3);

ax = [];
curr_frame = 0; % counter for switch location
for gg = 1:nFigs

    % PLOT DATA
    fig = getfig('',1);   
    set(fig, 'Visible', 'off') % hold off on drawing the figure until it is all plotted
    for ff = 1:maxExmp
        subplot(r,c,ff);
        curr_frame = curr_frame + 1;
        try  start_frame = confident_frames(curr_frame);
        catch
            if curr_frame>length(confident_frames)
                continue
            end
        end

        % Initialize fly skeletons for plotting later
        % hM = initFlySkeleton(ax(ff), mBaseColor, true, sz);
        % hF = initFlySkeleton(ax(ff), fBaseColor, true, sz);
        ax = gca;
        hM = initFlySkeleton(ax, mBaseColor, true, sz);
        hF = initFlySkeleton(ax, fBaseColor, true, sz);

        roi = allROI(curr_frame,:);% frames to plot
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
            updateFlySkeleton(hM, x, y, mColor);
            x =  fPos(roi(ii),:,1);
            y =  fPos(roi(ii),:,2);
            updateFlySkeleton(hF, x, y, fColor);
        end
        % if any(start_frame==confident_frames) || any(stop_frame==confident_frames)
        %     titleString{ff} = 'Confident';
        % else 
        %     titleString{ff} = '';
        % end
    end

    % format figure: 
    formatFig(fig, blkbgd,[r,c]);
     for ff = 1:maxExmp
        subplot(r,c,ff);
        axis equal tight
        set(gca, 'xcolor', 'none', 'ycolor', 'none')
        if curr_frame>length(confident_frames)
            continue
        end
        title(num2str(confident_frames(curr_frame)), 'color', foreColor,'fontsize', 10)
        % find the size and build a rectagle 'frame'
        ax = axis; % [xmin xmax ymin ymax]
        rectangle('Position',[ax(1) ax(3) ax(2)-ax(1) ax(4)-ax(3)],...
              'EdgeColor',foreColor,'LineWidth',0.5, 'Curvature',[0.1 0.1],...
              'linestyle', ':')
     end
     set(fig, 'Visible', 'on')
     drawnow
    % save_figure(fig, [saveDir fly(trial).name  ' speed error trials ' num2str(gg)],fig_type, true,false, '-r100');
end



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








































































