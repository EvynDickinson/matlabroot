
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

fig = figure; hold on
    plot(D,'color', Color('vaporwavepurple'))
    plot(test, 'color', Color('vaporwaveyellow'))
    h_line(50,'vaporwaveblue', '-',2)
    v_line([B(1),B(end)], 'metrored','-',2)
    xlabel('frame number')
    ylabel('fly speed (mm/s)')
    formatFig(fig, blkbgd);
    title(titleString,'Color',foreColor)
    
    save_figure(fig,[saveDir  ' ' fly(trial).name  ' speed over time trial ' num2str(maxtrial)]);


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


%% ANALYSIS :  can we pair skip frames together?

trial = 8;

% how many frames above 50mm/s?





























































