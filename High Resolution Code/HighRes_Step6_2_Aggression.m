
%% ANALYSIS: Determine male position relative to female

% Adapted from 5_1 analysis to determine possible courtship positions

clearvars('-except',initial_var{:})
pix2mm = conversion(4).pix2mm; %updates 5.12.25

% Create aggression structure and pull position data relative to female fly from fly structure 
aggression = [];
for trial = 1:num.trials
    aggression(trial).name = fly(trial).name;
    aggression(trial).MfacingF = fly(trial).position;
end
initial_var{end+1} = 'aggression';


%% FIGURES: Ground truthing 
foreColor = formattingColors(blkbgd); % get background colors
Mcolor = Color('dodgerblue');
Fcolor = Color('deeppink');

% hard code trial number for ground truthing
trial = 1;

% Demo random male positions relative to female fly
zoom = [-250,250];

% adjust skip size based on total number of points
displayNum = 50; 
xM = fly(trial).mX;
yM = fly(trial).mY;
xF = fly(trial).fX;
yF = fly(trial).fY;

plotLoc = round(linspace(1,length(xM),displayNum)); % equally spaced frames throughout the experiment

% FIGURE
fig = getfig('Random selection of all male positions relative to female',1,[1075 871]);
hold on

% plot male coordinates for head and body 
x = xM(plotLoc,[body.head,body.center]);
y = yM(plotLoc,[body.head,body.center]);
plot(x',y','color',Mcolor)
scatter(x(:,1),y(:,1),15,Mcolor,"filled","^") % arrow head on male

x = xF(plotLoc,[body.head,body.center]);
y = yF(plotLoc,[body.head,body.center]);
plot(x',y','color',Fcolor, 'LineWidth', 2)

% format figure
axis  equal square
h_line(0,'gray',':',2)
v_line(0,'grey',':',2)
formatFig(fig,blkbgd);
set(gca,'XColor','none','YColor','none')

% ____________________________________________________________________________________________________

% Visualize likely and unlikely positions 
likelyFrames = find(aggression(trial).MfacingF.all_likely); % find likely positions
rand_frames_finder = round(linspace(1,length(likelyFrames),displayNum)); % pull equally spaced frames from likely total
yah_frames = likelyFrames(rand_frames_finder); % get frame locations for random likely 

unlikelyFrames = find(aggression(trial).MfacingF.unlikely); % find unlikely positions
rand_frames_finder = round(linspace(1,length(unlikelyFrames),displayNum)); % pull equally spaced frames from unlikely total
nah_frames = unlikelyFrames(rand_frames_finder); % get frame locations for random unlikely

% FIGURE
fig = getfig('Random selection of all male positions relative to female',1,[1075 871]);
hold on

% plot male coordinates for head and body 
% unlikely positions
x = xM(nah_frames,[body.head,body.center]);
y = yM(nah_frames,[body.head,body.center]);
plot(x',y','color',Color('red'))
scatter(x(:,1),y(:,1),15,Color('red'),"filled","^") % arrow head on male
% likely positions
x = xM(yah_frames,[body.head,body.center]);
y = yM(yah_frames,[body.head,body.center]);
plot(x',y','color',Color('limegreen'))
scatter(x(:,1),y(:,1),15,Color('limegreen'),"filled","^") % arrow head on male

% plot female coordinates for head and body
x = xF(plotLoc,[body.head,body.center]);
y = yF(plotLoc,[body.head,body.center]);
plot(x',y','color',foreColor,'LineWidth', 2)

% format figure
axis  equal square
h_line(0,'gray',':',2)
v_line(0,'grey',':',2)
formatFig(fig,blkbgd);
set(gca,'XColor','none','YColor','none')


%% ANALYSIS: Calculate M aggressive wing extension

clearvars('-except',initial_var{:})

Lwing = 1;
Rwing = 2; 

ag_loc = []; % instances of aggression
min_wa = 50:5:90;


for trial = 1:num.trials
    facing = aggression(trial).MfacingF.all_likely;
    for angle = 1:length(min_wa)
        % Pull wing angles equal or greater than extension minimum for L and R
        wa_cutoff = min_wa(angle); % minimum wing extension angle for aggression
        wing_ext = ((fly(trial).data(M).wingangle(:,Lwing)) >= wa_cutoff & (fly(trial).data(M).wingangle(:,Rwing)) >= wa_cutoff); % both wings must be > 90 deg
        wing_ext = wing_ext & facing; % wing angle and facing F requirements met
        if isempty(wing_ext)
            disp(['no instances of aggression in ' aggression(trial).name])
        end
        
        % Each value subtracted by the value before it (1 = ext starts, -1 = ext stops, 0 = no state change)
        a = diff(wing_ext); 
        % Add the first extension value to the list to account for the starting condition
        b = [wing_ext(1); a]; 
        % Locations in wing_ext where extension period starts/end
        ext_start = find(b == 1); 
        ext_stop = find(b == -1);
        % If wing ext doesn't stop by end, add stop location at end of ext_stop (loc = length of experiment value)
        if wing_ext(end)
            ext_stop(end + 1) = length(fly(trial).time);
        end
        % Calculate the length of each wing ext bout
        ext_dur = ext_stop - ext_start;
        % Find where wing ext lasts longer than 1sec
        dur_loc = find(ext_dur > fps);
        
        % Create new aggression matrix with only true wing ext for bouts longer than 1sec
        mt = false(size(fly(trial).time));
        for i = 1:length(dur_loc)
            ii = dur_loc(i);
            mt(ext_start(ii):ext_stop(ii)) = true;
        end
        aggression(trial).one_sec.(['frames_' num2str(wa_cutoff)]) = mt;
        aggression(trial).all_time.(['frames_' num2str(wa_cutoff)]) = wing_ext;
        aggression(trial).instances.(['frames_' num2str(wa_cutoff)]) = ...
            find(aggression(trial).all_time.(['frames_' num2str(wa_cutoff)]));
    end
end

initial_var{end+1} = 'min_wa';

%% FIGURE: Histogram of agressive wing angles

clearvars('-except',initial_var{:})

plotdata = [];
y = [];

% Concatenate all instances in each wing angle across all trials
for angle = 1:length(min_wa)
    loc = [];
    for trial = 1:num.trials    
        wa_cutoff = min_wa(angle);
        idx = ['frames_' num2str(wa_cutoff)];
        loc = autoCat(loc,aggression(trial).instances.(idx),true,true);
    end
    % Calculate number of instances in for each angle
    y = [y,size(loc,1)];
end

% FIGURE
fig = getfig;
x = min_wa; % 50:5:90
% Plot bar plot
bar(x,y,'FaceColor', Color('FireBrick'))

% Format figure
formatFig(fig,blkbgd);
xlabel('Angles of individual mirrored wings')
ylabel('Total instances across trials')

save_figure(fig,[figDir 'number of aggressive wing angle instances'],fig_type);

