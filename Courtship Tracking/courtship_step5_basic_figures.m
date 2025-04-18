%% 
% pos (frame, body points, XY)
% body points (head, center, abdomen, right wing, left wing)

%% (DONE) Determine the conversion between pixels to mm for any new video configuration

% vidpath = "S:\Evyn\DATA\Courtship Videos\09.26.2024\Berlin_courtship_F_LRR_caviar_ramp\compiled_video_2.avi";
% movieInfo = VideoReader(vidpath); %read in video
% demoImg = (read(movieInfo,1));
% well_loc = readPoints(demoImg,4); % click on center positions of the wells in the arena
% 
% % distance between well 1 and 3
% d = [];
% well_1 = well_loc(:,1);
% well_3 = well_loc(:,3);
% d = [d;sum((well_1-well_3).^2).^0.5];
% % distance between well 2 and 4
% well_1 = well_loc(:,2);
% well_3 = well_loc(:,4);
% d = [d;sum((well_1-well_3).^2).^0.5];
% 
% pixelsbetweenwells = mean(d); %pixels
% actualdistance =  31.2; %mm
% pix2mm = actualdistance/pixelsbetweenwells; % multiplier

%% TODOs 
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

%% LOAD data
clear; clc;
baseFolder = [getDataPath(6,0),'Trial Data/'];
trialDir = selectFolder(baseFolder); 
baseDir = [baseFolder, trialDir{:} '/']; % full folder directory for that trial
figDir = [baseDir,'Figures/']; 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

load([baseDir, 'basic data.mat']) % load the parameters and temp table
disp('data loaded')

% Experiment parameters
nvids = parameters.nVids; % number of videos
fps = parameters.FPS;
time = T.time;
blkbnd = true;
fig_type = '-png';
[foreColor,backColor] = formattingColors(blkbnd); % get background colors
pix2mm = 0.0289; % calculated on the new back setup 11/7/24

% Create variable to store body points
body = [];
for i = 1:length(parameters.node_names)
    body.(parameters.node_names{i}) = i;
end

M = 1; % male fly index number
F = 2; % female fly index number

% Initial variables
initial_var = who; % who = all variables created so far
initial_var{end+1} = 'initial_var';
initial_var{end+1} = 'well';

disp_fig = false; % display baseline figures?
initial_var{end+1} = 'disp_fig';

%% ANALYSIS: Extract calculated variables
clearvars('-except',initial_var{:})

% Inter-fly-distance from the fly's center point
x1 = m.pos(:,2,1); % x location for male center
y1 = m.pos(:,2,2);
x2 = f.pos(:,2,1);
y2 = f.pos(:,2,2);
% calculate interfly distance
IFD = (sqrt((x1-x2).^2 + (y1-y2).^2)).*pix2mm;
% add interfly distance to the data table
T.IFD = IFD; 

% Fly speed
speed = [];
D = (sqrt((x1(1:end-1)-x1(2:end)).^2 + (y1(1:end-1)-y1(2:end)).^2)).*pix2mm; % male speed
m.speed = [0;(D./(1/fps))];
D = (sqrt((x2(1:end-1)-x2(2:end)).^2 + (y2(1:end-1)-y2(2:end)).^2)).*pix2mm; % female speed
f.speed = [0; (D./(1/fps))];

%% ANALYSIS: Screen for frames with funky wing positions
clearvars('-except',initial_var{:})

% postions: 1-head, 2-center, 3-abdomen, 4-left wing, 5-right wing

% Create variable to hold X and Y data for both male and female
data = [];
data(M).rawX = m.pos(:,:,1);
data(M).rawY = m.pos(:,:,2);
data(F).rawX = f.pos(:,:,1);
data(F).rawY = f.pos(:,:,2);
initial_var{end+1} = 'data';

% For male and female, shift and rotate frame
for sex = 1:2
    % for each loop/sex, switch x and y to appropriate values
    switch sex
        case 1
            sex_type = 'male';
            x = m.pos(:,:,1);
            y = m.pos(:,:,2);
        case 2
            sex_type = 'female';
            x = f.pos(:,:,1);
            y = f.pos(:,:,2);
    end

    % 1) Shift each frame so that the center is at the origin (0,0)
    x_offset = x(:,2); % x values for center
    y_offset = y(:,2);
    % subtract the offset from the x values for each body point
    X = x-x_offset; % new shifted x data set
    Y = y-y_offset;
    
    % 2) Rotate each frame so that the head and center are aligned on the Y axis
    newX = nan([length(X),5]); % establish variable to hold new rotated x data set
    newY = nan([length(Y),5]);
    for frame = 1:length(X)
        % identify coordinates for each body point in that frame
        points = [X(frame,:)',Y(frame,:)'];
        % rotate points so head and center align with y axis
        rotatedPoints = rotateToVertical(points, 2,1,false);
        % load new rotated points into variable
        newX(frame,:) = rotatedPoints(:,1);
        newY(frame,:) = rotatedPoints(:,2);
    end
    
    % Plot body point positions for all frames 
    SZ = 5; r = 1; c = 2;
    if disp_fig
        fig = getfig(sex_type,1,[1058 781]);
        % before corrections plot
        subplot(r,c,1)
            hold on
            CList = {'teal', 'red', 'grey', 'blue', 'gold'};
            % plot each body point using new shifted and rotated data
            for i = 1:5
                scatter(newX(:,i),newY(:,i),SZ, Color(CList{i}),'filled')
            end
            % format figure
            axis square equal
            set(gca, 'xcolor', 'none', 'ycolor','none')
            title('Before Correction')
    end
    % 3) Eliminate frames that are absolutely incorrect
    failstate = newY(:,1)<0; % state where head is below the center
    failstate = [failstate, newY(:,3:5)>0]; % add state where abdomen or wings are above the center point (on head side)
    % identify which frames have a failstate (where there are nonzeros)
    loc = any(failstate,2); 
    % display number of frames to be eliminated for that fly
    disp([sex_type ' fly elimination points: ' num2str(sum(loc))])
    % replace elimination frames with nans
    newX(loc,:) = nan;
    newY(loc,:) = nan;
    
    % 4) Flip R and L wing for mirror image ones
    flipstate = newX(:,5)<newX(:,3) & newX(:,5)<newX(:,4); % state where right wing is to the left of the abdomen and L wing
    flipstate = [flipstate, newX(:,4)>newX(:,3) & newX(:,4)>newX(:,5)]; % add state where left wing is to the right of the abdomen and R wing
    % identify which frames have a failstate (where there are nonzeros)
    loc = any(flipstate,2);
    % display number of wing points to be flipped for that fly
    disp([sex_type ' fly flip wing points: ' num2str(sum(loc))])
    % flip x values for wing points needing flipped
    wing_L = newX(loc,5);
    wing_R = newX(loc,4);
    newX(loc,4:5) = [wing_L,wing_R];
    
    % 5) Eliminate frames with both wings on one side of the abdomen
    flipstate = newX(:,5)<newX(:,3) & newX(:,4)<newX(:,3); % state where wings are to the left of the abdomen
    flipstate = [flipstate, newX(:,5)>newX(:,3) & newX(:,4)>newX(:,3)]; % add state where wings are to the right of the abdomen 
    % identify which frames have a failstate (where there are nonzeros)
    loc = any(flipstate,2);
    % display number of frames to be eliminated for that fly
    disp([sex_type ' fly wing-abdo elimination points: ' num2str(sum(loc))])
    % replace elimination frames with nans
    newX(loc,:) = nan;
    newY(loc,:) = nan;
    % Save the cleaned, adjusted data
    data(sex).x = newX;
    data(sex).y = newY;
    
    if disp_fig
        % Plot fly positions after adjustments: 
        % after corrections plot
        subplot(r,c,2)
            hold on
            CList = {'teal', 'red', 'grey', 'blue', 'gold'};
            % plot each body point using cleaned data
            for i = 1:5
                scatter(newX(:,i),newY(:,i),SZ, Color(CList{i}),'filled')
            end
            % format figure
            axis square equal
            set(gca, 'xcolor','none', 'ycolor','none')
            title('After Correction')

        % Save figures
        save_figure(fig,[figDir  sex_type ' wing position correction'], fig_type);
    end
end


if disp_fig
    % Plot male and female body points on the same figure
     fig = getfig('Male vs Female',1,[1058 781]);
     for sex = 1:2
         switch sex
            case 1
                sex_type = 'Male';
            case 2
                sex_type = 'Female';
        end
        subplot(r,c,sex)
        hold on
            CList = {'teal', 'red', 'grey', 'blue', 'gold'};
            % plot each body point
            for i = 1:5
                scatter(data(sex).x(:,i),data(sex).y(:,i),SZ, Color(CList{i}),'filled')
            end
            % format figure
            axis square equal
            set(gca, 'xcolor','none', 'ycolor','none')
            title(sex_type)
     end
    fig = matchAxis(fig,true);
    
    % Save figure
    save_figure(fig,[figDir  'M vs F wing position scatter'], fig_type);
end


%% FIGURE: Compare wing angles within M and between M and F
clearvars('-except',initial_var{:})

% % Calculate wing angles for male and female
% wing = [];
% for sex = 1:2
%     % for each loop/wing, switch x and y to appropriate values
%     for w = 1:2
%         switch w 
%             case 1
%                 P1 = [data(sex).x(:,4),data(sex).y(:,4)]; % left
%             case 2
%                 P1 = [data(sex).x(:,5),data(sex).y(:,5)]; % right
%         end
% 
%         P2 = [data(sex).x(:,2),data(sex).y(:,2)]; % center
%         P3 = [data(sex).x(:,3),data(sex).y(:,3)]; % abdomen
% 
%         % 1: Calculate vectors from P2 to P1 and from P2 to P3
%         v1 = P1 - P2;  % Nx2 matrix, vector from P2 to P1 for each time step
%         v2 = P3 - P2;  % Nx2 matrix, vector from P2 to P3 for each time step
% 
%         % 2: Calculate the dot product of v1 and v2 for each time step
%         dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);
% 
%         % 3: Calculate the magnitudes of v1 and v2 for each time step
%         mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2);
%         mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2);
% 
%         % 4: Calculate the cosine of the angle for each time step
%         cosTheta = dotProduct ./ (mag_v1 .* mag_v2);
% 
%         % 5: Compute the angle in radians, and then convert to degrees
%         anglesRadians = acos(cosTheta);  % Nx1 vector of angles in radians
%         anglesDegrees = rad2deg(anglesRadians);  % convert to degrees
% 
%         % Save wing angle values as variable
%         data(sex).wingangle(:,w) = anglesDegrees;
%     end
% 
%      % identify where there are nans for wing angle for either wing
%      loc = any(isnan(data(sex).wingangle),2);
%      % create wing spread variable by summing wingangles
%      data(sex).wingspread = sum(data(sex).wingangle,2);
%      % replace wing spread values with nans where wing angles have nans
%      data(sex).wingspread(loc,:) = nan;
% end

% Compare male L and R wing angles
fig = getfig('',true, [1032 300]); 
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
fig = getfig('',true, [1032 300]); 
hold on
    plot(time, data(M).wingspread,'color', Color('dodgerblue'),'linewidth', 1) % male wing spread
    plot(time, data(F).wingspread,'color', Color('deeppink'),'linewidth', 1) % female wing spread
% format figure
xlabel('time (s)')
ylabel('wing angle (\circ)')
formatFig(fig, blkbnd);
save_figure(fig,[figDir 'M and F wing angles'],fig_type);



%% FIGURE: Visualize body angles over time
clearvars('-except',initial_var{:})

% x = 1;
% y = 2;
% head = 1;
% center = 2;
% % positions of M and F head and center
% P1 = [m.pos(:,center,x),m.pos(:,center,y)]; % male center
% P2 = [m.pos(:,head,x),m.pos(:,head,y)]; % male head
% P3 = [f.pos(:,center,x),f.pos(:,center,y)]; % female center
% P4 = [f.pos(:,head,x),f.pos(:,head,y)]; % female head
% 
% % 1: Calculate body vectors
% v1 = P1 - P2;  % Nx2 matrix, vector for male body vector
% v2 = P3 - P4;  % Nx2 matrix, vector for female body vector
% 
% % 2: Calculate the dot product of v1 and v2 for each time step
% dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);
% 
% % 3) Compute the magnitudes of the vectors
% mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
% mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 
% 
% % 4) Calculate the cosine of the angle
% cosTheta = dotProduct ./ (mag_v1 .* mag_v2);
% 
% % 5) Compute the angle in radians and convert to degrees
% angleRadians = acos(cosTheta);  % angle in radians
% angleDegrees = rad2deg(angleRadians);  % convert to degrees
% data(1).mfbodyangle = angleDegrees;
% data(2).mfbodyangle = angleDegrees;

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




%% FIGURES: Male position relative to female
clearvars('-except',initial_var{:})
% center all points to female fly
% align all points to female fly heading (center = 0)

% data(M).color = Color('dodgerblue');
% data(F).color = Color('deeppink');
% 
% % ------------------------- 1) Shift each frame to the origin for female fly ------------------------- 
% x_offset = data(F).rawX(:,body.center); % x values for female center
% y_offset = data(F).rawY(:,body.center);
% mx = data(M).rawX-x_offset; % subtract the offset from the x values for each body point
% my = data(M).rawY-y_offset;
% fx = data(F).rawX-x_offset;
% fy = data(F).rawY-y_offset;
% 
% % ------------------------- 2) Rotate each set of points to head and center on the Y axis ------------------------- 
% [mX,mY,fX,fY] = deal(nan(size(mx)));
% for frame = 1:length(mX)
%     % identify coordinates for each body point in that frame - female only
%     fpoints = [fx(frame,:)',fy(frame,:)']; 
%     % rotate points so head and center align with y axis
%     [rotatedPoints, R] = rotateToVertical(fpoints, body.center,body.head,false);
%     % load new rotated points into variable
%     fX(frame,:) = rotatedPoints(:,1);
%     fY(frame,:) = rotatedPoints(:,2);
%     % identify coordinates for each body point in that frame - male only
%     mpoints = [mx(frame,:)',my(frame,:)'];
%     % rotate points in relation to female (using female R)
%     mpoints = (R * mpoints')';
%     % load new rotated points into variable
%     mX(frame,:) = mpoints(:,1);
%     mY(frame,:) = mpoints(:,2);
% end
% initial_var{end+1} = 'mX';
% initial_var{end+1} = 'mY';
% initial_var{end+1} = 'fX';
% initial_var{end+1} = 'fY';
% 
% % Calulate the angle between fly bodies (both - and +)
% theta = data(M).mfbodyangle; 
% test = theta<90; % when is male less than 90 deg from female

% ------------------------- 3) Demo selected angles in male positions relative to female fly ------------------------- 
zoom = [-250,250];
skip = 20;

fig = getfig('',1,[1075 871]);
hold on

% plot male coordinates for head and body within test constraints
x = mX(test,[body.head,body.center]);
y = mY(test,[body.head,body.center]);
% only plot every [skip value] points to reduce volume
x = x(1:skip:end,:);
y = y(1:skip:end,:);
plot(x',y','color',data(M).color)
scatter(x(:,1),y(:,1),15,data(M).color,"filled","^") % arrow head on male

% plot female coordiates for head, body, and abdomen
x = fX(1:skip:end,[body.head,body.center,body.abdomen]);
y = fY(1:skip:end,[body.head,body.center,body.abdomen]);
plot(x',y','color',data(F).color, 'LineWidth', 2)

% format figure
axis  equal square
h_line(0,'gray',':',2)
v_line(0,'grey',':',2)
xlim(zoom)
ylim(zoom)
formatFig(fig,blkbnd);
set(gca,'XColor','none','YColor','none')

save_figure(fig,[figDir ' all male fly body positions relative to female'],'-png');

% ------------------------- 4) Establish likely and unlikely courtship positions ------------------------- 

% % Determine fly body length (~2.5mm)
% BL = 2.5/pix2mm; % in pixels
% 
% % Body position location rules
% r.a = mX(:,body.center) >= 0; % everything to right of y axis
% r.b = mX(:,body.center) <= 0; % everything to left of y axis
% r.c = mY(:,body.center) >= 0; % everything above x axis
% r.d = mY(:,body.center) <= 0; % everything below x axis
% q1 = r.b & r.c; % quadrant one occupancy
% q2 = r.a & r.c; % quadrant two occupancy
% q3 = r.d & r.b; % quadrant three occupancy
% q4 = r.a & r.d; % quadrant four occupancy
% 
% % Heading direction rules
% r.e = mY(:,body.head) >= mY(:,body.center); % heading direction facing north
% r.f = mY(:,body.head) <= mY(:,body.center); % heading direction facing south
% r.g = mX(:,body.head) >= mX(:,body.center); % heading direction facing east
% r.h= mX(:,body.head) <= mX(:,body.center); % heading direction facing west
% ne = r.e & r.g; % north east direction
% nw = r.e & r.h; % north west direction
% se = r.f & r.g; % south east direction
% sw = r.f & r.h; % south west direction
% 
% 
% % Likely rules
% BLidx = abs(mX(:,body.center)) <= (5*BL) & abs(mY(:,body.center)) <= (5*BL); % limits likely roi's to 5 body lengths
% % 1-4)
% roi1 = q1 & se & BLidx; % quadrant 1, southeast
% roi2 = q2 & sw & BLidx; % quadrant 2, southwest
% roi3 = q3 & ne & BLidx; % quadrant 3, northeast
% roi4 = q4 & nw & BLidx; % quadrant 4, northwest
% 
% % Gray/maybe rules
% r.i = theta < 90;
% r.j = theta > 90;
% deg = [45, 20, 10];
% 
%     % 5) quadrant 1, southwest
%     roi5 = [];
%     for i = 1:3 % each angle
%         % body center is within appropriate BL, body angle is between [deg] and 180, and inside quadrant 1
%         dummy = abs(mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q1 & sw;
%         roi5(:,i) = dummy;
%     end
%     roi5 = any(roi5,2);
% 
%     % 6) quadrant 2, southeast
%     roi6 = [];
%     for i = 1:3
%         dummy = abs(mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q2 & se;
%         roi6(:,i) = dummy;
%     end
%     roi6 = any(roi6,2);
% 
%     % 7) quadrant 3, northwest
%     roi7 = [];
%     for i = 1:3
%         dummy = abs(mX(:,body.center)) <= i*BL & theta < deg(i) & q3 & nw;
%         roi7(:,i) = dummy;
%     end
%     roi7 = any(roi7,2);
% 
%     % 8) quadrant 4, northeast
%     roi8 = [];
%     for i = 1:3
%         dummy = abs(mX(:,body.center)) <= i*BL & theta < deg(i) & q4 & ne;
%         roi8(:,i) = dummy;
%     end
%     roi8 = any(roi8,2);
% 
%     % 9) quadrant 1, northeast
%     roi9 = [];
%     for i = 1:3
%         dummy = abs(mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q1 & ne;
%         roi9(:,i) = dummy;
%     end
%     roi9 = any(roi9,2);
% 
%     % 10) quadrant 2, northwest
%     roi10 = [];
%     for i = 1:3
%         dummy = abs(mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q2 & nw;
%         roi10(:,i) = dummy;
%     end
%     roi10 = any(roi10,2);
% 
%     % 11) quadrant 3, southeast
%     roi11 = [];
%     for i = 1:3
%         dummy = abs(mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q3 & se;
%         roi11(:,i) = dummy;
%     end
%     roi11 = any(roi11,2);
% 
%     % 12) quadrant 4, southwest
%     roi12 = [];
%     for i = 1:3
%         dummy = abs(mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q4 & sw;
%         roi12(:,i) = dummy;
%     end
%     roi12 = any(roi12,2);
% 
% % ROI's that are likely courtship positions
% likely = roi1 | roi2 | roi3 | roi4;
% % ROI's that are maybe courtship positions
% gray = roi5 | roi6 | roi7 | roi8; % exceptions moving towards X axis
% gray2 = roi9 | roi10 | roi11 | roi12; %exceptions moving towards Y axis
% 
% % Saving ROIs
% position = [];
% position.L1 = roi1;
% position.L2 = roi2;
% position.L3 = roi3;
% position.L4 = roi4;
% position.GX1 = roi5;
% position.GX2 = roi6;
% position.GX3 = roi7;
% position.GX4 = roi8;
% position.GY1 = roi9;
% position.GY2 = roi10;
% position.GY3 = roi11;
% position.GY4 = roi12;
% 
% 
% % Create variable for likely and maybe courtship positions
% loc = likely | gray | gray2;
% T.courtposition = loc;  % yes = green, no = red on graph vv
% likelyidx = find(loc);
% unlikelyidx = find(~loc);
% position.likely = likely;
% position.GX = gray;
% position.GY = gray2;
% position.all_likely = loc;
% position.unlikely = ~loc;
% initial_var{end+1} = 'position';

% ------------------------- 5) Visualize likely and unlikely courtship positions ------------------------- 

zoom = [-250,250];

% pull point locations that will be plotted
skip = 20;

likelyidx = likelyidx(1:skip:end);
unlikelyidx = unlikelyidx(1:skip:end);
allidx = [likelyidx; unlikelyidx];

fig = getfig('',1,[1075 871]);
hold on

% % plot all fly positions unlikely courtship male fly body positions
% kolor = Color('grey');
% x = mX(allidx,[body.head,body.center]);
% y = mY(allidx,[body.head,body.center]);
% plot(x',y','color',kolor)
% scatter(x(:,1),y(:,1),15,kolor,"filled","^")
% % axis equal square
% formatFig(fig,blkbnd)
% set(gca,'XColor','none','YColor','none')

% plot the unlikely courtship male fly body positions
kolor = Color('red');
x = mX(unlikelyidx,[body.head,body.center]);
y = mY(unlikelyidx,[body.head,body.center]);
plot(x',y','color',kolor)
scatter(x(:,1),y(:,1),15,kolor,"filled","^")

% plot the likely courtship male fly body positions
kolor = Color('limegreen');
x = mX(likelyidx,[body.head,body.center]);
y = mY(likelyidx,[body.head,body.center]);
plot(x',y','color',kolor)
scatter(x(:,1),y(:,1),15,kolor,"filled","^")

% plot female head, center, and abdoment
x = fX(1:skip:end,[body.head,body.center,body.abdomen]);
y = fY(1:skip:end,[body.head,body.center,body.abdomen]);
plot(x',y','color',foreColor, 'LineWidth', 2)
% xlim([-2000,1500]); ylim([-2000,2000])

% format figure
axis  equal square
h_line(0,'gray',':',2)
v_line(0,'grey',':',2)
xlim(zoom)
ylim(zoom)
formatFig(fig,blkbnd);
set(gca,'XColor','none','YColor','none')

formatFig(fig,blkbnd);
rectangle('Position',[zoom(1),zoom(1) sum(abs(zoom)) sum(abs(zoom))],'edgecolor',foreColor,'linewidth', 1)

save_figure(fig,[figDir 'likely and unlikely male body pos. relative to female'],'-png');
 


%% FIGURE: Distance to food histogram
clearvars('-except',initial_var{:})

% % determine if the well outlines already exist
% well_file = [baseDir 'well locations.mat'];
% if ~exist(well_file, 'file')    
%     % pull up picture to find food well
%     vidpath = [getDataPath(6, 2), parameters.date, '/', parameters.videoName, '/compiled_video_1.avi'];
%     % vidpath = '/Volumes/OnTheGoData/Courtship Videos/09.26.2024/Berlin_courtship_F_LRR_caviar_ramp/compiled_video_1.avi';
%     movieInfo = VideoReader(vidpath); %read in video
%     demoImg = (read(movieInfo,T.vidFrame(1)));
%     img = imadjust(demoImg,[72/255, 215/255]); % adjust the contrast of the image
% 
%     % save food well location
%     txt = {'12', '3', '6', '9'};
%     well = struct; % initialize the new well structure
%     fig = getfig('');
%         imshow(img)
%         for i = 1:4 
%             h = warndlg(['Outline the ' txt{i} ' oclock well']);
%             uiwait(h)
%             roi = drawcircle; % manually add in the circle over the food well
%             well.radius(i) = roi.Radius;
%             well.center(i,:) = roi.Center;
%         end
%     save_figure(fig,[figDir 'well outlines'],'-pdf',0,1, '-r100');
% 
%     well.R = mean(well.radius)*pix2mm;
% 
%     % select the well with the food
%     fig = getfig('');
%         imshow(img); hold on
%         for i = 1:4
%             viscircles(well.center(i,:), well.radius(i),'color',Color('white'));
%             % drawcircle('Center',[well.center(i,:)],'Radius',well.radius(i),'StripeColor','red','Color',Color('grey'));
%         end
%         title('Click the food well','FontSize',18)
%         [xi, yi] = crosshairs(1,{'black','black','yellow','yellow'});      % get a point
%         % find the well with the selected data point
%         % distance between each well center and the selected point
%         [~, index] = min((sqrt((xi-well.center(:,1)).^2 + (yi-well.center(:,2)).^2)).*pix2mm);
% 
%         % plot the selected point
%         scatter(xi,yi,60,"white",'filled')
%         scatter(xi,yi,30,"red",'filled')
%         % change the color of the well outline to highlight selection
%         viscircles(well.center(index,:), well.radius(index),'color',Color('red'));
%     if strcmp(questdlg('okay food well selection?'),'Yes')
%         well.food_idx = index;
%         well.food = well.center(index,:);
%         close(fig)
%     else
%         warndlg('Rerun the food well identification')
%         return
%     end
%     % SAVE DATA?
%     if strcmp(questdlg('save well locations?'),'Yes')
%         save(well_file,'well')
%     else 
%     end
% else
%     load(well_file,'well')
%     disp('Loaded prior well locations')
% end
% 
% % Center of the arena
% WC = well.center';
% N = [];
% x1 = WC(1,1:2:4);
% y1 = WC(2,1:2:4);
% x2 = WC(1,2:2:4);
% y2 = WC(2,2:2:4);
% [xi,yi] = polyxpoly(x1,y1,x2,y2);
% well.center(5,:) = [xi,yi];
% 
% % calculate distance to food (from fly head)
% x1 = m.pos(:,1,1); % x location for male center
% y1 = m.pos(:,1,2);
% x2 = f.pos(:,1,1);
% y2 = f.pos(:,1,2);
% c1 = well.center(1);
% c2 = well.center(2);
% 
% m.dist2food = (sqrt((x1-c2).^2 + (y1-c2).^2)).*pix2mm;
% f.dist2food = (sqrt((c1-x2).^2 + (c1-y2).^2)).*pix2mm;
% T.dist2food = [m.dist2food, f.dist2food];

fig = figure; 
    histogram(T.dist2food(:,M),'FaceColor',data(M).color,'FaceAlpha',0.8)
    hold on
    histogram(T.dist2food(:,F),'FaceColor',data(F).color,'FaceAlpha',0.8)
    xlabel('distance to food (mm)')
    formatFig(fig, blkbnd);
    set(gca,'ycolor', 'none')

% flies ON food
T.FlyOnFood = T.dist2food<=well.R; % fly head must be within the food circle
disp('Flies on food: M & F')
sum(T.FlyOnFood)

%% ANALYSIS: Determine fly occupancy in ring, quadrant, food circle  
clearvars('-except',initial_var{:})

% Max Distance from Center of Arena : 29.612mm = 30mm radius
% ArenaArea = 2827.43;
R = 30; %mm
innerR = R*sqrt(3/4); % radius of the inner 50% occupancy space R*sqrt(1/2)
dist_from_edge = (R - innerR); % distance acceptable for start of outer 25%
maxR = R*sqrt(0.1); % radius of a circle occupying 10% of the arena

% find distance from center for each fly center point: 
xi = well.center(5,1);
yi = well.center(5,2);

% eccentricity, outter ring occupancy of the flies
for sex = 1:2
    x_loc = data(sex).rawX(:,body.center); % x position of the fly
    y_loc = data(sex).rawY(:,body.center); % y position of the fly
    
    %distance from center of arena
    D = sqrt(((x_loc-xi).^2 + (y_loc-yi).^2)).*pix2mm; 
    data(sex).eccentricity = D;

    % outter ring occupancy
    data(sex).OutterRing = D<=R & D>=innerR; % find the locations that are between edge and inner R

    % quadrant occupancy 
    center = well.center(5,:);
    % r = data(exp).data(trial).data.r;
    foodWell = well.food; %data(exp).T.foodLoc(trial);

    % Adjust the X and Y coordinates relative to the new center
    adjustedX = x_loc - center(1);
    adjustedY = y_loc - center(2);
    
    % Initialize matrix to hold quadrant classification (same size as input matrices)
    quadrantMatrix = zeros(size(x_loc));
    
    % Define quadrant masks based on the new center
    Q = [];
    Q(1).Mask = (adjustedY > adjustedX) & (adjustedY <= -adjustedX);  % Top
    Q(2).Mask = (adjustedY <= adjustedX) & (adjustedY <= -adjustedX); % Bottom
    Q(3).Mask = (adjustedY <= adjustedX) & (adjustedY > -adjustedX);  % Left
    Q(4).Mask = (adjustedY > adjustedX) & (adjustedY > -adjustedX);   % Right
    
    % Determine the well locations (which determine quadrant assignment):
    adjusted_wx = foodWell(1) - center(1);
    adjusted_wy = foodWell(2) - center(2);
    
    idx_loc = false(1,4);
    % Find the food quadrant (find location with the food well coordinates included)
    idx_loc(1) = (adjusted_wy > adjusted_wx) & (adjusted_wy <= -adjusted_wx);  % top
    idx_loc(2) = (adjusted_wy <= adjusted_wx) & (adjusted_wy <= -adjusted_wx); % right
    idx_loc(3) = (adjusted_wy <= adjusted_wx) & (adjusted_wy > -adjusted_wx);  % bottom
    idx_loc(4) = (adjusted_wy > adjusted_wx) & (adjusted_wy > -adjusted_wx);   % left
    quad_loc = find(idx_loc);

    fly_loc = ~isnan(x_loc); %gives logical for all fly positions in the position matrix
    foodQuad = Q(quad_loc).Mask & fly_loc; % flies in food quad
    
    data(sex).foodQuad = foodQuad;
    
    % 10% Food Circle
    D = sqrt(((x_loc-foodWell(1)).^2 + (y_loc-foodWell(2)).^2)); % distance from center of food well
    D = D.*pix2mm;
    data(sex).foodcircle = D<=maxR ; % flies within 10% circle around the food

end

% quick plot of when the two flies are in the various regions over time: 


%% ANALYSIS: Determine temperature bins and directions
clearvars('-except',initial_var{:})

initial_var{end+1} = 'tRate';
temp_file = [baseDir 'temp regions.mat'];
if ~exist(temp_file, 'file') 
    
    % manually select the different time points for now TODO automate them
    x = []; 
    fig = getfig; 
    plot(time, T.temperature,'color', 'k')
    title('click the start or the ramp, bottom, and top')
    for i = 1:3
        [xi, ~] = crosshairs(1,{'black','black','red','red'});      % get a point
        [~,xidx] = min(abs(time-xi));
        x(i) = xidx(1);
        v_line(time(x(i)),'r','--',1)
    end
    % find the x-time value for each time period
    tRate = struct;
    tRate(1).idx = [1, x(1)-1];
    tRate(1).name = 'start hold';
    tRate(1).color = Color('grey');
    tRate(2).idx = [x(1), x(2)-1];
    tRate(2).name = 'cooling';
    tRate(2).color = Color('dodgerblue');
    tRate(3).idx = [x(2), x(3)-1];
    tRate(3).name = 'warming';
    tRate(3).color = Color('red');
    tRate(4).idx = [ x(3),T.frame(end)];
    tRate(4).name = 'end hold';
    tRate(4).color = Color('grey');
    
    % plot the regions onto the graph
    hold on
    ylims = ylim;
    for i = 1:4
        roi = time(tRate(i).idx);
        % h = rectangle('Position', [roi(1), ylims(1), roi(2), ylims(2)], 'FaceColor', tRate(i).color);
        h = rectangle('Position', [roi(1), ylims(1), diff(roi),diff(ylims)], 'FaceColor', tRate(i).color);
    end
    plot(time, smooth(T.temperature,fps*10,'moving'),'color', 'k','linewidth', 5)
    xlabel('time (min)')
    ylabel('temperature (\circC)')
    formatFig(fig,false)
    title('temperature conditions')
    
     if strcmp(questdlg('save temperature regions?'),'Yes')
         save_figure(fig, [figDir, 'temperature regions'],'-pdf',1);
         save(temp_file,'tRate')
         disp('saved temperature regions')
    else 
     end
else
    load(temp_file,'tRate')
    disp('loaded temperature regions')
end

% Set up logicals for each of the temperature regions
MT = false(size(time));
%cooling
idx = find(strcmp('cooling',{tRate(:).name}));
T.cooling = MT;
T.cooling(tRate(idx).idx(1):tRate(idx).idx(2)) = true;
%warming
idx = find(strcmp('warming',{tRate(:).name}));
T.warming = MT;
T.warming(tRate(idx).idx(1):tRate(idx).idx(2)) = true;
%hold
idx = find(strcmp('start hold',{tRate(:).name}));
T.hold = MT;
T.hold(tRate(idx).idx(1):tRate(idx).idx(2)) = true;
idx = find(strcmp('end hold',{tRate(:).name}));
T.hold(tRate(idx).idx(1):tRate(idx).idx(2)) = true;

% % Find temperature bins: 
tP = getTempTurnPoints(parameters.protocol);

   



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
        ylabel('temp (\circC')
    subplot(r,c,sb(2).idx); hold on
        for sex = 1:2
            % plot the smoothed average distance
            plot(time, smooth(T.dist2food(:,sex),sSpan, 'moving'), 'color', data(sex).color,'linewidth', lw)
        end
        y_base = rangeLine(fig, 5, false);
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



%% FIGURE: Plot body position during wing extension
% wing angle > 60 deg
% L wing angle in appropriate quadrants
% R wing angle in appropriate quandrants
% matrix for yes courtship wing extension, no not courtship
% diff to see if extension lasts > 1sec

clearvars('-except',initial_var{:})

% Lwing = [];
% Rwing = [];
% % Determine which positions require which wing to be extended
% L_items = {'L1', 'L4', 'GX1', 'GX4', 'GY2', 'GY3'};
% R_items = {'L2', 'L3', 'GX2', 'GX3', 'GY1', 'GY4'};
% % Identify if male is in an appropriate position for each wing direction across each item
% for i = 1:length(L_items)
%     Lwing = [Lwing, position.(L_items{i})];
%     Rwing = [Rwing, position.(R_items{i})];
% end
% % Condense to identify if male is in any of the appropriate positions
% Lwing = any(Lwing,2);
% Rwing = any(Rwing,2);
% 
% % Pull wing angles equal or greater than extension minimum for L and R
% wa_cutoff = 50; % minimum wing extension angle for courtship
% wing_ext = (Lwing & (m.wing.angle(:,1) >= wa_cutoff)) | (Rwing & (m.wing.angle(:,2) >= wa_cutoff)); % wing must be facing the female fly
% % Each value subtracted by the value before it (1 = ext starts, -1 = ext stops, 0 = no state change)
% a = diff(wing_ext); 
% % Add the first extension value to the list to account for the starting condition
% b = [wing_ext(1); a]; 
% % Locations in wing_ext where extension period starts/end
% ext_start = find(b == 1); 
% ext_stop = find(b == -1);
% % If wing ext doesn't stop by end, add stop location at end of ext_stop (loc = length of experiment value)
% if wing_ext(end)
%     ext_stop(end + 1) = length(time);
% end
% % Calculate the length of each wing ext bout
% ext_dur = ext_stop - ext_start;
% % Find where wing ext lasts longer than 1sec
% dur_loc = find(ext_dur > fps);
% 
% % Create new courtship matrix with only true wing ext for bouts longer than 1sec
% mt = false(size(time));
% for i = 1:length(dur_loc)
%     ii = dur_loc(i);
%     mt(ext_start(ii):ext_stop(ii)) = true;
% end
% T.wing_ext = mt;
% T.wing_ext_all = wing_ext;

% demo images of wing extension:

frames = find(T.wing_ext==1);
endidx = (find(diff(frames)>1)); % where the does first 'extension' bout end?
try 
    roi = frames(1):frames(endidx(1));
catch
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
% 
% xlim(xlims)
% ylim(ylims)
% % 
% xlims = xlim;
% ylims = ylim;

%% ANALYSIS: Chase identification
% < 120 deg area behind female x
% facing female x
% 7mm between M center and F center x
% female speed > 0 x
% if all are 1 then courtship
% diff to see if chasing lasts > 2sec x

clearvars('-except',initial_var{:})

x = 1;
y = 2;
head = 1;
center = 2;
% positions of M head and F head and center
P1 = [m.pos(:,head,x),m.pos(:,head,y)]; % male head
P2 = [f.pos(:,center,x),f.pos(:,center,y)]; % female center
P3 = [f.pos(:,head,x),f.pos(:,head,y)]; % female head

% 1) Calculate body vectors
v1 = P3 - P1;  % Nx2 matrix, vector for female head to male head
v2 = P3 - P2;  % Nx2 matrix, vector for female head to female center

% 2) Calculate the dot product of v1 and v2 for each time step
dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);

% 3) Compute the magnitudes of the vectors
mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 

% 4) Calculate the cosine of the angle
cosTheta = dotProduct ./ (mag_v1 .* mag_v2);

% 5) Compute the angle in radians and convert to degrees
angleRadians = acos(cosTheta);  % angle in radians
angleDegrees = rad2deg(angleRadians);  % convert to degrees
mfpos_angle = angleDegrees;

% Identify when male position angle is less than 60 degrees from female
pos_angle = abs(mfpos_angle) <= 60;

% Identify when male is facing female
facing = [];
m_items = {'L3', 'L4', 'GX3', 'GX4', 'GY3', 'GY4'};
% Identify if male is in an appropriate position for each wing direction across each item
for i = 1:length(m_items)
    facing = [facing, position.(m_items{i})];
end
facing = any(facing,2);

% Identify when male is behind female AND facing her
mbehindf = (facing & pos_angle);

% Identify when male is within 7mm of female
close_dist = T.IFD <= 7; % mm

% Identify when female is moving
fmoving = f.speed >= 0.1; % min speed up for debate

chase = (mbehindf & close_dist & fmoving);
a = diff(chase); 
% Add the first chase value to the list to account for the starting condition
b = [chase(1); a]; 
% Locations in chase where chasing period starts/end
ch_start = find(b == 1); 
ch_stop = find(b == -1);
% If chasing doesn't stop by end, add stop location at end of ch_stop (loc = length of experiment value)
if chase(end)
    ch_stop(end + 1) = length(time);
end
% Calculate the length of each chasing bout
ch_dur = ch_stop - ch_start;
% Find where chasing lasts longer than 2sec
dur_loc = find(ch_dur > (2*fps));

% Create new courtship matrix with only true chasing bouts longer than 2sec
mt = false(size(time));
if isempty(dur_loc)
    m.chaseroi = [];
else
    for i = 1:length(dur_loc)
        ii = dur_loc(i);
        mt(ch_start(ii):ch_stop(ii)) = true;
        m.chaseroi(i,:) = [ch_start(ii), ch_stop(ii)];
    end
end
T.court_chase = mt; % time restriction 2 seconds
T.chase_all = chase; % NO time limit



%% FIGURE: M body positions during chase
% Pull point locations that will be plotted
skip = 20;
zoom = [-250,250];

% screening = close_dist;

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
save_figure(fig,[figDir 'chase positions M fly'],fig_type,1,0);

%% FIGURE: Chase overlaid on arena with temp timecourse

switch questdlg(['Show all ' num2str(size(m.chaseroi,1)) 'instances of chase?'])
    case 'Yes'
    case 'No'
        return
    case 'Cancel'
        return
end
% Subplots
r = 5;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2:5;

% For each bout of chase, plot M and F movement on arena image
for i = 1:size(m.chaseroi,1)
    % Frame numbers at start and end of each chase bout
    plotroi = m.chaseroi(i,1):m.chaseroi(i,2);
    % Identify which frame (last in each bout) and video to pull 
    frame = m.chaseroi(i,2);
    vidnum = T.vidNums(frame);
    
    % Pull and read in video
    vidpath = [getDataPath(6, 2), parameters.date, '\', parameters.videoName, '\compiled_video_', num2str(vidnum), '.avi'];
    movieInfo = VideoReader(vidpath);
    demoImg = (read(movieInfo,T.vidFrame(frame)));
    img = imadjust(demoImg,[72/255, 180/255]);

    % Scatter point size and linewidth
    sz = 10;
    lw = 1;
    % Plot body position at every __ frame
    frame_skip = 1;

    % Create x and y limit variable for zoomed figures
    [xlimits, ylimits] = deal([]);

    % Plot body positions over chase bout
    fig = getfig(' ', true, [759 900]); 
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
        xlabel('time (s)')
        ylabel('\circC')
        xlim([0,time(end)])
    % 2) Body positions overlaid on arena image 
    subplot(r, c, sb(2).idx)
        % Plot arena image
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
    
    % Format figure
    formatFig(fig, blkbnd, [r,c], sb);
    % lineROI = drawline(gca);
    
    % Save figure
    save_figure(fig,[figDir 'Chase_', num2str(i), ' ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type,false, false);
    
    % ---------------------------------------- Zoom in on arena and display skeletons -----------------------------------
        % Zoom in on the flies
        xlimits(1) = min([x1; x2; m.pos(frame, :,1)'; f.pos(frame, :,1)']);
        xlimits(2) = max([x1; x2; m.pos(frame, :,1)'; f.pos(frame, :,1)']);
        ylimits(1) = min([y1;y2; m.pos(frame, :,2)'; f.pos(frame, :,2)']);
        ylimits(2) = max([y1;y2; m.pos(frame, :,2)'; f.pos(frame, :,2)']);
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
    save_figure(fig,[figDir 'Chase_zoom_', num2str(i), ' ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type);
end

%% FIGURES: Visualize chasing overlaid on arena image, with side zoom on other parameters 
clearvars('-except',initial_var{:})

% Subplots
r = 7;
c = 9;
sb(1).idx = 2:4; %  temperature
sb(2).idx = [19:23, 28:32, 37:41, 46:50, 55:59]; % arena image
sb(3).idx = 6:9; % zoom in temp 
sb(4).idx = [15:18, 24:27]; % distance between flies
sb(5).idx = [33:36, 42:45]; % speed correlation
sb(6).idx = [51:54, 60:63]; % male wing angles

timebuff = 3; % time buffer before and after chase (in seconds)
timebuff = timebuff/60; % set to minutes

% Plot body movements and timecourses for each chase bout
for i = 1:size(m.chaseroi,1)
    % Establish x limits for non-temp timecourses
    xlimit = [time(m.chaseroi(i,1))-timebuff,time(m.chaseroi(i,2))+timebuff];

    % Frame numbers at start and end of each chase bout
    plotroi = m.chaseroi(i,1):m.chaseroi(i,2);
    % Identify which frame (last in each bout) and video to pull
    frame = m.chaseroi(i,2);
    vidnum = T.vidNums(frame);
    
    vidpath = [getDataPath(6, 2), parameters.date, '\', parameters.videoName, '\compiled_video_', num2str(vidnum), '.avi'];
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
        plot(time,data(M).wingangle(:,1),'color', Color('Dodgerblue'),'LineWidth', lw)
        plot(time,data(M).wingangle(:,2),'color', Color('Cyan'),'LineWidth', lw)
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
    save_figure(fig,[figDir 'Chase Bout_', num2str(i), ' from ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type,false, false);

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
    save_figure(fig,[figDir 'Chase Bout Zoom_', num2str(i), ' from ', num2str(time(ROI(1))) ' to '  num2str(time(ROI(end)))], fig_type);
end

%% ANALYSIS: Circling behavior
% TODO: add some visualization image for the periods of circling
% M head within 3mm? [head_dist]
% M facing female [position.likely]
% M velocity constant within 1sec [const_var]
% F not moving
% if all are 1 then courtship
% diff to see if circling lasts > 1sec

clearvars('-except',initial_var{:})

% Distance between male head and female
x1 = m.pos(:,1,1); % x location for male head
y1 = m.pos(:,1,2);
x2 = f.pos(:,2,1); % x location for female center
y2 = f.pos(:,2,2);
% Calculate male head - female center distance
d = ((sqrt((x1-x2).^2 + (y1-y2).^2)).*pix2mm);
head_dist = d <= 3; % mm
y = nan(size(head_dist));
y(head_dist) = 1;

% (for fun) how much time did the male fly spend really close to the female fly?
percent_time_within_3mm = (sum(head_dist)/length(time))*100;
disp(['Male spend ' num2str(percent_time_within_3mm) '% of experiment close to female'])

% position.all_likely

% female fly is not moving
f_speed_cut = (f.speed<=0.2);

% Determine if semiconstant male velocity
% 0.5sec smooth, then 2 sec std
binwidth = fps; %  1 second bins to look for avg speed
n = length(m.speed) - binwidth; % number of iterations
offset = ceil(binwidth * 0.5);
smoothspeed = smooth(m.speed,floor(0.5*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.5*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.3*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.2*fps),'moving');
sp_var = nan(size(time)); % initialize the speed variability variable
for i = 1:n % loop through each bin
    roi = i:i+binwidth; % frames from i to one second later
    sp_var(i+offset) = std(smoothspeed(roi));
end

testlim = 0.7;
const_var = sp_var>testlim;
ok_var = m.speed;
ok_var(const_var) = nan;% all data points with acceptable variance

% -- View the threshold for consistent speed --
% r = 2;
% c = 1;
% 
% fig = figure('Position',[1 10 1478 845]);
% subplot(r,c,1); hold on
%     plot(time,m.speed,'color', Color('red'))
%     plot(time, ok_var, 'color', foreColor)
% subplot(r,c,2)
%     plot(time,sp_var,'color', foreColor)
% % formatting
% xlimits = [1,2];
% subplot(r,c,1)
%     xlim(xlimits)
%     ylabel('male fly speed (mm/s)')
% subplot(r,c,2)
%     xlim(xlimits)
%     h_line(testlim, 'r','-',1)
%     ylabel('speed variance')
% formatFig(fig,blkbnd,[r,c])

% Full selection criteria: 
V = position.likely & head_dist & const_var & f_speed_cut;


%  --- Find periods longer than 1 second ---
a = diff(V); % when does the speed switch between stability and instability 
% Add the first chase value to the list to account for the starting condition
b = [V(1); a]; 
% Find when stability periods starts/end
v_start = find(b == 1); 
v_stop = find(b == -1);
% If speed stability doesn't stop by end, add stop location at end of v_stop (loc = length of experiment value)
if V(end)
    v_stop(end + 1) = length(time);
end
% Calculate the length of each speed bout
v_dur = v_stop - v_start;
% Find where speed stability lasts longer than 1sec
dur_loc = find(v_dur > fps);
% save the constant speed regions into a new metric
constant_velocity = false(size(time));
for i = 1:length(dur_loc)
    idx = dur_loc(i);
    constant_velocity(v_start(idx): v_stop(idx)) = true;
end

T.circling_all = V; % when the male fly is circling the female (no time restriction)
T.circling_1sec = constant_velocity; % when circling is longer than 1 second


%% FIGURE: (TODO) M body positions relative to female during circling
% Pull point locations that will be plotted
skip = 20;
zoom = [-250,250];

% find the longest circling behavior example: 
[~,idx] = max(smooth(T.circling_all,30,'moving'));
roi = idx-fps:idx+fps;
figure; plot(roi,T.circling_all(roi))
roi = 94845:94845+fps;
% 
% frames = find(T.wing_ext==1);
% endidx = (find(diff(frames)>1)); % where the does first 'extension' bout end?
% roi = frames(1):frames(endidx(1));

fig = getfig('',1); hold on
% plot the female fly
for ff = 1:length(roi)
    frame = roi(ff);
    for sex = 1:2
        x = data(sex).rawX(frame,:);
        y = data(sex).rawY(frame,:);
        kolor = data(sex).color;
        if T.circling_all(frame) && sex==M
            plotFlySkeleton(fig, x,y,foreColor,false);
        else
           plotFlySkeleton(fig, x,y,kolor,false);
        end
        scatter(x,y, 15, Color('grey'),'filled')
        scatter(x(body.head),y(body.head), 35, Color('yellow'),'filled', 'Marker','^')

    end
end
formatFig(fig,blkbnd);
set(gca, 'xcolor', 'none', 'ycolor', 'none')

save_figure(fig, [figDir, 'Wing extension example 1'],'-png');

xlim(xlims)
ylim(ylims)
% 



% screening = close_dist;

fig = getfig('',1,[1075 871]);
hold on
% Plot all male body positions
kolor = Color('grey');
x = mX(1:skip:end,[body.head,body.center]);
y = mY(1:skip:end,[body.head,body.center]);
plot(x',y','color',kolor)
scatter(x(:,1),y(:,1),15,kolor,"filled","^")
% axis equal square
formatFig(fig,blkbnd)
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
formatFig(fig,blkbnd)
set(gca,'XColor','none','YColor','none')

formatFig(fig,blkbnd)
% rectangle('Position',[zoom(1),zoom(1) sum(abs(zoom)) sum(abs(zoom))],'edgecolor',foreColor,'linewidth', 1)
% save_figure(fig, 'G:\My Drive\Jeanne Lab\Presentations\Data Presentation 12.6.2024\All Chase Positions',fig_type);

% Save figure
save_figure(fig,[figDir 'Circling positions M fly'],fig_type);

%% FIGURE: CI & behavior components
clearvars('-except',initial_var{:})

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


%% FIGURE: Fly turning over time 
clearvars('-except',initial_var{:})
% for sex = 1:2
%     % extract the x and y head and center positions to calculcate the slope of the body line
%     x = data(sex).rawX(:,body.head:body.center);
%     y = data(sex).rawY(:,body.head:body.center);
% 
%     % zero the flys center (actually unneccessary) 
%     X = x - x(:,2);
%     Y = y - y(:,2);
%     % slope for each point in time
%     slope = (Y(:,2)-Y(:,1))./(X(:,2)-X(:,1));
%     m1 = slope(1:end-1);   % 'past' heading
%     m2 = slope(2:end);      % 'current' heading
% 
%     % calculate the angle between the two slopes 
%     theta = atan((m1-m2)./(1+(m1.*m2)));
%     theta = rad2deg(theta);
%     data(sex).turning = [nan; theta].*(fps);
% end

fig = getfig('',1); hold on
for sex = 1:2
    plot(time,smooth(data(sex).turning,fps,'moving'),'color', data(sex).color)
end
h_line(0,'grey','--')
xlabel('time (min)')
ylabel('turning (\circ/s)')

%% FIGURE: Plot fly positions for given frames
for frame = 1:20

% plot fly positions: 
fig = getfig('',1,[560 420]);

    hold on
    % for frame = frames
        plotFlySkeleton(fig, data(1).rawX(frame, :),data(1).rawY(frame, :),Color('dodgerblue'),true);
        plotFlySkeleton(fig, data(2).rawX(frame, :),data(2).rawY(frame, :),Color('deeppink'),true);
    % end
    axis square equal
set(gca, 'xcolor', backColor,'ycolor',backColor)
title([num2str(data(1).mfbodyangle(frame))])
end

%% FIGURE: Sleep
clearvars('-except',initial_var{:})

% bout = 5*60*parameters.FPS;
% dummy = [];
% 
% % Extract sleep bouts from position data
% for sex = 1:2
%     switch sex
%         case 1
%             x = m.pos(:,2,1); % male center
%         case 2
%             x = f.pos(:,2,1); % female center
%     end
%     % Calculate difference between all x values
%     x_diff = diff(x); 
%     % Identify when position is not changing
%     u = abs(x_diff)<= 1;
%     % Each value subtracted by the value before it (1 = ext starts, -1 = ext stops, 0 = no state change)
%     a = diff(u);
%     % Add the first position value to the list to account for the starting condition
%     b = [u(1); a]; 
%     % Frames where 'position-no-change' period starts/end
%     slp_start = find(b == 1); 
%     slp_stop = find(b == -1);
%     % If sleep doesn't stop by end, add stop location at end of slp_stop
%     if u(end)
%         slp_stop(end + 1) = length(time);
%     end
%     % Calculate the length of each 'position-no-change' bout
%     slp_dur = slp_stop - slp_start;
%     % Find where bout lasts longer than 5min (when is sleep)
%     slp_loc = find(slp_dur > bout);
% 
%     % Create dummy matrix with only true sleep bouts
%     mt = false(size(time));
%     for i = 1:length(slp_loc)
%         ii = slp_loc(i);
%         mt(slp_start(ii):slp_stop(ii)) = true;
%     end
%     dummy(sex).sleep = mt;
% end
% 
% % Save sleep data
% m.sleep = dummy(1).sleep;
% f.sleep = dummy(2).sleep;

% Figure
lw = 2;

fig = getfig;
hold on
    sSpan = fps; % single second smoothing
    y1 = smooth(m.sleep,sSpan,'moving'); % male
    y2 = smooth(f.sleep,sSpan,'moving'); % female
    % plot(time, y1, 'color', Color('dodgerblue'),'LineWidth', lw)
    plot(time, y2, 'color', Color('deeppink'),'LineWidth', lw)
    plot(time, y1, 'color', Color('dodgerblue'),'LineWidth', lw)
    ylabel('sleep(mm/s)')
    xlabel('time(min)')

formatFig(fig);
ylim([0,2])



%%





















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
    % xlabel('time (s)')
    ylabel('inter-fly distance (mm)')

% FLY SPEED
subplot(r,c,3); hold on 
    sSpan = fps; %single second smoothing
    y1 = smooth(m.speed,sSpan,'moving'); % male
    y2 = smooth(f.speed,sSpan,'moving'); % female
    plot(time, y1, 'color', Color('dodgerblue'),'LineWidth', lw)
    plot(time, y2, 'color', Color('deeppink'),'LineWidth', lw)

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

    ylabel('speed (mm/s)')

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

xlimit = [0 64];
for i = 1:4
    subplot(r,c,i)
    set(gca, 'xcolor', 'none')
    xlim(xlimit)
end
subplot(r,c,5)
xlim(xlimit)


save_figure(fig,[figDir, 'Full timecourse zoom in'], fig_type);



%% FIGURE: Distance, speed, and speed correcation between flies over full video course
% use the center body point to determine between-fly distance
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

save_figure(fig,[figDir 'speed and distance timecourse'], fig_type);


%% FIGURE: Plot fly positions for a given frame
frame = 1;

% plot fly positions: 
fig = getfig('',1,[560 420]);
    hold on
    plotFlySkeleton(fig, m.pos(frame, :,1),m.pos(frame, :,2),Color('dodgerblue'),true);
    plotFlySkeleton(fig, f.pos(frame, :,1),f.pos(frame, :,2),Color('deeppink'),true);
    axis square equal
set(gca, 'xcolor', backColor,'ycolor',backColor)


%% FIGURE: Distance, speed, and speed correcation between flies over full video course
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
