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

%% ANALYSIS:Screen for frames with funky wing positions
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

%% ANALYSIS: Calculate M and F wing angles

clearvars('-except',initial_var{:})

% Calculate wing angles for male and female
wing = [];
for sex = 1:2
    % for each loop/wing, switch x and y to appropriate values
    for w = 1:2
        switch w 
            case 1
                P1 = [data(sex).x(:,4),data(sex).y(:,4)]; % left
            case 2
                P1 = [data(sex).x(:,5),data(sex).y(:,5)]; % right
        end
        
        P2 = [data(sex).x(:,2),data(sex).y(:,2)]; % center
        P3 = [data(sex).x(:,3),data(sex).y(:,3)]; % abdomen
        
        % 1: Calculate vectors from P2 to P1 and from P2 to P3
        v1 = P1 - P2;  % Nx2 matrix, vector from P2 to P1 for each time step
        v2 = P3 - P2;  % Nx2 matrix, vector from P2 to P3 for each time step
        
        % 2: Calculate the dot product of v1 and v2 for each time step
        dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);
        
        % 3: Calculate the magnitudes of v1 and v2 for each time step
        mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2);
        mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2);
        
        % 4: Calculate the cosine of the angle for each time step
        cosTheta = dotProduct ./ (mag_v1 .* mag_v2);
        
        % 5: Compute the angle in radians, and then convert to degrees
        anglesRadians = acos(cosTheta);  % Nx1 vector of angles in radians
        anglesDegrees = rad2deg(anglesRadians);  % convert to degrees
        
        % Save wing angle values as variable
        data(sex).wingangle(:,w) = anglesDegrees;
    end

     % identify where there are nans for wing angle for either wing
     loc = any(isnan(data(sex).wingangle),2);
     % create wing spread variable by summing wingangles
     data(sex).wingspread = sum(data(sex).wingangle,2);
     % replace wing spread values with nans where wing angles have nans
     data(sex).wingspread(loc,:) = nan;
end

%% ANALYSIS: Calculate the body angle between the flies

clearvars('-except',initial_var{:})

x = 1;
y = 2;
head = 1;
center = 2;
% positions of M and F head and center
P1 = [m.pos(:,center,x),m.pos(:,center,y)]; % male center
P2 = [m.pos(:,head,x),m.pos(:,head,y)]; % male head
P3 = [f.pos(:,center,x),f.pos(:,center,y)]; % female center
P4 = [f.pos(:,head,x),f.pos(:,head,y)]; % female head

% 1: Calculate body vectors
v1 = P1 - P2;  % Nx2 matrix, vector for male body vector
v2 = P3 - P4;  % Nx2 matrix, vector for female body vector

% 2: Calculate the dot product of v1 and v2 for each time step
dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);

% 3) Compute the magnitudes of the vectors
mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 

% 4) Calculate the cosine of the angle
cosTheta = dotProduct ./ (mag_v1 .* mag_v2);

% 5) Compute the angle in radians and convert to degrees
angleRadians = acos(cosTheta);  % angle in radians
angleDegrees = rad2deg(angleRadians);  % convert to degrees
data(1).mfbodyangle = angleDegrees;
data(2).mfbodyangle = angleDegrees;

%% ANALYSIS: Determine male position relative to female

clearvars('-except',initial_var{:})
% center all points to female fly
% align all points to female fly heading (center = 0)

data(M).color = Color('dodgerblue');
data(F).color = Color('deeppink');

% ------------------------- 1) Shift each frame to the origin for female fly ------------------------- 
x_offset = data(F).rawX(:,body.center); % x values for female center
y_offset = data(F).rawY(:,body.center);
mx = data(M).rawX-x_offset; % subtract the offset from the x values for each body point
my = data(M).rawY-y_offset;
fx = data(F).rawX-x_offset;
fy = data(F).rawY-y_offset;

% ------------------------- 2) Rotate each set of points to head and center on the Y axis ------------------------- 
[mX,mY,fX,fY] = deal(nan(size(mx)));
for frame = 1:length(mX)
    % identify coordinates for each body point in that frame - female only
    fpoints = [fx(frame,:)',fy(frame,:)']; 
    % rotate points so head and center align with y axis
    [rotatedPoints, R] = rotateToVertical(fpoints, body.center,body.head,false);
    % load new rotated points into variable
    fX(frame,:) = rotatedPoints(:,1);
    fY(frame,:) = rotatedPoints(:,2);
    % identify coordinates for each body point in that frame - male only
    mpoints = [mx(frame,:)',my(frame,:)'];
    % rotate points in relation to female (using female R)
    mpoints = (R * mpoints')';
    % load new rotated points into variable
    mX(frame,:) = mpoints(:,1);
    mY(frame,:) = mpoints(:,2);
end
initial_var{end+1} = 'mX';
initial_var{end+1} = 'mY';
initial_var{end+1} = 'fX';
initial_var{end+1} = 'fY';

% Calulate the angle between fly bodies (both - and +)
theta = data(M).mfbodyangle; 
test = theta<90; % when is male less than 90 deg from female


% ------------------------- 4) Establish likely and unlikely courtship positions ------------------------- 

% Determine fly body length (~2.5mm)
BL = 2.5/pix2mm; % in pixels

% Body position location rules
r.a = mX(:,body.center) >= 0; % everything to right of y axis
r.b = mX(:,body.center) <= 0; % everything to left of y axis
r.c = mY(:,body.center) >= 0; % everything above x axis
r.d = mY(:,body.center) <= 0; % everything below x axis
q1 = r.b & r.c; % quadrant one occupancy
q2 = r.a & r.c; % quadrant two occupancy
q3 = r.d & r.b; % quadrant three occupancy
q4 = r.a & r.d; % quadrant four occupancy

% Heading direction rules
r.e = mY(:,body.head) >= mY(:,body.center); % heading direction facing north
r.f = mY(:,body.head) <= mY(:,body.center); % heading direction facing south
r.g = mX(:,body.head) >= mX(:,body.center); % heading direction facing east
r.h= mX(:,body.head) <= mX(:,body.center); % heading direction facing west
ne = r.e & r.g; % north east direction
nw = r.e & r.h; % north west direction
se = r.f & r.g; % south east direction
sw = r.f & r.h; % south west direction


% Likely rules
BLidx = abs(mX(:,body.center)) <= (5*BL) & abs(mY(:,body.center)) <= (5*BL); % limits likely roi's to 5 body lengths
% 1-4)
roi1 = q1 & se & BLidx; % quadrant 1, southeast
roi2 = q2 & sw & BLidx; % quadrant 2, southwest
roi3 = q3 & ne & BLidx; % quadrant 3, northeast
roi4 = q4 & nw & BLidx; % quadrant 4, northwest

% Gray/maybe rules
r.i = theta < 90;
r.j = theta > 90;
deg = [45, 20, 10];

    % 5) quadrant 1, southwest
    roi5 = [];
    for i = 1:3 % each angle
        % body center is within appropriate BL, body angle is between [deg] and 180, and inside quadrant 1
        dummy = abs(mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q1 & sw;
        roi5(:,i) = dummy;
    end
    roi5 = any(roi5,2);
    
    % 6) quadrant 2, southeast
    roi6 = [];
    for i = 1:3
        dummy = abs(mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q2 & se;
        roi6(:,i) = dummy;
    end
    roi6 = any(roi6,2);
    
    % 7) quadrant 3, northwest
    roi7 = [];
    for i = 1:3
        dummy = abs(mX(:,body.center)) <= i*BL & theta < deg(i) & q3 & nw;
        roi7(:,i) = dummy;
    end
    roi7 = any(roi7,2);
    
    % 8) quadrant 4, northeast
    roi8 = [];
    for i = 1:3
        dummy = abs(mX(:,body.center)) <= i*BL & theta < deg(i) & q4 & ne;
        roi8(:,i) = dummy;
    end
    roi8 = any(roi8,2);
    
    % 9) quadrant 1, northeast
    roi9 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q1 & ne;
        roi9(:,i) = dummy;
    end
    roi9 = any(roi9,2);
    
    % 10) quadrant 2, northwest
    roi10 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q2 & nw;
        roi10(:,i) = dummy;
    end
    roi10 = any(roi10,2);
    
    % 11) quadrant 3, southeast
    roi11 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q3 & se;
        roi11(:,i) = dummy;
    end
    roi11 = any(roi11,2);
    
    % 12) quadrant 4, southwest
    roi12 = [];
    for i = 1:3
        dummy = abs(mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q4 & sw;
        roi12(:,i) = dummy;
    end
    roi12 = any(roi12,2);

% ROI's that are likely courtship positions
likely = roi1 | roi2 | roi3 | roi4;
% ROI's that are maybe courtship positions
gray = roi5 | roi6 | roi7 | roi8; % exceptions moving towards X axis
gray2 = roi9 | roi10 | roi11 | roi12; %exceptions moving towards Y axis

% Saving ROIs
position = [];
position.L1 = roi1;
position.L2 = roi2;
position.L3 = roi3;
position.L4 = roi4;
position.GX1 = roi5;
position.GX2 = roi6;
position.GX3 = roi7;
position.GX4 = roi8;
position.GY1 = roi9;
position.GY2 = roi10;
position.GY3 = roi11;
position.GY4 = roi12;


% Create variable for likely and maybe courtship positions
loc = likely | gray | gray2;
T.courtposition = loc;  % yes = green, no = red on graph vv
likelyidx = find(loc);
unlikelyidx = find(~loc);
position.likely = likely;
position.GX = gray;
position.GY = gray2;
position.all_likely = loc;
position.unlikely = ~loc;
initial_var{end+1} = 'position';

%% ANALYSIS: Identify food well and calulate distance to food

% determine if the well outlines already exist
well_file = [baseDir 'well locations.mat'];
if ~exist(well_file, 'file')    
    % pull up picture to find food well
    vidpath = [getDataPath(6, 2), parameters.date, '/', parameters.videoName, '/compiled_video_1.avi'];
    % vidpath = '/Volumes/OnTheGoData/Courtship Videos/09.26.2024/Berlin_courtship_F_LRR_caviar_ramp/compiled_video_1.avi';
    movieInfo = VideoReader(vidpath); %read in video
    demoImg = (read(movieInfo,T.vidFrame(1)));
    img = imadjust(demoImg,[72/255, 215/255]); % adjust the contrast of the image
    
    % save food well location
    txt = {'12', '3', '6', '9'};
    well = struct; % initialize the new well structure
    fig = getfig('');
        imshow(img)
        for i = 1:4 
            h = warndlg(['Outline the ' txt{i} ' oclock well']);
            uiwait(h)
            roi = drawcircle; % manually add in the circle over the food well
            well.radius(i) = roi.Radius;
            well.center(i,:) = roi.Center;
        end
    save_figure(fig,[figDir 'well outlines'],'-pdf',0,1, '-r100');

    well.R = mean(well.radius)*pix2mm;

    % select the well with the food
    fig = getfig('');
        imshow(img); hold on
        for i = 1:4
            viscircles(well.center(i,:), well.radius(i),'color',Color('white'));
            % drawcircle('Center',[well.center(i,:)],'Radius',well.radius(i),'StripeColor','red','Color',Color('grey'));
        end
        title('Click the food well','FontSize',18)
        [xi, yi] = crosshairs(1,{'black','black','yellow','yellow'});      % get a point
        % find the well with the selected data point
        % distance between each well center and the selected point
        [~, index] = min((sqrt((xi-well.center(:,1)).^2 + (yi-well.center(:,2)).^2)).*pix2mm);
        
        % plot the selected point
        scatter(xi,yi,60,"white",'filled')
        scatter(xi,yi,30,"red",'filled')
        % change the color of the well outline to highlight selection
        viscircles(well.center(index,:), well.radius(index),'color',Color('red'));
    if strcmp(questdlg('okay food well selection?'),'Yes')
        well.food_idx = index;
        well.food = well.center(index,:);
        close(fig)
    else
        warndlg('Rerun the food well identification')
        return
    end
    % SAVE DATA?
    if strcmp(questdlg('save well locations?'),'Yes')
        save(well_file,'well')
    else 
    end
else
    load(well_file,'well')
    disp('Loaded prior well locations')
end
     
% Center of the arena
WC = well.center';
N = [];
x1 = WC(1,1:2:4);
y1 = WC(2,1:2:4);
x2 = WC(1,2:2:4);
y2 = WC(2,2:2:4);
[xi,yi] = polyxpoly(x1,y1,x2,y2);
well.center(5,:) = [xi,yi];

% calculate distance to food (from fly head)
x1 = m.pos(:,1,1); % x location for male center
y1 = m.pos(:,1,2);
x2 = f.pos(:,1,1);
y2 = f.pos(:,1,2);
c1 = well.center(1);
c2 = well.center(2);

m.dist2food = (sqrt((x1-c2).^2 + (y1-c2).^2)).*pix2mm;
f.dist2food = (sqrt((c1-x2).^2 + (c1-y2).^2)).*pix2mm;
T.dist2food = [m.dist2food, f.dist2food];

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

%% ANALYSIS: Calculate male wing extension

clearvars('-except',initial_var{:})

Lwing = [];
Rwing = [];
% Determine which positions require which wing to be extended
L_items = {'L1', 'L4', 'GX1', 'GX4', 'GY2', 'GY3'};
R_items = {'L2', 'L3', 'GX2', 'GX3', 'GY1', 'GY4'};
% Identify if male is in an appropriate position for each wing direction across each item
for i = 1:length(L_items)
    Lwing = [Lwing, position.(L_items{i})];
    Rwing = [Rwing, position.(R_items{i})];
end
% Condense to identify if male is in any of the appropriate positions
Lwing = any(Lwing,2);
Rwing = any(Rwing,2);

% Pull wing angles equal or greater than extension minimum for L and R
wa_cutoff = 50; % minimum wing extension angle for courtship
wing_ext = (Lwing & (m.wing.angle(:,1) >= wa_cutoff)) | (Rwing & (m.wing.angle(:,2) >= wa_cutoff)); % wing must be facing the female fly
% Each value subtracted by the value before it (1 = ext starts, -1 = ext stops, 0 = no state change)
a = diff(wing_ext); 
% Add the first extension value to the list to account for the starting condition
b = [wing_ext(1); a]; 
% Locations in wing_ext where extension period starts/end
ext_start = find(b == 1); 
ext_stop = find(b == -1);
% If wing ext doesn't stop by end, add stop location at end of ext_stop (loc = length of experiment value)
if wing_ext(end)
    ext_stop(end + 1) = length(time);
end
% Calculate the length of each wing ext bout
ext_dur = ext_stop - ext_start;
% Find where wing ext lasts longer than 1sec
dur_loc = find(ext_dur > fps);

% Create new courtship matrix with only true wing ext for bouts longer than 1sec
mt = false(size(time));
for i = 1:length(dur_loc)
    ii = dur_loc(i);
    mt(ext_start(ii):ext_stop(ii)) = true;
end
T.wing_ext = mt;
T.wing_ext_all = wing_ext;

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

%% ANALYSIS: Fly turning 

clearvars('-except',initial_var{:})
for sex = 1:2
    % extract the x and y head and center positions to calculcate the slope of the body line
    x = data(sex).rawX(:,body.head:body.center);
    y = data(sex).rawY(:,body.head:body.center);
    
    % zero the flys center (actually unneccessary) 
    X = x - x(:,2);
    Y = y - y(:,2);
    % slope for each point in time
    slope = (Y(:,2)-Y(:,1))./(X(:,2)-X(:,1));
    m1 = slope(1:end-1);   % 'past' heading
    m2 = slope(2:end);      % 'current' heading
    
    % calculate the angle between the two slopes 
    theta = atan((m1-m2)./(1+(m1.*m2)));
    theta = rad2deg(theta);
    data(sex).turning = [nan; theta].*(fps);
end

%% ANALYSIS: Extract sleep data

clearvars('-except',initial_var{:})

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

%%