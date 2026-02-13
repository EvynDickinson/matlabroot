
%% ANALYSIS: Determine male position relative to female

% Adapted from 5_1 analysis to determine possible courtship positions

clearvars('-except',initial_var{:})
pix2mm = conversion(4).pix2mm; %updates 5.12.25

aggression = [];
for trial = 1:num.trials

    % % ------------------------- 1) Shift each frame to the origin for female fly ------------------------- 
    % x_offset = fly(trial).data(F).rawX(:,body.center); % x values for female center
    % y_offset = fly(trial).data(F).rawY(:,body.center);
    % mx = fly(trial).data(M).rawX-x_offset; % subtract the offset from the x values for each body point
    % my = fly(trial).data(M).rawY-y_offset;
    % fx = fly(trial).data(F).rawX-x_offset;
    % fy = fly(trial).data(F).rawY-y_offset;
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

    aggression(trial).name = fly(trial).name;

    % Calulate the angle between fly bodies (both - and +)
    theta = fly(trial).data(M).mfbodyangle; 
    test = theta<90; % when is male less than 90 deg from female
    
    % ------------------------- 3) Establish likely and unlikely courtship positions ------------------------- 
    
    % Determine fly body length (~2.5mm)
    BL = 2.5*pix2mm; % in pixels
    
    r=[];
    q = [];
    % Body position location rules
    r.a = fly(trial).mX(:,body.center) >= 0; % everything to right of y axis
    r.b = fly(trial).mX(:,body.center) <= 0; % everything to left of y axis
    r.c = fly(trial).mY(:,body.center) >= 0; % everything above x axis
    r.d = fly(trial).mY(:,body.center) <= 0; % everything below x axis
    q1 = r.b & r.c; % quadrant one occupancy
    q2 = r.a & r.c; % quadrant two occupancy
    q3 = r.d & r.b; % quadrant three occupancy
    q4 = r.a & r.d; % quadrant four occupancy
    
    % Heading direction rules
    r.e = fly(trial).mY(:,body.head) >= fly(trial).mY(:,body.center); % heading direction facing north
    r.f = fly(trial).mY(:,body.head) <= fly(trial).mY(:,body.center); % heading direction facing south
    r.g = fly(trial).mX(:,body.head) >= fly(trial).mX(:,body.center); % heading direction facing east
    r.h= fly(trial).mX(:,body.head) <= fly(trial).mX(:,body.center); % heading direction facing west
    ne = r.e & r.g; % north east direction
    nw = r.e & r.h; % north west direction
    se = r.f & r.g; % south east direction
    sw = r.f & r.h; % south west direction
    
    
    % Likely rules
    BLidx = abs(fly(trial).mX(:,body.center)) <= (5*BL) & abs(fly(trial).mY(:,body.center)) <= (5*BL); % limits likely roi's to 5 body lengths
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
            dummy = abs(fly(trial).mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q1 & sw;
            roi5(:,i) = dummy;
        end
        roi5 = any(roi5,2);
        
        % 6) quadrant 2, southeast
        roi6 = [];
        for i = 1:3
            dummy = abs(fly(trial).mX(:,body.center)) <= i*BL & theta > (180 - deg(i)) & q2 & se;
            roi6(:,i) = dummy;
        end
        roi6 = any(roi6,2);
        
        % 7) quadrant 3, northwest
        roi7 = [];
        for i = 1:3
            dummy = abs(fly(trial).mX(:,body.center)) <= i*BL & theta < deg(i) & q3 & nw;
            roi7(:,i) = dummy;
        end
        roi7 = any(roi7,2);
        
        % 8) quadrant 4, northeast
        roi8 = [];
        for i = 1:3
            dummy = abs(fly(trial).mX(:,body.center)) <= i*BL & theta < deg(i) & q4 & ne;
            roi8(:,i) = dummy;
        end
        roi8 = any(roi8,2);
        
        % 9) quadrant 1, northeast
        roi9 = [];
        for i = 1:3
            dummy = abs(fly(trial).mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q1 & ne;
            roi9(:,i) = dummy;
        end
        roi9 = any(roi9,2);
        
        % 10) quadrant 2, northwest
        roi10 = [];
        for i = 1:3
            dummy = abs(fly(trial).mY(:,body.center)) <= i*BL & theta < 90 & theta > (90 - deg(i)) & q2 & nw;
            roi10(:,i) = dummy;
        end
        roi10 = any(roi10,2);
        
        % 11) quadrant 3, southeast
        roi11 = [];
        for i = 1:3
            dummy = abs(fly(trial).mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q3 & se;
            roi11(:,i) = dummy;
        end
        roi11 = any(roi11,2);
        
        % 12) quadrant 4, southwest
        roi12 = [];
        for i = 1:3
            dummy = abs(fly(trial).mY(:,body.center)) <= i*BL & theta > 90 & theta < (90 + deg(i)) & q4 & sw;
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
    % T.courtposition = loc;  % yes = green, no = red on graph vv
    likelyidx = find(loc);
    unlikelyidx = find(~loc);
    position.likely = likely;
    position.GX = gray;
    position.GY = gray2;
    position.all_likely = loc;
    position.unlikely = ~loc;
    % Add positions where M is facing F to aggression structure
    aggression(trial).MfacingF = position;
end
initial_var{end+1} = 'aggression';


%%
foreColor = formattingColors(blkbgd); % get background colors
Mcolor = Color('dodgerblue');
Fcolor = Color('deeppink');

 % 4.1) Demo selected angles in male positions relative to female fly
        zoom = [-250,250];
        % adjust skip size based on total number of points
        displayNum = 50; 
        if sum(test)<displayNum
            displayNum = sum(test);
        end
        plotLoc = randi(sum(test), [displayNum,1]);
        
        fig = getfig('Random selection of all male positions relative to female',1,[1075 871]);
        hold on
        
        % plot male coordinates for head and body within test constraints
        x = fly(1).mX(test,[body.head,body.center]);
        y = fly(1).mY(test,[body.head,body.center]);
        % only plot every [skip value] points to reduce volume
        x = x(plotLoc,:);
        y = y(plotLoc,:);
        plot(x',y','color',Mcolor)
        scatter(x(:,1),y(:,1),15,Mcolor,"filled","^") % arrow head on male
        
        % plot female coordiates for head, body, and abdomen
        x = fly(1).fX(plotLoc,[body.head,body.center,body.abdomen]);
        y = fly(1).fY(plotLoc,[body.head,body.center,body.abdomen]);
        plot(x',y','color',Fcolor, 'LineWidth', 2)
        
        % format figure
        axis  equal square
        h_line(0,'gray',':',2)
        v_line(0,'grey',':',2)
        xlim(zoom)
        ylim(zoom)
        formatFig(fig,blkbgd);
        set(gca,'XColor','none','YColor','none')
         

        %---------------------------------------------------------------------------------------------------------

        % 4.2) Visualize likely and unlikely courtship positions 
        zoom = [-250,250];
        skip = 20;

        
        likelyidx = aggression(1).MfacingF.likely(1:skip:end);
        unlikelyidx = aggression(1).MfacingF.unlikely(1:skip:end);
        allidx = [likelyidx; unlikelyidx];
        
        fig = getfig('Likely and unlikely courtship positions',1,[1075 871]);
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
        x = fly(1).mX(unlikelyidx,[body.head,body.center]);
        y = fly(1).mY(unlikelyidx,[body.head,body.center]);
        plot(x',y','color',kolor)
        scatter(x(:,1),y(:,1),15,kolor,"filled","^")
        
        % plot the likely courtship male fly body positions
        kolor = Color('limegreen');
        x = fly(1).mX(likelyidx,[body.head,body.center]);
        y = fly(1).mY(likelyidx,[body.head,body.center]);
        plot(x',y','color',kolor)
        scatter(x(:,1),y(:,1),15,kolor,"filled","^")
        
        % plot female head, center, and abdoment
        x = fly(1).fX(1:skip:end,[body.head,body.center,body.abdomen]);
        y = fly(1).fY(1:skip:end,[body.head,body.center,body.abdomen]);
        plot(x',y','color',foreColor, 'LineWidth', 2)
        % xlim([-2000,1500]); ylim([-2000,2000])
        
        % format figure
        axis  equal square
        h_line(0,'gray',':',2)
        v_line(0,'grey',':',2)
        xlim(zoom)
        ylim(zoom)
        formatFig(fig,blkbgd);
        set(gca,'XColor','none','YColor','none')
        
        formatFig(fig,blkbgd);
        rectangle('Position',[zoom(1),zoom(1) sum(abs(zoom)) sum(abs(zoom))]...
            ,'edgecolor',foreColor,'linewidth', 1)

%% Calculate M wing extension

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