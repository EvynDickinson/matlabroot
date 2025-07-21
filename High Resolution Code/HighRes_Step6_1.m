
% TODO: 
% R =  25.6;
% Adjust the outer ring occupancy to represent the accessible space in the arena

%% Load grouped data to compare across flies and ramps

clear; close all; clc
baseDir = getDataPath(6,0);
baseFolder = [baseDir,'Trial Data/'];

% Find files that can be run
[excelfile, Excel, xlFile] = load_HighResExperiments;

if strcmp(questdlg('Load premade data structure?'),'Yes')
    FileNames = dir([baseDir 'grouped/']);
    FileNames = {FileNames(3:end).name};
    fileIdx = listdlg('ListString', FileNames,'ListSize',[350,450],'promptstring', 'Select data group to load');
    if isempty(fileIdx)
        disp('No trials selected, build a new one')
        return
    else
        load([baseDir 'grouped/' FileNames{fileIdx} '/GroupData.mat'])
        % set up the data paths for the current computer configuration: 
        baseDir = getDataPath(6,0,'Select the location for saving new data and figures:');
        baseFolder = [baseDir,'Trial Data/'];
        % make or assign the base grouped data folder
        groupDir = createFolder([baseDir 'grouped/' FileNames{fileIdx} '/']);
        % make or assign the figure folder in the grouped data folder
        figDir = createFolder([groupDir 'Figures/']);

        disp('Data structure finished loading')
        return % downloading the data
    end
else
    % Find trials that have nothing in their Proofed Column & have been tracked: 
    switch questdlg('Load options from excel?')
        case 'Yes'
                loc = cellfun(@isnan,excelfile(2:end,Excel.groupready));
                loc = ~loc;
                rownums = find(loc)+1; 
                eligible_files = excelfile([false;loc],[Excel.date, Excel.expID, Excel.ramp]);
                FileNames = format_eligible_files(eligible_files);
                
                fileIdx = listdlg('ListString', FileNames,'ListSize',[350,450],'promptstring', 'Select trial to process');
                if isempty(fileIdx)
                    disp('No trials selected')
                    return
                else
                    selectedFiles = {};
                    for i = 1:length(fileIdx)
                        selectedFiles{i} = [eligible_files{fileIdx(i),1} '_' eligible_files{fileIdx(i),2}];
                    end
                end
    
                
                % pull the list of dates and arenas to be loaded
                % TODO: need some way to pair these by group and ramp 
                % dateStr = eligible_files(fileIdx,1);
                % trialStr = eligible_files(fileIdx,2); 
    
                % for now, load all the trials as 'independent' trials
        case 'No'
              fileList = dir(baseFolder);
              fileIdx = listdlg('ListString', {fileList(4:end).name},'ListSize',[350,450],'promptstring', 'Select trial to process');
                if isempty(fileIdx)
                    disp('No trials selected')
                    return
                end
                fileIdx = fileIdx + 3; % offset to account for the skipped presentation of the '.' files above
                selectedFiles = {fileList(fileIdx).name};
    
        case 'Cancel'
            return
    end

    field_list = {'T','data','f', 'fX', 'fY', 'fps', 'm', 'mX', 'mY', 'nvids', 'parameters', 'pix2mm', 'position', 'tRate', 'time', 'well'};
    fly = struct;
    fill_idx = 1;
    for i = 1:length(fileIdx)
        try dummy = load([baseFolder selectedFiles{i} '/post-5.1.1 data.mat']);
            fly(fill_idx).name = selectedFiles{i};
            for ii = 1:length(field_list)
                fly(fill_idx).(field_list{ii}) = dummy.(field_list{ii});
            end
            fill_idx = fill_idx +1;
        catch; disp(['Missing data: ' selectedFiles{i}])
        end
    end
    
    F = 2;
    M = 1;
    body = dummy.body;
    
    blkbgd = false;
    
    num.trials = length(fly);
    
    clear i ii dummy rownums ans fileIdx field_list FileNames loc trialStr dateStr fileList selectedFiles eligible_files
    initial_var = who; 
    initial_var{end+1} = 'initial_var';
    initial_var{end+1} = 'data';

end

% make or assign the figure folder in the grouped data folder
% temporary option:
switch questdlg('Will this data be saved to an existing data structure?')
    case 'Yes'
        fileList = dir([baseDir 'grouped/']);
        fileIdx = listdlg('ListString', {fileList(3:end).name},'ListSize',[350,450],'promptstring', 'Select group:');
        if isempty(fileIdx)
            disp('No group selected')
            return
        end
        fileIdx = fileIdx + 2; % offset to account for the skipped presentation of the '.' files above
        groupName = fileList(fileIdx).name;

    case 'No'
        str = inputdlg('Write group name:');
        groupName = str{:};
    case 'Cancel'
        return
end

% make or assign the base grouped data folder
groupDir = createFolder([baseDir 'grouped/' groupName '/']);
% make or assign the figure folder in the grouped data folder
figDir = createFolder([groupDir 'Figures/']);

initial_var{end+1} = 'groupDir';
initial_var{end+1} = 'figDir';
initial_var{end+1} = 'groupName';
initial_var{end+1} = 'conversion';

conversion = getConversion; % this pulls in the standard pix2mm values for trials

fig_type = '-png';
initial_var{end+1} = 'fig_type';

clearvars('-except',initial_var{:})

%% Update the data structures to have appropriate distances for the arena 
% Fixes any old data that uses incorrect distance measurments from faulty
% pix2mm conversion AS OF May 12, 2025 ESD

% TODO  : 5.12.25
% check if there is a 'plate update' flag for a given trial ... if not, it
% will process the data and replace the current distance related
% information with the distance updated information 


% figDir = '/Users/evyndickinson/Desktop/Berlin LTS caviar/';
% % save([figDir 'GroupData.mat'])

% baseDir = '/Users/evyndickinson/Documents/Jeanne Lab/Courtship Videos/';
% baseFolder  = [baseDir 'grouped/Berlin F LRR 25-17 caviar MF/'];
% figDir = [baseFolder 'Figures/'];

%% Find the temperature alignment across trials...
% TODO -- slow down the figure progression so that you can see all of them
% and approve etc.
[foreColor, ~] = formattingColors(blkbgd); %get background colors

alignmentDir = createFolder([figDir, 'Alignment Figures/']);

fig = getfig('',1); hold on
    for i = 1:num.trials
        x = fly(i).time;
        y = fly(i).T.temperature;
        y = smooth(y,90,'moving');
        plot(x,y)
        z = [fly(i).tRate(:).idx];
        v_line(x(z(2:2:end)),'r')
    end
    xlabel('time (min)')
    ylabel('temp (\circC)')
    formatFig(fig);
    title([groupName ' | n = ' num2str(num.trials)])
save_figure(fig, [alignmentDir 'raw temp alignment'],fig_type);

% update for the correct temp protocol names: 
temp_protocols = [];
for i = 1:num.trials
    try temp_protocols{i} = fly(i).parameters.protocol;
    catch
        if contains(fly(i).parameters.expID,'courtship_F_LRR')
            temp_protocols{i} =  'courtship_F_LRR_25-17';
        else 
            temp_protocols{i} =  ' ';
        end
    end
end
disp(temp_protocols')


% create a new matrix to look at trial alignment times
switch fly(1).parameters.protocol
    case 'high_res_LTS_35-15'
        idx = [1 3 5 7 9 10]; % LTS alignment points
        midpoint_idx = 4; % bottom trough of full sweep
    case {'courtship_F_LRR_25-17','high_res_F_LRR_25-17'}
        idx = [1 3 5 7 8]; % F LRR alignment points (beginning of each temp region)
        midpoint_idx = 3; % center of ramp
end

% get a matrix of all the transition points between the temp roi types
% (e.g. cooling vs warming vs cooling again) 
MT = nan(num.trials, length(idx));
for i = 1:num.trials
    z = [fly(i).tRate(:).idx]; % concatenate all the ROI transition points across the sections
    MT(i,:) = z(idx);
end

% find the ramp durations to make sure that everything is equally correct,
% even if not perfectly time-synced
xtick_lab = {};
buff = 0.3;
y = diff(MT,1,2);
rampDur = mean(y(:,2:3),"all"); % mean ramp duration for our 
y = y./(fly(1).fps*60);

fig = figure; 
    hold on
    for i = 1:length(idx)-1
        x = i*ones(num.trials,1);
        scatter(x, y(:,i),35, fly(1).tRate(i).color,'filled', 'XJitter','density')
        plot([i-buff, i+buff], [mean(y(:,i)), mean(y(:,i))],'color', foreColor)
        xtick_lab{i} = fly(1).tRate(i).name;
    end
    ylabel('period duration (min)')
    set(gca, 'xtick', 1:length(xtick_lab), 'XTickLabel', xtick_lab)
    title([groupName ' | n = ' num2str(num.trials)])
    formatFig(fig, blkbgd);
save_figure(fig, [alignmentDir 'raw section durations'],fig_type);

% Shift trials to align to the ramp trough:
% preallocate empty space for all the variables that need to be time-shifted across trials
max_len = size(fly(1).T,1); % set up a buffered time (assume all trials are the same length) 
optA =  nan(max_len,num.trials); % single data column for each trial
optB = nan(max_len,2,num.trials); % two data columns per trial
categories = {'temperature', 'frame', 'IFD', 'cooling', 'warming', 'hold', 'wing_ext', 'wing_ext_all',...
              'court_chase', 'chase_all', 'circling_all', 'circling_1sec', 'CI'}; % from fly.T
sub_cat = {'eccentricity', 'OutterRing', 'foodQuad', 'foodcircle', 'turning'}; % from fly.data


% DATA STRUCT is time aligned across trials!!!! AKA use this data moving
% forward in the analyses
data = struct; 
for i = 1:length(categories)
    data.(categories{i}) = optA;
end
[data.dist2food, data.FlyOnFood] = deal(optB);
for i = 1:length(sub_cat)
    data.(sub_cat{i}) = optB;
end

% fill in data and shift to appropriate place...
middlePoint = median(MT(:,midpoint_idx)); % this is the new 'center point' in time for the aligned trials to minimize lost data across trials
middlePoint = int32(middlePoint); % make sure the number is indexable

sexlist = {'m', 'f'};
MT_pos = nan([max_len, 5, num.trials]);
[data.x_loc(1).pos, data.x_loc(2).pos, data.x_loc(1).pos, data.y_loc(2).pos]  = deal(MT_pos);

for i = 1:num.trials
    curr_center = MT(i,midpoint_idx); %current index of first temp ramp trough for this specific trial
    shift_num = middlePoint-curr_center; %index offset value
    if shift_num<0 %shift data points to the left 
        curr_idx = abs(shift_num)+1:max_len; % index of the current data to take
        new_idx = 1:max_len-abs(shift_num);
    elseif shift_num>0 % shift data points to the right
        new_idx = shift_num+1:max_len;
        curr_idx = 1:max_len-shift_num;
    elseif shift_num==0
        new_idx = 1:max_len;
        curr_idx = 1:max_len;
    else 
        disp('unknown error')
        return
    end

    loc = false(max_len, 1);
    loc(curr_idx) = true;
    % update the single data line trials from fly.T
    for j = 1:length(categories)
        data.(categories{j})(new_idx,i) = fly(i).T.(categories{j})(curr_idx);
    end
    % update the sex specific trial data lines from fly.data
    for sex = 1:2
        for j = 1:length(sub_cat)
            data.(sub_cat{j})(new_idx,sex,i) = fly(i).data(sex).(sub_cat{j})(curr_idx);
        end
        data.dist2food(new_idx,sex,i) = fly(i).T.dist2food(curr_idx,sex);
        data.FlyOnFood(new_idx,sex,i) = fly(i).T.FlyOnFood(curr_idx,sex);
        % update body position locations to time-shifted alignment
        data.x_loc(sex).pos(new_idx,:,i) = fly(i).(sexlist{sex}).pos(curr_idx,:,1);
        data.y_loc(sex).pos(new_idx,:,i) = fly(i).(sexlist{sex}).pos(curr_idx,:,2);
    end
    data.sleep(new_idx,M,i) = fly(i).m.sleep(curr_idx);
    data.sleep(new_idx,F,i) = fly(i).f.sleep(curr_idx);
end

data.time = fly(1).time;
data.color(1,:) = fly(1).data(1).color; % male color
data.color(2,:) = fly(1).data(2).color; % female color
data.middlepoint = middlePoint;

% TODO: adjust this for protocols that are not mirror image parallel
switch temp_protocols{1}

    case 'high_res_LTS_35-15'
        % save new indexes for plotting regions: 
        a = (diff(data.cooling)); % use this to look for transition periods in and out of cooling        
        [r,~] = find(a==1); %  find the time points where cooling begins
        c1 = median(r(1:2:end)); % cooling start 1
        c2 = median(r(2:2:end)); % cooling start 2 
        c3 =  size(a,1); % end of cooling (end of experiment) 

        a = (diff(data.warming));
        [r,~] = find(a==1); % when warming section two starts
        h1 = median(r(1:2:end)); % warming start 1
        h2 = median(r(2:2:end)); % warming start 2 
        
        % find compile the start and stops for the cooling or warming periods
        data.warming_idx = [h1,c1-1; h2,c2-1];
        data.cooling_idx = [c1,h2-1; c2,c3];
     
    case 'courtship_F_LRR_25-17' % TODO : 7/10 update this to work with the new ramp strucutres... 
        % save new indexes for plotting regions: 
        a = (diff(data.cooling));
        data.cooling_idx_all(:,1) = find(a==1);  % start of cooling
        data.cooling_idx_all(:,2) = find(a==-1); % end of cooling (middlepoint)
        a = (diff(data.warming));
        data.warming_idx_all(:,1) = find(a==1);  % start of warming
        data.warming_idx_all(:,2) = find(a==-1); % end of warming (middlepoint)
        % universal timings: 
        data.cooling_idx = int32(median(data.cooling_idx_all));
        data.warming_idx = median(data.warming_idx_all);
        data.warming_idx(1) = data.warming_idx(1)+1; % offset the start of warming by 1 frame
        data.warming_idx = int32(data.warming_idx);
end 


% test to see how the alignment worked for the heating and cooling periods
fig = getfig('',0);
    subplot(1,2,1)
    plot(data.cooling)
    v_line(middlePoint,'red', '-',2)
    v_line(middlePoint,'yellow', '--',2)
    title('cooling')
    xlabel('frame')
    set(gca, 'ytick', [0,1],'yticklabel', {'Off', 'On'})
    ylim([-0.1,1.1])
    subplot(1,2,2)
    plot(data.warming)
    v_line(middlePoint,'red', '-',2)
    v_line(middlePoint,'yellow', '--',2)
    xlabel('frame')
    set(gca, 'ytick', [0,1],'yticklabel', {'Off', 'On'})
    ylim([-0.1,1.1])
    title('warming')
    formatFig(fig, blkbgd,[1,2]);
save_figure(fig, [alignmentDir 'cooling and warming alignment'],fig_type);


% show the alignment across temperature
LW  = 2;
fig = getfig('', 1);
    plot(data.temperature)
    y = rangeLine(fig,5,true);
    hold on
    % plot warming regions from the warming indx
    for i = 1:size(data.warming_idx,1)
        plot([data.warming_idx(i,1),data.warming_idx(i,2)],[y,y],'color', 'r','linewidth', LW)
    end
    % plot cooling regions from the cooling indx
    for i = 1:size(data.cooling_idx,1)
        plot([data.cooling_idx(i,1),data.cooling_idx(i,2)],[y,y],'color', Color('dodgerblue'),'linewidth', LW)
    end
    v_line((data.cooling_idx(:)),'grey', '--',0.5)
    v_line((data.warming_idx(:)),'grey', '--',0.5)
    xlabel('time')
    ylabel('temperature (\circC)')
    formatFig(fig,blkbgd);
    title([groupName ' | n = ' num2str(num.trials)])
save_figure(fig, [alignmentDir 'final temperature alignment'],fig_type);

disp('next:')

%% Create temp bin groups and a universal temperature ramp:
clearvars('-except',initial_var{:})

% universal temperature profile: 
data.temp = smooth(mean(data.temperature,2, 'omitnan'),5*30,'moving');

% create logical indexes for heat/cool at specific temperature bins
temp_bins = floor(min(data.temp)):0.5:ceil(max(data.temp)); % 0.5 deg temperature bins
idx = discretize(data.temp, temp_bins);
[c_idx, h_idx] = deal(false(size(data.temp)));
for i = 1:size(data.warming_idx,1)
    h_idx(data.warming_idx(i,1):data.warming_idx(i,2)) = true;
end
for i = 1:size(data.cooling_idx,1)
    c_idx(data.cooling_idx(i,1):data.cooling_idx(i,2)) = true;
end

tempbin = [];
tempbin.c_idx = c_idx; % times when cooling is happening
tempbin.h_idx = h_idx; % times when warming is happening
tempbin.temps = temp_bins;
[tempbin.all, tempbin.cooling, tempbin.warming] = deal(false(length(data.temp),length(temp_bins)));

for i = 1:length(temp_bins)
    tempbin.all(:,i) = (idx==i); % both warming and cooling
    tempbin.cooling(:,i) = (idx==i) & c_idx; % cooling
    tempbin.warming(:,i) = (idx==i) & h_idx; % cooling
end
        
data.tempbin = tempbin;
disp('next:')

%% Update the distance to food metrics that were incorrect before 1.23.25

% conversion = getConversion;
pix2mm = conversion(4).pix2mm;

for i = 1:num.trials
    % calculate distance to food (from fly head)
    x1 = fly(i).m.pos(:,body.center,1);  % x location for male center
    y1 = fly(i).m.pos(:,body.center,2); % y location for male center
    x2 = fly(i).f.pos(:,body.center,1); % x location for female center
    y2 = fly(i).f.pos(:,body.center,2); % y location for female center
    c1 = fly(i).well.food(1); % food well center
    c2 = fly(i).well.food(2); % food well center
    
    M_dist  = (sqrt((x1-c1).^2 + (y1-c2).^2))./fly(i).pix2mm;
    F_dist = (sqrt((x2-c1).^2 + (y2-c2).^2))./fly(i).pix2mm;
   
    data.dist2food(:,M,i) = M_dist;
    data.dist2food(:,F,i) = F_dist;

    fly(i).m.dist2food = M_dist;
    fly(i).f.dist2food = F_dist;
    fly(i).T.dist2food = [M_dist, F_dist];

end
disp('updated distance to food data structures')

%% Analysis: create spatial regions screens for ROIs of interest in each trial arena
% TODO: working here 7/16

%% ANALYSIS: Create ROImasks -- logical masks for all the regions in the arena for each trial
clearvars('-except',initial_var{:})
% Add a new structure that contains the masks for each of the regions of
% interest for the flies in this experiment which can be used to filter
% other properties later on...like speed in the outer ring etc.

ct = 4; % condition type for the current high resolution experiment trials

% Initialize empty mask structure:
initial_var{end+1} = 'ROImask';
initial_var{end+1} = 'quadOrder';
ROImask = struct; 
fields = {'all','ring', 'inner75'}; 
food_fields = {'fullquad','innerquad','quadring','circle10','circle7','circle5'};
quadOrder = {'food','right','opp','left'};
for trial = 1:num.trials
    % universal regions:
    for i = 1:length(fields)
        ROImask.(fields{i})(trial).m = [];
    end
    % quadrant regions:
    for i = 1:length(food_fields)
        for q = 1:4
            ROImask.(food_fields{i}).(quadOrder{q})(trial).m = [];
        end
    end
end

exp = 1;
% find locations of flies within each region of the arena
for trial = 1:num.trials
    % pull parameters for this trial:
    center = fly(trial).well.center(5,:);
    % M & F fly body positions
    x1 = data.x_loc(M).pos(:,body.center,trial); % x location for male center
    y1 = data.y_loc(M).pos(:,body.center,trial);  % y location for male center
    x2 = data.x_loc(F).pos(:,body.center,trial); % x location for female center
    y2 = data.y_loc(F).pos(:,body.center,trial);  % y location for female center

    % matrix with both flies: male first then female -- but this way we can run the processing in parallel
    x_loc = [x1,x2]; 
    y_loc = [y1, y2]; 
    
    % ct = data(exp).con_type(trial); % experiment lens configuration
    pix2mm = conversion(ct).pix2mm;
    R = conversion(ct).R;
    circle75 = conversion(ct).circle75; % defines the distance to the inside edge of the outer ring
    circle10 = conversion(ct).circle10;
    circle7 = conversion(ct).circle7;
    circle5 = conversion(ct).circle5;
    foodWell = fly(trial).well.food_idx;  % which well contains the food TODO: update this to account for no food trials...
    % foodWell = data(exp).T.foodLoc(trial);

    % Adjust the X and Y coordinates relative to new center of (0,0)
    adjustedX = x_loc - center(1); 
    adjustedY = y_loc - center(2);

    % Pull distance measures:  
    D = hypot(adjustedX, adjustedY)./pix2mm; % distance to center in mm

    % add screen here for regions in the temp protocol that are to be excluded:
    % tp = getTempTurnPoints(data(exp).T.TempProtocol{trial});
    % all_points = [tp.DownROI,tp.UpROI,tp.HoldROI];
    % loc = false(size(x_loc));
    % loc(all_points,:) = true;
    % if tp.holdexp
        loc = true(size(x1)); % start with all temporal locations being allowable
    % end

    % Define base location logicals
    fly_loc = ~isnan(x_loc) & (D <= R) & loc; % data points that are valid flies

    innerQuad = (D < circle75) & fly_loc; %logical of all points that have a fly within the inner circle
    adjustedX(~fly_loc) = nan;
    adjustedY(~fly_loc) = nan;
    Q = findQuadLocation(adjustedX,adjustedY);

    % sum of total flies in the three large regions
    nTotal = sum(fly_loc,2); % total across time
    nTotalInnerROI = sum(innerQuad,2); % this gives the bottom of the 'inner distribution' fraction
    nTotalOuterROI = nTotal-nTotalInnerROI; % this gives the bottom of the 'outer distribution' fraction
    ROImask(exp).all(trial).nflies = nTotal;
    ROImask(exp).ring(trial).nflies = nTotalOuterROI;
    ROImask(exp).inner75(trial).nflies = nTotalInnerROI;

    % ------- Find the flies in the outer ring and inner region -------
    % Define quadrant masks based on the new center & excluding flies that are in the outer ring:
    % food_fields = {'circle10','circle7','circle5'};
    ROImask(exp).all(trial).m = fly_loc; % ALL FLIES IN FULL ARENA (with bounds)
    ROImask(exp).ring(trial).m = fly_loc & ~innerQuad; % ALL FLIES WITHIN THE OUTER RING
    ROImask(exp).inner75(trial).m = fly_loc & innerQuad; % ALL FLIES IN INNER CIRCLE

    % Find the food quadrant (find location with the food well coordinates included)
    adjusted_wx = fly(trial).well.center(foodWell,1) - center(1); % adjusted well position
    adjusted_wy = fly(trial).well.center(foodWell,2) - center(2); 
    well_locations = findQuadLocation(adjusted_wx,adjusted_wy);
    quad_loc = find([well_locations(:).Mask]); % quadrant idx that has food

    % opposition Matrix: orientation of the quadrants such that loc 1 is the food quad and 
    % then it goes quad right, opposite food quad, and finally quad left of the food quad
    switch quad_loc 
        case 1
            opLoc = [1 4 3 2]; 
        case 2
            opLoc = [2 1 4 3];
        case 3
            opLoc = [3 2 1 4];
        case 4
            opLoc = [4 3 2 1];
    end

    % Find the quadrant that each well belongs to
    well_opt_x = fly(trial).well.center(1:4,1) - center(1); % adjusted well position
    well_opt_y = fly(trial).well.center(1:4,2) - center(2);
    well_quads = findQuadLocation(well_opt_x,well_opt_y);

    for i = 1:4 % for quadrant type (food, R, opp, L)
        idx = opLoc(i); % rearrange order based on trial food location
        well_idx = find(well_quads(idx).Mask);
        % well_idx = find(well_quad_loc==idx); % index of the well that falls into this quadrant of interest
        ROImask(exp).fullquad.(quadOrder{i})(trial).m = Q(idx).Mask; % full quadrant
        ROImask(exp).innerquad.(quadOrder{i})(trial).m = Q(idx).Mask & innerQuad; % quadrant inside inner circle
        ROImask(exp).quadring.(quadOrder{i})(trial).m = Q(idx).Mask & ~innerQuad; % quadrant portion of outer ring

        % well circle related distances
        dX = adjustedX - well_opt_x(well_idx);% well_idx
        dY = adjustedY - well_opt_y(well_idx);% well_idx
        well_D = hypot(dX,dY)./pix2mm;
        ROImask(exp).circle10.(quadOrder{(i)})(trial).m = well_D <= circle10; % within 10% area of well
        ROImask(exp).circle7.(quadOrder{(i)})(trial).m = well_D <= circle7; % within 7% area of well
        ROImask(exp).circle5.(quadOrder{(i)})(trial).m = well_D <= circle5; % within 5% area of well
    end
end


%% ANALYSIS:  Secondary organization of the ROImask data grouped by sex and not trial
% TODO: 7.18 create a new structure that has the ROImask data grouped by
% sex and not trial that can more easily be used with the 'data' structure
clearvars('-except',initial_var{:})
initial_var{end+1} = 'flyROImask';
% initiate empty structures
flyROImask = struct;
fields = {'all','ring', 'inner75'}; 
food_fields = {'fullquad','innerquad','quadring','circle10','circle7','circle5'};
quadOrder = {'food','right','opp','left'};
MT = nan([size(ROImask.all(1).m,1), num.trials]);
% blanks for the full arena regions
for i = 1:length(fields)
    flyROImask.(fields{i})(M).m = MT;
    flyROImask.(fields{i})(F).m = MT;
end
% blanks for the areas with four regions
for i = 1:length(food_fields)
    for r = 1:length(quadOrder)
        flyROImask.(food_fields{i})(M).(quadOrder{r}).m = MT;
        flyROImask.(food_fields{i})(F).(quadOrder{r}).m = MT;
    end
end

% Pull and restructure data from the original ROImask (by trial) structure:
for trial = 1:num.trials
    for i = 1:length(fields)
        flyROImask.(fields{i})(M).m(:,trial) = ROImask.(fields{i})(trial).m(:,M);
        flyROImask.(fields{i})(F).m(:,trial) = ROImask.(fields{i})(trial).m(:,F);
    end
    
    for i = 1:length(food_fields)
        for r = 1:length(quadOrder)
            flyROImask.(food_fields{i})(M).(quadOrder{r}).m(:,trial) = ROImask.(food_fields{i}).(quadOrder{r})(trial).m(:,M);
            flyROImask.(food_fields{i})(F).(quadOrder{r}).m(:,trial) = ROImask.(food_fields{i}).(quadOrder{r})(trial).m(:,F);
        end
    end
end

% Get the summary numbers from each region across the flies: 
for i = 1:length(fields)
    flyROImask.(fields{i})(M).nflies = sum(flyROImask.(fields{i})(M).m,2,'omitnan');
    flyROImask.(fields{i})(F).nflies = sum(flyROImask.(fields{i})(F).m,2,'omitnan');
end

% demo figure: 
% CRAZY FIGURE (keeeeeep tattoo fodder)
% fig = getfig('',1);
%     plot(flyROImask.inner75.nflies,'color',Color('dodgerblue'))
%     hold on
%     plot(flyROImask.ring.nflies,'color',Color('orange'))
% formatFig(fig, false);
% set(gca,'xcolor', 'none', 'ycolor', 'none')
% save_figure(fig, [figDir 'tattoo fodder'],'-pdf')

% CRAZY FIGURE (keeeeeep tattoo fodder)
% fig = getfig('',1);
%     plot(flyROImask.inner75(F).nflies,'color',Color('dodgerblue'))
%     hold on
%     plot(flyROImask.ring(F).nflies,'color',Color('orange'))
% formatFig(fig, true);
% set(gca,'xcolor', 'none', 'ycolor', 'none')
% save_figure(fig, [figDir 'tattoo fodder'],'-pdf')

%% ANALYSIS: Find occupancy in the different regions 
clearvars('-except',initial_var{:})
% TODO (7/17): set this up to have a pooled group across all the male
% flies and one for all the female flies i.e. not combined across the sexes
regionList = {'fullquad','innerquad', 'quadring', 'circle10', 'circle7', 'circle5'};


for sex = 1:2
    nFull = flyROImask.all(sex).nflies; % total fly count in arena
    nInner = flyROImask.inner75(sex).nflies; % fly count in the inner 75% region
    nRing = flyROImask.ring(sex).nflies; % fly count in the outer ring
    
    % ring & inner -- each get compared to the full arena compliment:
    n = getPercentFlies(flyROImask.ring(sex).m,nFull);
    flyROImask.ring(sex).avg = n;
    
    n = getPercentFlies(flyROImask.inner75(sex).m,nFull);
    flyROImask.inner75(sex).avg = n;

    % for each of the quadrant-related regions:
    for rr = 1:length(regionList)
        for q = 1:4 % each of the quadrants
            temp = flyROImask.(regionList{rr})(sex).(quadOrder{q}).m;
            a = mean(temp,2,'omitnan');
            flyROImask.(regionList{rr})(sex).(quadOrder{q}).avg = getPercentFlies(temp,nFull);

            if strcmp(regionList{rr},'innerquad')
                flyROImask.(regionList{rr})(sex).(quadOrder{q}).partial_avg = getPercentFlies(temp,nInner);
            elseif strcmp(regionList{rr},'quadring')
                flyROImask.(regionList{rr})(sex).(quadOrder{q}).partial_avg = getPercentFlies(temp,nRing);
            end
        end
    end
end


% TODO (7.18): update this to work with static temperature holds...
all_regions = [regionList, 'ring', 'inner75'];
for sex = 1:2
    
    temps = data.tempbin.temps; % binned temp groups
    nTemp = length(temps);
    cIDX = data.tempbin.cooling; % logical locations for each temp bin across time
    hIDX = data.tempbin.warming;

    for i = 1:length(all_regions)

        switch all_regions{i}
            case {'ring','inner75'}
                nQ = 1; % number of quadrants ...
                nP = 1; % number of partial region percentages to check
            case {'innerquad','quadring'}
                nQ = 4;
                nP = 2;
            case {'fullquad','circle10','circle7','circle5'}
                nQ = 4;
                nP = 1;
        end
        for p = 1:nP % for the number of full or partial percentages to compare
            for q = 1:nQ % each of the quadrants
                % pull data for the right quadrant (if there are quadrants)
                if nQ>1
                    baseY = flyROImask.(all_regions{i})(sex).(quadOrder{q});
                else 
                    baseY = flyROImask.(all_regions{i})(sex);
                end
                % pull the right type of data to run
                if p==1
                    y = baseY.avg;
                else
                    y = baseY.partial_avg;
                end
                % extract the temp and rate dependent information for each region: 
                [c_avg, h_avg, c_std, h_std] = deal(nan([nTemp,1])); 
                for tt = 1:nTemp %  for each temperature bin
                    c_loc = cIDX(:,tt); % logical vector where the temp rate and temp match
                    h_loc = hIDX(:,tt);
                    % find the avg and err within the temp bin: 
                    c_avg(tt) = mean(y(c_loc),'omitnan');
                    h_avg(tt) = mean(y(h_loc),'omitnan');
                    c_std(tt) = std(y(c_loc),'omitnan');
                    h_avg(tt) = std(y(h_loc),'omitnan');
                end
                % save back to the larger structure:                         
                switch all_regions{i}
                    case {'ring','inner75'}
                        flyROImask.(all_regions{i})(sex).increasing.avg = h_avg;
                        flyROImask.(all_regions{i})(sex).increasing.std = h_avg;
                        flyROImask.(all_regions{i})(sex).decreasing.avg = c_avg;
                        flyROImask.(all_regions{i})(sex).decreasing.std = h_std;
                        flyROImask.(all_regions{i})(sex).temps = temps;
                    case {'fullquad','innerquad','quadring','circle10','circle7','circle5'}
                        if p == 1 %aka full avg and not based on partial region occupancy
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).increasing.avg = h_avg;
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).increasing.std = h_std;
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).decreasing.avg = c_avg;
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).decreasing.std = c_std;
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).temps = temps;
                        else
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).increasing.p_avg = h_avg;
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).increasing.p_std = h_std;
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).decreasing.p_avg = c_avg;
                            flyROImask.(all_regions{i})(sex).(quadOrder{q}).decreasing.p_std = c_std;
                        end
                end
            end
        end
    end
end


%% Create new comparisons of distances for the group

% % TODO (1/28) give a demo image of the ring size for each area on an image still
% innerR = conversion(4).circle75;
% 
% for i = 1:num.trials
%     % % find distance from center for each fly center point: 
%     % xi = fly(i).well.center(5,1);
%     % yi = fly(i).well.center(5,2);
%     % % eccentricity, outer ring occupancy of the flies
%     % xm = fly(i).m.pos(:,body.center,1);
%     % ym = fly(i).m.pos(:,body.center,2);
%     % xf = fly(i).f.pos(:,body.center,1);
%     % yf = fly(i).f.pos(:,body.center,2);
%     % % distance from center of arena
%     % Dm = sqrt(((xm-xi).^2 + (ym-yi).^2)).*fly(i).pix2mm; 
%     % Df = sqrt(((xf-xi).^2 + (yf-yi).^2)).*fly(i).pix2mm; 
%     for sex = 1:2
%         D = fly(i).data(sex).eccentricity;
%         % outer ring occupancy
%         % data.OutterRing(:,sex,i) = D>=innerR; % find the locations that are greater than the inner 75%
%         fly(i).data(sex).OutterRing = D>=innerR; % find the locations that are greater than the inner 75%
%     end
% end
% 
% disp('updated outer ring to food data structures')
% 
% % TODO: need to update the behavior matrix now with these... redo and add
% this to 5.1 and reprocess the data on that front...

%% TODO: create a save data in a structure thing here so that we can save figures etc to an idea
clearvars('-except',initial_var{:})

% save the existing data (if desired): 
clearvars('-except',initial_var{:})
switch questdlg('Save data to grouped data structure?')
    case 'Yes'
        save([groupDir 'GroupData.mat'],'-v7.3');
        disp('Data saved')
end













