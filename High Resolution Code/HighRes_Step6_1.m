
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
optA =  nan(max_len,num.trials); %single data column for each trial
optB = nan(max_len,2,num.trials); %two data columns per trial
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
        % TODO!! 3.23.25 work on adjusting this to work for the LTS
        % protocol -- need to get the code to work dynamically for there
        % being two different warming periods and only one hold period 
        
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

% find locations of flies within each region of the arena
for trial = 1:num.trials
    % pull parameters for this trial:
    center = fly(trial).well.center(5,:);
    % pull fly center positions for the fly body location...
    


    x_loc = data(exp).data(trial).data.x_loc; % all fly X positions in the arena
    y_loc = data(exp).data(trial).data.y_loc; % all fly Y positions in the arena
    ct = data(exp).con_type(trial); % experiment lens configuration
    pix2mm = conversion(ct).pix2mm;
    R = conversion(ct).R;
    circle75 = conversion(ct).circle75; % defines the distance to the inside edge of the outer ring
    circle10 = conversion(ct).circle10;
    circle7 = conversion(ct).circle7;
    circle5 = conversion(ct).circle5;
    foodWell = data(exp).T.foodLoc(trial);

    % Adjust the X and Y coordinates relative to new center of (0,0)
    adjustedX = x_loc - center(1); 
    adjustedY = y_loc - center(2);

    % Pull distance measures:  
    D = hypot(adjustedX, adjustedY)./pix2mm; % distance to center in mm

    % add screen here for regions in the temp protocol that are to be excluded:
    tp = getTempTurnPoints(data(exp).T.TempProtocol{trial});
    all_points = [tp.DownROI,tp.UpROI,tp.HoldROI];
    loc = false(size(x_loc));
    loc(all_points,:) = true;
    if tp.holdexp
        loc = true(size(x_loc));
    end

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
    adjusted_wx = data(exp).data(trial).data.wellcenters(1,foodWell) - center(1); % adjusted well position
    adjusted_wy = data(exp).data(trial).data.wellcenters(2,foodWell) - center(2); 
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
    well_opt_x = data(exp).data(trial).data.wellcenters(1,1:4) - center(1); % adjusted well position
    well_opt_y = data(exp).data(trial).data.wellcenters(2,1:4) - center(2);
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



%% ANALYSIS: Find occupancy in the different regions 
clearvars('-except',initial_var{:})

regionList = {'fullquad','innerquad', 'quadring', 'circle10', 'circle7', 'circle5'};

% initialize empty variables
[ring.all,inner75.all] = deal([]);
region = struct;
for rr = 1:length(regionList)
    for q = 1:4 % each of the quadrants
        region(rr).(quadOrder{q}).all = [];
        region(rr).(quadOrder{q}).all_info = {'percent relative to all the flies in the whole arena'};
        if strcmp(regionList{rr},'innerquad') || strcmp(regionList{rr},'quadring')
            region(rr).(quadOrder{q}).partial = [];
            region(rr).(quadOrder{q}).partial_info = {'percent relative to only the flies in the four quadrants of this space'};
        end
    end
end

% Fill structures with all the trials processed data
for trial = 1:num.trial(exp)
    nFull = ROImask(exp).all(trial).nflies;
    nInner = ROImask(exp).inner75(trial).nflies;
    nRing = ROImask(exp).ring(trial).nflies;

    % ring & inner -- each get compared to the full arena compliment:
    n = getPercentFlies(ROImask(exp).ring(trial).m,nFull);
    ring.all = autoCat(ring.all,n,false);

    n = getPercentFlies(ROImask(exp).inner75(trial).m,nFull);
    inner75.all = autoCat(inner75.all,n,false);

    % fig = figure; hold on
    %     plot(inner75.all,'color', Color('teal'))
    %     plot(ring.all,'color', Color('gold'))
    %     plot(inner75.all + ring.all,'color', Color('white'))
    %     formatFig(fig,true);
    %     ylabel('% flies')
    %     legend({'inner', 'ring','total'},'textcolor', 'w','box', 'off');
    %     set(gca, 'xcolor', 'none')
    %     ylim([0,100])

    % for each of the quadrant-related regions:
    for rr = 1:length(regionList)
        for q = 1:4 % each of the quadrants
            temp = ROImask(exp).(regionList{rr}).(quadOrder{q})(trial).m;
            n = getPercentFlies(temp,nFull);
            region(rr).(quadOrder{q}).all = autoCat(region(rr).(quadOrder{q}).all,n,false);
            if strcmp(regionList{rr},'innerquad')
                n = getPercentFlies(temp,nInner);
                region(rr).(quadOrder{q}).partial = autoCat(region(rr).(quadOrder{q}).partial,n,false);
            elseif strcmp(regionList{rr},'quadring')
                n = getPercentFlies(temp,nRing);
                region(rr).(quadOrder{q}).partial = autoCat(region(rr).(quadOrder{q}).partial,n,false);
            end
        end
    end
end

% Pull out the averages and errors across the trials: 
ring.avg = mean(ring.all,2,'omitnan');
ring.std = std(ring.all,0,2,'omitnan');
inner75.avg = mean(inner75.all,2,'omitnan');
inner75.std = std(inner75.all,0,2,'omitnan');
for rr = 1:length(regionList)
    for q = 1:4 % each of the quadrants
        region(rr).(quadOrder{q}).avg = mean(region(rr).(quadOrder{q}).all,2,'omitnan');
        region(rr).(quadOrder{q}).std = std(region(rr).(quadOrder{q}).all,0,2,'omitnan');
        if strcmp(regionList{rr},'innerquad') | strcmp(regionList{rr},'quadring')
            region(rr).(quadOrder{q}).partial_avg = mean(region(rr).(quadOrder{q}).partial,2,'omitnan');
            region(rr).(quadOrder{q}).parital_std = std(region(rr).(quadOrder{q}).partial,0,2,'omitnan');
        end
    end
end   

% Assign the data to the grouped structure
grouped(exp).ring = ring;
grouped(exp).inner75 = inner75;
for rr = 1:length(regionList)
    grouped(exp).(regionList{rr}) = region(rr);
end


% % TODO 5.30.25: demo figure (to be moved to 4.2 something else...)
% r = 4;
% c = 1;
% sb(1).idx = 1;
% sb(2).idx = 2:4;
% lw = 2;
% kolor = {'gold', 'grey', 'white', 'grey'};
% x_lim = [0,700];
% 
% fig = getfig('',1);
% set(fig, 'windowstyle', 'docked');
% subplot(r,c,sb(1).idx);
%     x = grouped(exp).time;
%     plot(x,grouped(exp).temp,'color','w', 'linewidth',lw)
%     ylabel('(\circC)')
%     xlim(x_lim)
% subplot(r,c,sb(2).idx)
%     y_all = [];
%     hold on
%     for q = 1:4
%         y = smooth(grouped(exp).quadring.(quadOrder{q}).partial_avg,180,'moving');
%         plot(x,y,'color',Color(kolor{q}),'linewidth', lw)
%         y_all = [y_all, y];
%     end
%     plot(x, sum(y_all,2),'color', 'r')
%     xlim(x_lim)
%     ylabel('quad ring occupancy (%)')
%     xlabel('time (min)')
%     ylim([0, 100])
% formatFig(fig,true,[r,c],sb);
% subplot(r,c,sb(1).idx);
% set(gca,'xcolor', 'none');
% subplot(r,c,sb(2).idx);
% legend(quadOrder, 'textcolor', 'w', 'box', 'off');
% save_figure(fig, [saveDir 'Figures/' grouped(exp).name ' fly quadring occupancy over time'],'-png');


% Pool the data for heating and cooling together across the
% different regions for temp-tuning curve comparisons 
all_regions = [regionList, 'ring', 'inner75'];
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;

    for rr = 1:length(all_regions) % each type of region (e.g. outer ring, food circle etc)
        switch all_regions{rr}
            case {'ring','inner75'}
                nQ = 1; % quadrant number...
                nP = 1; % partial region percentages to check
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
                    baseY = grouped(exp).(all_regions{rr}).(quadOrder{q});
                else 
                    baseY = grouped(exp).(all_regions{rr});
                end
                % pull the right type of data to run
                if p==1
                    y = baseY.all;
                else
                    y = baseY.partial;
                end
                % initialize empty structures for the variables
                [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
                % [rawC_in, rawH_in, rawC_all, rawH_all] = deal(struct);
                % [rawC_in(1:4).raw, rawH_in(1:4).raw, rawC_all(1:5).raw, rawH_all(1:5).raw] = deal(nan(nTemp,num.trial(exp)));

                % Update the averages for the preset temperature bins 
                for t = 1:nTemp
                    % frame indexes for this temp bin
                    c_frames = locs(cIdx,t).frames; % cooling frames
                    h_frames = locs(hIdx,t).frames; % heating frames
                    if all(isnan(c_frames)) || all(isnan(h_frames)) % if no cooling or heating for this temp bin
                        continue
                    end
                    % pull frames associated with this temp and temp rate 
                    raw_c(t,:) = mean(y(c_frames,:),1,'omitnan'); 
                    raw_h(t,:) = mean(y(h_frames,:),1,'omitnan');
                end
                % find the avg and err and save to group structure
                h_avg = mean(raw_h, 2, 'omitnan');
                h_err = std(raw_h, 0, 2, 'omitnan');
                c_avg = mean(raw_c, 2, 'omitnan');
                c_err = std(raw_c, 0, 2, 'omitnan');


                switch all_regions{rr}
                    case {'ring','inner75'}    
                        grouped(exp).(all_regions{rr}).increasing.raw = raw_h;
                        grouped(exp).(all_regions{rr}).increasing.avg = mean(raw_h, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).decreasing.raw = raw_c;
                        grouped(exp).(all_regions{rr}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).temps = temps;
                    case {'fullquad','innerquad','quadring'}
                        if p == 1
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.raw = raw_h;
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.avg = mean(raw_h, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.raw = raw_c;
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).temps = temps;
                        else
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.p_raw = raw_h;
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.p_avg = mean(raw_h, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.p_std = std(raw_h, 0, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.p_raw = raw_c;
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.p_avg = mean(raw_c, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.p_std = std(raw_c, 0, 2, 'omitnan');
                        end
                    case {'circle10','circle7','circle5'}
                        grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.raw = raw_h;
                        grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.avg = mean(raw_h, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(quadOrder{q}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.raw = raw_c;
                        grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(quadOrder{q}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(quadOrder{q}).temps = temps;
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













