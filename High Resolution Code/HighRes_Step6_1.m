
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
    case {'courtship_F_LRR_25-17','high_res_F_LRR_25-17'}
        idx = [1 3 5 7 8]; % F LRR alignment points (beginning of each temp region)
end

MT = nan(num.trials, length(idx));
for i = 1:num.trials
    z = [fly(i).tRate(:).idx];
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
for i = 1:4
    x = i*ones(num.trials,1);
    scatter(x, y(:,i),35, fly(1).tRate(i).color,'filled', 'XJitter','density')
    plot([i-buff, i+buff], [mean(y(:,i)), mean(y(:,i))],'color', foreColor)
    xtick_lab{i} = fly(1).tRate(i).name;
end
ylabel('period duration (min)')
set(gca, 'xtick', 1:4, 'XTickLabel', xtick_lab)
title([groupName ' | n = ' num2str(num.trials)])
formatFig(fig, blkbgd);
save_figure(fig, [alignmentDir 'raw section durations'],fig_type);

% Shift trials to align to the ramp trough:
% preallocate empty space for all the variables that need to be time-shifted across trials
max_len = size(fly(1).T,1); % set up a buffered time
optA =  nan(max_len,num.trials); %single data column for each trial
optB = nan(max_len,2,num.trials); %two data columns per trial
categories = {'temperature', 'frame', 'IFD', 'cooling', 'warming', 'hold', 'wing_ext', 'wing_ext_all',...
              'court_chase', 'chase_all', 'circling_all', 'circling_1sec', 'CI'}; % from fly.T
sub_cat = {'eccentricity', 'OutterRing', 'foodQuad', 'foodcircle', 'turning'}; % from fly.data
data = struct;
for i = 1:length(categories)
    data.(categories{i}) = optA;
end
[data.dist2food, data.FlyOnFood] = deal(optB);
for i = 1:length(sub_cat)
    data.(sub_cat{i}) = optB;
end

% fill in data and shift to appropriate place...
middlePoint = median(MT(:,3)); % this is the new 'center point' in time for the aligned trials to minimize lost data across trials
middlePoint = int32(middlePoint); % make sure the number is indexable

for i = 1:num.trials
    curr_center = fly(i).tRate(3).idx(1); %current index of first temp ramp trough 
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
        [data.cooling_idx_all(:,1),~,~] = find(a==1);  % start of cooling
        [data.cooling_idx_all(:,2),~,~] = find(a==-1); % end of cooling (middlepoint)
        a = (diff(data.warming));
        [data.warming_idx_all(:,1),~,~] = find(a==1);  % start of warming
        [data.warming_idx_all(:,2),~,~] = find(a==-1); % end of warming (middlepoint)

    case 'courtship_F_LRR_25-17'
        % save new indexes for plotting regions: 
        a = (diff(data.cooling));
        [data.cooling_idx_all(:,1),~,~] = find(a==1);  % start of cooling
        [data.cooling_idx_all(:,2),~,~] = find(a==-1); % end of cooling (middlepoint)
        a = (diff(data.warming));
        [data.warming_idx_all(:,1),~,~] = find(a==1);  % start of warming
        [data.warming_idx_all(:,2),~,~] = find(a==-1); % end of warming (middlepoint)
        % universal timings: 
        data.cooling_idx = int32(median(data.cooling_idx_all));
        data.warming_idx = median(data.warming_idx_all);
        data.warming_idx(1) = data.warming_idx(1)+1; % offset the start of warming by 1 frame
        data.warming_idx = int32(data.warming_idx);
end 


% test to see how the alignment worked for the heating and cooling periods
fig = getfig('',1);
subplot(1,2,1)
plot(data.cooling)
v_line(middlePoint,'teal', '--',2)
title('cooling')
xlabel('frame')
set(gca, 'ytick', [0,1],'yticklabel', {'Off', 'On'})
ylim([-0.1,1.1])
subplot(1,2,2)
plot(data.warming)
v_line(middlePoint,'teal', '--',2)
xlabel('frame')
set(gca, 'ytick', [0,1],'yticklabel', {'Off', 'On'})
ylim([-0.1,1.1])
title('warming')
formatFig(fig, blkbgd,[1,2]);
save_figure(fig, [alignmentDir 'cooling and warming alignment'],fig_type);

[r,~,~] = find(a==1); 

% show the alignment across temperature
fig = getfig('', 1);
plot(data.temperature)
hold on
% v_line(middlePoint,'teal', '--',2)
v_line((data.cooling_idx),'dodgerblue', '--',2)
v_line((data.warming_idx),'red', '--',2)
xlabel('time')
ylabel('temperature (\circC)')
formatFig(fig,blkbgd);
title([groupName ' | n = ' num2str(num.trials)])
save_figure(fig, [alignmentDir 'final temperature alignment'],fig_type);

disp('next:')

%% TODO: create temp bin groups here and a universal temperature ramp:
clearvars('-except',initial_var{:})

% universal temperature profile: 
data.temp = smooth(mean(data.temperature,2, 'omitnan'),5*30,'moving');

temp_bins = floor(min(data.temp)):0.5:ceil(max(data.temp)); % 0.5 deg temperature bins
idx = discretize(data.temp, temp_bins);
[c_idx, h_idx] = deal(false(size(data.temp)));
c_idx(data.cooling_idx(1):data.cooling_idx(2)) = true;
h_idx(data.warming_idx(1):data.warming_idx(2)) = true;
tempbin = [];
tempbin.c_idx = c_idx;
tempbin.h_idx = h_idx;
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

%% Create new comparisons of distances for the group

% TODO (1/28) give a demo image of the ring size for each area on an image still

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
%         data.OutterRing(:,sex,i) = D>=innerR; % find the locations that are greater than the inner 75%
%         fly(i).data(sex).OutterRing = D>=innerR; % find the locations that are greater than the inner 75%
%     end
% end
% 
% disp('updated outer ring to food data structures')

% TODO: need to update the behavior matrix now with these... redo and add
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













