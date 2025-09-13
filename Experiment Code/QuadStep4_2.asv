% RUN FULL SECTION
% DataTransferGUI

% SELECT AND LOAD DATA FOR ANAYLSIS
% CAVEATS:
% ** THIS ASSUMES ALL LOADED GROUPS HAVE THE SAME WITHIN-GROUPING TEMPERATURE PROTOCOLS
% *** DOESN'T WORK FOR TEMP PROTOCOLS WITH MORE THAN 1 HEATING AND COOLING TEMP RATE

%% Select data groups to compare
% add matlabroot folder to the directory path
% addpath(genpath('C:\matlabroot'));
warning off % turn off the warnings about the oversized figures & how they might be slow to save
format shortG
clear; close all; clc
paths = getPathNames; % get the appropriate file path names
baseFolder = getDataPath(5,0,'Select where you want to find the grouped data structures');
% baseFolder = getCloudPath;
structFolder = [baseFolder paths.grouped_trials];
UpdatedFlag = false;

% Load preliminary data structure or determine which trials to load:
switch questdlg('Load existing data?','Quad Step 4 data processing','Yes','No','Cancel','Yes')
    case 'Cancel'
        return
    case 'Yes' % LOAD PRE-EXISTING DATA
            % find and select list of possible experimental groups:
            list_dirs = dir([baseFolder paths.group_comparision]);  list_dirs = {list_dirs(:).name};
            list_dirs(1:2) = [];
            [dirIdx, v] = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[500,700]);
            if v
                expGroup = list_dirs{dirIdx}; %name of experiment groups selected
            else
                disp('No group selected')
                return
            end
            % Find the group name of the existing data to load
            saveDir = [baseFolder paths.group_comparision expGroup '/'];
            filePath = [saveDir expGroup ' data.mat'];
            
            % Load exisiting data structure
            temp = load(filePath);
            % Pull  variables from structure (prevents errors in the local saving path for the current computer)
            data = temp.data;
            expNames = [];
            initial_vars = temp.initial_vars;
            num = struct; 
            num.exp = length(data);
            for i = 1:num.exp
                expNames{i} = data(i).ExpGroup;
            end
            clear temp v
    
            % Determine if any new groups need to be added/removed to structure by user
            list_dirs = dir(structFolder);  list_dirs = {list_dirs(:).name}; list_dirs(1:2) = [];
            includedIdx = nan(1,num.exp); [removeFlag, skipcheckFlag] = deal(false(1,num.exp));
            for i = 1:num.exp
                loc = find(strcmp(expNames{i},list_dirs));
                if ~isempty(loc)
                    includedIdx(i) = loc;
                else
                    disp(['Cannot find the following data structure: ' expNames{i}])
                    % Decide --  if the file doesn't exist, mark the data for deletion??
                    strg = ['WARNING: "' expNames{i} '" not found in "Data Structures." What do you want to do with this group?'];
                    switch questdlg(strg, ' ', 'Remove', 'Keep', 'Cancel','Remove')
                        case 'Remove'
                            removeFlag(i) = true;
                        case 'Keep'
                            skipcheckFlag(i) = true;
                        case 'Cancel'
                            return
                    end
                end
            end
            includedIdx(isnan(includedIdx)) = [];
            %Check that the current selection of experiments is okay by showing the included groups
            prompt_string = {'Select the groups for the structure: ', expGroup};
            expIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'multiple','ListSize',[450,600], ...
                            'InitialValue',includedIdx,'PromptString',prompt_string); clear prompt_string     

            % List of included data for comparison & updates
            dataList = struct;
            for i = 1:length(data)
                dataList(i).T = data(i).T;
                dataList(i).name = expNames{i};
                dataList(i).remove = false;
                if removeFlag(i) 
                    dataList(i).remove = true;
                else
                    dataList(i).remove = false;
                end
                 if skipcheckFlag(i) 
                    dataList(i).skipCheck = true;
                    else
                    dataList(i).skipCheck = false;
                end
            end
    
            % Determine if new data was added to the list: 
            added_data = setdiff(expIdx,includedIdx); 
            if ~isempty(added_data)
                %check that the data structures for these are all updated (from step 3.1)
                idx = size(dataList,2) + 1;
                for i = 1:length(added_data)
                    dataList(idx).name = list_dirs{added_data(i)};
                    dataList(idx).remove = false;
                    dataList(idx).skipCheck = false;
                    idx = idx +1;
                end
            end
    
            % Determine if data was removed: TODO
            remove_data = setdiff(includedIdx,expIdx);
            if ~isempty(remove_data)
                for i = 1:length(remove_data)
                     loc = strcmp(list_dirs{remove_data(i)},{dataList(:).name});
                     dataList(loc).remove = true;
                end
            end

    case 'No' % CREATE NEW DATA STRUCTURE
            UpdatedFlag = true;
            % Select processed data structures to compare:
            list_dirs = dir(structFolder);  list_dirs = {list_dirs(:).name};  list_dirs(1:2) = [];
            expIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'multiple','ListSize',[450,600]);
            expNames = list_dirs(expIdx); %name of experiment groups selected
            num.exp = length(expIdx);  %number of groups selected
           
            % Load selected experiment data groups
            for i = 1:num.exp

                % get field list for loading data:
                dummy = load([structFolder expNames{i} '/' expNames{i} ' post 3.2 data.mat']);
                if isfield(dummy,'autoRun')
                    dummy = rmfield(dummy,'autoRun');
                end
                %field not present in dummy, so add blank
                if i > 1
                    unmatched_in_dummy = setdiff(fields(data),fields(dummy));
                else
                    unmatched_in_dummy = [];
                end
                if ~isempty(unmatched_in_dummy)
                    for ii = 1:length(unmatched_in_dummy)
                        dummy.(unmatched_in_dummy{ii}) = [];
                    end
                end
                if ~isfield(dummy,'hold_exp') % 1/24/24 updates for hold temperature experiments
                      dummy.hold_exp = false; % account for new data structures
                      dummy.temp_protocol = dummy.T.TempProtocol{1};
                end
                if i == 1
                    data = dummy;
                else
                    data(i) = dummy;
                end
            end
    
            clear list_dirs expIdx dirIdx dummy
            % Set up base variables
            initial_vars = who;
            initial_vars = [initial_vars(:); 'initial_vars'; 'grouped'; 'expGroup'; 'saveDir'; 'mat';'expOrder'];
            initial_vars = unique(initial_vars);
    
            % List of included data for comparison & updates
            dataList = struct;
            for i = 1:length(data)
                dataList(i).T = data(i).T;
                dataList(i).name = expNames{i};
                dataList(i).remove = false;
                dataList(i).skipCheck = false;
            end
end

% CHECK FOR DATA UPDATES OF INCLUDED STRUCTURES        
% Find the files within each data set and compare to the existing base 
% data structure for missing data & then updata any missing data if possible
for i = 1:size(dataList,2)
    if dataList(i).remove || dataList(i).skipCheck % dont update data that will be deleted
        dataList(i).rebuild = false;
        dataList(i).add = false;
        dataList(i).extradata = false;      
        continue
    end
    % check the included data structure with the data list in the folder 
    [a,b,c] = matchDataStructure(dataList(i).name, dataList(i).T, structFolder);
    dataList(i).rebuild = a; % data in the struct folder doesn't match the data list, must rebuild from 3.1
    dataList(i).add = b;    % add/update data to the current structure
    dataList(i).extradata = c; % more data in the struct than in the data list           
end; clear a b c

% temp = data;

% Update the grouped structure according to the dataList
for i = 1:size(dataList, 2)
    exp_name = dataList(i).name;
    % remove unwanted data sets
        if dataList(i).remove
            loc = find(strcmp(exp_name,{data(:).ExpGroup}));
            data(loc) = [];
            UpdatedFlag = true;
            continue
        end
    % Load or reload existing datasets
        if any([dataList(i).rebuild, dataList(i).extradata])
            disp([exp_name ' needs to be rebuilt from Step 3.2'])
            switch questdlg(['"' exp_name ''' is not up-to-date in this structure. Continue anyway?'])
                case 'No'
                        disp('Check for other rebuilds in the structure')
                        return
                case 'Cancel'
                    return
            end
        end
        % ADD or UPDATE data structures
        if dataList(i).add
            disp(['Updating data for ' exp_name])
            UpdatedFlag = true;
            loc = find(strcmp(exp_name,{data(:).ExpGroup})); 
            dummy = load([structFolder exp_name '/' exp_name ' post 3.2 data.mat']);
             
            % check for structure alignment and potential mis-matches:
             if ~isfield(dummy,'hold_exp') % 1/24/24 compatability for hold temperature experiments
                   dummy.hold_exp = false; % account for new data structures
                   dummy.temp_protocol = dummy.T.TempProtocol{1};
             end

             struct_fields = fields(data);
             new_fields = fields(dummy);
             missing_in_dummy = setdiff(struct_fields, new_fields);
             if ~isempty(missing_in_dummy)
                 for xx = 1:length(missing_in_dummy)
                     dummy.(missing_in_dummy{xx}) = [];
                     disp(['Adding "' missing_in_dummy{xx} '" field to ' dummy.ExpGroup])
                 end
             end
             missing_in_data = setdiff(new_fields, struct_fields);
             if ~isempty(missing_in_data)
                 for xx = 1:length(missing_in_data)
                     data(1).(missing_in_data{xx}) = [];
                     disp(['Adding "' missing_in_data{xx} '" field to ' dummy.ExpGroup])
                 end
             end

             % determine location of data in structure or add to structure if not there: 
             if isempty(loc)
                 loc = length(data)+1;
             end
             data(loc) = dummy;
        end

        % Check and add field to preexisting data that needs this:
        loc = find(strcmp(exp_name,{data(:).ExpGroup}));
        a = false;
        try a = isstring(data(loc).temp_protocol);
        catch            
        end
        if ~a
            data(loc).temp_protocol = data(loc).T.TempProtocol{1};
        end
        a = false;
        try a = islogical(data(loc).hold_exp);
        catch            
        end
        if ~a
            data(loc).hold_exp = false; % account for new data structures
        end
end 
clear dummy missing_in_data missing_in_dummy struct_fields new_fields

% clear blank space **this might differ for a trial that ends with a removal???
i = 1; idx = 0; maxIdx = size(data,2);
while i <= size(data,2) 
    idx = idx + 1;
    % Check if the current element should be removed
    if isempty(data(i).ExpGroup)
        % Remove the blank slots
        data(i) = [];
    else
        % Keep group and move to next
        i = i + 1;
    end
    % Break the loop if idx exceeds the original max index
    if idx > maxIdx
        break;
    end
end

num.exp = size(data,2);

%check for temperature protocol assignments (make back-compatible)
for i = 1:num.exp
    if ~isfield(data(i),'temp_protocol')
        data(i).temp_protocol = data(i).T.TempProtocol{1};
    end
    if isempty(data(i).temp_protocol)
        data(i).temp_protocol = data(i).T.TempProtocol{1};
    end
end

%check for temperature hold notation (make back-compatible)
for i = 1:num.exp
    if ~isfield(data(i),'hold_exp')
        data(i).hold_exp = false;
    end
    if isempty(data(i).hold_exp)
        data(i).hold_exp = false;
    end
end

% Save data / make new grouped data folder
result = cellfun(@isnumeric, initial_vars);
initial_vars(result) = [];
initial_vars{end+1} = 'UpdatedFlag';
clearvars('-except',initial_vars{:})
result = cellfun(@(x) strcmpi(x, 'UpdatedFlag'), initial_vars);
if result(end)==1
    initial_vars  = initial_vars(1:end-1);
end

% List of included data for comparison & updates
expNames = [];
dataList = struct;
for i = 1:length(data)
    dataList(i).T = data(i).T;
    dataList(i).name = data(i).ExpGroup;
    expNames{i} = data(i).ExpGroup;
end

paths = getPathNames;

if UpdatedFlag
    switch questdlg('Select data saving format:','','new structure','existing structure', 'cancel','existing structure')
    case 'new structure'
        expGroup = char(inputdlg('Structure name:'));
        saveDir = [baseFolder paths.group_comparision expGroup '/'];
        if ~exist(saveDir,'dir')
            mkdir(saveDir);
        end 
        save([saveDir expGroup ' data.mat'],'-v7.3');
        disp([expGroup ' saved'])
    case 'existing structure'
        list_dirs = dir([baseFolder paths.group_comparision]);
        list_dirs = {list_dirs(:).name};
        list_dirs(1:2) = [];
        if exist('expGroup','var')
            guessLoc = find(strcmp(expGroup,list_dirs)); % TODO check if there is an existing folder and then offer that as the base selection
            dirIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[350,500], 'InitialValue',guessLoc);
        else
            dirIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[350,500]);
        end
        expGroup = list_dirs{dirIdx}; %name of experiment groups selected
        saveDir = [baseFolder paths.group_comparision expGroup '/'];
        save([saveDir expGroup ' data.mat'],'-v7.3');
        disp([expGroup ' saved'])
    case 'cancel'
        return
    end
else % append the dataList structure to the existing file if nothing else changed
    saveDir = [baseFolder paths.group_comparision expGroup '/'];
    % save([saveDir expGroup ' data.mat'], 'dataList','-append');
end

% Find plate identity for all the trials in the experiment
% pull excel information 
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;

% extract the arena plate information for each trial: 
for i = 1:num.exp
    num.trial(i) = data(i).ntrials;
    data(i).plate = nan(num.trial(i),1);
    % find each trial in the excel file and extract the plate ID
    T = data(i).T;
    for trial = 1:num.trial(i)
        date_str = T.Date{trial};
        arena_str = T.Arena{trial};
        exp_str = T.ExperimentID{trial};
        loc1 = (strcmp(excelfile(:,Excel.date),date_str));
        loc2 = (strcmp(excelfile(:,Excel.arena),arena_str));
        loc3 = (strcmp(excelfile(:,Excel.expID),exp_str));
        loc = (loc1 & loc2 & loc3);
        data(i).plate(trial) = excelfile{loc, Excel.plate};
        [~, con_type] = getConversion(date_str, data(i).plate(trial),1);
        data(i).con_type(trial) = con_type;
    end
end

disp(expGroup)
disp(expNames')

figDir = createFolder([saveDir, 'Figures/']);
initial_vars{end+1} = 'figDir';

%% ANALYSIS: organize data for each group
clearvars('-except',initial_vars{:})
fig_type = '-png'; 
blkbgd = true;
initial_vars = [initial_vars(:); 'initial_vars'; 'grouped'; 'expGroup'; 'saveDir'; 'mat';'expOrder'; 'fig_type';'f2m';'conversion';'blkbgd']; % changed pix2mm to conversion
initial_vars = unique(initial_vars);
% f2m = 3*60; %fps*min = number of frames in a minute -- this is now held
% in the getTempRate function
grouped = struct;

% Color selections
switch expGroup
   case 'Berlin LTS 15-35 caviar vs empty'
        expOrder = 1:2;
        colors = {'DodgerBlue', 'Gray'};
     case 'Berlin LTS 15-35 plate 1 vs plate 2'
        expOrder = 1:4;
        colors = {'Tomato', 'Dodgerblue', 'Peachpuff', 'Powderblue'};  
    case 'Berlin LTS 15-35 intact vs no antenna no food'
        expOrder = 1:2;
        colors = {'dodgerblue', 'peachpuff'};
    case 'Berlin F LRR 25-17 caviar intact vs wax vs hold'
        expOrder = 1:4;
        colors = {'Dodgerblue', 'Tomato', 'Grey', 'Grey'};
    case 'Berlin temp rate caviar'
        expOrder = [5, 3, 2, 1, 4]; % slow to fast
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};
    case 'Berlin LTS 15-35 no food mechanical removal comparisons'
        expOrder = [1,2,3]; % berlin, no antenna, no arista? 
        colors = {'Grey', 'MediumSlateBlue', 'Gold'};
    case 'Berlin temperature holds' % includes food and empty trials
        expOrder = 1:num.exp;
        kolor1 = Color('Blue', 'LightGrey', 5); % ultimately 5 (with all temp trials added)
        kolor2 = Color('LightGrey', 'Red', 5);  % ultimately 5 (with all temp trials added)
        colors = nan([18, 3]); % ultimately 18x3 (with all temp trials added)
        colors(1:2:end,:) = [kolor1(1:end-1,:); kolor2];
        colors(2:2:end,:) = [kolor1(1:end-1,:); kolor2];
        % figure; scatter(1:num.exp, 1:num.exp, 50, colors(1:num.exp,:), 'filled'); % color test
    case {'TrpA1-Gal4 x UAS-Kir2.1 LTS 15-35 no food','TrpA1-gal4 LTS 15-35 no food'}
        expOrder = 1:num.exp; % UAS control, GAL4 control, GAL4>UAS
        colors = {'Peachpuff', 'Powderblue','Magenta'};
    case {'TrpA1-Gal4 x UAS-Kir2.1_A1 LTS 15-35 caviar','TrpA1-Gal4 x UAS-Kir2.1_A1 LTS 15-35 caviar plate 1'}
        expOrder = 1:num.exp; % UAS control, GAL4 control, GAL4>UAS
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
    case 'TrpA1-gal4 x Kir2.1 no antenna LTS 15-35 no food'
        expOrder = 1:2; % intact vs antenna-less
        colors = {'Grey', 'DodgerBlue'};
    case 'Berlin vs UAS-Kir2.1 caviar background comparison'
        expOrder = 1:num.exp;
        colors = {'turquoise', 'DarkOrchid'};
    case 'IR25a-gal4 x Kir2.1 and controls LTS 15-25 caviar'
        expOrder = 1:3;
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
    case 'R77C10-gal4 x Kir2.1 and controls F LRR 17-25 caviar'
        expOrder = 1:3; % gal4 control, UAS control, 
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
    case 'IR40a_TM2-gal4 x Kir2.1 and controls F LRR 17-25 caviar'
        expOrder = 1:3; % gal4 control, UAS control, 
        colors = {'LightPink', 'HotPink', 'DodgerBlue'};
    case 'IR40a_TM2-gal4 x Kir2.1 no arista LTS 15-35 empty'
        expOrder = 1:3;
        colors = {'magenta', 'grey', 'dodgerblue'};

end

% % Reset colors if needed
% for i = 1:num.exp
%     exp = expOrder(i);
%     grouped(exp).color = Color(colors{i});
% end

if ~exist('colors','var')
    expOrder = 1:num.exp;
    if num.exp<=7
        colors = {'DarkOrchid','Gold','dodgerblue','turquoise','lime','red','Orange'};
    else
        cMap = colormap(turbo(num.exp+2));
        cMap(1,:) = []; cMap(num.exp+1,:) = []; %remove ends since they are too dark for black background
    end
end
% Display the experiment 'order'
disp('Experiment order:')
for exp = 1:num.exp
    i = expOrder(exp);
    disp(expNames{i})
end

% % update color map without running the rest of the data: 
% for i = 1:num.exp
%     grouped(expOrder(i)).color = Color(colors{(i)});
% end

for i = 1:num.exp % FOR EACH DATA GROUP
    % GENERAL
    if ~isfield(data,'fps')
        data(i).fps = 3;
    end
    grouped(i).name = data(i).ExpGroup;
    if exist('colors','var')
        if iscell(colors)
            grouped(expOrder(i)).color = Color(colors{(i)});
        else
            grouped(expOrder(i)).color = colors(i,:);
        end
    else 
        grouped(expOrder(i)).color = cMap(i,:);
    end
    % TIME COURSE DATA
    [time,temp,speed,distance] = deal([]);
    for trial = 1:num.trial(i)
        time = autoCat(time,data(i).data(trial).occupancy.time,false,true);
        temp = autoCat(temp,data(i).data(trial).occupancy.temp,false,true);
        speed = autoCat(speed,data(i).data(trial).speed.avg,false,true);
    % movement = autoCat(speed,data(i).data(trial).speed.avg,false,true);
        distance = autoCat(distance,data(i).data(trial).dist2wells(:,data(i).T.foodLoc(trial)),false,true);
    end
    grouped(i).time = mean(time,2,'omitnan');
    grouped(i).temp = mean(temp,2,'omitnan');
    grouped(i).speed.all = speed;
    grouped(i).speed.avg = mean(speed,2,'omitnan');
    grouped(i).speed.err = std(speed,0,2,'omitnan')/sqrt(num.trial(i));
    grouped(i).speed.zscore.all = zscore(grouped(i).speed.all);
    grouped(i).speed.zscore.avg = mean(grouped(i).speed.zscore.all,2,'omitnan');
    grouped(i).dist.all = distance;
    grouped(i).dist.avg = mean(distance,2,'omitnan');
    grouped(i).dist.err = std(distance,0,2,'omitnan')/sqrt(num.trial(i));
    grouped(i).dist.zscore.all = zscore(grouped(i).dist.all);
    grouped(i).dist.zscore.avg = mean(grouped(i).dist.zscore.all,2,'omitnan');

    % AVG POSITION BINNED BY TEMP
    for trial = 1:num.trial(i)
        tempList = data(i).G(trial).TR.temps;
        y = data(i).G(trial).TR.data(:,3);
        for tt = 1:length(tempList)
            grouped(i).dist.tempBinned(trial,tt) = mean(y(data(i).G(trial).TR.tempIdx==tt),'omitnan');
        end
        grouped(i).dist.tempList(trial,:) = tempList;
    end
    grouped(i).dist.distavgbytemp = [grouped(i).dist.tempList(1,:)',...
                                     mean(grouped(i).dist.tempBinned,1,'omitnan')'];
    grouped(i).dist.distavgbytemp_err = [grouped(i).dist.tempList(1,:)',...
                                         (std(grouped(i).dist.tempBinned,0,'omitnan')/sqrt(num.trial(i)))'];

    % BINNED
    [tempRates,decreasing,increasing,temperatures] = deal([]);
    for trial = 1:num.trial(i)
        % Account for multiple numbers of rate trials
        rates = data(i).G(1).TR.rates;
        nRates(i) = size(rates,2);
        if nRates(i)==2
            blankdata = data(i).G(trial).TR.dist_mat.avg;
        elseif nRates == 3
            blankdata = data(i).G(trial).TR.dist_mat.avg;
            blankdata(rates==0,:) = [];
            rates(rates==0) = [];
        else
            warndlg('Temp protocol has too many rates for this analysis')
            return
        end
        downIdx = find(rates<0);
        upIdx = find(rates>0);
        decreasing(:,trial) = blankdata(downIdx,:);
        increasing(:,trial) = blankdata(upIdx,:);
        tempRates = autoCat(tempRates,rates,true,true);
        temperatures(:,trial) = data(i).G(trial).TR.temps;
    end
    grouped(i).increasing.temps = median(temperatures,2);
    grouped(i).increasing.all = increasing;
    grouped(i).increasing.avg = mean(increasing,2,'omitnan');
    grouped(i).increasing.zscore.all = zscore(increasing,0,'omitnan');
    grouped(i).increasing.zscore.avg = mean(grouped(i).increasing.zscore.all,2,'omitnan');
    grouped(i).increasing.err = std(increasing,0,2,'omitnan')/sqrt(num.trial(i));
    grouped(i).increasing.rate = median(tempRates(:,upIdx));
    grouped(i).decreasing.temps = median(temperatures,2);
    grouped(i).decreasing.all = decreasing;
    grouped(i).decreasing.avg = mean(decreasing,2,'omitnan');
    grouped(i).decreasing.zscore.all = zscore(decreasing,0,'omitnan');
    grouped(i).decreasing.zscore.avg = mean(grouped(i).decreasing.zscore.all,2,'omitnan');
    grouped(i).decreasing.err = std(decreasing,0,2,'omitnan')/sqrt(num.trial(i));
    grouped(i).decreasing.rate = median(tempRates(:,downIdx));

end

% REAL TEMP TIME COURSE FOR EMPTY TRIALS
for exp = 1:num.exp
    if data(exp).hold_exp
        str_fill = 'real_temp';
    else 
        str_fill = 'temp';
    end
    temp = [];
    for trial = 1:num.trial(exp)
        temp = autoCat(temp,data(exp).data(trial).occupancy.(str_fill),false,true);
    end
    grouped(exp).real_temp = mean(temp,2,'omitnan');
end


% RAMP-TO-RAMP COMPARISONS
binWidth = 0.5; %temp bin increment
mat = struct;
% Extract ramp by ramp information from each trial
for i = 1:num.exp
    for trial = 1:num.trial(i)
        tempPoints = getTempTurnPoints(data(i).T.TempProtocol{trial});
        foodLoc = data(i).T.foodLoc(trial);
        threshLow = tempPoints.threshLow;
        threshHigh = tempPoints.threshHigh;
        dist = data(i).data(trial).occupancy.dist2wells(:,foodLoc);
        temp = data(i).data(trial).occupancy.temp;
        tempBins = floor(threshLow):binWidth:ceil(threshHigh);
        if tempBins(end)<ceil(threshHigh)
            tempBins(end+1) = ceil(threshHigh) + binSpace;
        end

        for idx = 1:tempPoints.nDown
            % COOLING
            downROI = tempPoints.down(idx,1):tempPoints.down(idx,2);
            mat(i).cooling(idx).dist(:,trial) = dist(downROI);
            mat(i).cooling(idx).temp(:,trial) = temp(downROI);
            mat(i).cooling(idx).tempIdx(:,trial) = discretize(temp(downROI),tempBins);

            % HEATING
            upROI = tempPoints.up(idx,1):tempPoints.up(idx,2);
            mat(i).heating(idx).dist(:,trial) = dist(upROI);
            mat(i).heating(idx).temp(:,trial) = temp(upROI);
            mat(i).heating(idx).tempIdx(:,trial) = discretize(temp(upROI),tempBins);

            % HYSTERESIS MEASURE
            for n = 1:length(tempBins)-1 %pull avg temp for this particular cooling ramp
                %cooling
                C_loc = mat(i).cooling(idx).tempIdx(:,trial)==n;
                mat(i).cooling(idx).bintemp(n,trial) = mean(mat(i).cooling(idx).dist(C_loc,trial),'omitnan');
                %heating
                H_loc = mat(i).heating(idx).tempIdx(:,trial)==n;
                mat(i).heating(idx).bintemp(n,trial) = mean(mat(i).heating(idx).dist(H_loc,trial),'omitnan');
            end
            mat(i).hysteresis(idx).all(:,trial) = mat(i).cooling(idx).bintemp(:,trial)-mat(i).heating(idx).bintemp(:,trial);
            mat(i).hysteresis(idx).sum(trial) = sum(mat(i).hysteresis(idx).all(:,trial),'omitnan');
            mat(i).hysteresis(idx).tempbins = tempBins;
            mat(i).distHist(:,trial,idx) = mat(i).hysteresis(idx).all(:,trial); %hyst by temp bin
            mat(i).cumHist(trial,idx) = mat(i).hysteresis(idx).sum(trial); %cumulative hysteresis
        end
    end
end

conversion = getConversion;

%% ANALYSIS: Find the quadrants with the highest and lowest occupancy over the exp for null comparison
% as defined by the 7% space around the well (aka using the 'occupancy' measure
clearvars('-except',initial_vars{:})

for exp = 1:num.exp
    grouped(exp).occ_idx = (nan([num.trial(exp),2])); % initialize an empty occupancy index
    for trial = 1:num.trial(exp) 
        quad_avg = mean(data(exp).data(trial).data.occ_P,1,'omitnan'); % pull the percent occupancy for each food well
        [~, low_idx] = min(quad_avg); % which quad was lowest occ on avg
        [~, high_idx] = max(quad_avg); % which quad was highest occ on avg
        grouped(exp).occ_idx(trial,1) = low_idx; % register into the index the quad identity for lowest occupancy 
        grouped(exp).occ_idx(trial,2) = high_idx; % register into the index the quad identity for highest occupancy 
    end
    % add a logical for if this is an empty trial (make it easier for switching code sections later
    if strcmp(data(exp).foodNames, 'Movement')
        data(exp).emptytrial = true;
    else
        data(exp).emptytrial = false;
    end
end

%% ANALYSIS: normalize fly position within arena to common orientation (food location) 
clearvars('-except',initial_vars{:})
OG_Orientation = datetime('10.20.2023','InputFormat','MM.dd.yyyy'); % camera & lens change accounting

for type = 1:3 % assigned food well, lowest occupied, highest occupied
    switch type
        case 1 % food assigned position 
            field_label = 'position';
        case 2 % lowest occupancy position aligned
            field_label = 'position_low';
        case 3 % highest occupancy position aligned
            field_label = 'position_high';
    end

    for i = 1:num.exp
      % BINNED
      [tempRates,decreasing,increasing,temperatures] = deal([]);
      % Pull temperature and rate information for the temp protocol
      temp_rates = data(i).G(1).TR.rates;
      temp_list = data(i).G(1).TR.temps;
      nrates = data(i).G(1).TR.nRates;
      ntemps = data(i).G(1).TR.nTemps;
      rate_idx = data(i).G(1).TR.rateIdx;
      temp_idx = data(i).G(1).TR.tempIdx;
      loc_mat = struct;
    
      % find frame index for each temp bin & rate of change
      for tt = 1:ntemps
        for rr = 1:nrates
            rateAligned = rate_idx==rr;
            tempAligned = temp_idx==tt;
            loc = find(rateAligned & tempAligned);
            if isempty(loc)
                loc = nan;
            end
            loc_mat(rr,tt).frames = loc;
            loc_mat(rr,tt).x = []; % set empty space for appending
            loc_mat(rr,tt).y = [];
        end
      end
      wellXY = [];
      for trial = 1:num.trial(i)
        trial_date = datetime(data(i).T.Date{trial},'InputFormat','MM.dd.yyyy');
        
        % get arena information
         switch type
             case 1 % food assigned position 
                well_loc = data(i).T.foodLoc(trial);
             case 2 % lowest occupancy position aligned
                well_loc = grouped(i).occ_idx(trial,1);
             case 3 % highest occupancy position aligned
                well_loc = grouped(i).occ_idx(trial,2);
        end
        wells = data(i).data(trial).data.wellcenters;
    
        % % find offset to make the food well the origin
        % x_offset = wells(1,well_loc);
        % y_offset = wells(2,well_loc);
        % wells_x = wells(1,:)-x_offset;
        % wells_y = wells(2,:)-y_offset;
        % 
        % X = data(i).data(trial).data.x_loc;
        % Y = data(i).data(trial).data.y_loc;
        % X = X-x_offset;
        % Y = Y-y_offset;
        % % save the position normalized data into the grouped structure or something
    
        % use the arena center as rotation center rather than the wells
        arena_center_x = mean(wells(1,:));
        arena_center_y = mean(wells(2,:));
        
        % offset all data by arena center
        X = data(i).data(trial).data.x_loc - arena_center_x;
        Y = data(i).data(trial).data.y_loc - arena_center_y;
        
        wells_x = wells(1,:) - arena_center_x;
        wells_y = wells(2,:) - arena_center_y;
    
        [WELLS,x_data,y_data] = deal([]);
        if trial_date > OG_Orientation %new camera orientation
            % Rotate to correct orientation
            switch well_loc
                case 1
                    x_data = X;
                    y_data = Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
                case 2
                    x_data = Y;
                    y_data = -X;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 3
                    x_data = X;
                    y_data = -Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 4
                    x_data = -Y;
                    y_data = X;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
            end
        else % Rotate to correct orientation with older camera arrangement            
            switch well_loc
                case 1
                    x_data = Y;
                    y_data = -X;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 2
                    x_data = X;
                    y_data = -Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 3
                    x_data = -Y;
                    y_data = X;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
                case 4
                    x_data = X;
                    y_data = Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
            end
    
        end
        
        wellXY.x(:,trial) = WELLS(:,1);
        wellXY.y(:,trial) = WELLS(:,2);
    
        % For each temp and rate, pool all (now normalized) fly positions
        for tt = 1:ntemps
          for rr = 1:nrates
            frame_idx = loc_mat(rr,tt).frames;
            % shift positions for food well as origin
            if ~isnan(frame_idx)
                x = x_data(frame_idx,:);
                y = y_data(frame_idx,:);
                
                % save data into structure:
                loc_mat(rr,tt).data(trial).pos = [x(:),y(:)];
                loc_mat(rr,tt).x = [loc_mat(rr,tt).x; x(:)];
                loc_mat(rr,tt).y = [loc_mat(rr,tt).y; y(:)];
            else
                loc_mat(rr,tt).data(trial).pos = [nan nan];
                loc_mat(rr,tt).x = nan;
                loc_mat(rr,tt).y = nan;
            end
          end
        end
        grouped(i).(field_label).trial(trial).x = x_data; % this is the full x data for each trial reoriented! 
        grouped(i).(field_label).trial(trial).y = y_data; % this is the full y data for each trial reoriented! 
      end
      
      % save trial data into broad structure
      grouped(i).(field_label).loc = loc_mat;
      grouped(i).(field_label).temp_rates = temp_rates;
      grouped(i).(field_label).temp_list = temp_list; %
      grouped(i).(field_label).well_pos = wellXY; %position of the wells in the newly oriented arena
    end

end


%% Visual test of arena alignment and overlay ... make sure the alignment actually worked! 

field_label = 'position_high';

[c, r] = subplot_numbers(num.exp);
fig = getfig(field_label,1);
for exp = 1:num.exp
    subplot(r,c,exp)
    % plot the well positions as a '+'
    barX1 = grouped(exp).(field_label).well_pos.x([1,3],:);
    barX2 = grouped(exp).(field_label).well_pos.x([2,4],:);
    barY1 = grouped(exp).(field_label).well_pos.y([1,3],:);
    barY2 = grouped(exp).(field_label).well_pos.y([2,4],:);
    centreX = grouped(exp).(field_label).well_pos.x(5,:);
    centreY = grouped(exp).(field_label).well_pos.y(5,:);
    hold on
        plot(barX1, barY1,'Marker','.','MarkerSize',15)
        plot(barX2, barY2,'Marker','.','MarkerSize',15)
        scatter(centreX, centreY,35,'filled')
        title(grouped(exp).name,'color', 'w')
end
formatFig(fig, true,[r,c]);

%% ANALYSIS: calculate group occupancy for food circle
% this is fairly legacy at this point (8/4/25) and the new quad occupancy
% should be used over this metric moving forward. This uses the 7% circle
% data
clearvars('-except',initial_vars{:})

% Pull the data together: 
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
    all_occ = [];
    % combine all the occupancy for a given temp bin across trials:
    idx = 0;
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        else idx = idx + 1;
        end
        for trial = 1:num.trial(exp)
            % pull locations for temp and speed bins
            well_loc = data(exp).T.foodLoc(trial);
            y = data(exp).data(trial).occupancy.occ(:,well_loc);
            raw_c(t,trial) = mean(y(c_frames),'omitnan');
            raw_h(t,trial) = mean(y(h_frames),'omitnan');
            if idx == 1
                all_occ = autoCat(all_occ, y,false);
            end
        end
    end
    
    % find the avg and err and save to group structure
    grouped(exp).occ.increasing.raw = raw_h;
    grouped(exp).occ.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).occ.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).occ.decreasing.raw = raw_c;
    grouped(exp).occ.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).occ.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).occ.temps = temps;
    grouped(exp).occ.all = all_occ;
    grouped(exp).occ.avg = mean(all_occ,2,'omitnan');
    grouped(exp).occ.std = std(all_occ,0,2,'omitnan');
    grouped(exp).occ.warning = '10% ring : this has been adjusted for each plate size 5.25.25';
end

%% ANALYSIS: pull binned speed
clearvars('-except',initial_vars{:})

% Pull the data together:
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp

    % combine all the occupancy for a given temp bin across trials:
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        end
        % pull locations for temp and speed bins
        raw_c(t,:) = mean(grouped(exp).speed.all(c_frames,:),1,'omitnan');
        raw_h(t,:) = mean(grouped(exp).speed.all(h_frames,:),1,'omitnan');
    end

    % find the avg and err and save to group structure
    grouped(exp).speed.increasing.raw = raw_h;
    grouped(exp).speed.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).speed.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).speed.decreasing.raw = raw_c;
    grouped(exp).speed.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).speed.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).speed.temps = temps;
end

% Add eccentricity data
for i = 1:num.exp
    ecent = [];
    for trial = 1:num.trial(i)
        y = data(i).data(trial).data.occupancy.eccentricity(:,1);
        ecent = autoCat(ecent,y,false);
    end
    grouped(i).ecent.all = ecent;
end

%% ANALYSIS: Calculate flies within the outer 25% of the arena ring of the region [w/ fig of spatial distributions]
% clearvars('-except',initial_vars{:})
% conversion = getConversion;  
% 
% % % Max Distance from Center of Arena 
% % % PLATE 1: 
% % R1 = 30; %mm for plate 1
% % num.R1 = R1;
% % innerR1 = R1*sqrt(3/4); % radius of the inner 25% occupancy space R*sqrt(1/2)
% % dist_from_edge1 = (R1 - innerR1);
% % maxR1 = R1*sqrt(0.1); % radius of a circle occupying 10% of the arena
% % % PLATE 2: 
% % R2 = 25.6;%mm for plate 2
% % num.R2 = R2;
% % innerR2 = R2*sqrt(3/4); % radius of the inner 25% occupancy space R*sqrt(1/2)
% % dist_from_edge2 = (R2 - innerR2); % distance acceptable for start of outer 25%
% % maxR2 = R2*sqrt(0.1); % radius of a circle occupying 10% of the arena
% 
% % Find the percent of the flies that are in the outer ring
% for exp = 1:num.exp
%     counts = []; ring_per = [];
%     for trial = 1:num.trial(exp)
%         conType = data(exp).con_type(trial); % flag for the specif plate and experiment type
%         R = conversion(conType).R; %outer reaches of the arena 
%         innerR = conversion(conType).circle75; % inner limits of a 25% outer ring line
% 
%         x = data(exp).data(trial).data.x_loc; % x locations for the entire experiment
%         y = data(exp).data(trial).data.y_loc; % x locations for the entire experiment
%         centre = data(exp).data(trial).data.centre; %distance to center of arena 
%         dX = (x-centre(1)); % x-difference from center
%         dY = (y-centre(2)); % y-difference from center
%         D = hypot(dX,dY);  % distance from center of arena
%         D = D./conversion(conType).pix2mm; % in mm distance to center of the arena 
%         loc = D<=R & D>=innerR; % find the locations that are between edge (R) and inner R
%         ringCount = sum(loc,2);
%         counts = autoCat(counts, ringCount, false); %count #flies in the outer ring
%         ring_per = autoCat(ring_per,(ringCount./data(exp).T.NumFlies(trial)).*100,false); % convert to percent & combine
%     end
%     % pool the data
%     grouped(exp).ring.counts = counts;
%     grouped(exp).ring.all = ring_per;
%     grouped(exp).ring.percent = ring_per; 
%     grouped(exp).ring.avg = mean(ring_per,2,'omitnan');
% end
% 
% % Pull the data together: 
% for exp = 1:num.exp
%     temps = grouped(exp).position.temp_list; % pre-binned temperatures
%     nTemp = length(temps);
%     rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
%     cIdx = find(rates<0); %cooling index
%     hIdx = find(rates>0); %heating index
%     locs = grouped(exp).position.loc;
%     [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
%     all_ring = [];
% 
%     % Update the averages for the classic temperature bins 
%     for t = 1:nTemp
%         % cooling frames for this temp
%         c_frames = locs(cIdx,t).frames;
%         h_frames = locs(hIdx,t).frames;
%         if all(isnan(c_frames)) || all(isnan(h_frames))
%             continue
%         end
%         raw_c(t,:) = mean(grouped(exp).ring.percent(c_frames,:),1,'omitnan');
%         raw_h(t,:) = mean(grouped(exp).ring.percent(h_frames,:),1,'omitnan');
%     end
% 
%     % find the avg and err and save to group structure
%     grouped(exp).ring.increasing.raw = raw_h;
%     grouped(exp).ring.increasing.avg = mean(raw_h, 2, 'omitnan');
%     grouped(exp).ring.increasing.err = std(raw_h, 0, 2, 'omitnan');
%     grouped(exp).ring.decreasing.raw = raw_c;
%     grouped(exp).ring.decreasing.avg = mean(raw_c, 2, 'omitnan');
%     grouped(exp).ring.decreasing.err = std(raw_c, 0, 2, 'omitnan');
%     grouped(exp).ring.temps = temps;
% end
% 
% % Visual comparision of the fly positions within the arena for each trial
% skipstep = 100; lw = 0.5;
% for exp = 1:num.exp
%     r = floor(sqrt(num.trial(exp)));
%     c = ceil(num.trial(exp)/r);
%     fig = getfig(grouped(exp).name,1); hold on
%     for trial = num.trial(exp):-1:1 %1:num.trial(exp)
%         subplot(r,c,trial); hold on
%         x = grouped(exp).position.trial(trial).x;
%         y = grouped(exp).position.trial(trial).y;
%         X = x(1:skipstep:end,:); X = X(:);
%         Y = y(1:skipstep:end,:); Y = Y(:);
%         conType = data(exp).con_type(trial);
%         switch conType
%             case {1,2} % old plate
%                 kolor = Color('yellow');
%                 sz = 1;
%             case 3 % new plate
%                 kolor = Color('magenta');
%                 sz = 2;
%         end
%         R = conversion(conType).R*conversion(conType).pix2mm; % outer limits of the arena radius for comparison
%         innerR = conversion(conType).circle75*conversion(conType).pix2mm;
%         scatter(X,Y,sz,kolor,'filled'); % draw 'all' the fly positions (some skipped for size constraints)
%         scatter(0, 0, 20,Color('gold')) % draw food location: (this is centered and aligned data to the food well at (0,0)       
%         cX = grouped(exp).position.well_pos.x(5,trial);
%         cY = grouped(exp).position.well_pos.y(5,trial);
%         axis square equal
%         viscircles([cX, cY],R,'Color','white','LineWidth',lw); % outer expanse of arena accessiblity 
%         viscircles([cX, cY],innerR,'Color','white','LineWidth',lw); % outer expanse of arena accessiblity 
%         formatFig(fig,blkbgd, [r,c]);
%         for i = 1:num.trial(exp)
%             subplot(r,c,i)
%             set(gca, 'xcolor', 'none', 'ycolor', 'none')
%         end
%     end
%     figDir = createFolder([saveDir 'Figures/']);
%     save_figure(fig, [figDir, grouped(exp).name ' fly spatial distribution across trials'],fig_type);
% end
% 
% 
% % % % This wil work if we have already accounted for the subtle differences in the
% % % distance measures due to the differences in the plate and lense configurations over
% % % time
% % % Can also look at the actual max distances from the center as an empirical value:
% % % (using eccentricity and making a histogram but separated by plate number) 
% % for exp = 1:num.exp
% %     g1 = data(exp).plate==1;
% %     g2 = data(exp).plate==2;
% % 
% %     fig = getfig('',1);
% %         y = grouped(exp).ecent.all(:,g1);
% %         histogram(y(:),'FaceColor',Color('magenta'))
% %         yyaxis right
% %         y = grouped(exp).ecent.all(:,g2);
% %         histogram(y(:),'FaceColor',Color('dodgerblue'))
% % 
% %         formatFig(fig, blkbgd)
% %         set(gca, 'ycolor', 'none')
% %         yyaxis left 
% %         set(gca, 'ycolor', 'none', 'box', 'off')
% %         xlabel('distance from arena center (mm)')
% % 
% %         % plot the 'elimination' lines
% %         v_line(innerR1,'magenta', '--')
% %         v_line(innerR2,'dodgerblue', '--')
% %         title(data(exp).ExpGroup,'color', Color('grey'))
% %         xlim([5,30])
% % 
% %         save_figure(fig, [saveDir, data(exp).ExpGroup ' eccentricity by plate types'],fig_type);
% % end
% 
% % % Compare the distance between the well centers across the two plate types since we
% % % know that they are equal in real life -- is there a subtle difference in the cam
% % % distance setting?
% % [plate1 plate2] = deal([]);
% % for exp = 1:num.exp
% %     for trial = 1:num.trial(exp)
% %         WC = (data(exp).data(trial).data.wellcenters);
% %         inter_dist = [pdist([WC(:,1),WC(:,3)]','euclidean'), pdist([WC(:,2),WC(:,4)]','euclidean')];
% %         switch data(exp).plate(trial)
% %             case 1
% %                 plate1 = [plate1, inter_dist];
% %             case 2
% %                 plate2 = [plate2, inter_dist];
% %         end
% %     end
% % end
% % disp(['Plate 1 mean: ' num2str(mean(plate1))])
% % disp(['Plate 2 mean: ' num2str(mean(plate2))])
% % 
% % bins = 460:5:480;
% % fig = figure;
% % hold on
% % histogram(plate1, bins, 'facecolor', Color('magenta'))
% % yyaxis right
% % histogram(plate2, bins,'facecolor', Color('dodgerblue'))
% % formatFig(fig, blkbgd)
% %         set(gca, 'ycolor', 'none')
% %         yyaxis left 
% %         set(gca, 'ycolor', 'none', 'box', 'off')
% %         xlabel('distance between food wells (pix)')
% % save_figure(fig, [saveDir, 'dist between wells by plate type'],fig_type);


%% ANALYSIS: Cluster eccentricity temperature
clearvars('-except',initial_vars{:})

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
            % Divided by heating / cooling
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            grouped(i).ecent.(g_name)(temp,1) = mean(mean(grouped(i).ecent.all(loc,:),2,'omitnan'),'omitnan'); %avg
            grouped(i).ecent.(g_name)(temp,2) = std(mean(grouped(i).ecent.all(loc,:),1,'omitnan'),'omitnan');%./num.trial(i); %err
        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        grouped(i).ecent.temp_all(temp,1) = mean(mean(grouped(i).ecent.all(loc,:),2,'omitnan'),'omitnan'); %avg
        grouped(i).ecent.temp_all(temp,2) = std(mean(grouped(i).ecent.all(loc,:),1,'omitnan'),'omitnan')./num.trial(i);% %err
    end
    grouped(i).ecent.temps = temps;
end


%% ANALYSIS: Create ROImasks -- logical masks for all the regions in the arena for each trial
clearvars('-except',initial_vars{:})
% Add a new structure that contains the masks for each of the regions of
% interest for the flies in this experiment which can be used to filter
% other properties later on...like speed in the outer ring etc.

% Initialize empty mask structure:
initial_vars{end+1} = 'ROImask';
initial_vars{end+1} = 'quadOrder';
ROImask = struct; 
fields = {'all', 'ring', 'inner75'}; 
food_fields = {'fullquad', 'innerquad', 'quadring', 'circle10', 'circle7', 'circle5'};
quadOrder = {'food', 'right', 'opp', 'left'}; 
for exp = 1:num.exp
    for trial = 1:num.trial(exp)
        % universal regions:
        for i = 1:length(fields)
            ROImask(exp).(fields{i})(trial).m = [];
        end
        % quadrant regions:
        for i = 1:length(food_fields)
            for q = 1:4
                ROImask(exp).(food_fields{i}).(quadOrder{q})(trial).m = [];
            end
        end
    end
end

for exp = 1:num.exp
    % find locations of flies within each region of the arena
    for trial = 1:num.trial(exp)
        % pull parameters for this trial:
        center = data(exp).data(trial).data.centre; % arena center
        x_loc = data(exp).data(trial).data.x_loc; % all fly X positions in the arena
        y_loc = data(exp).data(trial).data.y_loc; % all fly Y positions in the arena
        ct = data(exp).con_type(trial); % experiment lens configuration
        pix2mm = conversion(ct).pix2mm; % conversion to mm for this configuration
        R = conversion(ct).R; % arena radius
        circle75 = conversion(ct).circle75; % defines the distance to the inside edge of the outer ring
        circle10 = conversion(ct).circle10;
        circle7 = conversion(ct).circle7;
        circle5 = conversion(ct).circle5;
        foodWell = data(exp).T.foodLoc(trial); % location of food or randomly assigned region
        lowWell = grouped(exp).occ_idx(trial,1); % location of lowest occupancy region
        highWell = grouped(exp).occ_idx(trial,2); % location of highest occupancy region

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

        for CC = 1:3 % for food, low, high quad 
            switch CC
                case 1 % food well occupancy 
                    selectedWell = foodWell;
                    all_quads = true; % run the data for all the other quadrants
                    subfields = quadOrder;
                case 2 % low quad occupancy
                    selectedWell = lowWell;
                    all_quads = false; % don't run data on all other quadrants
                    subfields = {'low'};
                case 3 % high quad occupancy
                    selectedWell = highWell;
                    all_quads = false; % don't run data on all other quadrants
                    subfields = {'high'};
            end
            % Find the food quadrant (find location with the food well coordinates included)
            adjusted_wx = data(exp).data(trial).data.wellcenters(1,selectedWell) - center(1); % adjusted well position
            adjusted_wy = data(exp).data(trial).data.wellcenters(2,selectedWell) - center(2); 
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
    
            for i = 1:length(subfields) % for quadrant type (food, R, opp, L)
                idx = opLoc(i); % rearrange order based on trial food location
                well_idx = find(well_quads(idx).Mask);
                % well_idx = find(well_quad_loc==idx); % index of the well that falls into this quadrant of interest
                ROImask(exp).fullquad.(subfields{i})(trial).m = Q(idx).Mask; % full quadrant
                ROImask(exp).innerquad.(subfields{i})(trial).m = Q(idx).Mask & innerQuad; % quadrant inside inner circle
                ROImask(exp).quadring.(subfields{i})(trial).m = Q(idx).Mask & ~innerQuad; % quadrant portion of outer ring
    
                % well circle related distances
                dX = adjustedX - well_opt_x(well_idx);% well_idx
                dY = adjustedY - well_opt_y(well_idx);% well_idx
                well_D = hypot(dX,dY)./pix2mm;
                ROImask(exp).circle10.(subfields{(i)})(trial).m = well_D <= circle10; % within 10% area of well
                ROImask(exp).circle7.(subfields{(i)})(trial).m = well_D <= circle7; % within 7% area of well
                ROImask(exp).circle5.(subfields{(i)})(trial).m = well_D <= circle5; % within 5% area of well
            end
        end
    end
end

%% FIGURE: Visual demonstration of all the arena regions....
if strcmp('Yes', questdlg('Show demo image of arena regions of interest?','','Yes', 'No','Cance', 'No'))
    clearvars('-except',initial_vars{:})
    exp = 1;
    trial = 1;
    r = 4;
    c = 7;
    k = Color('Dodgerblue');
    s = 1;
    ct = data(exp).con_type(trial);
    R = conversion(ct).R*conversion(ct).pix2mm;
    center = data(exp).data(trial).data.centre;
    n = 28;
    foreColor = formattingColors(blkbgd);
    
    % regions: 
    roi = false([size(ROImask(exp).all(trial).m),n]);
    roi(:,:,1) = true(size(ROImask(exp).all(trial).m)); % all
    roi(:,:,8) = ROImask(exp).all(trial).m; % all screened
    roi(:,:,15) = ROImask(exp).ring(trial).m; % ring
    roi(:,:,22) = ROImask(exp).inner75(trial).m; % ring
    regionList = {'fullquad', 'innerquad','quadring', 'circle10','circle7','circle5'};
    for t = 2:length(regionList)+1
        q = 1;
        for i = t:7:n
            roi(:,:,i) = ROImask(exp).(regionList{t-1}).(quadOrder{q})(trial).m;
            q = q+1;
        end
    end
    
    fig = getfig;
    x = data(exp).data(trial).data.x_loc;
    y = data(exp).data(trial).data.y_loc;
    FW = data(exp).data(trial).data.wellcenters(:,data(exp).T.foodLoc(trial));
    % ALL points no filter
    for i = 1:n
        subplot(r,c,i); hold on
        scatter(x(roi(:,:,i)),y(roi(:,:,i)),s,k,'filled')
        viscircles(center',R,'color',foreColor);
        scatter(FW(1),FW(2),30, foreColor)
    end
    formatFig(fig,false,[r,c]);
    for i = 1:n
        subplot(r,c,i)
        set(gca,'xcolor', 'none', 'ycolor','none');
        axis square equal
    end
    save_figure(fig, [saveDir 'Figures/ROI outlines'],'-pdf',true,false);
    save_figure(fig, [saveDir 'Figures/ROI outlines'],'-png',true);
end

%% ANALYSIS: Extract occupancy data from the different regions 
clearvars('-except',initial_vars{:})

regionList = {'fullquad','innerquad', 'quadring', 'circle10', 'circle7', 'circle5'};
ext_quadOrder = [quadOrder, 'low', 'high'];

for exp = 1:num.exp
    % initialize empty variables
    [ring.all,inner75.all] = deal([]);
    region = struct;
    for rr = 1:length(regionList)
        for q = 1:length(ext_quadOrder) % each of the quadrants + low + high
            region(rr).(ext_quadOrder{q}).all = [];
            region(rr).(ext_quadOrder{q}).all_info = {'percent relative to all the flies in the whole arena'};
            if strcmp(regionList{rr},'innerquad') || strcmp(regionList{rr},'quadring')
                region(rr).(ext_quadOrder{q}).partial = [];
                region(rr).(ext_quadOrder{q}).partial_info = {'percent relative to only the flies in the four quadrants of this space'};
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
            for q = 1:length(ext_quadOrder) % each of the quadrants
                temp = ROImask(exp).(regionList{rr}).(ext_quadOrder{q})(trial).m;
                n = getPercentFlies(temp,nFull);
                region(rr).(ext_quadOrder{q}).all = autoCat(region(rr).(ext_quadOrder{q}).all,n,false);
                if strcmp(regionList{rr},'innerquad')
                    n = getPercentFlies(temp,nInner);
                    region(rr).(ext_quadOrder{q}).partial = autoCat(region(rr).(ext_quadOrder{q}).partial,n,false);
                elseif strcmp(regionList{rr},'quadring')
                    n = getPercentFlies(temp,nRing);
                    region(rr).(ext_quadOrder{q}).partial = autoCat(region(rr).(ext_quadOrder{q}).partial,n,false);
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
        for q = 1:length(ext_quadOrder) % each of the quadrants
            region(rr).(ext_quadOrder{q}).avg = mean(region(rr).(ext_quadOrder{q}).all,2,'omitnan');
            region(rr).(ext_quadOrder{q}).std = std(region(rr).(ext_quadOrder{q}).all,0,2,'omitnan');
            if strcmp(regionList{rr},'innerquad') || strcmp(regionList{rr},'quadring')
                region(rr).(ext_quadOrder{q}).partial_avg = mean(region(rr).(ext_quadOrder{q}).partial,2,'omitnan');
                region(rr).(ext_quadOrder{q}).parital_std = std(region(rr).(ext_quadOrder{q}).partial,0,2,'omitnan');
            end
        end
    end   
    
    % Assign the data to the grouped structure
    grouped(exp).ring = ring;
    grouped(exp).inner75 = inner75;
    for rr = 1:length(regionList)
        grouped(exp).(regionList{rr}) = region(rr);
    end
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
                nQ = length(ext_quadOrder);
                nP = 2;
            case {'fullquad','circle10','circle7','circle5'}
                nQ = length(ext_quadOrder);
                nP = 1;
        end
        for p = 1:nP % for the number of full or partial percentages to compare
            for q = 1:nQ % each of the quadrants
                % pull data for the right quadrant (if there are quadrants)
                if nQ>1
                    baseY = grouped(exp).(all_regions{rr}).(ext_quadOrder{q});
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
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.raw = raw_h;
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.avg = mean(raw_h, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.raw = raw_c;
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).temps = temps;
                        else
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.p_raw = raw_h;
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.p_avg = mean(raw_h, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.p_std = std(raw_h, 0, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.p_raw = raw_c;
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.p_avg = mean(raw_c, 2, 'omitnan');
                            grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.p_std = std(raw_c, 0, 2, 'omitnan');
                        end
                    case {'circle10','circle7','circle5'}
                        grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.raw = raw_h;
                        grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.avg = mean(raw_h, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.raw = raw_c;
                        grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                        grouped(exp).(all_regions{rr}).(ext_quadOrder{q}).temps = temps;
                end
            end
        end
    end
end

%% Demo test of the actual occupancy of the arena to see if the spatial constraints fit:
% clearvars('-except',initial_vars{:})
% R = 30; %mm
% innerR = R*sqrt(3/4); % radius of the inner 50% occupancy space R*sqrt(1/2)
% dist_from_edge = (R - innerR);
% 
% % 1) pull up image of the arena from one of the videos OR just go blank screen
% % 2) overlay all the positions of the flies during the experiments
% fig = getfig('',1); hold on
%     exp = 1;
%     trial = 1;
%     centre = data(exp).data(trial).data.centre; %distance to center of arena 
%     x = data(exp).data(trial).data.x_loc; % x locations for the entire experiment
%     y = data(exp).data(trial).data.y_loc; % x locations for the entire experiment
%     scatter(x(:),y(:),2,'w','filled')
%     formatFig(fig,true)
% 
%     viscircles(centre',R*pix2mm,'color', 'r','LineWidth',0.25);
%     viscircles(centre',innerR*pix2mm,'color', 'r','LineWidth',0.25);


%% ANALYSIS: Flies on Food
clearvars('-except',initial_vars{:})
maxR = 3; % 2.5 mm radius of the physical well -- give 0.5mm buffer zone outside well since body mark is body center not head
% TODO: (8.4.25) update this to include the null low and high well data

% Find the percent of the flies that are on the food 
for exp = 1:num.exp
    counts = []; ring_per = [];
    for trial = 1:num.trial(exp)
        pix2mm = conversion(data(exp).con_type(trial)).pix2mm;
        x = data(exp).data(trial).data.x_loc; % x locations for the entire experiment
        y = data(exp).data(trial).data.y_loc; % x locations for the entire experiment
        % find food well center and distance of flies to that point
        wellLoc  = data(exp).T.foodLoc(trial);
        centre = data(exp).data(trial).data.wellcenters(:,wellLoc);
        dX = x-centre(1);
        dY = y-centre(2);
        D = hypot(dX,dY); % distance from center of arena
        D = D./pix2mm;
        loc = D<=maxR ; % flies within 10% circle around the food
        circleCount = sum(loc,2);
        counts = autoCat(counts, circleCount, false); % count #flies in the outer ring
        ring_per = autoCat(ring_per,(circleCount./data(exp).T.NumFlies(trial)).*100,false); % convert to percent & combine
    end
    % pool the data
    grouped(exp).fliesonfood.counts = counts;
    grouped(exp).fliesonfood.all = ring_per;
    grouped(exp).fliesonfood.percent = ring_per; 
    grouped(exp).fliesonfood.avg = mean(ring_per,2,'omitnan');
end

% Pull the data together: 
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h,count_c,count_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
    
    % Update the averages for the classic temperature bins 
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        end
        raw_c(t,:) = mean(grouped(exp).fliesonfood.percent(c_frames,:),1,'omitnan');
        raw_h(t,:) = mean(grouped(exp).fliesonfood.percent(h_frames,:),1,'omitnan');
        count_c(t,:) = mean(grouped(exp).fliesonfood.counts(c_frames,:),1,'omitnan');
        count_h(t,:) = mean(grouped(exp).fliesonfood.counts(h_frames,:),1,'omitnan');

    end

    % find the avg and err and save to group structure
    grouped(exp).fliesonfood.increasing.raw = raw_h;
    grouped(exp).fliesonfood.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).fliesonfood.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).fliesonfood.decreasing.raw = raw_c;
    grouped(exp).fliesonfood.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).fliesonfood.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).fliesonfood.increasing.count_raw = count_h;
    grouped(exp).fliesonfood.increasing.count_avg = mean(count_h, 2, 'omitnan');
    grouped(exp).fliesonfood.increasing.count_err = std(count_h, 0, 2, 'omitnan');
    grouped(exp).fliesonfood.decreasing.count_raw = count_c;
    grouped(exp).fliesonfood.decreasing.count_avg = mean(count_c, 2, 'omitnan');
    grouped(exp).fliesonfood.decreasing.count_err = std(count_c, 0, 2, 'omitnan');
    grouped(exp).fliesonfood.temps = temps;
end

%% ANALYSIS: Generate and save the sleep quanitification
clearvars('-except',initial_vars{:})
paths = getPathNames;
nbins = 50;

% Create sleep data for unprocessed files (trial by trial)
for i = 1:num.exp
    
    TP = getTempTurnPoints(data(i).temp_protocol);
    fps = TP.fps;

    % How long does a fly need to be still to count as 'sleep'
    min_duration = 5*fps*60; % 5 mins * 3fps*60sec = data point number that must be met 'unmoving'
    
    T = data(i).T;
    for trial = 1:num.trial(i)
        trial_ID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
        sleep_file = [baseFolder paths.single_trial '/' trial_ID  '/' T.ExperimentID{trial} ' sleeping data v2.mat'];  
        
        if ~exist(sleep_file,"file")
            sleepData = struct;
            %preallocate for speed and space
            trial_length = length(data(i).data(trial).occupancy.time);
            [N,frameCount] = deal(nan(nbins,nbins,trial_length));
            sleepingCount = zeros(trial_length,1);

            % Set axis limits for the selected arena for bins that will dictate single fly sizes
            x = data(i).data(trial).data.centre(1);
            y = data(i).data(trial).data.centre(2);
            r = data(i).data(trial).data.r;
            xlimit = [x-(r+50),x+(r+50)];
            ylimit = [y-(r+50),y+50+r];

            % find the 'auto bin' lines (aka the bin edge coordinates)
            xedge = linspace(xlimit(1),xlimit(2),nbins+1); 
            yedge = linspace(ylimit(1),ylimit(2),nbins+1);

            % pull the fly locations during the trial
            x_loc = data(i).data(trial).data.x_loc;
            y_loc = data(i).data(trial).data.y_loc;
            % find the bins in which flies were present for this frame and add to the giant 'N' structure
            for frame = 1:trial_length
                X = x_loc(frame,:); % x-coordinates for all flies on this camera frame
                X(isnan(X)) = []; % removes any nan locations from the list
                Y = y_loc(frame,:); 
                Y(isnan(Y)) = [];
                N(:,:,frame) = histcounts2(X,Y,xedge,yedge); 
            end
            % (at this point, N is a large matrix where we will look for flies that
            % stay in the same x-y space bin for longer than 5 minutes (the def of
            % sleep in flies))

            % find grid space that have continuous occupation for more than min_duration
            frameCount(:,:,1) = N(:,:,1); % frameCount is the running count for # of frames that have a fly
            for frame = 2:trial_length
                currFrame = N(:,:,frame); % current frame locations
                resetLoc = currFrame==0; % locations that do not have flies and thus need a count reset

                tempCount = frameCount(:,:,frame-1)+currFrame; % add 1 to the frame count from previous frame count
                tempCount(resetLoc) = 0; % reset counts for spots with no flies

                frameCount(:,:,frame) = tempCount; % add current count into the saving structure
            end

            % ---- Vectorize the data (find the flies that are sleeping....) -----

            % pull coordinates of the food well for distance capture later
            foodWellLoc = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));
            c1 = foodWellLoc(1); % x-coordinate for the center of the food well
            c2 = foodWellLoc(2); % y-coordinate for the center of the food well
            
            % create empty matrixes for the x and y positions of the sleeping flies
            sleeping = struct;
            [sleeping.X, sleeping.Y, sleeping.all_distance] = deal(nan(trial_length,data(i).T.NumFlies(trial))); 
            sleeping.sleepNum = zeros(trial_length,1);
            [sleeping.dist_avg, sleeping.dist_err] = deal(nan(trial_length,1));

            % assign data by frame
            for frame = 1:trial_length
                frame_data = frameCount(:,:,frame) > min_duration; % puts a 'T' for any bin location that has had a stationary fly for >5 mins
                binLoc = find(frame_data>0); % vectorize the location % TODO HERE: resume sweep
                
                % Find the coordinates of the sleeping flies bins from the discretized data
                y_row = ceil(binLoc/nbins);
                x_row = rem(binLoc-1,nbins)+1;
                x_position = (xedge(x_row) + xedge(x_row+1))/2; %x-location is the middle of the x-bins
                y_position = (yedge(y_row) + yedge(y_row+1))/2; %y-location is the middle of the y-bins
                
                % add position data to the matrix:
                if ~isempty(binLoc)
                    % number of flies sleeping
                    sleepNum = length(x_position);
                    sleeping.sleepNum(frame) = sleepNum;
                    % location of sleeping flies
                    sleeping.X(frame,1:sleepNum) = x_position;
                    sleeping.Y(frame,1:sleepNum) = y_position;
                    % distance to food...
                    temp_dist = sqrt((x_position-c1).^2 + (y_position-c2).^2)./pix2mm;
                    sleeping.all_distance(frame,1:sleepNum) = temp_dist;
                    % average distance:
                    sleeping.dist_avg(frame) = mean(temp_dist);
                    sleeping.dist_err(frame) = std(temp_dist);
                end
            end

            save(sleep_file,'sleeping','-v7.3'); 
        end   
        disp([num2str(i) ' | ' num2str(trial)])
        clear N preallocate frameCount sleepingCount sleepLoc resetLoc tempCount
    end
    disp(['Done exp ' expNames{i}])
end
clearvars('-except',initial_vars{:})



%% ANALYSIS: Load previously created sleep data files and process data:
paths = getPathNames;
sleep = struct;
for i = 1:num.exp
    T = data(i).T;
    for trial = 1:num.trial(i)
        trial_ID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
        sleep_file = [baseFolder paths.single_trial '/' trial_ID  '/' T.ExperimentID{trial} ' sleeping data v2.mat'];  

        if exist(sleep_file,"file")
            load(sleep_file,'sleeping');
            sleep(i).trial(trial) = sleeping;
            clear sleeping
        else
            h = warndlg(['Warning: missing sleep data for ' expNames{i}]);
            uiwait(h)
        end
    end
end
disp('Loaded all sleep data')

% Process and prep the data for further analysis
for i = 1:num.exp
    [sleep(i).num, sleep(i).fract_sleep] = deal(zeros(length(grouped(i).temp),num.trial(i)));
    sleep(i).distance = nan(length(grouped(i).temp),num.trial(i));
    for trial = 1:num.trial(i)
        %number of sleeping flies
        inputdata = sleep(i).trial(trial).sleepNum;
        sleep(i).num(1:length(inputdata),trial) = inputdata;
        sleep(i).fract_sleep(1:length(inputdata),trial) = inputdata/data(i).T.NumFlies(trial);
        %distance to food for sleeping flies
        inputdata = sleep(i).trial(trial).dist_avg;
        sleep(i).distance(1:length(inputdata),trial) = inputdata;
    end
    sleep(i).sleepfract_avg = mean(sleep(i).fract_sleep,2,'omitnan');
    sleep(i).sleepfract_err = std(sleep(i).fract_sleep,0,2,'omitnan');
end

% Cluster the sleeping flies by temperature
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
            % fraction of flies sleeping
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            sleep(i).(g_name)(temp,1) = mean(mean(sleep(i).fract_sleep(loc,:),2,'omitnan'),'omitnan'); %avg 
            sleep(i).(g_name)(temp,2) = std(mean(sleep(i).fract_sleep(loc,:),1,'omitnan'),'omitnan');%./num.trial(i); %err
            % distance of sleeping flies
            sleep(i).([g_name '_dist'])(temp,1) = mean(mean(sleep(i).distance(loc,:),2,'omitnan'),'omitnan');
            sleep(i).([g_name '_dist'])(temp,2) = std(mean(sleep(i).distance(loc,:),1,'omitnan'),'omitnan');
        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        sleep(i).temp_all(temp,1) = mean(mean(sleep(i).fract_sleep(loc,:),2,'omitnan'),'omitnan'); %avg 
        sleep(i).temp_all(temp,2) = std(mean(sleep(i).fract_sleep(loc,:),1,'omitnan'),'omitnan');% %err
        %distance
        sleep(i).tempBinDist(temp,1) = mean(mean(sleep(i).distance(loc,:),2,'omitnan'),'omitnan');
        sleep(i).tempBinDist(temp,2) = std(mean(sleep(i).distance(loc,:),1,'omitnan'),'omitnan');
    end
    sleep(i).temps = temps;
end

% Thermal threat quantification and avg sleep quantity
for i = 1:num.exp
    tempProto = data(i).temp_protocol;
    if data(i).hold_exp
        temp_protocol  = 'linear_ramp_F_25-17';
    else
        temp_protocol = data(i).temp_protocol;
    end
    tPoints = getTempTurnPoints(temp_protocol); 
    fps = tPoints.fps; 

    if strcmp(data(i).temp_protocol,'Large_temp_sweep_15_35') ||...
       strcmp(data(i).temp_protocol,'Large_temp_sweep_15_35_FPS6') 
       sleep(i).avg_quant = nan;
       sleep(i).thermalThreat = nan;
       disp(['Skipping thermal threat for ' grouped(i).name])
       continue
    end

    demoRamp_idx = tPoints.down(1,1):tPoints.up(1,2); % points for a full ramp down and up

    % Thermal threat
    temp_ramp = grouped(i).temp(demoRamp_idx);
    thermalThreat = (sum(25-temp_ramp)/fps)/1000; % quant of time not at 25 *only punishes cold though...

    % Sleep duration

    nFlies = zeros(1,num.trial(i));
    nRamps = size(tPoints.up,1);
    for ramp = 1:nRamps
        rampROI = tPoints.down(ramp,1):tPoints.up(ramp,2);
        nFlies = nFlies + sum(sleep(i).num(rampROI,:),1);
    end
    avg_sleepDuration = ((nFlies/nRamps)./(data(i).T.NumFlies'))./fps; % avg fly sleep duration (in sec) during a single temp ramp
    % save data to sleep structure
    sleep(i).avg_quant = avg_sleepDuration;
    sleep(i).thermalThreat = thermalThreat;
end

initial_vars{end+1} = 'sleep';
initial_vars{end+1} = 'fps';
clearvars('-except',initial_vars{:})

%ANALYSIS: Sleep duration & start and stop of sleep
% Avg sleep duration
for i = 1:num.exp
    timing = struct;
    [sleepON,sleepOFF,sleepLength] = deal([]);
    for trial = 1:num.trial(i)
        dist_avg = sleep(i).trial(trial).all_distance;
        loc = diff(dist_avg)==0; % find frames where flies are stationary (where the distance to food does not change)
        duration_matrix = zeros(size(loc));
        duration_matrix = [zeros(1,size(loc,2)); duration_matrix]; %account for the first step where there is no sleep
        old_row = loc(1,:);
        for frame = 2:size(loc,1) % add frame count for each next moment of sleep
            curr_row = loc(frame,:);
            addLoc = ~(curr_row==0);
            duration_matrix(frame,addLoc) = old_row(addLoc) + curr_row(addLoc);
            old_row = duration_matrix(frame,:);
        end

        sleepDurationMatrix = duration_matrix./(fps*60); %sleep duration in minutes
        timing(trial).sleepduration = sleepDurationMatrix;
        % what is the max sleep length for these segments? 
        sleepDurationMatrix(sleepDurationMatrix==0)=nan;
%         fig = getfig('',1);
%         hold on
%         time = data(i).data(trial).occupancy.time;
%         for trial = 1:num.trial(i)
%             scatter(time, sleepDurationMatrix(:,trial),10)
%         end
%         xlabel('time (min)')
%         ylabel('sleep duration (min)')
%         formatFig(fig,true);
%         title(expNames{i},'color','w')
%         xlim([0,1000])

        % Find the 'end' of sleep duration and stop location (time)
        startLoc = (diff(diff(duration_matrix)))<0; % Find the start of sleep (time)
        stopLoc = (diff(duration_matrix))<0;
        total_duration = duration_matrix(stopLoc)/(fps*60); %duration in minutes

        [sleepStoppedFrame, sleepStartFrame] = deal([]); %when in time (frame) did the end of sleep occur?
        for fly = 1:size(stopLoc,2)
              sleepStoppedFrame = autoCat(sleepStoppedFrame,find(stopLoc(:,fly)),false);
              sleepStartFrame = autoCat(sleepStartFrame,find(startLoc(:,fly)),false);
        end
        % save individual start / stop sleep timing (for later alignment
        % with distance and location within the arena)
        sleep(i).trial(trial).sleepON = sleepStartFrame;
        sleep(i).trial(trial).sleepOFF = sleepStoppedFrame;
        sleep(i).trial(trial).sleepLoc = startLoc;

        % save trial information into larger matrix
        temp = sleepStartFrame(:);
        temp(isnan(temp)) = [];
        sleepON = autoCat(sleepON,temp,false);
        temp = sleepStoppedFrame(:);
        temp(isnan(temp)) = [];
        sleepOFF = autoCat(sleepOFF,temp,false);
        sleepLength = autoCat(sleepLength,total_duration,false);

    end

    % save into sleep structure
    sleep(i).sleepON = sleepON;
    sleep(i).sleepOFF = sleepOFF;
    sleep(i).sleepLength = sleepLength;
    sleep(i).durationMatrix = timing;

end

if ~exist([saveDir 'Sleep'],'dir')
    mkdir([saveDir 'Sleep'])
end

clearvars('-except',initial_vars{:})

%% Pull sleep data into the grouped structure format
for exp = 1:num.exp
    % add simple sleep data to grouped struct
    grouped(exp).sleep.counts = sleep(exp).num;
    grouped(exp).sleep.all = sleep(exp).fract_sleep*100;
    grouped(exp).sleep.percent = sleep(exp).fract_sleep*100; 
    grouped(exp).sleep.avg = sleep(exp).sleepfract_avg*100;
end

% Extract trial lines for warming and cooling
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
    
    % Update the averages for the classic temperature bins 
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        end
        raw_c(t,:) = mean(grouped(exp).sleep.percent(c_frames,:),1,'omitnan');
        raw_h(t,:) = mean(grouped(exp).sleep.percent(h_frames,:),1,'omitnan');
    end

    % find the avg and err and save to group structure
    grouped(exp).sleep.increasing.raw = raw_h;
    grouped(exp).sleep.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).sleep.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).sleep.decreasing.raw = raw_c;
    grouped(exp).sleep.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).sleep.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).sleep.temps = temps;
end



%% Sleep percent in quadrants, ring, and food circle...
% clearvars('-except',initial_vars{:})
% % R = 30; % mm
% % innerR = R*sqrt(3/4); % (old: radius of the inner 50% occupancy space R*sqrt(1/2))
% % dist_from_edge = (R - innerR);
% 
% for exp = 1:num.exp
%     quad_occ = [];
%     counts = [];
%     ring_per = [];
%     for trial = 1:num.trial(exp)
%         center = data(exp).data(trial).data.centre;
%         % find sleeping locations in the arena:
%         x_loc = sleep(exp).trial(trial).X;
%         y_loc = sleep(exp).trial(trial).Y;
% 
%         % QUADRANT:
%         r = data(exp).data(trial).data.r;
%         foodWell = data(exp).T.foodLoc(trial);
%         % Adjust the X and Y coordinates relative to the new center
%         adjustedX = x_loc - center(1);
%         adjustedY = y_loc - center(2);
%         % Initialize matrix to hold quadrant classification (same size as input matrices)
%         quadrantMatrix = zeros(size(x_loc));
%         % Define quadrant masks based on the new center
%         Q = [];
%         Q(1).Mask = (adjustedY > adjustedX) & (adjustedY <= -adjustedX);  % Top
%         Q(2).Mask = (adjustedY <= adjustedX) & (adjustedY <= -adjustedX); % Bottom
%         Q(3).Mask = (adjustedY <= adjustedX) & (adjustedY > -adjustedX);  % Left
%         Q(4).Mask = (adjustedY > adjustedX) & (adjustedY > -adjustedX);   % Right
%         % Determine the well locations (which determine quadrant assignment):
%         adjusted_wx = data(exp).data(trial).data.wellcenters(1,foodWell) - center(1);
%         adjusted_wy = data(exp).data(trial).data.wellcenters(2,foodWell) - center(2);   
%         idx_loc = false(1,4);
%         % Find the food quadrant (find location with the food well coordinates included)
%         idx_loc(1) = (adjusted_wy > adjusted_wx) & (adjusted_wy <= -adjusted_wx);  % top
%         idx_loc(2) = (adjusted_wy <= adjusted_wx) & (adjusted_wy <= -adjusted_wx); % right
%         idx_loc(3) = (adjusted_wy <= adjusted_wx) & (adjusted_wy > -adjusted_wx);  % bottom
%         idx_loc(4) = (adjusted_wy > adjusted_wx) & (adjusted_wy > -adjusted_wx);   % left
%         quad_loc = find(idx_loc);
%         fly_loc = ~isnan(x_loc); %gives logical for all fly positions in the position matrix
%         foodQuad = Q(quad_loc).Mask & fly_loc; % flies in food quad
%         nflies = sum(~isnan(x_loc),2);
%         y = (sum(foodQuad,2)./nflies).*100;
%         quad_occ = autoCat(quad_occ, y,false);
% 
%         % OUTER RING
%         D = sqrt(((x_loc-center(1)).^2 + (y_loc-center(2)).^2)); %distance from center of arena
%         D = D./pix2mm;
%         loc = D<=R & D>=innerR; % find the locations that are between edge and inner R
%         ringCount = sum(loc,2);
%         counts = autoCat(counts, ringCount, false); %count #flies in the outer ring
%         ring_per = autoCat(ring_per,(ringCount./nflies).*100,false); % convert to percent & combine
%     end
% 
%     % save quad occ to structure
%     grouped(exp).sleep.quad.all = quad_occ;
%     grouped(exp).sleep.quad.avg = mean(quad_occ,2,'omitnan');
%     grouped(exp).sleep.quad.std = std(quad_occ,0,2,'omitnan');
% 
%     % save ring occ to struct
%     grouped(exp).sleep.ring.counts = counts;
%     grouped(exp).sleep.ring.all = ring_per;
%     grouped(exp).sleep.ring.percent = ring_per; 
%     grouped(exp).sleep.ring.avg = mean(ring_per,2,'omitnan');
% end

%% SLEEPING analysis -- find the locations of the sleeping flies and create a sleep mask
clearvars('-except',initial_vars{:})
% Initialize empty sleep mask structure:
initial_vars{end+1} = 'sROImask';
sROImask = struct;  

% Initialize empty mask structure:
fields = {'all','ring', 'inner75'}; 
food_fields = {'fullquad','innerquad','quadring','circle10','circle7','circle5'};
ext_quadOrder = [quadOrder, 'low', 'high']; 

for exp = 1:num.exp
    for trial = 1:num.trial(exp)
        % universal regions:
        for i = 1:length(fields)
            sROImask(exp).(fields{i})(trial).m = [];
        end
        % quadrant regions + low/high regions:
        for i = 1:length(food_fields)
            for q = 1:length(ext_quadOrder)
                sROImask(exp).(food_fields{i}).(ext_quadOrder{q})(trial).m = [];
            end
        end
    end
end

for exp = 1:num.exp
    % find locations of flies within each region of the arena
    for trial = 1:num.trial(exp)
        % pull parameters for this trial:
        center = data(exp).data(trial).data.centre;
        % TODO : working here 5.30 updating the structure to work for sleep
        % location data
        x_loc = sleep(exp).trial(trial).X;
        y_loc = sleep(exp).trial(trial).Y;

        ct = data(exp).con_type(trial); % experiment lens configuration
        pix2mm = conversion(ct).pix2mm;
        R = conversion(ct).R;
        circle75 = conversion(ct).circle75; % defines the distance to the inside edge of the outer ring
        circle10 = conversion(ct).circle10;
        circle7 = conversion(ct).circle7;
        circle5 = conversion(ct).circle5;
        foodWell = data(exp).T.foodLoc(trial);
        lowWell = grouped(exp).occ_idx(trial, 1); % low well occupancy region
        highWell = grouped(exp).occ_idx(trial, 2); % low well occupancy region

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
        sROImask(exp).all(trial).nflies = nTotal;
        sROImask(exp).ring(trial).nflies = nTotalOuterROI;
        sROImask(exp).inner75(trial).nflies = nTotalInnerROI;

        % ------- Find the flies in the outer ring and inner region -------
        % Define quadrant masks based on the new center & excluding flies that are in the outer ring:
        % food_fields = {'circle10','circle7','circle5'};
        sROImask(exp).all(trial).m = fly_loc; % ALL FLIES IN FULL ARENA (with bounds)
        sROImask(exp).ring(trial).m = fly_loc & ~innerQuad; % ALL FLIES WITHIN THE OUTER RING
        sROImask(exp).inner75(trial).m = fly_loc & innerQuad; % ALL FLIES IN INNER CIRCLE

        for CC = 1:3
            switch CC
                case 1 % food well
                    subfields = quadOrder;
                    well_selection = foodWell;
                case 2 % low occupancy well
                    subfields = {'low'};
                    well_selection = lowWell;
                case 3 % high occupancy well
                    subfields = {'high'};
                    well_selection = highWell;
            end

            % Find the food quadrant (find location with the food well coordinates included)
            adjusted_wx = data(exp).data(trial).data.wellcenters(1,well_selection) - center(1); % adjusted well position
            adjusted_wy = data(exp).data(trial).data.wellcenters(2,well_selection) - center(2); 
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
    
            for i = 1:length(subfields) % for quadrant type (food, R, opp, L)
                idx = opLoc(i); % rearrange order based on trial food location
                well_idx = find(well_quads(idx).Mask);
                % well_idx = find(well_quad_loc==idx); % index of the well that falls into this quadrant of interest
                sROImask(exp).fullquad.(subfields{i})(trial).m = Q(idx).Mask; % full quadrant
                sROImask(exp).innerquad.(subfields{i})(trial).m = Q(idx).Mask & innerQuad; % quadrant inside inner circle
                sROImask(exp).quadring.(subfields{i})(trial).m = Q(idx).Mask & ~innerQuad; % quadrant portion of outer ring
    
                % well circle related distances
                dX = adjustedX - well_opt_x(well_idx);% well_idx
                dY = adjustedY - well_opt_y(well_idx);% well_idx
                well_D = hypot(dX,dY)./pix2mm;
                sROImask(exp).circle10.(subfields{(i)})(trial).m = well_D <= circle10; % within 10% area of well
                sROImask(exp).circle7.(subfields{(i)})(trial).m = well_D <= circle7; % within 7% area of well
                sROImask(exp).circle5.(subfields{(i)})(trial).m = well_D <= circle5; % within 5% area of well
            end
        end
    end
end


%% SLEEP ANALYSIS: Find fly sleeping locations within the regions of interest (takes about 20 seconds)
clearvars('-except',initial_vars{:})

regionList = {'fullquad','innerquad', 'quadring', 'circle10', 'circle7', 'circle5'};
groupType = {'allflies','sleepflies', 'withinROI'};
ext_quadOrder = [quadOrder, 'low', 'high'];

disp('Processing sleep locations...')
tic
for exp = 1:num.exp
    % initialize empty variables
    [ring.allflies.all,ring.sleepflies.all, ring.withinROI.all] = deal([]);
    [inner75.allflies.all,inner75.sleepflies.all, inner75.withinROI.all] = deal([]);
    [fullarena.allflies.all,fullarena.sleepflies.all] = deal([]);

    region = struct;
    for rr = 1:length(regionList)
        for q = 1:length(ext_quadOrder) % each of the quadrants
            region(rr).(ext_quadOrder{q}).allflies.all = [];
            region(rr).(ext_quadOrder{q}).allflies.all_info = {['percent of sleeping flies in the region relative ...' ...
                                        'to all the flies (not just sleeping flies) in the whole arena']};
            region(rr).(ext_quadOrder{q}).sleepflies.all = [];
            region(rr).(ext_quadOrder{q}).sleepflies.all_info = {'in region percent of the total flies that are sleeping'};
            region(rr).(ext_quadOrder{q}).withinROI.all = [];
            region(rr).(ext_quadOrder{q}).withinROI.all_info = {'in region sleep percent of the total flies in that region'};
            if strcmp(regionList{rr},'innerquad') || strcmp(regionList{rr},'quadring')
                region(rr).(ext_quadOrder{q}).allflies.partial = [];
                region(rr).(ext_quadOrder{q}).allflies.partial_info = {'percent relative to only the total flies in the four quadrants of this space'};
                region(rr).(ext_quadOrder{q}).sleepflies.partial = [];
                region(rr).(ext_quadOrder{q}).sleepflies.partial_info = {'percent relative to only the sleeping flies in the four quadrants of this space'};
            end
        end
    end
    
    % Fill structures with all the trials processed data
    for trial = 1:num.trial(exp)
        % % locations of flies that are sleeping in this trial:

        % Total numbers within each type of group: (denomenators) 
        nFull = ROImask(exp).all(trial).nflies;
        nInner = ROImask(exp).inner75(trial).nflies;
        nRing = ROImask(exp).ring(trial).nflies;
        snFull = sROImask(exp).all(trial).nflies; % s for sleeping
        snInner = sROImask(exp).inner75(trial).nflies;
        snRing = sROImask(exp).ring(trial).nflies;

        % sleeping fly percentages in the arena as a whole: 
        temp = sROImask(exp).all(trial).m;
        n = getPercentFlies(temp,nFull);  %total sleeping out of all flies 
        fullarena.allflies.all = autoCat(fullarena.allflies.all,n,false);
        n = getPercentFlies(temp,snFull); % total sleeping out of sleeping (should be 100 here...)
        fullarena.sleepflies.all = autoCat(fullarena.sleepflies.all,n,false);
    
        % sleeping fly percentages in the outer ring 
        temp = sROImask(exp).ring(trial).m; %logical indx of flies sleeping in the outer ring
        n = getPercentFlies(temp,nFull);
        ring.allflies.all = autoCat(ring.allflies.all,n,false);
        n = getPercentFlies(temp,snFull);
        ring.sleepflies.all = autoCat(ring.sleepflies.all,n,false);
        n = getPercentFlies(temp,nRing);
        ring.withinROI.all = autoCat(ring.withinROI.all,n,false);

        % sleeping flies in the inner circle
        temp = sROImask(exp).inner75(trial).m;
        n = getPercentFlies(temp,nFull);
        inner75.allflies.all = autoCat(inner75.allflies.all,n,false);
        n = getPercentFlies(temp,snFull);
        inner75.sleepflies.all = autoCat(inner75.sleepflies.all,n,false);
        n = getPercentFlies(temp,nInner);
        inner75.withinROI.all = autoCat(inner75.withinROI.all,n,false);

        % fig = figure; set(fig,'windowstyle', 'docked'); hold on
        % y = mean(inner75.sleepflies.all,2,'omitnan');
        % plot(y)
        % y = mean(ring.sleepflies.all,2,'omitnan');
        % plot(y)

        % fig = figure; set(fig,'windowstyle', 'docked');hold on
        %     plot(inner75.allflies.all,'color', Color('teal'))
        %     plot(ring.allflies.all,'color', Color('gold'))
        %     plot(inner75.all + ring.all,'color', Color('white'))
        %     formatFig(fig,true);
        %     ylabel('% flies')
        %     legend({'inner', 'ring','total'},'textcolor', 'w','box', 'off');
        %     set(gca, 'xcolor', 'none')
        %     ylim([0,100])

        % for each of the quadrant-related regions:
        for rr = 1:length(regionList)
            for q = 1:length(ext_quadOrder) % each of the quadrants
                temp = sROImask(exp).(regionList{rr}).(ext_quadOrder{q})(trial).m;
                for g = 1:3 % for each percentage type (e.g. all flies, within sleeping flies or within the region
                    switch g
                        case 1 % all flies
                            n = getPercentFlies(temp,nFull);
                        case 2 % within sleeping flies
                            n = getPercentFlies(temp,snFull);
                        case 3 % within the flies of that region
                            nArea = sum(ROImask(exp).(regionList{rr}).(ext_quadOrder{q})(trial).m,2);
                            n = getPercentFlies(temp,nArea);
                    end
                    
                    region(rr).(ext_quadOrder{q}).(groupType{g}).all = autoCat(region(rr).(ext_quadOrder{q}).(groupType{g}).all,n,false);
                end
                if strcmp(regionList{rr},'innerquad')
                    n = getPercentFlies(temp,nInner);
                    region(rr).(ext_quadOrder{q}).allflies.partial = autoCat(region(rr).(ext_quadOrder{q}).allflies.partial,n,false);
                    n = getPercentFlies(temp,snInner);
                    region(rr).(ext_quadOrder{q}).sleepflies.partial = autoCat(region(rr).(ext_quadOrder{q}).sleepflies.partial,n,false);
                elseif strcmp(regionList{rr},'quadring')
                    n = getPercentFlies(temp,nRing);
                    region(rr).(ext_quadOrder{q}).allflies.partial = autoCat(region(rr).(ext_quadOrder{q}).allflies.partial,n,false);
                    n = getPercentFlies(temp,snRing);
                    region(rr).(ext_quadOrder{q}).sleepflies.partial = autoCat(region(rr).(ext_quadOrder{q}).sleepflies.partial,n,false);
                end
            end
        end
    end

    % Pull out the averages and errors across the trials: 
    for g = 1:3
        ring.(groupType{g}).avg = mean(ring.(groupType{g}).all,2,'omitnan');
        ring.(groupType{g}).std = std(ring.(groupType{g}).all,0,2,'omitnan');
        inner75.(groupType{g}).avg = mean(inner75.(groupType{g}).all,2,'omitnan');
        inner75.(groupType{g}).std = std(inner75.(groupType{g}).all,0,2,'omitnan');
    end
    for rr = 1:length(regionList)
        for q = 1:length(ext_quadOrder) % each of the quadrants
            for g = 1:3
                region(rr).(ext_quadOrder{q}).(groupType{g}).avg = mean(region(rr).(ext_quadOrder{q}).(groupType{g}).all,2,'omitnan');
                region(rr).(ext_quadOrder{q}).(groupType{g}).std = std(region(rr).(ext_quadOrder{q}).(groupType{g}).all,0,2,'omitnan');
            end
            if strcmp(regionList{rr},'innerquad') || strcmp(regionList{rr},'quadring')
                region(rr).(ext_quadOrder{q}).allflies.partial_avg = mean(region(rr).(ext_quadOrder{q}).allflies.partial,2,'omitnan');
                region(rr).(ext_quadOrder{q}).allflies.parital_std = std(region(rr).(ext_quadOrder{q}).allflies.partial,0,2,'omitnan');
                region(rr).(ext_quadOrder{q}).sleepflies.partial_avg = mean(region(rr).(ext_quadOrder{q}).sleepflies.partial,2,'omitnan');
                region(rr).(ext_quadOrder{q}).sleepflies.parital_std = std(region(rr).(ext_quadOrder{q}).sleepflies.partial,0,2,'omitnan');
            end
        end
    end   
    
    % Assign the data to the grouped structure
    grouped(exp).sleep.ring = ring;
    grouped(exp).sleep.inner75 = inner75;
    grouped(exp).sleep.fullarena = fullarena;
    for rr = 1:length(regionList)
        grouped(exp).sleep.(regionList{rr}) = region(rr);
    end
end
toc
disp('...done')

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
all_regions = [regionList, 'ring', 'inner75', 'fullarena'];
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;

    for rr = 1:length(all_regions) % each type of region (e.g. outer ring, food circle etc)
      for g = 1:3
        switch all_regions{rr}
            case {'ring','inner75','fullarena'}
                nQ = 1; % quadrant number...
                nP = 1; % partial region percentages to check
            case {'innerquad','quadring'}
                nQ = length(ext_quadOrder);
                nP = 2;
            case {'fullquad','circle10','circle7','circle5'}
                nQ = length(ext_quadOrder);
                nP = 1;
        end
        for p = 1:nP % for the number of full or partial percentages to compare
            for q = 1:nQ % each of the quadrants
                % skip within region for the partial regino percentages
                if p==2 && g==3
                    continue
                end
                % skip within ROI for full arena (this doesn't exist)
                if g==3 && strcmp('fullarena',  all_regions{rr})
                    continue
                end
                % pull data for the right quadrant (if there are quadrants)
                if nQ>1
                    baseY = grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g});
                else 
                    baseY = grouped(exp).sleep.(all_regions{rr}).(groupType{g});
                end
                % pull the right type of data to run
                if p==1
                    y = baseY.all;
                else
                    y = baseY.partial;
                end
                
                % initialize empty structures for the variables
                [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp

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
                    case {'ring','inner75','fullarena'}    
                        grouped(exp).sleep.(all_regions{rr}).(groupType{g}).increasing.raw = raw_h;
                        grouped(exp).sleep.(all_regions{rr}).(groupType{g}).increasing.avg = mean(raw_h, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(groupType{g}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(groupType{g}).decreasing.raw = raw_c;
                        grouped(exp).sleep.(all_regions{rr}).(groupType{g}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(groupType{g}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(groupType{g}).temps = temps;
                    case {'fullquad','innerquad','quadring'}
                        if p == 1
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.raw = raw_h;
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.avg = mean(raw_h, 2, 'omitnan');
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.raw = raw_c;
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).temps = temps;
                        else
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.p_raw = raw_h;
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.p_avg = mean(raw_h, 2, 'omitnan');
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.p_std = std(raw_h, 0, 2, 'omitnan');
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.p_raw = raw_c;
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.p_avg = mean(raw_c, 2, 'omitnan');
                            grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.p_std = std(raw_c, 0, 2, 'omitnan');
                        end
                    case {'circle10','circle7','circle5'}
                        grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.raw = raw_h;
                        grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.avg = mean(raw_h, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).increasing.std = std(raw_h, 0, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.raw = raw_c;
                        grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.avg = mean(raw_c, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).decreasing.std = std(raw_c, 0, 2, 'omitnan');
                        grouped(exp).sleep.(all_regions{rr}).(ext_quadOrder{q}).(groupType{g}).temps = temps;
                end
            end
        end
      end
    end
end


%% ANALYSIS: Event aligned signals  -- food occupancy
% TODO: update this to work with the new occupancy structure format!
% 5.30.25

% TODO: ERROR -- this structure does not work for LTS because the two ramps
% are not identical and thus cannot be averaged in the same way ... 
clearvars('-except',initial_vars{:})
data_types = {'occ', 'dist', 'speed', 'ring','fullquad', 'innerquad','circle10'};

% account for any temperature hold groups
fictive_proto = false;
if any([data(:).hold_exp])
     result = questdlg(['There are temperature hold experiments, '...
         'do you want to assign them a fictive temperature protocol for '...
         'temperature event aligned comparisons?']);
    if strcmp(result, 'Yes') % offer a selection from the other trials within the grouping
         % Import and Load Temperature Protocol Options
        drivepath = getCloudPath;
        drivepath = drivepath(1:end-5);
        xlFile = [drivepath 'Quad Bowl Experiments.xlsx'];
        [~,~,excelfile] = xlsread(xlFile, 'Temp Protocols'); %#ok<XLSRD> % suppress the warning
        protocol_options = excelfile(:,1);
        idx  = listdlg('PromptString','Select fictive temp protocol','ListString',protocol_options, 'ListSize',[225, 500]);
        fictive_protocol = protocol_options{idx};
        fictive_proto = true;
        % OLD VERSON: only select from temp protocols within the grouped data sets
        % protocol_options = unique({data(:).temp_protocol});
        % idx  = listdlg('PromptString','Select fictive temp protocol','ListString',protocol_options, 'ListSize',[175, 150]);
        % fictive_protocol = protocol_options{idx};
        % fictive_proto = true;
    end
end

for i = 1:num.exp
    if fictive_proto && data(i).hold_exp
        grouped(i).aligned.temp_protocol = fictive_protocol;
        grouped(i).fictivetemp = fictive_protocol;
    else
        grouped(i).aligned.temp_protocol = data(i).temp_protocol;
        grouped(i).fictivetemp = data(i).temp_protocol;
    end
    tp = getTempTurnPoints(grouped(i).aligned.temp_protocol);
    
    % skip this alignment if this is the LTS since there isn't symmetry for
    % this type of temperature alignment
    if strcmp(grouped(i).fictivetemp,'Large_temp_sweep_15_35')
        continue
    end

    sectIDX = false(1,3);
    sections = {'decreasing', 'increasing', 'holding'};
    if tp.nDown >= 1
       sectIDX(1) = true;
    end
    if tp.nUp >= 1
        sectIDX(2) = true;
    end
    if tp.nHold >= 1
        sectIDX(3) = true;
    end
    % if tp.nHold>0 % account for protocols with no holding period
    %     sections = {'decreasing','increasing','holding'};
    % else
    %     sections = {'decreasing','increasing'};
    % end
    for ss = 1:3
        if ~sectIDX(ss) % skip groups that don't have a section
            continue
        end
        switch sections{ss}
            case 'increasing'
                tpBin = 'up';
                nrr = tp.nUp;
            case 'decreasing'
                tpBin = 'down';
                nrr = tp.nDown;
            case 'holding'
                tpBin = 'hold';
                nrr = tp.nHold;
        end
        [speed, occ, dist, ring, temperature,fullquad, innerquad, circle10] = deal([]);
        duration = min(tp.(tpBin)(:,2)-tp.(tpBin)(:,1)); % get shortest ramp period
        time = 1:duration;
        time = time./(tp.fps*60);
        % 'fullquad', 'innerquad','circle10'
        for rr = 1:nrr % for each ramp...
            % ROI and temperaure:
            ROI = tp.(tpBin)(rr,1):tp.(tpBin)(rr,1)+duration;
            temperature(:,rr) = grouped(i).temp(ROI);
            % food occupancy:
            occ(:,:,rr) = grouped(i).occ.all(ROI,:).*100;
            % distance to food: 
            dist(:,:,rr) = grouped(i).dist.all(ROI,:);
            % speed
            speed(:,:,rr) = grouped(i).speed.all(ROI,:);
            % ring occupancy
            ring(:,:,rr) = grouped(i).ring.all(ROI,:);
            % full quadrant occupancy
            fullquad(:,:,rr) = grouped(i).fullquad.food.all(ROI,:);
            % inner quadrant occupancy
            innerquad(:,:,rr) = grouped(i).innerquad.food.all(ROI,:);
            % 10% circle occupancy
            circle10(:,:,rr) = grouped(i).circle10.food.all(ROI,:);
        end
        % loop through all the data types to extract
        for type = 1:length(data_types)
            switch data_types{type}
                case 'dist' % DISTANCE: 
                    y_all = dist;
                case 'occ' % FOOD OCCUPANCY
                    y_all = occ;
                case 'speed' % SPEED
                    y_all = speed;
                case 'ring' % RING OCCUPANCY
                    y_all = ring;
                case 'fullquad'
                    y_all = fullquad;
                case 'innerquad'
                    y_all = innerquad;
                case 'circle10'
                    y_all = circle10;
            end
            % calculate the cross-ramp means and normalization
            y = mean(y_all,3);        
            y_norm = y_all-mean(y_all(1:10,:,:),'omitnan'); %normalize to zero distance
            y_norm = mean(y_norm,3); %this is now the avg over the different temp ramps within a trial
            % add to the grouped data
            grouped(i).aligned.(sections{ss}).(data_types{type}).all = y; 
            grouped(i).aligned.(sections{ss}).(data_types{type}).mean = mean(y,2,'omitnan');
            grouped(i).aligned.(sections{ss}).(data_types{type}).std = std(y,0,2,'omitnan');
            grouped(i).aligned.(sections{ss}).(data_types{type}).norm.all = y_norm; 
            grouped(i).aligned.(sections{ss}).(data_types{type}).norm.mean = mean(y_norm,2,'omitnan');
            grouped(i).aligned.(sections{ss}).(data_types{type}).norm.std = std(y_norm,0,2,'omitnan');
            grouped(i).aligned.(sections{ss}).temp = mean(temperature,2,'omitnan');
            grouped(i).aligned.(sections{ss}).time = time;
        end
    end
end

close all

%% EXPORT: save specific parameter data from groupings for quick load and comparison later
% [takes a couple minutes to run]
clearvars('-except',initial_vars{:})

% find out if they already exist
save_list = {'ring', 'inner75', 'fullquad', 'innerquad', 'quadring', 'circle10', 'circle7', 'circle5', 'fliesonfood', 'sleep'};
if isfile([saveDir, save_list{1} ' data.mat'])
    add_str = 'ps, data already exists';
else 
    add_str = 'ps, no existing data found';
end

if strcmp(questdlg({'Save subfields into separate data files for later use?'; add_str}),'Yes')
    tic
    for param = 1:length(save_list) % save each field as its own 
        y = struct;
        for exp = 1:num.exp
            y(exp).(save_list{param}) = grouped(exp).(save_list{param});
            y(exp).name = grouped(exp).name;
            y(exp).color = grouped(exp).color;
        end
        % Save the data from the single field to its own file: 
        save([saveDir save_list{param} ' data.mat'], 'y','-v7.3'); 
        disp(save_list{param})
    end
    toc
end

disp('All finished')

%% Food vs no food control trial information....
% TODO: 6.15 : update this to somehow account for the control vs food
% trials to have some kind of paired indicator value...


% dummy = [];
% plotData = [];
% for i = 1:4
%     a = grouped(np).fullquad.(quadOrder{i}).all;
%     dummy(i,:) = sum(a,1,'omitnan');
%     plotData(:,:,i) = a;
% end
% % find min and max occupancy quadrants: 
% [~, lowerIDX] = min(dummy);
% [~, upperIDX] = max(dummy);
% 
% [minOcc,maxOcc] = deal([]);
% for i = 1:num.trial(np)
%     minOcc(:,i) = squeeze(plotData(:,i,lowerIDX(i)));
%     maxOcc(:,i) = squeeze(plotData(:,i,upperIDX(i)));
% end
%     y_err = smooth(mean((minOcc-maxOcc)./2,2,'omitnan'),sSpan,'moving');
%     y_avg = smooth(mean([minOcc,maxOcc],2,'omitnan'),sSpan, 'moving');


%% VISUAL CHECK: how do the fly locations align with the plate ROIs

% conversion = getConversion;
% 
% dX = [];
% dY = [];
% for exp = 1:num.exp
%     for trial = 1:num.trial(exp)
%         dX = [dX; data(exp).data(trial).data.x_loc(:)];
%         dY = [dY; data(exp).data(trial).data.y_loc(:)];
%     end
% end
% 
% dX(isnan(dX)) = [];
% dY(isnan(dY)) = [];
% 
% fig = getfig('',1,[800,800]);
%     scatter(dX, dY, 1, 'w')
%     formatFig(fig,true)
% 
%     % overlay outer rim circles for each trial 
%     hold on
%     for exp = 1:num.exp
%         for trial = 1:num.trial(exp)
%             center = data(exp).data(trial).data.centre;
%             viscircles(center',conversion(2).R*pix2mm,'color', 'r','LineWidth',0.25);
%         end
%     end
% 
% % histogram to see why they are all not aligning??? 



%% old colors

%  case 'Berlin LTS 15-35 caviar vs empty'
%         expOrder = 1:2;
%         colors = {'DodgerBlue', 'Gray'};
%      case 'Berlin F Hold vs LRR 25-17 caviar'
%         expOrder = [4, 1:3, 5];
%         colors = {'DeepSkyBlue','DimGrey','Gray','DarkGrey','Gainsboro'};
% 
%     case 'Berlin vs UAS-Kir2.1 backcross F LRR 25-17 caviar'
%         expOrder = [1,4:6,2,3];
%         colors = {'DarkOrchid','LightSkyBlue','Blue','DeepSkyBlue','Yellow','Orange'};
% 
% % ---- WILD TYPE COMPARISONS ----
%     case 'WT linear recovery caviar'
%         expOrder = [];
%         expList = {'Berlin', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
%         % expList = {'Swedish', 'Berlin WT', 'OregonR','CantonS', 'Malawi', 'Zimbabwe'};
%         colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','Grey'};
%         for ii = 1:num.exp
%             try expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar']));
%             catch
%                 expOrder(ii) = find(strcmp(expNames,[expList{ii} ' S LRR 23-15 caviar']));
%             end
%         end
%     case 'WT linear recovery no food'
%         expOrder = [];
%         expList = {'Berlin WT', 'CantonS ', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'};
%         colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','Grey'};
%         for ii = 1:num.exp
%             expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp no food']));
%         end
%     case 'All linear recovery ramp'
%         expOrder = [];
%         expList = {'Berlin WT','CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe','IsoD1','W1118'};
%         colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','Grey','LightSalmon','Blue','DeepPink'};
%         for ii = 1:num.exp
%             expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar']));
%         end
%     case 'WT linear recovery caviar 27-19'
%         expOrder = [];
% %         expList = {'Berlin WT', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
% %         colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
%         expList = {'Berlin WT', 'OregonR', 'SwedishC', 'Malawi', 'Zimbabwe'}; %desired exp order
%         colors = {'DarkOrchid','LimeGreen','Red','Gold','Grey'};
%         for ii = 1:num.exp
%             expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar 27-19'])  | ...
%                                 strcmp(expNames,[expList{ii} ' linear recovery ramp 27-19 caviar']));
%         end
%     case 'WT linear recovery caviar 25-17'
%         expOrder = [];
% %         expList = {'Berlin WT', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
% %         colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
%         expList = {'Berlin WT', 'OregonR', 'SwedishC', 'Malawi', 'Zimbabwe'}; %desired exp order
%         colors = {'DarkOrchid','LimeGreen','Red','Gold','White'};
%         for ii = 1:num.exp
%             expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar 25-17'])  | ...
%                                 strcmp(expNames,[expList{ii} ' linear recovery ramp 25-17 caviar']));
%         end
%  % ---- FOOD VS NO FOOD CONTROLS ----
%     case 'Berlin linear recovery ramp food vs no food'
%         expOrder = [2,1];
%         colors = {'black','dodgerblue'}; %
%     case 'Berlin LRR 25-17 food vs no food'
%         expOrder = [1,2];
%         colors = {'white', 'gold'};
%     case 'Berlin F LRR 25-17 caviar vs no food'
%         expOrder = [1,2];
%         colors = {'teal', 'gold'};   
%     case 'Berlin S LRR 25-17 caviar vs no food'
%         expOrder = [2,1];
%         colors = {'white', 'gold'};
%     case 'Berlin S LRR 25-17 food vs no food'
%         expOrder = [2,1];
%         colors = {'white', 'gold'};
%     case 'Berlin LRR 25-17 new arena vs OG caviar'
%         expOrder = [1,2];
%         colors = {'teal', 'white'};
%     case 'Berlin F 25-17 food vs no food'
%         expOrder = [1,2];
%         colors = {'blueviolet', 'white'};
%     case 'Berlin giant ramp food vs no food'
%         expOrder = [1,2];
%         colors = {'grey','DarkOrchid'};
%     % case 'Berlin LTS 15-35 caviar vs empty'
%     %     expOrder = [2,1];
%     %     colors = {'Black','Teal'};
%     case 'CantonS LRR food vs no food'
%         expOrder = [1,2];
%         colors = {'white','DarkOrchid'};
%     case 'Berlin F hold 25C food vs no food'
%         expOrder = 1:2;
%         colors = {'slategrey', 'white'};
%     case 'Berlin LRR Temp-Shift Trials caviar vs empty'
%         expOrder = [2 3 5 6 4 1];
%         colors = {'Blue','DodgerBlue','LightSkyBlue','Plum','MediumPurple','BlueViolet'};
% 
% % ---- FOOD GROUP COMPARISIONS ----
%     case 'Berlin F LRR 25-17 ACV'
%         expOrder = 1:num.exp;
%         colors = {'White','magenta','dodgerblue','Orange','DarkOrchid'};
%     % case 'Berlin F LRR 25-17'
%     %     expOrder = [6,1,2,3,4,5];
%     %     colors = {'White','grey','yellow','orange', 'red','dodgerblue'};
%     case 'Berlin LRR 25-17 different food comp'
%         expOrder = [1,2];
%         colors = {'gold','DarkOrchid'};
%     case 'Berlin S LRR 23-15 food raised comparison'
%         expOrder = [1,2];
%         colors = {'gold','Dodgerblue'};
%     case 'Berlin S LRR 25-17 food raised comparison'
%         expOrder = [1,2];
%         colors = {'gold','Dodgerblue'};
%     case 'Berlin S LRR 23-15 caviar'
%         expOrder = [1,3,2];
%         colors = {'white','gold', 'DarkOrchid'};
%     case 'Berlin S LRR 25-17 caviar'
%         expOrder = [1,3,2]; % combined, old, new
%         colors = {'white','gold', 'DarkOrchid'};
%     case 'Berlin S LRR 27-19 caviar'
%         expOrder = [1,3,2]; % combined, old, new
%         colors = {'white','gold', 'DarkOrchid'};
%     % ---- SEX COMPARISONS ----
%     case 'Zimbabwe sex comparison'
%         expOrder = 1:num.exp;
%         colors = {'White','magenta','dodgerblue','Orange'};
%     case 'Berlin sex comparison'
%         expOrder = 1:num.exp;
%         colors = {'White','magenta','dodgerblue','Orange'};
% 
%     % ---- TEMP RATE COMPARISONS ----
%     case 'Berlln LRR 25-17 temprate comp caviar'
%         expOrder = 4:-1:1; 
%         colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};
%     case 'Berlin LRR 25-17 temprate caviar'
%         expOrder = [3 2 1 4 5]; % slow to fast
%         colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};
%     case 'Berlin LRR temprate comp'
%         expOrder = 4:-1:1; % slow to fast
%         colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue'};
%     case 'Berlin LRR 25-17 temp rate no food'
%         expOrder = [5 3 2 1 4]; % slow to fast
%         colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};
%     case 'Berlin LRR temprate comp no food'
%         expOrder = 4:-1:1; % slow to fast
%         colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue'};
%     case 'Zimbabwe LRR caviar temprate comparison'
%         expOrder = [3,4,2,1];
%         colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue'};
% 
%     % ---- TEMP SHIFT COMPARISONS ----
%     case 'Swedish LRR tempshift comp cavair full'
%         expOrder = 1:4; %lowest to highest
%         colors = {'dodgerblue','powderblue','peachpuff','tomato'};
%     case 'Swedish LRR temp shift caviar'
%         expOrder = [1,2,3]; %lowest to highest temp
%         colors = {'dodgerblue','gold','red'};
%     case 'Malawi LRR temp shift caviar'
%         expOrder = 1:3; %lowest to highest temp
%         colors = {'dodgerblue','gold','red'};
%     case 'OregonR LRR temp shift comp caviar'
%         expOrder = 1:3; %lowest to highest
%         colors = {'dodgerblue','gold','red'};
%     case 'Zimbabwe temp shifted linear recovery ramp caviar'
%         expOrder = [3,1,2]; %lowest to highest
%         colors = {'dodgerblue','gold','red'};
%     case 'Berlin temp shifted comp linear recovery ramp caviar'
%         expOrder = [3,1,2]; %lowest to highest
%         colors = {'dodgerblue','gold','red'};
%     case 'Berlin S LRR temp shift'
%         expOrder = [1,2,3]; %lowest to highest
%         colors = {'dodgerblue','gold','red'};
%     case 'Berlin S LRR TempShift caviar'
%         expOrder = 1:3;
%         colors = {'dodgerblue','gold','red'};
%     case 'Berlin LRR tempshift no food'
%         expOrder = [3 1 2];
%         colors = {'dodgerblue','gold','red'};
%     case 'LRR 27-19 no food comparison'
%         expOrder = 1:3;
%         colors = {'DarkOrchid', 'LimeGreen', 'Red'};
%     case 'Berlin WT mechanical manipulation linear recovery ramp caviar'
%         expOrder = [1,2]; %lowest to highest
%         colors = {'dodgerblue','peachpuff'};
%     case 'W1118 genetic background comp'
%         expOrder = [1,2,3]; %lowest to highest
%         colors = {'deeppink','dodgerblue','greenyellow'};
%     case 'R60H12-gal4 S LRR 25-17 food vs no food'
%         expOrder = 1:2;
%         colors = {'white', 'dodgerblue'};
%     case 'Berlin vs UAS-Chrimson-TM2 S LRR 25-17 caviar'
%         expOrder  = 1:2;
%         colors  = {'dodgerblue', 'yellow'};
%     case 'Berlin S LRR 25_17 caviar Cedar vs Collge'
%         expOrder  = 1:2;
%         colors  = {'dodgerblue', 'yellow'};
%     case 'Berlin F LRR caviar Cedar vs College'
%         expOrder  = 1:2;
%         colors  = {'dodgerblue', 'yellow'};
%     % ===  controls ======
%     case 'Berlin S LRR 23-15 all controls'
%         expOrder = [2 1 4 3]; %NF-NT; F-NT; NF-T; F-T
%         colors = {'lightslategray', 'white', 'lightpink', 'deeppink'};
%     case 'Berlin F LRR 25-17 controls and caviar'
%         expOrder = [4 3 2 1]; %NF-NT; F-NT; NF-T; F-T
%         colors = {'white', 'lightslategray', 'lightpink', 'deeppink'};
%     case 'Berlin Temp Holds caviar'
%         expOrder = 1:5;
%         colors = {'Blue', 'Dodgerblue', 'Cyan', 'lightcyan', 'white'}; %cold to warm
%     case 'Berlin Temp Holds no food'
%         expOrder = 1:4;
%         colors = { 'Blue','Dodgerblue', 'Cyan', 'lightcyan'}; %cold to warm
%     case 'Berlin 23C Hold food vs no food'
%         expOrder = [2,1];
%         colors = {'white', 'aquamarine'};
%     case 'Berlin S LRR 25-17 caviar vs temp holds'
%         colors = {'gold','lightcyan', 'dodgerblue'};
%         expOrder = 3:-1:1;
%     case 'Berlin S LRR 25-17 no food vs temp holds'
%         colors = {'gold','lightcyan', 'dodgerblue'};
%         expOrder = 3:-1:1;
%     case 'Berlin F LRR 25-17 food vs odor'
%         expOrder = 1:num.exp;
%         colors = {'DarkOrchid','Gold','dodgerblue','turquoise','lime','red','pink','Orange'};
% 
%         % ======  Sensory Component Experiments ======
%     case 'Berlin F LRR 25-17 antenna arista intact comparisons'
%         expOrder = 1:4; % empty, no antenna, no arista, intact
%         colors = { 'white','lightslategray', 'lightpink', 'deeppink'};
%     case 'Berlin F LRR 25-17 olfaction contribution'
%         expOrder = [1,2,5,4,3];
%         colors = {'gold','grey','dodgerblue','darkorchid','lime'};
%     case 'Berlin F LRR 25-17 sensory components'
%         expOrder = [2 3 4 5 1]; %full, water, waxed, 25, 17
%         colors = {'black','blue','orange', 'magenta', 'pink'};
%         % ============= SILENCING EXPERIMENTS ===============
%     case 'R77C10-gal4 F LRR caviar comparisons'
%         expOrder = 1:3; %
%         colors = {'dodgerblue','Gold','magenta'};
%     case 'IR25a-gal4;TM2 F LRR caviar comparisons'
%         expOrder = [1,3,2];
%         colors = {'darkorchid','gold', 'turquoise'};
%     case 'IR25a-gal4 F LRR caviar comparisons'
%         expOrder = [1,3,2,4,5,6]; %tm2, tm2 control, tb6b, tm6b control, kir control, test line
%         colors = {'cyan', 'dodgerblue', 'lightpink', 'deeppink', 'lightcyan', 'gold'};
% 



























   







