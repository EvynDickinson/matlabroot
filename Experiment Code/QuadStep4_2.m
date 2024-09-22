% RUN FULL SECTION
% DataTransferGUI

% SELECT AND LOAD DATA FOR ANAYLSIS
% CAVEATS:
% ** THIS ASSUMES ALL LOADED GROUPS HAVE THE SAME WITHIN-GROUPING TEMPERATURE PROTOCOLS
% *** DOESN'T WORK FOR TEMP PROTOCOLS WITH MORE THAN 1 HEATING AND COOLING TEMP RATE


%% Select data groups to compare
% add matlabroot folder to the directory path
% addpath(genpath('C:\matlabroot'));

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
            [dirIdx, v] = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[400,700]);
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
            expIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'multiple','ListSize',[300,450], ...
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
            expIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'multiple','ListSize',[300,450]);
            expNames = list_dirs(expIdx); %name of experiment groups selected
            num.exp = length(expIdx);  %number of groups selected
           
            % Load selected experiment data groups
            for i = 1:num.exp

                % get field list for loading data:
                dummy = load([structFolder expNames{i} '/' expNames{i} ' post 3.1 data.mat']);
                
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
            disp([exp_name ' needs to be rebuilt from Step 3.1'])
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
            dummy = load([structFolder exp_name '/' exp_name ' post 3.1 data.mat']);
             
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
        guessLoc = strcmp(expGroup,list_dirs); % TODO check if there is an existing folder and then offer that as the base selection
        dirIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[300,450], 'InitialValue',includedIdx,guessLoc);
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

disp(expNames')

%% ANALYSIS: organize data for each group
clearvars('-except',initial_vars{:})
fig_type = '-pdf'; 
blkbgd = false;
initial_vars = [initial_vars(:); 'initial_vars'; 'grouped'; 'expGroup'; 'saveDir'; 'mat';'expOrder'; 'fig_type';'f2m';'pix2mm';'blkbgd'];
initial_vars = unique(initial_vars);
% f2m = 3*60; %fps*min = number of frames in a minute -- this is now held
% in the getTempRate function
pix2mm = 12.8; % pixels per mm of real space
grouped = struct;

% Color selections
switch expGroup
    case 'Berlin LTS 15-35 caviar vs empty'
        expOrder = 1:2;
        colors = {'Gray', 'DodgerBlue'};
     case 'Berlin F Hold vs LRR 25-17 caviar'
        expOrder = [4, 1:3, 5];
        colors = {'DeepSkyBlue','DimGrey','Gray','DarkGrey','Gainsboro'};

    case 'Berlin vs UAS-Kir2.1 backcross F LRR 25-17 caviar'
        expOrder = [1,4:6,2,3];
        colors = {'DarkOrchid','LightSkyBlue','Blue','DeepSkyBlue','Yellow','Orange'};

% ---- WILD TYPE COMPARISONS ----
    case 'WT linear recovery caviar'
        expOrder = [];
        expList = {'Berlin', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
        % expList = {'Swedish', 'Berlin WT', 'OregonR','CantonS', 'Malawi', 'Zimbabwe'};
        colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','Grey'};
        for ii = 1:num.exp
            try expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar']));
            catch
                expOrder(ii) = find(strcmp(expNames,[expList{ii} ' S LRR 23-15 caviar']));
            end
        end
    case 'WT linear recovery no food'
        expOrder = [];
        expList = {'Berlin WT', 'CantonS ', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'};
        colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','Grey'};
        for ii = 1:num.exp
            expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp no food']));
        end
    case 'All linear recovery ramp'
        expOrder = [];
        expList = {'Berlin WT','CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe','IsoD1','W1118'};
        colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','Grey','LightSalmon','Blue','DeepPink'};
        for ii = 1:num.exp
            expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar']));
        end
    case 'WT linear recovery caviar 27-19'
        expOrder = [];
%         expList = {'Berlin WT', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
%         colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
        expList = {'Berlin WT', 'OregonR', 'SwedishC', 'Malawi', 'Zimbabwe'}; %desired exp order
        colors = {'DarkOrchid','LimeGreen','Red','Gold','Grey'};
        for ii = 1:num.exp
            expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar 27-19'])  | ...
                                strcmp(expNames,[expList{ii} ' linear recovery ramp 27-19 caviar']));
        end
    case 'WT linear recovery caviar 25-17'
        expOrder = [];
%         expList = {'Berlin WT', 'CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'}; %desired exp order
%         colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
        expList = {'Berlin WT', 'OregonR', 'SwedishC', 'Malawi', 'Zimbabwe'}; %desired exp order
        colors = {'DarkOrchid','LimeGreen','Red','Gold','White'};
        for ii = 1:num.exp
            expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar 25-17'])  | ...
                                strcmp(expNames,[expList{ii} ' linear recovery ramp 25-17 caviar']));
        end
 % ---- FOOD VS NO FOOD CONTROLS ----
    case 'Berlin linear recovery ramp food vs no food'
        expOrder = [2,1];
        colors = {'black','dodgerblue'}; %
    case 'Berlin LRR 25-17 food vs no food'
        expOrder = [1,2];
        colors = {'white', 'gold'};
    case 'Berlin S LRR 25-17 caviar vs no food'
        expOrder = [2,1];
        colors = {'white', 'gold'};
    case 'Berlin S LRR 25-17 food vs no food'
        expOrder = [2,1];
        colors = {'white', 'gold'};
    case 'Berlin LRR 25-17 new arena vs OG caviar'
        expOrder = [1,2];
        colors = {'teal', 'white'};
    case 'Berlin F 25-17 food vs no food'
        expOrder = [1,2];
        colors = {'blueviolet', 'white'};
    case 'Berlin giant ramp food vs no food'
        expOrder = [1,2];
        colors = {'grey','DarkOrchid'};
    case 'Berlin LTS 15-35 caviar vs empty'
        expOrder = [2,1];
        colors = {'Black','Teal'};
    case 'CantonS LRR food vs no food'
        expOrder = [1,2];
        colors = {'white','DarkOrchid'};
    case 'Berlin F hold 25C food vs no food'
        expOrder = 1:2;
        colors = {'slategrey', 'white'};

% ---- FOOD GROUP COMPARISIONS ----
    case 'Berlin F LRR 25-17 ACV'
        expOrder = 1:num.exp;
        colors = {'White','magenta','dodgerblue','Orange','DarkOrchid'};
    % case 'Berlin F LRR 25-17'
    %     expOrder = [6,1,2,3,4,5];
    %     colors = {'White','grey','yellow','orange', 'red','dodgerblue'};
    case 'Berlin LRR 25-17 different food comp'
        expOrder = [1,2];
        colors = {'gold','DarkOrchid'};
    case 'Berlin S LRR 23-15 food raised comparison'
        expOrder = [1,2];
        colors = {'gold','Dodgerblue'};
    case 'Berlin S LRR 25-17 food raised comparison'
        expOrder = [1,2];
        colors = {'gold','Dodgerblue'};
    case 'Berlin S LRR 23-15 caviar'
        expOrder = [1,3,2];
        colors = {'white','gold', 'DarkOrchid'};
    case 'Berlin S LRR 25-17 caviar'
        expOrder = [1,3,2]; % combined, old, new
        colors = {'white','gold', 'DarkOrchid'};
    case 'Berlin S LRR 27-19 caviar'
        expOrder = [1,3,2]; % combined, old, new
        colors = {'white','gold', 'DarkOrchid'};
    % ---- SEX COMPARISONS ----
    case 'Zimbabwe sex comparison'
        expOrder = 1:num.exp;
        colors = {'White','magenta','dodgerblue','Orange'};
    case 'Berlin sex comparison'
        expOrder = 1:num.exp;
        colors = {'White','magenta','dodgerblue','Orange'};

    % ---- TEMP RATE COMPARISONS ----
    case 'Berlln LRR 25-17 temprate comp caviar'
        expOrder = 4:-1:1; 
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};
    case 'Berlin LRR 25-17 temprate caviar'
        expOrder = [3 2 1 4 5]; % slow to fast
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};
    case 'Berlin LRR temprate comp'
        expOrder = 4:-1:1; % slow to fast
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue'};
    case 'Berlin LRR 25-17 temp rate no food'
        expOrder = [5 3 2 1 4]; % slow to fast
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue', 'dodgerblue'};
    case 'Berlin LRR temprate comp no food'
        expOrder = 4:-1:1; % slow to fast
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue'};
    case 'Zimbabwe LRR caviar temprate comparison'
        expOrder = [3,4,2,1];
        colors = {'Deeppink','Gold','MediumSpringGreen','mediumslateblue'};

    % ---- TEMP SHIFT COMPARISONS ----
    case 'Swedish LRR tempshift comp cavair full'
        expOrder = 1:4; %lowest to highest
        colors = {'dodgerblue','powderblue','peachpuff','tomato'};
    case 'Swedish LRR temp shift caviar'
        expOrder = [1,2,3]; %lowest to highest temp
        colors = {'dodgerblue','gold','red'};
    case 'Malawi LRR temp shift caviar'
        expOrder = 1:3; %lowest to highest temp
        colors = {'dodgerblue','gold','red'};
    case 'OregonR LRR temp shift comp caviar'
        expOrder = 1:3; %lowest to highest
        colors = {'dodgerblue','gold','red'};
    case 'Zimbabwe temp shifted linear recovery ramp caviar'
        expOrder = [3,1,2]; %lowest to highest
        colors = {'dodgerblue','gold','red'};
    case 'Berlin temp shifted comp linear recovery ramp caviar'
        expOrder = [3,1,2]; %lowest to highest
        colors = {'dodgerblue','gold','red'};
    case 'Berlin S LRR temp shift'
        expOrder = [1,2,3]; %lowest to highest
        colors = {'dodgerblue','gold','red'};
    case 'Berlin S LRR TempShift caviar'
        expOrder = 1:3;
        colors = {'dodgerblue','gold','red'};
    case 'Berlin LRR tempshift no food'
        expOrder = [3 1 2];
        colors = {'dodgerblue','gold','red'};
    case 'LRR 27-19 no food comparison'
        expOrder = 1:3;
        colors = {'DarkOrchid', 'LimeGreen', 'Red'};
    case 'Berlin WT mechanical manipulation linear recovery ramp caviar'
        expOrder = [1,2]; %lowest to highest
        colors = {'dodgerblue','peachpuff'};
    case 'W1118 genetic background comp'
        expOrder = [1,2,3]; %lowest to highest
        colors = {'deeppink','dodgerblue','greenyellow'};
    case 'R60H12-gal4 S LRR 25-17 food vs no food'
        expOrder = 1:2;
        colors = {'white', 'dodgerblue'};
    case 'Berlin vs UAS-Chrimson-TM2 S LRR 25-17 caviar'
        expOrder  = 1:2;
        colors  = {'dodgerblue', 'yellow'};
    case 'Berlin S LRR 25_17 caviar Cedar vs Collge'
        expOrder  = 1:2;
        colors  = {'dodgerblue', 'yellow'};
    case 'Berlin F LRR caviar Cedar vs College'
        expOrder  = 1:2;
        colors  = {'dodgerblue', 'yellow'};
    % ===  controls ======
    case 'Berlin S LRR 23-15 all controls'
        expOrder = [2 1 4 3]; %NF-NT; F-NT; NF-T; F-T
        colors = {'lightslategray', 'white', 'lightpink', 'deeppink'};
    case 'Berlin F LRR 25-17 controls and caviar'
        expOrder = [4 3 2 1]; %NF-NT; F-NT; NF-T; F-T
        colors = {'white', 'lightslategray', 'lightpink', 'deeppink'};
    case 'Berlin Temp Holds caviar'
        expOrder = 1:5;
        colors = {'Blue', 'Dodgerblue', 'Cyan', 'lightcyan', 'white'}; %cold to warm
    case 'Berlin Temp Holds no food'
        expOrder = 1:4;
        colors = { 'Blue','Dodgerblue', 'Cyan', 'lightcyan'}; %cold to warm
    case 'Berlin 23C Hold food vs no food'
        expOrder = [2,1];
        colors = {'white', 'aquamarine'};
    case 'Berlin S LRR 25-17 caviar vs temp holds'
        colors = {'gold','lightcyan', 'dodgerblue'};
        expOrder = 3:-1:1;
    case 'Berlin S LRR 25-17 no food vs temp holds'
        colors = {'gold','lightcyan', 'dodgerblue'};
        expOrder = 3:-1:1;
    case 'Berlin F LRR 25-17 food vs odor'
        expOrder = 1:num.exp;
        colors = {'DarkOrchid','Gold','dodgerblue','turquoise','lime','red','pink','Orange'};
 
        % ======  Sensory Component Experiments ======
    case 'Berlin F LRR 25-17 antenna arista intact comparisons'
        expOrder = 1:4; % empty, no antenna, no arista, intact
        colors = { 'white','lightslategray', 'lightpink', 'deeppink'};
    case 'Berlin F LRR 25-17 olfaction contribution'
        expOrder = [1,2,5,4,3];
        colors = {'gold','grey','dodgerblue','darkorchid','lime'};
    case 'Berlin F LRR 25-17 sensory components'
        expOrder = [2 3 4 5 1]; %full, water, waxed, 25, 17
        colors = {'black','blue','orange', 'magenta', 'pink'};
end

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

for i = 1:num.exp % FOR EACH DATA GROUP
    % GENERAL
    if ~isfield(data,'fps')
        data(i).fps = 3;
    end
    grouped(i).name = data(i).ExpGroup;
    if exist('colors','var')
        grouped(expOrder(i)).color = Color(colors{(i)});
    else 
        grouped(expOrder(i)).color = cMap(i,:);
    end
    % TIME COURSE DATA
    num.trial(i) = data(i).ntrials;
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


%% ANALYSIS: normalize fly position within arena 
clearvars('-except',initial_vars{:})
OG_Orientation = datetime('10.20.2023','InputFormat','MM.dd.yyyy');

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
        loc_mat(rr,tt).x = []; %set empty space for appending
        loc_mat(rr,tt).y = [];
    end
  end
  wellXY = [];
  for trial = 1:num.trial(i)
    trial_date = datetime(data(i).T.Date{trial},'InputFormat','MM.dd.yyyy');
    % get arena information
    well_loc = data(i).T.foodLoc(trial);
    wells = data(i).data(trial).data.wellcenters;
    % find offset to make the food well the origin
    x_offset = wells(1,well_loc);
    y_offset = wells(2,well_loc);
    wells_x = wells(1,:)-x_offset;
    wells_y = wells(2,:)-y_offset;

    X = data(i).data(trial).data.x_loc;
    Y = data(i).data(trial).data.y_loc;
    X = X-x_offset;
    Y = Y-y_offset;
    % save the position normalized data into the grouped structure or
    % something

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
  end
  
  % save trial data into broad structure
  grouped(i).position.loc = loc_mat;
  grouped(i).position.temp_rates = temp_rates;
  grouped(i).position.temp_list = temp_list;
  grouped(i).position.well_pos = wellXY;
end

%% ANALYSIS: calculate group occupancy
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

%% ANALYSIS: Calculate flies within the outer ring of the region
clearvars('-except',initial_vars{:})

R = data(1).data(1).data.r; % radius of the arena in pixels
innerR = R/sqrt(2); % radius of the inner 50% occupancy space
dist_from_edge = (R - innerR)/pix2mm;

% Find the percent of the flies that are in the outer ring
for exp = 1:num.exp
    counts = []; ring_per = [];
    for trial = 1:num.trial(exp)
        x = data(exp).data(trial).data.x_loc; % x locations for the entire experiment
        y = data(exp).data(trial).data.y_loc; % x locations for the entire experiment
        centre = data(exp).data(trial).data.centre; %distance to center of arena 
        D = sqrt(((x-centre(1)).^2 + (y-centre(2)).^2)); %distance from center of arena
        loc = D<=R & D>=innerR; % find the locations that are between edge and inner R
        ringCount = sum(loc,2);
        counts = autoCat(counts, ringCount, false); %count #flies in the outer ring
        ring_per = autoCat(ring_per,(ringCount./data(exp).T.NumFlies(trial)).*100,false); % convert to percent & combine
    end
    % pool the data
    grouped(exp).ring.all = counts;
    grouped(exp).ring.percent = ring_per; 
    grouped(exp).ring.avg = mean(ring_per,2,'omitnan');
end

% Pull the data together: 
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
    all_ring = [];
    
    % Update the averages for the classic temperature bins 
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        end
        raw_c(t,:) = mean(grouped(exp).ring.percent(c_frames,:),1,'omitnan');
        raw_h(t,:) = mean(grouped(exp).ring.percent(h_frames,:),1,'omitnan');
    end

    % find the avg and err and save to group structure
    grouped(exp).ring.increasing.raw = raw_h;
    grouped(exp).ring.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).ring.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).ring.decreasing.raw = raw_c;
    grouped(exp).ring.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).ring.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).ring.temps = temps;
end


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


%% ANALYSIS: Event aligned signals  -- food occupancy
clearvars('-except',initial_vars{:})
data_types = {'occ', 'dist', 'speed', 'ring'};

% account for any temperature hold groups
fictive_proto = false;
if any([data(:).hold_exp])
     result = questdlg(['There are temperature hold experiments, '...
         'do you want to assign them a fictive temperature protocol for '...
         'temperature event aligned comparisons?']);
    if strcmp(result, 'Yes') % offer a selection from the other trials within the grouping
        protocol_options = unique({data(:).temp_protocol});
        idx  = listdlg('PromptString','Select fictive temp protocol','ListString',protocol_options, 'ListSize',[175, 150]);
        fictive_protocol = protocol_options{idx};
        fictive_proto = true;
    end
end

for i = 1:num.exp
    if fictive_proto && data(i).hold_exp
        grouped(i).aligned.temp_protocol = fictive_protocol;
    else
        grouped(i).aligned.temp_protocol = data(i).temp_protocol;
    end
    tp = getTempTurnPoints(grouped(i).aligned.temp_protocol);

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
        [speed, occ,dist,ring,temperature] = deal([]);
        duration = min(tp.(tpBin)(:,2)-tp.(tpBin)(:,1)); %get shortest ramp period
        time = 1:duration;
        time = time./(tp.fps*60);
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
            ring(:,:,rr) = grouped(i).ring.percent(ROI,:);
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






%% ANALYSIS: Food quadrant occupation
clearvars('-except',initial_vars{:})

for exp = 1:num.exp
    quad_occ = [];
    for trial = 1:num.trial(exp)
        center = data(exp).data(trial).data.centre;
        x_loc = data(exp).data(trial).data.x_loc;
        y_loc = data(exp).data(trial).data.y_loc;
        r = data(exp).data(trial).data.r;
        foodWell = data(exp).T.foodLoc(trial);

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
        adjusted_wx = data(exp).data(trial).data.wellcenters(1,foodWell) - center(1);
        adjusted_wy = data(exp).data(trial).data.wellcenters(2,foodWell) - center(2);
        
        idx_loc = false(1,4);
        % Find the food quadrant (find location with the food well coordinates included)
        idx_loc(1) = (adjusted_wy > adjusted_wx) & (adjusted_wy <= -adjusted_wx);  % top
        idx_loc(2) = (adjusted_wy <= adjusted_wx) & (adjusted_wy <= -adjusted_wx); % right
        idx_loc(3) = (adjusted_wy <= adjusted_wx) & (adjusted_wy > -adjusted_wx);  % bottom
        idx_loc(4) = (adjusted_wy > adjusted_wx) & (adjusted_wy > -adjusted_wx);   % left
        quad_loc = find(idx_loc);

        fly_loc = ~isnan(x_loc); %gives logical for all fly positions in the position matrix
        foodQuad = Q(quad_loc).Mask & fly_loc; % flies in food quad
        nflies = data(exp).T.NumFlies(trial);
        y = (sum(foodQuad,2)./nflies).*100;
        quad_occ = autoCat(quad_occ, y,false);

        %  A = quad_occ;
        %  B = data(exp).data(trial).data.occupancy.dist2wells(:,foodWell);
        %  C = data(exp).data(trial).data.occupancy.occ(:,foodWell)*100;
        % 
        % figure; hold on
        % plot(A,'color', Color('red'));
        % plot(C,'color', 'k')
        % yyaxis right 
        % plot(B, 'color', 'y')
        % R = corrcoef(A,B);
        % disp(['Distance: ' num2str(R(2))])
        %  R = corrcoef(A,C);
        % disp(['ROI: ' num2str(R(2))])
    
    end

    % save the data to a larger structure
    grouped(exp).quadrant.all = quad_occ;
    grouped(exp).quadrant.avg = mean(quad_occ,2,'omitnan');
    grouped(exp).quadrant.std = std(quad_occ,0,2,'omitnan');
end
          
% Pool the data for heating and cooling together: 
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
    all_quad = [];
    
    % Update the averages for the classic temperature bins 
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        end
        raw_c(t,:) = mean(grouped(exp).quadrant.all(c_frames,:),1,'omitnan');
        raw_h(t,:) = mean(grouped(exp).quadrant.all(h_frames,:),1,'omitnan');
    end

    % find the avg and err and save to group structure
    grouped(exp).quadrant.increasing.raw = raw_h;
    grouped(exp).quadrant.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).quadrant.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).quadrant.decreasing.raw = raw_c;
    grouped(exp).quadrant.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).quadrant.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).quadrant.temps = temps;
end