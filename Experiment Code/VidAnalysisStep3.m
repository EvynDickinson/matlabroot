


%% Load in multiple flies that are grouped together:
% baseFolder = getCloudPath;
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
baseFolder = getCloudPath;

% Select structure to load:
structure_names = unique(excelfile(3:end, Excel.structure));
ExpGroup = structure_names{listdlg('ListString', structure_names, 'ListSize', [300, 400])};

% Pull structure details:
loc = find(strcmpi(ExpGroup, excelfile(:,Excel.structure))); %rows with sel exp.
% Find the number of trials within the structure
[trial, rowloc] = deal([]);
for tt = loc
    a = cell2mat(excelfile(tt, Excel.structurenum));
    if ischar(a) %no value for the structure
    else %has a number in the structure
        trial = [trial, a]; %find total number of flies
        rowloc = [rowloc, tt]; %get the rownumber for each fly
    end
end
ntrial = max(trial); clear flynum bb a tt
% Error Check
if sum(ntrial) == 0
    warndlg('Check name of flies in Excel File, no matching files found')
    return
end

% Load data from each trial
for trial = 1:ntrial
    % print the experiments as they are loaded
    trialName = excelfile{rowloc(trial), Excel.expID};
    trialDate = excelfile{rowloc(trial), Excel.date};
    fprintf('\n Trial: \n'); disp([trialDate,'  ' trialName])

    % build the path for the trial data
    dirc = [baseFolder, trialDate, '\analysis\' trialName];
    nvids = length(dir([dirc, '*.mat']));
    
    % load data
    load()

        
end