
figure;
exp = 1;
trial = 1;
y = ones(size(sROImask(exp).all(trial).m));
y = y.*[1:size(y,2)];
y(~sROImask(exp).all(trial).m) = nan;

plot(y)

%% Load data
clear; clc;
baseFolder = getDataPath(2,0);

% Load excel file
[excelfile, Excel, xlFile] = load_SurvivalQuadCounts;

% Find the experiments that have been counted
loc = cellfun(@isnan,excelfile(2:end,Excel.counted));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.expID, Excel.arena, Excel.counted]);
FileNames = format_eligible_files(eligible_files);

% Select files 
fileIdx = listdlg('ListString', FileNames,'ListSize',[350,450],'promptstring', 'Select trial to process');
if isempty(fileIdx)
    disp('No trials selected')
    return
end

% Establish date, trial, and base directories
dateDir = eligible_files(fileIdx(1),1);
trialDir = eligible_files(fileIdx(1),2); 
baseDir = [baseFolder, dateDir{:} '\', trialDir{:} '\'];

data = [];
nfiles = length(fileIdx);

for i = 1:nfiles
    data = excelfile(i)