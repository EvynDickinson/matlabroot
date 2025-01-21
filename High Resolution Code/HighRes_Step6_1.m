

%% Load grouped data to compare across flies and ramps

clear; close all; clc
baseFolder = [getDataPath(6,0),'Trial Data/'];

% Find files that can be run
[excelfile, Excel, xlFile] = load_HighResExperiments;

% Find trials that have nothing in their Proofed Column & have been tracked: 
loc = cellfun(@isnan,excelfile(2:end,Excel.groupready));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.expID, Excel.ramp]);
FileNames = format_eligible_files(eligible_files);

fileIdx = listdlg('ListString', FileNames,'ListSize',[350,450],'promptstring', 'Select trial to process');
if isempty(fileIdx)
    disp('No trials selected')
    return
end
% pull the list of dates and arenas to be loaded
% TODO: need some way to pair these by group and ramp 
% dateStr = eligible_files(fileIdx,1);
% trialStr = eligible_files(fileIdx,2); 

% for now, load all the trials as 'independent' trials

field_list = {'T','data','f', 'fX', 'fY', 'fps', 'm', 'mX', 'mY', 'nvids', 'parameters', 'pix2mm', 'position', 'tRate', 'time', 'well'};
fly = struct;
for i = 1:length(fileIdx)
    dateStr = eligible_files{i,1};
    trialStr = eligible_files{i,2}; 
    dummy = load([baseFolder dateStr '_' trialStr '/post-5.1 data.mat']);
    fly(i).name = [dateStr '_' trialStr];
    for ii = 1:length(field_list)
        fly(i).(field_list{ii}) = dummy.(field_list{ii});
    end
end

F = 2;
M = 1;
body = dummy.body;

num.trials = length(fly);

clear i ii dummy rownums ans fileIdx field_list FileNames loc trialStr dateStr
initial_var = who; initial_var{end+1} = 'initial_var';
clearvars('-except',initial_var{:})




























