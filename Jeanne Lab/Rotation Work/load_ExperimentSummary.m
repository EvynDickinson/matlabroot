
function [excelfile, Excel, xlFile] = load_ExperimentSummary
% [excelfile, Excel, xlFile] = load_ExperimentSummary
% Loads the excel file data into the structure 'raw' and finds the col headers
% Outputs:
% 'excelfile' [excel file: fly summary 2]
% 'Excel' [structure with locations of column headers]
% 'xlFile' ''E:\Evyn\Experiment Summary.xlsx''
%
%
% ES Dickinson, Yale 2020


% [baseName, folder] = uigetfile('*xlsx', 'Select the Fly Summary 2 file');
% xlFile = fullfile(folder, baseName);
xlFile = 'E:\Evyn\Experiment Summary.xlsx';
%load excel sheet data
[~,~,excelfile] = xlsread(xlFile,'Sheet2');

% Excel column for the headers:
Excel.headers = excelfile(1,:); %xltitles
Excel.date = find(strcmpi('Date',Excel.headers) == 1);
Excel.flynum = find(strcmpi('Fly Num',Excel.headers) == 1);

Excel.genetic_cross = find(strcmpi('Genotype',Excel.headers) == 1);
Excel.sex = find(strcmpi('Sex',Excel.headers) == 1);
Excel.age = find(strcmpi('Age',Excel.headers) == 1);
Excel.experiment = find(strcmpi('Experiement',Excel.headers) == 1);

Excel.concentration = find(strcmpi('Concentration',Excel.headers) == 1);
Excel.odor = find(strcmpi('odor',Excel.headers) == 1);
Excel.odor_dur = find(strcmpi('odor dur',Excel.headers) == 1);
Excel.shock_dur = find(strcmpi('shock dur',Excel.headers) == 1);
Excel.stim_delay = find(strcmpi('shock delay',Excel.headers) == 1);
Excel.stim_delay = find(strcmpi('odor delay',Excel.headers) == 1);
Excel.trail_duration = find(strcmpi('trail duration',Excel.headers) == 1);
Excel.folder_name = find(strcmpi('folder name',Excel.headers) == 1);
Excel.tag = find(strcmpi('tag',Excel.headers) == 1);
Excel.session = find(strcmpi('session',Excel.headers) == 1);
Excel.time = find(strcmpi('time',Excel.headers) == 1);
end

