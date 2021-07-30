
function [excelfile, Excel, xlFile] = load_flysummary
% [raw, Excel, xlFile] = load_flysummary
% Loads the excel file data into the structure 'raw' and finds the col headers
% Outputs:
% 'excelfile' [excel file: fly summary 2]
% 'Excel' [structure with locations of column headers]
% 'xlFile' 'G:\My Drive\Data\FicTrac Raw Data\Fly Summary 2.xlsx'
%
%
% ES Dickinson, University of Washington, Dec 2018


% [baseName, folder] = uigetfile('*xlsx', 'Select the Fly Summary 2 file');
% xlFile = fullfile(folder, baseName);

xlFile = 'G:\My Drive\Tuthill Lab Shared\Evyn\Data\FicTrac Raw Data\Fly Summary 2.xlsx';
% 'G:\My Drive\Data\FicTrac Raw Data\Fly Summary 2.xlsx';      
% xlFile = '/Volumes/Evyn SSD/Evyn UW work/Fly Summary 2.xlsx';
%load excel sheet data
[~,~,excelfile] = xlsread(xlFile);

% Excel column for the headers:
Excel.headers = excelfile(1,:); %xltitles
Excel.date = find(strcmpi('Date',Excel.headers) == 1);
Excel.flynum = find(strcmpi('Fly #',Excel.headers) == 1);
Excel.structurenum = find(strcmpi('Fly(#)',Excel.headers) == 1);
Excel.new_struct_name = find(strcmpi('Structure Name',Excel.headers) == 1);
Excel.genetic_cross = find(strcmpi('Fly Line',Excel.headers) == 1);
Excel.structure = find(strcmpi('Structure Name', Excel.headers)==1);
Excel.distance = find(strcmpi('Total Distance', Excel.headers)==1);
Excel.frames = find(strcmpi('Frames Aligned', Excel.headers)==1);

end

