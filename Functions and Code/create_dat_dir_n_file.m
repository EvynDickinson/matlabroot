function [blank_data_file, param] = create_dat_dir_n_file(param)
% 
% [blank_data_file, param] = create_dat_dir_n_file(param)
% 
% Create or check for the existance of all the data
% folders|directories, opens a blank file for data
% to fill during the experiment
% 
% ES Dickinson, University of Washington, 2019


% Get information for today's experiment
FormatOut = 'mmddyyyy';
Today = datestr(datetime,FormatOut);
param.matlab_data_file = [Today '_fly' param.fly_num];
param.folder_date = Num2Month(date);    % today's date


% create folder for today's experiment
param.save_location_matlab = [param.directory param.folder_date];
if ~exist(param.save_location_matlab, 'dir') 
    mkdir(param.save_location_matlab);
    fprintf('\n Created directory: \n')
    disp(param.save_location_matlab)
end

% create folder for fly number:
param.fly_dir = [param.save_location_matlab '\Fly ' param.fly_num '\'];
if ~exist(param.fly_dir, 'dir')
    mkdir(param.fly_dir);
    fprintf(['\n Added folder for fly ' param.fly_num '\n'])
end

% create folder for fly data 'fictrac' etc.
param.data_dir = [param.fly_dir, 'FicTrac Data\'];
if ~exist(param.data_dir, 'dir')
    mkdir(param.data_dir);
    fprintf(['\n Added data folder\n'])
end


% Create folder for videos for each fly 
param.Basler_folder_name = [param.fly_dir '\Uncompressed Video\'];
if ~exist(param.Basler_folder_name, 'dir')
    mkdir(param.Basler_folder_name);
    fprintf('\n Added uncompressed video folder \n')
end

% Specify the data files
filename = 'test_'; now = datetime('now','TimeZone','local');
formatOut = 'yyyymmddHHMMSS'; fileroot = 'C:\matlabroot\data\';
%opens this data file to write into it
blank_data_file = fopen([fileroot filename datestr(now,formatOut)],'w'); 


end




