
function data_org(Date)
%   DATA_ORG copies data from a single session, Date, into 
%   file locations on the D: drive and G: google drive
%   It also compresses and copies video files from the session
% 
%   date_folder = name of date folder as string, e.g. '20181003'
% 
%   ES Dickinson, University of Washington, 12/2018

% Directories to use
d_dir = ['D:\Evyn Data Files\', Date]; %computer storage, d drive
e_dir = ['E:\FicTrac Raw Data\', Date]; %initial data location, e drive
g_dir = ['G:\My Drive\Data\FicTrac Raw Data\', Date]; %google drive

% Inventory number of flies run per day
list_files = dir([e_dir '\Video']);
for ii = 3:length(list_files)
    fly_num{ii-2} = list_files(ii).name; %names of the flies for the day
end
clear list_files

% Create new folder system in the two saving directories
mkdir(d_dir) %make directory
mkdir(g_dir) %make directory

% Copy files over to the new 
for ii = 1:length(fly_num)
%-------------------------------------------------------------------------%
%Analysis directory
    analysis_dir1 = [d_dir '\' fly_num{ii} '\' Analysis];
    mkdir(analysis_dir1)
    analysis_dir2 = [g_dir '\' fly_num{ii} '\' Analysis];
    mkdir(analysis_dir2)
    %search for the fly moniker
    a = fly_num{ii}(5:7);
    search_list = ['fly' a];  
    %find the data files with specific fly moniker in the common directory
    list_files = dir([e_dir '\Analysis' ['\*' search_list '*']]);
    %move the data files to their respective d and g drive analysis folders
    for nn = 1:length(list_files)
        copyfile([list_files(nn).folder '\' list_files(nn).name], analysis_dir1)
        copyfile([list_files(nn).folder '\' list_files(nn).name], analysis_dir2)
    end
    clear list_files a search_list  
%-------------------------------------------------------------------------%    
%FicTrac Data directory
    data_dir1 = [d_dir '\' fly_num{ii} '\' FicTrac Data];
    mkdir(data_dir1)
    data_dir2 = [g_dir '\' fly_num{ii} '\' FicTrac Data];
    mkdir(data_dir2)
    %search for the fly moniker
    a = fly_num{ii}(5:7);
    search_list = ['fly' a];
    %find the data files with specific fly moniker in the common directory
    list_files = dir([e_dir ['\*' search_list '*']]);
    %move the data files to their respective d and g drive analysis folders
    for nn = 1:length(list_files)
        copyfile([list_files(nn).folder '\' list_files(nn).name], data_dir1)
        copyfile([list_files(nn).folder '\' list_files(nn).name], data_dir1)
    end
    clear list_files a search_list 

%-------------------------------------------------------------------------%
%Raw Video directory -- make the folders but don't put anything into them yet
    vid_dir = [d_dir '\' fly_num{ii} '\' Raw Video];
    mkdir(vid_dir)
end

% compress and copy the videos into the D-drive
compress_avi(Date)

% copy the video files from the data drive into google drive
for ii = 1:length(fly_num)
    %Raw Video directories
    vidin_dir = [d_dir '\' fly_num{ii} '\' Raw Video];
    vidout_dir = [g_dir '\' fly_num{ii}];
    %copy the videos over to google drive
    copyfile(vidin_dir, vidout_dir)
end
clc
fprintf(['\n Successfully compressed, transferred and copied \n the data from ' Date ' to: \n'])
disp(d_drive)
disp(g_drive)
fprintf('\n Erase data from the temporary data collection drive (E:)\n')
end





