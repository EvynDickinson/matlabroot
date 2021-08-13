clear
%% 
%raw data folder:
date_today = strrep(datestr(datetime,'mm-dd-yyyy'),'-','.');
start_dir = ['C:\Users\jeannelab\Documents\Evyn\DATA\'];

%get base folder pathway
switch getenv('COMPUTERNAME')
    case 'DENALI'
        baseFolder = 'E:\My Drive\Jeanne Lab\DATA\';
    case 'TOGIAK'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';  
end

%select folder date   
%TODO update this to iterate through multiple folders if need be
list_dirs = dir(start_dir);
for i = 3:length(list_dirs)
    folderNames{i-2} = list_dirs(i).name;
end
folderNames = ['Today', folderNames];

indx = listdlg('ListString', folderNames, 'SelectionMode', 'Single');
if strcmpi(folderNames{indx}, 'Today')==true
    dir_sel = date_today;
else
    dir_sel = folderNames{indx+1};
end
folder = fullfile(start_dir, dir_sel);

% Move folders to google drive:
copyfile(folder, [baseFolder dir_sel])

fprintf('done')    