
function [folder, subFolders] = cropCheck(folder, arena)
% [folder, subFolders] = cropCheck(folder, arena)
% Pull image of cropped video to compare to full size image
% ES Dickinson
% Yale University 2021

if nargin==0
    % Select the complete experiments to process
    [~, folder] = getCloudPath(1); 
    list_dirs = dir(folder); 
    dirFlags = [list_dirs.isdir]; % logical vector for directories.
    
    % Extract only those that are directories.
    subFolders = list_dirs(dirFlags);
    subFolders = {subFolders(3:end).name};
    % select the arena
    arena = subFolders{listdlg('ListString', subFolders, 'SelectionMode', 'Single')};
end
% select the video
vidList = dir([folder '\' arena '\*.avi']);
vidList = {vidList(:).name};
video = vidList{listdlg('ListString', vidList, 'SelectionMode', 'Single',...
                        'PromptString', 'Select vid for image comparison')};

% Pull the full video and the cropped vid
OGvid = [folder '\' video];
movieInfo = VideoReader(OGvid); %read in video
OGimg = read(movieInfo,1);
CROPvid = [folder '\' arena '\' video];
movieInfo = VideoReader(CROPvid); %read in video
CROPimg = read(movieInfo,1);

% Visually compare the full and cropped video
fig = figure; set(fig, 'pos', [233 239 1454 733], 'color', 'k')
imshowpair(OGimg, CROPimg, 'montage')
title(arena, 'Color', 'w')
annotation('textbox',[.16 .47 .5 .3],'String',' B','FitBoxToText','on', 'color', 'w');
annotation('textbox',[.16 .035 .5 .3],'String',' A','FitBoxToText','on', 'color', 'w');
annotation('textbox',[.38 .47 .5 .3],'String',' D','FitBoxToText','on', 'color', 'w');
annotation('textbox',[.38 .035 .5 .3],'String',' C','FitBoxToText','on', 'color', 'w');
% crop label
annotation('textbox',[.59 .47 .5 .3],'String', [' ' arena],'FitBoxToText','on', 'color', 'w');

end