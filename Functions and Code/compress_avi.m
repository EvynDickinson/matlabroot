

function compress_avi(date_folder)
%COMPRESS_AVI Compresses avi files 
% 
%
% 
%   date_folder = name of date folder as string, e.g. '20181003'
% 
%   CJ Dallmann, University of Washington, 10/2018
%   Edited and adjusted by ES Dickinson, UW, 12/2018

% date_folder = '10.31.18';

% Set input directory
dir_input =  ['E:\FicTrac Raw Data\', date_folder, '\Video\'];

% Create date folder & set output directory
dir_root = 'D:\Evyn Data Files\';

% Get list of fly folders
list_flies = dir([dir_input,'*Fly*']);

% Get total number of videos in all sessions
nVideos = 0;
for iSession = 1:length(list_flies)
    list_videos = dir([dir_input,list_flies(iSession).name,'\*.avi']);
    nVideos = nVideos+length(list_videos);
end

% Create waitbar
h = waitbar(0,['Compressing videos from ',date_folder]);

% Iterate through sessions
nVideo = 1;
for iSession = 1:length(list_flies)  %per each fly session
    % Get list of videos in current session
    list_videos = dir([dir_input,list_flies(iSession).name,'\*.avi']);

    % Designate session folder in output directory
    dir_output = ['D:\Evyn Data Files\' date_folder '\' list_flies(iSession).name '\Raw Video\']; 
    
    % Iterate through videos
    for iVideo = 1:length(list_videos)
        % Create objects to read and write the video, and open the AVI file
        % for writing
        disp(['Compressing ',list_videos(iVideo).name,' ...'])
        reader = VideoReader([dir_input,list_flies(iSession).name,'\',list_videos(iVideo).name]);        
        writer = VideoWriter([dir_output,list_videos(iVideo).name],'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
        writer.FrameRate = reader.FrameRate;
        writer.Quality = 75;
        open(writer);
        
        % Read and write each frame
        while hasFrame(reader)
            img = readFrame(reader);
            writeVideo(writer,img);
        end
        writeVideo(writer,img);
        close(writer);           
        disp('Done!')
        
        % Update waitbar
        waitbar(nVideo/nVideos,h)
        nVideo = nVideo+1;
    end
end
close(h)
clc

