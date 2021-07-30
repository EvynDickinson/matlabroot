


function compress_vid_n_plots(fly)

% Set input directory
dir_input =  ['D:\Evyn Data Files\' Folder_date '\Fly ' fly.parameters.fly_num '\Analysis\'];   

% Create date folder & set output directory
dir_root = [dir_input, 'Condition Videos\'];
if ~isfolder(dir_root)
    mkdir(dir_input, 'Condition Videos')
end

% Get list of fly folders
list_flies = dir([dir_input,['*' fly.parameters.matlab_data_file ' R*']]);

% Get total number of videos in all sessions
nVideos = 0;
for iSession = 1:length(list_flies)
    list_videos = dir([dir_input,list_flies(iSession).name]);
    nVideos = nVideos+length(list_videos);
end

% Create waitbar
h = waitbar(0,['Compressing videos from ', fly.parameters.matlab_data_file]);
% Get list of videos in current session
list_videos = dir([dir_input,['*' fly.parameters.matlab_data_file ' R*']]);
% Designate session folder in output directory
dir_output = dir_root; 

% Iterate through videos
for iVideo = 1:length(list_videos)
    % Create objects to read and write the video, and open the AVI file
    % for writing
    disp(['Compressing ',list_videos(iVideo).name,' ...'])
    reader = VideoReader([dir_input,list_flies(iVideo).name]);        
    writer = VideoWriter([dir_output,list_videos(iVideo).name],'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
    writer.FrameRate = reader.FrameRate;
    writer.Quality = 60;
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

close(h)
clc

end

