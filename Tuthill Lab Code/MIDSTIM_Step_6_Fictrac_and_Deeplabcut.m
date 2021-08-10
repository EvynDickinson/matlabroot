

close all; clear all
 
rep = 3;
cond = 4;

% Select fly data:
folder_date = '11.26.18'; % date of folder with videos/data acq.
Folder_date = folder_date; % date of folder with videos/data acq.
pathway = ['D:\Evyn Data Files\' folder_date]; % directory
list_flies = dir(pathway);
for ii = 3:length(list_flies)
    fly_num{ii-2} = list_flies(ii).name; %names of the flies for the day
end
VideoDir = 'G:\My Drive\Data\FicTrac Raw Data\';
%Select the fly subsets to analyze
indx = listdlg('ListString', fly_num, 'SelectionMode', 'Multiple');


for kk = 1:length(indx)
    ifly = indx(kk);    
    % Load the data:   (maybe make this a function...) 
    filepath = [pathway '\' fly_num{ifly} '\Analysis\'];
    folder_date(regexp(folder_date,'[.]'))=[];
    fly_ID = [folder_date(1:end-2), '2018_fly' fly_num{ifly}(end-2:end)];
    %load the combined fly data
    load([filepath fly_ID])
    fprintf(['\n Loaded: ' fly_ID '\n'])
    
    % add address/locations to the fly parameters
    fly.parameters.date = Folder_date;
    fly.parameters.vid_dir = VideoDir;
    loc = [VideoDir Folder_date '\Fly ' fly.parameters.fly_num '\'];
    fly.parameters.vid_folder = loc;
    fly.parameters.vid_3d = [loc 'videos-3d\'];
    fly.parameters.pose_2d = [loc, 'pose-2d-filtered\'];
    
    % % THROW DATA TO THE FUNCTION THAT WILL GENERATE THE PLOT % %
    num = NUM(fly.parameters);
   
%     combined_vid_n_plots(fly, cond, rep) 
end

%%  Compress the videos

% Set input directory
dir_input = uigetdir;


% Create date folder & set output directory
dir_root = [dir_input, '\Condition Videos\'];
if ~isfolder(dir_root)
    mkdir(dir_input, 'Condition Videos')
end
% 
% % Get list of fly folders
% list_flies = dir([dir_input,['*' fly.parameters.matlab_data_file ' R*']]);

% % Get total number of videos in all sessions
% nVideos = 0;
% for iSession = 1:length(list_flies)
%     list_videos = dir([dir_input,list_flies(iSession).name]);
%     nVideos = nVideos+length(list_videos);
% end

% % Create waitbar
% h = waitbar(0,['Compressing videos from ', fly.parameters.matlab_data_file]);
% Get list of videos in current session
list_videos = dir([dir_input, '\*.avi']);
% Designate session folder in output directory
dir_output = dir_root; 

% Iterate through videos
for iVideo = 1:length(list_videos)
    % Create objects to read and write the video, and open the AVI file
    % for writing
    disp(['Compressing ',list_videos(iVideo).name,' ...'])
    reader = VideoReader([dir_input '\' list_videos(iVideo).name]);        
    writer = VideoWriter([dir_output '\' list_videos(iVideo).name(1:end-4)],'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
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

%     % Update waitbar
%     waitbar(nVideo/nVideos,h)
%     nVideo = nVideo+1;
end


clc









