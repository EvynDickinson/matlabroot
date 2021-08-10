
% clear background
clc; close all; clear all
origState = warning;
warning('off')   
fps = 300;
buffer = 10; %pixel buffer around the edge of the video
box_size = 25; %length of light indicator box's sides
% root_directory = 'E:\Basler Trig\Behavior Phenotypes/';
root_directory = 'E:\Basler Trig\INTERESTING EFFECTS/Backwards Walking/';
% Load the Behaviors Notebook from Excel  
[~, ~, vids] = load_flybehaviors;
num.vids = size(vids,2);

%% ---
% root_directory = 'E:\Basler Trig\INTERESTING EFFECTS/Backwards Walking/';
% num.vids = 1;
% vids(1).start_frame = 1;
% vids(1).duration  = 599;
% vids(1).end_frame       = start_frame+frame_duration;
% vids(1).replay_fps      = 30;
% 
% 
% 
% pathroot = 'E:\Basler Trig\INTERESTING EFFECTS/';
% selpath  = uigetdir(pathroot, 'Select the video folder');
% % Get the video names|frames for each camera:
% FilePath = dir([selpath, '\*C*']); 
% 
%        
%  
% 
% 
% ProcessedImages.light_length = 1;
% 
% flyVideo = VideoReader(video_name);
% 
% 
% num.cam = length(FilePath);
% for ii = 1:num.cam
%     % Read video information and get name:
%     cam.(['name_' Alphabet(ii)]) =  [FilePath(ii).folder, '\', FilePath(ii).name];
%     cam.(Alphabet(ii)) = VideoReader(cam.(['name_' Alphabet(ii)]));
%     % Load frames as a single B&W image matrix:
%     video_name = cam.(['name_' Alphabet(ii)]); 
%     [~, ProcessedImages.(Alphabet(ii))] = ...
%         convert_vid_to_pixel_values(video_name, 1, 0);
%     fprintf(['Camera ' num2str(ii) ' loaded \n'])
% end
% 
% % determine the light length:
% switch str2double(COND)
%     case {1, 8, 15, 22}
%         light_length = 0.0;
%     case {2, 9, 16, 23}
%         light_length = 0.03;
%     case {3, 10, 17, 24}
%         light_length = 0.06;
%     case {4, 11, 18, 25}
%         light_length = 0.09;
%     case {5, 12, 19, 26}
%         light_length = 0.18;
%     case {6, 13, 20, 27}
%         light_length = 0.32;
%     case {7, 14, 21, 28}
%         light_length = 0.72;
% end
% 
% ProcessedImages.light_length = light_length;
% %---

%% Iterate Through the Videos:
for ivid = 1:num.vids
    % convert this to work in the loop format
    % Video Parameters: %These Must Be Selected Each Time%
    start_frame     = vids(ivid).start_frame;
    frame_duration  = vids(ivid).duration;
    end_frame       = start_frame+frame_duration;
    frame_rate      = vids(ivid).replay_fps;

    % Select the regions/load the videos
    [ProcessedImages, FilePath] = select_video_file(vids, ivid); 
    light_length = ProcessedImages.light_length;
    
    % Resize the images to the same height
    [ProcessedImages, height, width] = scale_and_size(ProcessedImages);
    
    temp.row1 = width.F+width.E+width.D; %top row max pixel size
    temp.row2 = width.A+width.B+width.C; %bottom row max pixel size
    width.VID = max([temp.row1,temp.row2]) + (4*buffer); %total video width
    height.VID = (height.A*2) + (3*buffer); %total video height
    clear temp
    
    % BLANK FIGURE create a black video frame:
    image_size = [height.VID, width.VID];

    % Create images of text to superimpose on videos
    %cross:
    cross = vids(ivid).genetic_cross; cross = strrep(cross,'_','-');
    ProcessedImages.cross_image = generate_image(cross, 0.8, -90);
        height.cross = ProcessedImages.cross_image.height; 
        width.cross = ProcessedImages.cross_image.width;
    
    %camera labels:
    labels = {'R', 'L', 'T'};
    for ii = 1:3
        ProcessedImages.([labels{ii} '_image']) = generate_image(labels{ii}, 1.1);
        height.(labels{ii}) = ProcessedImages.([labels{ii} '_image']).height;
        width.((labels{ii})) = ProcessedImages.([labels{ii} '_image']).width;
    end
    
    % Generate pixel data for the light indicator box
    ProcessedImages = LIB_generation(ProcessedImages, box_size);
    width.LIB = box_size;
    height.LIB = box_size;
    
    % Create the pixel ranges for each camera | LIB | label
    [H_zone, W_zone] = GenerateCamPositions(buffer, height, width);
    
    % Set blank image as background
    ProcessedImages.blank_frame = uint8(zeros(image_size)); 
    
%     % Preview a frame:
%     frame_num = 0.5*300 + 5; %light should be on at this point
%     buffer_frame = generate_single_frame(frame_num, ProcessedImages, H_zone, W_zone);
%     fig = figure;
%     set(fig, 'Name', 'Preview Only')
%     imshow(buffer_frame);
%     disp('Press a key to continue!') 
%     pause;
%     close(fig)

    % Saving Info:
    analysis_directory = [root_directory, cross '/'];
    if ~exist(analysis_directory, 'dir')
        mkdir(analysis_directory);
    end
    addpath(analysis_directory)
    temp.a = length(ProcessedImages.flydate);
    if temp.a == 5
        ProcessedImages.flydate = ['0' ProcessedImages.flydate];
    elseif temp.a == 4 %both a single for month and day
        ProcessedImages.flydate = ['0' ProcessedImages.flydate(1) '0' ProcessedImages.flydate(2:4)];
    end
    clear temp 
    fly_ID = [ProcessedImages.flydate(1:4) '20' ProcessedImages.flydate(5:6)...
              '_fly' vids(ivid).flynum ' R' num2str(vids(ivid).rep)...
              'C' num2str(vids(ivid).cond)];
    video_name = [fly_ID ' ' num2str(frame_rate) '-fps'];
    
    % Make the Video:
    vid_name = [analysis_directory video_name ' uncompressed.avi'];
    fig = figure; set(fig, 'pos',[300,190,1321,776], 'color','k');
    currAxes.Visible = 'off';
    hold on
    v = VideoWriter(vid_name, 'Uncompressed AVI');
    v.FrameRate = frame_rate;
    open(v);
    for datapoint = [1, start_frame+2:end_frame+1]
        set(fig, 'color','k');
        buffer_frame = generate_single_frame(datapoint, ProcessedImages, H_zone, W_zone);
        if datapoint == 1
          buffer_frame = uint8(zeros(image_size));  
        end
        %display the iamge
        imshow(buffer_frame, 'border', 'tight')
        currAxes.Visible = 'off';
        % Collect image info for movie file
        f = getframe(fig);  
        writeVideo(v, f)
        clf('reset')
    end
    close(v)
    close all 
    fprintf('\n Video Saved! \n')  

    % Compress the video
    new_vid = [analysis_directory video_name];
    reader = VideoReader(vid_name);        
    writer = VideoWriter(new_vid, 'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
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
    
    % Write back into Excel the Video Name and Change Symbol from 'y' to
    % 'D' for done
    excelinfo.vid_name = video_name;
    excelinfo.row_num = vids(ivid).excel_idx;
    load_flybehaviors(excelinfo);
    
    fprintf(['\n Finished Video ' num2str(ivid) '/' num2str(num.vids) '\n \n'])

end
fprintf('\n DONE! \n')
warning(origState)
clear all







