
% clear background
clc; close all; clear all
tic
origState = warning; 
warning('off') 
fps = 300;
buffer = 10; %pixel buffer around the edge of the video
box_size = 25; %length of light indicator box's sides
root_directory = 'E:\Basler Trig\Behavior Phenotypes/';


% Load the Behaviors Notebook from Excel  
[~, ~, vids] = load_flybehavior_RA;
num.vids = size(vids,2); 

% Iterate Through the Videos:
for ivid = 1:num.vids
    % Video Parameters: %These Must Be Selected Each Time%
    start_frame     = vids(ivid).start_frame; 
    frame_duration  = vids(ivid).duration;
    end_frame       = start_frame+frame_duration;
    frame_rate      = vids(ivid).replay_fps;
    for icond = [1:7, 15:21]
        vids(ivid).cond = icond;
        % Select the regions/load the videos
        [ProcessedImages, FilePath] = select_video_file_RA(vids, ivid); 
        light_length = ProcessedImages.light_length;

        % Resize the images to the same height
        [ProcessedImages, height, width] = scale_and_size_RA(ProcessedImages);

        width.VID = (width.A*3) + (4*buffer); %total video width
        height.VID = (height.A*2) + (5*buffer); %total video height

        % BLANK FIGURE create a black video frame:
        image_size = [height.VID, width.VID];

        % Create images of text to superimpose on video s
        %cross:
        cross = vids(ivid).genetic_cross; cross = strrep(cross,'_','-');
        name.a =  strsplit(FilePath(ivid).name);
        name.b = strsplit(name.a{4}, '-');
        desc = [cross ' -- ' name.b{1} ' ' name.b{3} ' sec'];

        ProcessedImages.cross_image = generate_image(desc, 0.8);
            height.cross = ProcessedImages.cross_image.height; 
            width.cross = ProcessedImages.cross_image.width;

        %camera labels:
        labels = {'ccw-2', 'ccw-3', 'ccw-1', 'cw-3', 'cw-1', 'cw-2'};
        for ii = 1:length(labels)
            ProcessedImages.([Alphabet(ii) '_label']) = generate_image(labels{ii}, 1.1);
            height.([Alphabet(ii) '_label']) = ProcessedImages.([Alphabet(ii) '_label']).height;
            width.([Alphabet(ii) '_label']) = ProcessedImages.([Alphabet(ii) '_label']).width;
        end

        % Generate pixel data for the light indicator box
        ProcessedImages = LIB_generation(ProcessedImages, box_size);
        width.LIB = box_size;
        height.LIB = box_size;

        % Create the pixel ranges for each camera | LIB | label
        [H_zone, W_zone] = GenerateCamPositions_RA(buffer, height, width);

        % Set blank image as background
        ProcessedImages.blank_frame = uint8(zeros(image_size)); 

        % Preview a frame:
    %     frame_num = 0.5*300 + 5; %light should be on at this point
    %     buffer_frame = generate_single_frame_RA(frame_num, ProcessedImages, H_zone, W_zone);
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

        fly_ID = [name.a{1} ' ' name.b{1} '-' name.b{3} ' ' name.a{3}]; clear temp      
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
            buffer_frame = generate_single_frame_RA(datapoint, ProcessedImages, H_zone, W_zone);
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


    fprintf(['\n Finished cond ' num2str(icond) '\n \n'])
        
    end
    fprintf(['\n Finished Video ' num2str(ivid) '/' num2str(num.vids) '\n \n'])
        % Write back into Excel the Video Name and Change Symbol from 'y' to
        % 'D' for done
        video_name = [fly_ID(1:15) ' Cam-' vids(ivid).cam ' ' num2str(frame_rate) '-fps'];
        excelinfo.vid_name = video_name;
        excelinfo.row_num = vids(ivid).excel_idx;
        load_flybehavior_RA(excelinfo);
end
fprintf('\n DONE! \n')
warning(origState)

toc





