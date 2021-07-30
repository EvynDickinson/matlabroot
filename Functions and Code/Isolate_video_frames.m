
%% Create a video to isolate a small segment of a larger video and create
% a quilt-style video with all 6 camera angles of the same frames
clear all
close all
clc

%% Video Parameters: %These Must Be Selected Each Time%
fps             = 300;

start_frame     = 100;
frame_duration  = 200;
end_frame       = start_frame+frame_duration;
frame_rate      = 30;

buffer          = 10; %pixel buffer around the edge of the video

% Select the regions/load the videos
[ProcessedImages, FilePath] = select_video_file;
light_length = ProcessedImages.light_length;

%% Scale the videos
[ProcessedImages, height, width] = scale_and_size(ProcessedImages);
 
%% POSITIONING of videos

temp.row1 = width.F+width.E+width.D; %top row max pixel size
temp.row2 = width.A+width.B+width.C; %bottom row max pixel size
Vid_width = max([temp.row1,temp.row2]) + (4*buffer); %total video width
Vid_height = (FinalHeight*2) + (3*buffer); %total video height

% BLANK FIGURE create a black video frame:
image_size = [Vid_height, Vid_width];


% TEXT LABELS:
% insert cross text:
txt = 'Name of the fly cross sits here';
imtext = 255*imcomplement(text2im(txt)); %convert to black and white rgb
crosstext = imrotate(uint8(imtext),-90);
crosstext = imresize(crosstext, 0.8);
cross_height = size(crosstext,1);
cross_width = size(crosstext,2);
H_zone.cross_txt = (buffer):(buffer + cross_height)-1;
W_zone.cross_txt = (Vid_width - buffer - cross_width):(Vid_width - buffer)-1;

% 'R' label:
txt = 'R';
imtext = 255*imcomplement(text2im(txt)); %convert to black and white rgb
Rtxt = uint8(imtext);
Rtxt = imresize(Rtxt, 1.1);
R_height = size(Rtxt,1);
R_width = size(Rtxt,2);

% 'L' label:
txt = 'L';
imtext = 255*imcomplement(text2im(txt)); %convert to black and white rgb
Ltxt = uint8(imtext);
Ltxt = imresize(Ltxt, 1.1);
L_height = size(Ltxt,1);
L_width = size(Ltxt,2);

% 'T' label:
txt = 'T';
imtext = 255*imcomplement(text2im(txt)); %convert to black and white rgb
Ttxt = uint8(imtext);
Ttxt = imresize(Ttxt, 1.1);
T_height = size(Ttxt,1);
T_width = size(Ttxt,2);
% figure;imshow(Ttxt);

% ATTACHING LIGHT INDICATOR TO THE TOP RIGHT CORNER OF CAMERA B
% % Light indicator box (LIB):
% on_val = 255; %white color
% box_size = 25;
% 
% LIB_Height = box_size;
% LIB_Width = box_size;
% LIB_on = on_val*ones(box_size, box_size);
% LIB_off = zeros(box_size, box_size);
% 
% % Make a dummy light on/off seqeunece:
% total_frames = numel(ProcessedImages.A)+3;
% light_on = (0.5*fps)+2;
% light_off = light_on + (light_length*fps)+2;
% for iframes = 1:total_frames
%     if iframes <= light_on || iframes >= light_off
%         temp.data = uint8(LIB_off);
%     elseif iframes > light_on && iframes < light_off
%         temp.data = uint8(LIB_on);
%     end
%     ProcessedImages.LIB(iframes).data = temp.data;
% end
% % Generate pixel data for the light indicator box
ProcessedImages = LIB_generation(ProcessedImages, box_size);

%% Positions: 
% Order: [EFD];[CAB]

% F CAMERA
H_zone.F = (buffer):(buffer + height.F)-1;
W_zone.F = (2*buffer + width.E):(2*buffer + width.E) + width.F-1;
H_zone.F_text = buffer:buffer + R_height-1;
W_zone.F_text = (2*buffer + width.E):(2*buffer + width.E)+R_width-1;

% E CAMERA
H_zone.E = (buffer):(buffer + height.E)-1;
W_zone.E = (buffer):(buffer + width.E)-1;
H_zone.E_text = buffer:buffer + R_height-1;
W_zone.E_text = buffer:buffer + R_width-1;
H_zone.E_LI = H_zone.E(1):H_zone.E(1)+LIB_Height-1;
W_zone.E_LI = W_zone.E(end)-LIB_Width:W_zone.E(end)-1;

% D CAMERA
H_zone.D = (buffer):(buffer + height.D)-1;
W_zone.D = (W_zone.F(end) + buffer):(W_zone.F(end) + buffer + width.D)-1;
H_zone.D_text = buffer:buffer + T_height-1;
W_zone.D_text = (W_zone.F(end) + buffer):(W_zone.F(end) + buffer) + T_width-1;
H_zone.D_LI = H_zone.D(1):H_zone.D(1)+LIB_Height-1;
W_zone.D_LI = W_zone.D(end)-LIB_Width:W_zone.D(end)-1;

% A CAMERA
H_zone.A = (Vid_height - buffer - height.A):(Vid_height - buffer)-1;
W_zone.A = (2*buffer + width.C):(2*buffer + width.C) + width.A-1;
H_zone.A_text = (Vid_height - buffer - height.A):(Vid_height - buffer - height.A)+L_height-1;
W_zone.A_text = (2*buffer + width.C):(2*buffer + width.C)+L_width-1;
H_zone.A_LI = H_zone.A(1):H_zone.A(1)+LIB_Height-1;
W_zone.A_LI = W_zone.A(end)-LIB_Width:W_zone.A(end)-1;

% C CAMERA
H_zone.C = (Vid_height - buffer - height.C):(Vid_height - buffer)-1;
W_zone.C = (buffer):(buffer + width.C)-1;
H_zone.C_text = (Vid_height - buffer - height.C):(Vid_height - buffer - height.C)+L_height-1;
W_zone.C_text = buffer:buffer+L_width-1;
H_zone.C_LI = H_zone.C(1):H_zone.C(1)+LIB_Height-1;
W_zone.C_LI = W_zone.C(end)-LIB_Width:W_zone.C(end)-1;

% B CAMERA
H_zone.B = (Vid_height - buffer - height.B):(Vid_height - buffer)-1;
W_zone.B = (W_zone.A(end) + buffer):(W_zone.A(end) + buffer + width.B)-1;
H_zone.B_text = (Vid_height - buffer - height.B):(Vid_height - buffer - height.B)+L_height-1;
W_zone.B_text = (W_zone.A(end) + buffer):(W_zone.A(end) + buffer)+L_width-1;
H_zone.B_LI = H_zone.B(1):H_zone.B(1)+LIB_Height-1;
W_zone.B_LI = W_zone.D(end)-LIB_Width:W_zone.D(end)-1;

for icam = 1:6
    H_data = H_zone.(Alphabet(icam));
    W_data = W_zone.(Alphabet(icam));
    H_zone.([Alphabet(icam) '_LI']) = H_data(1):H_data(1)+LIB_Height-1;
    W_zone.([Alphabet(icam) '_LI']) = W_data(end)-LIB_Width:W_data(end)-1;
end
H_zone.B_LI = H_zone.B(1):H_zone.B(1)+LIB_Height-1;
W_zone.B_LI = W_zone.D(end)-LIB_Width:W_zone.D(end)-1;


%% Preview a frame:

test_frame = 200;
buffer_frame = uint8(zeros(image_size)); 
% buffer_frame(LIB_Height, LIB_Width) = ProcessedImages.LIB(test_frame).data;
for icam = 1:6
    %image
    buffer_frame(H_zone.(Alphabet(icam)), W_zone.(Alphabet(icam))) = ...
        ProcessedImages.(Alphabet(icam))(test_frame).data;
    %text
    if icam <= 3 
        text_image = Ltxt;
    elseif icam >= 5
        text_image = Rtxt;
    else
        text_image = Ttxt;
    end
    buffer_frame(H_zone.([Alphabet(icam) '_text']), W_zone.([Alphabet(icam) '_text'])) = text_image;
    % LIB
    buffer_frame(H_zone.([Alphabet(icam) '_LI']), W_zone.([Alphabet(icam) '_LI'])) = ...
        ProcessedImages.LIB(test_frame).data;
end

 

buffer_frame(H_zone.cross_txt, W_zone.cross_txt) = crosstext;

fig = figure;
I = imshow(buffer_frame);


%% Saving Info:
root_directory = 'E:\Basler Trig\Behavior Phenotypes/';
cross = select_cross;
cross = strrep(cross,'_','-');

% fly cross text:
imtext = 255*imcomplement(text2im(cross)); %convert to black and white rgb
crosstext = imresize(imrotate(uint8(imtext),-90), 0.8);
cross_height = size(crosstext,1);
cross_width = size(crosstext,2);
H_zone.cross_txt = (buffer):(buffer + cross_height)-1;
W_zone.cross_txt = (Vid_width - buffer - cross_width):(Vid_width - buffer)-1;

% create folder for fly cross:
analysis_directory = [root_directory, cross '/'];
if ~exist(analysis_directory, 'dir')
    mkdir(analysis_directory);
end
video_name = [FilePath(1).name(1:21) ' ' num2str(frame_rate) '-fps'];
addpath(analysis_directory)

%% Make the Video:
origState = warning;
warning('off')

fig = figure; set(fig, 'pos',[300,190,1321,776], 'color','k');
currAxes.Visible = 'off';
hold on
v = VideoWriter([analysis_directory video_name ' uncompressed.avi'], 'Uncompressed AVI');
v.FrameRate = frame_rate;
open(v);
for datapoint = [1, start_frame+2:end_frame+1]
    set(fig, 'color','k');
    
    
    buffer_frame = uint8(zeros(image_size));
    %set the image
    for icam = 1:6
        %image
        buffer_frame(H_zone.(Alphabet(icam)), W_zone.(Alphabet(icam))) = ...
            ProcessedImages.(Alphabet(icam))(datapoint).data;
        %text labels
        if icam <= 3 
            text_image = Ltxt;
        elseif icam >= 5
            text_image = Rtxt;
        else
            text_image = Ttxt;
        end
        buffer_frame(H_zone.([Alphabet(icam) '_text']), W_zone.([Alphabet(icam) '_text'])) = text_image;
        % LIB
        buffer_frame(H_zone.([Alphabet(icam) '_LI']), W_zone.([Alphabet(icam) '_LI'])) = ...
                ProcessedImages.LIB(datapoint).data;
        
    end
%     % LIB
%     buffer_frame(H_zone.LIB, W_zone.LIB) = ProcessedImages.LIB(datapoint).data;
    %Cross:
    buffer_frame(H_zone.cross_txt, W_zone.cross_txt) = crosstext;
    % for the first frame, just black pixels
    if datapoint == 1
      buffer_frame = uint8(zeros(image_size));  
    end
    %display the 
    imshow(buffer_frame, 'border', 'tight')
    currAxes.Visible = 'off';
    % Collect image info for movie file
    f = getframe(fig);
    writeVideo(v, f)
    clf('reset')
end
close(v)
close all
warning(origState)
fprintf('\n Video Saved! \n')

%% Compress the video
old_vid = [analysis_directory video_name ' uncompressed.avi'];
new_vid = [analysis_directory video_name];

reader = VideoReader(old_vid);        
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
disp(' Done!')







    
    
    
    
    
    
    