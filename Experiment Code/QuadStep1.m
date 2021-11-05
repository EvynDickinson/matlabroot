% LOAD VIDEOS FOR AN EXPERIMENT AND CROP THEM INTO THE FOUR INDEPENDENT
% QUADRANTS / ARENAS FOR TRACKING AND POST PROCESSING
clear; clc  

%% Select data to split:
%get base folder pathway
[baseFolder, folder] = getCloudPath(2); 

% Select the complete experiments to process
list_dirs = dir([baseFolder folder, '\*.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
expName = expName(1:end-1);
clear expNames
% Pull fly summary sheet information on selected experiment
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;

% Make a new video folder:
videosDir = [baseFolder folder '\Split Videos\'];
if ~isfolder(videosDir); mkdir(videosDir); end

%% Assign quadrants to split
%Load vid template and manually align vid split
params = load([baseFolder folder '\' expName ' parameters.mat']);
nvids = params.parameters.numVids;
movieInfo = VideoReader([baseFolder folder '\' expName '_1.avi']); %read in video
demoImg = rgb2gray(read(movieInfo,1));
% figure; imshow(demoImg)
% [J,rect] = imcrop(demoImg);

% use fixed crop size for all images and videos
cropWidth = 946;
width = 1920;
% RGB = insertShape(demoImg,'rectangle',[0,0,cropWidth, cropWidth]);
% imshow(RGB)
% cropping indexes
roiA = [width-cropWidth, width, 1, cropWidth];
roiB = [1, cropWidth, 1, cropWidth];
roiC = [width-cropWidth, width, width-cropWidth, width];
roiD = [1, cropWidth, width-cropWidth, width];

% crop the image (will be used for later processing)
imgA = demoImg(roiA(1):roiA(2), roiA(3):roiA(4));
imgB = demoImg(roiB(1):roiB(2), roiB(3):roiB(4));
imgC = demoImg(roiC(1):roiC(2), roiC(3):roiC(4));
imgD = demoImg(roiD(1):roiD(2), roiD(3):roiD(4));

% Preview selecion:
nrow = 2; ncol = 4;
fig = figure; set(fig, 'pos', [267 233 959 420], 'color', 'k');
hold on;
subplot(nrow, ncol, 1) %Arena B
imshow(imgB)
subplot(nrow, ncol, 2) %Arena D
imshow(imgD)
subplot(nrow, ncol, 5) %Arena A
imshow(imgA)
subplot(nrow, ncol, 6) %Arena C
imshow(imgC)
subplot(nrow, ncol, [3:4,7:8]) %Arena B
imshow(demoImg)
save_figure(fig, [baseFolder folder '\' expName ' Video Cropping'], '-pdf')

switch questdlg('Acceptable quadrants?')
    case 'Yes'
    case 'No'
        return
    case 'Cancel'
        return
end

%% Split the videos
switch questdlg('Group AB or Group CD?','', 'AB', 'CD','Cancel') %
    case 'AB'
        tic
        for vid = 1:nvids
            movieInfo = VideoReader([baseFolder folder '\' expName '_' num2str(vid) '.avi']); %read in video
            nframes = movieInfo.duration;
            h = waitbar(0,['Cropping video ' num2str(vid) '/' num2str(nvids)]);
            % initiate arena A writer
            videoA = [baseFolder folder '\Split Videos\' expName 'A_' num2str(vid) '.avi'];
            vidA = VideoWriter(videoA, 'Motion JPEG AVI');
            vidA.Quality = 75;
            vidA.FrameRate = 6;
            open(vidA);
            % initiate arena B writer
            videoB = [baseFolder folder '\Split Videos\' expName 'B_' num2str(vid) '.avi'];
            vidB = VideoWriter(videoB, 'Motion JPEG AVI');
            vidB.Quality = 75;
            vidB.FrameRate = 6;
            open(vidB);
            fig = figure; set(fig, 'color', 'k','visible', 'off');

            for frame = 1:nframes
                %Extract and crop frames
                demoImg = rgb2gray(read(movieInfo,frame));
                imgA = demoImg(roiA(1):roiA(2), roiA(3):roiA(4));
                imgB = demoImg(roiB(1):roiB(2), roiB(3):roiB(4));

                %Save A frames 
                imshow(imgA)
                f = getframe(fig);
                writeVideo(vidA, f)  

                %Save B frames
                imshow(imgB)
                f = getframe(fig);
                writeVideo(vidB, f)  

                % update the waitbar every 10 frames
                if rem(frame,10) == 0
                    waitbar(frame/nframes,h)
                end   
            end
            close(vidA); close(vidB)
             close(h); close(fig)   
        end
        toc
    case 'CD'
        tic
        for vid = 1:nvids
            movieInfo = VideoReader([baseFolder folder '\' expName '_' num2str(vid) '.avi']); %read in video
            nframes = movieInfo.duration;
            h = waitbar(0,['Cropping video ' num2str(vid) '/' num2str(nvids)]);

            % initiate arena C writer
            videoC = [baseFolder folder '\Split Videos\' expName 'C_' num2str(vid) '.avi'];
            vidC = VideoWriter(videoC, 'Motion JPEG AVI');
            vidC.Quality = 75;
            vidC.FrameRate = 6;
            open(vidC);
            % initiate arena D writer
            videoD = [baseFolder folder '\Split Videos\' expName 'D_' num2str(vid) '.avi'];
            vidD = VideoWriter(videoD, 'Motion JPEG AVI');
            vidD.Quality = 75;
            vidD.FrameRate = 6;
            open(vidD);

            fig = figure; set(fig, 'color', 'k','visible', 'off');
            for frame = 1:nframes
                %Extract and crop frames
                demoImg = rgb2gray(read(movieInfo,frame));
                imgC = demoImg(roiC(1):roiC(2), roiC(3):roiC(4));
                imgD = demoImg(roiD(1):roiD(2), roiD(3):roiD(4));
  
                %Save C frames
                imshow(imgC)
                f = getframe(fig);
                writeVideo(vidC, f)        
    
                %Save D frames
                imshow(imgD)
                f = getframe(fig);
                writeVideo(vidD, f)  

                % update the waitbar every 10 frames
                if rem(frame,10) == 0
                    waitbar(frame/nframes,h)
                end   
            end
            close(vidC); close(vidD)
             close(h); close(fig)   
        end
        toc
end

