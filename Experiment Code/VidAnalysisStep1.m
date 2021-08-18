
clear; close all

%% Select videos to load and analyze

%get base folder pathway
switch getenv('COMPUTERNAME')
    case 'DENALI'
        baseFolder = 'E:\My Drive\Jeanne Lab\DATA\';
    case 'TOGIAK'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
end

% Select the date of the experiment to analyze
dirc = dir(baseFolder);
dirc = flip(dirc(find(~cellfun(@isdir,{dirc(:).name}))));
folderNames = ['Today', {dirc(:).name}];
indx = listdlg('ListString', folderNames, 'SelectionMode', 'Single');
if strcmpi(folderNames{indx}, 'Today')==true
    dir_sel = strrep(datestr(datetime,'mm-dd-yyyy'),'-','.');
else
    dir_sel = folderNames{indx};
end
folder = fullfile(baseFolder, dir_sel);

% Select the complete experiments to process
list_dirs = dir([folder, '\*.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
indx = listdlg('ListString', expNames, 'SelectionMode', 'Multiple');

% Make new analyzed file directory
analysisDir = [folder '\analysis\']; 
if ~isfolder(analysisDir); mkdir(analysisDir); end


%% Process videos
warning off

for Exp = 1:length(indx)
    
    % load matlab file for experiment
    expData = load([folder '\' expNames{indx(Exp)} 'dataMat.mat']);
    vidName = expData.videoNames;
    
    % load and process videos:
    for vid = 1:expData.num.vids
        
        % video parameters
        vidPath = fullfile(folder, [vidName '_' num2str(vid) '.avi']);
        movieInfo = VideoReader(vidPath); %read in video

        % extract parameters from video
        nframes = movieInfo.Duration * movieInfo.FrameRate;
        height = movieInfo.Height;
        width = movieInfo.Width;
        occCumSum = zeros(height, width); %setup blank occupation
        n_tot = nframes;
        
        % save processing params for this experiment:
        nWells = 1;
        pixel_thresh = 105;
        cluster_thresh = 4;
        if vid == 1
            % demo image
            demoImg = rgb2gray(read(movieInfo,1));
            demoAdj = imopen(demoImg>pixel_thresh,strel('disk',cluster_thresh));
            f = figure; imshowpair(demoImg, demoAdj,'montage'); title('Thresholding test')
            uiwait(f)

            % COLORTHRESHOLD
            switch questdlg('Use default image threshold?')
               case 'No' 
                   f = figure; imtool(demoImg); uiwait(f)
                   pixel_thresh = str2double(cell2mat(inputdlg('Light pixel threshold?')));
                   demoAdj = imopen(demoImg>pixel_thresh,strel('disk',4));
                   f = figure; imshowpair(demoImg, demoAdj,'montage'); uiwait(f)
                   switch questdlg('New threshold acceptable?')
                        case 'No' 
                            return
                   end
               case 'Cancel'
                   return
           end
           % MASK OPTIONS
           switch questdlg('Load previous arena mask?')
               case 'Yes'
                   [file,path] = uigetfile('*.mat');
                   a = load(fullfile(path, file));
                   mask = a.mask;
               case 'No'
                   mask = drawArenaMask(demoImg,nWells);
                   save([analysisDir 'Mask ' vidName], 'mask');
               case 'Cancel'
                   return
           end
           % Test all image processing:
           demoAdj = imopen(demoImg>pixel_thresh,strel('disk',4));
           demoAdj(mask) = 0;
           f = figure; imshowpair(demoImg, demoAdj,'montage'); title('Thresholding test')
           uiwait(f)
           switch questdlg('Processing acceptable?')
                case 'No' 
                    return
           end
           clear file path demoImg demoAdj f a
        end
        
        % process video frames
        h = waitbar(0,...
        ['reading in frames from vid ' num2str(vid) '/' num2str(expData.num.vids)]);
        for i = 1:n_tot
            % Image processing:
            Img = processOccImg(read(movieInfo,i), pixel_thresh, cluster_thresh);
            Img(mask) = 0; % mask for arena size
            
            % Update the occupancy count
            currOcc = gather(Img/sum(sum(Img))); %probabilty for this frame (1 across total image)
            occCumSum = (occCumSum + currOcc);
            waitbar(i/n_tot,h)
        end
        close(h)
        
        % Group the video and temp data:
%         videoData....
        
        
    end
    
    
    
end


fig = figure; set(fig, 'color', 'k')
h = heatmap(gather(occCumSum)); 
colormap hot;
set(h,'gridvisible', 'off')
ax = gca;
ax.XDisplayLabels = nan(size(ax.XDisplayData));
ax.YDisplayLabels = nan(size(ax.YDisplayData));
set(ax,'FontColor', 'w', 'FontName', 'Arial');
title('Spatial occupation probability');


%% Set up basics with a GPU:

% gpu = gpuDevice();

for n = 1:length(indx)
    
    vid = indx(n); %select video file from list
    vidPath = fullfile(folder, vidNames{vid}); %full video path
    movieInfo = VideoReader(vidPath); %read in video
    
    %extract parameters from video
    nframes = movieInfo.Duration * movieInfo.FrameRate;
    height = movieInfo.Height;
    width = movieInfo.Width;
    occCumSum = zeros(height, width); %setup blank occupation
    n_tot = nframes;
    
%     currOcc = nan(height, width, nframes);
    
    h = waitbar(0,...
        ['reading in frames from vid ' num2str(n) '/' num2str(length(indx))]);
    % Read and process the videos frame-by-frame
    tic
    for i = 1:n_tot
        
     % IMAGE PROCESSING: TODO move to a unique function
        IMoriginal = read(movieInfo,i); %convert to greyscale
%         imshow(IMoriginal)
        imGPUoriginal = gpuArray(IMoriginal);
        imGPUgray = rgb2gray(imGPUoriginal); % convert to greyscale
%         imtool(gather(imGPUgray)) % TODO: make this a dynamic option for
%         each chunk of data
        imflyGPU = imGPUgray>90; %mask for everything that is a fly color or lighter        
%         imshow(imflyGPU)       
        imFlyMask = imopen(imflyGPU,strel('disk',5));
%         imshow(imFlyMask)
        % Mask the probability image:
        Img = imFlyMask;
        Img(mask) = 0;
%         imshow(Img)
        
                   
        

        currOcc = gather(Img/sum(sum(Img))); %probabilty for this frame (1 across total image)
        occCumSum = (occCumSum + currOcc);
        
%         fig = figure; set(fig, 'color', 'k')
%         h = heatmap(gather(currOcc)); 
%         colormap hot;
%         set(h,'gridvisible', 'off')
%         ax = gca;
%         ax.XDisplayLabels = nan(size(ax.XDisplayData));
%         ax.YDisplayLabels = nan(size(ax.YDisplayData));
%         set(ax,'FontColor', 'w', 'FontName', 'Arial');
%         title('Spatial occupation probability');

        waitbar(i/n_tot,h)
    end
    close(h)
    toc
    clear IMoriginal imGPUoriginal imGPUgray imFlyMask Img imflyGPU
    % calculate the probability of a fly being in a given location
    prob = occCumSum / i;
        
%         OccSum = (OccSum + gather(currOcc))/i;
%     
%     % save specific vid data
%     data.tot(:,:,n) = tot_occ;
%     data.prob(:,:,n) = currProb;

end
fprintf('\n All loaded! \n')  







fig = figure; set(fig, 'color', 'k')
h = heatmap(prob); 
colormap hot;
set(h,'gridvisible', 'off')
ax = gca;
ax.XDisplayLabels = nan(size(ax.XDisplayData));
ax.YDisplayLabels = nan(size(ax.YDisplayData));
set(ax,'FontColor', 'w', 'FontName', 'Arial');
title('Spatial occupation probability');



























