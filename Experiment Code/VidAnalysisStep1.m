% PROCESS EACH VIDEOO & EXTRACT THE PIXEL COUNT PROBABILITY
clear; close all
closepreview
%% Select videos to load and analyze

%get base folder pathway
[baseFolder, folder] = getCloudPath(1);

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

    % load the temperature log for the experiment
    tempLog = readmatrix([folder '\' vidName '_RampLog']);

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
        for ii = 1:4; wellList{ii} = expData.params.(['well_' num2str(ii)]); end 
        nWells = 4-sum(strcmpi(wellList, 'Empty'));

        pixel_thresh = 60; %105;
        cluster_thresh = 5;
        if vid == 1 % define and check image processing parameters

            % demo image
            demoImg = read(movieInfo,1);
            demoAdj = processOccImg(demoImg, false, pixel_thresh, cluster_thresh);
            f = figure; imshowpair(demoImg, demoAdj,'montage'); title('Thresholding test')
            
            % COLORTHRESHOLD
            switch questdlg('Use default image threshold?')
                case 'Yes'
                    close(f)
                case 'No'
                    close(f)
                    badThresh = true;
                    while badThresh
                       testImg = processOccImg(demoImg, false, pixel_thresh, cluster_thresh,1);
                       f = imtool(testImg); uiwait(f)
                       pixel_thresh = str2double(cell2mat(inputdlg('Light pixel threshold?')));
                       demoAdj = processOccImg(demoImg, false, pixel_thresh, cluster_thresh);
                       f = figure; imshowpair(demoImg, demoAdj,'montage')
                       switch questdlg('New threshold acceptable?')
                           case 'Yes'
                               close(f)
                               badThresh = false;
                           case 'No'
                               close(f)
                           case 'Cancel'
                               return
                       end
                    end
                case 'Cancel'
                   return
            end
           
           % MASK OPTIONS
           switch questdlg('Load previous arena mask?')
               case 'Yes'
                   [file,path] = uigetfile([folder '\analysis\*mask*.mat']);
                   load(fullfile(path, file));
               case 'No'
                   mask = drawArenaMask(demoImg,nWells);
                   save([analysisDir 'Mask ' vidName], 'mask');
               case 'Cancel'
                   return
           end
           % Test all image processing:
           demoAdj = processOccImg(demoImg, false, pixel_thresh, cluster_thresh);
           demoAdj(mask) = 0;
           f = figure; imshowpair(demoImg, demoAdj,'montage'); title('Full processing preview')
           switch questdlg('Processing acceptable?')
               case 'Yes'
                   params.pixel_thresh = pixel_thresh;
                   params.cluster_thresh = cluster_thresh;
                   close(f)
               case 'No'
                   return
               case 'Cancel'
                   return
           end
           clear file path demoImg demoAdj f a badThresh
        end

        % process video frames
        h = waitbar(0,...
        ['reading in frames from vid ' num2str(vid) '/' num2str(expData.num.vids)]);
        for i = 1:n_tot
            % Image processing:
            Img = processOccImg(read(movieInfo,i), true, pixel_thresh, cluster_thresh);
            Img(mask) = 0; % mask out arena and food wells

            % Update the occupancy count
            currOcc = gather(Img/sum(sum(Img))); %probabilty for this frame (1 across total image)
            occCumSum = (occCumSum + currOcc);
            waitbar(i/n_tot,h)
        end
        close(h)

        % add data summary into data structure
        videoData.currOcc = currOcc;
        videoData.occ_Prob = occCumSum;

        % Group the video and temp data:
        videoData.temprange = [expData.tempLogStart(vid,4), expData.tempLogEnd(vid,4)];

        % Extract the appropriate chunk of the data log
        roi = [expData.tempLogStart(vid,3), expData.tempLogEnd(vid,3)];
        loc = tempLog(:,1)>=roi(1) & tempLog(:,1)<=roi(2);
        videoData.tempLog = tempLog(loc,:);
        fprintf(['\n' num2str(length(tempLog(loc,:)))])
        videoData.name = [vidName '_' num2str(vid)];

        % Save the VideoData into the analysis folder
        params = expData.params;
        save([analysisDir videoData.name ' data'], 'videoData', 'params')
    end
    % TODO: add an analyzed checkbox to the excel sheet here
end

warning on

% fig = figure; set(fig, 'color', 'k')
% h = heatmap(gather(occCumSum));
% colormap hot;
% set(h,'gridvisible', 'off')
% ax = gca;
% ax.XDisplayLabels = nan(size(ax.XDisplayData));
% ax.YDisplayLabels = nan(size(ax.YDisplayData));
% set(ax,'FontColor', 'w', 'FontName', 'Arial');
% title('Spatial occupation probability');



















