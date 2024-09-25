

%% Pull experiment values from the acquired data
% 
% rootDir = getDataPath(5, 2, 'Select location for data');
% paths = getPathNames;
% dateDir = selectFolder([rootDir, paths.courtship]);
% rampName = selectFolder([rootDir paths.courtship  dateDir{1}]);
% baseDir = [rootDir paths.courtship  dateDir{1} '\' rampName{:} '\'];
% 
% fileList = dir([baseDir 'file*.mat']);
% tic
% temp = [];
% for i = 1:length(fileList)
%     temp(i).data = load([baseDir, 'file' num2str(i) '.mat'],'timestamp');
% end
% toc
% 
% 
% %% Plot the fragment lengths
% 
% % find difference in time between the video saves: 
% timestamps = [];
% for i = 1:length(temp)
%     [h,m,s] = hms(temp(i).data.timestamp);
%     timestamps(i) = seconds(duration(h,m,s));
% 
% end
% 
% y = diff(timestamps);
% 
% fig = getfig('',false,[1123 601]); 
% scatter(1:length(y),y,50,Color('black'),'filled')
% h_line(5,'red','--',1.5)
% h_line(mean(y),'k','-',1)
% formatFig(fig)
% xlabel('Video Number')
% ylabel('Video record and save time (s)')
% 
% save_figure(fig,[baseDir, 'Video fragment length'],'-png')
% 
% %% convert single files to videos
% % % This is something that can be parallized and run in multiples! %TODO
% tic
% 
% % Loop through all the video files 
% % baseFolder = 'F:\Evyn\DATA\09.23.2024\output_1\';
% for i = 1:length(fileList)
%     % lad data matrix
%     data = load([baseDir fileList(i).name],'data');
% 
%     % convert matrix to video file
%     v = VideoWriter([baseDir  fileList(i).name(1:end-3) 'avi'],'Motion JPEG AVI');
%     v.Quality = 95;
%     open(v)
%     writeVideo(v, data.data)
%     close(v)
%     disp(['Finished ' fileList(i).name])
% end
% toc


%% Parallelized version: 
clc; clear
% TODO -- update this to read from the Courtship Experiments.xlsx file
rootDir = getDataPath(5, 2, 'Select location for data');
paths = getPathNames;
dateDir = selectFolder([rootDir, paths.courtship]);
rampName = selectFolder([rootDir paths.courtship  dateDir{1}]);
A = strsplit(rampName{:}, '_');
expName = strjoin(A(1:end-1), '_');

baseDir = [rootDir paths.courtship  dateDir{1} '\' rampName{:} '\'];

% load parameters:
load([baseDir, expName ' alignment table'],'parameters')

nVids = parameters.nVids;
vROI = parameters.vROI;

% % Create a 'compile' list
% compileLength = 3; %  minutes
% compile_size = floor((compileLength*60)/parameters.fragmentduration);
% vid_rois = 0:compile_size:parameters.numFrag;
% vid_rois(1) = 1;
% if vid_rois(end) < parameters.numFrag
%     vid_rois = [vid_rois, parameters.numFrag];
% end
% nVids = length(vid_rois)-1; % how many compiled videos there will be
% vROI = [vid_rois(1), vid_rois(2)];
% for i = 2:nVids
%     vROI(i,:) = [vid_rois(i)+1, vid_rois(i+1)];
% end


% Start parallel pool if it's not already running
if isempty(gcp('nocreate'))
    parpool;
end

tic
%  loop for parallel processing
iStart = vROI(:,1);
iEnd = vROI(:,2);
N = (iEnd-iStart)+1;
hz = parameters.FPS;

parfor vid = 1:nVids
    % Each compiled video is processed independently in parallel
    vidPath = [baseDir 'compiled_video_' num2str(vid)];
    % Convert matrix to video file
    v = VideoWriter([vidPath '.avi'], 'Motion JPEG AVI');
    v.Quality = 95;
    v.FrameRate = hz;
    open(v);

    for ii = 1:N(vid)
        ROI = iStart(vid):iEnd(vid);
        i = ROI(ii);
        % Load data matrix
        data = load([baseDir, 'file' num2str(i) '.mat']); 
        writeVideo(v, data.data);
        disp(['Finished ' 'file' num2str(i) '.mat']);
    end
    
    close(v);
end
disp('Finished all videos');
toc







































