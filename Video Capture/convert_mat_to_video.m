
% clear
% 
% tic
% 
% % load matrix data file
% baseFolder = 'F:\Evyn\DATA\09.23.2024\output_1\';
% load([baseFolder 'file1.mat']);
% 
% % convert matrix to video file
% v = VideoWriter([baseFolder 'file1.avi'],'Motion JPEG AVI');
% v.Quality = 95;
% open(v)
% writeVideo(v, data)
% close(v)
% 
% toc



%% How long can we write to buffer without crapping the memory?

rootDir = getDataPath(5, 2, 'Select location for data');
paths = getPathNames;
dateDir = selectFolder([rootDir, paths.courtship]);
rampName = selectFolder([rootDir paths.courtship  dateDir{1}]);
baseDir = [rootDir paths.courtship  dateDir{1} '\' rampName{:} '\'];
% baseFolder = [rootDir paths.courtship  dateDir{1} '\Video Testing\'];

% fold = 3;
% baseDir = [baseFolder 'output_' num2str(fold) '/'];
fileList = dir([baseDir 'file*.mat']);
 

% find difference in time between the video saves: 
timestamps = [];
for i = 1:length(temp)
    [h,m,s] = hms(temp(i).data.timestamp);
    timestamps(i) = seconds(duration(h,m,s));
end

y = diff(timestamps);

fig = getfig('',false,[1123 601]); 
scatter(1:length(y),y,50,Color('black'),'filled')
h_line(5,'red','--',1.5)
h_line(mean(y),'k','-',1)
formatFig(fig)
xlabel('Video Number')
ylabel('Video record and save time (s)')

save_figure(fig, [baseDir 'Video fragment length'],'-png');

%% convert files to videos
% This is something that can be parallized and run in multiples! %TODO
tic

% Loop through all the video files 
baseFolder = 'F:\Evyn\Courtship Tracking\09.26.2024\Berlin_courtship_F_LRR_caviar_lighttest_1\';
fileList = dir([baseFolder 'file*.mat']);
for i = 1:1%length(fileList)
    % lad data matrix
    data = load([baseFolder fileList(i).name],'data');
    
    % convert matrix to video file
    v = VideoWriter([baseFolder  fileList(i).name(1:end-3) 'avi'],'Motion JPEG AVI');
    v.Quality = 95;
    open(v)
    writeVideo(v, data.data)
    close(v)
    disp(['Finished ' fileList(i).name])
end
toc

%% convert files to a single video file (or an 8-minute segment)

newVidLength = 4; % four minutes
% need to know the total video length & the video frequency & individual
% vid lengths
tot_vids = length(fileList);
vid_length = 5;
compile_size = (newVidLength*60)/vid_length;
vid_rois = 1:compile_size:tot_vids;
if vid_rois(end)<tot_vids
    vid_rois = [vid_rois, tot_vids];
end
nVids = length(vid_rois)-1; % how many compiled videos there will be
vROI = [vid_rois(1),vid_rois(2)];
for i = 2:nVids
    vROI(i,:) = [vid_rois(i)+1, vid_rois(i+1)];
end

% This is something that can be parallized and run in multiples! %TODO
tic
for vid = 1:nVids
    % Loop through all the video files 
    vidPath = [baseDir  'compiled video_' num2str(vid)];
    parameters = {};
    % convert matrix to video file
    v = VideoWriter([vidPath '.avi'],'Motion JPEG AVI');
    v.Quality = 95;
    open(v)
    for i = vROI(vid,1):vROI(vid,2)
        % lad data matrix
        data = load([baseDir, 'file' num2str(i) '.mat']'); 
        parameters{end+1} = data.timestamp;
        writeVideo(v, data.data)
        disp(['Finished ' 'file' num2str(i) '.mat'])
    end
    close(v)
    save([vidPath ' parameters.mat'],'parameters')
end
disp('Finished all videos')
toc







































