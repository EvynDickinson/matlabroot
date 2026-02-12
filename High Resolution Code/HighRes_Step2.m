%    moveRawFiles

%% Parallelized version: 
clc; clear
% TODO -- update this to read from t
% he HighRes Experiments.xlsx file
% specifically: have it pick up the files that have 'R' for ready and then
% turn them to a 'Y' when it's done
% TODO: have this copy the tracking files into each of the folders! 

rootDir = getDataPath(5, 0, 'Select data location');
paths = getPathNames;
dateDir = selectFolder([rootDir, paths.courtship]);
rampName = selectFolder([rootDir paths.courtship  dateDir{1}],true,'Select Folders to Compile');

% option to change the save folder to a secondary location and not back to
% the server: 
if strcmp(questdlg('Do you want to save videos back to the same start folder?'), 'No')
    disp('select the location where you want to create a new folder to save the compiled videos')
    endDir = [uigetdir '/'];
else
    endDir = [rootDir paths.courtship];
end

startDir = [rootDir paths.courtship];

initial_var = who;
initial_var{end+1} = 'ramp';
initial_var{end+1} = 'initial_var';


%% convert files to a single video file (or an 8-minute segment)
disp('COMPILING IN PROGRESS')

for ramp = 1:size(rampName,2)
    clearvars('-except',initial_var{:})
    
    expName = rampName{ramp};

    baseDir = [startDir  dateDir{1} '\' expName '\'];
    vidBase = [endDir dateDir{1} '\' expName '\'];
    if ~exist(vidBase, 'dir')
        mkdir(vidBase)
    end
    
    % load parameters:
    load([baseDir, expName ' alignment table.mat'],'parameters')
    nVids = parameters.nVids;
    vROI = parameters.vROI;
    
    % find videos that have already been compiled and extract them from the processing list
    present_vid_dir = dir([vidBase 'compiled_video_*.avi']);
    nameList = {present_vid_dir(:).name};
    run_idx = true([nVids,1]);
    if ~isempty(nameList) % only run this if there are already compiled videos
        for i = 1:length(nameList)
            dummy = strsplit(nameList{i},{'_','.'});
            vidNum = str2double(dummy{3});
            run_idx(vidNum) = false;
        end
    end

    % Start parallel pool if it's not already running
    if isempty(gcp('nocreate'))
        if strcmpi(getenv('COMPUTERNAME'), 'SLEEPINGGIANT')
            parpool(22);
        else
            parpool;
        end
    end
    
    %  loop for parallel processing
    iStart = vROI(:,1); % file number that starts the compiling
    iEnd = vROI(:,2); % file number that ends the compiling list
    N = (iEnd-iStart)+1; % total number of files to compile
    hz = parameters.FPS; 
    videoList = (find(run_idx))';
    nVids = length(videoList);
    
    tic
    parfor idx = 1:nVids
        vid = videoList(idx);
        % Each compiled video is processed independently in parallel
        vidPath = [vidBase 'compiled_video_' num2str(vid)];
        % Convert matrix to video file
        v = VideoWriter([vidPath '.avi'], 'Motion JPEG AVI');
        v.Quality = 95;
        v.FrameRate = hz;
        open(v);
        % loop through all the files and compile them to the new video file
        for ii = 1:N(vid)
            ROI = iStart(vid):iEnd(vid); % current 'file' to compile
            i = ROI(ii);
            % Load data matrix
            data = load([baseDir, 'file' num2str(i) '.mat']); 
            writeVideo(v, data.data);
            % disp(['Finished ' 'file' num2str(i) '.mat']);
        end    
        close(v);
    end
    disp(['Finished ' expName ' videos'])
    disp([num2str(ramp) ' / ' num2str(size(rampName,2))])
    toc

end

ans = questdlg('All videos finished processing','','Okay','Okay');











