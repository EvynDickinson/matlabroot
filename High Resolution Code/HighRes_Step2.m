



%% Parallelized version: 
clc; clear
% TODO -- update this to read from the Courtship Experiments.xlsx file
rootDir = getDataPath(5, 2, 'Select Date of data to process');
paths = getPathNames;
dateDir = selectFolder([rootDir, paths.courtship]);
rampName = selectFolder([rootDir paths.courtship  dateDir{1}],true,'Select Folders to Compile');

initial_var = who;
initial_var{end+1} = 'ramp';
initial_var{end+1} = 'initial_var';


%% convert files to a single video file (or an 8-minute segment)

for ramp = 1:size(rampName,2)
    clearvars('-except',initial_var{:})
    
    expName = rampName{ramp};
    baseDir = [rootDir paths.courtship  dateDir{1} '\' expName '\'];
    
    % load parameters:
    load([baseDir, expName ' alignment table.mat'],'parameters')
    nVids = parameters.nVids;
    vROI = parameters.vROI;
    
    % Start parallel pool if it's not already running
    if isempty(gcp('nocreate'))
        parpool(11);
    end
    
    %  loop for parallel processing
    iStart = vROI(:,1);
    iEnd = vROI(:,2);
    N = (iEnd-iStart)+1;
    hz = parameters.FPS;
    
    tic
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
    disp(['Finished ' expName ' videos'])
    disp([num2str(ramp) ' / ' num2str(size(rampName,2))])
    toc

end

ans = questdlg('All videos finished processing','','Okay','Okay');











