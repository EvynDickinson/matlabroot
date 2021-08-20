
%get base folder pathway
baseFolder = getCloudPath;

params.well_1 = 'Plant';
params.well_2 = 'Empty';
params.well_3 = 'Empty';
params.well_4 = 'Empty';

params.genotype = select_cross; % choose the fly genotypes
params.protocol = select_protocol;  % choose experiment protocol

% move params into structure:
params.date = dirName;
params.expID = videoNames;
params.num = num;


writeExptoExcel(params) % save the data to the experiement list

% resave with same name: 
save([baseFolder dirName '\' videoNames 'dataMat'])


%% 
baseFolder = getCloudPath;

for i = 1:num.vids
    filename = [baseFolder dirName '\analysis\' videoNames '_' num2str(i) ' data'];
    load(filename);
    
    params.well_1 = 'Plant';
    
    save(filename, 'videoData', 'params')
end


%% Fix the temp log information being added...
clear
% manually load the exp data file
% Fix the MAIN file
[~,baseFolder] = getCloudPath(1);

tempLogStart = [(1:num.vids)', tempLogStart];
save([baseFolder, '\', videoNames, 'dataMat.mat'])

% load each vid file & adjust the temp log:
tempLog = readmatrix([baseFolder '\' videoNames '_RampLog']);
roi = [tempLogStart(:,3), tempLogEnd(:,3)];
filepath = [baseFolder, '\analysis\' videoNames '_'];

for i = 1:num.vids
    load([filepath num2str(i) ' data.mat']);
    
    % update file
    loc = tempLog(:,1)>=roi(i,1) & tempLog(:,1)<=roi(i,2);
    videoData.tempLog = tempLog(loc,:);
    
    % save file
    save([filepath num2str(i) ' data.mat'], 'params', 'videoData')
end




























