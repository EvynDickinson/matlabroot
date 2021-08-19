
% VidAnalysisStep2 --> group flies and plot graphs of temp & occupancy
% follows VidAnalysisStep1
clear
%% Select data to group

% TODO: add dynamic reading of experiment xl file 
[~, folder] = getCloudPath(1);
dirc = [folder '\analysis'];

% offer exp. options:
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
% ... TODO ...

expName = 'PlantFoodTest_N1'; % TODO make this dynamic
num.vids = size(dir([dirc '\' expName '*.mat']),1);
% Load occupation data from each video:
for i = 1:num.vids
    filename = [dirc '\' expName '_' num2str(i) ' data.mat'];
    load(filename); % should load: params & videoData
    
    % Extract information to add to grouped set
    temp(i) = mean(videoData.tempLog(:,2));
    data(:,:,i) = videoData.occ_Prob;
    % -----------------------------------------
end


figure;
plot(temp)

length(videoData.tempLog)
