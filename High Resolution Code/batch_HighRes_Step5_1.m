

%% LOAD data
clear; clc;
path = getDataPath(6,0);
baseFolder = [path,'Trial Data/'];
trials_list = dir(baseFolder);
trials_list = {trials_list(3:end).name}; % full list of trials to process

processed_log = {};

for this_trial = 1:length(trials_list)
    trialDir = [baseFolder, trials_list{this_trial} '/'];
    
    if ~isfile([trialDir 'post-5.1.1 data.mat'])
        disp(trials_list{this_trial})
    end
end


for this_trial = 1:length(trials_list)
    trialDir = [baseFolder, trials_list{this_trial} '/'];

    if isfile([trialDir 'post-5.1.1 data.mat'])
        continue
    % else
    %     break
    end

    try runBatch_HighRes_Step5_1(trialDir);
    catch 
        continue
    end
    disp([num2str(this_trial) '/' num2str(length(trials_list)) ' ' trials_list{this_trial}])
    processed_log{this_trial} =  trials_list{this_trial};
end

% [processed_log', trials_list']

