


%% load in spike time data
clear

spikedatafolder = 'JJ222';

serverPath  = 'S:\Maddie\';
% serverPath  = '/Volumes/shared/Maddie/'; % maddies laptop

basepath = [serverPath 'data from ephys rig\'];
fullpath = [basepath spikedatafolder,'/'];
foldercontents = dir([fullpath,'AllSpikeTimes*.mat']); %search for only files that start with AllSpikeTimes
nTrials = length(foldercontents);
saveDir = [basepath 'spike rate data/']; % folder path for processd spike rate data (needs to match the folder in Finder)

% odor parameters
odorON = 4; % odor start time 
odorOFF = 6; % odor end time

% Initalize new variable to put the spike data in
spiketimeData = nan(1,nTrials);

for ii = 1:nTrials
    
    % pull out the names of the data files in the folder 
    trial_name = foldercontents(ii).name;
    % pull out the trial number from the name
    splitName = strsplit(trial_name,'_');
    trial_number = str2double(splitName{end}(1:end-4));
    % full data location pathway
    filePath = [fullpath trial_name];
    
    % load data from the pathway
    % TODO 11.20: figure out why it is skipping data and not finding the 'allSpikeTimeData' variable    
    try tempData = load(filePath);
        data = tempData.allSpikeTimeData.spikeTimesPrecise;
    catch
        disp(['skipping ' trial_name])
        continue
    end
    % Find avg spike rate during odor period
    odorSpikes = data>=odorON & data<=odorOFF;
    spikeRate = sum(odorSpikes)/(odorOFF-odorON);

    % save spike rate into spike time data for this trial
    spiketimeData(trial_number) = spikeRate;

    % update on processed files
    disp(['Finished ' num2str(ii) '/' num2str(nTrials) ' trial ' num2str(trial_number)])
end
disp('All data loaded')

%% Load in Kristyn imaging data
imaging_path =[serverPath 'Kristyn Data/Mean_PNresponse_duringodor.mat']; % hard coded from laptop
imaging_data = load(imaging_path); % load in the data from the above path

% DM1 row 10
% VC1 row 24
% VC2 row 25
% VA4 row 19
% DP1m = row 16
% Vl2p = row 27

% find the row in the imaging data that corresponds to the specific cell type
switch spikedatafolder
    case 'MA15'
        imaging_row = 10;
        cell_type = 'DM1';
        % Spike trials that correspond to the following odor order
        trial_loc = [1,5; 6,10; 11,15; 16,20; 21,25;...
            26,30; 31,35; 36,40; 41,45; 46,50; 51,55; 56,60; 61,65];
        % spike trials order
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But'};
        legend_gap = 1;
    case 'MA16'
        imaging_row = 19;
        cell_type = 'VA4';
        trial_loc = [82,86; 77,81; 72,76; 22,26; 27,31; 32,36; 37,41; 42,46; 47,51; 52,56; 57,61; 62,66;];
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But'};
        legend_gap = 1;
    case 'MA17'
        imaging_row = 19;
        cell_type = 'VA4';
        trial_loc = [1,5; 6,10; 13,17; 18,22; 23,27; 28,32; 33,37; 40,44; 45,49; 50,54; 55,59; 60,66];
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But'};
        legend_gap = 1;
    case 'MA18'
        imaging_row = 19;
        cell_type = 'VA4';
        trial_loc = [9,13; 14,18; 19,23; 24,28; 29,33; 34,38];
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc'};
        legend_gap = 1;
    case 'MA19'
        imaging_row = 10;
        cell_type = 'DM1';
        trial_loc = [39,43; 44,48; 49,53];
        spike_order = {'ACV', 'AA', 'EA'};
        legend_gap = 1;
     case 'MA20'
        imaging_row = 10;
        cell_type = 'DM1';
        trial_loc = [2,6; 7,11; 12,16; 17,21];
        spike_order = {'ACV', 'AA', 'MA', 'ISOPENT'};
        legend_gap = 1;
     case 'MA21'
        imaging_row = 10;
        cell_type = 'DM1';
        trial_loc = [2,6; 7,11; 12,16];
        spike_order = {'ACV', 'AA', 'MA'};
        legend_gap = 1;
    case 'MA22'
        imaging_row = 10;
        cell_type = 'DM1';
        trial_loc = [79,83; 23,27; 68,72; 28,32; 33,37; 38,42; 43,47; 48,52; 53,57; 58,62; 51,55; 63,67];
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But'};
        legend_gap = 3;
    case 'MA23'
        imaging_row = 19;
        cell_type = 'VA4';
        trial_loc = [2,6; 7,11; 12,16; 17,21; 22,26; 27,31; 32,36; 37,41; 42,46; 47,51; 52,56; 57,61];
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But'};
        legend_gap = 3;
    case 'MA25'
        imaging_row = 10;
        cell_type = 'DM1';
        trial_loc = [1,5; 6,10; 12,16; 17,21; 22,26; 27,32; 63,67; 33,36; 38,42; 43,47; 48,52; 53,57];
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But'};
        legend_gap = 3;
    case 'JJ221'
        imaging_row = 27;
        cell_type = 'Vl2p';
        trial_loc = [24,28; 29,43];
        spike_order = {'ACV', 'AA'};
        legend_gap = 1;
    case 'JJ222'
        imaging_row = 27;
        cell_type = 'V12p';
        trial_loc = [91,95; 96,97; 81,85; 37,43; 59,63; 44,48; 49,53; 54,58; 86,90; 64,70; 71,75; 76,80; 16,26; 27,31];
        spike_order = {'ACV', 'AA', '2BUT', 'EA', 'PentAc', 'EthBut', 'MA', 'ISOPENT', '2,3-But', 'EHex', 'BENZ', '3OCT','ACV','AA'};
        legend_gap = 0.5;
    case 'JJ223'
        imaging_row = 16;
        cell_type = 'DP1m';
        trial_loc = [6,10; 11,20; 21,25; 26,30; 31,35];
        spike_order = {'ACV', 'AA', '2BUT', 'EA', 'PentAc'};
        legend_gap = 3;
end


% extract the specific data for this cell from the loaded imaging data
df_f = [imaging_data.Mean_PNresponse_duringodor{imaging_row,:}]; % pull the data into a matrix 

% Kristyn's odor order: 
dff_order = {'ACV','AA', 'BENZ', 'ISOPENT', 'EA', 'EthBut', 'PentAc', '2BUT', '3OCT', 'EHex', '2,3-But', 'MA'};
nOdors = length(dff_order); % number of odors in the imaging data

% reorder spike data to match the df_f data and find the average for each
avg_spike_rate = nan([1,nOdors]); % initialize empty variable for averages
for ii = 1:length(dff_order)
    match_odor = dff_order{ii}; % what odor in kristyns data you are matching
    idx = find(strcmp(match_odor, spike_order)); % find the location of the search odor in the spike rate odor order data
    if isempty(idx)
        continue
    end
    % pull out the trial numbers that match the odor of the correct order
    idx_trials = [];
    for odor = 1:length(idx)
        odor_loc = idx(odor);
        new_trials = trial_loc(odor_loc,1):trial_loc(odor_loc,2);
        idx_trials = [idx_trials, new_trials]; % find the trial numbers that correspond to this odor
    end
    avg_spike_rate(ii) = mean(spiketimeData(idx_trials),'omitnan'); % find the average spike rate for the odor trials
end

%% Save the processed data

save([saveDir, spikedatafolder '_avg_spike_rate.mat'], 'avg_spike_rate', 'dff_order', 'df_f', 'spikedatafolder');
disp('Saved data file')

%% Plot the comparison figure:


cmap = turbo(nOdors); % select colors for each odor from a color map
 
fig = figure; set(fig, 'color', 'w'); hold on
for ii = 1:nOdors
    scatter(df_f(ii),avg_spike_rate(ii),105,cmap(ii,:),"filled")
end
% format the figure
ylabel('mean spike rate (hz)')
xlabel('mean df/F')
set(gca, 'fontsize',18,'fontname', 'calibri')
% find the x limits
xlims = xlim;
xlim([xlims(1), xlims(end)+legend_gap])
legend(dff_order,'box', 'off','FontSize',15)
title_str =(['df/F vs spike rate in ' cell_type]);
title(title_str)
 

