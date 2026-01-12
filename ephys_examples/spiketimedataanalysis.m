
% This script extracts spike rate averages during the odor delivery period
% for all trials related to a given cell

% the output is a matlab file with the spike average data organized by the
% odor order that Kristyn used in her imaging data so that the spike data
% can be compared to the imaging data by odor

%% load in spike time data
clear

% input the name of the cell that will be processed:
spikedatafolder = 'MA22';

% organize pathways to data folders
serverPath = '/Volumes/shared/Maddie/'; % pathway to maddie's folder on the server
basepath = [serverPath 'data from ephys rig/']; % pathway to the folder that holds the processed cell data
fullpath = [basepath spikedatafolder,'/']; % combined specific data folder and base path to get full folder pathway
foldercontents = dir([fullpath,'AllSpikeTimes*.mat']); %search for only files that start with AllSpikeTimes in the selected folder
nTrials = length(foldercontents); % total number of spike files to load and process
saveDir = [basepath 'spike rate data/']; % saved OUTPUT folder path for processd spike rate data (needs to match the folder in Finder)

% odor timing parameters
odorON = 4; % odor start time 
odorOFF = 6; % odor end time

% Initalize new variable to load with the spike trial data
spiketimeData = nan(1,nTrials);

for i = 1:nTrials % repeat for each of the odor trials in the folder
    
    % pull out the names of the data files in the folder 
    trial_name = foldercontents(i).name;
    % pull out the trial number from the name
    splitName = strsplit(trial_name,'_');
    trial_number = str2double(splitName{end}(1:end-4)); % this keeps the trial number aligned by actual number not just the loaded in order
    % full data location pathway
    filePath = [fullpath trial_name];
    
    % load data from the pathway 
    try tempData = load(filePath);
        data = tempData.allSpikeTimeData.spikeTimesPrecise; % each spike time is the timepoint at which it occurred 
    catch
        disp(['skipping ' trial_name])
        continue
    end
    % Find avg spike rate during odor period
    odorSpikes = data>=odorON & data<=odorOFF; % logical vector of all the spikes within the given time frame
    spikeRate = sum(odorSpikes)/(odorOFF-odorON); %total number of spikes divided by the duration of the odor stim 

    % save spike rate into spike time data for this trial
    spiketimeData(trial_number) = spikeRate;

    % update on processed files
    disp(['Finished ' num2str(i) '/' num2str(nTrials) ' trial ' num2str(trial_number)])
end
disp('All data loaded')

%% Load in Kristyn imaging data
imaging_path = [serverPath 'Kristyn Data/Mean_PNresponse_duringodor.mat']; % hard coded from laptop
imaging_data = load(imaging_path); % load in the data from the above path


% find the row in the imaging data that corresponds to the specific cell type
switch spikedatafolder
    case 'MA22'
        cell_type = 'DM1';
        % Spike trials that correspond to the following odor order
        trial_loc = [79,83; 23,27; 68,72; 28,32; 33,37;...
            38,42; 43,47; 48,52; 53,57; 58,62; 51,55; 63,67; 74,78];
        % spike trials order
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But','H20'};
        legend_gap = 3;
    case '2025-11-18'
        cell_type = 'DM1';
        % Spike trials that correspond to the following odor order
        trial_loc = [1,5; 6,10; 11,15; 16,20; 21,25;...
            26,30; 31,35; 36,40; 41,45; 46,50; 51,55; 56,60; 61,65];
        % spike trials order
        spike_order = {'ACV', 'AA', 'EA', 'MA', 'ISOPENT', 'PentAc', 'EthBut', 'EHex', 'BENZ', '3OCT', '2BUT', '2,3-But'};
        legend_gap = 1;
     case '2025-12-02'
        cell_type = 'DM1';
        % Spike trials that correspond to the following odor order
        trial_loc = [39,43; 44,48; 49,53];
        % spike trials order
        spike_order = {'ACV', 'AA', 'EA'};
        legend_gap = 1;
     case '2025-12-18'
        cell_type = 'DM1';
        % Spike trials that correspond to the following odor order
        trial_loc = [2,6; 7,11; 12,16; 17,21];
        % spike trials order
        spike_order = {'ACV', 'AA', 'MA', 'ISOPENT'};
        legend_gap = 1;
     case 'MA21'
        cell_type = 'DM1';
        % Spike trials that correspond to the following odor order
        trial_loc = [2,6; 7,11; 12,16];
        % spike trials order
        spike_order = {'ACV', 'AA', 'MA'};
        legend_gap = 1;
end

% auto select the imaging row based on Kristyn's data organization for the given cell type
switch cell_type
    case 'DM1'
        imaging_row = 10;
    case 'VC1'
        imaging_row = 24;
    case 'VC2'
        imaging_row = 25;
    case 'VA4'
        imaging_row = 19;
end


% extract the specific data for this cell from the loaded imaging data
df_f = [imaging_data.Mean_PNresponse_duringodor{imaging_row,:}]; % pull the data into a matrix 

% Kristyn's odor order: 
dff_order = {'ACV', 'AA', 'BENZ', 'ISOPENT', 'EA', 'EthBut', 'PentAc', '2BUT', '3OCT', 'EHex', '2,3-But', 'MA'};
nOdors = length(dff_order); % number of odors in the imaging data

% reorder Maddie's spike data to match Kristyn's df_f data & find the average for each
avg_spike_rate = nan([1,nOdors]); % initialize empty variable for averages in the size of total number of odors
for i = 1:length(dff_order)
    match_odor = dff_order{i}; % find the target odor in Kristyns odor order to match
    idx = find(strcmp(match_odor, spike_order)); % find the location of the search odor in the spike rate odor order data
    if isempty(idx)
        continue
    end
    % pull out the trial numbers that match the odor of the correct order
    idx_trials = trial_loc(idx,1):trial_loc(idx,2); % find the trial numbers that correspond to this odor
    avg_spike_rate(i) = mean(spiketimeData(idx_trials),'omitnan'); % find the average spike rate for the odor trials
end


%% Save the processed data

save([saveDir, spikedatafolder '_avg_spike_rate.mat'], 'avg_spike_rate', 'dff_order', 'df_f', 'spikedatafolder');
disp(['Saved ' spikedatafolder ' data file'])

% TODO: do you want to save the name with the putative cell type as well
% for more easy selection moving forward? 

%% Plot the comparison figure:


cmap = turbo(nOdors); % select colors for each odor from a color map
 
fig = figure; set(fig, 'color', 'w'); hold on
for i = 1:nOdors
    scatter(df_f(i),avg_spike_rate(i),105,cmap(i,:),"filled")
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

% save the figure as a jpg in the spike timing data folder
saveas(fig, [fullpath title_str '.jpg']);
disp(['Saved figure ' title_str])

% ask to close the figure if desired
if strcmp(questdlg('Figure saved. Close figure?'),'Yes')
    close(fig)
end

