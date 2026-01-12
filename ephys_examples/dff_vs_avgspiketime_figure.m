
%% Load in data

% organize pathways to data folders
serverPath = '/Volumes/shared/Maddie/'; % pathway to maddie's folder on the server
basepath = [basepath 'spike rate data/']; % pathway to the folder that holds the processed and extracted spike rate data

% selecte the cell type that you are interested in?
cell_type = 'DM1';

foldercontents = dir([fullpath,'*' cell_type '.mat']); %search for only files that start with AllSpikeTimes in the selected folder
nCells = length(foldercontents); % total number of cells for this type

% TODO -- working here! 
for i = 1:nCells
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

end

% save([saveDir, spikedatafolder '_avg_spike_rate.mat'], 'avg_spike_rate', 'dff_order', 'df_f', 'spikedatafolder');


%% Find the average and error in the spike rate data for each odor

%% FIGURE: plot the df_f data vs spike rate for EACH cell (one subplot for each)

%% FIGURE: plot the avgs of the df_f vs spike rate for ALL cells





%% FIGURE: plot the df_f data compared to the spike rate data for each
% 
% cmap = turbo(nOdors); % select colors for each odor from a color map
% 
% fig = figure; set(fig, 'color', 'w')
% scatter(df_f,avg_spike_rate,75,cmap,"filled")
% % format the figure
% ylabel('mean spike rate (hz)')
% xlabel('mean df/F')
% set(gca, 'fontsize',20)
% legend(dff_order,'box', 'off')





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
 
