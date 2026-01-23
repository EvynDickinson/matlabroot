
clear
clc

%load the data - already avg spike rate by dff for each cell type 
%plot each dff v spike rate hold on thing where it puts them all on the same plot - color coded by cell type not odor (legend)
%make the fit line for each cell type 
%receptor type - ors = DM1 and VA4, irs = DP1m and Vl2p
%second plot = make two fit lines, separated by receptor type

%% Load in spike rate data
%rootpath = '/Volumes/shared/Maddie/';
rootpath = 'S:\Maddie\';
basepath = [rootpath 'data from ephys rig/'];
fullpath = [basepath 'spike rate data/']; % folder with extracted spike rate data
foldercontents = dir([fullpath,'*avg_spike_rate.mat']); % find all the cells avg data
 
imaging_path =[rootpath 'Kristyn Data/Mean_PNresponse_duringodor.mat']; % hard coded from laptop
imaging_data = load(imaging_path); % load in the data from the above path

glomerulus_ID = {'DM1', 'VA4', 'DP1m', 'Vl2p'};
color_legend = {'k', 'k', 'r', 'r'};
marker_type = [true false true false]; %determines whether or not the circle is filled in the plot 
line_type = {'-', '--', '-', '--'}; %determines if line is solid or dashed 

%receptor types, 1 = OR, 2 = IR
receptor_typestr = {'OR', 'IR'};
spike_rate_avg = [];
dff_avg = [];
receptor_list = [];

for ii = 1:length(glomerulus_ID)
    PN_name = glomerulus_ID{ii};

    % pull cell type specific information
    switch PN_name 
        case 'DM1'
            imaging_row = 10;
            cell_idx = {'MA15', 'MA22', 'MA25'};
            receptor_type = 1;
        case 'VA4'
            imaging_row = 19;
            cell_idx = {'MA16', 'MA17', 'MA23'}; 
            receptor_type = 1;
        case 'Vl2p'
            imaging_row = 27;
            cell_idx = {'JJ221', 'JJ222'};
            receptor_type = 2;
        case 'DP1m'
            imaging_row = 16;
            cell_idx = {'JJ223'};
            receptor_type = 2;
    end

    receptor_list(ii) = receptor_type;

    % Find spike rate files that match the selected data type:
    spike_rate = []; 
    for jj = 1:length(cell_idx)
        fileName = [cell_idx{jj} '_avg_spike_rate.mat'];
        dummy = load([fullpath fileName]);
        spike_rate = [spike_rate; dummy.avg_spike_rate];
    end
    odor_names = dummy.dff_order;
    clear dummy

   %finding avg spike rate for each cell type per odor  
   spike_rate_avg = [spike_rate_avg; mean(spike_rate,1,'omitnan')]; 
   %finding the avg df/f for each cell type per odor 
   df_f = [imaging_data.Mean_PNresponse_duringodor{imaging_row,:}]; % pull the data into a matrix 
   dff_avg = [dff_avg; df_f];

end 

%% find the fit for each cell type 

%plotting parameters
uni_x = [-0.5, 5];
uni_y = [-20,150];

legend_str = {};
%make figure 
fig = figure; set(fig, 'color', 'w'); hold on
    for ii = 1:length(glomerulus_ID)
        x = dff_avg(ii,:);
        y = spike_rate_avg(ii,:);
        %scatter plotting 
        if marker_type(ii)
            scatter(x,y,35,color_legend{ii},"filled")
        else 
            scatter(x,y,35,color_legend{ii})
        end 
        legend_str{end+1} = glomerulus_ID{ii};
        % make line of best fit for each PN type 
        mdl = fitlm(x, y, 'Intercept', false);
        ypred = predict(mdl,uni_x');
        plot(uni_x,ypred,'color',color_legend{ii},'LineStyle',line_type{ii})
        legend_str{end+1} = [glomerulus_ID{ii} ' fit'];
    end 

% format the figure
title_str =(['df/F vs spike rate']);
title(title_str)
ylabel('mean spike rate (hz)')
xlabel('mean df/F')
set(gca, 'fontsize',18,'fontname', 'calibri')
ylim([uni_y])

%make legend
legend(legend_str,"Location","northwest")



%% Collapse across ORs and IRs for avg lines


%plotting parameters
color_legend = {'k', 'r'};

legend_str = {};
%make figure 
fig = figure; set(fig, 'color', 'w'); hold on
    
    for ii = 1:length(receptor_typestr)
        loc = receptor_list == ii;
        x = dff_avg(loc,:);
        x_avg = mean(x,1,'omitnan');
        y = spike_rate_avg(loc,:);
        y_avg = mean(y,1,'omitnan');
        %scatter plotting 
        scatter(x(:),y(:),35,color_legend{ii},"filled")
        legend_str{end+1} = receptor_typestr{ii};
        % make line of best fit for each PN type 
        mdl = fitlm(x_avg, y_avg, 'Intercept', false);
        ypred = predict(mdl,uni_x');
        plot(uni_x,ypred,'color',color_legend{ii})
        legend_str{end+1} = [receptor_typestr{ii} ' fit'];
    end 

% format the figure
title_str =(['df/F vs spike rate']);
title(title_str)
ylabel('mean spike rate (hz)')
xlabel('mean df/F')
set(gca, 'fontsize',18,'fontname', 'calibri')
ylim(uni_y)

%make legend
legend(legend_str,"Location","northwest")


% naming scheme:  dff_vs_spikerate_byreceptortype_1_23_26











%%
clear
clc

% collapse data across the ORs and IRs

%% Load in spike rate data
%rootpath = '/Volumes/shared/Maddie/';
rootpath = 'S:\Maddie\';
basepath = [rootpath 'data from ephys rig/'];
fullpath = [basepath 'spike rate data/']; % folder with extracted spike rate data
foldercontents = dir([fullpath,'*avg_spike_rate.mat']); % find all the cells avg data
 
imaging_path =[rootpath 'Kristyn Data/Mean_PNresponse_duringodor.mat']; % hard coded from laptop
imaging_data = load(imaging_path); % load in the data from the above path

glomerulus_ID = {'DM1', 'VA4', 'DP1m', 'Vl2p'};
color_legend = {'k', 'k', 'r', 'r'};
marker_type = [true false true false]; %determines whether or not the circle is filled in the plot 
line_type = {'-', '--', '-', '--'}; %determines if line is solid or dashed 

spike_rate_avg = [];
dff_avg = [];

for ii = 1:length(glomerulus_ID)
    PN_name = glomerulus_ID{ii};

    % pull cell type specific information
    switch PN_name 
        case 'DM1'
            imaging_row = 10;
            cell_idx = {'MA15', 'MA22', 'MA25'};
        case 'VA4'
            imaging_row = 19;
            cell_idx = {'MA16', 'MA17', 'MA23'}; 
        case 'Vl2p'
            imaging_row = 27;
            cell_idx = {'JJ221', 'JJ222'};
        case 'DP1m'
            imaging_row = 16;
            cell_idx = {'JJ223'};
    end

    % Find spike rate files that match the selected data type:
    spike_rate = []; 
    for jj = 1:length(cell_idx)
        fileName = [cell_idx{jj} '_avg_spike_rate.mat'];
        dummy = load([fullpath fileName]);
        spike_rate = [spike_rate; dummy.avg_spike_rate];
    end
    odor_names = dummy.dff_order;
    clear dummy

   %finding avg spike rate for each cell type per odor  
   spike_rate_avg = [spike_rate_avg; mean(spike_rate,1,'omitnan')]; 
   %finding the avg df/f for each cell type per odor 
   df_f = [imaging_data.Mean_PNresponse_duringodor{imaging_row,:}]; % pull the data into a matrix 
   dff_avg = [dff_avg; df_f];

end 

%% find the fit for each cell type 

%plotting parameters
uni_x = [-0.5, 5];
uni_y = [-20,120];

legend_str = {};
%make figure 
fig = figure; set(fig, 'color', 'w'); hold on
    for ii = 1:length(glomerulus_ID)
        x = dff_avg(ii,:);
        y = spike_rate_avg(ii,:);
        %scatter plotting 
        if marker_type(ii)
            scatter(x,y,35,color_legend{ii},"filled")
        else 
            scatter(x,y,35,color_legend{ii})
        end 
        legend_str{end+1} = glomerulus_ID{ii};
        % make line of best fit for each PN type 
        mdl = fitlm(x, y, 'Intercept', false);
        ypred = predict(mdl,uni_x');
        plot(uni_x,ypred,'color',color_legend{ii},'LineStyle',line_type{ii})
        legend_str{end+1} = [glomerulus_ID{ii} ' fit'];
    end 

% format the figure
title_str =(['df/F vs spike rate']);
title(title_str)
ylabel('mean spike rate (hz)')
xlabel('mean df/F')
set(gca, 'fontsize',18,'fontname', 'calibri')

%make legend
legend(legend_str,"Location","northwest")



