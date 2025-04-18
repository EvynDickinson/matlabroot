clear
baseFolder = getCloudPath;

% Load waxed fly responses
data_dir = [baseFolder 'Imaging/Waxed antenna 2P imaging Kristyn/'];
load([data_dir 'meanDFF_ACV_heavywax_antennaeandpalps_new.mat']) % load avg data

% load([data_dir 'meanDFF_ACV_heavy_wax_antennaeandpalps.mat']) % load avg data
load([data_dir 'glomeruli_labels.mat']) % load labels for the glomeruli

% Load control fly responses
control = load([baseFolder 'Imaging/Intact 2P imaging ACV Kristyn/meanDFF_ACV.mat']);
intact = control.meanDFF;

n = length(glom_ID);
hz = 4.22; %check this with Kristyn (was it really ~12 sec?)
odor_on = 4;
odor_dur = 2;

%% SKIP :  load provided data and pull the glom lables
% 
% % get handles for each subplot on the graph
% subplot_handles = findobj(gcf, 'type', 'axes');
% 
% n = length(subplot_handles);
% 
% % pull names into cell format
% glom_ID = cell(size(subplot_handles));
% 
% for i = 1:n
%     glom_ID{i} = get(subplot_handles(i), 'Title').String;
% end
% 
% glom_ID = flip(glom_ID);
% 
% save([data_dir 'glomeruli_labels.mat'],'glom_ID')

%% FIGURE: plot the mean df/f for each glomeruli
%figure parameters
blkbgd = false;
fig_type = '-pdf';
[foreColor,backColor] = formattingColors(blkbgd); %get background colors
kolor = Color('red');
LW = 1.5;
col = 7; 
row = 6;

%plot figure
new_fig = getfig('',1);
xlimits = []; skipIdx = [];
for i = 1:n
    subplot(row,col,i)
    hold on
    
    % control flies
    y1 = intact{1,i};
    len = size(y1,1);
    x = linspace(0,len/hz,len);
    plot(x,y1,'Color', foreColor,'linewidth',LW)


    % waxed flies
    y = meanDFF{1,i};
    len = length(y);
    x = linspace(0,len/hz,len);
    plot(x,y,'Color', kolor,'linewidth',LW)

    %

    % add glomerulus identity
    title(glom_ID{i},'FontSize',12,'Color',foreColor)
    % get info for normalizing
    xlimits(i,:) = xlim;
    ylimits(i,:) = ylim;
    if all(all(isnan(y))) && all(all(isnan(y1)))
        skipIdx(i) = true;
    else
        skipIdx(i) = false;
    end
end

X_norm = [min(min(xlimits)), max(max(xlimits))];
Y_norm = [min(min(ylimits)), max(max(ylimits))];


formatFig(new_fig,blkbgd,[row,col]);
xRect = [odor_on,odor_on, odor_on+odor_dur, odor_on+odor_dur];
yRect = [Y_norm(1), Y_norm(2), Y_norm(2), Y_norm(1)];
for i = 1:n
    subplot(row,col,i)

    xlim(X_norm)
    xlim([1,X_norm(2)]) % cut out the start with the photobleach response
    ylim(Y_norm)
    set(gca,'xcolor', backColor,'xtick', [],'ycolor',backColor,'ytick',[]);
    if ~skipIdx(i)
        % Draw the odor stim rectangle
        hold on
        % Create the rectangle as a semi-transparent patch
        patch(xRect, yRect, foreColor, 'FaceAlpha', 0.3,'EdgeColor','none');
    end
    % color the first column 
    if ismember(i,1:col:n)
        set(gca,'ycolor',foreColor,'ytick',Y_norm);
        ylabel('df/F')
    end
end

% save_figure(new_fig,[data_dir 'HeavyWax_ACV_meanDFF'],fig_type);


save_figure(new_fig,[data_dir 'HeavyWax_ACV_meanDFF ALL'],fig_type);

%%

 X_norm = [0   14.2180];
 Y_norm = [-1, 10];


% Example plot
x = linspace(0, 2*pi, 100);
y = sin(x);
plot(x, y);
hold on; % Keep the plot for adding the rectangle

% Coordinates for the rectangle (specify lower left corner, width, and height)
 rectangle('Position', [odor_on, Y_norm(1), odor_dur, diff(Y_norm)], ...
              'FaceColor', foreColor, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

xRect = [odor_on,odor_on, odor_on+odor_dur, odor_on+odor_dur];
yRect = [Y_norm(1), Y_norm(2), Y_norm(2), Y_norm(1)];

% Create the rectangle as a semi-transparent patch
patch(xRect, yRect, foreColor, 'FaceAlpha', 0.3);

hold off; % Release the plot hold


















