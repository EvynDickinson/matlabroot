
%% Select and load video
clearvars
% select video file:
cd 'C:\Users\evynd\Desktop\HeadlessFliesDLC\tracking\off_ball-flies\training_videos\videos-labeled\';

[baseFileName, folderName, FilterIndex] = uigetfile('*.mp4');
movieFullFileName = fullfile(folderName, baseFileName);
v = VideoReader(movieFullFileName);

%% Display Frames (six frames w/ pixel boxes)
fig = getfig; idx = 1;
frameOpt = randi([1,600],1,6);
% select the frame to display
for frameNum = frameOpt
    subplot(2,3,idx)
%     frameNum = 300;
    frame = read(v,frameNum);

    % display the frame
    hold all; 
    imshow(frame)

    invisibleboundary_X = 120;
    invisibleboundary_Y = 40;

    % draw a pixel box of a given size
    imrect(gca, [0, 0, invisibleboundary_X, invisibleboundary_Y]);
    xlabel(['Frame: ' num2str(frameNum)])
    idx = idx+1;
end

%% Reset settings

% reset the current folder to matlabroot
cd 'C:\matlabroot'
close all

%% Look at placement of the original points (all overlaid) from VIA output (not converted from Pierres)
clearvars

cd 'C:\Users\evynd\Desktop\HeadlessFliesDLC\tracking\off_ball-flies\labeled-data\';

[baseFileName, folderName, ~] = uigetfile('*.csv');
csvFile = fullfile(folderName, baseFileName);

metaData = readtable(csvFile);

A = metaData(:,6); %xy coord column
a = table2cell(A); %convert to a cell

B = metaData(:,5); %xy coord column
b = table2cell(B); %convert to a cell


tic
for ii = 1:length(a)
    % find the coordinates
    splitdata = strsplit(a{ii}, {'"', '{', '}', ',', ':'});

    output(ii,1) = str2double(splitdata{5});
    output(ii,2) = str2double(splitdata{7});
    % find the point index
    output(ii,3) = b{ii};
end
toc 

% %cam A
% xmax = 528;
% ymax = 450;
% % cam B
% xmax = 464;
% ymax = 355;
% % cam C
% xmax = 832;
% ymax = 544;
% % cam D
xmax = 488;
ymax = 521;
% % cam E
% xmax = 560;
% ymax = 425;
% cam F
% xmax = 464;
% ymax = 344;

clearvars -except metaData output csvFile ymax xmax

% check for missing data points:
% figure; plot(output(:,3)) 
% 2277/30

close all
% split the data based on the label location:

%% Plot all the points from all frames and the boundary lines
% scatter plot the points:

% boundaries
invisibleboundary_X = 70;
invisibleboundary_Y = 50;


imgStart = 1;
imgEnd = 50;

ROI = [(imgStart-1)*30+1:imgEnd*30];

figure; %make the x and y lims the image crop size - taken from config.yaml:
xlim([0, xmax])
ylim([0, ymax])
hold all
scatter(output(ROI,1), output(ROI,2), 30, 'b', 'filled')
vline(invisibleboundary_X, 'k')
hline(invisibleboundary_Y, 'k')


%% Look at 9 figs at a time -- zoom on the region

startNum = 92 ;

fig = getfig('',1); idx = 1;
for ii = startNum:startNum+8
    subplot(3,3,idx)
    imgStart = ii;
    imgEnd = ii;

    ROI = [(imgStart-1)*30+1:imgEnd*30];

    %make the x and y lims the image crop size - taken from config.yaml:
    xlim([0, invisibleboundary_X+50])
    ylim([0, invisibleboundary_Y+50])
    hold all
    scatter(output(ROI,1), output(ROI,2))
    vline(invisibleboundary_X, 'k')
    hline(invisibleboundary_Y, 'k')
    idx = idx+1;
    set(gca, 'YTickLabel', [], 'XTickLabel', [])
    xlabel(['# ' num2str(ii)]) 
end


%%  LOOK AT DATA POINTS FROM PIERRES DLC CONVERT VIA SCRIPT

clearvars

cd 'C:\Users\evynd\Desktop\HeadlessFliesDLC\tracking\off_ball-flies\labeled-data\';

[baseFileName, folderName, FilterIndex] = uigetfile('*.csv');
csvFile = fullfile(folderName, baseFileName);

metaData = readtable(csvFile);

A = metaData(:,6); %xy coord column
a = table2cell(A); %convert to a cell


tic
for iframe = 1:length(a)
    % find the coordinates
    splitdata = strsplit(a{iframe}, {'"', '{', '}', ',', ':'}); %separate the data
    for ii = 2:length(splitdata) %pull the point locations
        output(iframe,ii-1) = str2double(splitdata{ii});
    end
end

% CamName = 'Cam A';
% xmax = 528;
% ymax = 450;

% CamName = 'Cam B';
% xmax = 464;
% ymax = 355;

% CamName = 'Cam C';
% xmax = 832;
% ymax = 544;

CamName = 'Cam D';
xmax = 488;
ymax = 521;

% CamName = 'Cam E';
% xmax = 560;
% ymax = 425;

% CamName = 'Cam F';
% xmax = 464;
% ymax = 344;

clearvars -except metaData output csvFile ymax xmax CamName

% check for missing data points:
% figure; plot(output(:,3)) 
% 2277/30

close all
% split the data based on the label location:


% visualize the data points: 
invisibleboundary_X = 70;
invisibleboundary_Y = 50;

sz = 30;
fig = getfig('', 1);
xlim([0,xmax])
ylim([0,ymax])
hold all

% plot outlier square:
plot([0,invisibleboundary_X],[invisibleboundary_Y,invisibleboundary_Y],...
    'linewidth', 2, 'color', 'r')
plot([invisibleboundary_X,invisibleboundary_X],[0,invisibleboundary_Y],...
    'linewidth', 2, 'color', 'r')


for frameNum = 90:100
    for joints = 1:size(output,2)/2
%         x = output(frameNum,joints);
%         y = output(frameNum,joints+1);
        x = flipud(output(frameNum,joints));
        y = (output(frameNum,joints+1));
        
        scatter(x,y, sz, 'g', 'filled')

    end
end

title(CamName)
toc 


% save_figure(fig, CamName)


%%  OVERLAY BOTH DLC AND VIA OUTPUTS TO COMPARE
% check for missing data points:
% figure; plot(VIAoutput(:,3)) 
% ceil(839/30)

% missing frames: 
% missing_frame = ceil(find((diff(VIAoutput(:,3))<1 & diff(VIAoutput(:,3))>-29)==true)/30);
% if sum(missing_frame>0)
%     fprintf('\n Missing label(s) in frame:')
%     disp(missing_frame)
% end

%% VIA and DLC_ConvertVIA overlay
clearvars
base = 'C:\Users\evynd\Desktop\HeadlessFliesDLC\tracking\off_ball-flies\labeled-data\';

CamChoice = {'A','B', 'C', 'D', 'E', 'F'};
choice = listdlg('ListString', CamChoice, 'PromptString','Temperature?',...
                'SelectionMode', 'Single', 'ListSize', [80 150]);
Cam = CamChoice{choice}; 
clear choice

% FolderName = [base 'offball_flies_2020-06-18_11-27_random_' Cam '\']; %this needs to be updated to match the folders
FolderName = [base 'offball_flies_2020-07-13_19-59_bad_' Cam '\']; %this needs to be updated to match the folders
nframes = 50;

% Load the OG VIA file data
VIAdata = readtable([FolderName, 'via_export_csv.csv']);

A = VIAdata(:,6); %xy coord column
a = table2cell(A); %convert to a cell
B = VIAdata(:,5); %xy coord column
b = table2cell(B); %convert to a cell

for ii = 1:length(a)
    % find the coordinates
    splitdata = strsplit(a{ii}, {'"', '{', '}', ',', ':'});

    VIAoutput(ii,1) = str2double(splitdata{5});
    VIAoutput(ii,2) = str2double(splitdata{7});
    % find the point index
    VIAoutput(ii,3) = b{ii};
end; clear a A b B

missing_frame = ceil(find((diff(VIAoutput(:,3))<1 & diff(VIAoutput(:,3))>-29)==true)/30);
if sum(missing_frame>0)
    fprintf('\n Missing label(s) in frame:')
    disp(missing_frame)
    return
end

% Load the DLC label data
DLCdata = readtable([FolderName, 'CollectedData_TuthillLab.csv']);

A = DLCdata(:,6); %xy coord column
a = table2cell(A); %convert to a cell

for iframe = 1:length(a)
    % find the coordinates
%     splitdata = strsplit(a{iframe}, {'"', '{', '}', ',', ':'}); %separate the data
    splitdata = strsplit(a{iframe}, ','); %separate the data
    for ii = 2:length(splitdata) %pull the point locations
        DLCoutput(iframe,ii-1) = str2double(splitdata{ii});
    end
end
loc = nansum(DLCoutput)==0;
if sum(loc>0)
    DLCoutput(:,loc) = [];
    fprintf('\n Deleted extra data column\n')
end

% Visualize the data: 
sz = 30;

fig = getfig('',1);
hold all

imgStart = 1;
imgEnd = nframes;
imageShow = 1;

% load a demo frame into the space
fileList = dir(fullfile(FolderName, '*.png'));
imshow([FolderName, fileList(imageShow).name]); hold all


% invisible boundary lines:
invisibleboundary_X = 70;
invisibleboundary_Y = 50;
plot([0,invisibleboundary_X],[invisibleboundary_Y,invisibleboundary_Y],...
    'linewidth', 1, 'color', 'r')
plot([invisibleboundary_X,invisibleboundary_X],[0,invisibleboundary_Y],...
    'linewidth', 1, 'color', 'r')


% Overlay the VIA labels:
ROI = (imgStart-1)*30+1:imgEnd*30;
scatter(VIAoutput(ROI,1), VIAoutput(ROI,2), sz, 'b', 'filled')


% Overlay the DLC_ConvertVIA labels: 
for ii = imgStart:imgEnd
    raw = DLCoutput(ii,:);
    x_idx = (1:2:length(raw));
    x = raw(x_idx);
    y_idx = (2:2:length(raw));
    y = raw(y_idx);
    scatter(x,y, sz+20, 'y')
    raw = [];
end


save_figure(fig, ['C:\matlabroot\DLC\Demo Images\Cam ' Cam ' bad VIA points']);
% close all



%%  which points have a nan in the DLC but not the OG via?
%restructure the VIA points into DLC format:
idx = 1;
for iframe = 1:100
    ROI = idx:idx+29;
    
    x_loc = 1:2:60;
    y_loc = 2:2:60;
    
    VIA_rearranged(iframe,x_loc) = VIAoutput(ROI,1);
    VIA_rearranged(iframe,y_loc) = VIAoutput(ROI,2);
    
    idx = ROI(end)+1;
end

sz = 30;

fig = getfig('',1);
hold all

imgStart = 1;
imgEnd = 1;
imageShow = 1;

% load a demo frame into the space
fileList = dir(fullfile(FolderName, '*.png'));
imshow([FolderName, fileList(imageShow).name]); hold all


% invisible boundary lines:
invisibleboundary_X = 70;
invisibleboundary_Y = 50;
plot([0,invisibleboundary_X],[invisibleboundary_Y,invisibleboundary_Y],...
    'linewidth', 1, 'color', 'r')
plot([invisibleboundary_X,invisibleboundary_X],[0,invisibleboundary_Y],...
    'linewidth', 1, 'color', 'r')


% Overlay the DLC_ConvertVIA labels: 
for ii = imgStart:imgEnd
    raw = VIA_rearranged(ii,:);
    x_idx = (1:2:length(raw));
    x = raw(x_idx);
    y_idx = (2:2:length(raw));
    y = raw(y_idx);
    scatter(x,y, sz+40, 'm', 'filled')
    raw = [];
end

% Overlay the VIA labels:
ROI = (imgStart-1)*30+1:imgEnd*30;
scatter(VIAoutput(ROI,1), VIAoutput(ROI,2), sz, 'b', 'filled')



% Overlay the DLC_ConvertVIA labels: WHY DOESNT IT WORK??
for ii = imgStart:imgEnd
    adj = 44;
    raw = DLCoutput(ii,:);
    x_idx = (1:2:adj);
%     x_idx = (1:2:length(raw));
    x = raw(x_idx);
    y_idx = (2:2:adj);
%     y_idx = (2:2:length(raw));
    y = raw(y_idx);
    scatter(x,y, sz-10, 'y', 'filled')
    raw = [];
end





