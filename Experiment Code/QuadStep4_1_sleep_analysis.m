
% Find 'sleep' points in data

% set a buffer circle
% create a 'running' count for occupancy of that circle 


X = data(1).data(1).data.x_loc;
Y = data(1).data(1).data.y_loc;

X()




%% grid spacing???

i = 1;
trial = 1;
vid = 1;
frame = 1;

% pull info for the first trial:
dataDate = data(i).T.Date{trial};
vid_name = data(i).T.ExperimentID{trial};
vidDir = [baseFolder dataDate '/' vid_name '_'];
videoPath = [vidDir num2str(vid) '.avi'];
movieInfo = VideoReader(videoPath); %read in video

% Set axis limits for the selected arena
x = data(i).data(trial).data.centre(1);
y = data(i).data(trial).data.centre(2);
r = data(i).data(trial).data.r;
xlimit = [x-(r+50),x+(r+50)];
ylimit = [y-(r+50),y+50+r];

% Plot image of video frame
fig = figure; set(fig,'pos',[-1030 279 772 1009],'color','k');
currentImg = rgb2gray(read(movieInfo,frame));
imshow(currentImg)
xlim(xlimit); ylim(ylimit);   

% find the 'auto bin' lines
nbins = 90;
xedge = linspace(xlimit(1),xlimit(2),nbins+1);
yedge = linspace(ylimit(1),ylimit(2),nbins+1);

% Plot the bin outline edges:
h_line(yedge,'yellow','-',0.25) .n 
v_line(xedge,'yellow','-',0.25)

%% 




% MOVEMENT:
xmin = arenaData(arena).centre(1)-r; % arena spacing divide by pixel size
xmax = arenaData(arena).centre(1)+r;
ymin = arenaData(arena).centre(2)-r;
ymax = arenaData(arena).centre(2)+r;
nbins = 100;
xedge = linspace(xmin,xmax,nbins+1);
yedge = linspace(ymin,ymax,nbins+1);
N = [];
for frame = 1:T.frame(end)
    X = arenaData(arena).x_loc(frame,:); X(isnan(X)) = [];
    Y = arenaData(arena).y_loc(frame,:); Y(isnan(Y)) = [];
    N(:,:,frame) = histcounts2(X,Y,xedge,yedge);
end
binDiff = diff(N,1,3);
binDiff(binDiff<0) = 0;
movement = squeeze(sum(sum(binDiff,1),2));
binSizeinMM = ((2*r)/nbins)/pix2mm;
movement = movement./binSizeinMM;
movement = movement*3; % convert to per second WITH 3fps videos
avg_movement = movement./T.flyCount(1:end-1,arena);% normalize for the number of flies actually tracked on the page:
% update data structures
occupancy.movement = avg_movement;
occupancy.baseMovement = movement;