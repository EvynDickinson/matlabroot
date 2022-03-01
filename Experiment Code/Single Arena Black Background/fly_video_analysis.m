clear

%% Select folder and video to load

%get base folder pathway
switch getenv('COMPUTERNAME')
    case 'DENALI'
        baseFolder = 'E:\My Drive\Jeanne Lab\DATA\';
    case 'TOGIAK'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
end

%select folder date
list_dirs = dir(baseFolder);
for i = 3:length(list_dirs)
    folderNames{i-2} = list_dirs(i).name;
end
folderNames = ['Today', folderNames];

%Select the flies to analyze
indx = listdlg('ListString', folderNames, 'SelectionMode', 'Single');
if strcmpi(folderNames{indx}, 'Today')==true
    dir_sel = strrep(datestr(datetime,'mm-dd-yyyy'),'-','.');
else
    dir_sel = folderNames{indx};
end
    
folder = fullfile(baseFolder, dir_sel);

%Get list of available videos
list_dirs = dir([folder, '\*.avi']); %only videos
for i = 1:length(list_dirs)
    vidNames{i} = list_dirs(i).name;
end
indx = listdlg('ListString', vidNames, 'SelectionMode', 'Single');
vidPath = fullfile(folder, vidNames{indx});

vid = 1;

%% Load video -- save cumulative pixel counts

%Get list of available videos
list_dirs = dir([folder, '\*.avi']); %only videos
for i = 1:length(list_dirs)
    vidNames{i} = list_dirs(i).name;
end
indx = listdlg('ListString', vidNames, 'SelectionMode', 'Multiple');

data = [];
for n = 1:length(indx)
    
    vid = indx(n); %select video file from list
    vidPath = fullfile(folder, vidNames{vid}); %full video path
    movieInfo = VideoReader(vidPath); %read in video
    
    %extract parameters from video
    nframes = movieInfo.Duration * movieInfo.FrameRate;
    height = movieInfo.Height;
    width = movieInfo.Width;
    tot_occ = zeros(height, width); %setup blank occupation
    n_tot = nframes;
        
    h = waitbar(0,...
        ['reading in frames from vid ' num2str(n) '/' num2str(length(indx))]);
    % Read and process the videos frame-by-frame
    for i = 1:n_tot
        % Read in current frame from raw video
        frame = (rgb2gray(read(movieInfo,i))); %convert to greyscale
        BWframe = imbinarize(im2double(frame),... % threshold black and white
                  'adaptive','ForegroundPolarity','bright','Sensitivity',0.3);
        [currProb, tot_occ] = (getOccProb(tot_occ,BWframe,i)); % spatial probabilty 
        waitbar(i/n_tot,h)
    end
    close(h)

    % save specific vid data
    data.tot(:,:,n) = tot_occ;
    data.prob(:,:,n) = currProb;

end
fprintf('\n All loaded! \n')  

 
%% save the data! 
save([folder, '\occupational prob data'])



%% Basic visualization of the flies density over the total video duration

% get the total across all loaded videos
tot_occupancy = sum(data.tot,3);
tot_probabilty = sum(data.prob,3)/length(indx);
  

% --- log spatial distribution ---
fig = getfig; set(fig, 'pos', [50 604 1141 346], 'color', 'k'); 
subplot(1,2,1)
    h = heatmap(tot_probabilty); 
    colormap hot;
    set(h,'gridvisible', 'off')
    ax = gca;
    ax.XDisplayLabels = nan(size(ax.XDisplayData));
    ax.YDisplayLabels = nan(size(ax.YDisplayData));
    ax.ColorScaling = 'log';
    set(ax,'FontColor', 'w', 'FontName', 'Arial');
    title('Spatial occupation probability (log scale)');
subplot(1,2,2)
    h = heatmap(tot_probabilty); 
    colormap hot;
    set(h,'gridvisible', 'off')
    ax = gca;
    ax.XDisplayLabels = nan(size(ax.XDisplayData));
    ax.YDisplayLabels = nan(size(ax.YDisplayData));
    set(ax,'FontColor', 'w', 'FontName', 'Arial');
    title('Spatial occupation probability');
    
    
    save_figure(fig, [folder, '\occupation probability'], '-png')


    

%% Generate a circular mask for the bowl: % do this early on?
figure;
imshow(frame)
roi = drawcircle;
centre = roi.Center;
radius = roi.Radius;

% Define coordinates and radius
x1 = centre(1);
y1 = centre(2);

% Generate grid with binary mask representing the circle. Credit to Jonas for original code.
[xx,yy] = ndgrid((1:height)-y1,(1:width)-x1);
mask = (xx.^2 + yy.^2>radius^2);


fig = figure; set(fig, 'color', 'k')
imshow(Im)   
save_figure(fig, [folder, '\arena masked'], '-png')
 
fig = figure; set(fig, 'color', 'k')
imshow(frame)   
save_figure(fig, [folder, '\arena NOTmasked'], '-png')


% Mask the probability image:
Img = tot_probabilty;
Img(mask) = 0;


% Binned pixels
fig = getfig; set(fig, 'pos', [50 50 1101 900], 'color', 'k'); 
    h = heatmap(imresize(Img,[20 20])); 
    colormap hot;
    set(h,'gridvisible', 'off')
    ax = gca;
    ax.XDisplayLabels = nan(size(ax.XDisplayData));
    ax.YDisplayLabels = nan(size(ax.YDisplayData));
%     ax.ColorScaling = 'log';
    set(ax,'FontColor', 'w', 'FontName', 'Arial', 'FontSize', 20);
    title('Spatial occupation probability (log scale)');
    
save_figure(fig, [folder, '\binned occupation probability'], '-png')


%% Load selected video and write processed video
% movieFullFileName = fullfile(folder, 'test_vid_1.avi');
% Check to see that it exists.
if ~exist(vidPath, 'file')
  strErrorMessage = sprintf('File not found:\n%s\nYou can choose a new one, or cancel', movieFullFileName);
  response = questdlg(strErrorMessage, 'File not found', 'OK - choose a new movie.', 'Cancel', 'OK - choose a new movie.');
  if strcmpi(response, 'OK - choose a new movie.')
    [baseFileName, folderName, FilterIndex] = uigetfile('*.avi');
    if ~isequal(baseFileName, 0)
      movieFullFileName = fullfile(folderName, baseFileName);
    else
      return;
    end
  else
    return;
  end
end

movieInfo = VideoReader(vidPath);
nframes = movieInfo.Duration * movieInfo.FrameRate;
height = movieInfo.Height;
width = movieInfo.Width;

% pull in specs for labeling images
txt = OccParams(vidPath, 5, 50);

%New video information
vid_name = [vidPath(1:end-4) ' feature extraction'];
FrameRate = 20; %playback rate
% Set up recording video parameters
set(0,'DefaultFigureVisible','off');
fig = getfig; 
    hold on
    v = VideoWriter(vid_name, 'Motion JPEG AVI');
    v.Quality = 75;
    v.FrameRate = FrameRate;
    open(v);
    set(fig, 'Color', 'k')   
    % add scale bar
    imshow(zeros(txt.frameSize))
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset') 
 

%load all the data for the video & early process images:
n_tot = nframes;
tot_occ = (zeros(height, width));

tic
h = waitbar(0,'reading & writing frames');
for i = 1:n_tot 
    % Read in current frame from raw video
    frame = (rgb2gray(read(movieInfo,i)));
    BWframe = imbinarize(im2double(frame),...
              'adaptive','ForegroundPolarity','bright','Sensitivity',0.3);
    [currProb, tot_occ] = (getOccProb(tot_occ,BWframe,i)); % spatial probabilty
    currProb(currProb==1) = 0; %filter out light particles
    % Build an image frame
    newFrame = OccBuildFrame(im2double(frame), BWframe, currProb, txt);
    % Collect image info for movie file
    set(fig, 'Color', 'k')
    imshow(newFrame)
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')
    
    waitbar(i/n_tot,h)
end
close(h)
close(v)
toc 

fprintf('\n Video Saved! \n')  
  
set(0,'DefaultFigureVisible','on');
    

%% Grid population density test

% 1) select a video
% 2) load frame
% 3) 

Z = tot_occ;
gridProb  = imresize(tot_occ, [10 10]);


figure; h = heatmap(gridProb); 
set(h,'gridvisible', 'off')
ax = gca;
% ax.ColorScaling = 'log';


figure;
imshow(gridProb)

figure;
imshow(tot_occ)

range(range(tot_occ))

X = 1:height;
Y = 1:width;


figure;
surf(Y,X,Z)


figure; imagesc(tot_occ)




%% Build video frames with side by side of raw & processed images

% space buffers
MT = zeros(height, width);
line_size = 5;
vert = ones(height,line_size);
hori = ones(line_size,(width*3)+(line_size*4));

% build a colorbar
barWidth = 50;
horimini = ones(line_size,barWidth+(line_size*2));
buffZone = zeros(height+(2*line_size), barWidth);
clims = [1,0];
cmap = repmat(linspace(clims(1),clims(2), height)',1,barWidth);
cbar = [horimini; vert,cmap,vert; horimini];

% set up text locations:
offsets = [2,4]; %top&bottom, right&left
blankFrame = ([vert, MT, vert, MT,vert, MT, vert]);
blankFrame = [hori;blankFrame; hori];    
blankFrame = [blankFrame,buffZone,cbar];
blankFrame = [zeros(barWidth*offsets(1),size(blankFrame,2)); blankFrame]; % top buffer
blankFrame = [zeros(size(blankFrame,1), barWidth),blankFrame];%right side buffer
blankFrame = [blankFrame, zeros(size(blankFrame,1), offsets(2)*barWidth)]; % left side buffer
blankFrame = [blankFrame; zeros(barWidth,size(blankFrame,2))]; % bottom buffer
frameSize = size(blankFrame);

endloc = frameSize(2)-barWidth*(offsets(2)-1);
pos = [barWidth, barWidth+line_size+round(width/2);...
        barWidth, barWidth+2*line_size+round(width*1.5);...
        barWidth, barWidth+3*line_size+round(width*2.5);...
        2*barWidth, endloc;...
        frameSize(1)-barWidth, endloc];
loc = [pos(:,2), pos(:,1)];
boxcolors = repmat([0,0,0],[5,1]);
str = {'raw', 'binary', 'occupation probabilty', '1', '0'};


% generate each frame
for n = 1:nframes
    newFrame = ([vert, mov(:,:,n), vert, BWvid(:,:,n),vert, probOccupation(:,:,n), vert]);
    newFrame = [hori;newFrame; hori];    
    newFrame = [newFrame,buffZone,cbar];
    newFrame = [zeros(barWidth*offsets(1),size(newFrame,2)); newFrame]; % top buffer
    newFrame = [zeros(size(newFrame,1), barWidth),newFrame];%right side buffer
    newFrame = [newFrame, zeros(size(newFrame,1), offsets(2)*barWidth)]; % left side buffer
    newFrame = [newFrame; zeros(barWidth,size(newFrame,2))]; % bottom buffer

    % insert title and text...
    newFrame = (insertText(newFrame, loc, str,'FontSize',40,'BoxColor',...
        boxcolors,'BoxOpacity',1,'TextColor','white', 'anchorpoint', 'center'));

    title_loc = [size(vidFrame,2)-barWidth, barWidth+round(height/2)];
    newFrame = (insertText(newFrame, title_loc, 'Prob.','FontSize',40,'BoxColor',...
        'black','BoxOpacity',1,'TextColor','white', 'anchorpoint', 'center'));
    movie(n).frame = newFrame;
end

% 
% fig = figure; set(fig, 'color', 'k');
% imshow(newFrame, 'Border', 'tight')

%% Write demo vid
set(0,'DefaultFigureVisible','off');

vid_name = [vidPath(1:end-4) ' feature extraction 3'];
FrameRate = 20;

h = waitbar(0,'Writing Video');
fig = getfig; 
hold on
v = VideoWriter(vid_name, 'Motion JPEG AVI');
v.Quality = 90;
v.FrameRate = FrameRate;
open(v);
set(fig, 'Color', 'k')   
    % add scale bar
    imshow(zeros(frameSize))
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset') 
    
for n = 1:nframes
    set(fig, 'Color', 'k')
    
    % add scale bar
    imshow(movie(n).frame)
    
    % Collect image info for movie file
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')
    waitbar(n/nframes,h)
end
close(v)
close(h)
close all
fprintf('\n Video Saved! \n')  

   
set(0,'DefaultFigureVisible','on');
    
%% Pull the cumulative probabilty for each frame

%make a demo vid with raw vid, binary vid, heatmap


% Place the three images side by side
n = 198; %frame number
fig = getfig;
set(fig, 'Color', 'k')
% extract the three frames
rawframe = mov(:,:,n);
editframe = BWvid(:,:,n);
occProb = sum(BWvid(:,:,1:n),3)/n;

vidFrame = ([vert, rawframe, vert, editframe,vert, occProb, vert]);
    vidFrame = [hori;vidFrame; hori];

% add scale bar
imshow(vidFrame)
    cb = colorbar;
        %label and aesthetics
        cb.FontName = 'Arial';
        cb.FontSize = 10;
        cb.Color = 'w';
        cb.Label.String = 'Occupation probabilty/pixel';
    % TITLES:
    dim = [0.1860 0.93 0.0343 0.0535];
    a = annotation('textbox',dim,'String','Raw','FitBoxToText','on');
    set(a, 'color', 'w', 'EdgeColor', 'k')
    dim2 = [0.4568 0.93 0.0343 0.0535];
    a = annotation('textbox',dim2,'String','Binary','FitBoxToText','on');
    set(a, 'color', 'w', 'EdgeColor', 'k')
    dim3 = [0.7451 0.93 0.0343 0.0535];
    a = annotation('textbox',dim3,'String','Occupation','FitBoxToText','on');
    set(a, 'color', 'w', 'EdgeColor', 'k')

%% Write demo vid
set(0,'DefaultFigureVisible','off');

vid_name = [vidPath(1:end-4) ' feature extraction 2'];
FrameRate = 20;
line_size = 5;
vert = ones(height,line_size);
hori = ones(line_size,(width*3)+(line_size*4));

fig = getfig; 
hold on
v = VideoWriter(vid_name, 'Motion JPEG AVI');
v.Quality = 90;
v.FrameRate = FrameRate;
open(v);
set(fig, 'Color', 'k')
    % extract the three frames
    MT = zeros(height, width);
    editframe = BWvid(:,:,1);
    occProb = sum(BWvid(:,:,1:n),3)/n;
    vidFrame = ([vert, MT, vert, MT,vert, MT, vert]);
    vidFrame = [hori;vidFrame; hori];
    buff = zeros(size(vidFrame));
    % add scale bar
    imshow(buff)
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset') 
    
for n = 1:15
    
    set(fig, 'Color', 'k')
    % extract the three frames
    rawframe = mov(:,:,n);
    editframe = BWvid(:,:,n);
    occProb = sum(BWvid(:,:,1:n),3)/n;
    vidFrame = ([vert, rawframe, vert, editframe,vert, occProb, vert]);
    vidFrame = [hori;vidFrame; hori];
   
    % add scale bar
    imshow(vidFrame)
    cb = colorbar;
        %label and aesthetics
        cb.FontName = 'Arial';
        cb.FontSize = 10;
        cb.Color = 'w';
        cb.Label.String = 'Occupation probabilty/pixel';
    % TITLES:
    dim = [0.1860 0.93 0.0343 0.0535];
    a = annotation('textbox',dim,'String','Raw','FitBoxToText','on');
    set(a, 'color', 'w', 'EdgeColor', 'k')
    dim2 = [0.4568 0.93 0.0343 0.0535];
    a = annotation('textbox',dim2,'String','Binary','FitBoxToText','on');
    set(a, 'color', 'w', 'EdgeColor', 'k')
    dim3 = [0.7451 0.93 0.0343 0.0535];
    a = annotation('textbox',dim3,'String','Occupation','FitBoxToText','on');
    set(a, 'color', 'w', 'EdgeColor', 'k')

    % Collect image info for movie file
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')  
end

close(v)
close all
fprintf('\n Video Saved! \n')  

   
set(0,'DefaultFigureVisible','on');

%%


imshowpair(mov(:,:,1), BWvid(:,:,1), 'montage')

testframe = emptyImg - currFrame;

thresh = 2*std(reshape(testframe, numel(testframe),1));
testframe(testframe>thresh) = 1;
imshow(testframe)

figure; histogram(testframe); vline(thresh)

imshowpair(currFrame, testframe, 'montage') % look at raw and filtered image
imshow(testframe)



%non moving flies will count as background...need to just take a background picture...
background = repmat(emptyImg,[1,1,nframes]);
backSubImg = mov-background;

%sum each pixel's value over time:
total_val = sum(backSubImg,3);
minval = min(min(total_val));
maxval = max(max(total_val));

h = heatmap(total_val-minval);


% TODO : keep working on this *simple* version of pixel tracking...










