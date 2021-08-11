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
folderNames = [folderNames, 'Today'];

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

%% Load background video
movieFullFileName = fullfile(folder, 'background_vid.avi');
% Check to see that it exists.
if ~exist(movieFullFileName, 'file')
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

movieInfo = VideoReader(movieFullFileName);
nframes = movieInfo.Duration * movieInfo.FrameRate;
height = movieInfo.Height;
width = movieInfo.Width;
mov = zeros(height, width, nframes);

%load all the data for the video:
for i = 1:nframes
    mov(:,:,i) = im2double(read(movieInfo, i));
end

% Pull the median pixel value in the background images
emptyImg = median(mov,3);
imshow(emptyImg)

backgroundImg = emptyImg - mov(:,:,1);
imshowpair(emptyImg, backgroundImg, 'montage')


%% Load selected video
movieFullFileName = fullfile(folder, 'test_vid_1.avi');
% Check to see that it exists.
if ~exist(movieFullFileName, 'file')
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

movieInfo = VideoReader(movieFullFileName);
nframes = movieInfo.Duration * movieInfo.FrameRate;
height = movieInfo.Height;
width = movieInfo.Width;
mov = zeros(height, width, nframes);

%load all the data for the video:
while hasFrame(movieInfo)
    mov(:,:,i) = im2double(readFrame(movieInfo));
end

% How do the pixel values change across time for a random selection?



% Can you identify the flies from the background image?
currFrame = mov(:,:,1);
imshowpair(emptyImg, currFrame, 'diff')

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










