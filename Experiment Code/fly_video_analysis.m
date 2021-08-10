
baseFolder = ['E:\My Drive\Jeanne Lab\DATA\8.10.2021'];



%% Load new video

folder = fullfile(baseFolder); % TODO add the date or other dynamic title bits here later
movieFullFileName = fullfile(folder, 'test_vid.avi');
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
tic
for i = 1:nframes
    mov(:,:,i) = im2double(read(movieInfo, i));
end
toc


imshow(img)



% try plotting the mediam greyscale value for each pixel...are the flies
% invisible now?

emptyImg = median(mov,3);
%non moving flies will count as background...need to just take a background picture...

background = repmat(emptyImg,[1,1,nframes]);

backSubImg = mov-background;


imshow(median(backSubImg,3)) 

% TODO : keep working on this *simple* version of pixel tracking...












