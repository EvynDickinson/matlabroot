function hf = vidplayback(VideoName, FrameRate, starttime, stoptime, zoom)
% hf = vidplayback(VideoName, FrameRate, starttime, stoptime)
% Min input is VideoName, complete with address to find it
% Defaults: 
% FrameRate = 50;
% starttime = 0.3;
% stoptime = 0.5;

switch nargin
    case 1
        FrameRate = 50;
        starttime = 0.3;
        stoptime = 1.5;
    case 2
        starttime = 0.3;
        stoptime = 1.5;
    case 3
        stoptime = 1.5;
end

% % Create a VideoReader to use for reading frames from the file.
flyVideo = VideoReader(VideoName);
flyVideo.CurrentTime = starttime;
stptime = stoptime;
duration = round((stptime-flyVideo.CurrentTime)*300+1);
% 
% height = flyVideo.Height*zoom;
% width = flyVideo.Width*zoom;

mov(duration) = struct('cdata',zeros(flyVideo.Height,flyVideo.Width,3,'uint8'),...
    'colormap',[]);

% mov(duration) = struct('cdata',zeros(height, width, 3, 'uint8'), 'colormap',[]);

hf = figure;
set(hf,'position',[151,433 flyVideo.Width flyVideo.Height]);
set(hf, 'Name', VideoName(end-45:end))

if nargin < 5
% NO SCALING
    ii = 1;
    currAxes = axes;
    while flyVideo.CurrentTime <= stptime
       mov(ii).cdata = readFrame(flyVideo);
       image(mov(ii).cdata, 'Parent', currAxes);
       currAxes.Visible = 'off';
       pause(1/FrameRate);
       ii = ii+1;
    end
else
% SCALING
    ii = 1;
    currAxes = axes;
    while flyVideo.CurrentTime <= stptime
       pre = readFrame(flyVideo);
       post = imresize(pre, zoom);
       mov(ii).cdata = post;
       image(mov(ii).cdata, 'Parent', currAxes);
       currAxes.Visible = 'off';
       pause(1/FrameRate);
       ii = ii+1;
    end
end


