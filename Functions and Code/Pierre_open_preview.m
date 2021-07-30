

%% parameters
% param_fname = 'calibration.toml';
% param_fname = 'recorder.toml';
param_fname = 'imaging_recording.toml';

%% load libraries
addpath C:/Users/Tuthill/Pierre/builds/matlab-toml

%% read the parameters
params = toml.read(param_fname);

%% dates
today = char(datetime('now','TimeZone','local', 'Format','y-MM-dd'));
currtime = char(datetime('now','TimeZone','local', 'Format','y-MM-dd--HH-mm-ss'));

outfolder = fullfile( params.DestFolder, today ) ;

if ~exist(outfolder,'dir')
    mkdir(outfolder)
end

fprintf("%s %s\n", today, currtime);

%% make backup of parameter file in session
[path, basename, ext] = fileparts(param_fname);
copyfile(param_fname, fullfile(outfolder, [basename ext]));

%% setup the cameras based on parameters

% get the cameras
imaqreset;
info = imaqhwinfo('gentl');

% setup defaults
default_present = isfield(params.cameras, 'default');

if default_present
    default_params = params.cameras.default;
else
    default_params = struct();
end

% setup cam names
cam_names = fieldnames(params.cameras);
cam_names = cam_names(~ismember(cam_names, 'default'));


vids = cell(0);


for cnum=1:length(cam_names)
    cam_name = cam_names{cnum};

    cam_params = mergestructs(default_params, params.cameras.(cam_name));

    % TODO: get the device id here from device info and cam_params.DeviceName
    device_id = -1;
    for dnum=1:length(info.DeviceInfo)
       if info.DeviceInfo(dnum).DeviceName == cam_params.DeviceName
            device_id = info.DeviceIDs{dnum};
       end
    end
    
    if device_id == -1
        error('Camera %s not found: %s', cam_name, cam_params.DeviceName)
    else
        fprintf('Camera %s found: %s\n', cam_name, cam_params.DeviceName);
    end
    
    % setup video
    vid = videoinput('gentl', device_id, cam_params.Format);
    triggerconfig(vid, 'hardware');  % set trigger to come from hardware

    vid.FramesPerTrigger = params.FPS * params.Length;

    f=fieldnames(cam_params.video);
    for i=1:length(f)
        vid.(f{i}) = cam_params.video.(f{i});
    end

%     disp(vid.ROIPosition) 
    
    % setup source
    src = getselectedsource(vid);
    src.AcquisitionFrameRate = params.FPS;

    f=fieldnames(cam_params.source);
    for i=1:length(f)
        src.(f{i}) = cam_params.source.(f{i});
    end
    
    % handle line options
    lines=fieldnames(cam_params.lines);
    for lnum=1:length(lines)
        line=lines{lnum};
        src.LineSelector = line;
        f=fieldnames(cam_params.lines.(line));
        for i=1:length(f)
            src.(f{i}) = cam_params.lines.(line).(f{i});
        end
    end
    
    % handle trigger options
    triggers=fieldnames(cam_params.triggers);
    for tnum=1:length(triggers)
        trigger=triggers{tnum};
        src.TriggerSelector = trigger;
        f=fieldnames(cam_params.triggers.(trigger));
        for i=1:length(f)
            src.(f{i}) = cam_params.triggers.(trigger).(f{i});
        end
    end

%     diskLogger = struct();
    outname = fullfile(outfolder, ...
        sprintf('%s_%s_%s.avi', params.Prefix, currtime, cam_name));
%     diskLogger = VideoWriter(outname, 'Motion JPEG AVI');
    diskLogger = VideoWriter(outname, cam_params.Compression);
    diskLogger.FrameRate = params.FPS; 

    if isfield(cam_params, 'Quality')
        diskLogger.Quality = cam_params.Quality; 
    end
    vid.DiskLogger = diskLogger;  

    vids{cnum} = vid;

end


%% setup the daq trigger

% % LED|Basler session
s_vid_light = daq.createSession('ni');
s_vid_light.Rate = 10000;
s_vid_light.DurationInSeconds = params.Length;

% % add analog output channels for LED|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler output

s_vid_light.addAnalogInputChannel('Dev1', 'ai0', 'Voltage');
s_vid_light.addAnalogInputChannel('Dev1', 'ai1', 'Voltage');
s_vid_light.addAnalogInputChannel('Dev1', 'ai2', 'Voltage');


basler_volts = 9;
basler_outsig = zeros(params.Length*s_vid_light.Rate, 1);
basler_rate = round(s_vid_light.Rate/params.FPS);
basler_outsig(1:basler_rate:end) = basler_volts;
% basler_outsig(100) = basler_volts;

queueOutputData(s_vid_light, basler_outsig)


%% do the recording

% flush the data
for cnum=1:length(vids)
    flushdata(vids{cnum});
end

closepreview;
close('all');

fig = figure('Renderer', 'painters', 'Position', [200 300 1500 500]);
 
for cnum=1:length(vids)

    %fig.Name = cam_names{cnum};
    % ax = axes(fig)
    ax = subplot_tight(1, length(vids), cnum); 
    vidRes = get(vids{cnum}, 'ROIPosition');
    nBands = get(vids{cnum}, 'NumberOfBands');
    im = image(ax, zeros(vidRes(4), vidRes(3), nBands) );
    axis(ax,'image');
    preview(vids{cnum}, im);
    title(cam_names{cnum});
    h = gca;
    h.Visible = 'On';
end

fprintf("starting recording...\n");

% start recording
for cnum=1:length(vids)
    start(vids{cnum});
end

% start trigger
data = s_vid_light.startForeground();
% s_vid_light.startBackground();

% wait 
% pause(params.Length);

fprintf("done recording!\n");

% stop recording
for cnum=1:length(vids)
    stop(vids{cnum});
end


% make sure all the video is saved
for cnum=1:length(vids)
    while (vids{cnum}.FramesAcquired > vids{cnum}.DiskLoggerFrameCount)
        pause(.1)
    end
end


closepreview;
close('all');

imaqreset;