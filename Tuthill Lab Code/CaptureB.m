

%-----------------Capture B ---------------------
FramesPerTrigger = 300;
exposure_time = 2000;


% Cam Specs B: %395
cam_spec.num = 2;
cam_spec.ROI_positioning = [398 337 228 178]; 
cam_spec.ROI_film = [268 284 458 355]; 
cam_spec.ROI = [268 284 458 355];
cam_spec.Gain = 5.9242674326438092;
cam_spec.Gamma = 0.92999267578125;
% cam_spec.ExposureTime = 2000;
cam_spec.ExposureTime = exposure_time;
cam_spec.FramesPerTrigger = FramesPerTrigger;


% Setup source and B_vid objects for Basler Camera
B_vid = videoinput('gentl', cam_spec.num, 'Mono8');
triggerconfig(B_vid, 'hardware');                 % set trigger to come from hardware i.e. line in
B_vid.LoggingMode = 'disk'; %'disk&memory';
B_vid.FramesPerTrigger = cam_spec.FramesPerTrigger;
B_vid.ROIPosition = cam_spec.ROI;            % fly-only video acq.
Basler_src = B_vid.Source;
% Basler_src = getselectedsource(B_vid);
Basler_src.TriggerMode = 'Off';

% Set up for input wires from Basler Cam
Basler_src.LineSelector = 'Line4';             % brings up settings for line4
Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
Basler_src.LineInverter = 'False';             % should be 'False'
Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
Basler_src.LineSelector = 'Line3';             % brings up settings for line3
Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
Basler_src.TriggerMode = 'Off';
Basler_src.LineSelector = 'Line3';             % brings up settings for line3
Basler_src.TriggerActivation = 'RisingEdge';
Basler_src.TriggerMode = 'On';
Basler_src.GainAuto = 'Off';
Basler_src.ExposureTime = cam_spec.ExposureTime;                % Exposure setting for Basler
Basler_src.Gain = cam_spec.Gain;          % Lighting gain
Basler_src.Gamma = cam_spec.Gamma;         % White enhancement on video

fprintf(['\n Camera ' num2str(cam_spec.num) ' configuration completed \n'])
