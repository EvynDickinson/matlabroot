
%-----------------Capture A ---------------------

FramesPerTrigger = 300;
exposure_time = 2000;


% Cam Specs A: %400
cam_spec.num = 4;
cam_spec.ROI_positioning = [146, 201, 224, 185];
cam_spec.ROI_film = [0, 105, 516 450]; %[36 165 447 365]; 
cam_spec.ROI = [0, 105, 516 450];
% cam_spec..fighandle = [];
cam_spec.Gain = 6.334203998982872;
cam_spec.Gamma = 1;
% cam_spec..ExposureTime = 2000;
cam_spec.ExposureTime = exposure_time + 300;
cam_spec.FramesPerTrigger = FramesPerTrigger;


% Setup source and A_vid objects for Basler Camera
A_vid = videoinput('gentl', cam_spec.num, 'Mono8');
triggerconfig(A_vid, 'hardware');                 % set trigger to come from hardware i.e. line in
A_vid.LoggingMode = 'disk'; %'disk&memory';
A_vid.FramesPerTrigger = cam_spec.FramesPerTrigger;
A_vid.ROIPosition = cam_spec.ROI;            % fly-only video acq.
Basler_src = A_vid.Source;
% Basler_src = getselectedsource(A_vid);
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




tic
parfor icam = 1:6
   cam(icam).vid =  testparallel(icam);
   
 
end
toc
