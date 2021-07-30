
function [vid, Basler_src] = initiate_camera(cam_spec)
% 
% [vid, Basler_src] = load_cam_params(cam_spec)
% 
% Load the camera specs for a camera and set up the
% trigger and parameters
% 
% Inputs:
% 'cam_spec' [structure with num, ROI, Gain, Gamma, Exposure
%             for the desired camera linking]
%             
% Outputs:
% 'vid' [video input structure]
% 'Basler_src' [source for the Basler cam with controls]
% 
% ES Dickinson, University of Washington, 2019


% Setup source and vid objects for Basler Camera
vid = videoinput('gentl', cam_spec.num, 'Mono8');
triggerconfig(vid, 'hardware');                 % set trigger to come from hardware i.e. line in
vid.LoggingMode = 'disk'; %'disk&memory';
vid.FramesPerTrigger = cam_spec.FramesPerTrigger;
vid.ROIPosition = cam_spec.ROI;            % fly-only video acq.
Basler_src = vid.Source;
% Basler_src = getselectedsource(vid);
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

end

