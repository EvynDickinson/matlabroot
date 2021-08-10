


function [vid, cam_spec] = testparallel(camnum)

exposure_time = 2000;
FramesPerTrigger = 300;
fprintf(['\n Cam: ' num2str(camnum)])
switch camnum
    case 4
    % Cam Specs A: %400
    cam_spec.num = 4;
    cam_spec.ROI_positioning = [146, 201, 224, 185];
    cam_spec.ROI_film = [0, 105, 516 450]; %[36 165 447 365]; 
    cam_spec.Gain = 6.334203998982872;
    cam_spec.Gamma = 1;
    cam_spec.ExposureTime = exposure_time + 300;
    
    case 2
    % Cam Specs B: %395
    cam_spec.num = 2;
    cam_spec.ROI_positioning = [398 337 228 178]; 
    cam_spec.ROI_film = [268 284 458 355]; %[279 310 408 297]; 
    cam_spec.Gain = 5.9242674326438092;
    cam_spec.Gamma = 0.92999267578125;
    cam_spec.ExposureTime = exposure_time;
    
    case 3
    % Cam Specs C: %083
    cam_spec.num = 3;
    cam_spec.ROI_positioning = [283 98 354 293]; 
    cam_spec.ROI_film = [12 53 820 544]; %[92 98 677 448];   
    cam_spec.Gain = 8.9931973119681032;  
    cam_spec.Gamma = 0.6869964599609375;
    cam_spec.ExposureTime = exposure_time;
    
    case 6
    % Cam Specs D: %483
    cam_spec.num = 6; 
    %cam_spec.ROI_film = [138 25 380 461];  
    cam_spec.ROI_film = [105 0 434 521];
    cam_spec.ROI_positioning = [48 0 537 543]; %[265 134 135 237]; 
    cam_spec.Gain = 4.3728721948228992;  
    cam_spec.Gamma = 0.75;
    cam_spec.ExposureTime = exposure_time;
    
    case 5
    % Cam Specs E: %635
    cam_spec.num = 5; 
    cam_spec.ROI_positioning = [206 252 229 131];
    cam_spec.ROI_film = [61 142 550 425]; %[117 240 424 272];
    cam_spec.Gain = 8.9931973119681032;  
    cam_spec.Gamma = 0.875;
    cam_spec.ExposureTime = exposure_time;
    
    case 1
    % Cam Specs F: %652
    cam_spec.num = 1;
    cam_spec.ROI_positioning = [448 361 191 130];
    cam_spec.ROI_film = [323 279 450 344]; %[371 334 355 260]; 
    cam_spec.Gain = 9.876538964450301;  
    cam_spec.Gamma = 1.00;
    cam_spec.ExposureTime = exposure_time;

end

cam_spec.FramesPerTrigger = FramesPerTrigger;



% Setup source and vid objects for Basler Camera
vid = videoinput('gentl', camnum, 'Mono8');
triggerconfig(vid, 'hardware');                 % set trigger to come from hardware i.e. line in
% vid.LoggingMode =  'disk&memory'; %'disk';

vid.ROIPosition = cam_spec.ROI_film;            % fly-only video acq.
Basler_src = getselectedsource(vid);
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

cam_spec.FramesAcquiredFcnCount = 2;

fprintf(['\n Camera ' num2str(cam_spec.num) ' configuration completed \n'])

