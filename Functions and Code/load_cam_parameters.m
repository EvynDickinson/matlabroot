

function cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger)
% 
% cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger)
% 
% Load camera parameters to faciliate setting up 
% 'initiate_cameras' function
%     
% ES Dickinson
% 


% Camera ROIs for calibration: 
full_spec = [0 0 832 632];
for icam = 1:6
    cam_spec.(['Cam' Alphabet(icam)]).ROI_FULL = full_spec;
end

exposure_time = 2000;
camNum = num.camNums;

% Cam Specs A: %400
cam_spec.CamA.num = camNum(1); %4
cam_spec.CamA.ROI_positioning = [146, 201, 224, 185];
cam_spec.CamA.ROI_film = [0, 105, 516 450]; %[36 165 447 365]; 
% cam_spec.CamA.fighandle = [];
cam_spec.CamA.Gain = 6.334203998982872;
cam_spec.CamA.Gamma = 1;
% cam_spec.CamA.ExposureTime = 2000;
cam_spec.CamA.ExposureTime = exposure_time + 300;
cam_spec.CamA.FramesPerTrigger = FramesPerTrigger;

% Cam Specs B: %395
cam_spec.CamB.num = camNum(2);%2;
cam_spec.CamB.ROI_positioning = [398 337 228 178]; 
cam_spec.CamB.ROI_film = [268 284 458 355]; %[279 310 408 297]; 
cam_spec.CamB.Gain = 5.9242674326438092;
cam_spec.CamB.Gamma = 0.92999267578125;
% cam_spec.CamB.ExposureTime = 2000;
cam_spec.CamB.ExposureTime = exposure_time;
cam_spec.CamB.FramesPerTrigger = FramesPerTrigger;

% Cam Specs C: %083
cam_spec.CamC.num = camNum(3);%3;
cam_spec.CamC.ROI_positioning = [283 98 354 293]; 
cam_spec.CamC.ROI_film = [12 53 820 544]; %[92 98 677 448];   
cam_spec.CamC.Gain = 8.9931973119681032;  
cam_spec.CamC.Gamma = 0.6869964599609375;
% cam_spec.CamC.ExposureTime = 2000;
cam_spec.CamC.ExposureTime = exposure_time;
cam_spec.CamC.FramesPerTrigger = FramesPerTrigger;

% Cam Specs D: %483
cam_spec.CamD.num = camNum(4);%6; 
%cam_spec.CamD.ROI_film = [138 25 380 461];  
cam_spec.CamD.ROI_film = [105 0 434 521];
cam_spec.CamD.ROI_positioning = [48 0 537 543]; %[265 134 135 237]; 
cam_spec.CamD.Gain = 4.3728721948228992;  
cam_spec.CamD.Gamma = 0.75;
% cam_spec.CamD.ExposureTime = 2000;
cam_spec.CamD.ExposureTime = exposure_time;
cam_spec.CamD.FramesPerTrigger = FramesPerTrigger;

% Cam Specs E: %635
cam_spec.CamE.num = camNum(5);%5; 
cam_spec.CamE.ROI_positioning = [206 252 229 131];
cam_spec.CamE.ROI_film = [61 142 550 425]; %[117 240 424 272];
cam_spec.CamE.Gain = 8.9931973119681032;  
cam_spec.CamE.Gamma = 0.875;
% cam_spec.CamE.ExposureTime = 2000;
cam_spec.CamE.ExposureTime = exposure_time;
cam_spec.CamE.FramesPerTrigger = FramesPerTrigger;

% Cam Specs F: %652
cam_spec.CamF.num = camNum(6);%1;
cam_spec.CamF.ROI_positioning = [448 361 191 130];
cam_spec.CamF.ROI_film = [323 279 450 344]; %[371 334 355 260]; 
cam_spec.CamF.Gain = 9.876538964450301;  
cam_spec.CamF.Gamma = 1.00;
% cam_spec.CamF.ExposureTime = 2000;
cam_spec.CamF.ExposureTime = exposure_time;
cam_spec.CamF.FramesPerTrigger = FramesPerTrigger;







% % Camera ROIs for calibration: 
% full_spec = [0 0 832 632];
% for icam = 1:num.cams
%     cam_spec.(['Cam' Alphabet(icam)]).ROI_FULL = full_spec;
% end
% 
% exposure_time = 2000;
% 
% 
% % Cam Specs A: %400
% cam_spec.CamA.num = 4;
% cam_spec.CamA.ROI_positioning = [146, 201, 224, 185];
% cam_spec.CamA.ROI_film = [0, 105, 516 450]; %[36 165 447 365]; 
% % cam_spec.CamA.fighandle = [];
% cam_spec.CamA.Gain = 6.334203998982872;
% cam_spec.CamA.Gamma = 1;
% % cam_spec.CamA.ExposureTime = 2000;
% cam_spec.CamA.ExposureTime = exposure_time + 300;
% cam_spec.CamA.FramesPerTrigger = FramesPerTrigger;
% 
% % Cam Specs B: %395
% cam_spec.CamB.num = 2;
% cam_spec.CamB.ROI_positioning = [398 337 228 178]; 
% cam_spec.CamB.ROI_film = [268 284 458 355]; %[279 310 408 297]; 
% cam_spec.CamB.Gain = 5.9242674326438092;
% cam_spec.CamB.Gamma = 0.92999267578125;
% % cam_spec.CamB.ExposureTime = 2000;
% cam_spec.CamB.ExposureTime = exposure_time;
% cam_spec.CamB.FramesPerTrigger = FramesPerTrigger;
% 
% % Cam Specs C: %083
% cam_spec.CamC.num = 3;
% cam_spec.CamC.ROI_positioning = [283 98 354 293]; 
% cam_spec.CamC.ROI_film = [12 53 820 544]; %[92 98 677 448];   
% cam_spec.CamC.Gain = 8.9931973119681032;  
% cam_spec.CamC.Gamma = 0.6869964599609375;
% % cam_spec.CamC.ExposureTime = 2000;
% cam_spec.CamC.ExposureTime = exposure_time;
% cam_spec.CamC.FramesPerTrigger = FramesPerTrigger;
% 
% % Cam Specs D: %483
% cam_spec.CamD.num = 6; 
% %cam_spec.CamD.ROI_film = [138 25 380 461];  
% cam_spec.CamD.ROI_film = [105 0 434 521];
% cam_spec.CamD.ROI_positioning = [48 0 537 543]; %[265 134 135 237]; 
% cam_spec.CamD.Gain = 4.3728721948228992;  
% cam_spec.CamD.Gamma = 0.75;
% % cam_spec.CamD.ExposureTime = 2000;
% cam_spec.CamD.ExposureTime = exposure_time;
% cam_spec.CamD.FramesPerTrigger = FramesPerTrigger;
% 
% % Cam Specs E: %635
% cam_spec.CamE.num = 5; 
% cam_spec.CamE.ROI_positioning = [206 252 229 131];
% cam_spec.CamE.ROI_film = [61 142 550 425]; %[117 240 424 272];
% cam_spec.CamE.Gain = 8.9931973119681032;  
% cam_spec.CamE.Gamma = 0.875;
% % cam_spec.CamE.ExposureTime = 2000;
% cam_spec.CamE.ExposureTime = exposure_time;
% cam_spec.CamE.FramesPerTrigger = FramesPerTrigger;
% 
% % Cam Specs F: %652
% cam_spec.CamF.num = 1;
% cam_spec.CamF.ROI_positioning = [448 361 191 130];
% cam_spec.CamF.ROI_film = [323 279 450 344]; %[371 334 355 260]; 
% cam_spec.CamF.Gain = 9.876538964450301;  
% cam_spec.CamF.Gamma = 1.00;
% % cam_spec.CamF.ExposureTime = 2000;
% cam_spec.CamF.ExposureTime = exposure_time;
% cam_spec.CamF.FramesPerTrigger = FramesPerTrigger;
% 



end