
cam_spec.Basler_fps = 300;
cam_spec.light_length = 0.5; 
cam_spec.basler_length = cam_spec.light_length+2.5;
cam_spec.basler_delay = 0.5;    %sec of basler running before light comes on  
num.cams = 6;
FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger);
for icam = 1:num.cams
    cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_film; %ROI_positioning;
end 
% Grab Cameras and Initiate Set up:
[A_vid, A_Basler_src] = initiate_camera(cam_spec.CamA);
[B_vid, B_Basler_src] = initiate_camera(cam_spec.CamB);
[C_vid, C_Basler_src] = initiate_camera(cam_spec.CamC);
[D_vid, D_Basler_src] = initiate_camera(cam_spec.CamD);
[E_vid, E_Basler_src] = initiate_camera(cam_spec.CamE);
[F_vid, F_Basler_src] = initiate_camera(cam_spec.CamF);


A_Basler_src.TriggerMode = 'Off'; % preview(A_vid)
B_Basler_src.TriggerMode = 'Off'; %preview(B_vid)
C_Basler_src.TriggerMode = 'Off';% preview(C_vid)
D_Basler_src.TriggerMode = 'Off'; %preview(D_vid)
E_Basler_src.TriggerMode = 'Off'; %preview(E_vid)
F_Basler_src.TriggerMode = 'Off'; %preview(F_vid)

z = 1.2; %zoom ratio for previews

% load into windows:
AFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam A 400 gentl-4');
n = cam_spec.CamA.ROI_film(4);
m = cam_spec.CamA.ROI_film(3);
hImage = image(zeros(n, m, 1));
AFig.Position = [738 610 m*z n*z];
preview(A_vid, hImage);


BFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam B 395 gentl-2');
n = cam_spec.CamB.ROI_film(4);
m = cam_spec.CamB.ROI_film(3);
hImage = image(zeros(n, m, 1));
BFig.Position = [1338 633 m*(z+0.2) n*(z+0.2)];
preview(B_vid, hImage);


CFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam C 083 gentl-3');
n = cam_spec.CamC.ROI_film(4);
m = cam_spec.CamC.ROI_film(3);
hImage = image(zeros(n, m, 1));
CFig.Position = [37 601 m n];
preview(C_vid, hImage);


DFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam D 483 gentl-6');
n = cam_spec.CamD.ROI_film(4);
m = cam_spec.CamD.ROI_film(3);
hImage = image(zeros(n, m, 1));
DFig.Position = [894 98 m n];
preview(D_vid, hImage);


EFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam E 635 gentl-5');
n = cam_spec.CamE.ROI_film(4);
m = cam_spec.CamE.ROI_film(3);
hImage = image(zeros(n, m, 1));
EFig.Position = [449 285 m n];
preview(E_vid, hImage);


FFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam F 652 gentl-1');
n = cam_spec.CamF.ROI_film(4);
m = cam_spec.CamF.ROI_film(3);        
hImage = image(zeros(n, m, 1));
FFig.Position = [10 245 m*z  n*z];
preview(F_vid, hImage);



%%

% Load config file with all figures
load([handles.settings.dir_scripts,'/trial_viewer_config.mat']);

% Adjust main window
pos_fig.main = handles.figure1.Position;

% Iterate through other open figures
names_fig  = fieldnames(handles.fig);
for iFigure = 1:numel(fieldnames(handles.fig))
   if isfield(pos_fig,names_fig{iFigure})
       if isvalid(handles.fig.(names_fig{iFigure}))
           pos_fig.(names_fig{iFigure}) = handles.fig.(names_fig{iFigure}).Position;
       end
   end
end

% Save config file
save([handles.settings.dir_scripts,'trial_viewer_config.mat'],'pos_fig')
fprintf(['\t Display settings saved as ',handles.settings.dir_scripts,'trial_viewer_config.mat \n'])

% Update handles structure
guidata(hObject, handles);


% GUI.handles position = [1293 43 621 496];


