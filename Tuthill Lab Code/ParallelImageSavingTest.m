


imaqreset


fileName = 'E:\Basler Trig\TestFolder\Test';

p = gcp('nocreate'); %check for current parallel pool i.e. cameras or something
if isempty(p)
    p = parpool(6);
end

tic
A_vid = testparallel(4);
B_vid = testparallel(2);
C_vid = testparallel(3);
D_vid = testparallel(6);
E_vid = testparallel(5);
F_vid = testparallel(1);
toc


tic
for icam = 1:6
    f(icam) = parfeval(@testparallel,1,icam);
end

for idx = 1:6
    [completedIdx,value] = fetchNext(f);
    switch completedIdx
        case 1
            F_vid = value;
        case 2
            B_vid = value;
        case 3
            C_vid = value;
        case 4
            A_vid = value;
        case 5
            E_vid = value;
        case 6
            D_vid = value;
    end
    fprintf('Got result with index: %d.\n', completedIdx);
end
toc



tic
a = parfeval(@testparallel,1,4);
b = parfeval(@testparallel,1,2);
c = parfeval(@testparallel,1,3);
d = parfeval(@testparallel,1,5);
e = parfeval(@testparallel,1,6);
f = parfeval(@testparallel,1,1);


A_vid = fetchOutputs(a); % Block until a completes
B_vid = fetchOutputs(b);
C_vid = fetchOutputs(c);
D_vid = fetchOutputs(d);
E_vid = fetchOutputs(e);
F_vid = fetchOutputs(f);
toc

 
A_vid = fetchOutputs(a);





for idx = 1:6
  % fetchNext blocks until next results are available.
  [completedIdx,value] = fetchNext(f);
  magicResults{completedIdx} = value;
  fprintf('Got result with index: %d.\n', completedIdx);
end




%%

% v = videoinput('gige', 1, 'Mono8');
v = videoinput('gentl', 1, 'Mono8');
s = v.Source;
% % Set up for input wires from Basler Cam
% Basler_src.LineSelector = 'Line4';             % brings up settings for line4
% Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
% Basler_src.LineInverter = 'False';             % should be 'False'
% Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
% Basler_src.LineSelector = 'Line3';             % brings up settings for line3
% Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
% Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
% Basler_src.TriggerMode = 'Off';
% Basler_src.LineSelector = 'Line3';             % brings up settings for line3
% Basler_src.TriggerActivation = 'RisingEdge';
% Basler_src.TriggerMode = 'On';
% Basler_src.GainAuto = 'Off';
% Basler_src.ExposureTime = cam_spec.ExposureTime;                % Exposure setting for Basler
% Basler_src.Gain = cam_spec.Gain;          % Lighting gain
% Basler_src.Gamma = cam_spec.Gamma;         % White enhancement on video




% Set optimum camera streaming parameters PacketSize and PacketDelay
% Refer to http://www.mathworks.com/matlabcentral/answers/91834
% s.PacketSize = 9000;
% s.PacketDelay = 

v.FramesPerTrigger = 300;

v.FramesAcquiredFcnCount = 2;
v.FramesAcquiredFcn = {@saveImageData, fileName};

v.TimerFcn = @checkFuturesErrors;
v.TimerPeriod = 1;

tic
start(v);
pause(1)
% 
% hImg = preview(v);
% hFig = ancestor(hImg, 'figure');

% Run until preview window is manually closed
% waitfor(hFig);

stop(v);
delete(v)
toc
%%


% To request multiple evaluations, use a loop.
for idx = 1:10
  f(idx) = parfeval(p,@magic,1,idx); % Square size determined by idx
end

% Collect the results as they become available.
magicResults = cell(1,10);
for idx = 1:10
  % fetchNext blocks until next results are available.
  [completedIdx,value] = fetchNext(f);
  magicResults{completedIdx} = value;
  fprintf('Got result with index: %d.\n', completedIdx);
end









num.cams = 6;
param = load_cam_parameters(param, num, FramesPerTrigger);
for ii = 1:num.cams
    param.(['Cam' Alphabet(ii)]).ROI = param.(['Cam' Alphabet(ii)]).ROI_film;
end
param.num_cams = num.cams;



tic
for icam = 1:6
vid = testparallel(icam);
end
toc


f = parfeval(@testparallel,1,1);



p = gcp('nocreate'); %check for current parallel pool i.e. cameras or something

% To request multiple evaluations, use a loop.
for icam = 1:6
  f(icam) = parfeval(@testparallel,1,icam); 
end

f(icam) = parfeval(@testparallel,1,icam);
f(icam) = parfeval(@testparallel,1,icam);
f(icam) = parfeval(@testparallel,1,icam);
f(icam) = parfeval(@testparallel,1,icam);
f(icam) = parfeval(@testparallel,1,icam);
f(icam) = parfeval(@testparallel,1,icam);

% Collect the results as they become available.
magicResults = cell(1,10);
for idx = 1:10
  % fetchNext blocks until next results are available.
  [completedIdx,value] = fetchNext(f);
  magicResults{completedIdx} = value;
  fprintf('Got result with index: %d.\n', completedIdx);
end



parallel.FevalFuture



