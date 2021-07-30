
function param = load_parameters
% 
% param = load_parameters
% 
% Load the param for the experiement for timing
% arena stim, analoge channels, vid and light, etc.
% 
% ES Dickinson, University of Washington, 2019


% Genetic Cross Information:
param.cross = select_cross;


param.directory = 'E:\FicTrac Raw Data\';
% param.directory = 'F:\FicTrac Raw Data\';

param.num_reps =      3;  
param.vid_on =        1;    % 1 for on, 0 for off
param.light_length =  [0, 0.03, 0.06, 0.09, 0.18, 0.36, 0.720];  % LED trigger time for optigentic flies
% param.light_length =  [0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09];  % LED trigger time for optigentic flies
param.LED_light =  'red';

%start|end timing:
param.start_pause =   20;   % length of time to run tracking without stimulus at experiment start
param.end_pause =     5;     % pause after all stimulus

%light and video info
param.LED_intensity = 4; %1.8;     % Intensity of the LED light


param.basler_length = 2;     % Basler trigger length in seconds
param.basler_delay =  0.5;   % How long before laser to start the cameras
param.Basler_fps =    300;   % Desired fps from the basler cam
param.fps =           30;    % Pt. Grey camera fps
param.basler_samplerate = 10000; %sampling frequency of session
 
%arena stimiulus
param.Pattern = {'str', 'rot'}; %arena patterns
param.Dir = {'cw', 'ccw'};     %arena directions
param.num.pats =      2;     % number of patterns
param.dir =       [1 0];     % 'direction' of rotation for pattern
param.speeds = [68 -68];     % stimulus group D
param.width =         4;     % width of stripe in the arena (by pixel)
param.stripe_Pat =    2;     % wiggle pattern
param.grating_Pat =   3;     % dark stripe, bright background-rotating pattern


%testing light effect stimulus
param.test.Pattern = {'str', 'str'}; %arena patterns
param.test.Dir = {'cw', 'ccw'};     %arena directions
param.test.num.pats =      2;     % number of patterns
param.test.dir =       [1 0];     % 'direction' of rotation for pattern
param.test.speeds = [68 -68];     % stimulus group D
param.test.width =         4;     % width of stripe in the arena (by pixel)
param.test.stripe_Pat =    2;     % wiggle pattern
param.test.grating_Pat =   2;     % dark stripe, bright background-rotating pattern

%timing info
param.CL_time =            3;
param.OL_time =            2; %time before and time after the laser onset
param.control_time =       3;  
param.writetodiskpause =   5;

% Analogue channels
param.num_conds = param.num.pats*length(param.light_length)*length(param.dir);
param.num_inchans = 12; %former 9
param.chan_labels = {'trigger', 'x pos', 'y pos', 'cond sig', 'FicTrac heading',...
                          'LED trigger', 'A Cam', 'B Cam', 'C Cam', 'D Cam', 'E Cam', 'F Cam'};

% Temperature:
temp_opt = {'<75','76-79', '80-83', '84-87', '88-90', '>90'};
choice = listdlg('ListString', temp_opt, 'PromptString','Temperature?',...
                'SelectionMode', 'Single', 'ListSize', [80 150]);
% choice = menu('Temperature?', temp_opt);
param.temp = temp_opt(choice);

% Humidity: 
humidity_opt = {'<45','46-49', '50-53', '54-57', '58-60', '61-64', '65-67', '68-70', '71-74','>75'};
choice = listdlg('ListString', humidity_opt, 'PromptString','Humidity?',...
                'SelectionMode', 'Single', 'ListSize', [80 150]);
% choice = menu('Humidity?', humidity_opt);
param.humidity = humidity_opt(choice);

param.sex = questdlg('Select fly sex:', 'Fly Sex', 'Male', 'Female', 'Female');

param.body = questdlg('Select starting condition:', '', 'intact', 'headless', 'intact');
param.loading = questdlg('Select loading condition:', '', 'onball', 'offball', 'onball');

% Camera Parameters:
FramesPerTrigger = param.Basler_fps*param.basler_length;


num.camNums = [4, 2, 3, 6, 5, 1];
num.cams = 6;
param = load_cam_parameters(param, num, FramesPerTrigger);
for ii = 1:num.cams
    param.(['Cam' Alphabet(ii)]).ROI = param.(['Cam' Alphabet(ii)]).ROI_film;
end
param.num_cams = num.cams;

% Laser param
param.laser_freq = 1200; %Hz of laser
param.laser_ratio = 0.5; % 1:num ratio of light off
param.basler_volts = 9; %basler camera signal volts

fprintf('\n Parameters loaded \n') 


end




