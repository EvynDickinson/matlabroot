

%% Queue blue light laser signal:
fs = 10E3;
LED_length = 5;

% Load DAQ session:
s = daq.createSession('ni');
s.addAnalogOutputChannel('Dev1',[0:3] , 'Voltage'); %LED,ODOR,SHOCK,CAM
s.Rate = fs;
disp('DAQ connected')

%make LED stimulus vector
LED_volt = 5;
outsig = LED_volt*ones(fs*LED_length,1);
outsig(end,1) = 0;

%make camera trig vector

cam_sig = 5*ones(fs*LED_length,1);
cam_sig(end,1) = 0;

% dummy channels for odor, shock, cam:
dummy = zeros(fs*LED_length,1);
data = [outsig, dummy, dummy, cam_sig];


%load data to LED:
queueOutputData(s, data) 
startBackground(s) 


s.queueOutputData(data);
x = s.startForeground();

% great -- LED is just like the old laser -- need to end on a the 0 output to
% turn off the laser
   

preview(vid)
closepreview(vid)
imaqreset

%% set trial parameters
for k = 1:numtrials
    startingtrial = 1; 
    todaysDate = datestr(datetime('now'), 'yyyy-mm-dd');
    n = startingtrial-1+k;
    data(n).date = todaysDate;                            % experiment date
    data(n).expnumber = expnumber;                        % experiment number
    data(n).trial = n;                                    % trial number
    % data(n).sampleTime = clock;
    data(n).acquisition_filename = mfilename('fullpath'); % saves name of mfile that generated data
    % sampling rates
    data(n).sampratein = 10E3;                            % input sample rate
    data(n).samprateout = fs;                             % output sample rate becomes input rate as well when both input and output present
    data(n).trialduration = trialduration;                % trial duration
    data(n).ODORstimulustime = ODORstimulustime;
    data(n).LEDstart = 1.5;         % LED triggers on 1.5 seconds into trial
    data(n).camerastart = 2.0;      % camera triggers on 2.0 seconds into trial     
    
    %make odor stimulus vector
    trailingdurationODOR = trialduration - stimulusduration - ODORstimulustime;
    data(n).ODORstimulusVector = [zeros(1,fs*ODORstimulustime)...
                                  5*ones(1,fs*stimulusduration)...
                                  zeros(1,fs*trailingdurationODOR)]';
    
    %make electric shock stimulus vector
    trailingdurationSHOCK = trialduration - stimulusduration - SHOCKstimulustime;   
    data(n).SHOCKTrigger = [zeros(1,fs*SHOCKstimulustime)...
                            5*ones(1,fs*stimulusduration)...
                            zeros(1,fs*trailingdurationSHOCK)]';

    %Make vector to trigger the camera to start recording.
    % here, camera triggers at second 2 of the 10 second length trial
    cameraon_time = data(n).trialduration*fs - data(n).camerastart*fs - 1;
    data(n).cameraTrigger = [zeros(1,data(n).camerastart*fs) 5*ones(1,cameraon_time) zeros(1,1)]';
    
    %Make vector to turn on the light source shortly before starting to
    %record.
    % LED triggers at second 1.5
    LEDon_time = data(n).trialduration*fs - data(n).LEDstart*fs - 1;
    data(n).LEDTrigger = [zeros(1,data(n).LEDstart*fs) 5*ones(1,LEDon_time) zeros(1,1)]';
    
    dt = 1/data(n).sampratein;
    timevec = dt:dt:data(n).trialduration;
    %   CHANNEL SET-UP
    %   0 VAR OUT
    %   1 Im  Through conditioner - gain 10X filtering 2K
    %   2 10Vm Through conditioner - gain 1 X filtering 2K
    %   3 Mode
    %   4 Gain
    s = daq.createSession('ni');
    s.addAnalogInputChannel('Dev1',[0:2],'Voltage');
    for i=1:3
        s.Channels(1,i).InputType = 'SingleEnded';
    end
    s.addAnalogOutputChannel('Dev1',[0:3] , 'Voltage');
    s.Rate = fs;
    disp('DAQ connected')
    %Wiring:
    %        AO0 = LED trigger (command to turn light on)
    %        AO1 = ODOR stimulus (command to odor delivery device)
    %        AO2 = SHOCK stimulus (command to trigger pulse generator)
    %        AO3 = camera trigger (command to start acquisition from camera)
    s.queueOutputData([data(n).LEDTrigger data(n).ODORstimulusVector data(n).SHOCKTrigger  data(n).cameraTrigger]);
    x = s.startForeground();
    
%     fullData(:,k) = x(:,1); % 2/19 changed 2 to 1

end