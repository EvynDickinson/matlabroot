function data = AcquireTrace_shock(expnumber,numtrials,trialduration,ODORstimulustime,SHOCKstimulustime,ODORdur,SHOCKdur)
% data = AcquireTrace_shock(expnumber,numtrials,trialduration,ODORstimulustime,SHOCKstimulustime,ODORdur,SHOCKdur)
% expnumber = 1;
% trialduration = 10; %in seconds;
% numtrials = 1;
% ODORstimulustime % time when the odor stim starts into trial (s)
% SHOCKstimulustime % when the shock starts (s);
% ODOR_dur % duration of odor stimulus (s)
% SHOCK_dur % duration of shock (s) -- this is actually set by the AM systems pulse
    
% shorthand
fs = 10E3;
fullData = zeros(trialduration*fs, numtrials);

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
    data(n).LEDstart = 1.0;         % LED triggers on 1.5 seconds into trial
    data(n).camerastart = 2.0;      % camera triggers on @ 2.0 seconds     
    data(n).SHOCK_dur = SHOCKdur;
    data(n).ODOR_dur = ODORdur;
    
    %make odor stimulus vector
    trailingdurationODOR = trialduration - ODORdur - ODORstimulustime;
    data(n).ODORstimulusVector = [zeros(1,fs*ODORstimulustime)...
                                  5*ones(1,fs*ODORdur)...
                                  zeros(1,fs*trailingdurationODOR)]';
    
    %make electric shock stimulus vector
    trailingdurationSHOCK = trialduration - SHOCKdur - SHOCKstimulustime;   
    data(n).SHOCKTrigger = [zeros(1,fs*SHOCKstimulustime)...
                            5*ones(1,fs*SHOCKdur)...
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
    
    output_signal = [data(n).LEDTrigger data(n).ODORstimulusVector data(n).SHOCKTrigger  data(n).cameraTrigger];
    
%     figure; plot(output_signal); ylim([-0.5,5.5])
%     legend({'LED', 'Odor', 'Shock', 'Camera'})
    
    
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
    s.queueOutputData(output_signal);
    x = s.startForeground();
    
%     fullData(:,k) = x(:,1); % 2/19 changed 2 to 1

end

% average across rows
% data = mean(fullData, 2);

%plot(timevec, lowpass(x(:,2),150))
% figure(1)
% plot(timevec, data)
% figure(2)
% plot(timevec, sgolayfilt(data,1,17)) % random fit from the MATLAB curve smoothing website
% figure(3)
% plot(timevec, sgolayfilt(data,1,121))
% figure(4)
% plot(timevec, lowpass(data,0.5)) 
% figure(5)
% windowSize = 150; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% plot(timevec, filter(b,a,sgolayfilt(data,1,71))) % plot filtered PID waveform (channel 2) and stimulus 1 
% title(['Odor stimulus waveform (PID measurement)'])

% save([basepath,'/WCwaveform_' data(n).date,'_E',num2str(expnumber)],'data');