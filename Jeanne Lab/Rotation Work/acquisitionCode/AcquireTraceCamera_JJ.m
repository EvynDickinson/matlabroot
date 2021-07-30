function data = AcquireTraceCamera_JJ(expnumber,numtrials,trialduration,stimulustime,stimulusduration,bursts)
% expnumber = 1;
% trialduration = 10; %in seconds;
% numtrials = 1;
% stimulustime = 4; %in seconds;
% stimulusduration = 2.0; %in seconds;

fullData = zeros(100000, numtrials);

%% set trial parameters
for k = 1:numtrials
    startingtrial = 1; 
    todaysDate = datestr(datetime('now'), 'yyyy-mm-dd');
    n = startingtrial-1+k;
    data(n).date = todaysDate;                                 % experiment date
    data(n).expnumber = expnumber;                          % experiment number
    data(n).trial = n;                                        % trial number
    % data(n).sampleTime = clock;
    data(n).acquisition_filename = mfilename('fullpath');    % saves name of mfile that generated data
    % sampling rates
    data(n).sampratein = 10000;                              % input sample rate
    data(n).samprateout = 10000;                           % output sample rate becomes input rate as well when both input and output present
    data(n).trialduration = trialduration;                            % trial duration
    data(n).stimulustime = stimulustime;
    data(n).LEDstart = 1.5;         % LED triggers on 1.5 seconds into trial
    data(n).camerastart = 2.0;      % camera triggers on 2.0 seconds into trial     
    
    %make column vector for use as master8 trigger
    trailingdurationSTIM = trialduration - stimulusduration - stimulustime;
    trailingdurationCOMMAND = trialduration - 2;
    
    % create vector for odor stimulus pattern
    switch bursts  % time from stim start to trial end is 60000 samples
        case '250' % single 250 ms pulse
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,2500) zeros(1,57500)]';
        case '500' % single 500 ms pulse
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,55000)]'; 
        case '2000' % single 2000 ms pulse
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,data(n).samprateout*stimulusduration) zeros(1,data(n).samprateout*trailingdurationSTIM)]';
        case '500-250-500' % two pulses, length 500, ISI of 250 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,2500) 5*ones(1,5000) zeros(1,47500)]';
        case '500-500-500' % two pulses, length 500, ISI of 500 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,5000) 5*ones(1,5000) zeros(1,45000)]';
        case '500-1000-500' % two pulses, length 500, ISI of 1000 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,10000) 5*ones(1,5000) zeros(1,40000)]';
        case '2000-250-500' % two pulses, first 2000, second 500, ISI of 250 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,20000) zeros(1,2500) 5*ones(1,5000) zeros(1,32500)]';
        case '2000-500-500' % two pulses, first 2000, second 500, ISI of 500 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,20000) zeros(1,5000) 5*ones(1,5000) zeros(1,30000)]';
        case '2000-1000-500' % two pulses, first 2000, second 500, ISI of 1000 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,20000) zeros(1,10000) 5*ones(1,5000) zeros(1,25000)]';
        case '500-250-500-250-500' % three pulses, length 500, ISI of 250 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,2500) 5*ones(1,5000) zeros(1,2500) ...
                                        5*ones(1,5000) zeros(1,data(n).samprateout*trailingdurationSTIM)]';
        case '500-500-500-500-500' % three pulses, length 500, ISI of 500 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,5000) 5*ones(1,5000) zeros(1,5000) ...
                                        5*ones(1,5000) zeros(1,35000)]';   
        case '500-100-500-100-500' % three pulses, length 500, ISI of 100 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,1000) 5*ones(1,5000) zeros(1,1000) ...
                                        5*ones(1,5000) zeros(1,43000)]';
        case '500-250-500-250-500-250-500' % four pulses, length 500, ISI of 250 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,5000) zeros(1,2500) 5*ones(1,5000) zeros(1,2500) ...
                                        5*ones(1,5000) zeros(1,2500) 5*ones(1,5000) zeros(1,32500)]';
        case '250-250-250-250-250-250-250' % four pulses, length 250, ISI of 250 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,2500) zeros(1,2500) 5*ones(1,2500) zeros(1,2500) ...
                                        5*ones(1,2500) zeros(1,2500) 5*ones(1,2500) zeros(1,42500)]';
        case '100-100-100-100-100-100-100' % four pulses, length 100, ISI of 100 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,1000) zeros(1,1000) 5*ones(1,1000) zeros(1,1000) ...
                                        5*ones(1,1000) zeros(1,1000) 5*ones(1,1000) zeros(1,53000)]';
        case '50-50-50-50-50-50-50' % four pulses, length 50, ISI of 50 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,500) zeros(1,500) 5*ones(1,500) zeros(1,500) ...
                                        5*ones(1,500) zeros(1,500) 5*ones(1,500) zeros(1,56500)]';
        case '20-20-20-20-20-20-20' % four pulses, length 20, ISI of 20 ms
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,200) zeros(1,200) 5*ones(1,200) zeros(1,200) ...
                                        5*ones(1,200) zeros(1,200) 5*ones(1,200) zeros(1,58600)]';
        otherwise % single 2.0 s pulse
            data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,data(n).samprateout*stimulusduration) zeros(1,data(n).samprateout*trailingdurationSTIM)]';
    end
    
    % here, camera triggers at second 2 of the 10 second length trial
    cameraon_time = data(n).trialduration*data(n).samprateout - data(n).camerastart*data(n).samprateout - 1;
    data(n).cameraTrigger = [zeros(1,data(n).camerastart*data(n).samprateout) 5*ones(1,cameraon_time) zeros(1,1)]';
    % LED triggers at second 1.5
    LEDon_time = data(n).trialduration*data(n).samprateout - data(n).LEDstart*data(n).samprateout - 1;
    data(n).LEDTrigger = [zeros(1,data(n).LEDstart*data(n).samprateout) 5*ones(1,LEDon_time) zeros(1,1)]';
    
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
    s.addAnalogOutputChannel('Dev1',[1:3] , 'Voltage');
    s.Rate = data(n).samprateout;
    
    s.queueOutputData([data(n).stimulusVector  data(n).LEDTrigger  data(n).cameraTrigger]);
    x = s.startForeground();
    
    fullData(:,k) = x(:,1); % 2/19 changed 2 to 1
    
end

% average across rows
data = mean(fullData, 2);

%plot(timevec, lowpass(x(:,2),150))
figure(1)
plot(timevec, data)
figure(2)
plot(timevec, sgolayfilt(data,1,17)) % random fit from the MATLAB curve smoothing website
figure(3)
plot(timevec, sgolayfilt(data,1,121))
figure(4)
plot(timevec, lowpass(data,0.5)) 
figure(5)
windowSize = 150; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
plot(timevec, filter(b,a,sgolayfilt(data,1,71))) % plot filtered PID waveform (channel 2) and stimulus 1 
title(['Odor stimulus waveform (PID measurement)'])

% save([basepath,'/WCwaveform_' data(n).date,'_E',num2str(expnumber)],'data');