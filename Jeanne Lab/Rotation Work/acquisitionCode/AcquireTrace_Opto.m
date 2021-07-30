function data = AcquireTrace_Opto(expnumber,numtrials,trialduration,stimulustime,stimulusduration,OptoStim,mode,VGain,IGain)
%expnumber = 1;
%trialduration = 6; %in seconds;
%numtrials = 1;
%stimulustime = 4; %in seconds;
%stimulusduration = 1; %in seconds;
%mode = 'IC';
%VGain = 10;
%IGain = 10;


%Set date appropriately:
todaysDate = datestr(datetime('now'), 'yyyy-mm-dd');
if ~isdir(['E:\Jamie\data\' todaysDate])
    mkdir(['E:\Jamie\data\' todaysDate]);
end
basepath = ['E:\Jamie\data\' todaysDate];

D = dir([basepath '\WCwaveform_',todaysDate,'_E',num2str(expnumber),'.mat']);
if isempty(D)
    % if no saved data exists then this is the first trial
    startingtrial=1;
else
    %load current data file
    load(['E:\Jamie\Data\' todaysDate,'\WCwaveform_',todaysDate,'_E',num2str(expnumber),'.mat']','data');
    startingtrial = length(data)+1;
end



%% set trial parameters
for k = 1:numtrials
    n = startingtrial-1+k;
    data(n).date = todaysDate;                                 % experiment date
    data(n).expnumber = expnumber;                          % experiment number
    data(n).trial = n;                                        % trial number
    % data(n).sampleTime = clock;
    data(n).acquisition_filename = mfilename('fullpath');    %saves name of mfile that generated data
    % sampling rates
    data(n).sampratein = 10000;                              % input sample rate
    data(n).samprateout = 10000;                           % output sample rate becomes input rate as well when both input and output present
    data(n).trialduration = trialduration;                            % trial duration
    data(n).stimulustime = stimulustime;
    data(n).variableGain1 = NaN;                             %Amplifier 1 alpha
    data(n).ImGain1 = IGain;
    data(n).VmGain1 = VGain;
    data(n).mode = mode;
    
    data(n).testVoltage = 0/20; %10mV/(10mV/V) in mV (for voltage clamp)
    data(n).testCurrent = -5/400; %10pA/(400pA/V) in pA (for current clamp)
    
    %make column vector for use as master8 trigger
    trailingdurationSTIM = trialduration - stimulusduration - stimulustime;
    trailingdurationOPTOSTIM = trialduration - OptoStim.startTime - OptoStim.duration;
    trailingdurationCOMMAND = trialduration - 2;
    
    data(n).stimulusVector = [zeros(1,data(n).samprateout*stimulustime) 5*ones(1,data(n).samprateout*stimulusduration) zeros(1,data(n).samprateout*trailingdurationSTIM)]';
    data(n).OptoStimVector = [zeros(1,data(n).samprateout*OptoStim.startTime) OptoStim.amplitude*ones(1,data(n).samprateout*OptoStim.duration) zeros(1,data(n).samprateout*trailingdurationOPTOSTIM)]';

    if strcmp(data(n).mode, 'VC')
        data(n).commandVector = [zeros(1,data(n).samprateout*1) data(n).testVoltage*ones(1,data(n).samprateout*1) zeros(1,data(n).samprateout*trailingdurationCOMMAND)]';
    else
        data(n).commandVector = [zeros(1,data(n).samprateout*1) data(n).testCurrent*ones(1,data(n).samprateout*1) zeros(1,data(n).samprateout*trailingdurationCOMMAND)]';
    end
    
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
    s.addAnalogOutputChannel('Dev1', [0:2] , 'Voltage');
    s.Rate = data(n).samprateout;
    
    s.queueOutputData([data(n).commandVector data(n).stimulusVector data(n).OptoStimVector]);
    x = s.startForeground();
    
    if strcmp(data(n).mode, 'VC')
        %for voltage clamp...
        %This is 10Vm
        voltage1 = x(:,2)/(data(n).VmGain1*10)*1000;  %+data(n).variableOffset1;
        %This is .5Im
        current1 = x(:,1)/(data(n).ImGain1*.5)*1000;  %+data(n).ImOffset1;
    else
        %for current clamp
        voltage1 = x(:,1)/(data(n).VmGain1*10)*1000;  %+data(n).variableOffset1;
        current1 = x(:,2)/(data(n).ImGain1*.5)*1000;  %+data(n).ImOffset1;
    end
    
    figure(1)
    if strcmp(data(n).mode, 'VC')
        subplot(8,1,[1 2 3])
        plot(timevec,current1)
        title(['Experiment: ' todaysDate '_E',num2str(expnumber) '. Trial number ' num2str(n)],'Interpreter', 'none')
        ylabel('current (pA)')
        set(gca, 'XTickLabel', [])
        box off
        
        subplot(8,1,[4:6])
        plot(timevec,voltage1)
        ylabel('voltage (pA)')
        box off
        set(gca, 'XTickLabel', [])
        
        subplot(8,1,7)
        plot(timevec,data(n).stimulusVector)
        box off
        set(gca, 'XTickLabel', [])
        ylabel({'stimulus','voltage (mV)'})
        
        subplot(8,1,8)
        plot(timevec,data(n).commandVector)
        box off
        ylabel({'command','voltage (mV)'})
        xlabel('Time (s)')
        set(gcf, 'Position', [20 210 1430 760])
        
    else
        subplot(8,1,[1 2 3])
        plot(timevec,voltage1)
        title(['Experiment: ' todaysDate '_E',num2str(expnumber) '. Trial number ' num2str(n)],'Interpreter', 'none')
        ylabel('voltage (mV)')
        set(gca, 'XTickLabel', [])
        box off
        
        subplot(8,1,[4:6])
        plot(timevec,current1)
        ylabel('current (pA)')
        box off
        set(gca, 'XTickLabel', [])
        
        subplot(8,1,7)
        plot(timevec,data(n).stimulusVector)
        box off
        set(gca, 'XTickLabel', [])
        ylabel({'stimulus','voltage (mV)'})
        
        subplot(8,1,8)
        plot(timevec,data(n).commandVector)
        box off
        ylabel({'command','voltage (mV)'})
        xlabel('Time (s)')
        set(gcf, 'Position', [20 210 1430 760])        
    end
    
    %Measure whole-cell parameters (current clamp only):
    if strcmp(data(n).mode, 'IC')
        data(n).VRest = mean(voltage1(1:9999));
        
        voltageDifference = (mean(voltage1(15000:19999)) - mean(voltage1(1:9999)))*1e-3; %in Volts.
        currentDifference = (mean(current1(15000:19999)) - mean(current1(1:9999)))*1e-12; %in Amps.
        data(n).Rin = (voltageDifference/currentDifference)/1e6; %In MOhms
    end
    
    figure(2)
        for i = 1:n
            if length(data(i).Rin) == 0;
                allRin(i) = NaN;
                allVrest(i) = NaN;
            else
                allRin(i) = data(i).Rin;
                allVrest(i) = data(i).VRest;
            end
        end
        subplot(2,1,1)
        plot(1:n, allRin, '.')
        ylabel('R_{in} (M\Omega)')
        subplot(2,1,2)
        plot(1:n,allVrest,'.')
        ylabel('V_{rest} (mV)')
        xlabel('Trial num.')
        set(gcf, 'Position', [1460 700 440 270])
    
    
    %% save data(n)
    save([basepath,'/Raw_WCwaveform_' data(n).date,'_E',num2str(expnumber),'_',num2str(n)],'current1','voltage1');
end

save([basepath,'/WCwaveform_' data(n).date,'_E',num2str(expnumber)],'data');