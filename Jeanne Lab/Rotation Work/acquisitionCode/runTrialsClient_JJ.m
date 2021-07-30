expnumber = 1;
numtrials = 1;
trialduration = 8;
stimulustime = 4;
stimulusduration = 4;
VGain = 10;
IGain = 10;
data = AcquireTrace(expnumber,numtrials,trialduration,stimulustime,stimulusduration,'VC',VGain,IGain);


data = AcquireTrace(expnumber,numtrials,trialduration,stimulustime,stimulusduration,'IC',VGain,IGain);



data = AcquireTraceCamera_JJ(expnumber,numtrials,trialduration,stimulustime,stimulusduration,bursts)



% OptoStim.startTime = 4;
% OptoStim.duration = .5;
% OptoStim.duty = 1;
% OptoStim.amplitude = .04;
% OptoStim.amplitude = 0;
% data = AcquireTrace_Opto(1,1,8,4,2,OptoStim,'IC',10,10);


