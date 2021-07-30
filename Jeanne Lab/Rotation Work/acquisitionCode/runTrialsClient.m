%% Former code from Jamie

expnumber = 1;
numtrials = 1;
trialduration = 8;
stimulustime = 4;
stimulusduration = 0;
VGain = 10;
IGain = 10;
data = AcquireTrace(expnumber,numtrials,trialduration,stimulustime,stimulusduration,'VC',VGain,IGain);

data = AcquireTrace(expnumber,numtrials,trialduration,stimulustime,stimulusduration,'IC',VGain,IGain);



% OptoStim.startTime = 4;
% OptoStim.duration = .5;
% OptoStim.duty = 1;
% OptoStim.amplitude = .04;
% OptoStim.amplitude = 0;
% data = AcquireTrace_Opto(1,1,8,4,2,OptoStim,'IC',10,10);



%% Run odor and blue LED for GCamp6s
clear
% save parameters:

% configuration parameters: 
expnumber = 1;       % Fly number
numtrials = 1;       % Num of repetitions
ODOR_start = 10;     % odor application start time (s)
ODOR_dur = 1;        % odur duration (s)
SHOCK_start = 10.5;  % shock trigger start time (s)
SHOCK_dur = 2;       % shock duration (s)
stim_length = 2;     % length (s) of shock&odor stimulus
postStim_dur = 20;   % time after stim off to keep recording

trial_duration = ODOR_start+stim_length+postStim_dur;
disp(['Num frames: ' num2str((trial_duration-2)/3*100)])
disp(['Video record time: ' num2str(trial_duration-2)])

% could add a visual cue here of the timing

%% Run experiment
data = AcquireTrace_shock(expnumber,numtrials,trial_duration,ODOR_start,SHOCK_start,ODOR_dur,SHOCK_dur);
disp(['Done exp ' num2str(expnumber)])






%
% % Add comment into 'Notes' section in Excel
% answer = questdlg('Add comment on fly behavior?', 'Notes', 'Yes', 'No', 'No');
% a = strcmpi('Notes:', Excel.headers); %column number
% NoteRange = [Alphabet(a) num2str(end_pos)]; %cell address
% switch answer
% case 'Yes'
%     newNote = inputdlg('Note for Excel');
%     Notes = {[param.body ', ' param.loading ', ' newNote]};
%     xlswrite(xlFile, Notes, 'Sheet1', NoteRange);
%     fprintf('\n Note added to Excel file \n')
% case 'No'
%     Notes = {[param.body ', ' param.loading]};
%     xlswrite(xlFile, Notes, 'Sheet1', NoteRange);
%     fprintf('\n No comments on behavior \n')
% end



