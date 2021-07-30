% Data acquisition script for both odor and shock aversive learning


% %  DAQ_gui
% param.fly_num = '5.1';
% FileList = dir([dir_root '\*.tif']);
% session = length(FileList)+1;
% AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
%                           param.ODOR_start,param.SHOCK_start,...
%                           param.ODOR_dur,param.SHOCK_dur);
% % write information to excel on trial
% write_OdorShockTrial_excel(param, 'ScentBlast', '1%', '2-butanone', num2str(session));
% % write_OdorShockTrial_excel(param, 'ScentBlast', '1%', '3-octanol', num2str(session));
% % write_OdorShockTrial_excel(param, 'ScentBlast', '0%', 'empty', num2str(session));


%% Parameters
clear 

param = load_learningParams;
param.fly_num = '5';
num.emptyTrials = 1;    %empty bottle trial for start and finish
num.solventTrials = 2;  %number of solvent trials
num.odors = 2;          %number of odors outside to test *does not include paraffin oil
num.baseTrials = 5;     %number of pre-training session for each odor
num.trainingTrials = 3; %number of conditioning trials
num.testTrials = 5;     %number of test trials post conditioning
param.sex = 'F';
param.folder_name = param.date;
param.tag = 'mpImage';
dir_root = ['E:\Evyn\' param.date];
if ~exist(dir_root,'dir')
    mkdir(dir_root)
end
h = warndlg('Confirm camera save folder');
uiwait(h)

% configuration parameters: 
param.ODOR_start = 10;     % odor application start time (s)
param.ODOR_dur = 1;        % odur duration (s)
param.SHOCK_start = 10.5;  % shock trigger start time (s)
param.SHOCK_dur = 2;       % shock duration (s)
param.postStim_dur = 15;   % time after stim off to keep recording
stim_end = max([(param.ODOR_dur+param.ODOR_start),(param.SHOCK_dur+param.SHOCK_start)]);

param.trial_duration = ceil(stim_end + param.postStim_dur);
param.postStim_dur = param.trial_duration-stim_end;

disp(['Goal Num frames: ' num2str(floor((param.trial_duration-2)/3*100-1))])
disp(['Video record time: ' num2str(param.trial_duration-2) ' sec'])

h = warndlg({'Check:'; '1) shock OFF'; '2) Frame count matches'; '3) All Excel files are closed'});
uiwait(h)

%% Baseline acquisition
disp('Begin CONTROL phase')
% ------- Empty bottle control -------
h = warndlg('Load empty bottle');
uiwait(h)
% Confirm that the session number & saving number match:
FileList = dir([dir_root '\*.tif']);
session = length(FileList)+1;
param.start_session = session;
disp(['Video name: ' param.tag ' ' num2str(session)]);
disp('Loaded empty bottle')
answer = questdlg('Loaded empty bottle and camera ready?');
switch answer
    case 'Cancel'
       return
end
data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                          param.ODOR_start,param.SHOCK_start,...
                          param.ODOR_dur,param.SHOCK_dur);
% write information to excel on trial
write_OdorShockTrial_excel(param, 'Empty', param.solvent_conc, 'empty', num2str(session));


% ------- Solvent control -------
h = warndlg(['Load odor: ' param.solvent]);
uiwait(h)
for n = 1:num.solventTrials
    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    session = session+1;
    param.start_session = session;

    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.solvent])
    answer = questdlg(['Loaded ' param.solvent ' and camera ready?']);
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'Control', param.solvent_conc, param.solvent, num2str(session));
end

% ------- Conditioned Stimulus -------
h = warndlg(['Load odor: ' param.csPOS_odor]);
uiwait(h)
for n = 1:num.baseTrials
    session = session+1;
    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    if ~((length(FileList)+1)==session)
        warndlg('Video number and session don''t match')
    end
    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.csPOS_odor])
    answer = questdlg({['Loaded ' param.csPOS_odor ' and camera ready?'];...
                       ['Video ' num2str(n) '/' num2str(num.baseTrials)]});
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'Control', param.csPOS_odor_conc, param.csPOS_odor, num2str(session));
end

% ------- Conditioned Stimulus -------
h = warndlg(['Load odor: ' param.csNEG_odor]);
uiwait(h)
for n = 1:num.baseTrials
    session = session+1;
    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    if ~((length(FileList)+1)==session)
        warndlg('Video number and session don''t match')
    end
    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.csNEG_odor])
    answer = questdlg({['Loaded ' param.csNEG_odor ' and camera ready?'];...
                       ['Video ' num2str(n) '/' num2str(num.baseTrials)]});
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'Control', param.csNEG_odor_conc, param.csNEG_odor, num2str(session));
end

%% TRAINING 
disp('Begin TRAINING phase')
h = warndlg('TURN ON SHOCK AM-SYSTEMS');
uiwait(h)
% ------- Conditioned Stimulus -------
h = warndlg(['Load odor: ' param.csPOS_odor]);
uiwait(h)
for n = 1:num.trainingTrials
    session = session+1;
    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    if ~((length(FileList)+1)==session)
        warndlg('Video number and session don''t match')
    end
    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.csPOS_odor])
    answer = questdlg({['Loaded ' param.csPOS_odor ' and camera ready?'];...
                       ['Video ' num2str(n) '/' num2str(num.trainingTrials)]});
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'CS+ training', param.csPOS_odor_conc, param.csPOS_odor, num2str(session));
end
h = warndlg('TURN OFF SHOCK AM-SYSTEMS');
uiwait(h)


% ------- Unconditioned Stimulus -------
h = warndlg(['Load odor: ' param.csNEG_odor]);
uiwait(h)
for n = 1:num.trainingTrials
    session = session+1;
    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    if ~((length(FileList)+1)==session)
        warndlg('Video number and session don''t match')
    end
    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.csNEG_odor])
    answer = questdlg({['Loaded ' param.csNEG_odor ' and camera ready?'];...
                        ['Video ' num2str(n) '/' num2str(num.trainingTrials)]});
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'CS- training', param.csNEG_odor_conc, param.csNEG_odor, num2str(session));
end

%% TEST
disp('Begin TESTING phase')

% ------- Conditioned Stimulus -------
h = warndlg(['Load odor: ' param.csPOS_odor]);
uiwait(h)
for n = 1:num.testTrials
    session = session+1;
    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    if ~((length(FileList)+1)==session)
        warndlg('Video number and session don''t match')
    end
    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.csPOS_odor])
    answer = questdlg({['Loaded ' param.csPOS_odor ' and camera ready?'];...
                       ['Video ' num2str(n) '/' num2str(num.testTrials)]});
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'Test', param.csPOS_odor_conc, param.csPOS_odor, num2str(session));
end

% ------- Conditioned Stimulus -------
h = warndlg(['Load odor: ' param.csNEG_odor]);
uiwait(h)
for n = 1:num.testTrials
    session = session+1;
    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    if ~((length(FileList)+1)==session)
        warndlg('Video number and session don''t match')
    end
    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.csNEG_odor])
    answer = questdlg({['Loaded ' param.csNEG_odor ' and camera ready?'];
                       ['Video ' num2str(n) '/' num2str(num.testTrials)]});
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'Test', param.csNEG_odor_conc, param.csNEG_odor, num2str(session));
end

% ------- Solvent -------
h = warndlg(['Load odor: ' param.solvent]);
uiwait(h)
for n = 1:num.solventTrials
    session = session + 1;
    param.final_session = session;

    % Confirm that the session number & saving number match:
    FileList = dir([dir_root '\*.tif']);
    if ~((length(FileList)+1)==session)
        warndlg('Video number and session don''t match')
    end
    disp(['Video name: ' param.tag ' ' num2str(session)]);
    disp(['Loaded odor: ' param.solvent])
    answer = questdlg(['Loaded ' param.solvent ' and camera ready?']);
    switch answer
        case 'Cancel'
           return
    end
    data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                              param.ODOR_start,param.SHOCK_start,...
                              param.ODOR_dur,param.SHOCK_dur);
    % write information to excel on trial
    write_OdorShockTrial_excel(param, 'Test', param.solvent_conc, param.solvent, num2str(session));
end

% ------- Empty bottle control -------
h = warndlg('Load empty bottle');
uiwait(h)
% Confirm that the session number & saving number match:
FileList = dir([dir_root '\*.tif']);
session = session+1;
param.start_session = session;
disp(['Video name: ' param.tag ' ' num2str(session)]);
disp('Loaded empty bottle')
answer = questdlg('Loaded empty bottle and camera ready?');
switch answer
    case 'Cancel'
       return
end
data(session) = AcquireTrace_shock(param.fly_num,1,param.trial_duration,...
                          param.ODOR_start,param.SHOCK_start,...
                          param.ODOR_dur,param.SHOCK_dur);
% write information to excel on trial
write_OdorShockTrial_excel(param, 'Empty', param.solvent_conc, 'empty', num2str(session));

disp('DONE!')

% Save data to folder:

flyID = generate_flyID(param.date, param.fly_num);
save([dir_root '\' flyID '_data'])

disp('Data saved!')















