
%% Histogram of experiment start times

% load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;

% make unique trial ID
expstarttime = []
expstarttime_sec = []
ntrials = size(T,1)

for trial = 1:ntrials
    unqID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
    explist = excelfile(:,Excel.trialID);
    rownumber = find(strcmp(explist, unqID));
    t = excelfile{rownumber,Excel.starttime};
    % create contingency plan for non-duration values
    % if...
    tsec = (t*24*60*60);
    expstarttime_sec(trial) = tsec;
    t = seconds(t*24*60*60);
    t.Format = 'hh:mm:ss';
    expstarttime{trial} = t;
end

