

% data groupings: 
% 10.29.2021 -- 2.8.2022 : date>arena>analysis (:)   [black background video copy] 
% these are before line 149 in the folder



%%  Copy data from post 2.8.22
clear
clc

baseFolder = getCloudPath;
trialDataBase = getBasePath;

%load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;
trial_ID_list = excelfile(:,Excel.trialID);
rownums = [];
% rownums = 2:1550;
for i = 2:length(trial_ID_list)
    if ~ischar(trial_ID_list{i})
        rownums = [rownums; i];
    end
end

OG_Orientation = datetime('02.20.2022','InputFormat','MM.dd.yyyy');

for i = 1: length(rownums)

    [success, success_1, success_2] = deal(true);

    % trial ID structure: [date_expname_arena]
    trial_date = excelfile{rownums(i),Excel.date};
    trial_expName = excelfile{rownums(i),Excel.expID};
    trial_arena = excelfile{rownums(i),Excel.arena};
    
    % generate a unique trial ID:
    trialID = [trial_date '_' trial_expName '_' trial_arena];

    % check the date
    try currTime = datetime(trial_date,'InputFormat','MM.dd.yyyy');
    catch 
        continue
    end
    if  currTime < OG_Orientation %new data organization
        % disp(['Skipped: ' trialID])
        continue
    end

    % make a new folder if it doesn't exist for the new items: 
    startDir = [baseFolder trial_date '/Arena '  trial_arena '/'];
    endDir = [trialDataBase trialID '\'];
    if ~exist(endDir, 'dir')
        mkdir(endDir)
    end

    % check for more than 1 experiment (which would require subdividing the data files individually)
    tempDir = dir([startDir '*timecourse data.mat']);
    if size(tempDir,1)>1 % more than 1 experiment for the day
        % split the data by the trial:
        
        % timecourse data:
        try
        success_1 = copyfile([startDir trial_expName trial_arena ' timecourse data.mat'],endDir);
        catch
            warndlg(['No timecourse data found for trial: ' startDir trial_expName trial_arena])
            disp(i)
            break
        end
        if ~success_1
            warndlg(['Copy error on ' startDir trial_expName trial_arena])
        end

        % speed data:
        try
        success_2 = copyfile([startDir trial_expName ' speed data.mat'],endDir);
        catch
            warndlg(['No speed data found for trial: ' startDir trial_expName trial_arena])
            disp(i)
            break
        end
         if ~success_2
            warndlg(['Copy error on ' startDir trial_expName trial_arena])
         end

        % sleep data
        try
        copyfile([startDir trial_expName ' sleeping data.mat'],endDir);
        disp(['Copied sleep data for: ' startDir trial_expName])
        catch
            disp(['no sleeping data found for: ' trial_expName])
            % break
        end

        % figures:
        copyfile([startDir '*' trial_expName ' *.pdf'],endDir)
        copyfile([startDir '*' trial_expName ' *.png'],endDir)
 
    else % COPY WHOLE FOLDER, ONLY 1 TRIAL CONTAINED

        success = copyfile(startDir, endDir);
        if ~success
            disp(['Error on copying full directory: ' startDir])
        end
    end
    
    % update the excel sheet that the file has been copied
    if all([success, success_1, success_2])
        writecell({trialID},XL,'Sheet','Exp List','Range',[Alphabet(Excel.trialID) num2str(rownums(i))]);
    end
    disp(['Finished ' trialID])
end


%%  Copy data from PRE 2.8.22
clear
clc

baseFolder = getCloudPath;
trialDataBase = getBasePath;

%load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;
trial_ID_list = excelfile(:,Excel.trialID);
rownums = [];
% rownums = 2:1550;
for i = 2:length(trial_ID_list)
    if ~ischar(trial_ID_list{i})
        rownums = [rownums; i];
    end
end

OG_Orientation = datetime('02.20.2022','InputFormat','MM.dd.yyyy');

for i = 1: length(rownums)

    [success, success_1, success_2] = deal(true);

    % trial ID structure: [date_expname_arena]
    trial_date = excelfile{rownums(i),Excel.date};
    trial_expName = excelfile{rownums(i),Excel.expID};
    trial_arena = excelfile{rownums(i),Excel.arena};
    
    % generate a unique trial ID:
    trialID = [trial_date '_' trial_expName '_' trial_arena];

    % check the date
    try currTime = datetime(trial_date,'InputFormat','MM.dd.yyyy');
    catch 
        continue
    end
    if  currTime > OG_Orientation %new data organization
        % disp(['Skipped: ' trialID])
        continue
    end

    % make a new folder if it doesn't exist for the new items: 
    startDir = [baseFolder trial_date '/Arena '  trial_arena '/analysis/'];
    endDir = [trialDataBase trialID '\'];
    if ~exist(endDir, 'dir')
        mkdir(endDir)
    end

    % check for more than 1 experiment (which would require subdividing the data files individually)
    tempDir = dir([startDir '*timecourse data.mat']);
    if size(tempDir,1)>1 % more than 1 experiment for the day
        % split the data by the trial:
        
        % timecourse data:
        try
        success_1 = copyfile([startDir trial_expName trial_arena ' timecourse data.mat'],endDir);
        catch
            warndlg(['No timecourse data found for trial: ' startDir trial_expName trial_arena])
            disp(i)
            break
        end
        if ~success_1
            warndlg(['Copy error on ' startDir trial_expName trial_arena])
        end

        % preformed data:
        try
        copyfile([startDir trial_expName ' preformed data.mat'],endDir);
        catch
            disp(['No preformed data found for trial: ' startDir trial_expName trial_arena])
        end

        % falsepoints data
        try
        copyfile([startDir trial_expName ' FalsePoints.mat'],endDir);
        disp(['Copied falsepoints data for: ' startDir trial_expName])
        catch
            disp(['no falsepoints data found for: ' trial_expName])
            % break
        end

        % arena mask data
        try
        copyfile([startDir trial_expName trial_arena ' ArenaMask.mat'],endDir);
        disp(['Copied falsepoints data for: ' startDir trial_expName])
        catch
            disp(['no arena mask data found for: ' trial_expName])
            % break
        end

        % figures: 
        try copyfile([startDir '*' trial_expName ' *.pdf'],endDir)
        catch
             copyfile([startDir '*' trial_expName trial_arena ' *.pdf'],endDir)
        end
        try copyfile([startDir '*' trial_expName ' *.png'],endDir)
        catch
            copyfile([startDir '*' trial_expName trial_arena ' *.png'],endDir)
        end
 
    else % COPY WHOLE FOLDER, ONLY 1 TRIAL CONTAINED

        success = copyfile(startDir, endDir);
        if ~success
            disp(['Error on copying full directory: ' startDir])
        end
    end
    
    % update the excel sheet that the file has been copied
    if all([success, success_1, success_2])
        writecell({trialID},XL,'Sheet','Exp List','Range',[Alphabet(Excel.trialID) num2str(rownums(i))]);
    end
    disp(['Finished ' trialID])
end


