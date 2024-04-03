

function [rebuildFlag, addFlag, extradataFlag] = matchDataStructure(StructureName, test_T)
% [rebuildFlag, addFlag, extradataFlag] = matchDataStructure(structFolder, StructureName, T)
%
% rebuildFlag      -->  rerun Step3.1 to update the data structure for current data
% addFlag           --> data is missing from test structure but is present in current post3.1
%                                so load the existing data struct into the grouped data struct
% extradataFlag --> there are more trials in the test data structure than exist in
%                                 the 3.1 data structure and should be looked at 
% ESD 04.2024 Yale

    baseFolder = getCloudPath;
    structFolder = [baseFolder 'Data structures\'];

    [rebuildFlag, addFlag, extradataFlag] = deal(false);
    C_sel = 1:4; %selected columns to compare


% 1) Does the flielist in the structure folder match the 3.1 processed data file?
    T_truth = load([structFolder StructureName '\fileList.mat'],'T'); %loads 'T' which is the most up-to-date flieList
    Tpost3_1 = load([structFolder StructureName '\' StructureName ' post 3.1 data.mat'],'T');
    % Find rows unique to the grouped data or single data
    missingDataFiles = setdiff(T_truth.T(:,C_sel), Tpost3_1.T(:,C_sel), 'rows'); % files missing from the grouped data
    extraDataFiles = setdiff(Tpost3_1.T(:,C_sel), T_truth.T(:,C_sel), 'rows');   % extra files in the grouped but not single data list
    if ~isempty(missingDataFiles) || ~isempty(extraDataFiles)
        rebuildFlag = true;
    end

% 2) Does the test data structure match the ground truth?
    if nargin == 2
        if isempty(test_T)
            addFlag = true;
        else
            missingDataFiles = setdiff(T_truth.T(:,C_sel), test_T(:,C_sel), 'rows'); % files missing from the grouped data
            extraDataFiles = setdiff(test_T(:,C_sel), T_truth.T(:,C_sel), 'rows');   % extra files in the grouped but not single data list
             if ~isempty(missingDataFiles)
                addFlag = true;
             end
             if ~isempty(extraDataFiles) %more data in old group 
                extradataFlag = true;
             end

        end
    end

end