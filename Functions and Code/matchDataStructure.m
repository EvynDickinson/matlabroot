

function [rebuildFlag, addFlag, extradataFlag] = matchDataStructure(structFolder, StructureName, T)
% [rebuildFlag, addFlag, extradataFlag] = matchDataStructure(structFolder, StructureName, T)


    [rebuildFlag, addFlag, extradataFlag] = deal(false);
    
    
    % open the fileList and see if the Ns are matching
    checkList = load([structFolder StructureName '\fileList.mat']);
    checkTable = checkList.T(:,1:4); %single data
    testTable = T(:,1:4); %grouped data
    % Find rows unique to the grouped data or single data
    missingDataFiles = setdiff(checkTable, testTable, 'rows'); % files missing from the grouped data
    extraDataFiles = setdiff(testTable, checkTable, 'rows');   % extra files in the grouped but not single data list
    % if missing files, check that existing data structure in the single 
    % group actually has been updated and also includes those files:
    if ~isempty(missingDataFiles)
        checkList_2 = load([structFolder StructureName '\' StructureName ' post 3.1 data.mat'],'T');
        testTable = checkList_2.T(:,1:4);
        final_test = setdiff(testTable, checkTable, 'rows');   % extra files in the grouped but not single data  list
        if isempty(final_test) %if true the single structure has all the data and can be added as is
            % tag this structure for reloading into the grouped structure from the 3.1 data
            addFlag = true;
        else % this structure needs to be rebuilt from the 3.1 data processing step to include any new data files
            rebuildFlag = true; 
        end
    end
    % extra data in the grouped structure that requires a rebuild and manual eyes on the problem
    if ~isempty(extraDataFiles)
        extradataFlag = true;
        rebuildFlag = true;
    end


end