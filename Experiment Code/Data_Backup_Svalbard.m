% Data_BackUp_To_Server
clear; clc
 
% load excel file:   
[excelfile, Excel, XL] = load_QuadBowlExperiments;
loc = cellfun(@isnan,excelfile(2:end,Excel.processed));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.arena, Excel.expID,  Excel.backUp]);
loc1 = cellfun(@isnan,eligible_files(:,4));
c = cellfun(@string,eligible_files);
c(loc1,4) = ' ';
FileNames = join(c);
fileIdx = listdlg('ListString', FileNames,'ListSize',[300,450]);
%pull the list of dates and arenas to be 
List.date = eligible_files(fileIdx,1);
List.expID = eligible_files(fileIdx,3); 

% get base folder pathway
finishedFiles = [];
baseFolder = getCloudPath;
for ii = 1:length(fileIdx)
    if ~any(strcmp(finishedFiles,[List.date{ii} ' ' List.expID{ii}]))
        % COPY DATA TO SVALBARD
        from_folder = [baseFolder List.date{ii} '\'];
        to_folder = ['S:\Evyn\DATA\' List.date{ii} '\'];
        disp(['Backing up  ' List.date{ii}])
        copyfile(from_folder, to_folder) % Move folders to google drive:
    end
    finishedFiles{ii} = [List.date{ii} ' ' List.expID{ii}];
    % Write finished files to the excel master spreadsheet
    XLrow = rownums(fileIdx(ii));
    % write processed 'Y' into the excel sheet 

    try
        xlswrite(XL, {'T'}, 'Exp List', [Alphabet(Excel.backUp) num2str(XLrow)]);
    catch
        h = warndlg('Close Experiment Summary excel file and then close this warning box');
        uiwait(h)
        xlswrite(XL, {'T'}, 'Exp List', [Alphabet(Excel.backUp) num2str(XLrow)]);
    end

disp(['Finished ' FileNames(fileIdx(ii))])         
end

disp('Fullly uploaded');
