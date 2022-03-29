


%% -------- Find the files that haven't been analyized yet and run them -------------
clear; close all; clc
autoSave = true;
essentialfigs = true; 
excelWrite = true;

%load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;
loc = cellfun(@isnan,excelfile(2:end,Excel.numflies));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.arena, Excel.expID, Excel.processed]);
loc1 = cellfun(@isnan,eligible_files(:,4));
c = cellfun(@string,eligible_files);
c(loc1,4) = ' ';
FileNames = join(c);
fileIdx = listdlg('ListString', FileNames,'ListSize',[250,450]);
%pull the list of dates and arenas to be 
List.date = eligible_files(fileIdx,1);
List.expID = eligible_files(fileIdx,3); 

%get base folder pathway
finishedFiles = [];
baseFolder = getCloudPath;
for ii = 1:length(fileIdx)
    inputPath = [baseFolder List.date{ii} '/Analysis/' List.expID{ii} ' preformed data.mat'];
    if ~any(strcmp(finishedFiles,[List.date{ii} ' ' List.expID{ii}]))
        results = runQuadStep2_1(inputPath,autoSave,essentialfigs); % Run the basic figures
    end
    finishedFiles{ii} = [List.date{ii} ' ' List.expID{ii}];
    if excelWrite == true
        if strcmpi(results, 'Saved Data')
            XLrow = rownums(fileIdx(ii));
            % write number of flies into the excel sheet
            try
                xlswrite(XL, {'Y'}, 'Exp List', [Alphabet(Excel.processed) num2str(XLrow)]);
            catch
                h = warndlg('Close Experiment Summary excel file and then close this warning box');
                uiwait(h)
                xlswrite(XL, {'Y'}, 'Exp List', [Alphabet(Excel.processed) num2str(XLrow)]);
            end
        end
    end
disp(['Finished ' FileNames(fileIdx(ii))])         
end




   








