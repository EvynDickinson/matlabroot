function [names, structNum, structInfo]  = getExcelStructureNames(~)
% [names, structNum]  = getExcelStructureNames
% 
% Pull unique structure names from the load_FlyBowlExperiments
% Input
%     Any input triggers a structure selection choices that returns
%     the structInfo variable
% Returns
%     names --> cell array with unique names
%     structNum --> number of unique structures
%     structInfo --> rows for each selected structure name
%
% ES Dickinson, Yale University, Aug 2021

    % Load current excel fly summary data sheet:
    [excelfile, Excel] = load_FlyBowlExperiments;

    % get list of unique structure names:
    structName_loc = cellfun(@ischar, excelfile(:, Excel.structure));
    structure_names = unique(excelfile(structName_loc, Excel.structure));
    structure_names(strcmpi(structure_names,' ')) = []; %remove blank spaces for section headers
    structure_names(strcmpi(structure_names,'Structure')) = []; %remove blank spaces for section headers
    names = structure_names;
    
    % get number of structures
    structNum = length(names);

    % select structures for more detailed information
    if nargin == 1 
       choiceIndex = listdlg('ListString', names, 'ListSize', [300, 400],...
                             'PromptString', 'Select data structure(s)', 'SelectionMode', 'multiple');
       for ii = 1:length(choiceIndex)
           loc = choiceIndex(ii);
           structInfo(ii).StructName = names{loc};
           struct_loc = strcmpi(excelfile(:, Excel.structure),names{loc});
           expNum_loc = ~cellfun(@isempty, excelfile(:, Excel.structurenum));
           structInfo(ii).rowLoc = struct_loc & expNum_loc;
           structInfo(ii).rowNum = find(struct_loc & expNum_loc);
           structInfo(ii).numTrials = sum(struct_loc & expNum_loc);
       end
    
    else
        structInfo = [];
    end
    
end