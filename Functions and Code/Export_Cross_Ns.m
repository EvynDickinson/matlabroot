
clearvars; close all; clc 

% Load current excel fly summary data sheet:
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
sheet = 'Structure List';

[structure_names, nstructs, sInfo]  = getExcelStructureNames(true);

for ii = 1:nstructs
   % make counts
   nExp = sum(sInfo(ii).numTrials);
   nflies = sum(cell2mat(excelfile(sInfo(ii).rowLoc, Excel.numflies)));
   
   % add info to growing matrix
   outputMat{ii+1,1} = sInfo(ii).StructName;
   outputMat{ii+1,2} = num2str(nExp);
   outputMat{ii+1,3} = num2str(nflies);
end

outputMat{1,1} = 'Structure:';
outputMat{1,2} = 'Num Exp:';
outputMat{1,3} = 'Num Flies:';

% write numbers to Excel Sheet
writecell(outputMat, xlFile, 'Sheet', sheet);

%% TUTHILL LAB DATA STRUCTURE VERSION:





% 
% %% Define variables:
% 
% % minimum line in excel file to include in N value search
% strtpoint = 470; 
% % Save output name and file path:
% xlFile = 'C:\matlabroot\Fly Structure Numbers.xlsx'; % Path and Nam
% 
% %% Extract N values from Excel Sheet
% 
% % Structure List --> sheet name
% 
% % Pull data from fly summary:    
% [excelfile, Excel] = load_flysummary; 
% 
% % Pull out the structure names from the file:
% structure_names.excelfile = excelfile(strtpoint:end, Excel.new_struct_name);
% structure_names.numberslist = excelfile(strtpoint:end, Excel.structurenum);
% ind = 1;
% for ii = 1:length(structure_names.excelfile)
%     if ischar(structure_names.excelfile{ii})
%         structure_names.text{ind} = (structure_names.excelfile{ii});
%         structure_names.num(ind) = (structure_names.numberslist{ii});
%         ind = ind+1;
%     end
% end; clear ind
%   
% %find the unique structure names:
% [structure_names.unique, ia] = unique(structure_names.text);
% for ii = 1:length(ia)
%    crossName = structure_names.unique{ii};
%    Loc = find(strcmpi(crossName, structure_names.text) == 1);  
%    structure_names.finalnumbers{ii} = num2str(max(structure_names.num(Loc)));   
% end
% 
% % Build the output structure:
% output_Matrix{1,1} = 'Structures:';
% output_Matrix{1,2} = 'Fly Line Totals:';
% output_Matrix(2:(length(ia)+1),1) = structure_names.unique;
% output_Matrix(2:(length(ia)+1),2) = structure_names.finalnumbers;
% 
% %% Create excel sheet with fly line totals: 
% 
% writecell(output_Matrix, xlFile, 'Sheet', 'Structure Numbers')
% 
% 
% 
