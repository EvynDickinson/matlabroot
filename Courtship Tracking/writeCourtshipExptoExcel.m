

function rows = writeCourtshipExptoExcel(params)
% Write the experiment data into Excel master sheet
% single line for each arena, since they can function independently or
% grouped
% ES Dickinson

% load current excel file:
[excelfile, Excel, xlFile] = load_CourtshipExperiments;
isExcelFileOpen(xlFile); % test that excel file is not open before writing to it

%            headers: {1�32 cell}
%               date: 2 **
%              expID: 3 **
%           protocol: 4 **
%           genotype: 5 **
%       experimenter: 6 **
%            trialID: 1 
%          starttime: 7
%      zeitgebertime: 8
%        compiledvid: 9
%            tracked: 10
%              step1: 11
%          processed: 12
%          structure: 13
%       structurenum: 14
%                sex: 15 **
%      starved_hours: 16
%           daynight: 17 **
%           numflies: 18 **
%             well_1: 19 **
%             well_2: 20 **
%             well_3: 21 **
%             well_4: 22 **
%         ramplength: 23 **
%              ramps: 24 **
%                ITI: 25 ** 
%     fragmentlength: 26 **
%                FPS: 27 **
%             backUp: 28
%           facility: 30

   

sheet = 'Exp List'; % excel sheet to write into       
fprintf('\n Writing to excel...')

% Newest row to populate:
nrow = size(excelfile,1)+1;

base_paramList = {'date', 'expID', 'protocol','genotype','experimenter',...
                  'sex','daynight','numflies','well_1', 'well_2', 'well_3', 'well_4',...
                  'ramplength', 'ramps', 'ITI', 'fragmentlength', 'FPS'};

% TODO update here to check that the excel sheet is writable ...

rows = [];
for ramp = 1:params.ramps
    rows(ramp) = nrow;
    
    for i = 1:length(base_paramList)
        xlswrite(xlFile, {params.(base_paramList{i})}, sheet, [Alphabet(Excel.(base_paramList{i})),num2str(nrow)]);
    end
    xlswrite(xlFile, {ramp}, sheet, [Alphabet(Excel.ramp),num2str(nrow)]);

    nrow = nrow+1;
end
fprintf('done \n')


end
