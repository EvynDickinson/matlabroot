function rows = writeQuadExptoExcel(params)
% Write the experiment data into Excel master sheet
% single line for each arena, since they can function independently or
% grouped
% ES Dickinson

% load current excel file:
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
isExcelFileOpen(xlFile); % test that excel file is not open before writing to it

sheet = 'Exp List'; % excel sheet to write into       
fprintf('\n Writing to excel...')

% Newest row to populate:
nrow = size(excelfile,1)+1;

% params.well_1 
% params.well_2
% params.well_3 
% params.well_4 
% params.date 
% params.expID
% params.genotype
% params.protocol 
% params.num  
% params.arena
% params.(arena).sex
% params.(arena).starved_hours
% params.experimenter
% params.day_night

base_paramList = {'date', 'expID', 'protocol'};
arena_paramList = {'genotype','well_1', 'well_2', 'well_3', 'well_4', 'sex','starved_hours'};

rows = [];
for arena = 1:4
    rows(arena) = nrow;
    % experimenter:
    xlswrite(xlFile, {params.experimenter}, sheet, [Alphabet(Excel.experimenter),num2str(nrow)]);
    % day or night incubator:
    xlswrite(xlFile, {params.day_night}, sheet, [Alphabet(Excel.daynight),num2str(nrow)]);
    % arena:
    xlswrite(xlFile, {Alphabet(arena)}, sheet, [Alphabet(Excel.arena),num2str(nrow)]);
    % trial name:
    trial_name = [params.date '_' params.expID '_' Alphabet(arena)];
    xlswrite(xlFile, {trial_name}, sheet, [Alphabet(Excel.trialID),num2str(nrow)]);
    % shared data:
    for ii = 1:length(base_paramList)
        xlswrite(xlFile, {params.(base_paramList{ii})}, sheet, [Alphabet(Excel.(base_paramList{ii})),num2str(nrow)]);
    end
    % arena specific data:
    for ii = 1:length(arena_paramList)
        xlswrite(xlFile, {params.(['Arena' Alphabet(arena)]).(arena_paramList{ii})},...
                 sheet, [Alphabet(Excel.(arena_paramList{ii})),num2str(nrow)]);
    end
    nrow = nrow+1;
end
fprintf('done \n')


end
