
function writeExptoExcel(params)



% load current excel file:
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
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

paramList = {'date', 'expID', 'genotype', 'protocol',...
             'well_1', 'well_2', 'well_3', 'well_4',...
             'PF_Batch', 'YF_Batch'};
  
for ii = 1:length(paramList)
    xlswrite(xlFile, {params.(paramList{ii})}, sheet, [Alphabet(Excel.(paramList{ii})),num2str(nrow)]);
end
fprintf('done \n')

end