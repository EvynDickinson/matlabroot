function [excelfile, Excel, xlFile] = load_QuadBowlExperiments
% [excelfile, Excel, xlFile] = load_FlyBowlExperiments
% Loads the excel file data into the structure 'raw' and finds the col headers
% Outputs:
% 'excelfile' [excel file: fly summary 2]
% 'Excel' [structure with locations of column headers]
% 'xlFile' excel file path
%
%
% ES Dickinson, University of Washington, Dec 2018
% updated Yale Universitty, Aug 2021

basePath = getCloudPath;
baseFolder = basePath(1:end-5); % path to just the '/jeanne lab/' folder

% [baseName, folder] = uigetfile('*xlsx', 'Select the Fly Summary 2 file');
% xlFile = fullfile(folder, baseName);

xlFile = [baseFolder 'Quad Bowl Experiments.xlsx'];

%load excel sheet data
[~,~,excelfile] = xlsread(xlFile,'Exp List');
% excelfile = readcell(xlFile,'Sheet','Exp List'); %update to new format --
% but need to adjust other functions that call the 1st two data type
% options

% Excel column for the headers:
Excel.headers = excelfile(1,:); %xltitles
    
Excel.date = find(strcmpi('Date',Excel.headers) == 1);
Excel.expID = find(strcmpi('Experiment ID', Excel.headers)==1);
Excel.protocol = find(strcmpi('Temp Protocol',Excel.headers) == 1);
Excel.arena = find(strcmpi('Arena',Excel.headers) == 1);
Excel.processed = find(strcmpi('Step 2',Excel.headers) == 1);
Excel.structure = find(strcmpi('Structure', Excel.headers)==1);
Excel.structurenum = find(strcmpi('Exp Num',Excel.headers) == 1);
Excel.backUp = find(strcmpi('Backed Up',Excel.headers) == 1);

Excel.genotype = find(strcmpi('Genotype',Excel.headers) == 1);
Excel.numflies = find(strcmpi('Num Flies',Excel.headers) == 1);
Excel.well_1 = find(strcmpi('Well 1', Excel.headers)==1);
Excel.well_2 = find(strcmpi('Well 2', Excel.headers)==1);
Excel.well_3 = find(strcmpi('Well 3', Excel.headers)==1);
Excel.well_4 = find(strcmpi('Well 4', Excel.headers)==1);
Excel.sex = find(strcmpi('Sex', Excel.headers)==1);
Excel.starved_hours = find(strcmpi('Starved Hours', Excel.headers)==1);

Excel.tracked = find(strcmpi('Tracked', Excel.headers)==1);
Excel.starttime = find(strcmpi('Start time', Excel.headers)==1);
Excel.facility = find(strcmpi('Facility', Excel.headers)==1);
Excel.trialID = find(strcmpi('Trial_ID', Excel.headers)==1);
Excel.step1 = find(strcmpi('Step 1', Excel.headers)==1);

Excel.experimenter = find(strcmpi('Experimenter', Excel.headers)==1);

end


















