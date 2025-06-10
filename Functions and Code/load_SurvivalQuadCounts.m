
function [excelfile, Excel, xlFile] = load_SurvivalQuadCounts

basePath = getCloudPath;
baseFolder = basePath(1:end-5); % path to just the '/jeanne lab/' folder

xlFile = [baseFolder 'Survival Quad Bowl Counts.xlsx'];

% Load excel sheet data
[~,~,excelfile] = xlsread(xlFile,'Death counts');

% Excel row for the headers:
Excel.headers = excelfile(1,:); % xltitles

Excel.trialID = find(strcmpi('Trial ID',Excel.headers) == 1);
Excel.date = find(strcmpi('Date', Excel.headers)==1);
Excel.genotype = find(strcmpi('Genotype',Excel.headers) == 1);
Excel.expID = find(strcmpi('Exp ID', Excel.headers)==1);
Excel.numflies = find(strcmpi('Num',Excel.headers) == 1);
Excel.arena = find(strcmpi('Arena',Excel.headers) == 1);
Excel.videos = excelfile(1, 5:45);
Excel.counted = find(strcmpi('Counted',Excel.headers) == 1);
Excel.analyzed = find(strcmpi('Analyzed',Excel.headers) == 1);

end