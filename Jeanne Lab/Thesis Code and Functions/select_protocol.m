
function protocol = select_protocol
% protocol = select_protocol
% 
% Select the protocol from the list of possibilties
% contained within the 'Fly Bowl Experiments' excel 
% file, details found in the 'Temp Protocols' sheet
% 
% ES Dickinson, Yale University, Aug 2021

%get base folder pathway
switch getenv('COMPUTERNAME')
    case 'DENALI'
        baseFolder = 'E:\My Drive\Jeanne Lab\';
    case 'TOGIAK'
        baseFolder = 'G:\My Drive\Jeanne Lab\';
    case 'EVYNPC'
        baseFolder = 'G:\My Drive\Jeanne Lab\';
end

xlFile = [baseFolder 'Fly Bowl Experiments.xlsx'];

%load excel sheet data
[~,~,excelfile] = xlsread(xlFile, 'Temp Protocols');
protocolList = excelfile(:,1);

a = protocolList(listdlg('ListString', protocolList, 'PromptString',...
           'Select exp protocol', 'ListSize', [250, 400]));
protocol = a{:};


end