
%get base folder pathway
switch getenv('COMPUTERNAME')
    case 'DENALI'
        baseFolder = 'E:\My Drive\Jeanne Lab\DATA\';
    case 'TOGIAK'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'EVYNPC'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
end

params.well_1 = 'Plant';
params.well_2 = 'Empty';
params.well_3 = 'Empty';
params.well_4 = 'Empty';



params.genotype = select_cross; % choose the fly genotypes
params.protocol = select_protocol;  % choose experiment protocol

% move params into structure:
params.date = dirName;
params.expID = videoNames;
params.num = num;


writeExptoExcel(params) % save the data to the experiement list

% resave with same name: 
save([baseFolder dirName '\' videoNames 'dataMat'])


%% 
switch getenv('COMPUTERNAME')
    case 'DENALI'
        baseFolder = 'E:\My Drive\Jeanne Lab\DATA\';
    case 'TOGIAK'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'EVYNPC'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
end

for i = 1:num.vids
    filename = [baseFolder dirName '\analysis\' videoNames '_' num2str(i) ' data'];
    load(filename);
    
    params.well_1 = 'Plant';
    
    save(filename, 'videoData', 'params')
end


%% 

diff([tempLogStart(:,3), tempLogEnd(:,3)])


























