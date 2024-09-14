

%% Load tracking points for courtship

baseFolder = 'S:\Evyn\Courthship Videos\';
dataFile = 'labels.v000.analysis.h5';
filePath = [baseFolder, dataFile];

occupancy_matrix = h5read(filePath,'/track_occupancy');
tracks = h5read(filePath,'/tracks');

%% Create data labels and structure organization 

