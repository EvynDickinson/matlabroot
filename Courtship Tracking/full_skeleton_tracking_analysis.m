

%% Load tracking points for courtship
clear; clc;
baseFolder = '/Users/evyndickinson/Documents/Courtship Videos/';
dataFile = 'labels.v001.000_Courtship 0001.analysis.h5';
filePath = [baseFolder, dataFile];

occupancy_matrix = h5read(filePath,'/track_occupancy');
tracks = h5read(filePath,'/tracks');
node_names = h5read(filePath, '/node_names');
edge_names = h5read(filePath, '/edge_names');
track_names = h5read(filePath, '/track_names');


%% Create data labels and structure organization 

frame_loc = find(occupancy_matrix==1);
nframes = length(frame_loc);

% tracks organization = frames, body position, fly instance
fig = getfig('', false);

scatter()
frame_loc(1)


frame = frame_loc(1);

figure;
scatter(tracks(frame,:,1), tracks(frame,:,1))

