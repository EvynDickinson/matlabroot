
% https://sleap.ai/notebooks/Analysis_examples.html <-- python code for data example

filepath = 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\labels.v000.analysis.h5';
occupancy_matrix = h5read(filepath,'/track_occupancy');
tracks_matrix = h5read(filepath,'/tracks');


h5disp(filepath)


