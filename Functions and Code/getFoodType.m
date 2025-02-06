

function [food, wellLoc] = getFoodType(select_well)
% genotype = getGenotype(select_well);

folder = getCloudPath;
folder = folder(1:end-5);

file_id = fopen([folder 'FoodWellOptions.txt']);
A = fscanf(file_id, '%s');
foodList = strsplit(A,','); 
defaultIdx = find(strcmpi(foodList,'Caviar'));

idx = listdlg('PromptString', 'Select fly genotype:','ListString',foodList,'SelectionMode','single','InitialValue',defaultIdx,'ListSize',[225,450]);
food = foodList{idx};

wellOpt = {'1','2','3','4'};
if nargin>0 && select_well
    num = listdlg('PromptString', 'Food well?','ListString',wellOpt,'ListSize',[100,100]);
end
wellLoc = (wellOpt{num});