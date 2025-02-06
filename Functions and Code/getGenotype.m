
function genotype = getGenotype
% genotype = getGenotype;

folder = getCloudPath;
folder = folder(1:end-5);

file_id = fopen([folder 'cross.txt']);
A = fscanf(file_id, '%s');
genotypeList = strsplit(A,','); 
defaultIdx = find(strcmpi(genotypeList,'Berlin-WT'));

idx = listdlg('PromptString', 'Select fly genotype:','ListString',genotypeList,'SelectionMode','single','InitialValue',defaultIdx,'ListSize',[225,450]);
genotype = genotypeList{idx};
