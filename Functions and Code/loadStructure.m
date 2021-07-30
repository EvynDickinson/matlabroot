
function output = loadStructure
% Load Fly Structure (already built with "MIDSTIM_Step_4_groupflies.m")
% Select a structure to load
% variables to output
% output.structure_name = structure_name;
% output.fly = fly;
% output.group = group;
% output.num = num;
% output.traj = traj;



[filename, directory] = uigetfile('*.mat', 'Select your Fly Structure');
fullname = fullfile(directory, filename); 
load(fullname);
try fly = FLY; catch
end
if isfield(fly(1), 'parameters') 
    for kk = 1:length(fly)
        fly(kk).param = fly(kk).parameters;
    end
    fly = rmfield(fly,'parameters');
end
kk = 1;
parameters = fly(1).param;
fprintf(['\n Loaded ' filename '\n'])
setGlobalx(parameters)

% load labels & condition parameters
% [num, Fictrac, Matlab, labels] = NUM(fly(1).param);  
num = NUM(fly(1).param);  
num.fly = length(fly);
Type = {'Stim', 'Control'};

% Create a Figures Folder for this structure:
% figures_dir = [directory 'MN line figures\' filename(1:end-4) '\'];
figures_dir = [directory 'Interneuron Lines/' filename(1:end-4) '/'];
if ~isfolder(figures_dir)
    mkdir(figures_dir)
end 

cond_figures_dir = [figures_dir 'Conditions/'];
if ~isfolder(cond_figures_dir)
    mkdir(cond_figures_dir)
end 
clc
% ELIMINATE OUTLIERS
fly = eliminate_speed_outliers(fly, 5);

% variables to output
output.structure_name = structure_name;
output.fly = fly;
output.num = num;

% LOAD BEHAVIOR CLASSIFICATION DATA:
beep
switch questdlg('Load behavior classification data?')
    case 'Yes'
        load([directory, '/', structure_name,  ' behavior class'])
%         load([directory, '/behavior class/', structure_name,  ' behavior class'])
end

% TRAJECTORY PATHS AND HEATMAPS
ROI = 1:60;
beep
switch questdlg('Load former trajectory info?', 'Yes', 'No', 'Cancel', 'No')
    case 'Yes' 
        load([figures_dir, 'Trajectory Data'])
    case {'No', 'Cancel'}
        traj = [];
end
output.traj = traj;

% USE THE 'GROUP' STRUCTURE FROM 'Behavior_Categorization.m'

% Convert behvariors into numbers and then a filter:
for kk = 1:num.fly
    for cond = 1:num.conds
        for rep = 1:num.reps
            state = group(kk).behavior{cond, rep};
            phase = group(kk).phase{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            if ~ischar(phase)
                phase = cell2mat(phase);
            else
                PHASE = 0;
            end
            switch state
                case {'stationary', 'Stationary'}
                    STATE = 1;
                case {'walking', 'Walking'}
                    STATE = 2;
                case {'grooming', 'Grooming'}
                    STATE = 3;
                case {'other', 'Other'}
                    STATE = 4;
            end
            switch phase
                case {'stance', 'Stance'}
                    PHASE = 1;
                case {'swing', 'Swing'}
                    PHASE = 2;
                case {'-'}
                    PHASE = 0;
            end 
            group(kk).STATE(cond, rep) = STATE;
            group(kk).PHASE(cond, rep) = PHASE;
        end
    end
    group(kk).walking = (group(kk).STATE==2);
    group(kk).stationary = (group(kk).STATE==1);
    group(kk).grooming = (group(kk).STATE==3);
    group(kk).other = (group(kk).STATE==4);
    group(kk).stance = (group(kk).PHASE==1);
    group(kk).swing = (group(kk).PHASE==2);
end

output.group = group;


switch questdlg('Load joint angle info?', 'Yes', 'No', 'No') % JA_Data
    case 'Yes' 
        load([figures_dir, structure_name, ' joint angle data'])
    case {'No', 'Cancel'}
end


end
