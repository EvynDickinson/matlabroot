% 'DLC_Step1_loaddata.m'
% should be run after the data has been uploaded to Google drive and
% processed by Anipose
% Runs similarly to Step_3 for Fictrac data
% Step one analysis for DLC & anipose analyzed videos
% load and reorganize the angles & joint position data to fit the Fictrac data structure
% Do this for a single fly or a day of data

clearvars
clc
close all
warning('off')

% TODO section:
% Set base file paths:
fileroot = 'G:\My Drive\Evyn\Data\FicTrac Raw Data\'; %base folder for processed Anipose data
output_root = 'D:\Tuthill Lab Work\Anipose\';         %local PC save location
local_save = true;                                    %save data in the local PC folder?
online_save = false;                                  %save data to google drive folder?
% ----------------------
folder_date = '6.12.19';                              %folder to analyze
% ----------------------
% labels to match with anipose data 
Joints = {'CoFe', 'FeTi', 'TiTa'};
all_angles = {'L1_CF','L1_FTi','L1_TiTa',...
          'L2_CF','L2_FTi','L2_TiTa',...
          'L3_CF','L3_FTi','L3_TiTa',...
          'R1_CF','R1_FTi','R1_TiTa',...
          'R2_CF','R2_FTi','R2_TiTa',...
          'R3_CF','R3_FTi','R3_TiTa'};
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
leg_positions = {'body_wall','coxa', 'femur', 'tibia', 'tarsus'};

%% Select flies within the folder_date to load:
pathway = [fileroot folder_date]; % directory
list_flies = dir(pathway);
idx = 1;
for ii = 1:length(list_flies)
    a = list_flies(ii).name; %names of the flies for the day
    if strfind(a,'Fly')==1
        fly_num{idx} = a;
        idx = idx+1;
    else; continue;
    end
end; clear idx
%Select the flies to analyze
idx = listdlg('ListString', fly_num, 'SelectionMode', 'Multiple');
FliesToLoad = fly_num(idx);
nflies = length(FliesToLoad);
disp('Flies selected: ')
disp(FliesToLoad)
clear a ii idx list_flies


%% Load Anipose Data for each fly in 'FliesToLoad'

for ifly = 1:nflies
    anipose = [];
    % path to base dir for this fly
    flypath = [pathway '\' FliesToLoad{ifly} '\']; 
    flyID = generate_flyID(folder_date, FliesToLoad{ifly}(end-2:end));
    
    % load the param file for this fly to get num conds, reps etc.
    load([flypath, 'FicTrac Data/Parameters.mat']); 
    num = NUM(param);
    num.frames = param.basler_length*param.Basler_fps+1;

    % angle directory:
    angleDir = [flypath 'angles\'];
    % pose 3d directory:
    poseDir = [flypath 'pose-3d\'];
    % ---------------------------------------------------------------------
    % load the CSV files for each trial video
    for rep = 1:num.reps
        for cond = 1:num.conds
            prep = [ flyID ' R' num2str(rep) 'C' num2str(cond)...
                    '  ' param.conds_matrix(cond).label, '.csv'];
            % load angles files:
            metaData = [];    
            csvFile = [angleDir, prep];
            metaData = readtable(csvFile);
            angles.labels{cond,rep} = metaData.Properties.VariableNames;
            angles.raw{cond,rep} = table2array(metaData);
            % load pose-3d files:
            metaData = [];    
            csvFile = [poseDir, prep];
            metaData = readtable(csvFile);
            pose_3d.labels{cond,rep} = metaData.Properties.VariableNames;
            pose_3d.raw{cond,rep} = table2array(metaData);
            disp(prep)
        end
    end
    %TODO add a loading bar rather than a printed list of videos
    % DO NOT ERASE 'ANGLES' OR 'POSE_3D' or you have to load all the files again
    
    % ---------------------------------------------------------------------
    % Run basic data organization on the ANGLE csv data:
    
    % Joints and angles of interest:

    
    % extract and organize the angle information:
    % find the correct extraction data columns:
    num.joints = length(Joints);
    input(:,1) = angles.labels{1,1};
    %find the data columns for all angles
    for ii = 1:length(all_angles)
       angle_loc(ii) = find(strcmpi(all_angles(ii), input));
       label_check(ii,1) = input(angle_loc(ii));
    end
    angles.angle_loc = angle_loc; % from above
    angles.searchNames = all_angles; % save the joints search
    % extract the joints over time into the angles structure organized by
    % leg
    for cond = 1:num.conds
        for rep = 1:num.reps
           raw = angles.raw{cond,rep};
           % check that the number of frames is what was expected|replace with
           % nans
           if ~(size(raw,1)==num.frames) %frame numbers don't match=true
               raw = nan(num.frames,size(raw,2));
               fprintf(['\n Missing frames found in ' flyID ...
                        ' R' num2str(rep) 'C' num2str(cond) '\n'])
           end
           for leg = 1:6 %six legs on the fly
               ROI = leg*num.joints-(num.joints-1):leg*num.joints;
               loc = angle_loc(ROI);
               %extract the time series data:
               angles.Leg(leg).data(cond,rep).all = raw(:,loc); 
               a = angles.labels{cond,rep}(loc);
               angles.Leg(leg).labels = strrep(a,'_','-'); %swap underscores to dashes
               angles.Leg(leg).loc = loc;
           end
         end
    end
    % ---------------------------------------------------------------------
    % Run basic data organization on the POSE_3D csv data:
    pose_3d.leg_labels = leg_labels;
    for cond = 1:num.conds
     for rep = 1:num.reps
      for leg = 1:6 % for all legs
       for LP = 1:5 % leg position:
        pos = [leg_labels{leg} Alphabet(LP) '_x'];
        labels = pose_3d.labels{cond,rep};
        label_loc = find(strcmpi(pos, labels));
        raw = pose_3d.raw{cond,rep}(:,label_loc:label_loc+2);
        pose_3d.Leg(LP,leg).data{cond,rep}.all = raw;
       end
      end
     end
    end

    % ---------------------------------------------------------------------
    % Save the data - both Angles and Pose3d into the appropriate
    % folders/locations
    anipose.angles = angles;
    anipose.pose_3d = pose_3d;
    anipose.param = param;
    anipose.Joints = Joints;
    anipose.all_angles = all_angles;
    anipose.leg_labels = leg_labels;
    anipose.leg_positions = leg_positions;
    anipose.flyID = flyID;
    anipose.flypath = flypath;
    anipose.angleDir = angleDir;
    anipose.poseDir = poseDir;
    fly_file_name = [flyID '_AniposeData'];
    
    % Save 'anipose' structure into the respective locations:
    if local_save == true
        L_save = [output_root, folder_date,'\',FliesToLoad{ifly},'\Analysis'];
        if ~exist(L_save,'dir'); mkdir(L_save); end % make folder 
        save([L_save,'\', fly_file_name], 'anipose')
    end
    if online_save == true
        save([flypath, 'Analysis\' fly_file_name], 'anipose')
    end    
    
    fprintf(['\n Loaded: ' flyID '\n'])
end



warning('on')
disp('Done')









