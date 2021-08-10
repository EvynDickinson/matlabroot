 
clear all; close all; clc 
origState = warning;
warning('off')
strtpoint = 470; %line in excel file to start on
behavior_list = {'stationary', 'walking', 'grooming', 'other', 'repeat video!', 'wrong phase last time!'};
phase_list = {'stance', 'swing', 'swing', 'swing', 'WHOOPS'};
% Find files from Excel:    
[excelfile, Excel] = load_flysummary;
 
% Pull out the structure names from the file:
structure_names.excelfile = excelfile(strtpoint:end, Excel.structure);
ind = 1;
for ii = 1:length(structure_names.excelfile)
    if ischar(structure_names.excelfile{ii})
        structure_names.text{ind} = (structure_names.excelfile{ii});
        ind = ind+1;
    end
end; clear ind

%find the unique structure names:
structure_names.unique = unique(structure_names.text);
indx = listdlg('ListString', structure_names.unique, 'SelectionMode', 'Multiple', 'ListSize', [300, 400]);
% fig = getfig;
%for each structure selected, load the data corresponding to that fly
for ii = 1:length(indx)
   %FIND STRUCTURE DETAILS
    structure_name = structure_names.unique{indx(ii)};
    location = find(strcmpi(structure_name, excelfile(strtpoint:end,Excel.new_struct_name)))+strtpoint-1;
    bb = 1;
    for tt = 1:length(location)
        a = cell2mat(excelfile(location(tt), Excel.structurenum));
        if ischar(a) %no value for the structure
        else %has a number in the structure
            flynum(bb) = a; %find total number of flies
            Location.rownum(bb) = location(tt); %get the rownumber for each fly
            bb = bb+1; %increase the index for the number of flies
        end
    end
    num.flies = max(flynum); clear flynum bb a tt
    % Error Check
    if sum(num.flies) == 0
        warndlg('Check name of flies in Excel File, no matching files found')
        return
    end
    
    %LOAD DATA INTO MATLAB
    % begin loop to load data from each fly in 'fly.num'
    
%     group(1:num.flies).behavior = struct([]);
%     group(1:num.flies).phase = struct([]);
    
    
    for kk = 1:num.flies
        tic
        % Find Specific Fly & Load Data:
        fprintf('\n Fly: \n')
        disp([excelfile(Location.rownum(kk), Excel.date) excelfile(Location.rownum(kk), Excel.flynum)])
        
        % Date and number of the fly that matches the desired one
        %get row number of the fly(kk) 1_0 fly:
        %gives a logical index of NaN = 1
        NaN_index = cellfun(@(excelfile) any(isnan(excelfile(:,Excel.date))), excelfile(:,Excel.date));
        date_rownum = Location.rownum(kk);
        fly_date = string(excelfile(date_rownum, Excel.date));
        fly_num = string(excelfile(Location.rownum(kk), Excel.flynum));
        genetic_cross = string(excelfile(Location.rownum(kk), Excel.genetic_cross));
        datestring = dateconverter(fly_date);
        fly_ID = [datestring '_fly' cell2mat(fly_num)];

        % Load a preview of the first video:
        for irep = 1:3
            fprintf(['\nStarting ' fly_ID ' rep ' num2str(irep) '\nFly num: ' num2str(kk) '\n'])
            for icond = 1:28
                vid_root = ['D:\Evyn Data Files\' char(fly_date) '\Fly ' char(fly_num) '\Raw Video\'];
                vid_name = [vid_root, fly_ID, ' R' num2str(irep), 'C', num2str(icond), ' Cam-C*'];
                FilePath = dir(vid_name);
                vid_name = FilePath.name;
                FrameRate = 50;
                fig = vidplayback([vid_root, vid_name], FrameRate);
                
                % Behavior classifier:
                idx = listdlg('ListString', behavior_list, 'SelectionMode', 'Single');
                group(kk).behavior{icond, irep} = behavior_list{idx};
                switch idx
                    case num2cell(1:4) %normal selection
                    case 5 %replay the video
                    case 6 %correct phase of previous video
                end
                
                
                
                if idx == 5 % replay the video
                    close(fig)
                    fig = vidplayback([vid_root, vid_name], round(FrameRate/2));
                    idx = listdlg('ListString', behavior_list(1:4), 'SelectionMode', 'Single');
                end
                if idx == 6 % correct the phase from the previous video
                    if irep == 1 && icond == 1 %the last video was the former video/fly
                       warndlg(['Manually adjust this fly ' num2str(kk) ' cond 3, rep 28'])
                    elseif icond == 1  % cond wrapped repetitions
                        rep = 28;
                        cond = icond - 1;
                        fprintf(['\n adjusted ' group(kk).phase{cond, rep} ' to ']) %...
                        if strcmpi(group(kk).phase{cond, rep}, 'stance')==1
                        group(kk).phase{cond, rep} = 'swing';
                        else 
                            group(kk).phase{cond, rep} = 'stance';
                        end
                        fprintf([group(kk).phase{cond, rep} '\n'])
                        group(kk).behavior{icond, irep} = ...
                            behavior_list{listdlg('ListString', behavior_list, 'SelectionMode', 'Single')};
                    else % normal situations
                        rep = irep;
                        cond = icond-1;
                        fprintf(['\n adjusted ' group(kk).phase{cond, rep} ' to ']) %...
                        if strcmpi(group(kk).phase{cond, rep}, 'stance')==1
                        group(kk).phase{cond, rep} = 'swing';
                        else 
                            group(kk).phase{cond, rep} = 'stance';
                        end
                        fprintf([group(kk).phase{cond, rep} '\n'])
                        group(kk).behavior{icond, irep} = ...
                            behavior_list{listdlg('ListString', behavior_list, 'SelectionMode', 'Single')};
                    end
                end
                
                
                % Swing or Stance?
                group(kk).phase(icond, irep) = ...
                        phase_list(listdlg('ListString', phase_list, 'SelectionMode', 'Single'));
%                 group(kk).phase{icond, irep} = questdlg('Swing or stance?', '', 'swing', 'stance', 'WHOOPS', 'stance');
                if strcmpi(group(kk).phase{icond, irep},'WHOOPS')
                    group(kk).behavior{icond, irep} = ...
                        behavior_list(listdlg('ListString', behavior_list(1:4), 'SelectionMode', 'Single'));
                    group(kk).phase{icond, irep} = questdlg('Swing or stance?', '', 'swing', 'stance', 'stance');
                end
                    
                close(fig)
                fprintf([' Finished R' num2str(irep) 'C' num2str(icond) '\n'])
            end
        end
    save([structure_name ' behavior class'], 'group')
    behavior_classification = group(kk);
    save([FilePath.folder '\Behavior Classification'], 'behavior_classification')
    toc
    end
  
    
end %loop to the next structure
beep

warning(origState)
       
num.flies = length(group);
num.conds = 28;
num.reps = 3;
for kk = 1:num.flies
    for cond = 1:num.conds
        for rep = 1:num.reps
            state = group(kk).behavior{cond, rep};
            phase = group(kk).phase{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            if ~ischar(phase)
                phase = cell2mat(phase);
            end
            switch (state)
                case 'stationary'
                    STATE = 1;
                case 'walking'
                    STATE = 2;
                case 'grooming'
                    STATE = 3;
                case 'other'
                    STATE = 4;
            end
            switch phase
                case 'stance'
                    PHASE = 1;
                case 'swing'
                    PHASE = 2;
            end 
            group(kk).STATE(cond, rep) = STATE;
            group(kk).PHASE(cond, rep) = PHASE;
        end
    end
    group(kk).walking = (group(kk).STATE==2);
    group(kk).stationary = (group(kk).STATE==1);
    group(kk).stance = (group(kk).PHASE==1);
    group(kk).swing = (group(kk).PHASE==2);
end

save([structure_name ' behavior class'], 'group')
 


%% Error check -- look for unresolved vidoes
num.flies = length(group);
num.conds = 28;
num.reps = 3;


for kk = 1:num.flies
    for cond = 1:num.conds
        for rep = 1:num.reps
            state = group(kk).behavior{cond, rep};
            phase = group(kk).phase{cond, rep};
            if sum(strcmpi(state, '-')) > 0
                
                NaN_index = cellfun(@(excelfile) any(isnan(excelfile(:,Excel.date))), excelfile(:,Excel.date));
                date_rownum = Location.rownum(kk);

                fly_date = string(excelfile(date_rownum, Excel.date));
                fly_num = string(excelfile(Location.rownum(kk), Excel.flynum));
                genetic_cross = string(excelfile(Location.rownum(kk), Excel.genetic_cross));
                datestring = dateconverter(fly_date);
                fly_ID = [datestring '_fly' cell2mat(fly_num)];

                vid_root = ['D:\Evyn Data Files\' char(fly_date) '\Fly ' char(fly_num) '\Raw Video\'];
                vid_name = [vid_root, fly_ID, ' R' num2str(rep), 'C', num2str(cond), ' Cam-C*'];
                FilePath = dir(vid_name);
                vid_name = FilePath.name;
                FrameRate = 50;
                fig = vidplayback([vid_root, vid_name], FrameRate);
                 
                % Behavior classifier:
                idx = listdlg('ListString', behavior_list, 'SelectionMode', 'Single');
                group(kk).behavior{cond, rep} = behavior_list{idx};
                close(fig)
            end
            
            if sum(strcmpi(phase, '-')) > 0
                
                NaN_index = cellfun(@(excelfile) any(isnan(excelfile(:,Excel.date))), excelfile(:,Excel.date));
                date_rownum = Location.rownum(kk);

                fly_date = string(excelfile(date_rownum, Excel.date));
                fly_num = string(excelfile(Location.rownum(kk), Excel.flynum));
                genetic_cross = string(excelfile(Location.rownum(kk), Excel.genetic_cross));
                datestring = dateconverter(fly_date);
                fly_ID = [datestring '_fly' cell2mat(fly_num)];

                vid_root = ['D:\Evyn Data Files\' char(fly_date) '\Fly ' char(fly_num) '\Raw Video\'];
                vid_name = [vid_root, fly_ID, ' R' num2str(rep), 'C', num2str(cond), ' Cam-C*'];
                FilePath = dir(vid_name);
                vid_name = FilePath.name;
                FrameRate = 50;
                fig = vidplayback([vid_root, vid_name], FrameRate);
                
                group(kk).phase(cond, rep) = ...
                        phase_list(listdlg('ListString', phase_list, 'SelectionMode', 'Single'));
                close(fig)
            end
        end
    end
end
disp('No more repeat videos')
        

%% CORRECTIONS:
err = str2double(inputdlg('Number of edits'));
if isempty(err); return; end
for ii = 1:err % for all errors to fix:
    %find fly
    flynum = str2double(inputdlg('Fly number?'));
    irep =   str2double(inputdlg('Rep number?'));    % missing rep number
    icond = str2double(inputdlg('Condition number?'));    % missing cond number    
    fprintf(['\n Behavior-- ' group(flynum).behavior{icond, irep} '\n'...
             ' Phase-- ' group(flynum).phase{icond, irep}  '\n'])
    switch questdlg('correction in?', '', 'behavior', 'phase', 'both', 'phase')
        case 'behavior'
            group(flynum).behavior{icond, irep} = ...
                        behavior_list(listdlg('ListString', behavior_list(1:4), 'SelectionMode', 'Single'));
        case 'phase'
             group(flynum).phase{icond, irep} = ...
                        phase_list(listdlg('ListString', phase_list(1:2), 'SelectionMode', 'Single'));
        case 'both'
             group(flynum).behavior{icond, irep} = ...
                        behavior_list(listdlg('ListString', behavior_list(1:4), 'SelectionMode', 'Single'));
             group(flynum).phase{icond, irep} = ...
                        phase_list(listdlg('ListString', phase_list(1:2), 'SelectionMode', 'Single'));
    end
end

for kk = 1:num.flies
    for cond = 1:num.conds
        for rep = 1:num.reps
            state = group(kk).behavior{cond, rep};
            phase = group(kk).phase{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            if ~ischar(phase)
                phase = cell2mat(phase);
            end
            switch (state)
                case 'stationary'
                    STATE = 1;
                case 'walking'
                    STATE = 2;
                case 'grooming'
                    STATE = 3;
                case 'other'
                    STATE = 4;
            end
            switch phase
                case 'stance'
                    PHASE = 1;
                case 'swing'
                    PHASE = 2;
            end 
            group(kk).STATE(cond, rep) = STATE;
            group(kk).PHASE(cond, rep) = PHASE;
        end
    end
    group(kk).walking = (group(kk).STATE==2);
    group(kk).stationary = (group(kk).STATE==1);
    group(kk).stance = (group(kk).PHASE==1);
    group(kk).swing = (group(kk).PHASE==2);
end

save([structure_name ' behavior class'], 'group')






%                 tic
%                 [~, ProcessedImages] = convert_vid_to_pixel_values([vid_root, vid_name], 0);
%                 toc

%                 % PREVIEW 200ms BEFORE STIM
%                 fig = getfig;
%                 for datapoint = 97:157
%                     set(fig, 'color','k');
%                     imshow(ProcessedImages(datapoint).data, 'border', 'tight')
%                     currAxes.Visible = 'off';
%                     clf('reset')
%                 end
%                 close all
%                 % Behavior classifier:
%                 idx = listdlg('ListString', behavior_list, 'SelectionMode', 'Single');
%                 if idx == 5
%                     fig = getfig;
%                     for datapoint = 97:157
%                         set(fig, 'color','k');
%                         imshow(ProcessedImages(datapoint).data, 'border', 'tight')
%                         currAxes.Visible = 'off';
%                         clf('reset')
%                     end
%                     
%                     close all
%                     idx = listdlg('ListString', behavior_list, 'SelectionMode', 'Single');
%                 end







        
%         % Load analyzed data *can take a sec or two:
%         folder_name = cell2mat(fly_date);
%         try
%             directory = ['D:\Evyn Data Files\' folder_name '\Fly ' char(fly_num) '\Analysis\'];
%             load([directory fly_ID '.mat']);
%         catch
%             try
%                 directory = ['E:\FicTrac Raw Data\' folder_name '\Analysis\'];
%                 load([directory fly_ID '.mat']);
%             catch
%                 directory = ['C:\matlabroot\FicTrac Raw Data\' folder_name '\Analysis\'];
%                 load([directory fly_ID '.mat']);
%             end
%         end
%         fprintf(['\n Loaded fly data: ' fly_ID '\n'])
        
        
%         % Add a video category to the data
%         Sstrt = 1;Sstp = 44;Cstrt = 46;Cstp = 61;
% %         [num, Fictrac, Matlab, labels] = NUM(fly.param);
%         for cond = 1:28
%             for rep = 1:3
%                 fly.Video.data(cond, rep).raw = ...
%                     [fly.Control.data(cond,rep).raw(Cstrt:Cstp,:); fly.Stim.data(cond,rep).raw(Sstrt:Sstp,:)];
%                   
%             end
%             fly.Video.speed(cond).data = ...
%                 [fly.Control.speed(cond).data(Cstrt:Cstp,:); fly.Stim.speed(cond).data(Sstrt:Sstp,:)];
%             fly.Video.rotvelocity(cond).data = ...
%                 [fly.Control.rotvelocity(cond).data(Cstrt:Cstp,:); fly.Stim.rotvelocity(cond).data(Sstrt:Sstp,:)];
%         end
%         
%         % Save the updated fly video:
%         newfile = [directory, fly_ID];
%         save(newfile, 'fly');
%         
%         % Copy the folder to google drive:
%         % Destination:
%         rootpath = 'G:\My Drive\Data\FicTrac Raw Data\';
%         dest = [rootpath, fly.param.folder_date, '\Fly ', fly_ID(end-2:end), '\Analysis'];
%         copyfile(directory, dest)

%         fprintf([' Number: ' num2str(kk) '/' num2str(num.flies) '\n '])
  
    

    
    

