
%% Load Data from Matlab and FicTrac
close all; clear all; clc

origState = warning;
warning('off')

folder_date = '5.12.20'; % date of folder with videos/data acq.
Folder_date = folder_date; % date of folder with videos/data acq. 
pathway = ['F:\Evyn Data Files\' folder_date]; % directory
list_flies = dir(pathway);
for ii = 3:length(list_flies)
    fly_num{ii-2} = list_flies(ii).name; %names of the flies for the day
end
%Select the flies to analyze
indx = listdlg('ListString', fly_num, 'SelectionMode', 'Multiple');

for kk = 1:length(indx) % LOOP THROUGH THE ENTIRE PROGRAM FOR THESE FLIES
    ifly = indx(kk);    
    % Load the data:   (maybe make this a function...) 
    filepath = [pathway '\' fly_num{ifly} '\FicTrac Data\'];
    a = dir([filepath, '\*2020*.mat']);
    fly_ID = a(1).name(1:end-4); clear a 
    % 
    % folder_date(regexp(folder_date,'[.]'))=[];
    % fly_ID = [folder_date(1:end-2), '2019_fly' fly_num{ifly}(end-2:end)];

    %load the matlab data
    load([filepath fly_ID])
    %load the Fictrac data
    testdata = load([filepath fly_ID 'data_out.dat']);
    fprintf(['\n Loaded: ' fly_ID '\n'])

    % Update Variables on timing:  
    [num, Fictrac] = NUM(fly.param);
    % Add a change in heading column to fictrac data
    testdata = generate_change_in_heading(testdata);

    if ~isfield(param, 'cross')
        fly.param.cross = cell2mat(inputdlg('Cross name?'));
        param.cross = fly.param.cross;
    end
    
%     % %***ERROR CORRECTIONS FOR DATA BETWEEN: 11.16-11.28*** and fly1 from 2.6.19
%     %change name on conditions from start_loc to start_num for fictrac data 
%     conditions.fictracdata.Stim.start_num = conditions.fictracdata.Stim.start_loc; 
%     conditions.fictracdata.Stim.end_num = conditions.fictracdata.Stim.end_loc; 
%     conditions.fictracdata.Control.start_num = conditions.fictracdata.Control.start_loc;  
%     conditions.fictracdata.Control.end_num = conditions.fictracdata.Control.end_loc;
%     % adjust the control to be one more frame
    conditions.fictracdata.Control.start_num = conditions.fictracdata.Control.start_num-1;
%     % disp(conditions.fictracdata.Control.end_num - conditions.fictracdata.Control.start_num)

    % Find Fictrac timing (start, stop, stim, control) locations 
    conditions = assign_frames_fictrac(conditions, testdata, num);

    % Visual inspection for missing data
    fprintf('\n Control data points: \n');
    disp(conditions.fictracdata.Control.end_loc - conditions.fictracdata.Control.start_loc)
    fprintf('\n Stim data points: \n');
    disp(conditions.fictracdata.Stim.end_loc - conditions.fictracdata.Stim.start_loc)

    %-------------------------------------------------------------------------%
    %                           Error Checking Below                          %
    %-------------------------------------------------------------------------%
    %print missing frame information
    missing_values = alert_missing_frames(conditions, num);
    if missing_values(1).data(1) == 27272727
        %no missing frames
    else %missing data points
        num.missing = length(missing_values);
        for qq = 1:num.missing
            [testdata, conditions] = insert_averages_in_holes(testdata, conditions, param, missing_values(qq).data);
            conditions = assign_frames_fictrac(conditions, testdata, num);
            fprintf(['\n Round ' num2str(qq) ' of ' num2str(num.missing) ' missing numbers resolved \n'])    
        end
    end
    clear qq
    fly.testdata = testdata;
    fly.data = data;
    fly.conditions = conditions;
    %-------------------------------------------------------------------------%

    % create a combined data matrix with Fictrac and Matlab data
    %(the start/stop information for each cond is just the Fictrac info) 
    [fly.all_data, labels] = combine_fictrac_n_analogue(fly.testdata, fly.data);

    % pull out the data for each conditions from the 'all_data' matrix
    % & add fields for speed, velocity, & activity level
    fly = separate_conditions(fly);
    fly.param.fly_name = [fly.param.folder_date ' Fly ' param.fly_num(1) '-' param.fly_num(3)];

    % check/create folder for data analysis within raw data
    analysis_directory = [pathway '\Fly ' param.fly_num '\Analysis\'];
    if ~isfolder(analysis_directory)
        mkdir(analysis_directory)
        fprintf('\n Created ''Analysis'' folder \n')
    end

    % Flat path with condition subunits:
    [fig, fly.param.total_distance_traveled] = flat_path_fig(fly);
    export_fig(fig,[analysis_directory, param.matlab_data_file '-Flat Path and subunits.pdf'],...
                    '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
    fprintf('\n Figure saved: flat path trajectory \n')
    close all
    
    % Add a video category to the data for Pierre and Katie tracking
     Sstrt = 1;Sstp = 44;Cstrt = 46;Cstp = 61;
    for cond = 1:num.conds
        for rep = 1:num.reps
            fly.Video.data(cond, rep).raw = ...
                [fly.Control.data(cond,rep).raw(Cstrt:Cstp,:); fly.Stim.data(cond,rep).raw(Sstrt:Sstp,:)];

        end
            fly.Video.speed(cond).data = ...
                [fly.Control.speed(cond).data(Cstrt:Cstp,:); fly.Stim.speed(cond).data(Sstrt:Sstp,:)];
            fly.Video.rotvelocity(cond).data = ...
                [fly.Control.rotvelocity(cond).data(Cstrt:Cstp,:); fly.Stim.rotvelocity(cond).data(Sstrt:Sstp,:)];
    end
        

%     switch questdlg('Save file?')
%         case 'Yes'
            % Save the fly structure into analysis folder:
            save([analysis_directory fly_ID], 'fly')
            fprintf('\n Saved ''fly'' structure \n')
            % Write information to Excel
            excel_data_write(fly_ID, Folder_date, fly.param.total_distance_traveled)
            fprintf('\n Information written to Excel file. \n')
            fprintf(['\n Finished fly: ' fly_ID '\n'])
%         case 'No'
%             fprintf('\n No information saved \n')
%     end

    %-------------------------------------------------------------------------%
    %                               FIGURES
    %-------------------------------------------------------------------------%
    % switch questdlg('Make figures for the fly?')
    %     case 'No'
    %         fprintf('\n Ending program.\n')
    %         return
    %     case 'Yes'
    %         fprintf('\n Creating figures...\n')          
    % end

%     % Summary of one condition/stimulus 
%     for rep = 1:num.reps
%         for cond = 1:num.conds
%             fig = condition_summary_fig(fly, cond, rep);
%             fig_name{cond} = ['Condition Summary R' num2str(rep) 'C' num2str(cond)];
%             export_fig(fig,[analysis_directory fig_name{cond} '.pdf'],...
%                        '-pdf','-nocrop', '-r60' , '-painters', '-rgb');
% 
%         end
%         % append the pdfs into one big file
%     end
    
    
    % fprintf('\n Figure saved: flat path trajectory \n')
    close all
    % 
    % % Video of the plots coming in one frame at a time
    % render_condition_vid(fly, cond, rep)
    % 
    % % To create a video of the fly running while the data is plotting, run:
    % % 'Combined_vid_n_plots.m'

end
fprintf('\n Done!\n')

warning(origState)


