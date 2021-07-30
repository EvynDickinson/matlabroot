
function conditions = sort_stims(laser_trig, param, data)
% Sort laser trigger start points from order randomization
% Create structures for frame start and end locations and timing
% 
% Input: 
% 'laser_trig' [chronological start points for stimulus]
% 'param' [contains the randomization order]
% Output:
% 'conditions' [structure that contains the frame numbers and locations
%              of the start points for each conditon, sorted by rep]
% saves figure of the arena signals for each condition&rep
% 
% ES Dickinson, University of Washington, Dec 2018

%Laser_trigs.start = ; %% start points for the stimulus, chronologically
for jj = 1:param.num_reps
    randomization_order(:,jj) = param.conditions_rand(:,jj);
    cond_starts(randomization_order(:,jj),jj) = laser_trig.start(:,jj);
end

% create a structure with the condition start and stop location info
conditions.matlabdata.Stim.start_loc = cond_starts;
stim_length = (param.acq_in_fs*param.OL_time);
one_unit = round(param.acq_in_fs*(1/param.fps));
for jj = 1:param.num_reps
    conditions.matlabdata.Stim.end_loc(:,jj) = conditions.matlabdata.Stim.start_loc(:,jj)+stim_length;
    conditions.matlabdata.Control.start_loc(:,jj) = conditions.matlabdata.Stim.start_loc(:,jj)-stim_length;
    conditions.matlabdata.Control.end_loc(:,jj) = conditions.matlabdata.Stim.start_loc(:,jj)-one_unit;
end

% -------------Figure with analogue signals for the conditions------------%
fprintf('\n Saving figure...\n')
fig = figure; set(fig, 'pos',[10 10 1600 980], 'color','w'); hold all
figName = [param.matlab_data_file ' condition alignment'];
set(fig, 'Name', figName)
for ii = 1:param.num_conds
    subplot(4,7,ii)
    x = (-stim_length-2000:stim_length+4000); %% time points
    for jj = 1:param.num_reps
        strt = conditions.matlabdata.Control.start_loc(ii,jj)-2000;
        stp = conditions.matlabdata.Stim.end_loc(ii,jj)+4000;
        hold all;
        plot(x, data(8,strt:stp), 'b'); %basler sig
        plot(x, data(3,strt:stp), 'r'); %Arena x pos
        plot(x, data(7,strt:stp), 'g'); %Laser trigger
        plot(x, data(5,strt:stp), 'm'); %AO sign trigger
        stim_start = data(1,conditions.matlabdata.Stim.start_loc(ii,jj)); %for a vline with this timing info
        vline(0, 'k--')
    end
    axis tight
    title(['C' num2str(ii) ': ' num2str(param.conds_matrix(ii).opto) ' sec'])
    switch ii
        case 1
            ylabel('Str CW')
        case 8
            ylabel('Str CCW')
        case 15
            ylabel('Rot CW')
        case 22
            ylabel('Rot CCW')
    end
end
save_dir = [param.fly_dir 'FicTrac Data\'];
if ~isfolder(save_dir)
    mkdir(save_dir)
end
result = save_figure(fig, [save_dir figName], '-pdf');    

% ------------------------------------------------------------------------%

end