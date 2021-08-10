
%% Load Fly Structure (already built with "MIDSTIM_Step_4_groupflies.m")
clc; close all; clear all

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
fprintf(['\n Loaded ' filename '\n'])
parameters = fly(1).param;


% for kk = 1:7
%     fly(kk).param.cross = '81A06-gal4xgtACR1';
% end
% save(fullname)

% load labels & condition parameters
[num, Fictrac, Matlab, labels] = NUM(fly(1).param);  
num.fly = length(fly);
Type = {'Stim', 'Control'};

% Create a Figures Folder for this structure:
figures_dir = [directory 'Figures\' filename(1:end-4) '\'];
if ~isfolder(figures_dir)
    mkdir(figures_dir)
end 

cond_figures_dir = [directory 'Figures\' filename(1:end-4) '\Conditions\'];
if ~isfolder(cond_figures_dir)
    mkdir(cond_figures_dir)
end 

% ELIMINATE OUTLIERS
fly = eliminate_speed_outliers(fly, 5);

%% Velocity Average Combined across flies

for kk = 1:num.fly
    for tt = 1:2
        for cond = 1:num.conds
            raw = nanmean(fly(kk).(Type{tt}).speed(cond).data,2);
            speed.(Type{tt}).avg(cond).data(:,kk) = raw;
%             speed.(Type{tt}).avg(cond).err(:,kk) = raw;
        end
    end
end
   
for tt = 1:2
    clear raw
    idx = 0;
    for cond = 1:14
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).str(:,idx) = nanmean(raw,2);
    end
    speed.(Type{tt}).strerr = nanstd(speed.(Type{tt}).str,0,2)/sqrt(12);
    speed.(Type{tt}).stravg = nanmean(speed.(Type{tt}).str,2);
    
    clear raw
    idx = 0;
    for cond = 15:28
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).rot(:,idx) = nanmean(raw,2);
    end
    speed.(Type{tt}).roterr = nanstd(speed.(Type{tt}).rot,0,2)/sqrt(12);
    speed.(Type{tt}).rotavg = nanmean(speed.(Type{tt}).rot,2);
end
    
% Rotational Velocity

for kk = 1:num.fly
    for tt = 1:2
        for cond = 1:num.conds
            raw = nanmean(fly(kk).(Type{tt}).rotvelocity(cond).data,2);
            rotvelocity.(Type{tt}).avg(cond).data(:,kk) = raw;
        end
    end
end

for tt = 1:2
    clear raw
    idx = 0;
    for cond = 1:14
        idx = idx+1;
        raw = rotvelocity.(Type{tt}).avg(cond).data;
        
        if cond >=5 && cond <=14
            rotvelocity.(Type{tt}).str(:,idx) = nanmean(raw,2);
        else
            rotvelocity.(Type{tt}).str(:,idx) = nanmean(raw,2)*-1;
        end
       
    end
    rotvelocity.(Type{tt}).strerr = nanstd(rotvelocity.(Type{tt}).str,0,2)/sqrt(12);
    rotvelocity.(Type{tt}).stravg = nanmean(rotvelocity.(Type{tt}).str,2);
    
    clear raw
    idx = 0;
    for cond = 15:28
        idx = idx+1;
        raw = rotvelocity.(Type{tt}).avg(cond).data;
        
        if cond >=22 && cond <=28
            rotvelocity.(Type{tt}).rot(:,idx) = nanmean(raw,2);
        else
            rotvelocity.(Type{tt}).rot(:,idx) = nanmean(raw,2)*-1;
        end
    end
    rotvelocity.(Type{tt}).roterr = nanstd(rotvelocity.(Type{tt}).rot,0,2)/sqrt(12);
    rotvelocity.(Type{tt}).rotavg = nanmean(rotvelocity.(Type{tt}).rot,2);
end
       

fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
subplot(2,2,1) %str stim
hold all
% stim with light
    temp.x = (-60:60);
    temp.Y1 = speed.Control.stravg;
    temp.Y2 = speed.Stim.stravg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.strerr(1:end-1); speed.Stim.strerr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)
  
vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Speed during Walking')

subplot(2,2,2) %str stim
hold all
% stim with light
    temp.Y1 = speed.Control.rotavg;
    temp.Y2 = speed.Stim.rotavg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.roterr(1:end-1); speed.Stim.roterr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)

vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Speed during Turning')

% fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
subplot(2,2,3) %str stim
hold all
% stim with light
    temp.x = (-60:60);
    temp.Y1 = rotvelocity.Control.stravg;
    temp.Y2 = rotvelocity.Stim.stravg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [rotvelocity.Control.strerr(1:end-1); rotvelocity.Stim.strerr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)
  
vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Rot Vel (deg/s)')
title('Rotational Velocity during Walking')

subplot(2,2,4) %str stim
hold all
% stim with light
    temp.Y1 = rotvelocity.Control.rotavg;
    temp.Y2 = rotvelocity.Stim.rotavg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [rotvelocity.Control.roterr(1:end-1); rotvelocity.Stim.roterr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)

vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Rot velocity (deg/s)')
title({'Rotational Velocity during Turning'; [fly(1).param.cross ', N = ' num2str(num.fly)]})

% save_figure(fig, [figures_dir, filename(1:end-4), ' Walking and Turning Nof3'])

save_figure(fig, [figures_dir, filename(1:end-4), ' Walking and Turning'])


%% Figures of indivudual CW & CWW velocity and rotvelocity
%create label for combining CW and CCW data:
combined_label = create_combined_label(parameters);

% ROTATIONAL VELOCITY:
%cond placed as 2nd arg of MIDSTIM_(rot)velocity_fig returns cond & cond+7(CW&CCW)
%cond must be between 0-7;
[~, rotvelocity, rotvelocity_combined] = MIDSTIM_rotvelocity_fig(fly);

% SPEED | Velocity:
[~, velocity, velocity_combined] = MIDSTIM_velocity_fig(fly);

% Save figures for all conditions with CW % CCW data visualized together
fig_action = questdlg('Desired action for the rotational velocity figures? (''Esc'' for nothing)',...
        'Rotation Velocity', 'View Only', 'Save Individual', 'Save and Append', 'View Only');
if length(fig_action)>2 %anything but 'esc'
    %generate the figures  
    for cond = 1:7  
        figure_name(cond) = {[cond_figures_dir filename(1:end-4) ' rotational velocity combined cond ' num2str(cond) '.pdf']};
        fig(cond) = MIDSTIM_rotvel_and_speed_fig(fly, velocity, rotvelocity, cond);
    end
    %determine saving procedure
  switch fig_action
    case 'View Only'  
       uiwait(fig(1)) %wait to proceed until the last figure is closed
    case 'Save Individual'
      for cond = 1:7
        export_fig(fig(cond), figure_name{cond},'-pdf','-nocrop', '-r300' , '-painters', '-rgb');
        close all
      end
    case 'Save and Append'
      for cond = 7:-1:1
        export_fig(fig(cond), figure_name{cond},'-pdf','-nocrop', '-r300' , '-painters', '-rgb');
        close(fig(cond))
      end
      append_pdfs([figures_dir filename(1:end-4) ' rotational velocity combined.pdf'], figure_name{:})
      close all
      fprintf('\n Rotational velocity and speed figures appended!\n')
  end
end; clear figure_name fig_action

%% Visualize speed|rotational velocity averages without slow flies
% Remove Slow Flies & show filtered vs unfiltered
% Fig   
[activity_level, fig] = Activity_Level_Index(fly, 0.3, 0.5); %Min activity percent, min speed
save_figure(fig, [figures_dir filename(1:end-4) ' filtered vs unfiltered  speed'])

% Fig: Speed of all conditions combined
fig = MIDSTIM_filtered_fig(fly, activity_level, velocity_combined);
save_figure(fig, [figures_dir filename(1:end-4) ' filtered speed all traces'])
 
% Filtered by Speed Struct with avg and err:
velocity_filtered = MIDSTIM_filter_combined_data(velocity_combined, activity_level);
rotvelocity_filtered = MIDSTIM_filter_combined_data(rotvelocity_combined, activity_level);

%% Fourier Transform of Conditions Roll data: [mostly ignore this]
% Fig
cond = 5;
fig = MIDSTIM_make_fourier_transform(fly, cond);
save_figure(fig, [figures_dir filename(1:end-4) ' Roll Data Fourier Transform Cond ' num2str(cond)])

%% Analyze slope of effect and recovery after stimulus 
num.peak_select_buffer = 9;
% ***unfiltered data***
% rotvelocity -- all
[rotvelocity_analysis_unfiltered, fig] = MIDSTIM_analysis_effect_n_recovery(rotvelocity_combined, fly, num.peak_select_buffer);
save_figure(fig, [figures_dir filename(1:end-4) ' unfiltered rot velocity peak effect selection'])
% speed -- all
[velocity_analysis_unfiltered, fig] = MIDSTIM_analysis_effect_n_recovery(velocity_combined, fly, num.peak_select_buffer);
save_figure(fig, [figures_dir filename(1:end-4) ' unfiltered velocity peak effect selection'])

% ***filtered data***
%rotvelocity -- filtered
[rotvelocity_analysis, fig] = MIDSTIM_analysis_effect_n_recovery(rotvelocity_filtered, fly, num.peak_select_buffer);
save_figure(fig, [figures_dir filename(1:end-4) ' rot velocity peak effect selection'])
%speed -- filtered
[velocity_analysis, fig] = MIDSTIM_analysis_effect_n_recovery(velocity_filtered, fly, num.peak_select_buffer);
save_figure(fig, [figures_dir filename(1:end-4) ' velocity peak effect selection'])

% Basic Trial Tuning Curves
fig = MIDSTIM_basic_tuning_curves(rotvelocity_analysis, velocity_analysis, parameters);
save_figure(fig, [figures_dir filename(1:end-4) ' basic tuning curves'])

%% Single Trial Analysis - velocity
% Generate individual response (speed|velocity) analyses for each trial within each condition
for cond = 1:num.conds/2
   AD_velocity_analysis(cond) = MIDSTIM_trialsbytrial_analysis_effect_n_recovery...
                                (velocity_filtered, parameters, num.peak_select_buffer, cond);
end

% % Group the data into N bins based on their control steady state average
% num.bins = 6;
% [AD_velocity_single_trial_bins] = MIDSTIM_parameter_filter(AD_velocity_analysis, parameters, num.bins);
% 
% % Plot the single trial analysis tuning curves based on steady state
% % condition (3 figs: effect change, rate of response, rate of recovery)
% % replace the error bars with shaded regions
% fig = MIDSTIM_speed_based_tuning_curves_FIG(AD_velocity_analysis, AD_velocity_single_trial_bins, parameters);
% 
% Look at the area under the curves as a response analysis -- chunk out for
% various periods of time

% Overlay the average fly response to the various light lengths -- is the
% start the same? Do the reponses have similar trends? What about for speed
% based groups?

  
%% Single Trial Analysis - rotational velocity
% Generate individual response (speed|velocity) analyses for each trial within each condition
for cond = 1:num.conds/2
   AD_rotvel_analysis(cond) = MIDSTIM_trialsbytrial_analysis_effect_n_recovery...
                             (rotvelocity_filtered, parameters, num.peak_select_buffer, cond);
end

% % Group the data into N bins based on their control steady state average
% num.bins = 6;
% [AD_rotvel_single_trial_bins] = MIDSTIM_parameter_filter(AD_rotvel_analysis, parameters, num.bins);
% 
% % Plot the single trial analysis tuning curves based on steady state
% % condition (3 figs: effect change, rate of response, rate of recovery)
% % replace the error bars with shaded regions
% 
% %use speed filter rather than rotvelcity for splitting up the data
% fig = MIDSTIM_speed_based_tuning_curves_FIG(AD_rotvel_analysis, AD_velocity_single_trial_bins, parameters);

% Overlay the average fly response to the various light lengths -- is the
% start the same? Do the reponses have similar trends? What about for speed
% based groups? 

% Rotational velocity aligned
fig = MIDSTIM_convergent_parameters_graphs(AD_rotvel_analysis, parameters);
save_figure(fig(1), [figures_dir filename(1:end-4) ' ' fig(1).Name ' light activation aligned'])
save_figure(fig(2), [figures_dir filename(1:end-4) ' ' fig(2).Name ' light activation aligned'])

% Speed aligned
fig = MIDSTIM_convergent_parameters_graphs(AD_velocity_analysis, parameters);
save_figure(fig(1), [figures_dir filename(1:end-4) ' ' fig(1).Name ' light activation aligned'])
save_figure(fig(2), [figures_dir filename(1:end-4) ' ' fig(2).Name ' light activation aligned'])

%% Basic Tuning Curve with error

fig = MIDSTIM_tuning_curves_with_error(AD_velocity_analysis, AD_rotvel_analysis, parameters);
save_figure(fig, [figures_dir filename(1:end-4) ' rotvel time to peak response histogram']);
  

%%
%   save('Testing File 1072019')
%   load('Testing File 1072019')

%% Manual Selection of Effects:
% ** NEED TO COMBINE THE TWO DIRECTIONS HERE ***
% Find the fly average and then the grouped average
for tt = 1:2
  for cond = 1:num.conds
    % Add average of each fly into the 'Grouped Data' structures
    for kk = 1:num.fly
       temp.avg_speed = mean(fly(kk).(Type{tt}).speed(cond).data,2);
       GD_velocity(cond).(Type{tt})(:,kk) = temp.avg_speed; 
       temp.avg_rot = mean(fly(kk).(Type{tt}).rotvelocity(cond).data,2);
       GD_rotvelocity(cond).(Type{tt})(:,kk) = temp.avg_rot; 
    end
    % Average across all flies
    temp.velmean = mean(GD_velocity(cond).(Type{tt}),2);
    GD_velocity(cond).([Type{tt} '_avg']) = temp.velmean;
    temp.rotmean = mean(GD_rotvelocity(cond).(Type{tt}),2);
    GD_rotvelocity(cond).([Type{tt} '_avg']) = temp.rotmean;
    % Error across all flies
    temp.velerr = std(GD_velocity(cond).(Type{tt}),0,2)/sqrt(num.fly);
    GD_velocity(cond).([Type{tt} '_err']) = temp.velerr;
    temp.roterr = std(GD_rotvelocity(cond).(Type{tt}),0,2)/sqrt(num.fly);
    GD_rotvelocity(cond).([Type{tt} '_err']) = temp.roterr;
  end
end

%% Make a figure! Overlay Control and Laser condition

fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
    for cond = 1:num.conds
        subplot(4,7,cond)
        switch cond
            case num2cell(1:7)
                COND = 1;
            case num2cell(8:14)
                COND = 8;
            case num2cell(15:21)
                COND = 15;
            case num2cell(22:28)
                COND = 22;
        end
            temp.x = -num.stim_length:num.stim_length;
            temp.y = [GD_velocity(COND).Control_avg(1:end-1); GD_velocity(COND).Stim_avg];
            temp.err = [GD_velocity(COND).Control_err(1:end-1); GD_velocity(COND).Stim_err];
            hold all
            fill_data = error_fill(temp.x, temp.y, temp.err);
            h = fill(fill_data.X, fill_data.Y, get_color('grey'), 'EdgeColor','none');
                set(h, 'facealpha', 0.2)
            hh = plot(temp.x, temp.y, 'k', 'LineWidth', 2);
%          cond = 3;
            temp.x = -num.stim_length:num.stim_length;
            temp.y = [GD_velocity(cond).Control_avg(1:end-1); GD_velocity(cond).Stim_avg];
            temp.err = [GD_velocity(cond).Control_err(1:end-1); GD_velocity(cond).Stim_err];
            fill_data = error_fill(temp.x, temp.y, temp.err);
            h = fill(fill_data.X, fill_data.Y, get_color('lightgreen'), 'EdgeColor','none');
                set(h, 'facealpha', 0.2)
            hh = plot(temp.x, temp.y, 'g', 'LineWidth', 2);

        vline([0, parameters.conds_matrix(cond).opto*num.fps], 'k-')
    end

    
%% Average between all the conditions:

for tt = 1:2
  for cond = 1:num.conds
    % Add average of each fly into the 'Grouped Data' structures
    for kk = 1:num.fly
       temp.avg_speed = mean(fly(kk).(Type{tt}).speed(cond).data,2);
       GD_velocity(cond).(Type{tt})(:,kk) = temp.avg_speed; 
       temp.avg_rot = mean(fly(kk).(Type{tt}).rotvelocity(cond).data,2);
       GD_rotvelocity(cond).(Type{tt})(:,kk) = temp.avg_rot; 
    end
    % Average across all flies
    temp.velmean = mean(GD_velocity(cond).(Type{tt}),2);
    GD_velocity(cond).([Type{tt} '_avg']) = temp.velmean;
    temp.rotmean = mean(GD_rotvelocity(cond).(Type{tt}),2);
    GD_rotvelocity(cond).([Type{tt} '_avg']) = temp.rotmean;
    % Error across all flies
    temp.velerr = std(GD_velocity(cond).(Type{tt}),0,2)/sqrt(num.fly);
    GD_velocity(cond).([Type{tt} '_err']) = temp.velerr;
    temp.roterr = std(GD_rotvelocity(cond).(Type{tt}),0,2)/sqrt(num.fly);
    GD_rotvelocity(cond).([Type{tt} '_err']) = temp.roterr;
  end
end 

%% Overlay of all the tracks
% offset to zero
figure;
for cond = 1:14
    subplot(2, 7, cond)
    hold all
    for kk = 1:num.fly
        temp.offset1 = GD_velocity(cond+7).Stim(1,kk);
        temp.y1 = [GD_velocity(cond+7).Control(1:end-1,kk)-temp.offset1; GD_velocity(cond+7).Stim(:,kk)-temp.offset1];
        temp.offset = GD_velocity(cond).Stim(1,kk);
        temp.y = [GD_velocity(cond).Control(1:end-1,kk)-temp.offset; GD_velocity(cond).Stim(:,kk)-temp.offset];
        plot(temp.y, 'color', 'b', 'linewidth', 1)
        plot(temp.y1, 'color', 'm', 'linewidth', 1)
    end
    vline([60, 60+parameters.conds_matrix(cond).opto*num.fps], 'k-')
    title(combined_label{cond})
    hline(0, 'k-')
end





%% Figure Combining all the params:

for kk = 1:num.fly
    for tt = 1:2
        for cond = 1:num.conds
            raw = mean(fly(kk).(Type{tt}).speed(cond).data,2);
            speed.(Type{tt}).avg(cond).data(:,kk) = raw;
%             speed.(Type{tt}).avg(cond).err(:,kk) = raw;
        end
    end
end
for tt = 1:2
    clear raw
     % str stim control
    idx = 0;
    for cond = [1,8]
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).strnolight(:,idx) = mean(raw,2);
    end
    % rot stim control
    idx = 0;
    for cond = [15,22] 
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).rotnolight(:,idx) = mean(raw,2);
    end
    % str stim
    idx = 0;
    for cond = [2:7,9:14] 
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).strlight(:,idx) = mean(raw,2);
    end
    speed.(Type{tt}).errstrlight = std(speed.(Type{tt}).strlight,0,2)/sqrt(12);
    % rot stim
    idx = 0;
    for cond = [16:21,23:28] 
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).rotlight(:,idx) = mean(raw,2);
    end
    speed.(Type{tt}).errrotlight = std(speed.(Type{tt}).strlight,0,2)/sqrt(12);
end

%% plot the graphs:

fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
subplot(1,2,1) %str stim
hold all
% stim with light
temp.x = (-60:60);
    temp.Y1 = mean(speed.Control.strlight,2);
    temp.Y2 = mean(speed.Stim.strlight,2);
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.errstrlight(1:end-1); speed.Stim.errstrlight];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('red'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'r', 'linewidth', 1)
% control no light
    temp.y1 = mean(speed.Control.strnolight,2);
    temp.y2 = mean(speed.Stim.strnolight,2);
    temp.y = [temp.y1(1:end-1); temp.y2];
    temp.x = (-60:60);
    plot(temp.x, temp.y, 'color', 'k', 'linewidth', 1)    
vline([0,0.09*num.fps], 'k-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Walking')

subplot(1,2,2) %str stim
hold all
% stim with light
    temp.Y1 = mean(speed.Control.rotlight,2);
    temp.Y2 = mean(speed.Stim.rotlight,2);
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.errrotlight(1:end-1); speed.Stim.errrotlight];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('red'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'r', 'linewidth', 1)
% control no light
    temp.y1 = mean(speed.Control.rotnolight,2);
    temp.y2 = mean(speed.Stim.rotnolight,2);
    temp.y = [temp.y1(1:end-1); temp.y2];
    temp.x = (-60:60);
    plot(temp.x, temp.y, 'color', 'k', 'linewidth', 1)    
vline([0,0.09*num.fps], 'k-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Turning')




%% Velocity Average Combined across flies

for kk = 1:num.fly
    for tt = 1:2
        for cond = 1:num.conds
            raw = nanmean(fly(kk).(Type{tt}).speed(cond).data,2);
            speed.(Type{tt}).avg(cond).data(:,kk) = raw;
%             speed.(Type{tt}).avg(cond).err(:,kk) = raw;
        end
    end
end

for tt = 1:2
    clear raw
    idx = 0;
    for cond = 1:14
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).str(:,idx) = nanmean(raw,2);
    end
    speed.(Type{tt}).strerr = nanstd(speed.(Type{tt}).str,0,2)/sqrt(12);
    speed.(Type{tt}).stravg = nanmean(speed.(Type{tt}).str,2);
    
    clear raw
    idx = 0;
    for cond = 15:28
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).rot(:,idx) = nanmean(raw,2);
    end
    speed.(Type{tt}).roterr = nanstd(speed.(Type{tt}).rot,0,2)/sqrt(12);
    speed.(Type{tt}).rotavg = nanmean(speed.(Type{tt}).rot,2);
end
    
% Rotational Velocity

for kk = 1:num.fly
    for tt = 1:2
        for cond = 1:num.conds
            raw = nanmean(fly(kk).(Type{tt}).rotvelocity(cond).data,2);
            rotvelocity.(Type{tt}).avg(cond).data(:,kk) = raw;
        end
    end
end

for tt = 1:2
    clear raw
    idx = 0;
    for cond = 1:14
        idx = idx+1;
        raw = rotvelocity.(Type{tt}).avg(cond).data;
        
        if cond >=5 && cond <=14
            rotvelocity.(Type{tt}).str(:,idx) = nanmean(raw,2);
        else
            rotvelocity.(Type{tt}).str(:,idx) = nanmean(raw,2)*-1;
        end
       
    end
    rotvelocity.(Type{tt}).strerr = nanstd(rotvelocity.(Type{tt}).str,0,2)/sqrt(12);
    rotvelocity.(Type{tt}).stravg = nanmean(rotvelocity.(Type{tt}).str,2);
    
    clear raw
    idx = 0;
    for cond = 15:28
        idx = idx+1;
        raw = rotvelocity.(Type{tt}).avg(cond).data;
        
        if cond >=22 && cond <=28
            rotvelocity.(Type{tt}).rot(:,idx) = nanmean(raw,2);
        else
            rotvelocity.(Type{tt}).rot(:,idx) = nanmean(raw,2)*-1;
        end
    end
    rotvelocity.(Type{tt}).roterr = nanstd(rotvelocity.(Type{tt}).rot,0,2)/sqrt(12);
    rotvelocity.(Type{tt}).rotavg = nanmean(rotvelocity.(Type{tt}).rot,2);
end
       

fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
subplot(2,2,1) %str stim
hold all
% stim with light
    temp.x = (-60:60);
    temp.Y1 = speed.Control.stravg;
    temp.Y2 = speed.Stim.stravg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.strerr(1:end-1); speed.Stim.strerr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)
  
vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Speed during Walking')

subplot(2,2,2) %str stim
hold all
% stim with light
    temp.Y1 = speed.Control.rotavg;
    temp.Y2 = speed.Stim.rotavg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.roterr(1:end-1); speed.Stim.roterr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)

vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Speed during Turning')

% save_figure(fig, [figures_dir, filename(1:end-4), ' Walking and Turning Nof3'])







% fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
subplot(2,2,3) %str stim
hold all
% stim with light
    temp.x = (-60:60);
    temp.Y1 = rotvelocity.Control.stravg;
    temp.Y2 = rotvelocity.Stim.stravg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [rotvelocity.Control.strerr(1:end-1); rotvelocity.Stim.strerr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)
  
vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Rot Vel (deg/s)')
title('Rotational Velocity during Walking')

subplot(2,2,4) %str stim
hold all
% stim with light
    temp.Y1 = rotvelocity.Control.rotavg;
    temp.Y2 = rotvelocity.Stim.rotavg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [rotvelocity.Control.roterr(1:end-1); rotvelocity.Stim.roterr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)

vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Rot velocity (deg/s)')
title({'Rotational Velocity during Turning'; [fly(1).param.cross ', N = ' num2str(num.fly)]})

% save_figure(fig, [figures_dir, filename(1:end-4), ' Walking and Turning Nof3'])

save_figure(fig, [figures_dir, filename(1:end-4), ' Walking and Turning'])








