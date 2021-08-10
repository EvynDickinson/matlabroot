

%% Load Fly Structure (already built with "MIDSTIM_Step_4_groupflies.m")
clc; close all; clear all

newdata = 1;

switch newdata
    case 1 %create a new data set completely
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
        parameters = fly(1).param;
        fprintf(['\n Loaded ' filename '\n'])
        setGlobalx(parameters)
        % load labels & condition parameters
        [num, Fictrac, Matlab, labels] = NUM(fly(1).param);  
        num.fly = length(fly);
        Type = {'Stim', 'Control'};
        % Create a Figures Folder for this structure:
        figures_dir = ['C:\matlabroot\Sweta Paper Work\' filename(1:end-4) '\'];
        if ~isfolder(figures_dir)
            mkdir(figures_dir)
        end 
        % ELIMINATE OUTLIERS
        fly = eliminate_speed_outliers(fly, 5);
        % LOAD BEHAVIOR CLASSIFICATION DATA:
        load([directory, '\behavior class\', structure_name,  ' behavior class'])
        % TRAJECTORY PATHS AND HEATMAPS
        load(['C:\matlabroot\Interneuron Lines\', structure_name, '\Trajectory Data'])

    case 0
        % LOAD FROM CREATED DATA:
        ROI = 1:60;
        fly_cross = select_cross;
        figures_dir = ['C:\matlabroot\Sweta Paper Work\' fly_cross '\'];
        load([figures_dir, fly_cross ' speed time course figure'])
        Speed = InputData;
        load([figures_dir, fly_cross ' rotvel time course figure'])
        Rotvel = InputData;
        load([figures_dir, fly_cross ' Fly starts from stationary 720ms'])
        load([figures_dir, fly_cross ' Average Trajectory walking selective'])
        load(['C:\matlabroot\Interneuron Lines\', fly_cross, '\Trajectory Data'])
        load(['C:\matlabroot\behavior class\', fly_cross,  ' behavior class'])
        structure_name = fly_cross;
end

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


%% Create time course figures for speed and rotational velocity & save the data for the figure

param = {'speed', 'rotvelocity'};
for condition = 1:2 %speed|rotvel
    fig = getfig('', 1); hold all
    colorchoice = {'black', 'blue', 'red'};
    ii = 0;
    if condition == 1
        condchoice = [1,4,7];
        ADJ = 10; %convert to mm from cm
        ylabel({'Speed (mm/sec)'; ''})
        fig_name = [figures_dir, structure_name ' speed time course figure'];
        ylim([0, 15])
    else
        condchoice = [15,18,21];
        ADJ = 1;
        ylabel({'Rotational Velocity (deg/sec)'; ''})
        fig_name = [figures_dir, structure_name ' rotvel time course figure'];
        ylim([-75, 300])
    end
for icond = condchoice
    
    ii = ii+1;
    if newdata == 1 %loading new data
    InputData(ii).x = [];
    InputData(ii).y = [];
    InputData(ii).Ind = 1;
    for kk = 1:num.fly
        for tt = 1:2 %CW & CCW
            switch tt
                case 1
                    COND = icond;
                     if condition == 2
                        adj = -1; %change this for rot velocity
                    else 
                        adj = 1; 
                     end
                case 2
                    COND = icond+7;
                   adj = 1;
            end
            filter = group(kk).walking(COND,:);
            num.present = sum(filter)-1;
            if num.present < 0
            else
            y = [fly(kk).Control.(param{condition})(COND).data(2:end, filter);...
                     fly(kk).Stim.(param{condition})(COND).data(:, filter)];

            InputData(ii).y(InputData(ii).Ind:InputData(ii).Ind+num.present,:) = y'.*adj;
            InputData(ii).Ind = InputData(ii).Ind+num.present+1;
            end  
        end  
    end   
%     smoothing = 2;
    InputData(ii).x = (-2:1/num.fps:2)';
    InputData(ii).Color = Color(colorchoice{ii});
    InputData(ii).xavg = nanmean(InputData(ii).x,2);
    InputData(ii).yavg = nanmean(InputData(ii).y,1).*ADJ;
    InputData(ii).yerr = (nanstd(InputData(ii).y,1)/sqrt(num.fly)).*ADJ;
%         InputData(ii).yerr = sem(InputData(ii).y,1).*ADJ;
    InputData(ii).ymin = InputData(ii).yavg-InputData(ii).yerr;
    InputData(ii).ymax = InputData(ii).yavg+InputData(ii).yerr;
    else %using old data
        switch condition
            case 1 %speed
                InputData(ii) = Speed(ii);
            case 2 %rotvel
                InputData(ii) = Rotvel(ii);
        end
        
    end
    plot(InputData(ii).xavg, InputData(ii).yavg, 'LineWidth', 2, 'Color', InputData(ii).Color)
    plot(InputData(ii).xavg, InputData(ii).ymin, 'LineWidth', 1, 'Color', InputData(ii).Color)
    plot(InputData(ii).xavg, InputData(ii).ymax, 'LineWidth', 1, 'Color', InputData(ii).Color)
    set(gca,'TickDir','out');
    
end

vline([0, 0.09, 0.72], 'y')
xlabel({'';'Time (sec)' })
title({structure_name; param{condition};  'Blue short|black long'})
save(fig_name, 'InputData', 'parameters', 'num')
save_figure(fig, fig_name);

end

% clear adj ADJ adj2 ax ccw cw cond COND h ii ind idx LStyle p temp


%% Average Trajectory for walkers and non walkers: (finished)
if newdata == 1
    clear stim ctrl cw ccw
end
p.color_idx = Color('LimeGreen'); 
% p.color_idx = Color('Red'); 
BClass = {'Walking', 'Nonwalking'};
for classification = 0:0 %1 = selecting only for non walkers; 0=walkers
    p.wide = 3;
    p.thin = 1;
    p.dotwide = 40;
    p.dotthin = 20;
    p.cap = 0;
    x_lim = [-.6, 0.6]; 
    y_lim = [-.1, 2];

    fig = getfig('Destination Trials');
    % STIM LENGTH TRIAL:

%     for cond = 2:7 %WALKING
%         idx = cond-1;
%         stim(cond).x.data = []; stim(cond).y.data = [];
%         ctrl(cond).x.data = []; ctrl(cond).y.data = [];
%         ROI = ceil(parameters.light_length(cond)*parameters.fps);
%         for kk = 1:num.fly % concatenate all the data into one vector
%             for rep = 1:num.reps
%                 % stim X
%                     cw = traj(cond).x(kk, rep).data(ROI);
%                     if group(kk).walking(cond,rep)==classification
%     %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
%                         cw = NaN;
%                     end
%                     ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
%                     if group(kk).walking(cond+7,rep)==classification
%     %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
%                         ccw = NaN;
%                     end
%                 stim(cond).x.data = [stim(cond).x.data, cw, ccw];
%                 % stim Y
%                     cw = traj(cond).y(kk, rep).data(ROI);
%                     if group(kk).walking(cond,rep)==classification
%     %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
%                         cw = NaN;
%                     end
%                     ccw = traj(cond+7).y(kk, rep).data(ROI);
%                     if group(kk).walking(cond+7,rep)==classification
%     %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
%                         ccw = NaN;
%                     end
%                 stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
%                 % control X
%                 cntrlcond = 1;
%                     cw = traj(cntrlcond).x(kk, rep).data(ROI);
%                     if group(kk).walking(cntrlcond,rep)==classification
%     %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
%                         cw = NaN;
%                     end
%                     ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
%                     if group(kk).walking(cntrlcond+7,rep)==classification
%     %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
%                         ccw = NaN;
%                     end
%                 ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
%                 % control Y
%                     cw = traj(cntrlcond).y(kk, rep).data(ROI);
%                     if group(kk).walking(cntrlcond,rep)==classification
%     %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
%                         cw = NaN;
%                     end
%                     ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
%                     if group(kk).walking(cntrlcond+7,rep)==classification
%     %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
%                         ccw = NaN;
%                     end
%                 ctrl(cond).y.data = [ctrl(cond).y.data, cw, ccw]; 
%             end    
%         end    
%         % Stats:
%         for ax = ['x', 'y']
%             stim(cond).(ax).avg = nanmean(stim(cond).(ax).data);
%             ctrl(cond).(ax).avg = nanmean(ctrl(cond).(ax).data); 
%             stim(cond).(ax).err = sem(stim(cond).(ax).data,2);
%             ctrl(cond).(ax).err = sem(ctrl(cond).(ax).data,2);
%         end       
%         % PLOT THE DATA:     
%         figure(fig); subplot(1,2,1); hold all
%         %stim data -- wide bars and point in light dependent color
%            errorbar(stim(cond).x.avg, stim(cond).y.avg,... %coordinates
%                  stim(cond).y.err, stim(cond).y.err,... %vertical error bar
%                  stim(cond).x.err, stim(cond).x.err,... %horizontal error bar
%                  'LineStyle','none', 'color', p.color_idx,...
%                  'LineWidth', p.wide, 'CapSize', p.cap);
%            scatter(stim(cond).x.avg, stim(cond).y.avg, p.dotwide, p.color_idx, 'filled')
%         %control data -- wide bars|dot first in black and then thin in light dependent color
%                  % WIDE
%            errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
%                  ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
%                  ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
%                  'LineStyle','none', 'color', 'k',...
%                  'LineWidth', p.wide, 'CapSize', p.cap);
%            scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotwide, 'k', 'filled')
%     %              % THIN
%     %        errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
%     %              ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
%     %              ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
%     %              'LineStyle','none', 'color', p.color_idx(cond-1,:),...
%     %              'LineWidth', p.thin, 'CapSize', p.cap);
%     %        scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotthin, p.color_idx(idx,:), 'filled')     
% 
%         xlim(x_lim); ylim(y_lim)
%         vline(0, 'r-');hline(0, 'r-')
%         getaxes(fig, 10);
%     end
%     
    for cond = 16:21 %TURNING
        if newdata == 1
            idx = cond-15;
            stim(cond).x.data = []; stim(cond).y.data = [];
            ctrl(cond).x.data = []; ctrl(cond).y.data = [];
            ROI = ceil(parameters.light_length(idx+1)*parameters.fps);
            for kk = 1:num.fly % concatenate all the data into one vector
                for rep = 1:num.reps
                    % stim X
                        cw = traj(cond).x(kk, rep).data(ROI);
                        if group(kk).walking(cond,rep)==classification
        %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                            cw = NaN;
                        end
                        ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
                        if group(kk).walking(cond+7,rep)==classification
        %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                            ccw = NaN;
                        end
                    stim(cond).x.data = [stim(cond).x.data, cw, ccw];
                    % stim Y
                        cw = traj(cond).y(kk, rep).data(ROI);
                        if group(kk).walking(cond,rep)==classification
        %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                            cw = NaN;
                        end
                        ccw = traj(cond+7).y(kk, rep).data(ROI);
                        if group(kk).walking(cond+7,rep)==classification
        %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                            ccw = NaN;
                        end
                    stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
                    % control X
                        cntrlcond = 15;
                        cw = traj(cntrlcond).x(kk, rep).data(ROI);
                        if group(kk).walking(cntrlcond,rep)==classification
        %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                            cw = NaN;
                        end
                        ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
                        if group(kk).walking(cntrlcond+7,rep)==classification
        %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                            ccw = NaN;
                        end
                    ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
                    % control Y
                        cw = traj(cntrlcond).y(kk, rep).data(ROI);
                        if group(kk).walking(cntrlcond,rep)==classification
        %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                            cw = NaN;
                        end
                        ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
                        if group(kk).walking(cntrlcond+7,rep)==classification
        %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                            ccw = NaN;
                        end
                    ctrl(cond).y.data = [ctrl(cond).y.data, cw, ccw]; 
                end    
            end    
            % Stats:
            for ax = ['x', 'y']
                stim(cond).(ax).avg = nanmean(stim(cond).(ax).data);
                ctrl(cond).(ax).avg = nanmean(ctrl(cond).(ax).data); 
                stim(cond).(ax).err = sem(stim(cond).(ax).data,2);
                ctrl(cond).(ax).err = sem(ctrl(cond).(ax).data,2);
            end       
        % PLOT THE DATA:     
        end
        
        figure(fig); 
%         subplot(1,2,2); 
        hold all
        %stim data -- wide bars and point in light dependent color
           errorbar(stim(cond).x.avg, stim(cond).y.avg,... %coordinates
                 stim(cond).y.err, stim(cond).y.err,... %vertical error bar
                 stim(cond).x.err, stim(cond).x.err,... %horizontal error bar
                 'LineStyle','none', 'color', p.color_idx,...
                 'LineWidth', p.wide, 'CapSize', p.cap);
           scatter(stim(cond).x.avg, stim(cond).y.avg, p.dotwide, p.color_idx, 'filled')
        %control data -- wide bars|dot first in black and then thin in light dependent color
                 % WIDE
           errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
                 ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
                 ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
                 'LineStyle','none', 'color', 'k',...
                 'LineWidth', p.wide, 'CapSize', p.cap);
           scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotwide, 'k', 'filled')

        xlim(x_lim); ylim(y_lim)
        vline(0, 'r-');hline(0, 'r-')
        getaxes(fig, 10);
       
    end

    figure(fig)
%     subplot(1,2,1)
%         title({parameters.cross; 'Walking Cue'})
%         xlabel('position (cm)'); ylabel('position (cm)')
%     subplot(1,2,2)
        title({['Group: ' BClass{classification+1}]; [structure_name ' Turning Cue']})  
        set(gca,'TickDir','out');
        
figure_name = [figures_dir, structure_name ' Average Trajectory walking selective'];
    switch classification
        case 0
        save_figure(fig, figure_name);
        case 1
        save_figure(fig, [figures_dir, structure_name ' Average Trajectory nonwalking selective']);
    end
end
save(figure_name, 'stim', 'ctrl')


%% Stationary flies that start walking after activation|silencing
% save an individual image and then the data, then, once all four crosses
% are run, make the figure with all the appropriate data
min_speed = 0.3;
pointrange = 6:12; % 200ms period after 90ms post light length
num.divide = 6;

if newdata == 1
    for kk = 1:num.fly

    %control information
    data.control.speed(kk,1:num.divide) = nan; %fill the trial spots with NaN since moving flies won't count 

    idx = 0;
    for cond = [1,8] %control stimuli15,22
        for rep = 1:num.reps
            idx = idx + 1;
            if group(kk).stationary(cond,rep) == true
                data.control.speed(kk,idx) = mean(fly(kk).Stim.speed(cond).data(pointrange, rep));
            end
        end
    end
    data.control.stationary(kk) = (sum(~isnan(data.control.speed(kk,:))));
    data.control.idx(kk,:) = data.control.speed(kk,:)>=min_speed;
    data.control.fraction(kk) = sum(data.control.idx(kk,:))/data.control.stationary(kk)*100;
    data.control.increase(kk) = nanmean(data.control.speed(kk,data.control.idx(kk,:)));



    % stimulus information
    data.light.speed(kk,1:num.divide) = nan;
    idx = 0;
    for cond = [4,11] %,18,25
        for rep = 1:num.reps
            idx = idx + 1;
            if group(kk).stationary(cond,rep) == true
                data.light.speed(kk,idx) = mean(fly(kk).Stim.speed(cond).data(pointrange, rep));
            end
        end
    end
    data.light.stationary(kk) = (sum(~isnan(data.light.speed(kk,:))));
    data.light.idx(kk,:) = data.light.speed(kk,:)>=min_speed;
    data.light.fraction(kk) = sum(data.light.idx(kk,:))/data.light.stationary(kk)*100;
    data.light.increase(kk) = nanmean(data.light.speed(kk,data.light.idx(kk,:)));


    data.plot(kk,1) = data.control.fraction(kk);
    data.plot(kk,2) = data.light.fraction(kk);

    end

    % averages:
    data.control.fractionavg = nanmean(data.control.fraction);
    data.control.fractionerr = nanstd(data.control.fraction)/sqrt(num.fly);
    data.light.fractionavg = nanmean(data.light.fraction);
    data.light.fractionerr = nanstd(data.light.fraction)/sqrt(num.fly);

    data.control.increaseavg = nanmean(data.control.increase);
    data.control.increaseerr = nanstd(data.control.increase)/sqrt(num.fly);
    data.light.increaseavg = nanmean(data.light.increase);
    data.light.increaseerr = nanstd(data.light.increase)/sqrt(num.fly);
end

% data.light.increase(kk)
x1 = 0.5:(0.25/num.fly):0.75;
x2 = 1:(0.25/num.fly):1.25;
point_size = 100 ;
x3 = x1(end)+(0.25/num.fly*2);
x4 = x2(end)+(0.25/num.fly*2);

% plot the fraction of flies that started moving...
fig = getfig;
subplot(1,2,1); hold all
for kk = 1:num.fly
    plot([x1(kk), x2(kk)], data.plot(kk,:), 'color', Color('grey'))
    scatter(x1(kk),data.control.fraction(kk),point_size, 'filled', 'markerfacecolor', 'k')
    scatter(x2(kk),data.light.fraction(kk),point_size, 'filled', 'markerfacecolor', 'b')
end

errorbar(x3,(data.control.fractionavg),(data.control.fractionerr), 'LineWidth', 2, 'Color', 'g')
scatter(x3, data.control.fractionavg, point_size, 's', 'filled', 'markerfacecolor', 'g')
errorbar(x4,(data.light.fractionavg),(data.light.fractionerr), 'LineWidth', 2, 'Color', 'g')
scatter(x4, data.light.fractionavg, point_size, 's', 'filled', 'markerfacecolor', 'g')

xlim([0 2])
ylim([0,100])
title({'Fraction of times a fly started moving'; ''})
ylabel('% of trials')
xlabel('Control   |     0.09sec')
set(gca,'xticklabel',{[]})

% plot the increase in speed for those that started moving...
subplot(1,2,2)
hold all
for kk = 1:num.fly
    scatter(x1(kk),data.control.increase(kk)*10,point_size, 'filled', 'markerfacecolor', 'k')
    scatter(x2(kk),data.light.increase(kk)*10,point_size, 'filled', 'markerfacecolor', 'b')
end
errorbar(x3,(data.control.increaseavg*10),(data.control.increaseerr*10), 'LineWidth', 2, 'Color', 'g')
scatter(x3, data.control.increaseavg*10, point_size, 's', 'filled', 'markerfacecolor', 'g')
errorbar(x4,(data.light.increaseavg*10),(data.light.increaseerr*10), 'LineWidth', 2, 'Color', 'g')
scatter(x4, data.light.increaseavg*10, point_size, 's', 'filled', 'markerfacecolor', 'g')

xlim([0 2.5])
ylim([0,10])
title({'Average Speed increase following stimulus'; ''})
ylabel('Speed increase (mm/sec)')
xlabel('Control   |     0.09sec')
set(gca,'xticklabel',{[]})
set(gca,'TickDir','out');


figure_name = [figures_dir, structure_name, ' Fly starts from stationary'];
save(figure_name, 'data')

save_figure(fig, figure_name);

fprintf('\n DONE!')



%%
data = [];
newdata = 1;
% Stationary flies that start walking after activation|silencing
% save an individual image and then the data, then, once all four crosses
% are run, make the figure with all the appropriate data
type = {'control', 'light'};
min_speed = 0.3;
% pointrange = 4:19; % 500ms period after light offset
pointrange = 22:37; %long light
num.divide = 12; %(4 conds x 3 reps)
ind = 0;
if newdata == 1
    for kk = 1:num.fly
    %control information
    data.control.speed(kk,1:num.divide) = nan; %fill the trial spots with NaN since moving flies won't count 

    idx = 0;
    for cond = [1,8, 15,22] %control stimuli
        for rep = 1:num.reps
            idx = idx + 1;
            if group(kk).stationary(cond,rep) == true
                data.control.speed(kk,idx) = mean(fly(kk).Stim.speed(cond).data(pointrange, rep));
            end
        end
    end
    data.control.stationary(kk) = (sum(~isnan(data.control.speed(kk,:))));
    data.control.idx(kk,:) = data.control.speed(kk,:)>=min_speed;
    
    data.control.fraction(kk) = sum(data.control.idx(kk,:))/data.control.stationary(kk)*100;
    data.control.increase(kk) = nanmean(data.control.speed(kk,data.control.idx(kk,:)));



    % stimulus information
    data.light.speed(kk,1:num.divide) = nan;
    idx = 0;
    for cond = [7 14 21 28] %
        for rep = 1:num.reps
            idx = idx + 1;
            if group(kk).stationary(cond,rep) == true
                data.light.speed(kk,idx) = mean(fly(kk).Stim.speed(cond).data(pointrange, rep));
            end
        end
    end
    data.light.stationary(kk) = (sum(~isnan(data.light.speed(kk,:))));
    data.light.idx(kk,:) = data.light.speed(kk,:)>=min_speed;
    data.light.fraction(kk) = sum(data.light.idx(kk,:))/data.light.stationary(kk)*100;
    data.light.increase(kk) = nanmean(data.light.speed(kk,data.light.idx(kk,:)));


    data.plot(kk,1) = data.control.fraction(kk);
    data.plot(kk,2) = data.light.fraction(kk);

    end
    
    for tt = 1:2
        data.(type{tt}).totalstationary = sum(data.(type{tt}).stationary);
        data.(type{tt}).totalmoved = sum(sum(data.(type{tt}).idx));
        data.(type{tt}).percent = data.(type{tt}).totalmoved/data.(type{tt}).totalstationary*100;
        data.bargraph(tt) = data.(type{tt}).percent;
    end

end

fig = getfig;
bar([1,2], data.bargraph, 'BarWidth', 1)
ylim([0,100])
xlabel(['total stationary: ' num2str(data.control.totalstationary),...
        ' vs ' num2str(data.light.totalstationary)])
title(structure_name)
set(gca,'TickDir','out');
figure_name = [figures_dir, structure_name, ' Fly starts from stationary 720ms'];
save(figure_name, 'data')

save_figure(fig, figure_name);
    
    

beep





%% Create mean of means: time course figures for speed and rotational velocity & save the data for the figure
% 
% param = {'speed', 'rotvelocity'};
% for condition = 1:2 %speed|rotvel
%     fig = getfig('', 1); hold all
%     colorchoice = {'black', 'blue', 'red'};
%     ii = 0;
%     if condition == 1
%         condchoice = [1,4,7];
%         ADJ = 10; %convert to mm from cm
%         ylabel({'Speed (mm/sec)'; ''})
%         fig_name = [figures_dir, structure_name ' speed time course figure'];
%         ylim([0, 18])
%     else
%         condchoice = [15,18,21];
%         ADJ = 1;
%         ylabel({'Rotational Velocity (deg/sec)'; ''})
%         fig_name = [figures_dir, structure_name ' rotvel time course figure'];
%         ylim([-75, 300])
%     end
% for icond = condchoice
%     ii = ii+1;
%     if newdata == 1 %loading new data
%         InputData(ii).x = [];
%         InputData(ii).y = [];
%         InputData(ii).Ind = 1;
%         for kk = 1:num.fly
%             Ydata = nan(121,1);
%             for tt = 1:2 %CW & CCW
%                 switch tt
%                     case 1
%                         COND = icond;
%                          if condition == 2
%                             adj = -1; %change this for rot velocity
%                         else 
%                             adj = 1; 
%                          end
%                     case 2
%                         COND = icond+7;
%                        adj = 1;
%                 end
%                 filter = group(kk).walking(COND,:);
%                 %input all data
%                 y = [fly(kk).Control.(param{condition})(COND).data(2:end, :);...
%                          fly(kk).Stim.(param{condition})(COND).data(:, :)];
%                 y = y*adj;
%                 %set nonwalkers to NaN
%                 FILTER = ~filter;
%                 y(:,FILTER) = nan;
%                 Ydata = [Ydata, y];
%             end 
%             Yavg = nanmean(Ydata,2);
%             InputData(ii).y(:,kk) = Yavg;
%         end   
%     %     smoothing = 2;
%         InputData(ii).x = (-2:1/num.fps:2)';
%         InputData(ii).Color = Color(colorchoice{ii});
%         InputData(ii).xavg = nanmean(InputData(ii).x,2);
%         InputData(ii).yavg = nanmean(InputData(ii).y,2).*ADJ;
%         InputData(ii).yerr = (nanstd(InputData(ii).y,0,2)/sqrt(num.fly)).*ADJ;
% 
% %         InputData(ii).yerr = sem(InputData(ii).y,2).*ADJ;
%         InputData(ii).ymin = InputData(ii).yavg-InputData(ii).yerr;
%         InputData(ii).ymax = InputData(ii).yavg+InputData(ii).yerr;
%     else %using old data
%         switch condition
%             case 1 %speed
%                 InputData(ii) = Speed(ii);
%             case 2 %rotvel
%                 InputData(ii) = Rotvel(ii);
%         end
%     end
%     plot(InputData(ii).xavg, InputData(ii).yavg, 'LineWidth', 2, 'Color', InputData(ii).Color)
%     plot(InputData(ii).xavg, InputData(ii).ymin, 'LineWidth', 1, 'Color', InputData(ii).Color)
%     plot(InputData(ii).xavg, InputData(ii).ymax, 'LineWidth', 1, 'Color', InputData(ii).Color)
%     set(gca,'TickDir','out');
%     
% end
% 
% vline([0, 0.09, 0.72], 'y')
% xlabel({'';'Time (sec)' })
% title({structure_name; param{condition}})
% save(fig_name, 'InputData', 'parameters', 'num')
% save_figure(fig, fig_name);
% 
% end
% 
% % clear adj ADJ adj2 ax ccw cw cond COND h ii ind idx LStyle p temp