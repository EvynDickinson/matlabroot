
%% Load Fly Structure (already built with "MIDSTIM_Step_4_groupflies.m")
clc; close all; clear all

newdata = 0;
    
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
        beep
    case 0
        % LOAD FROM CREATED DATA:
        ROI = 1:60;
        fly_cross = select_cross;
        figures_dir = ['C:\matlabroot\Sweta Paper Work\' fly_cross '\'];
        load([figures_dir, fly_cross ' speed time course figure'])
        Speed = InputData;
        load([figures_dir, fly_cross ' rotvel time course figure'])
        Rotvel = InputData;
        load([figures_dir, fly_cross ' Fly starts from stationary'])
        stationary = data;
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
%
%%
data = [];


% Measuring the diff between light effect and behavior effects
% area under the curve diff between no light cond and light conditions --
% box plot or scatter for each fly
% average for each fly, then average of averages. 
for pp = 1:2
    switch pp
        case 1
            param = 'speed';
            InputData = Speed;
        case 2
            param = 'rotvelocity';
            InputData = Rotvel;
    end

    shorttime = [61:(61+round(.09*num.fps)+(.5*num.fps))];
    longtime = [61:(61+round(.72*num.fps)+(.5*num.fps))];
    timerange(1).data = shorttime;
    timerange(1).length = length(shorttime);
    timerange(2).data = longtime;
    timerange(2).length = length(longtime);
    %add data ranges to new structure
    LA(1).control = InputData(1).y(:, timerange(1).data); %short control
    LA(2).control = InputData(1).y(:, timerange(2).data); %long control
    LA(1).stim = InputData(2).y(:, timerange(1).data); %short light stim
    LA(2).stim = InputData(3).y(:, timerange(2).data); %long light stim
%     LA(1).control = InputData(1).y(timerange(1).data,:); %short control
%     LA(2).control = InputData(1).y(timerange(2).data,:); %long control
%     LA(1).stim = InputData(2).y(timerange(1).data,:); %short light stim
%     LA(2).stim = InputData(3).y(timerange(2).data,:); %long light stim
    % get stats for each structure

    for tt = 1:2
        LA(tt).controlavg = nanmean(sum(LA(tt).control));
        LA(tt).controlerr = nanstd(sum(LA(tt).control));
        LA(tt).stimavg = nanmean(sum(LA(tt).stim));
        LA(tt).stimerr = nanstd(sum(LA(tt).stim));
%         LA(tt).diff = LA(tt).stimavg/LA(tt).controlavg;
        LA(tt).diff = (LA(tt).controlavg-LA(tt).stimavg)/LA(tt).controlavg;
        LA(tt).differr = sqrt((LA(tt).controlerr/LA(tt).controlavg)^2+(LA(tt).stimerr/LA(tt).stimavg)^2);
    end

    fig = getfig;
    hold all
    for tt = 1:2
       scatter(tt, LA(tt).diff, 100, 'filled', 'k')
       errorbar(tt, LA(tt).diff , LA(tt).differr/num.fly, 'Color', 'k')
    end
    xlim([0, 3])
    ylim([-2, 2])
    xlabel('Short light, long light')
    title({structure_name; ['Change in ' param]; ''})
    ylabel(['fraction change in ' param])
    set(gca,'TickDir','out');
    figure_name = [figures_dir, structure_name, ' Percent Change in ' param];
    save(figure_name, 'LA')

    save_figure(fig, figure_name);
end

%% Create mean of means: time course figures for speed and rotational velocity & save the data for the figure

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
        ylim([0, 18])
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
            Ydata = nan(121,1);
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
                %input all data
                y = [fly(kk).Control.(param{condition})(COND).data(2:end, :);...
                         fly(kk).Stim.(param{condition})(COND).data(:, :)];
                y = y*adj;
                %set nonwalkers to NaN
                FILTER = ~filter;
                y(:,FILTER) = nan;
                Ydata = [Ydata, y];
            end 
            Yavg = nanmean(Ydata,2);
            InputData(ii).y(:,kk) = Yavg;
        end   
    %     smoothing = 2;
        InputData(ii).x = (-2:1/num.fps:2)';
        InputData(ii).Color = Color(colorchoice{ii});
        InputData(ii).xavg = nanmean(InputData(ii).x,2);
        InputData(ii).yavg = nanmean(InputData(ii).y,2).*ADJ;
        InputData(ii).yerr = (nanstd(InputData(ii).y,0,2)/sqrt(num.fly)).*ADJ;

%         InputData(ii).yerr = sem(InputData(ii).y,2).*ADJ;
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
title({structure_name; param{condition}})
save(fig_name, 'InputData', 'parameters', 'num')
save_figure(fig, fig_name);

end

clear adj ADJ adj2 ax ccw cw cond COND h ii ind idx LStyle p temp


%% Average Trajectory for walkers and non walkers:
if newdata == 1
    clear stim ctrl cw ccw
end
p.color_idx = Color('LimeGreen'); 
% p.color_idx = Color('Red'); 
BClass = {'Walking', 'Nonwalking'};
classification = 0; % selecting only for walkers
p.wide = 3;
p.thin = 1;
p.dotwide = 40;
p.dotthin = 20;
p.cap = 0;
x_lim = [-.6, 0.6]; 
y_lim = [-.1, 2];

fig = getfig('Destination Trials');


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
%                 stim(cond).x.data = [stim(cond).x.data, cw, ccw];
                stim(cond).x.data(kk,[rep, rep+3]) = [cw,ccw];
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
%                 stim(cond).y.data = [stim(cond).y.data, cw, ccw];  
                stim(cond).y.data(kk,[rep, rep+3]) = [cw,ccw];
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
%                 ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
                ctrl(cond).x.data(kk,[rep, rep+3]) = [cw,ccw];
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
%                 ctrl(cond).y.data = [ctrl(cond).y.data, cw, ccw]; 
                ctrl(cond).y.data(kk,[rep, rep+3]) = [cw,ccw];
            end   
        end    
        % Stats:
        for ax = ['x', 'y']
            stim(cond).(ax).allavg = nanmean(stim(cond).(ax).data,2);
            stim(cond).(ax).avg = nanmean(stim(cond).(ax).allavg);
            stim(cond).(ax).err = sem(stim(cond).(ax).allavg);
            
            ctrl(cond).(ax).allavg = nanmean(ctrl(cond).(ax).data,2);
            ctrl(cond).(ax).avg = nanmean(ctrl(cond).(ax).allavg);
            ctrl(cond).(ax).err = sem(ctrl(cond).(ax).allavg);
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
       
    figure(fig)
    title({['Group: ' BClass{classification+1}]; [structure_name ' Turning Cue']})  
    set(gca,'TickDir','out');
        



end

figure_name = [figures_dir, structure_name ' Average Trajectory walking selective'];
save_figure(fig, figure_name);
save(figure_name, 'stim', 'ctrl')











%%
newdata = 1;
% Stationary flies that start walking after activation|silencing
% save an individual image and then the data, then, once all four crosses
% are run, make the figure with all the appropriate data
type = {'control', 'light'};
min_speed = 0.3;
% pointrange = 4:19; % 500ms period after light offset
pointrange = 22:37; %long light
num.divide = 12;
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

%% Short light length starts from stationary
data = [];
% Stationary flies that start walking after activation|silencing
% save an individual image and then the data, then, once all four crosses
% are run, make the figure with all the appropriate data
type = {'control', 'light'};
min_speed = 0.3;
pointrange = 4:19; % 500ms period after light offset short light
% pointrange = 22:37; %long light

num.divide = 12;
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
    for cond = [4 11 18 25] %[7 14 21 28]
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
figure_name = [figures_dir, structure_name, ' Fly starts from stationary'];
save(figure_name, 'data')

save_figure(fig, figure_name);
    
    

beep




%% Measuring the diff between light effect and behavior effects
% area under the curve diff between no light cond and light conditions --
% box plot or scatter for each fly
% % average for each fly, then average of averages. 
% for pp = 1:2
%     switch pp
%         case 1
%             param = 'speed';
%             InputData = Speed;
%         case 2
%             param = 'rotvelocity';
%             InputData = Rotvel;
%     end
% 
%     shorttime = [61:(61+round(.09*num.fps)+(.5*num.fps))];
%     longtime = [61:(61+round(.72*num.fps)+(.5*num.fps))];
%     timerange(1).data = shorttime;
%     timerange(1).length = length(shorttime);
%     timerange(2).data = longtime;
%     timerange(2).length = length(longtime);
%     %add data ranges to new structure
%     LA(1).control = InputData(1).y(timerange(1).data,:); %short control
%     LA(2).control = InputData(1).y(timerange(2).data,:); %long control
%     LA(1).stim = InputData(2).y(timerange(1).data,:); %short light stim
%     LA(2).stim = InputData(3).y(timerange(2).data,:); %long light stim
%     % get stats for each structure
% 
%     for tt = 1:2
%         LA(tt).controlavg = nanmean(sum(LA(tt).control));
%         LA(tt).controlerr = nanstd(sum(LA(tt).control));
%         LA(tt).stimavg = nanmean(sum(LA(tt).stim));
%         LA(tt).stimerr = nanstd(sum(LA(tt).stim));
%         LA(tt).diff = LA(tt).stimavg/LA(tt).controlavg;
%         LA(tt).differr = sqrt((LA(tt).controlerr/LA(tt).controlavg)^2+(LA(tt).stimerr/LA(tt).stimavg)^2);
%     end
% 
%     fig = getfig;
%     hold all
%     for tt = 1:2
%        scatter(tt, LA(tt).diff, 100, 'filled', 'k')
%        errorbar(tt, LA(tt).diff , LA(tt).differr, 'Color', 'k')
%     end
%     xlim([0, 3])
%     ylim([-1, 4])
%     xlabel('Short light, long light')
%     title({structure_name; ['Change in ' param]; ''})
%     ylabel(['fraction change in ' param])
%     set(gca,'TickDir','out');
%     figure_name = [figures_dir, structure_name, ' Change in ' param];
%     save(figure_name, 'LA')
% 
%     save_figure(fig, figure_name);
% end








