

%% How to align and overlay a series of steps given a stimulus start
clearvars('-except',initial_vars{:})
% Try first with trials that are in stance at the oneset of the ...

% input: 
control_cond = [1,8];   %control conditions
stim_cond = [7,14];     %opto conditions 720ms
% stim_cond = [4, 11];     %opto conditions 90ms
ifly = 8;               %fly number


% find the desired trials with swing-stance fitting: 
[stimList, controlList] = deal([]);
idx = 0;
for ii = 1:length(control_cond)
    cond = control_cond(ii);
    for rep = 1:num.reps
        if ~isempty(angles(ifly).ball.Center{cond,rep})
            idx = idx+1;
            controlList(idx,:) = [ifly,rep,cond];
        end
    end
end
C_tot = size(controlList,1);
idx = 0;
for ii = 1:length(stim_cond)
    cond = stim_cond(ii);
    for rep = 1:num.reps
        if ~isempty(angles(ifly).ball.Center{cond,rep})
            idx = idx+1;
            stimList(idx,:) = [ifly,rep,cond];
        end
    end
end
S_tot = size(stimList,1);


%% Pull a specific group of joint traces for a given fly:
% clearvars('-except',initial_vars{:})

flyList = stimList;
% flyList = controlList;

% compare joint angles for a specific group:
splt = [1:2:6, 2:2:6];

fig = getfig('',1);
% group the data:
for leg = 1:6 % switch through legs
  for iJ = 2 % switch through joints
    for n = 1:size(flyList,1) % go through all the trials
        ifly = flyList(n,1); 
        rep = flyList(n,2);
        cond = flyList(n,3);
        x = (-param.basler_delay:1/fps:param.basler_length-param.basler_delay);
        light_on = 0;
        light_Loc = round((param.basler_delay)*fps);
        light_off = (param.conds_matrix(cond).opto);
        data(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        % are they in stance or swing at the stim start?
        SwSt = angles(ifly).stance(cond,rep).loc(leg,:);
        stance(n) = SwSt(light_Loc);
    end
    %change in joint angle:
    del_data = data-data(light_Loc,:);   
    % plot the traces and Avg trace:
    subplot(3,2,splt(leg)); hold all
    plot(x, data, 'color', Color('grey'), 'linewidth', 0.5, 'linestyle', ':')
    plot(x, mean(data,2), 'color', Color('black'), 'linewidth', 1)
    % labels: 
    ylim([0,180])
    vline([light_on, light_off], 'g')
  end
end


%% pull the avg joint angle during the control period & stimulus period:

[S_data, C_data, Stim, Control] = deal([]);
% ---input---
iJ = 2;
leg = 6;
disp(['Leg ' num2str(leg)])

% stim data
for n = 1:S_tot
    ifly = stimList(n,1); 
    rep = stimList(n,2);
    cond = stimList(n,3);
    Stim(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
end

% control:
for n = 1:C_tot
    ifly = controlList(n,1); 
    rep = controlList(n,2);
    cond = controlList(n,3);
    Control(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
end

% ROIs:
light_on = round((param.basler_delay)*fps);
light_off = round((param.basler_delay+param.conds_matrix(stimList(1,3)).opto)*fps);
preROI = 1:light_on;
lazROI = light_on+1:light_off;

% find the avg joint angles: 
S_data(:,1) = mean(Stim(preROI,:)); % stim control period
S_data(:,2) = mean(Stim(lazROI,:)); % stim laser period

C_data(:,1) = mean(Control(preROI,:)); % stim control period
C_data(:,2) = mean(Control(lazROI,:)); % stim laser period

% plot the avg angle for the two time periods:
SZ = 50;

fig = getfig('',1);
subplot(1,2,1) % control plots
hold on
scatter(ones(1,C_tot),C_data(:,1), SZ, 'k') % pre
scatter(2*ones(1,C_tot),C_data(:,2), SZ, 'k') % stim
x1 = 1;
x2 = 2;
y1 = C_data(:,1);
y2 = C_data(:,2);
for ii = 1:C_tot
    plot([x1,x2], [y1(ii), y2(ii)], 'color', 'k', 'linewidth', 1)
end

c = Color('teal');
scatter(ones(1,S_tot),S_data(:,1), SZ, c) % pre
scatter(2*ones(1,S_tot),S_data(:,2), SZ, c) % stim
x1 = 1;
x2 = 2;
y1 = S_data(:,1);
y2 = S_data(:,2);
for ii = 1:C_tot
    plot([x1,x2], [y1(ii), y2(ii)], 'color', c, 'linewidth', 1)
end
xlim([0,3])
xlabel('control --- stim')
ylabel(['Avg ' Joints{iJ} ' angle (\circ)'])

set1 = [num2str(param.conds_matrix(controlList(1,3)).opto) ' ms'];
set2 = [num2str(param.conds_matrix(stimList(1,3)).opto) ' ms'];
fig_title = [FilePath.locations{ifly,4} ': Opto ' set1 ' VS ' set2];
fig_title = strrep(fig_title, '_', '-');
title(fig_title)

% run a paired t-test on the before and during data:
fprintf('\n Joint angle stats: \n')
[~,p] = ttest(S_data(:,1), S_data(:,2));
disp(['Stim data p-val: ' num2str(p)])
[~,p] = ttest(C_data(:,1), C_data(:,2));
disp(['Control data p-val: ' num2str(p)])
if S_tot==C_tot
    [~,p] = ttest(C_data(:,1), S_data(:,1));
    disp(['Control vs stim preopto p-val: ' num2str(p)])
end

% Step frequency paired with above data:
[S_data, C_data, Stim, Control] = deal([]);
% ----------------------------- joint angles

% stim data
for n = 1:S_tot
    ifly = stimList(n,1); 
    rep = stimList(n,2);
    cond = stimList(n,3);
    Stim(:,n) = gait(ifly).Freq(cond,rep).f_avg;
end

% control:
for n = 1:C_tot
    ifly = controlList(n,1); 
    rep = controlList(n,2);
    cond = controlList(n,3);
    Control(:,n) =gait(ifly).Freq(cond,rep).f_avg;
end

% ROIs:
light_on = round((param.basler_delay)*fps);
light_off = round((param.basler_delay+param.conds_matrix(stimList(1,3)).opto)*fps);
preROI = 1:light_on;
lazROI = light_on+1:light_off;

% find the avg joint angles: 
S_data(:,1) = mean(Stim(preROI,:)); % stim control period
S_data(:,2) = mean(Stim(lazROI,:)); % stim laser period

C_data(:,1) = mean(Control(preROI,:)); % stim control period
C_data(:,2) = mean(Control(lazROI,:)); % stim laser period

% plot the avg angle for the two time periods:
SZ = 50;
figure(fig) % activate the figure from above
subplot(1,2,2) % control plots
hold on
scatter(ones(1,C_tot),C_data(:,1), SZ, 'k') % pre
scatter(2*ones(1,C_tot),C_data(:,2), SZ, 'k') % stim
x1 = 1;
x2 = 2;
y1 = C_data(:,1);
y2 = C_data(:,2);
for ii = 1:C_tot
    plot([x1,x2], [y1(ii), y2(ii)], 'color', 'k', 'linewidth', 1)
end

c = Color('teal');
scatter(ones(1,S_tot),S_data(:,1), SZ, c) % pre
scatter(2*ones(1,S_tot),S_data(:,2), SZ, c) % stim
x1 = 1;
x2 = 2;
y1 = S_data(:,1);
y2 = S_data(:,2);
for ii = 1:C_tot
    plot([x1,x2], [y1(ii), y2(ii)], 'color', c, 'linewidth', 1)
end
xlim([0,3])
xlabel('control --- stim')
ylabel('Step frequency (Hz)')

set1 = [num2str(param.conds_matrix(controlList(1,3)).opto) ' ms'];
set2 = [num2str(param.conds_matrix(stimList(1,3)).opto) ' ms'];
fig_title = [FilePath.locations{ifly,4} ': Opto ' set1 ' VS ' set2];
fig_title = strrep(fig_title, '_', '-');


save_figure(fig, [fig_dir, '\' FilePath.locations{ifly,4} '\' set1 ' VS ' set2 ' angle and frequency']);

% run a paired t-test on the before and during data:
fprintf('\n Step frequency stats: \n')
[~,p] = ttest(S_data(:,1), S_data(:,2));
disp(['Stim data p-val: ' num2str(p)])
[~,p] = ttest(C_data(:,1), C_data(:,2));
disp(['Control data p-val: ' num2str(p)])
if S_tot==C_tot
    [~,p] = ttest(C_data(:,1), S_data(:,1));
    disp(['Control vs stim preopto p-val: ' num2str(p)])
end

%% WORKING HERE -- FIND THE AVG STEP FREQUENCY FOR A GIVEN RANGE AND PLOT POINTS OF CHANGE


% Step frequency paired with above data:
[S_data, C_data, Stim, Control] = deal([]);
% ----------------------------- joint angles

% stim data
for n = 1:S_tot
    ifly = stimList(n,1); 
    rep = stimList(n,2);
    cond = stimList(n,3);
    Stim(:,n) = gait(ifly).Freq(cond,rep).f_avg;
end

% control:
for n = 1:C_tot
    ifly = controlList(n,1); 
    rep = controlList(n,2);
    cond = controlList(n,3);
    Control(:,n) =gait(ifly).Freq(cond,rep).f_avg;
end

% ROIs:
light_on = round((param.basler_delay)*fps);
light_off = round((param.basler_delay+param.conds_matrix(stimList(1,3)).opto)*fps);
preROI = 1:light_on;
lazROI = light_on+1:light_off;


% avg stp frequency 2nd have ...


%% Compare the avg joint angle across legs:

[S_data, C_data, Stim, Control] = deal([]);
% ---input---
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
iJ = 2; % joint

% ROIs:
light_on = round((param.basler_delay)*fps);
light_off = round((param.basler_delay+param.conds_matrix(stimList(1,3)).opto)*fps);
preROI = 1:light_on;
lazROI = light_on+1:light_off;


fig = getfig('',1); hold on
set(fig, 'color', 'k')
for leg = 1:6

    % stim data
    for n = 1:S_tot
        ifly = stimList(n,1); 
        rep = stimList(n,2);
        cond = stimList(n,3);
        Stim(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
    end
    % control:
    for n = 1:C_tot
        ifly = controlList(n,1); 
        rep = controlList(n,2);
        cond = controlList(n,3);
        Control(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
    end


    % find the avg joint angles: 
    S_data(:,1) = mean(Stim(preROI,:)); % stim control period
    S_data(:,2) = mean(Stim(lazROI,:)); % stim laser period

    C_data(:,1) = mean(Control(preROI,:)); % stim control period
    C_data(:,2) = mean(Control(lazROI,:)); % stim laser period

    % plot the average change in joint angle for control conditions:
    data = C_data(:,2)-C_data(:,1); % change in avg joint angle
    scatter((leg)*ones(length(data),1), data, 70, 'w', 'filled')
    % plot the avg line:
    plot([leg-0.2,leg+0.2], [median(data),median(data)], 'color', 'w', 'linewidth',3)
    
    % plot the average change in joint angle for stimulus conditions
    data = S_data(:,2)-S_data(:,1); % change in avg joint angle
    scatter((leg+0.2)*ones(length(data),1), data, 70, Color(leg_colors{leg}), 'filled')
    % plot the avg line:
    plot([leg,leg+0.4], [median(data),median(data)], 'color', Color(leg_colors{leg}), 'linewidth',3)
end

% labels etc
xlim([0,7])
hline(0, 'w:')
set(gca, 'color', 'k', 'YColor', 'w', 'XColor', 'w', 'TickDir', 'out')
xlabel('leg')
ylabel('\Delta angle (\circ)')

save_figure(fig, [fig_dir, '\' FilePath.locations{ifly,4} '\Change FeTi joint angles opto period'])


%% Compare the avg joint angle for diff legs between prestim and stim
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
fig = getfig('',1); hold on
set(fig, 'color', 'k')
for leg = 1:6
    [S_data, C_data, Stim, Control] = deal([]);
    % ---input---
    iJ = 2;

    disp(['Leg ' num2str(leg)])

    % stim data
    for n = 1:S_tot
        ifly = stimList(n,1); 
        rep = stimList(n,2);
        cond = stimList(n,3);
        Stim(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
    end

    % control:
    for n = 1:C_tot
        ifly = controlList(n,1); 
        rep = controlList(n,2);
        cond = controlList(n,3);
        Control(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
    end

    % ROIs:
    light_on = round((param.basler_delay)*fps);
    light_off = round((param.basler_delay+param.conds_matrix(stimList(1,3)).opto)*fps);
    preROI = 1:light_on;
    lazROI = light_on+1:light_off;

    % find the avg joint angles: 
    S_data(:,1) = mean(Stim(preROI,:)); % stim control period
    S_data(:,2) = mean(Stim(lazROI,:)); % stim laser period

    C_data(:,1) = mean(Control(preROI,:)); % stim control period
    C_data(:,2) = mean(Control(lazROI,:)); % stim laser period

    % plot the avg angle for the two time periods:
    SZ = 70;
    x1 = leg*2-1;
    x2 = leg*2;

    scatter(x1*ones(1,C_tot),C_data(:,1), SZ, 'w', 'filled') % pre
    scatter(x2*ones(1,C_tot),C_data(:,2), SZ, 'w', 'filled') % stim
    y1 = C_data(:,1);
    y2 = C_data(:,2);
    for ii = 1:C_tot
        plot([x1,x2], [y1(ii), y2(ii)], 'color', 'w', 'linewidth', 1)
    end
    %stim:
    x1 = leg*2-1;
    x2 = leg*2;
    c = Color(leg_colors{leg});
    scatter(x1*ones(1,S_tot),S_data(:,1), SZ, c, 'filled') % pre
    scatter(x2*ones(1,S_tot),S_data(:,2), SZ, c,'filled') % stim
    y1 = S_data(:,1);
    y2 = S_data(:,2);
    for ii = 1:C_tot
        plot([x1,x2], [y1(ii), y2(ii)], 'color', c, 'linewidth', 1)
    end
    % run a paired t-test on the before and during data:
    fprintf('\n Joint angle stats: \n')
    [~,p] = ttest(S_data(:,1), S_data(:,2));
    disp(['Stim data p-val: ' num2str(p)])
    [~,p] = ttest(C_data(:,1), C_data(:,2));
    disp(['Control data p-val: ' num2str(p)])
%     [~,p] = ttest(C_data(:,1), S_data(:,1));
%     disp(['Control vs stim preopto p-val: ' num2str(p)])

end
set(gca, 'color', 'k', 'YColor', 'w', 'XColor', 'w', 'TickDir', 'out')
xlim([0,13])
ylim([40,110])
xlabel('control --- stim')
ylabel(['Avg ' Joints{iJ} ' angle (\circ)'])

save_figure(fig, [fig_dir, '\' FilePath.locations{ifly,4} '\FeTi joint angle change all legs']);

%% Plot velocity profiles for a given fly:
ifly = 8;
% Load fictrac data for the given fly
switch questdlg('Load Fictrac data?')
    case 'Yes'
        f_root = [FilePath.locations{ifly,1}, FilePath.locations{ifly,2} '\' FilePath.locations{ifly,3}];
        load([f_root, '\Analysis\' FilePath.locations{ifly,4} '.mat'])
        clc
end

% extract the velocity information for the trials:
[S_data, C_data, S_freq, C_freq] = deal([]);
for n = 1:S_tot
    rep = stimList(n,2);
    cond = stimList(n,3);
    S_data(:,n) = fly.Video.speed(cond).data(:,rep);
    S_freq(:,n) = gait(ifly).Freq(cond,rep).f_avg;
end
for n = 1:C_tot
    rep = controlList(n,2);
    cond = controlList(n,3);
    C_data(:,n) = fly.Video.speed(cond).data(:,rep);
    C_freq(:,n) = gait(ifly).Freq(cond,rep).f_avg;
end

x = (-param.basler_delay:1/fly.param.fps:param.basler_length-param.basler_delay);
x2 = (-param.basler_delay:1/fps:param.basler_length-param.basler_delay);
C = Color('teal');
light_off = param.conds_matrix(stimList(1,3)).opto;



% plot the avg speed for the fly:
fig = getfig('',1);
% FLY SPEED PLOT
subplot(2,1,1)
    hold on
    ylim([0, 18])
    % stim (e.g. 720ms opto group)
    y = smooth(mean(S_data,2),3);
    err = sem(S_data,2);
    plot(x(2:end),y+err, 'linewidth',0.25, 'color', Color('teal'))
    plot(x(2:end),y-err, 'linewidth',0.25, 'color', Color('teal'))
    plot(x(2:end),y, 'linewidth',1, 'color', Color('teal'))

    % control (e.g. 0 ms opto group)
    y = smooth(mean(C_data,2),3);
    err = sem(C_data,2);
    plot(x(2:end),y+err, 'linewidth',0.25, 'color', 'k')
    plot(x(2:end),y-err, 'linewidth',0.25, 'color', 'k')
    plot(x(2:end),y, 'linewidth',1, 'color', 'k')

    % plot opto line:
    ylim([0,1.8])
    xlim([-0.4,1.4])
    y = rangeLine(fig);
    plot([0, light_off], [y,y], 'linewidth', 3, 'color', 'g')

    % labels:
    xlabel('time (s)')
    ylabel('speed (cm/s)')
    set(gca, 'TickDir', 'out')

% STEP FREQUENCY PLOT
subplot(2,1,2)
    hold on
    ylim([3,11])
    % stim (e.g. 720ms opto group)
    y = smooth(mean(S_freq,2),3);
    err = sem(S_freq,2);
    plot(x2,y+err, 'linewidth',0.25, 'color', Color('teal'))
    plot(x2,y-err, 'linewidth',0.25, 'color', Color('teal'))
    plot(x2,y, 'linewidth',1, 'color', Color('teal'))

    % control (e.g. 0 ms opto group)
    y = smooth(mean(C_freq,2),3);
    err = sem(C_freq,2);
    plot(x2,y+err, 'linewidth',0.25, 'color', 'k')
    plot(x2,y-err, 'linewidth',0.25, 'color', 'k')
    plot(x2,y, 'linewidth',1, 'color', 'k')

    % plot opto line:
    % ylim([0,1.8])
    xlim([-0.4,1.4]) %hide the edge effect for step frequency
    y = rangeLine(fig);
    plot([0, light_off], [y,y], 'linewidth', 3, 'color', 'g')
    
    % labels:
    xlabel('time (s)')
    ylabel('step frequency (Hz)')
    set(gca, 'TickDir', 'out')

% Save the figure:
save_figure(fig, [fig_dir, '\' FilePath.locations{ifly,4} '\' set1 ' VS ' set2 ' speed and step freq']);


%% Plot joint angle distributions -- to show which joints might be shifting their angle?
% WORKING HERE % PLOT THE INST. STEP FREQ IN SAME STYLE
leg = 1;

[S_data, C_data, Stim, Control] = deal([]);
% stim data
for n = 1:S_tot
    ifly = stimList(n,1); 
    rep = stimList(n,2);
    cond = stimList(n,3);
    Stim(n,:,:) = angles(ifly).Leg(leg).data(cond,rep).all; %(trial,joint,time)
end

% control:
for n = 1:C_tot
    ifly = controlList(n,1); 
    rep = controlList(n,2);
    cond = controlList(n,3);
    Control(n,:,:) = angles(ifly).Leg(leg).data(cond,rep).all;
end

% ROIs:
light_on = round((param.basler_delay)*fps);
light_off = round((param.basler_delay+param.conds_matrix(stimList(1,3)).opto)*fps);
preROI = 1:light_on;
lazROI = light_on+1:light_off;
% labels
set1 = [num2str(param.conds_matrix(controlList(1,3)).opto) '-ms'];
set2 = [num2str(param.conds_matrix(stimList(1,3)).opto) '-ms'];


% Plot the Pre-stim comparison: 3 rows-joints
Edges = 0:10:180;
fig = getfig('',1);
for iJ = 1:3 %extract joint angles:
    subplot(3,1,iJ); hold on
    %Stim region
    a = Stim(:,lazROI,iJ);
    h = histogram(a,Edges);
    h.FaceColor = Color('orangered'); 
    %pre region
    a = Control(:,lazROI,iJ);
    h = histogram(a,Edges);
    h.FaceColor = Color('white'); 
    % labels: 
    ylabel(Joints{iJ})
    xlim([0,180])
    set(gca, 'XTickLabel', [])
    set(gca, 'TickDir', 'out')
end 
set(gca, 'XTickLabel', (0:20:180))
xlabel('Joint angle (\circ)')
subplot(3,1,1)
title([set1 ' VS ' set2 ' laser Activation'])

save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\',set1,' VS ',set2,' joint angle histo Stim']);



% Plot the Pre-stim comparison: 3 rows-joints
Edges = 0:10:180;
fig = getfig('',1);
for iJ = 1:3 %extract joint angles:
    subplot(3,1,iJ); hold on
    %Stim region
    a = Stim(:,preROI,iJ);
    h = histogram(a,Edges);
    h.FaceColor = Color('orangered'); 
    %pre region
    a = Control(:,preROI,iJ);
    h = histogram(a,Edges);
    h.FaceColor = Color('white'); 
    % labels: 
    ylabel(Joints{iJ})
    xlim([0,180])
    set(gca, 'XTickLabel', [])
    set(gca, 'TickDir', 'out')
end 
set(gca, 'XTickLabel', (0:20:180))
xlabel('Joint angle (\circ)')
subplot(3,1,1)
title([set1 ' VS ' set2 ' control period'])

save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\',set1,' VS ',set2,' joint angle histo PreStim']);


%% Plot all the joint positions in space for a given trial (demo)
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
ifly = 8;
rep = 1; 
cond = 7;

fig = getfig('',1);
set(fig, 'color', 'k');
for LP = 3:5 %number of leg positions tracked

    % pull the tarsus position data: 
    for leg = 1:6
        % find the raw data from the pose 3d
        pos = [leg_labels{leg} Alphabet(LP) '_x'];
        labels = pose_3d(ifly).labels{cond,rep};
        label_loc = find(strcmpi(pos, labels));
        raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

        % plot the leg data:
        plot3(raw(:,1), raw(:,2), raw(:,3), 'linewidth', 1, 'color', Color(leg_colors{leg}))
        hold on
    end
end
grid on
% legend(leg_labels)
axis tight
box off
set(gca,'visible','off')

% optional sphere fit:
Radius = angles(ifly).ball.Radius(cond,rep);
Center = angles(ifly).ball.Center{cond,rep};
[x,y,z] = sphere;
s = surf(Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3));
set(s, 'FaceColor', Color('grey'))
alpha 0.8


axis vis3d % sets the aspect ratio for 3d rotation
ViewZ = [10,60; 50,60; 100,30; 180,30; 280, 10; 360, 10; 410, 10];
VidName = 'C:\Users\evynd\Desktop\All joints video';
OptionZ.FrameRate=10;OptionZ.Duration=10;OptionZ.Periodic=true;
CaptureFigVid(ViewZ, VidName, OptionZ) 


% save figure option:
save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\All joints overlaid'], '-png')

save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\All joints overlaid with sphere'], '-png')



%%  Single leg rotating trace of tarsus

leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
ifly = 3;
rep = 2; 
cond = 1; 
LP = 5; %number of leg positions tracked
leg = 1; %front left leg

fig = getfig('',1);
set(fig, 'color', 'k');

% find the raw data from the pose 3d
pos = [leg_labels{leg} Alphabet(LP) '_x'];
labels = pose_3d(ifly).labels{cond,rep};
label_loc = find(strcmpi(pos, labels));
raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

% plot the leg data:
plot3(raw(:,1), raw(:,2), raw(:,3), 'linewidth', 1, 'color', Color(leg_colors{leg}))
% labels and coloring preferences
set(gca, 'color', 'k', 'ZColor', 'w', 'YColor', 'w', 'XColor', 'w')
grid on
xlabel('x')
ylabel('y')
zlabel('z')

% save figure option:
save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\Single tarsus trace'], '-png');


% grid on
% % legend(leg_labels)
% axis tight
% box off
% set(gca,'visible','off')
% % 
% % optional sphere fit:
% Radius = angles(ifly).ball.Radius(cond,rep);
% Center = angles(ifly).ball.Center{cond,rep};
% [x,y,z] = sphere;
% s = surf(Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3));
% set(s, 'FaceColor', Color('grey'))
% alpha 0.8
% 
% % option to save as a video:
% axis vis3d % sets the aspect ratio for 3d rotation
% ViewZ = [10,60; 50,60; 100,30; 180,30; 280, 10; 360, 10; 410, 10];
% VidName = 'C:\Users\evynd\Desktop\Single joint video';
% OptionZ.FrameRate=10;OptionZ.Duration=10;OptionZ.Periodic=true;
% CaptureFigVid(ViewZ, VidName, OptionZ) 


% add demo points for swing / stance: 
% load STEP data from manually saving it within 'DLC_getStancePointsWALK'

leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
ifly = 3;
rep = 2; 
cond = 1; 
LP = 5; %number of leg positions tracked
leg = 1;

fig = getfig('',1);
set(fig, 'color', 'k');

% find the raw data from the pose 3d
pos = [leg_labels{leg} Alphabet(LP) '_x'];
labels = pose_3d(ifly).labels{cond,rep};
label_loc = find(strcmpi(pos, labels));
raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

% plot the leg data:
plot3(raw(:,1), raw(:,2), raw(:,3), 'linewidth', 1, 'color', Color(leg_colors{leg}))
% labels and coloring preferences
set(gca, 'color', 'k', 'ZColor', 'w', 'YColor', 'w', 'XColor', 'w')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
hold on

% plot the transition points:
% swing-stance driven
loc = find((diff(angles(ifly).stance(cond,rep).loc(leg,:))==0)==0);
scatter3(raw(loc,1), raw(loc,2), raw(loc,3), 100, 'pw', 'filled')


% save figure option:
save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\Single trace transition points'], '-png');


%% Plot the projected cloud points for fitting the sphere:

% add demo points for swing / stance: 
% load STEP data from manually saving it within 'DLC_getStancePointsWALK'
ifly = 3; 
n = 17; % need to match to ID column in angles.ballFit.ID

fig = getfig('',1);
set(fig, 'color', 'k');
% pull the data from the Angles structure:
raw = angles(ifly).ballFit(n).cloud;

% plot the leg data:
scatter3(raw(:,1), raw(:,2), raw(:,3), 100, Color('white'), 'filled')
% labels and coloring preferences
set(gca, 'color', 'k', 'ZColor', 'w', 'YColor', 'w', 'XColor', 'w')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
hold on

% svae the figure:
save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\Fitting Cloud points'], '-png');


% -------- Fit cloud with the projected sphere ----------------: 


% add demo points for swing / stance: 
% load STEP data from manually saving it within 'DLC_getStancePointsWALK'
ifly = 3; 
n = 17; % need to match to ID column in angles.ballFit.ID
cond = angles(ifly).ballFit(n).ID(2);
rep = angles(ifly).ballFit(n).ID(3);

fig = getfig('',1);
set(fig, 'color', 'k');
% pull the data from the Angles structure:
raw = angles(ifly).ballFit(n).cloud;

% plot the leg data:
scatter3(raw(:,1), raw(:,2), raw(:,3), 100, Color('white'), 'filled')
% labels and coloring preferences
set(gca, 'color', 'k', 'ZColor', 'w', 'YColor', 'w', 'XColor', 'w')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
hold on
% 
% axis tight
% box off
% set(gca,'visible','off')

% optional sphere fit:
Radius = angles(ifly).ball.Radius(cond,rep);
Center = angles(ifly).ball.Center{cond,rep};

offset = angles(ifly).ballFit(n).offset;
Center = Center-offset;

[x,y,z] = sphere;
s = surf(Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3));
set(s, 'FaceColor', Color('grey'))
alpha 0.8
axis vis3d % sets the aspect ratio for 3d rotation

save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\Fitting Cloud on sphere'], '-png');


%% Swing vs stance projection on sphere and off sphere:

% add demo points for swing / stance: 
% load STEP data from manually saving it within 'DLC_getStancePointsWALK'

leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
ifly = 8;
rep = 1; 
cond = 7; 
LP = 5; %number of leg positions tracked


fig = getfig('',1);
set(fig, 'color', 'k');

for leg = 1:6
% find the raw data from the pose 3d
pos = [leg_labels{leg} Alphabet(LP) '_x'];
labels = pose_3d(ifly).labels{cond,rep};
label_loc = find(strcmpi(pos, labels));
raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

% plot the leg data:
plot3(raw(:,1), raw(:,2), raw(:,3), 'linewidth', 1, 'color', Color(leg_colors{leg}))
% labels and coloring preferences
set(gca, 'color', 'k', 'ZColor', 'w', 'YColor', 'w', 'XColor', 'w')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
hold on
end

% save figure option:
save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\swing stance projections OFF sphere'], '-png');

% ---- Plot with the stance points hightlighted -------------

fig = getfig('',1);
set(fig, 'color', 'k');

for leg = 1:6
% find the raw data from the pose 3d
pos = [leg_labels{leg} Alphabet(LP) '_x'];
labels = pose_3d(ifly).labels{cond,rep};
label_loc = find(strcmpi(pos, labels));
raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

% plot the leg data:
plot3(raw(:,1), raw(:,2), raw(:,3), 'linewidth', 1, 'color', Color(leg_colors{leg}))
% labels and coloring preferences
set(gca, 'color', 'k', 'ZColor', 'w', 'YColor', 'w', 'XColor', 'w')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
hold on

% % plot the transition points:
% % swing-stance driven
% loc = angles(ifly).stance(cond,rep).loc(leg,:);
% raw(~loc,:) = [];

% scatter3(raw(:,1), raw(:,2), raw(:,3), 30, 'w', 'filled')

end

% save figure option:
save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\swing stance projections with points'], '-png');

% --- plot with sphere projected beneath ---------

fig = getfig('',1);
set(fig, 'color', 'k');

for leg = 1:6
% find the raw data from the pose 3d
pos = [leg_labels{leg} Alphabet(LP) '_x'];
labels = pose_3d(ifly).labels{cond,rep};
label_loc = find(strcmpi(pos, labels));
raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

% plot the leg data:
plot3(raw(:,1), raw(:,2), raw(:,3), 'linewidth', 1, 'color', Color(leg_colors{leg}))
% labels and coloring preferences
set(gca, 'color', 'k', 'ZColor', 'w', 'YColor', 'w', 'XColor', 'w')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
hold on
% 
% % plot the transition points:
% % swing-stance driven
% loc = angles(ifly).stance(cond,rep).loc(leg,:);
% raw(~loc,:) = [];
% 
% scatter3(raw(:,1), raw(:,2), raw(:,3), 30, 'w', 'filled')

end
% Sphere fit:
Radius = angles(ifly).ball.Radius(cond,rep);
Center = angles(ifly).ball.Center{cond,rep};

[x,y,z] = sphere;
s = surf(Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3));
set(s, 'FaceColor', Color('grey'))
alpha 0.8
axis vis3d % sets the aspect ratio for 3d rotation

axis tight
box off
set(gca,'visible','off')

% save figure option:
save_figure(fig, [fig_dir,'\',FilePath.locations{ifly,4},'\swing stance projections no high ON ball'], '-png');


% option to save as a video:
axis vis3d % sets the aspect ratio for 3d rotation
ViewZ = [10,60; 50,60; 100,30; 180,30; 280, 10; 360, 10; 410, 10];
VidName = 'C:\Users\evynd\Desktop\Swing stance projected on sphere';
OptionZ.FrameRate=10;OptionZ.Duration=10;OptionZ.Periodic=true;
CaptureFigVid(ViewZ, VidName, OptionZ) 







