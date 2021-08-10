
axsize = 14;
  
%% Rotational Velocity ratio
% find the diff between each one and control: 
for cond = 1:14
  for trial = 1:size(rotvelocity_combined(cond).Stim, 2)
    temp(cond).cntl(:,1) = rotvelocity_combined(1).stim_avg;
    temp(cond).controlraw(:,trial) = rotvelocity_combined(1).Control(:,trial);
    offset = temp(cond).cntl(1)-rotvelocity_combined(cond).Stim(1,trial);
    temp(cond).Stim(:,trial) = rotvelocity_combined(cond).Stim(:,trial)+offset;
    temp(cond).stimdiff(:,trial) = abs(temp(cond).cntl-temp(cond).Stim(:,trial));%
    [temp(cond).max(trial,1), temp(cond).max(trial,2)] = nanmax(temp(cond).stimdiff(:,trial));
%     [temp(cond).min(trial,1), temp(cond).min(trial,2)] = nanmin(temp(cond).stimdiff(:,trial));
    % difference in turning:
    temp(cond).cntrl_tot_turn = sum(temp(cond).cntl)/num.fps; %total amount turned during control stim
    temp(cond).stim_tot_turn(trial,1) = sum(rotvelocity_combined(cond).Stim(:,trial))/num.fps;
  end
end

for cond = 1:7 %straight stim
   templot(:,cond) = temp(cond).stim_tot_turn;   
   temp(cond).cntrl_turndiff = mean(templot(:,1));
   turndiff.data(cond) = mean(templot(:,cond))/temp(cond).cntrl_turndiff;
   a = templot(:,cond);
   turndiff.SEM(cond) = std(a)/length(a);
   %error for division:
   turndiff.err(cond) = abs(turndiff.data(cond))*sqrt((turndiff.SEM(cond)/mean(a))^2+(turndiff.SEM(1)/temp(1).cntrl_turndiff)^2);
end
for cond = 8:14 %rot stim
   templot(:,cond) = temp(cond).stim_tot_turn;   
   temp(cond).cntrl_turndiff = mean(templot(:,8));
   turndiff.data(cond) = mean(templot(:,cond))/temp(cond).cntrl_turndiff;
   a = templot(:,cond);
   turndiff.SEM(cond) = std(a)/length(a);
   turndiff.err(cond) = abs(turndiff.data(cond))*sqrt((turndiff.SEM(cond)/mean(a))^2+(turndiff.SEM(8)/temp(8).cntrl_turndiff)^2);
end

for cond = 1:14
    light_lengths(cond) = parameters.conds_matrix(cond).opto;
end


fig = figure; set(fig, 'pos',[15 50 1900 940], 'color','w'); 
    subplot(2,1,1)
    hold all
    x = 2:7;
    xdata = light_lengths(x);
    fill_data = error_fill(xdata, turndiff.data(x), turndiff.err(x));
    h = fill(fill_data.X, fill_data.Y, get_color('darkorange'), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(xdata, turndiff.data(x), 'color', get_color('darkorange'), 'LineWidth', 4)
    hline(1,'k:')
    ylabel('Fraction of control turns')
    xlabel('Activation Length')
    title({'Fraction of turns compared to control'; 'Walking'})
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')
    
    
    subplot(2,1,2)
    hold all
    x = 9:14;
    xdata = light_lengths(x);
    fill_data = error_fill(xdata, turndiff.data(x), turndiff.err(x));
    h = fill(fill_data.X, fill_data.Y, get_color('teal'), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(xdata, turndiff.data(x),  'color', get_color('teal'), 'LineWidth', 4)
    hline(1,'k:')
    ylabel('Fraction of control turns')
    xlabel('Activation Length')
    title({'Fraction of turns compared to control'; 'Turning'})
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')
save_figure(fig, [figures_dir filename(1:end-4) ' turn ratio tuning curve'])


%% Speed ratio:
% find the diff between each one and control: 
for cond = 1:14
  for trial = 1:size(rotvelocity_combined(cond).Stim, 2)
    temp(cond).cntl(:,1) = velocity_combined(1).stim_avg;
    temp(cond).controlraw(:,trial) = velocity_combined(1).Control(:,trial);
    offset = temp(cond).cntl(1)-velocity_combined(cond).Stim(1,trial);
    temp(cond).Stim(:,trial) = velocity_combined(cond).Stim(:,trial);% +offset
    temp(cond).stimdiff(:,trial) = abs(temp(cond).cntl-temp(cond).Stim(:,trial));%
    [temp(cond).max(trial,1), temp(cond).max(trial,2)] = nanmax(temp(cond).stimdiff(:,trial));
%     [temp(cond).min(trial,1), temp(cond).min(trial,2)] = nanmin(temp(cond).stimdiff(:,trial));
    % difference in turning:
    temp(cond).cntrl_tot_turn = nansum(temp(cond).cntl)/num.fps; %total amount turned during control stim
    temp(cond).stim_tot_turn(trial,1) = nansum(velocity_combined(cond).Stim(:,trial))/num.fps;
  end
end

for cond = 1:7 %straight stim
   templot(:,cond) = temp(cond).stim_tot_turn;   
   temp(cond).cntrl_turndiff = nanmean(templot(:,1));
   distdiff.data(cond) = nanmean(templot(:,cond))/temp(cond).cntrl_turndiff;
   a = templot(:,cond);
   distdiff.SEM(cond) = nanstd(a)/length(a);
   %error for division:
   distdiff.err(cond) = abs(distdiff.data(cond))*sqrt((distdiff.SEM(cond)/nanmean(a))^2+(distdiff.SEM(1)/temp(1).cntrl_turndiff)^2);
end
for cond = 8:14 %rot stim
   templot(:,cond) = temp(cond).stim_tot_turn;   
   temp(cond).cntrl_turndiff = nanmean(templot(:,8));
   distdiff.data(cond) = nanmean(templot(:,cond))/temp(cond).cntrl_turndiff;
   a = templot(:,cond);
   distdiff.SEM(cond) = nanstd(a)/length(a);
   distdiff.err(cond) = abs(distdiff.data(cond))*sqrt((distdiff.SEM(cond)/nanmean(a))^2+(distdiff.SEM(8)/temp(8).cntrl_turndiff)^2);
end

for cond = 1:14
    light_lengths(cond) = parameters.conds_matrix(cond).opto;
end


fig = figure; set(fig, 'pos',[15 50 1900 940], 'color','w'); 
    subplot(2,1,1)
    hold all
    x = 2:7;
    xdata = light_lengths(x);
    fill_data = error_fill(xdata, distdiff.data(x), distdiff.err(x));
    h = fill(fill_data.X, fill_data.Y, get_color('darkorange'), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(xdata, distdiff.data(x), 'color', get_color('darkorange'), 'LineWidth', 4)
    hline(1,'k:')
    ylabel('Fraction of control distance')
    xlabel('Activation Length')
    title({'Fraction of displacement compared to control'; 'Walking'})
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')

    subplot(2,1,2)
    hold all
    x = 9:14;
    xdata = light_lengths(x);
    fill_data = error_fill(xdata, distdiff.data(x), distdiff.err(x));
    h = fill(fill_data.X, fill_data.Y, get_color('teal'), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(xdata, distdiff.data(x), 'color', get_color('teal'), 'LineWidth', 4)
    hline(1,'k:')
    ylabel('Fraction of control distance')
    xlabel('Activation Length')
    title({'Fraction of displacement compared to control'; 'Turning'})
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')
    
save_figure(fig, [figures_dir filename(1:end-4) ' displacement ratio tuning curve'])


%% 
% Save the key data:
save_name = [parameters.cross ' tuning curve data'];
save(save_name, 'distdiff', 'turndiff')


%% Speed based Tuning curves for walking and turning:

WDdata = load([parameters.cross ' WD speed breakdown']);
TDdata = load([parameters.cross ' TD speed breakdown']);
   
fig = figure; set(fig, 'pos',[55 55 1600 700], 'color','w');
set(fig, 'Name', 'Walking and Turning Data')
    % SPEED WITH WALKING DATA (top)
    subplot(2,1,1)
    hold all
    color_output = get_color('darkorange', 'yellow', num.bins);
    for bin = 1:num.bins
        %remove NaN from data:
        nan_loc = isnan(WDdata.WD(bin).avg);
        Y = WDdata.WD(bin).avg;
        Y(nan_loc) = [];
        X = WDdata.x;
        X(nan_loc) = [];
        
        % plot the data
        if length(X) == length(WDdata.WD(bin).avg) && length(WDdata.WD(bin).avg) == length(WDdata.WD(bin).err)
            % adjust the start location to the control
            offset = Y(1);
            Y = Y-offset;
            fill_data = error_fill(X, WDdata.WD(bin).avg-offset, WDdata.WD(bin).err);
            h = fill(fill_data.X, fill_data.Y, color_output(bin,:), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        end
        hh = plot(X, Y);
        set(hh,'color', color_output(bin,:), 'LineWidth', 3);
        axis tight

    end
    hline(0, 'k-')
    ylabel('\Delta Speed (cm/s)')
    xlabel('Activation Length')
    title('Change in Speed during Walking', 'Color', 'k')
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')

    % ROT VELOCITY WITH TURNING DATA (bottom graph)
    subplot(2,1,2)
    hold all
    color_output = get_color('teal', 'lightcyan', num.bins);
    for bin = 1:num.bins
        %remove NaN from data:
        nan_loc = isnan(TDdata.TD(bin).avg);
        Y = TDdata.TD(bin).avg;
        Y(nan_loc) = [];
        X = TDdata.x;
        X(nan_loc) = [];
        
        % plot the data
        if length(X) == length(TDdata.TD(bin).avg) && length(TDdata.TD(bin).avg) == length(TDdata.TD(bin).err)
            % offset to 0 change at onset of laser
            offset = Y(1);
            Y = Y-offset;
            fill_data = error_fill(X, TDdata.TD(bin).avg-offset, TDdata.TD(bin).err);
            h = fill(fill_data.X, fill_data.Y, color_output(bin,:), 'EdgeColor','none');
            set(h, 'facealpha', 0.1)
        end
        hh = plot(X, Y);
        set(hh,'color', color_output(bin,:), 'LineWidth', 3);
    axis tight
    end
    hline(0, 'k-')
    ylabel('\Delta Rot Velocity (deg/s)')
    xlabel('Activation Length')
%     legend(h, legend_names, 'TextColor', 'k'); legend('boxoff');
    title('Change in rotational velocity during turning', 'Color', 'k')
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')
    
save_figure(fig, [figures_dir filename(1:end-4) ' speed based tuning curves'])





%% USE THIS FOR MANUAL SELECTION OF 'BEST OF' ROTATIONAL VELOCITY
clear temp x
cond = 11;
temp.Stim = rotvelocity_analysis(cond).stim_data.raw;
temp.Control = rotvelocity_analysis(cond).control_data.raw;
Type = {'Stim', 'Control'};
x.Stim = 0:60;
x.Control = -60:0;
maxii = size(rotvelocity_combined(cond).Stim,2);
% FIG of all the reps and for each fly laid out one by one
cond = 18;
rep = 0;
kk = 1; ind = 1;
fig = figure; set(fig, 'pos', [50, 50, 1400, 890], 'color', 'w')
for ii = 1:maxii
    % rep cycle through reps
    rep = rep+1;
    if rep > num.reps
        rep = 1;
        kk = kk+1; %go to the next fly
        if kk > num.fly
           cond = cond+7;
           ind = 2;
           kk = 1;
        end
    end
    subplot(6,7,ii)
    hold all
    for tt = 1:2
        plot(x.(Type{tt}), temp.(Type{tt}), 'k', 'LineWidth', 3) %control data
        switch ind
            case 1
                adjust = -1;
            case 2
                adjust = 1;
        end
        temp.fly = fly(kk).(Type{tt}).rotvelocity(cond).data(:,rep)*adjust;
        plot(x.(Type{tt}), temp.fly, 'color', 'r', 'LineWidth', 3)
    end
    title({fly(kk).parameters.fly_name;...
           ['cond: ' num2str(cond) ' rep: ' num2str(rep)]})  
    vline(parameters.conds_matrix(cond).opto*num.fps)
end

save_figure(fig, [figures_dir filename(1:end-4) ' rotvelocity examples all trials'])

% % look at individual selected fig:
% kk = 6;
% rep = 2;
% cond = 25;
% figure;
% hold all
% for tt = 1:2
%     plot(x.(Type{tt}), temp.(Type{tt}), 'k', 'LineWidth', 3)
%     temp.fly = fly(kk).(Type{tt}).rotvelocity(cond).data(:,rep);
%     plot(x.(Type{tt}), temp.fly, 'r', 'LineWidth', 3)
% end

%% USE THIS FOR MANUAL SELECTION OF 'BEST OF' SPEED
cond = 4;
temp.Stim = velocity_analysis(cond).stim_data.raw;
temp.Control = velocity_analysis(cond).control_data.raw;

% FIG of all the reps and for each fly laid out one by one
fig = figure; set(fig, 'pos', [50, 50, 1400, 890], 'color', 'w')
cond = 4;
rep = 0;
kk = 1; ind = 1;
for ii = 1:42
    % rep cycle through reps 
    rep = rep+1;
    if rep > num.reps
        rep = 1;
        kk = kk+1; %go to the next fly
        if kk > num.fly
           cond = cond+7;
           ind = 2;
           kk = 1;
        end
    end
    subplot(6,7,ii)
    hold all
    for tt = 1:2
        plot(x.(Type{tt}), temp.(Type{tt}), 'k', 'LineWidth', 3) %control data
        temp.fly = fly(kk).(Type{tt}).speed(cond).data(:,rep)*adjust;
        plot(x.(Type{tt}), temp.fly, 'color', 'r', 'LineWidth', 3)
    end
    title({fly(kk).parameters.fly_name;...
           ['cond: ' num2str(cond) ' rep: ' num2str(rep)]}) 
    vline(0,'k-')
end

save_figure(fig, [figures_dir filename(1:end-4) ' speed examples all trials'])

% % look at individual selected fig:
% kk = 1;
% rep = 1;
% cond = 18;
% figure;
% hold all
% for tt = 1:2
%     plot(x.(Type{tt}), temp.(Type{tt}), 'k', 'LineWidth', 3)
%     temp.fly = fly(kk).(Type{tt}).speed(cond).data(:,rep);
%     plot(x.(Type{tt}), temp.fly, 'r', 'LineWidth', 3)
% end


%% 

% Select fly lines to include
pathway = 'C:\matlabroot\'; % directory
a = dir([pathway, '*tuning curve data.mat']);
for ii = 1:length(a)
    fly_lines{ii} = a(ii).name; %names of the flies for the day
end
indx = listdlg('ListString', fly_lines, 'SelectionMode', 'Multiple');
% Load selected fly lines
for ii = 1:length(indx)
   data(ii) = load([pathway, fly_lines{indx(ii)}]);
   data(ii).name = fly_lines{indx(ii)};
end
% color choices for each fly line (dark = inactive, light = active)
colors_list = {'DarkOrange', 'LimeGreen', 'DarkGreen', 'SkyBlue', 'Blue'};
for ii = 1:length(colors_list)
    color_matrix(ii,:) = get_color(colors_list{ii});  
    data(ii).color = colors_list{ii};
end

% displacement {walking}
fig = figure; set(fig, 'pos',[15 50 1900 940], 'color','w'); 
    subplot(2,1,1) % 
    hold all
    x = 2:7;
    xdata = light_lengths(x);
    for ii = 1:length(fly_lines)
        fill_data = error_fill(xdata, data(ii).distdiff.data(x), data(ii).distdiff.err(x));
        h = fill(fill_data.X, fill_data.Y, color_matrix(ii,:), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(xdata, data(ii).distdiff.data(x), 'color', color_matrix(ii,:), 'LineWidth', 4)
    end
    hline(1,'k:')
    ylabel('Fraction of control distance')
    xlabel('Activation Length')
    title({'Fraction of displacement compared to control'; 'Walking'})
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')

    subplot(2,1,2)
    hold all
    x = 9:14;
    xdata = light_lengths(x);
    for ii = 1:length(fly_lines)
        fill_data = error_fill(xdata, data(ii).turndiff.data(x), data(ii).turndiff.err(x));
        h = fill(fill_data.X, fill_data.Y, color_matrix(ii,:), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(xdata, data(ii).turndiff.data(x),  'color', color_matrix(ii,:), 'LineWidth', 4)
    end
    hline(1,'k:')
    ylabel('Fraction of control turns')
    xlabel('Activation Length')
    title({'Fraction of turns compared to control'; 'Turning'})
    ax = gca;
    set(ax, 'FontWeight', 'Bold', 'FontSize', axsize, 'XColor', 'k', 'YColor', 'k')
    
save_figure(fig, 'Combined Fly Line Tuning Curves')





   
    
    
   













