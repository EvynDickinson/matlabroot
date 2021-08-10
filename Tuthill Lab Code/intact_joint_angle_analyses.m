




% JA_Data
% cell2mat(JA_Data(ii).CoFe)
% group all the data together, each fly is a column, all reps in the column
for kk = 1:13
    a = []; b = [];
    for irep = 1:3%length(JA_Data)
       a = [a; cell2mat(JA_Data(kk).CoFe(:,irep))];
       b = [b; cell2mat(JA_Data(kk).FeTi(:,irep))];  
       JA_Data(kk).CoFe_num(:,irep) = cell2mat(JA_Data(kk).CoFe(:,irep));
       JA_Data(kk).FeTi_num(:,irep) = cell2mat(JA_Data(kk).FeTi(:,irep));
    end
    JA.CoFe(:,kk) = a;
    JA.FeTi(:,kk) = b;
end



 
%% compare joint angles by activity: (use data from loaded intro to step 7)

% behavior = group
% sort the joint angles into groups by activity
behavior_list = {'walking', 'stationary', 'grooming', 'other'};
for kk = 1:num.fly
    for ii = 1:4
        temp = nan(num.conds, num.reps);
        temp(group(kk).(behavior_list{ii})) = JA_Data(kk).CoFe_num(group(kk).(behavior_list{ii}));  
        JA_Data(kk).(['CoFe_' behavior_list{ii}]) = temp;
        temp(group(kk).(behavior_list{ii})) = JA_Data(kk).FeTi_num(group(kk).(behavior_list{ii}));  
        JA_Data(kk).(['FeTi_' behavior_list{ii}]) = temp;
    end
end

% group all flies together
for ii = 1:4
    a = []; b = [];
    for kk = 1:num.fly
       a = [a; JA_Data(kk).(['CoFe_' behavior_list{ii}])];
       b = [b; JA_Data(kk).(['FeTi_' behavior_list{ii}])];  
    end
    JA.(behavior_list{ii}).CoFe = a;
    JA.(behavior_list{ii}).FeTi = b;
end

JA_label.name = {'acute', 'bent', 'extended', 'obtuse'};
JA_label.range = [0,70; 71,100; 101,150; 151,180];

% plot the joint angle distribution by activity
fig = getfig; 
for ii = 1:4
    subplot(2,2,ii)
    histogram(JA.(behavior_list{ii}).CoFe)
    hold on
    histogram(JA.(behavior_list{ii}).FeTi)
    title(['13B x CsChrimson : ' behavior_list{ii}])
    xlabel('Joint Angle (CoFe blue, FeTi red)')
    ylabel('Counts')
    vline([0; JA_label.range(:,2)], 'r--')
    fig = getaxes(fig, 15);
end
save_figure(fig, [figures_dir, structure_name, ' joint angle histograms by state'], '-png');


%% 


load(['C:\matlabroot\behavior class\' structure_name ' Group'])
% 
% % add a new list to the behavior data for the angle bins:
% 
% % Convert behvariors into numbers and then a filter:
% for kk = 1:num.fly
%     for cond = 1:num.conds
%         for rep = 1:num.reps
%             state = group(kk).behavior{cond, rep};
%             phase = group(kk).phase{cond, rep};
%             if ~ischar(state)
%                 state = cell2mat(state);
%             end
%             if ~ischar(phase)
%                 phase = cell2mat(phase);
%             end
%             switch state
%                 case {'stationary', 'Stationary'}
%                     STATE = 1;
%                 case {'walking', 'Walking'}
%                     STATE = 2;
%                 case {'grooming', 'Grooming'}
%                     STATE = 3;
%                 case {'other', 'Other'}
%                     STATE = 4;
%             end
%             switch phase
%                 case {'stance', 'Stance'}
%                     PHASE = 1;
%                 case {'swing', 'Swing'}
%                     PHASE = 2;
%             end 
%             group(kk).STATE(cond, rep) = STATE;
%             group(kk).PHASE(cond, rep) = PHASE;
%             % joint angle analysis:
%             
%         end
%     end
%     group(kk).walking = (group(kk).STATE==2);
%     group(kk).stationary = (group(kk).STATE==1);
%     group(kk).grooming = (group(kk).STATE==3);
%     group(kk).other = (group(kk).STATE==4);
%     group(kk).stance = (group(kk).PHASE==1);
%     group(kk).swing = (group(kk).PHASE==2);
% end
% 
% JA_label.name = {'acute', 'bent', 'extended', 'obtuse'};
% JA_label.range = [0,70; 71,100; 101,150; 151,180];
% clear a
% % % acute
% % 0, 70
% % % bent
% % 71, 100
% % % extended
% % 101, 150
% % % obtuse
% % 151, 180
% for kk = 1:num.fly
%     group(kk).angles = nan(num.conds, num.reps);
%     for ii = 1:4
%         low_lim = JA_label.range(ii,1);
%         high_lim = JA_label.range(ii,2);
%         
%         a = (JA_Data(kk).FeTi_num > low_lim & JA_Data(kk).FeTi_num <= high_lim);
%         group(kk).(JA_label.name{ii}) = a; 
%         group(kk).angles(a) = ii;
%     end
% end
% 
% % all labels named together in group structure:
% for kk = 1:num.fly
% %     group(kk).ANGLE = cell(num.conds, num.reps);
%     for cond = 1:num.conds
%         for rep = 1:num.reps
%            switch group(kk).angles(cond, rep)
%                case 1
%                    group(kk).ANGLE{cond, rep} = JA_label.name{1};
%                case 2
%                   group(kk).ANGLE{cond, rep} = JA_label.name{2};
%                case 3
%                    group(kk).ANGLE{cond, rep} = JA_label.name{3};
%                case 4
%                    group(kk).ANGLE{cond, rep} = JA_label.name{4};
%            end
%         end
%     end
% end
% 
% %reorganize the JA_Data before combining it with the group data
% [JA_Data.Coordinates] = JA_Data.angles; 
% JA_Data = orderfields(JA_Data,[1:0,14,1:13]); 
% JA_Data = rmfield(JA_Data,'angles');
% 
% for kk = 1:num.fly
%     group(kk).JA_Data = JA_Data(kk);
% end
% 
% save(['C:\matlabroot\behavior class\' structure_name ' Group'], 'group')




%% speed trajectories binned by joint angles...
color_list =   {'red', 'Thistle', 'Plum', 'Orchid', 'MediumPurple', 'BlueViolet', 'Indigo',...
                'red', 'Thistle', 'Plum', 'Orchid', 'MediumPurple', 'BlueViolet', 'Indigo',...
                'red', 'Thistle', 'Plum', 'Orchid', 'MediumPurple', 'BlueViolet', 'Indigo',...
                'red', 'Thistle', 'Plum', 'Orchid', 'MediumPurple', 'BlueViolet', 'Indigo'};
% JA_label.name = {'acute', 'bent', 'extended', 'obtuse'};
x = -2:1/30:2;
speed = struct;
fig = getfig;

for iangle = 1:4
    subplot(2,2,iangle)
    hold all
    % state = 'stationary';
    for cond = 1:28 %[2:7, 9:14, 16:21, 23:28]
        ii = 0;
        for kk = 1:num.fly
            idx = 0;
            for rep = 1:num.reps
                % filter for the appropriate behavior and joint angle
                if group(kk).stationary(cond,rep)==true && group(kk).(JA_label.name{iangle})(cond,rep)==true
                    temp = [fly(kk).Control.speed(cond).data(1:end-1,rep);...
                            fly(kk).Stim.speed(cond).data(:,rep)]';
                    plot(x, temp, 'color', Color(color_list{cond}))
                    idx = idx+1;
                    speed(cond).stationary(kk).(JA_label.name{iangle})(:,idx) = temp;
                end
            end
            if idx > 0
                ii = ii +1;
                a = mean(speed(cond).stationary(kk).(JA_label.name{iangle}),2);
                speed(cond).stationary(kk).([JA_label.name{iangle} '_avg']) = a;
                speed(cond).(['stationary_' JA_label.name{iangle}]).all(:,ii) = a;   
            end
        end
        if ii > 0
            b = speed(cond).(['stationary_' JA_label.name{iangle}]).all;
            speed(cond).(['stationary_' JA_label.name{iangle}]).avg = nanmean(b,2);
            speed(cond).(['stationary_' JA_label.name{iangle}]).err = sem(b,2);
        end
    end
    vline(0, 'k-')
    title((JA_label.name{iangle}))
end



fig = getfig;
for iangle = 1:3
    subplot(2,2,iangle)
    % combine cw & ccw:
    % acute data
    for cond = [1:7]
        IDX = 0;
        try
            b = speed(cond).(['stationary_' JA_label.name{iangle}]).all;
            IDX = 2;
        catch
            b = nan(length(x),1);
            IDX = 0;
        end
        try
            c = speed(cond+7).(['stationary_' JA_label.name{iangle}]).all;
            IDX = IDX + 2;
        catch
            c = nan(length(x),1);
            IDX = IDX + 0;
        end
         InputData(cond).x = x';
         InputData(cond).y = [b,c];
         InputData(cond).Ind = IDX;
         InputData(cond).Color = Color(color_list{cond});
    %     InputData(ii).Style = ':'; 
    end
    [fig, InputData] = TimeCourseFigure(InputData, fig, 3);
    ylim([0,1])
    vline(0, 'k-')
    hline(0.3, 'k-')
    title(['Avg speed - stationary flies - ' JA_label.name{iangle} ' angle'])
    fig = getaxes(fig, 15);  
end
save_figure(fig, [figures_dir, structure_name, ' stationary fly velocity profile by joint angle'], '-png');

 
% fill_data = error_fill(x, speed(cond).(['stationary_' JA_label.name{iangle}]).avg,...
%                         speed(cond).(['stationary_' JA_label.name{iangle}]).err);
% h = fill(fill_data.X, fill_data.Y, InputData(ii).Color, 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% plot(x, speed(cond).(['stationary_' JA_label.name{iangle}]).avg, 'LineStyle', LStyle, 'color', InputData(ii).Color,...
%  'LineWidth', 1);



% [fig, InputData] = TimeCourseFigure(InputData, fig, smoothing)
% InputData needs to include:
% x - x data == time data (e.g. -2:2)
% y - y data == speed or rotational velocity over the course of a trial
% Color == color values (rgb) for each type of data in InputData
% smoothing = number for moving average...1 = no smoothing
% EXAMPLE:
% cond = 1; %control for now
% coloropts = {'black', 'mediumblue', 'darkcyan', 'coral'}; %stationary, walking, grooming, other
% for ii = 1:4
%     InputData(ii).x = [];
%     InputData(ii).y = [];
%     InputData(ii).Ind = 1; % 1 means that data won't graph because there
%     isn't any data for that instance
%     InputData(ii).Color = Color(coloropts{ii});
%     InputData(ii).Style = ':'; 
% end
% 
% for kk = 1:num.fly
%   for rep = 1:num.reps
%     % filter:
%     filter = group(kk).STATE(cond,rep);
%     % data:
%     x = (-2:1/30:2)';
%     y = [fly(kk).Control.speed(cond).data(2:end, rep); fly(kk).Stim.speed(cond).data(:, rep)];
%     InputData(filter).x(:,InputData(filter).Ind) = x;
%     InputData(filter).y(:,InputData(filter).Ind) = y;
%     % set the next position
%     InputData(filter).Ind = InputData(filter).Ind+1; 
%   end
% end




