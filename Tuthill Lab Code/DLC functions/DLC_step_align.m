

%% How to align and overlay a series of steps given a stimulus start

% Try first with trials that are in stance at the oneset of the ...




%% Pull a specific group of joint traces for a given fly:
clearvars('-except',initial_vars{:})

% compare joint angles for a specific group:
flyList = [1,1,14; 1,2,7; 1,2,14; 1,3,7]; % fly,cond&rep list to include
ntotal = length(flyList);

x = (-param.basler_delay:1/fps:param.basler_length-param.basler_delay);
light_on = 0;
light_Loc = round((param.basler_delay)*fps);
light_off = (param.conds_matrix(cond).opto);
splt = [1:2:6, 2:2:6];

fig = getfig('',1);
% group the data:
for leg = 1:6 % switch through legs
  for iJ = 2 % switch through joints
    for n = 1:ntotal % go through all the trials
        ifly = flyList(n,1); 
        rep = flyList(n,2);
        cond = flyList(n,3);
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


%% align the steps from the first joint peak on a step: -- doesn't look good -- would need more trials
minpeakdist = 15;
minpeakh = 20;
light_Loc = round((param.basler_delay)*fps);

leg = 2; % left front leg
iJ = 2; % FeTi 

for n = 1:ntotal % go through all the trials
    ifly = flyList(n,1); 
    rep = flyList(n,2);
    cond = flyList(n,3);
    data(:,n) = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
    % are they in stance or swing at the stim start?
    SwSt = angles(ifly).stance(cond,rep).loc(leg,:);
    stance(n) = SwSt(light_Loc);
    % find the angle peaks 
    locs = peakseek(data(:,n), minpeakdist, minpeakh);
    locs(locs>light_Loc) = [];
    peakPoint(n) = locs(end); % closest peak to start of stim
end



figure; 
hold all
for n = 1:ntotal
    y = data(:,n);
    y = y-y(1);
    plot(y)
end


%% 
figure; hold all
plot(data(:,stance), 'color', Color('grey'), 'linewidth', 1)
plot(data(:,~stance), 'color', Color('black'), 'linewidth', 1)
vline([light_on, light_off], 'g')






% Single condition: all reps overlaid for each joint
% variables: 
t_delay = param.basler_delay*fps;
leg = 1; %1-LF, 2-LM, 3-LR, 4-RF, 5-RM, 6-RR
Joints = {'CoFe', 'FeTi', 'TiTa'};
CCW_offset = num.conds/4;
LW = 1; %linewidth

% ---- INPUT ------
ifly = 1; 
cond = 7;
% -----------------

light_on = param.conds_matrix(cond).opto;

% Plot the data: 
fig = getfig('',1);
for iJ = 1:num.joints %number of joints
    [CW_input, CCW_input] = deal([]);
    subplot(num.joints,1,iJ)
    hold all
    for rep = 1:num.reps
        % plot all CW trials of the condition
        CW_input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        CCW_input = angles(ifly).Leg(leg).data(cond+CCW_offset,rep).all(:,iJ);
        plot(x, CW_input, 'color', Color('teal'))
        plot(x, CCW_input, 'color', Color('maroon'))
    end
    ylabel({'Joint angle'; angles(1).Leg(1).labels{iJ}})  
    % plot the opto_stim
    yl = max(ylim); %y axis max number 
    plot([0, light_on], [yl,yl], 'linewidth', LW+2, 'color', Color(param.LED_light))
    fig_name = [FilePath.locations{ifly,4} ' ' param.conds_matrix(cond).label];
    fig_name = strrep(fig_name, '_', '-');
    fig_name = strrep(fig_name, 'cw-', '');
    if iJ == 1; title(fig_name); end
end
xlabel('time (s)')


save_figure(fig, [fig_dir, '\rep overlay ' fig_name]);

clearvars('-except',initial_vars{:})