

% processed data previously in 'DLC_behavior_predictions.m' and
% 'DLC_swing_stance_ID.m'
clearvars('-except',initial_vars{:})
var_list = [];
var_list = initial_vars;
% variables that need to be saved during the extraction for loops
new_vars = {'var_list';'cond'; 'ifly'; 'rep'; 'trial';...
            'controlROI'; 'LP'; 'idx'; 'ntotal'; 'trialLoc'};

var_list(end+1:end+length(new_vars)) = new_vars;
% clearvars('-except',var_list{:})

ifly = 3;
% find the trials with swing-stance fitting: 
trialLoc = [];
for rep = 1:num.reps
   a(:,1) = find(angles(ifly).ball.Radius(:,rep)>0);
   a(:,2) = rep*ones(length(a),1);
   trialLoc = [trialLoc; a];
   clear a
end
ntotal = length(trialLoc);

%% Pull Step frequency
LW = 1; % plot line width 
ifly = 3;

% for each trial, pull the step frequency over the full experiement:
% per leg: 
for n = 1:ntotal
    cond = trialLoc(n,1);
    rep = trialLoc(n,2);
    % pull the logical of swing|stance per frame
    raw = angles(ifly).stance(cond,rep).loc'; %flip orientation
    for leg = 1:6 %each leg
        input = raw(:,leg);
        stance_strt = find(diff(input)==1);
        swing_strt = find(diff(input)==-1);
        f = 1./([stance_strt(1); diff(stance_strt)]./fps); % step frequency at each stance start
        
        % extrapolate and fill the gaps to have a fps matrix: 
        Y = nan(param.basler_length*fps+1,1);
        X = [1; stance_strt; param.basler_length*fps+1];
        for ii = 1:length(X)-1
            if ii == length(X)-1; val = f(ii-1); 
            else val = f(ii); end
            strt = X(ii);
            fin = X(ii+1);
            Y(strt:fin-1) = val*ones(fin-strt,1);
        end
        Y(end) = Y(end-1);  
        y = smooth(Y,fps/10,'lowess'); 
        
        % save data into struct:
        freq(leg).stance_strt = stance_strt;
        freq(leg).swing_strt = swing_strt;
        freq(leg).freq = y;
        freq(leg).rawFreq = Y;
        freq(leg).inst_freq = [stance_strt,f];
        F_all(:,leg) = y;
    end
    % assign the frequency data to the Gait structure:
    gait(ifly).Freq(cond,rep).data = freq;
    gait(ifly).Freq(cond,rep).f_all = F_all;
    gait(ifly).Freq(cond,rep).f_avg = mean(F_all,2);
end

% clearvars('-except',var_list{:})


%% Plot step frequency for a single trial:
% select trial: 
n = 2;
cond = trialLoc(n,1);
rep = trialLoc(n,2);
x = -(param.basler_delay):1/fps:param.basler_length-(param.basler_delay);
light_on = (param.conds_matrix(cond).opto);
freq = [];

LW = 1;
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
% blue-front, red-mid, green-hind; dark-left, pastel-right

fig = getfig('',1); hold all
for leg = 1:6
    y = gait(ifly).Freq(cond,rep).data(leg).freq;
    plot(x, y, 'color', Color(leg_colors{leg}), 'linewidth', LW)
    freq(:,leg) = y;
end
% plot the average: 
plot(x, median(freq,2), 'color', 'k', 'linewidth', LW+1)
% laser region lines:
vline([0,light_on],'g')
% labels
xlabel('time (s)')
ylabel('Step frequency (hz)')
fly_ID = strrep([FilePath.locations{ifly,4}, ' R' num2str(rep) 'C' num2str(cond)], '_','-');
title({FilePath.cross; fly_ID})
legend([leg_labels, {'avg'}])

% save option
save_figure(fig, [fig_dir, '\' FilePath.locations{ifly,4} '\raw step frequency'], '-png');


% 
% %% Find Stride length IN PROGRESS 8-18-20
% clearvars('-except',var_list{:})
% 
% 
% LW = 1; % plot line width 
% ifly = 1;
% 
% % for each trial, pull the step frequency over the full experiement:
% for n = 1:ntotal
%     cond = trialLoc(n,1);
%     rep = trialLoc(n,2);
%     % pull the logical of swing|stance per frame
%     raw = angles(ifly).stance(cond,rep).loc'; %flip orientation
%     for leg = 1:6 %each leg
%         input = raw(:,leg);
%         stance_strt = find(diff(input)==1);
%         swing_strt = find(diff(input)==-1);
%         minT = min([length(stance_strt), length(swing_strt)]);
%         if stance_strt(1)<swing_strt(1) 
%             % fly is in swing at start of video
%             strideT = swing_strt(1:minT)-stance_strt(1:minT);
%         else
%             % fly is in stance at start of video (aka swing start is first
%             % transition in the video)
%             strideT = stance_strt(1:minT)-swing_strt(1:minT);
%         end
%         
% %         stride = 
%         
%         f = 1./([stance_strt(1); diff(stance_strt)]./fps); % step frequency at each stance start
%         
%         % extrapolate and fill the gaps to have a fps matrix: 
%         Y = nan(param.basler_length*fps+1,1);
%         X = [1; stance_strt; param.basler_length*fps+1];
%         for ii = 1:length(X)-1
%             if ii == length(X)-1; val = f(ii-1); 
%             else val = f(ii); end
%             strt = X(ii);
%             fin = X(ii+1);
%             Y(strt:fin-1) = val*ones(fin-strt,1);
%         end
%         Y(end) = Y(end-1);  
%         y = smooth(Y,fps/10,'lowess'); 
%         
%         % save data into struct:
%         freq(leg).stance_strt = stance_strt;
%         freq(leg).swing_strt = swing_strt;
%         freq(leg).freq = y;
%         freq(leg).rawFreq = Y;
%         freq(leg).inst_freq = [stance_strt,f];
%     end
%     % assign the frequency data to the Gait structure:
%     gait(ifly).Freq(cond,rep).data = freq;
% end
% 
% % clearvars('-except',var_list{:})
% 
% 
% 
% 
% %% create new data structure with calculated values: 
% % gait
% 



























