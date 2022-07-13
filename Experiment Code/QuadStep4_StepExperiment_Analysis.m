% This script looks at data files that had no food and used the step ramp
% temperature protocols

%% Compare movement and speed


% Pull plotting data:
for ii = 1:nfoods+1 %each food type
    % cut off the high and low ends of data (to clean):
    loc = food(ii).N(:,1)>threshHigh | food(ii).N(:,1)<threshLow;
    food(ii).N(loc,:) = [];
    % sort all the data by temperature:
    food(ii).N(:,3) = discretize(food(ii).N(:,1),t_roi);
    for tt = 1:length(t_roi)
        loc = food(ii).N(:,3)==tt;
        food(ii).avg(tt) = mean(food(ii).N(loc,2),1,'omitnan');
        food(ii).err(tt) = std(food(ii).N(loc,2),0,1,'omitnan')/sqrt(ntrials);
    end
end



% group speed over temperature
for trial = 1:ntrials
    % TODO...
end


row = 3;
col = 3;
sb(1).idx = [1,2];
sb(2).idx = [4,5,7,8];
sb(3).idx = [3,6,9];


fig = figure;
% temperature
subplot(row, col, sb(1).idx)
    x = data(1).occupancy.time;
    y = data(1).occupancy.time;
    plot(x,y,'color', 'w', 'LineWidth',1)
    ylabel('\circC')

% speed
subplot(row, col, sb(2).idx)
    hold on
    