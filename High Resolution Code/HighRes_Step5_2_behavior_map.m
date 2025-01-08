

%% Behavior probability map

% 1 = on food
% 2 = sleeping
% 3 = edge-occupancy
% 4 = courtship

behavior = nan([length(time),4]);

behavior(T.FlyOnFood(:,M),1) = 1;
behavior(m.sleep,2) = 2;
behavior(data(M).OutterRing,3) = 3;
behavior(T.CI,4) = 4;

figure;
for i = 1:4
    scatter(time, behavior(:,i))
end

sum(data(M).OutterRing)

% dummy data set: 
ntime = length(time);

behavior = nan([ntime,4]);
p = randperm(ntime);
loc = 1:floor(ntime*0.2):ntime;
for i = 1:4
    rloc = loc(i):loc(i+1);
    behavior(p(rloc),i) = i;
end

% 
figure; hold on
for i = 1:4
    scatter(time, behavior(:,i))
end


% what are the transitions? 

% 
% if sleep co-occurs with other behavior markers, sleep takes priority 
sleep_override = sum(~isnan(behavior(:,1:3)),2)>1 & ~isnan(behavior(:,2));
% TODO: set everything in these locations to just sleeeeeep

for i = 1:4













































