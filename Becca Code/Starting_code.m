load("S:\Evyn\Data structures\Berlin F LRR 25-17 caviar Becca Trials\Berlin F LRR 25-17 caviar Becca Trials post 3.1 data.mat")
clearvars('-except',initial_vars{:})



%% Figure of temp over time

pointsize=25

figure;
x=data(1).data.occupancy.time;
y=data(1).data.occupancy.temp;
% plot(x,y, 'Color', Color('purple'))
scatter(x,y,pointsize,Color('green'),'filled')
xlim([105,110])
ylim([22.4,22.5])

%% Figure of distance over time

figurehandle=figure
hold on

trials=[1:8]
clist=Color('blue','yellow',length(trials))
idx=0;

for trial=trials
    idx=idx+1
location=T.foodLoc(trial)


x=data(trial).data.occupancy.temp;
y=data(trial).data.occupancy.dist2wells(:,location);
scatter(x,y,7,clist(idx,:))
end 

formatFig(figurehandle,false)
xlabel('Temperature\circC')
ylabel('Distance\mm')

% trial=trials(2)
% location=T.foodLoc(trial)
% 
% 
% x=data(trial).data.occupancy.temp;
% y=data(trial).data.occupancy.dist2wells(:,location);
% scatter(x,y,15,Color('teal'),'filled')
