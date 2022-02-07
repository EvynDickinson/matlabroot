
sb(1).idx = 1;
sb(2).idx = 2;
subplot(nrows,ncols,sb(2).idx)

expelRange = std(flyCount)*2;
minCount = nflies-expelRange;
maxCount = nflies+expelRange;

goodLoc = (flyCount<maxCount & flyCount>minCount);


% OVERTRACKING OVER TIME
fig = figure;
subplot(2,1,1); hold on
scatter(occupancy.time, flyCount, 5, Color('white'))
plot(occupancy.time, smooth(flyCount,3), 'color',Color('teal'))
h_line(nflies, 'purple','-')
hline([minCount, maxCount], 'r')
xlim([0,occupancy.time(end)+10])
ylabel('Fly count')
xlabel('time (min)')
title('All data')


subplot(2,1,2); hold on
scatter(occupancy.time(goodLoc), flyCount(goodLoc), 5, Color('white'))
% plot(occupancy.time(goodLoc), smooth(flyCount(goodLoc),5), 'color',Color('yellow'))
scatter(occupancy.time(goodLoc), smooth(flyCount(goodLoc),180), 5, Color('yellow'))
xlim([0,occupancy.time(end)+10])
ylabel('Fly count')
xlabel('time (min)')
title('Data points excluding those two STDs from number of flies')

fig = formatFig(fig,true,[2,1]); 
    
save_figure(fig, [analysisDir expName arenaSel ' data cleanup comparison'], '-png');
    
    
%%
    

figure; scatter(1:length(y_loc), y_loc)
x_offset = (b(1)-x_loc);
y_offset = (b(2)-y_loc);

testY = [118.685438308716;127.649601440430;138.208458099365;NaN; NaN; 0;0];
testX = [470.133565979004;464.054027404785;460.798176879883; NaN; NaN; 0; 0];

euDist = sqrt((b(1)-testX).^2 + (b(2)-testY).^2); %euclidian distance from well center
    
figure;
plot(y_loc)
    
    
figure; scatter(1:length(occupancy.IFD), smoothdata(occupancy.IFD,sSpan))