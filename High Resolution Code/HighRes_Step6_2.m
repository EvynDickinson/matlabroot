
%% Find the temperature alignment across trials...
fig = getfig('',1); hold on
for i = 1:num.trials
    x = fly(i).time;
    y = fly(i).T.temperature;
    plot(x,y)
    z = [fly(i).tRate(:).idx];
    v_line(x(z(2:2:end)),'r')
end
xlabel('time (min)')
ylabel('temp (\circC)')
formatFig(fig);

%% Simple comparison across flies: distance to food over time

% compile the data: 
[pD(1).x, pD(2).x, pD(1).y,pD(2).y] = deal([]);
for i = 1:num.trials
    for sex = 1:2
        x = fly(i).time;
        y = fly(i).T.dist2food(:,sex);
        pD(sex).x = [pD(sex).x, x];
        pD(sex).y = [pD(sex).y, y];
    end
end

% plot the data:
lw = 1;
sSpan = 5*fly(1).fps; %  5 second smoothing

fig = getfig('',1); hold on
for sex = 1:2
    x = mean(pD(sex).x,2);
    y = mean(pD(sex).y,2);
    plot(x,smooth(y,5*60,'moving'),'color', fly(1).data(sex).color,'LineWidth',lw)
end

xlabel('time (min)')
ylabel('distance to food (mm)')
formatFig(fig);






















































































