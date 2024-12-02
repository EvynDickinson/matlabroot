

for i = 1:length(dur_loc)
    ii = dur_loc(i);
    m.extroi(i,:) = [ext_start(ii), ext_stop(ii)];
end




% Visualize periods of M wing extension

r = 2;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2;
lw = 1;

% xlimit = [2.8,3];
timebuff = 5; 
timebuff = timebuff/60;


    fig = getfig('',true, [1032 300]); 
    subplot(r,c,sb(1).idx); hold on 
        plot(time, data(M).wingangle(:,1),'color', Color('dodgerblue'),'linewidth', 1) % male wing spread
        plot(time, data(M).wingangle(:,2),'color', Color('cyan'),'linewidth', 1) % male wing spread
        h_line(50,'red','--',1)
        ylabel('M wing angle (\circ)')
        % xlim(xlimit)
    subplot(r,c,sb(2).idx)
        plot(time,T.IFD,'color', foreColor,'LineWidth', lw)
        xlabel('time (m)')
        ylabel('IFD (mm)')
        % xlim(xlimit)
    formatFig(fig, blkbnd, [r,c], sb);


    for i = 1:length(m.extroi)
        xlimit = [time(m.extroi(i,1))-timebuff,time(m.extroi(i,2))+timebuff];
        subplot(r,c,sb(1).idx)
        xlim(xlimit)     
        subplot(r,c,sb(2).idx)
        xlim(xlimit) 
        v_line(time(m.extroi(i,:)),'gold', ':',1)
        save_figure(fig,[figDir 'M wing extension roi ' num2str(xlimit(1)) ' to ' num2str(xlimit(2))],fig_type,1,0);
    end


