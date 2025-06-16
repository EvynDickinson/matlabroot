


%% Comparisons of food region measures: quadrants vs circles vs distances
clearvars('-except',initial_vars{:})
plot_err = true;
foreColor = formattingColors(blkbgd); %get background colors

r = 2;
c = 3;

lw = 2;
x_lim = [0,700];
sSpan = 180;
FA = 0.4;

regions = {'fullquad','innerquad','circle10', 'circle7', 'circle5'};
region_null = [25,18.75,10,7,5];

fig = getfig('',1);
set(fig, 'windowstyle', 'docked');
for exp = 1:num.exp
    kolor = grouped(exp).color;
    x = grouped(exp).time;
    % distance to food
    subplot(r,c,1); hold on
        y = smooth(grouped(exp).dist.avg,sSpan,'moving');
        y_err = smooth(grouped(exp).dist.err,sSpan,'moving');
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, FA);        
        plot(x,y,'color',kolor,'linewidth', lw)
        ylabel('distance to food (mm)')
        xlim(x_lim)
        set(gca,'ydir', 'reverse');
    % new occupancy measures
    for i = 2:r*c
        subplot(r,c,i); hold on
        y = smooth(grouped(exp).(regions{i-1}).food.avg,sSpan,'moving');
        y_err = smooth(grouped(exp).(regions{i-1}).food.std,sSpan,'moving');
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, FA);        
        plot(x,y,'color',kolor,'linewidth', lw)
        ylabel([regions{i-1} ' %'])
        xlim(x_lim)
        ylim([0,100])
        h_line(region_null(i-1),foreColor,'--',1)
    end
end
formatFig(fig,blkbgd,[r,c]);
for i = 1:3
    subplot(r,c,i); hold on
    set(gca,'xcolor', 'none')
end
for i = 4:6
    subplot(r,c,i); hold on
    set(gca,'xcolor', foreColor)
    xlabel('time (min)')
end

save_figure(fig, [figDir 'Food occupancy comparisons over time'],fig_type);

%% Correlation between distance to food and all the other measures of food attraction: 
clearvars('-except',initial_vars{:})
plot_err = true;
foreColor = formattingColors(blkbgd); %get background colors
sz = 35;
buff = 0.3;
buff2 = 0.4;


% correlation for each trial and plot as a scatter point for each trials
% and each condition type 
regions = {'fullquad','innerquad','circle10', 'circle7', 'circle5'};

for exp = 1:num.exp
    A = grouped(exp).dist.all;
    kolor = grouped(exp).color;
    
    fig = getfig('',1,[744 900]); hold on
        for i = 1:length(regions)
            B = grouped(exp).(regions{i}).food.all;
            y = nan(num.trial(exp),1);
            for trial = 1:num.trial(exp)
                [yA, yB] = rmnan(A(:,trial), B(:,trial), 0, true); % remove nans from data
                R = corrcoef(yA,yB);
                y(trial) = R(1,2);
            end
            y_avg = mean(y,'omitnan');
            x = shuffle_data(linspace(i-buff, i+buff,num.trial(exp)));
            scatter(x,y,sz,kolor,"filled","o")
            plot([i-buff2, i+buff2],[y_avg,y_avg],'color', kolor,'LineWidth',2)       
        end
    formatFig(fig,blkbgd);
    set(gca,'xtick', 1:length(regions),'xticklabel',regions,'ydir', 'reverse')
    ylabel('R (pearsons correlation)')

    save_figure(fig,[figDir, grouped(exp).name, ' ROI corr to distance'],fig_type);
end


%% Temperature-tuning curves for the different parameters: 

params = {'innerquad', 'circle7','ring'};
[~,~,foreColor] = formattingColors(blkbgd,true); %get background colors

r = 1; 
c = 3;
FA = 0.5;
lw = 1.5;

fig = getfig('', 1);
for exp = 1:num.exp
    for i = 1:length(params)
        switch params{i}
            case 'innerquad'
                base = grouped(exp).innerquad.food;
                null_line = 18.75;
            case 'circle7'
                base = grouped(exp).circle7.food;
                null_line = 7;
            case 'ring'
                base = grouped(exp).ring;
                null_line = 25;
        end

        subplot(r,c,i); hold on
        kolor = grouped(exp).color;
            % warming
            y = base.increasing.avg;
            y_err = base.increasing.std./sqrt(num.trial(exp));
            x = base.temps;
            plot_error_fills(true, x,y,y_err,kolor,fig_type,FA);
            plot(x,y,'color', kolor,'linestyle', '-','linewidth', lw);
            % cooling
            y = base.decreasing.avg;
            y_err = base.decreasing.std./sqrt(num.trial(exp));
            x = base.temps;
            plot_error_fills(true, x,y,y_err,kolor,fig_type,FA);
            plot(x,y,'color', kolor,'linestyle', '--','linewidth', lw);
            
            % formatting
            ylabel([params{i} ' occupancy (%)'])
            xlabel('temp (\circC)')
            h_line(null_line, foreColor);
    end
end
formatFig(fig, blkbgd,[r c]);

save_figure(fig,[figDir, 'innerquad circle7 and ring occ temp curves'],fig_type);


%% 
% clearvars('-except',initial_vars{:})
% plot_err = true;
% 
% r = 4;
% c = 1;
% sb(1).idx = 1;
% sb(2).idx = 2:4;
% lw = 2;
% kolor = {'gold', 'grey', 'white', 'grey'};
% x_lim = [0,700];
% 
% fig = getfig('',1);
% set(fig, 'windowstyle', 'docked');
% subplot(r,c,sb(1).idx);
%     x = grouped(exp).time;
%     plot(x,grouped(exp).temp,'color','w', 'linewidth',lw)
%     ylabel('(\circC)')
%     xlim(x_lim)
% subplot(r,c,sb(2).idx)
%     y_all = [];
%     hold on
%     for q = 1:4
%         y = smooth(grouped(exp).quadring.(quadOrder{q}).partial_avg,180,'moving');
%         plot(x,y,'color',Color(kolor{q}),'linewidth', lw)
%         y_all = [y_all, y];
%     end
%     plot(x, sum(y_all,2),'color', 'r')
%     xlim(x_lim)
%     ylabel('quad ring occupancy (%)')
%     xlabel('time (min)')
%     ylim([0, 100])
% formatFig(fig,true,[r,c],sb);
% subplot(r,c,sb(1).idx);
% set(gca,'xcolor', 'none');
% subplot(r,c,sb(2).idx);
% legend(quadOrder, 'textcolor', 'w', 'box', 'off');
% save_figure(fig, [saveDir 'Figures/' grouped(exp).name ' fly quadring occupancy over time'],'-png');

%% Drop origin:
exp = 2;

sSpan = 360;

[r,c] = deal(ceil(sqrt(num.trial(exp))));
y = (grouped(exp).dist.all);

fig = getfig('',1); 
for i = 1:num.trial(exp)
    subplot(r,c,i)
    y1 = y(:,i);
    y1 = smooth(y1,sSpan,'moving');
    y1 = smooth(y1,sSpan,'moving');
    plot(y1,'color', 'w')
    yyaxis right
    plot(grouped(exp).temp,'color', Color('gold'))
end
formatFig(fig, blkbgd,[r c]);
for i = 1:num.trial(exp)
subplot(r, c, i)
    yyaxis left
    ylabel('dist')
    set(gca,'ycolor', 'w','xcolor', 'w')
    yyaxis right
    set(gca,'ycolor', 'none')
end


%% 2D spatial histograms of arena locations....






    
    







































