

%% FIGURE: Temp-distance
% correlation ONLY during ramps

clearvars('-except',initial_vars{:})

[foreColor,backColor] = formattingColors(blkbgd);
ylimits = [-0.8,0.1];
autoLim = true;
corr_coef = [];
buff = 0.005;
SZ = 50;
LW = 1.5;

pixWidth = 60; % additional pixel size for image for each extra experiment group
figSize = [pixWidth + (pixWidth*num.exp),590];

% get correlation data
plotData = struct;
for exp = 1:num.exp
    i = expOrder(exp);
    disp(expNames{i})

    pooled_temp = [];
    pooled_dist = [];
    % get speed / distance information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp;
        y = grouped(i).dist.all(:,trial);
        pooled_temp = autoCat(pooled_temp,x,false);
        pooled_dist = autoCat(pooled_dist,y,false);
    end
    % screen out control periods (recovery holds and start/end of exp)
    tp = getTempTurnPoints(data(i).temp_protocol);
    loc = [tp.UpROI, tp.DownROI];
    loc = sort(loc);

    temp = pooled_temp(loc,:);
    dist = pooled_dist(loc,:);
    [rho,pval] = corr(temp,dist);

    plotData(exp).rho = rho(logical(eye(num.trial(i))));
    plotData(exp).pval = pval(logical(eye(num.trial(i))));
    plotData(exp).groupName = expNames(i);
    plotData(exp).color =  grouped(i).color;
    plotData(exp).tp = tp;
end
% correlation coefficients
t_rate_save = [];
fig = getfig('',true,figSize);
hold on
 for i = 1:num.exp
   disp(plotData(i).groupName)
   kolor = plotData(i).color;
   t_rate = abs(plotData(i).tp.rates(1));
   t_rate_save(i) = t_rate;
   disp([grouped(i).name '  '  num2str(t_rate)])
   xlow = t_rate-buff-.05;
   xhigh = t_rate+buff+.05;
   y = plotData(i).rho;
   y_avg = mean(plotData(i).rho);
   x = shuffle_data(linspace(t_rate-buff,t_rate+buff,length(y)));

   scatter(x,y,SZ,kolor,'filled')
   plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',LW)
 end
 % xlim([0.5,num.exp+.5])
 ylabel('temp-distance corr. coef.')
 h_line(0,foreColor,':',1)
 formatFig(fig,blkbgd);
 set(gca,'xtick', t_rate_save,'XTickLabel',[])
 if ~autoLim
     ylim(ylimits)
 end

% save figure
save([saveDir expGroup ' temprate aligned distance correlation ramps only'],'plotData');
save_figure(fig,[saveDir expGroup ' temprate aligned distance correlation ramps only'],'-png',true,false);
save_figure(fig,[saveDir expGroup ' temprate aligned distance correlation ramps only'],'-pdf',true,true);