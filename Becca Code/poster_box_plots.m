



%% Full experiment aligned
% 
clearvars('-except',initial_vars{:})

blkbgd = false;
foreColor = formattingColors(blkbgd);
buffer = 0.25;
SZ = 50;

x_start = [15 17 20 25 27 33 35]; % manually update this to include the temp for the holds

fig = getfig;
hold on
for ii = 1:num.exp
    i = expOrder(ii);
    % disp(i)
    tp = getTempTurnPoints(data(i).temp_protocol);
    fps = tp.fps;
    roi = [tp.DownROI, tp.UpROI, tp.HoldROI];
    roi_min = length(roi)/(fps*60);
    y = sleep(i).num(roi,:);
    y = y ./ data(i).T.NumFlies.'; % percent of flies sleeping
    y = (mean(y))*100;
    y_avg = median(y,'omitnan');
    x = shuffle_data(linspace(x_start(ii)-buffer,x_start(ii)+buffer,num.trial(i)));

    % save data into structure for loading later
    total_sleep(ii).x = x; total_sleep(ii).y = y; total_sleep(ii).y_avg = y_avg;
    total_sleep(ii).protocol = data(i).temp_protocol; total_sleep(ii).buff = buffer;
    total_sleep(ii).name = data(i).ExpGroup; total_sleep(ii).color = grouped(i).color;

    % plot data
    k = 'black';
    scatter(x,y,SZ,k,'filled')
    % plot([ii-buffer,ii+buffer],[y_avg,y_avg],'color',k,'linewidth',LW)
    boxchart(x_start(ii)*ones([length(y),1]), y',"BoxFaceColor",k,"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    % ylim([0,40])
end

%% Sleep - Fictive temperature aligned

clearvars('-except',initial_vars{:})
% initial_vars{end+1} = 'fakeTemp';

% fakeTemp = load('G:\My Drive\Jeanne Lab\LTS 15-35 temperature data.mat'); % hard coded for SG currently

x_start = [15 17 20 23 25 27 30 33 35]; % manually update this to include the temp for the holds

blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);
buffer = 0.25;
SZ = 75;

y_avg = [];
fig = getfig;
hold on
for ii = 1:num.exp
    i = expOrder(ii);
    % disp(i)
    tp = getTempTurnPoints(data(i).temp_protocol);
    fps = tp.fps;
    % look for time periods with the specific temperature at the hold period
    idx = find(fakeTemp.LTS_temp.temp_list == x_start(ii));
    roi = [fakeTemp.LTS_temp.loc(1, idx).frames; fakeTemp.LTS_temp.loc(2, idx).frames];

    roi_min = length(roi)/(fps*60); % time duration
    y = sleep(i).num(roi,:);
    y = y ./ data(i).T.NumFlies.'; % percent of flies sleeping
    y = (mean(y,'omitnan'))*100;
    y_avg = [y_avg,median(y,'omitnan')];
    x = shuffle_data(linspace(x_start(ii)-buffer,x_start(ii)+buffer,num.trial(i)));

    % plot data
    k = 'black';
    scatter(x,y,SZ,k,'filled')
    % plot([ii-buffer,ii+buffer],[y_avg,y_avg],'color',k,'linewidth',LW)
    boxchart(x_start(ii)*ones([length(y),1]), y',"BoxFaceColor",k,"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    % ylim([0,40])
end
plot(x_start,y_avg,'Color',k,'LineStyle','--','LineWidth',2)


save_figure(fig,[figDir ' sleep box plot'],'-pdf')

%% Ring occupancy - Fictive temperature aligned

clearvars('-except',initial_vars{:})
% initial_vars{end+1} = 'fakeTemp';

% fakeTemp = load('G:\My Drive\Jeanne Lab\LTS 15-35 temperature data.mat'); % hard coded for SG currently

x_start = [15 17 20 23 25 27 30 33 35]; % manually update this to include the temp for the holds

blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);
buffer = 0.25;
SZ = 75;

y_avg = [];
fig = getfig;
hold on
for ii = 1:num.exp
    i = expOrder(ii);
    % disp(i)
    tp = getTempTurnPoints(data(i).temp_protocol);
    fps = tp.fps;
    % look for time periods with the specific temperature at the hold period
    idx = find(fakeTemp.LTS_temp.temp_list == x_start(ii));
    roi = [fakeTemp.LTS_temp.loc(1, idx).frames; fakeTemp.LTS_temp.loc(2, idx).frames];

    roi_min = length(roi)/(fps*60); % time duration
    y = grouped(i).ring.all(roi,:);% percent of flies in ring
    y = mean(y,'omitnan');
    y_avg = [y_avg,median(y,'omitnan')];
    x = shuffle_data(linspace(x_start(ii)-buffer,x_start(ii)+buffer,num.trial(i)));

    % plot data
    k = 'black';
    scatter(x,y,SZ,k,'filled')
    % plot([ii-buffer,ii+buffer],[y_avg,y_avg],'color',k,'linewidth',LW)
    boxchart(x_start(ii)*ones([length(y),1]), y',"BoxFaceColor",k,"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    % ylim([0,40])
end
plot(x_start,y_avg,'Color',k,'LineStyle','--','LineWidth',2)


save_figure(fig,[figDir ' edge occ box plot'],'-pdf')

%% Speed - Fictive temperature aligned

clearvars('-except',initial_vars{:})
% initial_vars{end+1} = 'fakeTemp';

% fakeTemp = load('G:\My Drive\Jeanne Lab\LTS 15-35 temperature data.mat'); % hard coded for SG currently

x_start = [15 17 20 23 25 27 30 33 35]; % manually update this to include the temp for the holds

blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);
buffer = 0.25;
SZ = 75;

y_avg = [];
fig = getfig;
hold on
for ii = 1:num.exp
    i = expOrder(ii);
    % disp(i)
    tp = getTempTurnPoints(data(i).temp_protocol);
    fps = tp.fps;
    % look for time periods with the specific temperature at the hold period
    idx = find(fakeTemp.LTS_temp.temp_list == x_start(ii));
    roi = [fakeTemp.LTS_temp.loc(1, idx).frames; fakeTemp.LTS_temp.loc(2, idx).frames];

    roi_min = length(roi)/(fps*60); % time duration
    y = grouped(i).speed.all(roi,:);% speed
    y = mean(y,'omitnan');
    y_avg = [y_avg,median(y,'omitnan')];
    x = shuffle_data(linspace(x_start(ii)-buffer,x_start(ii)+buffer,num.trial(i)));

    % plot data
    k = 'black';
    scatter(x,y,SZ,k,'filled')
    % plot([ii-buffer,ii+buffer],[y_avg,y_avg],'color',k,'linewidth',LW)
    boxchart(x_start(ii)*ones([length(y),1]), y',"BoxFaceColor",k,"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    % ylim([0,40])
end
plot(x_start,y_avg,'Color',k,'LineStyle','--','LineWidth',2)


save_figure(fig,[figDir ' speed box plot'],'-pdf')

%% Food quadrant occpancy - Fictive temperature aligned

clearvars('-except',initial_vars{:})
% initial_vars{end+1} = 'fakeTemp';

% fakeTemp = load('G:\My Drive\Jeanne Lab\LTS 15-35 temperature data.mat'); % hard coded for SG currently

x_start = [15 17 20 23 25 27 30 33 35]; % manually update this to include the temp for the holds

blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);
buffer = 0.25;
SZ = 75;

y_avg = [];
fig = getfig;
hold on
for ii = 1:num.exp
    i = expOrder(ii);
    % disp(i)
    tp = getTempTurnPoints(data(i).temp_protocol);
    fps = tp.fps;
    % look for time periods with the specific temperature at the hold period
    idx = find(fakeTemp.LTS_temp.temp_list == x_start(ii));
    roi = [fakeTemp.LTS_temp.loc(1, idx).frames; fakeTemp.LTS_temp.loc(2, idx).frames];

    roi_min = length(roi)/(fps*60); % time duration
    y = grouped(i).fullquad.food.all(roi,:);% percent of flies in the food quadrant
    y = mean(y,'omitnan');
    y_avg = [y_avg,median(y,'omitnan')];
    x = shuffle_data(linspace(x_start(ii)-buffer,x_start(ii)+buffer,num.trial(i)));

    % plot data
    k = 'black';
    scatter(x,y,SZ,k,'filled')
    % plot([ii-buffer,ii+buffer],[y_avg,y_avg],'color',k,'linewidth',LW)
    boxchart(x_start(ii)*ones([length(y),1]), y',"BoxFaceColor",k,"BoxFaceAlpha",0.4,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    % ylim([0,40])
end
plot(x_start,y_avg,'Color',k,'LineStyle','--','LineWidth',2)


save_figure(fig,[figDir ' food occ box plot'],'-pdf')