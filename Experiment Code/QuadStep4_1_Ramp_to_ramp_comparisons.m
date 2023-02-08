
% QuadStep4_1_Ramp_to_ramp_comparisons

%% FIGURE: Ramp to ramp comparisons of food proximity
% Q: do we need a fourth ramp??
clearvars('-except',initial_vars{:})
CList = {'BlueViolet','MediumPurple','Plum','Thistle'};
buff = 0.2;
SZ = 50;

for i = 1:num.exp %experimental group

    % Pull data for this experimental group
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    ROIs = [tp.down(:,1),tp.up(:,2)];
    nRamps = size(ROIs,1);
    [avg_dist,avg_dist_err,dist_range,TD_corr] = deal([]);
    for ramp = 1:nRamps
        idx = ROIs(ramp,1):ROIs(ramp,2);
        temp = grouped(i).dist.all(idx,:);
        % avg distance
        avg_dist = autoCat(avg_dist,mean(temp,2,'omitnan'),false);
        avg_dist_err = autoCat(avg_dist_err,std(temp,0,2,'omitnan'),false);
        % temp range
        dist_range(:,ramp) = range(temp)';
        % temp_correlation
        temperature = grouped(i).temp(idx);
        TD_corr(:,ramp) = corr(temperature,temp)';
    end


    

    fig = figure; set(fig, 'pos',[116 189 1274 622]); 
    % Plot the ramp overlays:
    subplot(1,3,1); hold on
        for ramp = 1:nRamps

            y = avg_dist(:,ramp);
            y_err = avg_dist_err(:,ramp);
            x = 1:length(y); x = x./(3*60); %convert to minutes (fps*sec/min)

            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, Color(CList{ramp}), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
            plot(x,y,'linewidth',1,'color',Color(CList{ramp}))
        end
        set(gca,'YDir','reverse')
        ylabel('Proximity to food (mm)')
        xlabel('Time (min)')

    % Plot the ramp range:
    subplot(1,3,2); hold on
        for ramp = 1:nRamps
            x = shuffle_data(linspace(ramp-buff,ramp+buff,num.trial(i)));
            y = dist_range(:,ramp);
            scatter(x,y,SZ,Color(CList{ramp}),'filled')
            plot([ramp-buff,ramp+buff],[mean(y),mean(y)],'linewidth',1.5,'color', Color(CList{ramp}))
        end
        xlim([1-(buff*3),nRamps+(buff*3)])
        ylabel('Range traveled (mm)')
        set(gca,'xtick',1:nRamps)
        xlabel('Ramp')
    % Plot the temp-distance correlation:
    subplot(1,3,3); hold on
        for ramp = 1:nRamps
            x = shuffle_data(linspace(ramp-buff,ramp+buff,num.trial(i)));
            y = TD_corr(:,ramp);
            scatter(x,y,SZ,Color(CList{ramp}),'filled')
            plot([ramp-buff,ramp+buff],[mean(y),mean(y)],'linewidth',1.5,'color', Color(CList{ramp}))
        end
        xlim([1-(buff*3),nRamps+(buff*3)])
        ylabel('Temp-distance correlation')
        set(gca,'xtick',1:nRamps)
        xlabel('Ramp')

    formatFig(fig,true,[1,3]);

    save_figure(fig,[saveDir expNames{i} ' ramp-by-ramp overlay'],'-png',true);

end

%% FIGURE: is ramp 4 neccessary??
% How does the data with ramps 1-3 compare to with ramps 1-4??
clearvars('-except',initial_vars{:})

for i = 1:num.exp %experimental group

    % Pull data for this experimental group
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    ROIs = [tp.down(:,1),tp.up(:,2)];
    nRamps = size(ROIs,1);
    [avg_dist,avg_dist_err,dist_range,TD_corr] = deal([]);
    bins = tp.threshLow:0.5:tp.threshHigh;
    plotData = [];
    for g = 1:2 % Groups: 1-3 and full 1-4 ramps
        switch g
            case 1
                rampIdx = 1:4; 
                plotData(g).color = 'cyan';
            case 2
                rampIdx = 1:3;
                plotData(g).color = 'gold';
        end              
        for trial = 1:num.trial(i)
            temp = []; ramp_range = [];
            for ramp = rampIdx
                idx = ROIs(ramp,1):ROIs(ramp,2);
                y = data(i).G(trial).TR.data(idx,:);
                temp = autoCat(temp, y, true); % add all ramp data to other ramps
                ramp_range(ramp) = range(y(:,3)); % get the ramp-specific distance range
            end
            % find the avg temp, distance relationship for each trial:
            % [ignores heating and cooling differences]
            proximity = [];
            for bb = 1:length(bins)-1
                loc = temp(:,1)>=bins(bb) & temp(:,1)<bins(bb+1);
                proximity(bb) = mean(temp(loc,3),'omitnan');
            end
            % distance-temp avg
            plotData(g).prox(:,trial) = proximity;
            % correlation
            plotData(g).corr(trial) = corr(temp(:,1),temp(:,3));
            % range
            plotData(g).range(trial) = mean(ramp_range);
        end
    end

    % Overlay figure
    SZ = 50;
    LW = 2;

    fig = figure; 
    set(fig,'pos',[116 189 1274 622]); 
    % Average distance per temp bin
    subplot(1,3,1); hold on
        x = bins(1:end-1);
        for g = 1:2
            kolor = Color(plotData(g).color);
            y = mean(plotData(g).prox,2,'omitnan');
            y_err = std(plotData(g).prox,0,2,'omitnan');
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h, 'facealpha', 0.3)
            plot(x,y,'linewidth',LW,'color',kolor)
        end       
        xlabel('Temperature (\circC)')
        ylabel('Food proximity (mm)')
        set(gca,'ydir', 'reverse')
    % Distance range over ramps
    subplot(1,3,2); hold on
        x = [1,2];
        for trial = 1:num.trial(i)
            y = [plotData(1).range(trial),plotData(2).range(trial)];
            plot(x,y,'linewidth',LW,'Color','w')
        end
        for g = 1:2
            x = g*ones(1,num.trial(i));
            y = plotData(g).range;
            scatter(x,y,SZ,Color(plotData(g).color),'filled')
        end
        xlim([0.8,2.2])
        ylabel('range (mm)')
   % Distance range over ramps
    subplot(1,3,3); hold on
    x = [1,2];
        for trial = 1:num.trial(i)
            y = [plotData(1).corr(trial),plotData(2).corr(trial)];
            plot(x,y,'linewidth',LW,'Color','w')
        end
        for g = 1:2
            x = g*ones(1,num.trial(i));
            y = plotData(g).corr;
            scatter(x,y,SZ,Color(plotData(g).color),'filled')
        end
        ylabel('temp-distance correlation')
        xlim([0.8,2.2])

    formatFig(fig,true, [1,3]);
    subplot(1,3,2)
    set(gca,'xtick',[1,2],'xticklabel',{'full','reduced'})
    subplot(1,3,3)
    set(gca,'xtick',[1,2],'xticklabel',{'full','reduced'})

    save_figure(fig,[saveDir expNames{i} ' temp-prox ramp comparison'],'-png',true);

end

