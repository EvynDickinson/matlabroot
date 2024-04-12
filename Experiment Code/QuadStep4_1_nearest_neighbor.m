


%% !SLOW! ANALYSIS: Nearest neighbor analysis
% WARNING: TAKES A MINUTE OR TWO (COMPUTATION HEAVY)
clearvars('-except',initial_vars{:})

% Calculate nearest neighbor distance for each frame

for i = 1:num.exp
    NN_file = [baseFolder,'Data structures\',expNames{i},...
               '\',expNames{i},' nearest neighbor.mat'];
    if exist(NN_file,"file")
        load(NN_file,'NN')
        grouped(i).NN = NN;
    else
        NN = [];
        for trial = 1:num.trial(i)
            x = data(i).data(trial).data.x_loc;
            y = data(i).data(trial).data.y_loc;
            for frame = 1:size(x,1)
                D = pdist([x(frame,:)',y(frame,:)']);
                Z = squareform(D);
                loc = logical(eye(size(Z))); %remove self from list
                Z(loc) = nan;
                minD = min(Z);
                NN(frame) = mean(minD,'omitnan');
                %     % === demo figure ===
                %     [minD, Idx] = min(Z);
                %     fig = figure; hold on
                %     scatter(x(frame,:),y(frame,:),40,'w','filled')
                %     for ii = 1:length(Idx)
                %         plot([x(frame,ii);x(frame,Idx(ii))],[y(frame,ii);y(frame,Idx(ii))])
                %     end
                %     formatFig(fig,true);
                %     xlabel('X'); ylabel('Y')
                %     axis equal
                %     title(['Avg N.N. = ' num2str(mean(minD,'omitnan'))],'color','w')
                %     save_figure(fig,[saveDir expGroup ' nearest neighbor demo frame ' num2str(frame)],'-png');
            end
            % Save data into internal data structure
            grouped(i).NN.data(trial).all = NN./pix2mm; %get mm from pixel distance
        end
        % Save data into external dataset
        NN = grouped(i).NN;
        save(NN_file,'NN'); clear NN
    end
    disp(['Done exp ' num2str(i)])
end
disp('All finished')

for i = 1:num.exp
    NN_all = [];
    for trial = 1:num.trial(i)
        % avg across trials
        NN_all = autoCat(NN_all,grouped(i).NN.data(trial).all);
    end
    grouped(i).NN.avg = mean(NN_all,1,'omitnan');
    grouped(i).NN.err = std(NN_all,0,1,'omitnan');
    grouped(i).NN.all = NN_all;
end

% Calculate the avg nearest neighbor for each temperature
for i = 1:num.exp
    for trial = 1:num.trial(i)
        temps = unique(data(i).G(1).TR.temps);
        rateIdx = data(i).G(trial).TR.rateIdx;
        tempIdx = data(i).G(trial).TR.tempIdx;
        % find rate index
        heatRate = find(data(i).G(trial).TR.rates>0);
        coolRate = find(data(i).G(trial).TR.rates<0);
        holdRate = find(data(i).G(trial).TR.rates==0);
        for temp = 1:length(temps)
            % increasing rates:
            loc = rateIdx==heatRate & tempIdx==temp; %rate and temp align
            grouped(i).NN.increasing(trial,temp) = mean(grouped(i).NN.all(trial,loc),'omitnan');
            % decreasing rates:
            loc = rateIdx==coolRate & tempIdx==temp; %rate and temp align
            grouped(i).NN.decreasing(trial,temp) = mean(grouped(i).NN.all(trial,loc),'omitnan');
            % decreasing rates:
            loc = rateIdx==holdRate & tempIdx==temp; %rate and temp align
            grouped(i).NN.holding(trial,temp) = mean(grouped(i).NN.all(trial,loc),'omitnan');
        end
    end
    grouped(i).NN.temps = temps;
end



%% FIGURE: Nearest neighbor over time and by temp
clearvars('-except',initial_vars{:})
[~,backColor] = formattingColors(blkbgd); 

% ======== TIMECOURSE FIGURE =========
% set up figure aligments
r = 4; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8,10,11]; % nearest neighbor
sb(3).idx = [3,6,9,12]; % nearest neighbor distance by temperature
LW = 1;
sSpan = 180; 

fig = getfig('',true);  hold on
% TEMP
subplot(r,c,sb(1).idx); hold on
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;
    y = grouped(i).temp;
    plot(x,y,'LineWidth',LW,'Color',kolor)
end
ylabel('\circC')
% Nearest neighbor distance
subplot(r,c,sb(2).idx); hold on
for i = 1:num.exp
    x = grouped(i).time;
    y = mean(grouped(i).NN.all,1,'omitnan');
    y = smooth(y,sSpan,'moving');
    y(length(x)+1:end) = [];
    plot(x,y,'color', grouped(i).color,'linewidth',LW)
end
ylabel('Nearest neighbor distance (mm)')
xlabel('Time (min)')
% Nearest neighbor vs. temperature
subplot(r,c,sb(3).idx); hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    for tt = 1:2
        switch tt
            case 1
                y = grouped(i).NN.increasing;
                L_style = '-';
            case 2
                y = grouped(i).NN.decreasing;
                L_style = '--';
        end
        y = mean(y,1,'omitnan');
        x = grouped(i).NN.temps;
        loc = isnan(y);
        x(loc) = []; y(loc) = [];
        plot(x,y,'color',kolor,'linewidth',LW,'linestyle',L_style)
    end
end
ylabel('Nearest neighbor distance (mm)')
xlabel('Temp (\circC)')

% dataString{i} = grouped(i).name;
% legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 5)

formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca,'xcolor', backColor)

% save figure
save_figure(fig,[saveDir expGroup ' nearest neighbor timecourse'],fig_type);


% ======== DISTANCE TO FOOD VS N.N. DISTANCE FIGURE =========
fig = getfig('',true,[497 732]); hold on

for i = 1:num.exp
    kolor = grouped(i).color;
    for tt = 1:2
        switch tt
            case 1
                x = grouped(1).increasing.avg;
                y = grouped(i).NN.increasing;
                L_style = '-';
            case 2
                x = grouped(1).decreasing.avg;
                y = grouped(i).NN.decreasing;
                L_style = '--';
        end
        y = mean(y,1,'omitnan');
        loc = isnan(y) | isnan(x');
        x(loc) = []; y(loc) = [];
        plot(x,y,'color',kolor,'linewidth',LW,'linestyle',L_style)
    end
end
ylabel('Nearest neighbor distance (mm)')
xlabel('Distance to food (mm)')
formatFig(fig, blkbgd);

save_figure(fig,[saveDir expGroup ' nearest neighbor vs food distance'],fig_type);

%% FIGURE: Nearest neighbor-temp correlation

% Correlation between flies on food and distance to food:

clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
LW = 1.5;
corr_coef = [];
buff = 0.2;

fig_size = [80+(80*num.exp), 590];

% get correlation data
for i = 1:num.exp
    pooledData = [];
    % get speed / distance information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp; % temperature
        y = grouped(i).NN.all(trial,:)';       % nearest neighbor distance
        temp = [];
        temp = autoCat(temp,x,false);%temp for each trial within the exp.
        temp = autoCat(temp,y,false);%dist for each trial within the exp
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-distance correlation
        rho = corr(temp);
        corr_coef(i).all(trial) = rho(1,2);
        % save data for pooled comparison
        pooledData = [pooledData; temp];
    end
    % Pooled speed-distance correlation
    rho = corr(pooledData);
    corr_coef(i).group = rho(1,2);
end

% correlation coefficients
fig = getfig('',true,fig_size); hold on
hold on
 for ii = 1:num.exp
   i = expOrder(ii);
   kolor = grouped(i).color;
   xlow = ii-buff-0.1;
   xhigh = ii+buff+0.1;
   x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
   y = corr_coef(i).all;
   y_avg = mean(corr_coef(i).all);
   scatter(x,y,50,kolor,'filled')
   plot([xlow,xhigh],[corr_coef(i).group,corr_coef(i).group],'color',kolor,'linestyle',':','linewidth',LW)
   plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',LW)
 end
 xlim([0.5,num.exp+.5])
 ylabel('correlation between temp-clustering')
 h_line(0,foreColor,':',1)
 formatFig(fig,blkbgd);
 set(gca,'xcolor',backColor)

% save figure
save_figure(fig,[saveDir expGroup ' temp NearestNeighbor correlation'],fig_type);

