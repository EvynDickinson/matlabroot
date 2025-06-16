
fig_type = '-pdf';
blkbgd = false;


%% Occupancy 'null' distribution for no food trials
clearvars('-except',initial_vars{:})

% link each data set to it's 'null' data set: (make this more involved and
% automated later -- gui driven, maybe?)

null_pair = [1,2]; % first idx is the test trial and the second is the null for each row

% Plot out the quadrant data:

np = 2; %null-pair idx
exp = 1; % active trial
% find the quads with the highest and lowest occupancy over the course of
% the experiment:

sSpan = 180;

dummy = [];
plotData = [];
for i = 1:4
    a = grouped(np).fullquad.(quadOrder{i}).all;
    dummy(i,:) = sum(a,1,'omitnan');
    plotData(:,:,i) = a;
end
% find min and max occupancy quadrants: 
[~, lowerIDX] = min(dummy);
[~, upperIDX] = max(dummy);

[minOcc,maxOcc] = deal([]);
for i = 1:num.trial(np)
    minOcc(:,i) = squeeze(plotData(:,i,lowerIDX(i)));
    maxOcc(:,i) = squeeze(plotData(:,i,upperIDX(i)));
end
    y_err = smooth(mean((minOcc-maxOcc)./2,2,'omitnan'),sSpan,'moving');
    y_avg = smooth(mean([minOcc,maxOcc],2,'omitnan'),sSpan, 'moving');



fig = getfig('',1); 
    hold on
    % plot the null distribution data
    y_err = smooth(mean((minOcc-maxOcc)./2,2,'omitnan'),sSpan,'moving');
    y_avg = smooth(mean([minOcc,maxOcc],2,'omitnan'),sSpan, 'moving');
    kolor = grouped(np).color;
    time = grouped(np).time;
    y1 = smooth(mean(minOcc,2),sSpan, 'moving');
    y2 = smooth(mean(maxOcc,2),sSpan, 'moving');
    plot_error_fills(true, time, y_avg,y_err,kolor,fig_type,0.5);
    % plot(time,y1,'color',kolor)
    % plot(time,y2,'color',kolor)
    % plot the paired food trial on top:
    x = grouped(exp).time;
    kolor = grouped(exp).color;
    y = grouped(exp).fullquad.food.avg;
    y_err = grouped(exp).fullquad.food.std./sqrt(num.trial(exp));
    plot_error_fills(true, x, y, y_err, kolor, fig_type);
    plot(x,y,'color',kolor,'linewidth', 1)
    % formatting
    xlabel('time (min)')
    ylabel('food quadrant occupancy (%)')
    xlim([0 700])
    formatFig(fig,false);
    
save_figure(fig,[figDir, 'full quad occ over time'],fig_type);

%% FIGURE: STATIC caviar trials -- plot heatmap of fly position within arena
clearvars('-except',initial_vars{:})
save_path = createFolder([saveDir 'COM/']);
autoSave = true;

[foreColor] = formattingColors(blkbgd);

% Find the occupancy for each bin:

n = 26; % number of spatial bins
autoLim = false;
axis_limits = [0, 2];
start_time = 60; % when to start counting the behavior
duration_time = 60*4.5; % duration of time to include in the position summary

% Parameters: 
expList = 1:num.exp;
temp_list = [15 17 20 25 27 33 35]; % temps that we have temp hold data for...

% find the time point index: 
plotData = struct;
max_occ = [];
for exp = 1:length(temp_list) % could also do this as auto find of the avg temp for the trial...
    frames = start_time*(fps*60):(start_time+duration_time)*(fps*60);

    % find the center of the arena for this exp type
    Cx = mean(grouped(exp).position.well_pos.x(5,:)); %center X
    Cy = mean(grouped(exp).position.well_pos.y(5,:)); %center Y

    nflies = nan([n,n,num.trial(exp)]); % initialize the variable as zeros
    plotData(exp).wells = nan([2,num.trial(exp)]);
    for trial = 1:num.trial(exp)
        con_type = data(exp).con_type(trial);
        if any(con_type==[1 2]) % check if it matches plate 1 
            x = grouped(exp).position.trial(trial).x(frames,:);
            y = grouped(exp).position.trial(trial).y(frames,:);
            r = conversion(con_type).R*conversion(con_type).pix2mm; % pixel radius of this arena
            % find the edges of the spatial bins 
            x_edge = linspace(Cx-r,Cx+r,n);
            y_edge = linspace(Cy-r,Cy+r,n);
            % find which fly locations go to which spatial bin
            nanLoc = isnan(x)| isnan(y);
            x(nanLoc) = [];
            y(nanLoc) = []; % this also reorganizes the data as a vector
            xInd = discretize(x,x_edge);
            yInd = discretize(y,y_edge);

            % find the number of flies within each spatial bin:
            nflies = zeros(n); % initialize the variable as zeros
            for row = 1:n
                for col = 1:n
                    nflies(row,col,trial) = sum(yInd==row & xInd==col);
                end
            end
            
            xInd = discretize(0,x_edge);
            yInd = discretize(0,y_edge);
            plotData(exp).wells(:,trial) = [xInd,yInd];
        end
    end
    % find the occupancy across all the trials for this experiment group
    tot_flies = sum(nflies,3,'omitnan');
    tot_flies = (tot_flies./sum(sum(tot_flies))*100);

    plotData(exp).data = tot_flies;
    max_occ = max([max_occ,max(max(plotData(exp).data))]);
   disp(exp) 
end


disp(['Max occupancy: ' num2str(max_occ)])

square_unit = mean(diff(x_edge)); % pixel size for one bin
circ_r = r/square_unit; % arena radius in bin size
circ_X = discretize(Cx, x_edge);
circ_Y = discretize(Cy, y_edge);

% PLOT 
r = 1;
c = num.exp;
fig_W = 20 + (400*c);

fig = getfig('',false,[fig_W, 340]); 
    for i = expList
        subplot(r,c,i)
        hold on
        imagesc(plotData(i).data); hold on
        wellX = median(plotData(i).wells(1,:),'omitnan');
        wellY = median(plotData(i).wells(2,:),'omitnan');
        scatter(wellX, wellY,10,'r','filled')
        axis tight;
        axis square;
        % h = drawcircle('Center',[Cx,Cy],'Radius',conversion(1).R*conversion(1).pix2mm,'StripeColor',foreColor);
        v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
    end
    formatFig(fig, blkbgd,[r,c]);
    for i = expList
        subplot(r,c,i)
        set(gca,'XColor','none','Ycolor','none','XTick', [],'YTick', [])
        title(grouped(i).name,'color',foreColor,'fontsize', 12)
        % set(gca,'ColorScale','log')
        C = colorbar;
        C.Label.String = 'Occupancy Probability';
        C.Label.Color = foreColor;
        C.Color = foreColor;
        if autoLim
            clim([0,max_occ])
        else
            clim(axis_limits)
        end
    
        colormap(flipud(gray))
    end
% save the figure to a folder specific to that cohort?
save_figure(fig,[save_path '2D spatial distribution all experiments'], fig_type);

MaxOcc = 0.046287;


%% TODO: Make a structure for the caviar trial of this so they can be compared and then find the difference between them
clearvars('-except',initial_vars{:})
save_path = createFolder([saveDir 'COM/']);
autoSave = true;

[foreColor] = formattingColors(blkbgd);

temp_path = getCloudPath;
temp_path = temp_path(1:end-5);
temp = load([temp_path 'LTS 15-35 temp data.mat']); % pull in the fictive temp region information (and berlin caviar LTS data)
LTS = temp.LTS_temp; clear temp


% Find the occupancy for each bin:
n = 26; % number of spatial bins
autoLim = false;
axis_limits = [0, 2];
start_time = 60; % when to start counting the behavior
duration_time = 60*4.5; % duration of time to include in the position summary

% Parameters: 
expList = 1:num.exp;
temp_list = [15 17 20 25 27 33 35]; % temps that we have temp hold data for...

% find the time point index: 
plotData = struct;
max_occ = [];
for exp = 1:length(temp_list) % could also do this as auto find of the avg temp for the trial...
    frames = start_time*(fps*60):(start_time+duration_time)*(fps*60);

    % find the center of the arena for this exp type
    Cx = mean(grouped(exp).position.well_pos.x(5,:)); %center X
    Cy = mean(grouped(exp).position.well_pos.y(5,:)); %center Y

    nflies = nan([n,n,num.trial(exp)]); % initialize the variable as zeros
    plotData(exp).wells = nan([2,num.trial(exp)]);
    for trial = 1:num.trial(exp)
        con_type = data(exp).con_type(trial);
        if any(con_type==[1 2]) % check if it matches plate 1 
            x = grouped(exp).position.trial(trial).x(frames,:);
            y = grouped(exp).position.trial(trial).y(frames,:);
            r = conversion(con_type).R*conversion(con_type).pix2mm; % pixel radius of this arena
            % find the edges of the spatial bins 
            x_edge = linspace(Cx-r,Cx+r,n);
            y_edge = linspace(Cy-r,Cy+r,n);
            % find which fly locations go to which spatial bin
            nanLoc = isnan(x)| isnan(y);
            x(nanLoc) = [];
            y(nanLoc) = []; % this also reorganizes the data as a vector
            xInd = discretize(x,x_edge);
            yInd = discretize(y,y_edge);

            % find the number of flies within each spatial bin:
            nflies = zeros(n); % initialize the variable as zeros
            for row = 1:n
                for col = 1:n
                    nflies(row,col,trial) = sum(yInd==row & xInd==col);
                end
            end
            
            xInd = discretize(0,x_edge);
            yInd = discretize(0,y_edge);
            plotData(exp).wells(:,trial) = [xInd,yInd];
        end
    end
    % find the occupancy across all the trials for this experiment group
    tot_flies = sum(nflies,3,'omitnan');
    tot_flies = (tot_flies./sum(sum(tot_flies))*100);

    plotData(exp).data = tot_flies;
    max_occ = max([max_occ,max(max(plotData(exp).data))]);
   disp(exp) 
end


disp(['Max occupancy: ' num2str(max_occ)])

square_unit = mean(diff(x_edge)); % pixel size for one bin
circ_r = r/square_unit; % arena radius in bin size
circ_X = discretize(Cx, x_edge);
circ_Y = discretize(Cy, y_edge);

% PLOT 
r = 1;
c = num.exp;
fig_W = 20 + (400*c);

fig = getfig('',false,[fig_W, 340]); 
    for i = expList
        subplot(r,c,i)
        hold on
        imagesc(plotData(i).data); hold on
        wellX = median(plotData(i).wells(1,:),'omitnan');
        wellY = median(plotData(i).wells(2,:),'omitnan');
        scatter(wellX, wellY,10,'r','filled')
        axis tight;
        axis square;
        % h = drawcircle('Center',[Cx,Cy],'Radius',conversion(1).R*conversion(1).pix2mm,'StripeColor',foreColor);
        v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
    end
    formatFig(fig, blkbgd,[r,c]);
    for i = expList
        subplot(r,c,i)
        set(gca,'XColor','none','Ycolor','none','XTick', [],'YTick', [])
        title(grouped(i).name,'color',foreColor,'fontsize', 12)
        % set(gca,'ColorScale','log')
        C = colorbar;
        C.Label.String = 'Occupancy Probability';
        C.Label.Color = foreColor;
        C.Color = foreColor;
        if autoLim
            clim([0,max_occ])
        else
            clim(axis_limits)
        end
    
        colormap(flipud(gray))
    end
% save the figure to a folder specific to that cohort?
save_figure(fig,[save_path '2D spatial distribution all experiments'], fig_type);
% save([save_path '2D spatial distribution.mat'],'plotData', 'circ_r', ...
%     'circ_Y','circ_X', )

MaxOcc = 0.046287;




%% WORK IN PROGRESS -- static hold alignment to heating and cooling paired to fictive temp protocol
% %% FIGURE: STATIC caviar trials -- plot heatmap of fly position within arena
% clearvars('-except',initial_vars{:})
% save_path = createFolder([saveDir 'COM/']);
% autoSave = true;
% 
% [foreColor] = formattingColors(blkbgd);
% 
% % Find the occupancy for each bin:
% 
% 
% 
% n = 26; % number of spatial bins
% autoLim = false;
% axis_limits = [0, 0.01];
% 
% % Parameters: 
% expList = 1:num.exp;
% temp_list = [15 17 20 25 27 33 35]; % temps that we have temp hold data for...
% temp_match_proto = 'Large_temp_sweep_15_35'; % this is where we want to match the time points for the experiment
% 
% % Find the time points that correspond to the desired temp period
% tp = getTempTurnPoints(temp_match_proto);
% temp_path = getCloudPath;
% temp_path = temp_path(1:end-5);
% temp = load([temp_path 'LTS 15-35 temp data.mat']); % pull in the fictive temp region information (and berlin caviar LTS data)
% LTS = temp.LTS_temp; clear temp
% 
% % heating and cooling separated for each desired temp bin
% % find the time index locations for the periods where the behavior should
% % be compiled to match the fictive temp protocol
% 
% % grouped(exp).position.loc(rate, tempbin).data(trial).pos  
% 
% find the time point index: 
% plotData2 = struct;
% max_occ2 = [];
% for tt = 1:length(temp_list) % could also do this as auto find of the avg temp for the trial...
%     temp = temp_list(tt); % this is the temp for the temp hold
%     [~,idx] = min(abs(LTS.temp_list-temp)); %  bin index for this temp
%     nRates = length(LTS.temp_rates);
% 
%     % find the frame numbers for the selected temp
%     for rr = 1:nRates % for heating and cooling, respectively
%         frames = LTS.loc(rr,idx).frames;
%         % find the center of the arena for this exp type
%         Cx = mean(LTS.well_pos.x(5,:)); %center X
%         Cy = mean(LTS.well_pos.y(5,:)); %center Y
% 
%         nflies = nan([n,n,length(temp_list)]); % initialize the variable as zeros
% 
%         for trial = 1:length(temp_list)
%             con_type = data(exp).con_type(trial);
%             if any(con_type==[1 2]) % check if it matches plate 1 
%                 x = grouped(exp).position.trial(trial).x(frames,:);
%                 y = grouped(exp).position.trial(trial).y(frames,:);
%                 r = conversion(con_type).R*conversion(con_type).pix2mm; % pixel radius of this arena
%                 % find the edges of the spatial bins 
%                 x_edge = linspace(Cx-r,Cx+r,n);
%                 y_edge = linspace(Cy-r,Cy+r,n);
%                 % find which fly locations go to which spatial bin
%                 nanLoc = isnan(x)| isnan(y);
%                 x(nanLoc) = [];
%                 y(nanLoc) = []; % this also reorganizes the data as a vector
%                 xInd = discretize(x,x_edge);
%                 yInd = discretize(y,y_edge);
% 
%                 % find the number of flies within each spatial bin:
%                 nflies = zeros(n); % initialize the variable as zeros
%                 for row = 1:n
%                     for col = 1:n
%                         nflies(row,col,trial) = sum(yInd==row & xInd==col);
%                     end
%                 end
% 
%                 xInd = discretize(0,x_edge);
%                 yInd = discretize(0,y_edge);
%                 plotData(exp,rr).wells(:,trial) = [xInd,yInd];
%             end
%         end
%         % find the occupancy across all the trials for this experiment group
%         tot_flies = sum(nflies,3,'omitnan');
%         tot_flies = tot_flies./sum(sum(tot_flies));
% 
%         plotData(exp,rr).data = tot_flies;
%         max_occ = max([max_occ,max(max(plotData(exp,rr).data))]);
%     end
%     disp(exp)
% end

% disp(['Max occupancy: ' num2str(max_occ)])
% 
% % PLOT 
% fig_W = 20 + (400*nRates);
% 
% for i = expList
%     fig = getfig('',false,[fig_W, 340]); 
%     for rr = 1:nRates
%         subplot(1,nRates,rr)
%         hold on
%         imagesc(plotData(i,rr).data); hold on
%         scatter(plotData(i,rr).wells(1,:),plotData(i,rr).wells(2,:),10,'r','filled')
%         axis tight;
%         axis square;
%         % h = drawcircle('Center',[circ_X,circ_Y],'Radius',circ_r,'StripeColor',foreColor);
%         % v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
%     end
%     formatFig(fig, blkbgd,[1,nRates]);
%     for rr = 1:nRates
%         subplot(1,nRates,rr)
%         set(gca,'XColor','none','Ycolor','none','XTick', [],'YTick', [])
%         t_str = [num2str(LTS.temp_rates(rr)) '\circC/min | ' num2str(temp_list(i)) '\circC'];
%         title({grouped(i).name; t_str},'color',foreColor,'fontsize', 12)
% 
%         % set(gca,'ColorScale','log')
% 
%         c = colorbar;
%         c.Label.String = 'Occupancy Probability';
%         c.Label.Color = foreColor;
%         c.Color = foreColor;
%         if autoLim
%             clim([0,max_occ])
%         else
%             clim(axis_limits)
%         end
%     end
%     colormap(flipud(gray))
%     % save the figure to a folder specific to that cohort?
%     save_figure(fig,[save_path grouped(i).name ' 2D spatial distribution'], fig_type,autoSave,true);
% end
% 
% 
% MaxOcc = 0.046287;
% 
% % compare and load the dynamic temp curve occupancy data
% 
% % Why are the hold experiments not matching up in lenght with the LTS?
% 
% fig = figure;
% for i = 1:num.exp
%     subplot(num.exp,1,i); hold on
%     temp = [];
%     time = [];
%     for trial = 1:num.trial
%         temp = autoCat(temp,data(i).data(trial).data.occupancy.temp,false);
%         time = autoCat(time,data(i).data(trial).data.occupancy.time,false);
%     end
%     plot(data(i).data(trial).occupancy.time,data(i).data(trial).data.occupancy.temp)
% end
% 

%% Sleep by location for different hold trials TODO UPDATE TO NEW STRUCTURE FORMAT
clearvars('-except',initial_vars{:})

[foreColor,backColor] = formattingColors(blkbgd);
sz = 50;
r = 1; c = 2;
buff = 0.2;
LW = 2;

% expList = 1:num.exp;
expList = expOrder;

fig = getfig('',1,[1064 837]);
groupNames = []; h = []; p = [];
for idx = 1:length(expList)
    i = expList(idx);
    k = grouped(i).color;
    tp = getTempTurnPoints(grouped(i).fictivetemp);
    trange = [tp.hold(:); tp.down(:); tp.up(:)];
    x_roi = min(trange):max(trange);
    subplot(r,c,1) % quadrant
    hold on
    y = mean(grouped(i).sleep.quad.all(x_roi,:),'omitnan');
    x = idx*ones(length(y),1); 
    boxchart(x, y',"BoxFaceColor",k,"BoxFaceAlpha",1,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
    scatter(x,y,sz, k, 'filled', 'xjitter', 'density')

    % check if the distribution is significantly different from the null distribution
    [h(idx),p(idx)] = ttest(y,25);
    

    subplot(r,c,2) % ring
    hold on
    y = mean(grouped(i).sleep.ring.percent(x_roi,:),'omitnan');

     boxchart(x, y',"BoxFaceColor",k,"BoxFaceAlpha",1,'BoxMedianLineColor',foreColor,'MarkerColor','none',...
        'BoxEdgeColor',foreColor,'WhiskerLineColor',foreColor,'BoxWidth',0.75,'LineWidth',2,'MarkerStyle','none')
     scatter(i*ones(length(y),1),y,sz, k, 'filled', 'xjitter', 'density')
    groupNames{end+1} = grouped(i).name;
end

formatFig(fig, blkbgd,[r,c]);
title_str = {'quadrant', 'ring'};
for ii = 1:2
    subplot(r,c,ii)
    set(gca, 'xtick', 1:num.exp,'xticklabel',groupNames,'XTickLabelRotation',45)
    ylim([0,100])
    set(gca, 'ytick',0:20:100)
    ylabel('flies sleeping in region (%)')
    title(title_str{ii},'color', foreColor)
    h_line(25,'grey', '--',1.5)
end

save_figure(fig,[saveDir 'sleeping flies in quad and ring'],'-png',true,false);
save_figure(fig,[saveDir 'sleeping flies in quad and ring'],'-pdf',true,true);

%% FIGURE: occupancy of inner and outer regions over time and temp: 

%% FIGURE: Time Course for single parameter -- select your metric
% TODO add flies on food to this section
clearvars('-except',initial_vars{:})

plot_err = true;
autoLim = true;
xlim_auto = true; % change the time range for the x axis
time_limits = [0,900]; % time limits if manual control over x-axis range
nMax =  num.exp; 

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(true);

switch questdlg('Plot error?','','True','False', 'Cancel','True')
    case 'True'
        plot_err = true;
    case 'False'
        plot_err = false;
    case 'Cancel'
        return
    case ''
        return
end
if isempty(title_str)
    return
end
fig_dir = [saveDir, dir_end];
% set figure folder
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; % dependent var timecourse
sb(3).idx = 3:c:r*c; %dependent var temp tuning curve

LW = 0.75;
sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true);
for i = num.exp:-1:1
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor)

   % selected parameter time course
    subplot(r,c,sb(2).idx); hold on
        switch dType
            case 1 % single trial lines
                for trial = 1:num.trial(i)
                    if ext
                        y = smooth(grouped(i).(pName).food.all(:,trial),sSpan, 'moving')*scaler;
                    else
                        y = smooth(grouped(i).(pName).all(:,trial),sSpan, 'moving')*scaler;
                    end
                    plot(x,y,'LineWidth',LW,'Color',kolor)
                end
            case {2, 3} % avg line
                if ext 
                    y = smooth(grouped(i).(pName).food.avg,sSpan, 'moving')*scaler;
                else
                    y = smooth(grouped(i).(pName).avg,sSpan, 'moving')*scaler;
                end
                plot(x,y,'LineWidth',LW,'Color',kolor)
        end

    %temp vs dependent variable tuning curve
    subplot(r,c,sb(3).idx); hold on
    
     switch dType
         case 1 % single trial lines
            for trial = 1:num.trial(i)
                if strcmp(pName, 'dist')
                    x = grouped(i).(pName).distavgbytemp(:,1);
                    rawY = [grouped(i).increasing.all(:,trial),grouped(i).decreasing.all(:,trial)];
                else
                    if ext % subregions
                        x = grouped(i).(pName).food.temps;
                        rawY = [grouped(i).(pName).food.increasing.raw(:,trial),grouped(i).(pName).food.decreasing.raw(:,trial)];
                    else
                        x = grouped(i).(pName).temps;
                        rawY = [grouped(i).(pName).increasing.raw(:,trial),grouped(i).(pName).decreasing.raw(:,trial)];
                    end
                end
                y = mean(rawY,2,'omitnan')*scaler;
                plot(x,y,'color',kolor,'linewidth',1.25)
            end

         case 2 % avg lines (combined heating and cooling)
            if strcmp(pName, 'dist')
                x = grouped(i).(pName).distavgbytemp(:,1);
                rawY = [grouped(i).increasing.all,grouped(i).decreasing.all];
            else
                if ext 
                    x = grouped(i).(pName).food.temps;
                    rawY = [grouped(i).(pName).food.increasing.raw,grouped(i).(pName).food.decreasing.raw];
                else
                    x = grouped(i).(pName).temps;
                    rawY = [grouped(i).(pName).increasing.raw,grouped(i).(pName).decreasing.raw];
                end
            end
            y = mean(rawY,2,'omitnan')*scaler;
            y_err = (std(rawY,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
            plot(x,y,'color',kolor,'linewidth',1.25)
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);

         case 3 % separated heating and cooling
            if strcmp(pName, 'dist')
                x = grouped(i).(pName).distavgbytemp(:,1);
                YC = grouped(i).decreasing.all;
                YH = grouped(i).increasing.all;
            else
                if ext 
                    x = grouped(i).(pName).food.temps;
                    YC = grouped(i).(pName).food.decreasing.raw;
                    YH = grouped(i).(pName).food.increasing.raw;
                else
                     x = grouped(i).(pName).temps;
                    YC = grouped(i).(pName).decreasing.raw;
                    YH = grouped(i).(pName).increasing.raw;
                end
            end
            % cooling
            y = mean(YC,2,'omitnan')*scaler;
            y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
            plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '--')
            % heating
            y = mean(YH,2,'omitnan')*scaler;
            y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
            plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
            plot(x,y,'color',kolor,'linewidth',1.25,'linestyle', '-')          
     end
     dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);

% indeendent measure over time
subplot(r,c,sb(2).idx)
ylabel(y_lab)
xlabel('time (min)')
set(gca,'ydir',y_dir)
xlims = xlim;

% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",'none')
xlim(xlims)

% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel(y_lab)
xlabel('temp (\circC)')
% if ~autoLim
%     ylim(dt_lim)
% end
h_line(nullD,'grey',':',2) %36.2
set(gca,'ydir',y_dir)

if ~xlim_auto
    subplot(r,c,sb(1).idx)
    set(gca, 'xlim', time_limits)
    subplot(r,c,sb(2).idx)
    set(gca, 'xlim', time_limits)
    % ylim([0,100])
    % ylim([20,100])
end

% legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[fig_dir 'Timecourse summary ' title_str],fig_type);



























