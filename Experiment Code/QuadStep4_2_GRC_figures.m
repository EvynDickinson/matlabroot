
fig_type = '-pdf';
blkbgd = false;


%% update the color selection for the temp hold trials

% color_list = {'Blue', 'DodgerBlue', 'LightSkyBlue',  'Silver', 'DarkSalmon', 'Red','FireBrick'};
color_list = {'Peachpuff', 'Powderblue','Magenta'};
for i = 1:num.exp
    grouped(i).color = Color(color_list{i});
end


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
axis_limits = [0, 1.5];
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
r = 2;
c = num.exp;
fig_W = 20 + (400*c); % set this to match sizing for trials with heating and cooling....

fig = getfig('',false,[fig_W, 340*2]); 
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

% Save matrix occupancy data so that it can be compared to other trial
% types in the future more easily: 
for i = 1:length(temp_list)
    plotData(i).name = ['static ' num2str(temp_list(i)) 'C'];
end
save([save_path '2D occupancy matrix.mat'],'plotData');

%% FIGURE: DYNAMIC caviar trials -- plot heatmap of fly position within arena
clearvars('-except',initial_vars{:})
save_path = createFolder([saveDir 'COM/']);
autoSave = true;

foreColor = formattingColors(blkbgd);

n = 26; % number of spatial bins
autoLim = false;
axis_limits = [0, 1.5];

% Parameters: 
expList = 1:num.exp;
temp_list = [15 17 20 25 27 33 35]; % temps that we have temp hold data for...\
nTemps = length(temp_list);
plotData = struct;
max_occ = [];

% find the time point index: 
for exp = 1:num.exp % could also do this as auto find of the avg temp for the trial...
    % find the center of the arena for this exp type
    Cx = mean(grouped(exp).position.well_pos.x(5,:)); %center X
    Cy = mean(grouped(exp).position.well_pos.y(5,:)); %center Y

    nflies = nan([n,n,num.trial(exp)]); % initialize the variable as zeros
    plotData(exp).wells = nan([2,num.trial(exp)]);
    for temp = 1:nTemps
        
        % find the temp bin that corresponds to the selected temp:
        [~,idx] = min(abs(grouped(exp).position.temp_list-temp_list(temp))); % temp bin
        nRates = length(grouped(exp).position.temp_rates);
        
        for rr = 1:nRates   
            nflies = nan([n n num.trial(exp)]); 
            for trial = 1:num.trial(exp)
                    con_type = data(exp).con_type(trial);
                    % check if it matches plate 1 
                    if ~any(con_type==[1 2]) 
                        continue
                    end
                    % pull the trial specific data:
                    r = conversion(con_type).R*conversion(con_type).pix2mm; % pixel radius of this arena    
                    x = grouped(exp).position.loc(rr,idx).data(trial).pos(:,1); % x positions of the flies
                    y = grouped(exp).position.loc(rr,idx).data(trial).pos(:,2); % x positions of the flies

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
                    for row = 1:n
                        for col = 1:n
                            nflies(row,col,trial) = sum(yInd==row & xInd==col);
                        end
                    end
                
                xInd = discretize(0,x_edge);
                yInd = discretize(0,y_edge);
                plotData(exp).wells(:,trial) = [xInd,yInd];
            end

            % find the occupancy across all the trials for this experiment group
            tot_flies = sum(nflies,3,'omitnan');
            tot_flies = (tot_flies./sum(sum(tot_flies))*100);
        
            plotData(exp).occ(rr,temp).n = tot_flies;
            max_occ = max([max_occ,max(max(plotData(exp).occ(rr,temp).n))]);
        end
    end
    disp(exp)
end

disp(['Max occupancy: ' num2str(max_occ) '%'])

square_unit = mean(diff(x_edge)); % pixel size for one bin
circ_r = (conversion(1).R*conversion(1).pix2mm)/square_unit; % arena radius in bin size
circ_X = discretize(Cx, x_edge);
circ_Y = discretize(Cy, y_edge);

% PLOT 
for exp = 1:num.exp
    nRates = length(grouped(exp).position.temp_rates);
    r = 2;
    c = nTemps;
    fig_W = 20 + (400*c);
    
    fig = getfig('',false,[fig_W, (340*2)]); 
    i = 1;
    for rr = 1:nRates
        for temp = 1:nTemps
            subplot(r,c,i)
            hold on
            imagesc(plotData(exp).occ(rr,temp).n); 
            hold on
            wellX = median(plotData(exp).wells(1,:),'omitnan');
            wellY = median(plotData(exp).wells(2,:),'omitnan');
            scatter(wellX, wellY,10,'r','filled')
            title([num2str(temp_list(temp)) '\circC @ ' num2str(grouped(exp).position.temp_rates(rr)) '\circC/min'],'color', foreColor)
            axis tight square
            % axis square;
            v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
            i = i+1; % update the subplot location
        end
    end
    % FORMATTING: 
    formatFig(fig, blkbgd,[r,c]);
    for i = 1: nTemps*nRates
        subplot(r,c,i)
        set(gca,'XColor','none','Ycolor','none','XTick', [],'YTick', [])
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
    save_figure(fig,[save_path '2D spatial distribution ' grouped(exp).name], fig_type);

end

% Save matrix occupancy data so that it can be compared to other trial
% types in the future more easily: 
for exp = 1:num.exp
    for rr = 1:nRates
        for temp = 1:nTemps
            plotData(exp).occ(rr,temp).name = [num2str(temp_list(temp)) 'C at ' num2str(grouped(exp).position.temp_rates(rr)) 'C/min'];
        end
    end
end
save([save_path '2D occupancy matrix.mat'],'plotData');

%% FIGURE: HARD CODED STATIC VS DYNAMIC 2D SPATIAL DATA COMPARISON
% load in the static temp hold 2D occupancy data so that we can compare
% across the two conditions more directly....

% currently hard coded but can be automated later TODO 
StaticOcc = load("S:\Evyn\DATA\Grouped Data Structures\Berlin Temp Holds Caviar\COM\2D occupancy matrix.mat");

exp = 1;

% align food wells across the experiment types: 
for i = 1:length(StaticOcc.plotData)
    disp('Static:')
    disp(median(StaticOcc.plotData(i).wells,2,'omitnan'))
    disp('Dynamic:')
    disp(median(plotData(exp).wells,2,'omitnan'))
end
% so only the 25deg hold is shifted by 1 unit in each direction...
raw = StaticOcc.plotData(4).data;
MT = zeros([size(raw,1)+1,size(raw,2)+1]);
MT(2:end,2:end) = raw;
MT(end,:) = [];
MT(:,end) = [];
StaticOcc.plotData(4).data = MT;


% compare for each temp between static and dynamic temps:
diffMat = struct;
exp = 1;
for rr = 1:nRates
    for temp = 1:nTemps
        m = plotData(exp).occ(rr,temp).n-StaticOcc.plotData(temp).data;
        diffMat(rr,temp).n = m;
    end
end

r = 2;
c = nTemps;
fig_W = 20 + (400*c);

fig = getfig('',false,[fig_W, (340*2)]); 
i = 1;
for rr = 1:nRates
    for temp = 1:nTemps
        subplot(r,c,i)
        hold on
        imagesc(diffMat(rr,temp).n); 
        hold on
        wellX = median(plotData(exp).wells(1,:),'omitnan');
        wellY = median(plotData(exp).wells(2,:),'omitnan');
        scatter(wellX, wellY,10,'r','filled')
        title([num2str(temp_list(temp)) '\circC @ ' num2str(grouped(exp).position.temp_rates(rr)) '\circC/min'],'color', foreColor)
        axis tight square
        % axis square;
        v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
        i = i+1; % update the subplot location
    end
end


% Create custom colormap for positive and negative change between static and dynamic trials
n = 256;  % number of colors
z0 = -2.5; z1 = 0; z2 = 1.5; % color axis limits...
% Number of colors below and above zero
n_neg = round(n * abs(z0) / (z2 - z0));
n_pos = n - n_neg;
cmap_neg = [linspace(0,1,n_neg)', linspace(0,1,n_neg)', ones(n_neg,1)];% Blue to white (for negative)
cmap_pos = [ones(n_pos,1), linspace(1,0,n_pos)', linspace(1,0,n_pos)'];% White to red (for positive)
% Combined color map
cmap = [cmap_neg; cmap_pos];

% FORMATTING: 
formatFig(fig, blkbgd,[r,c]);
for i = 1: nTemps*nRates
    subplot(r,c,i)
    set(gca,'XColor','none','Ycolor','none','XTick', [],'YTick', [])
    C = colorbar;
    C.Label.String = '\Delta Occupancy';
    C.Label.Color = foreColor;
    C.Color = foreColor;
    clim([z0, z2])
    colormap(cmap)
end
% save the figure to a folder specific to that cohort?
save_figure(fig,[save_path '2D spatial distribution difference static vs dynamic ' grouped(exp).name], fig_type);





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




%% FIGURE: Scatter plot for specific time points that can be selected from a selected independent variable...
% TODO add flies on food to this section
clearvars('-except',initial_vars{:})

plot_err = true;
autoLim = true;
xlim_auto = true; % change the time range for the x axis
time_limits = [0,900]; % time limits if manual control over x-axis range
nMax =  num.exp; 

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);

% TODO: option for equal spacing or temp-dependent spacing

tempList = [15 20 25 30 35];
% tempList = [17 25];

% set figure folder
fig_dir = createFolder([saveDir, 'Figures/']);

figWidth = 150 + (100*length(tempList));
sz = 35;
buff = 0.5;
tempStep = 2; % gap between diff temp conditions
expStep = 1; % gap between different experiments within the same temp cluster
foreColor = formattingColors(blkbgd);

fig = getfig('',1,[figWidth, 420]);
    [tempLabel, tempMark, stats, group1, group2] = deal([]);

    x_idx = 0;
    hold on 
    for temp = 1:length(tempList)
        for exp = 1:num.exp
            kolor = grouped(exp).color;
            targetTemp = tempList(temp);
             if ext % subregions
                all_temps = grouped(exp).(pName).food.temps;
                [~, tIDX]  = min(abs(targetTemp-all_temps));
                x = all_temps(tIDX); % this is the actual temp that the data corresponds to
                rawY = [grouped(exp).(pName).food.increasing.raw(tIDX,:);grouped(exp).(pName).food.decreasing.raw(tIDX,:)];
            else
                all_temps = grouped(exp).(pName).temps;
                [~, tIDX]  = min(abs(targetTemp-all_temps));
                x = all_temps(tIDX); % this is the actual temp that the data corresponds to
                rawY = [grouped(exp).(pName).increasing.raw(tIDX,:);grouped(exp).(pName).decreasing.raw(tIDX,:)];
            end
            Y = mean(rawY,1,'omitnan');
            Y_avg = mean(Y,'omitnan');
            Y_err = std(Y,0,2,'omitnan');%/sqrt(num.trial(exp));
              
            % plot the data: 
            X = x_idx*ones([1,num.trial(exp)]);
            a = scatter(X,Y,sz,kolor,'filled', 'XJitter','density');
            plot([x_idx-buff,x_idx+buff],[Y_avg, Y_avg], 'color', foreColor, 'linewidth',2)
            errorbar(x_idx, Y_avg, Y_err, 'color', foreColor, 'LineWidth',1,'CapSize',2)

            % update plotting location tracking
            tempLabel = [tempLabel, x]; % keep track of the ACTUAL temp value at this location
            tempMark = [tempMark, x_idx]; % the plot location for the data
            x_idx = x_idx+expStep;

            stats = [stats; Y'];
            group1 = [group1; x*ones(size(Y'))];
            group2 = [group2; exp*ones(size(Y'))];

        end
        if temp < length(tempList)
            x_idx = x_idx+tempStep;
        end
    end

% formatting
formatFig(fig,blkbgd);
xlim([-tempStep,x_idx+1])
set(gca, 'XTick',tempMark,'xticklabel',cellstr(string(tempLabel)))
xlabel('temp (\circ)')
ylabel(y_lab)

save_figure(fig,[fig_dir, pName ' occ scatter'],fig_type);

% disp('manually update and run the stats')
% return
% 
% % Multiway ANOVA
% [p, tbl, stats_out] = anovan(stats, {group1, group2}, 'model', 'interaction', 'varnames', {'Temp','ExpGroup'});
% multcompare(stats_out, 'Dimension', 1)  % for temp
% multcompare(stats_out, 'Dimension', 2)  % for exp group
% 
% % WORKING HERE: TODO 6.17 STATS on the individual comparisons....
% % manually check which stats are important to you...
% groupA1B1 = stats(group1 == 15 & group2 == 1);
% groupA2B2 = stats(group1 == 15 & group2 == 2);
% [h,p] = ttest2(groupA1B1, groupA2B2);
% 
% [results,~,~,gnames] = multcompare(stats_out,"Dimension",[1 2]);

%% FIGURE: Scatter plot for specific time points separated by heating and cooling for each exp 
% TODO add flies on food to this section
clearvars('-except',initial_vars{:})

plot_err = true;
autoLim = true;
xlim_auto = true; % change the time range for the x axis
time_limits = [0,900]; % time limits if manual control over x-axis range
nMax =  num.exp; 

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);

tempList = [15 25 35];
% tempList = [17 25];

% set figure folder
fig_dir = createFolder([saveDir, 'Figures/']);

figWidth = 200 + (100*length(tempList));
sz = 50;
buff = 0.5;
tempStep = 2; % gap between diff temp conditions
expStep = 1; % gap between different experiments within the same temp cluster
foreColor = formattingColors(blkbgd);
r = 1;
c = 2; % heating and cooling
buff = 0.5;
edgeBuff = 3;
ylims = [0, 100];

timeROI = 1:data(1).fps * 3600 * 4.5 ; % time region to average over...4.5 hours

for exp = 1:num.exp
    if ext
        yBase = grouped(exp).(pName).food;
    else
        yBase = grouped(exp).(pName);
    end
        
    fig = getfig('',1,[figWidth, 420]);
    kolor = grouped(exp).color;
    for tt = 1:2 % cooling then decreasing
        switch tt
            case 1
                tType = 'decreasing';
                % kolor = Color('dodgerblue');
                tName = 'cooling';
            case 2
                tType = 'increasing';
                % kolor = Color('red');
                tName = 'warming';
        end

        % cooling data
        subplot(r,c,tt); hold on
        for temp = 1:length(tempList)
            targetTemp = tempList(temp);
            % pull data and seperate by dynamic or static 
            if data(exp).hold_exp % hold experiments
                y_all = mean(yBase.all(timeROI,:),1,'omitnan');
            else % dyunamic experiments
                 if ext % subregions
                    all_temps = grouped(exp).(pName).food.temps;
                    [~, tIDX]  = min(abs(targetTemp-all_temps));
                    x = all_temps(tIDX); % this is the actual temp that the data corresponds to
                    y_all = grouped(exp).(pName).food.(tType).raw(tIDX,:);
                else
                    all_temps = grouped(exp).(pName).temps;
                    [~, tIDX]  = min(abs(targetTemp-all_temps));
                    x = all_temps(tIDX); % this is the actual temp that the data corresponds to
                    y_all = grouped(exp).(pName).(tType).raw(tIDX,:);
                 end

                 % check that the target temp and selected temps match
                 if ~(targetTemp==x)
                     warndlg('target temp not aligned with selected temp')
                     return
                 end
            end

            Y = mean(y_all,1,'omitnan');
            Y_avg = mean(Y,'omitnan');
            Y_err = std(Y,0,2,'omitnan');%/sqrt(num.trial(exp));

            % plot the data: 
            X = shuffle_data(linspace(targetTemp-buff,targetTemp+buff,num.trial(exp)));
            a = scatter(X,Y,sz,kolor,'filled');
            plot([targetTemp-(3*buff),targetTemp+(3*buff)],[Y_avg, Y_avg], 'color', foreColor, 'linewidth',2)
            errorbar(targetTemp, Y_avg, Y_err, 'color', foreColor, 'LineWidth',1,'CapSize',2)
            
        end
        % Formatting
        title(tName)
        ylabel(pName)
        xlabel('temp (\circC) ')
        set(gca, 'xtick', tempList)
        xlim([min(tempList)-edgeBuff, max(tempList)+edgeBuff])
        ylim(ylims)
    end
    % formatting full set 
    formatFig(fig, blkbgd, [r c]);
    save_figure(fig,[fig_dir, grouped(exp).name ' ' pName ' limited temp occ scatter'],fig_type);
end

       
    

%% LTS food vs no food: temp tuning curves for plotting separeated by heating and cooling
clearvars('-except',initial_vars{:})
autoLims = true;
[title_str, param,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);

% link each data set to it's 'null' data set: (make this more involved and
% automated later -- gui driven, maybe?)

% Parameters: 
null_pair = [1,2]; % first idx is the test trial and the second is the null for each row
np = 2; %null-pair idx
exp = 1; % active trial

% find the quads with the highest and lowest occupancy over the course of the experiment:
if ext % there are quadrant specific features
    dummy = [];
    plotData = [];
    for i = 1:4
        a = grouped(np).innerquad.(quadOrder{i}).all;
        dummy(i,:) = sum(a,1,'omitnan');
        % plotData(:,:,i) = a;
    end
    % find min and max occupancy quadrants: 
    [~, lowerIDX] = min(dummy);
    [~, upperIDX] = max(dummy);

    % extract trial specific information for averaging
    [cool_low, cool_high, warm_low, warm_high] = deal([]);
    for trial = 1:num.trial(np)
        iLow = lowerIDX(trial); % quadrant with lowest occ throughout experiment
        iHigh = upperIDX(trial);  % quadrant with highest occ throughout the experiment 
        cool_low = [cool_low, grouped(np).(param).(quadOrder{iLow}).decreasing.raw(:,trial)];
        cool_high = [cool_high, grouped(np).(param).(quadOrder{iHigh}).decreasing.raw(:,trial)];
        warm_low = [warm_low, grouped(np).(param).(quadOrder{iLow}).increasing.raw(:,trial)];
        warm_high = [warm_high, grouped(np).(param).(quadOrder{iHigh}).increasing.raw(:,trial)];
    end
    cool_low = mean(cool_low,2,'omitnan');
    cool_high = mean(cool_high,2,'omitnan');
    warm_low = mean(warm_low,2,'omitnan');
    warm_high = mean(warm_high,2,'omitnan');
    % find the avg and 'err' from the low and high occupancy regions
    yC_err = mean((cool_high-cool_low)./2,2,'omitnan');
    yC_avg = mean([cool_low,cool_high],2,'omitnan');
    yW_err = mean((warm_high-warm_low)./2,2,'omitnan');
    yW_avg = mean([warm_low,warm_high],2,'omitnan');
end

r = 1;
c = 2;
if ext
    x = grouped(np).(param).food.temps; % temperature bins from low to high in this experiment type
else
    x = grouped(np).(param).temps;
end

LW = 2;
FA = 0.5;
xlims = [14, 36];
if autoLims
    ylims = [];
else
    ylims = [0 90];
end

fig = getfig('',1);
% cooling temperatures 
subplot(r,c,1); hold on
    % 'null' no food distribution 
    if ext
        plot_error_fills(true, x, yC_avg,yC_err,grouped(np).color,fig_type,FA);
        plot(x, yC_avg, 'color', grouped(np).color,'linewidth', LW)
    else
        y = grouped(np).(param).decreasing.avg;
        try y_err = grouped(np).(param).decreasing.std./sqrt(num.trial(np));
        catch
            y_err = grouped(np).(param).decreasing.err./sqrt(num.trial(np));
        end
        plot_error_fills(true, x, y, y_err, grouped(np).color,fig_type,FA);
        plot(x, y, 'color', grouped(np).color,'linewidth', LW)
    end
    % food distribution 
    if ext
        baseY = grouped(exp).(param).food;
    else
        baseY = grouped(exp).(param);
    end
    y = baseY.decreasing.avg;
    try y_err = baseY.decreasing.std./sqrt(num.trial(exp));
    catch
        y_err = baseY.decreasing.err./sqrt(num.trial(exp));
    end
    plot_error_fills(true, x, y, y_err, grouped(exp).color,fig_type,FA);
    plot(x, y, 'color', grouped(exp).color,'linewidth', LW)
    %formatting
    set(gca,'xdir', 'reverse')
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('cooling')
    if autoLims
        ylims = [ylims; ylim];
    end

% warming temperatures 
subplot(r,c,2); hold on
    % 'null' no food distribution 
    if ext
        plot_error_fills(true, x, yW_avg,yW_err,grouped(np).color,fig_type,FA);
        plot(x, yW_avg, 'color', grouped(np).color,'linewidth', LW)
    else
        y = grouped(np).(param).increasing.avg;
        try y_err = grouped(np).(param).increasing.std./sqrt(num.trial(np));
        catch 
            y_err = grouped(np).(param).increasing.err./sqrt(num.trial(np));
        end
        plot_error_fills(true, x, y, y_err, grouped(np).color,fig_type,FA);
        plot(x, y, 'color', grouped(np).color,'linewidth', LW)
    end
    % food distribution 
    if ext
        baseY = grouped(exp).(param).food;
    else
        baseY = grouped(exp).(param);
    end
    y = baseY.increasing.avg;
    try y_err = baseY.increasing.std./sqrt(num.trial(exp));
    catch
        y_err = baseY.increasing.err./sqrt(num.trial(exp));
    end
    plot_error_fills(true, x, y, y_err, grouped(exp).color,fig_type,FA);
    plot(x, y, 'color', grouped(exp).color,'linewidth', LW)
    % formatting
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('warming')       
    if autoLims
        ylims = [ylims; ylim];
    end
formatFig(fig, blkbgd,[r c]);
for i = 1:2
    subplot(r,c,i); hold on
    if autoLims
        ylim([min(ylims(:,1)),max(ylims(:,2))]);
    else 
        ylim(ylims)
    end
end

save_figure(fig,[figDir, param ' occ temp tuning curve separated heating and cooling'],fig_type);

%% LTS STATIC TEMP HOLDS occupancy temp tuning curves
clearvars('-except',initial_vars{:})
[title_str, param,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);
foreColor = formattingColors(blkbgd);
r = 1;
c = 2;
autoLims = true;

timeROI = 1:data(1).fps * 3600 * 5 ; % time region to average over...5 hours

plotData = [];
for exp = 1:num.exp
    if ext
        y_all = mean(grouped(exp).(param).food.all(timeROI,:),1,'omitnan');
    else
        y_all = mean(grouped(exp).(param).all(timeROI,:),1,'omitnan');
    end
    plotData = autoCat(plotData, y_all',false);
end

temp_list = [15 17 20 23 25 27 30 33]; % temps that we have temp hold data for...
y_avg = mean(plotData, 1,'omitnan');
y_sem = std(plotData,0,1,'omitnan')./sqrt(num.trial);

sz = 35;
buff = 0.4;
if autoLims
    ylims = [];
else
    ylims = [0 90];
end
xlims = [14, 36];

kolor = foreColor;
LW = 2;

fig = getfig('',1);
% cooling cooling orientation 
subplot(r,c,1); hold on
    for exp = 1:num.exp
        temp = temp_list(exp);
        y = plotData(:,exp);
        y(isnan(y)) = [];
        x = shuffle_data(linspace(temp-buff, temp+buff, length(y)));
        scatter(x,y,sz, Color('grey'), 'filled')
        scatter(temp, y_avg(exp), 70,foreColor,'filled', 'square')
        errorbar(temp, y_avg(exp), y_sem(exp),'Color', foreColor, 'LineWidth',LW)
    end
    plot(temp_list, y_avg,'color', foreColor, 'linewidth', LW, 'linestyle', ':')
    plot_error_fills(true, temp_list, y_avg, y_sem, foreColor,fig_type);
    %formatting
    set(gca,'xdir', 'reverse')
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('static')
    if autoLims
        ylims = [ylims; ylim];
    end

% warming orientation
subplot(r,c,2); hold on
    for exp = 1:num.exp
        temp = temp_list(exp);
        y = plotData(:,exp);
        y(isnan(y)) = [];
        x = shuffle_data(linspace(temp-buff, temp+buff, length(y)));
        scatter(x,y,sz, Color('grey'), 'filled')
        scatter(temp, y_avg(exp), 70,foreColor,'filled', 'square')
        errorbar(temp, y_avg(exp), y_sem(exp),'Color', foreColor, 'LineWidth',LW)
    end
    plot(temp_list, y_avg,'color', foreColor, 'linewidth', LW, 'linestyle', ':')
    plot_error_fills(true, temp_list, y_avg, y_sem, foreColor,fig_type);
    %formatting
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('static')
    if autoLims
        ylims = [ylims; ylim];
    end
formatFig(fig, blkbgd,[r c]);
for i = 1:2
    subplot(r,c,i); hold on
    if autoLims
        ylim([min(ylims(:,1)),max(ylims(:,2))]);
    else 
        ylim(ylims)
    end
end

save_figure(fig,[figDir, param ' occ temp tuning curve scatter'],fig_type);


% figure;
% hold on
% for exp = 1:num.exp
%     plot(grouped(exp).(param).food.avg(timeROI))
% end

%% Static vs dynamic F LRR & waxed comparision figures: 
clearvars('-except',initial_vars{:})

% scatter plot of 17 and 25 degrees for the holds vs waxed and intact trials

[title_str, param,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);
foreColor = formattingColors(blkbgd);
% r = 1;
% c = 2;
autoLims = true;

timeROI = 1:data(1).fps * 3600 * 4.5 ; % time region to average over...4.5 hours

% pull hold trial information 
hold_trials = [3, 4];
dyn_trials = [1, 2];

% 
% plotData = [];
% for i = 1:length(hold_trials)
%     exp = hold_trials(i);
%     if ext
%         y_all = mean(grouped(exp).(param).food.all(timeROI,:),1,'omitnan');
%     else
%         y_all = mean(grouped(exp).(param).all(timeROI,:),1,'omitnan');
%     end
%     plotData = autoCat(plotData, y_all',false);






plotData = [];
for exp = 1:num.exp
    if ext
        y_all = mean(grouped(exp).(param).food.all(timeROI,:),1,'omitnan');
    else
        y_all = mean(grouped(exp).(param).all(timeROI,:),1,'omitnan');
    end
    plotData = autoCat(plotData, y_all',false);
end

temp_list = [15 17 20 25 27 33 35]; % temps that we have temp hold data for...
y_avg = mean(plotData, 1,'omitnan');
y_sem = std(plotData,0,1,'omitnan')./sqrt(num.trial);

sz = 35;
buff = 0.4;
if autoLims
    ylims = [];
else
    ylims = [0 90];
end
xlims = [14, 36];

kolor = foreColor;
LW = 2;

fig = getfig('',1);
% cooling cooling orientation 
subplot(r,c,1); hold on
    for exp = 1:num.exp
        temp = temp_list(exp);
        y = plotData(:,exp);
        y(isnan(y)) = [];
        x = shuffle_data(linspace(temp-buff, temp+buff, length(y)));
        scatter(x,y,sz, Color('grey'), 'filled')
        scatter(temp, y_avg(exp), 70,foreColor,'filled', 'square')
        errorbar(temp, y_avg(exp), y_sem(exp),'Color', foreColor, 'LineWidth',LW)
    end
    plot(temp_list, y_avg,'color', foreColor, 'linewidth', LW, 'linestyle', ':')
    plot_error_fills(true, temp_list, y_avg, y_sem, foreColor,fig_type);
    %formatting
    set(gca,'xdir', 'reverse')
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('static')
    if autoLims
        ylims = [ylims; ylim];
    end

% warming orientation
subplot(r,c,2); hold on
    for exp = 1:num.exp
        temp = temp_list(exp);
        y = plotData(:,exp);
        y(isnan(y)) = [];
        x = shuffle_data(linspace(temp-buff, temp+buff, length(y)));
        scatter(x,y,sz, Color('grey'), 'filled')
        scatter(temp, y_avg(exp), 70,foreColor,'filled', 'square')
        errorbar(temp, y_avg(exp), y_sem(exp),'Color', foreColor, 'LineWidth',LW)
    end
    plot(temp_list, y_avg,'color', foreColor, 'linewidth', LW, 'linestyle', ':')
    plot_error_fills(true, temp_list, y_avg, y_sem, foreColor,fig_type);
    %formatting
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('static')
    if autoLims
        ylims = [ylims; ylim];
    end
formatFig(fig, blkbgd,[r c]);
for i = 1:2
    subplot(r,c,i); hold on
    if autoLims
        ylim([min(ylims(:,1)),max(ylims(:,2))]);
    else 
        ylim(ylims)
    end
end



%% Temp rate comparisons: temp tuning curves for plotting separeated by heating and cooling
clearvars('-except',initial_vars{:})
autoLims = false;
[title_str, param,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);

r = 1;
c = 2;

LW = 2;
FA = 0.5;
xlims = [15, 35];
if autoLims
    ylims = [];
else
    ylims = [0 60];
end

fig = getfig('',1);
% cooling temperatures 
subplot(r,c,1); hold on
    for i = 1:num.exp
        exp = expOrder(i);
        % food distribution 
        if ext
            baseY = grouped(exp).(param).food;
            x = grouped(exp).(param).food.temps;
        else
            baseY = grouped(exp).(param);
            x = grouped(exp).(param).temps;
        end
        y = baseY.decreasing.avg;
        try y_err = baseY.decreasing.std./sqrt(num.trial(exp));
        catch
            y_err = baseY.decreasing.err./sqrt(num.trial(exp));
        end
        plot_error_fills(true, x, y, y_err, grouped(exp).color,fig_type,FA);
        plot(x, y, 'color', grouped(exp).color,'linewidth', LW)
    end
    %formatting
    set(gca,'xdir', 'reverse')
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('cooling')
    if autoLims
        ylims = [ylims; ylim];
    end

% warming temperatures 
subplot(r,c,2); hold on
    for i = 1:num.exp
        exp = expOrder(i);
        % food distribution 
        if ext
            baseY = grouped(exp).(param).food;
            x = grouped(exp).(param).food.temps;
        else
            baseY = grouped(exp).(param);
            x = grouped(exp).(param).temps;
        end
        y = baseY.increasing.avg;
        try y_err = baseY.increasing.std./sqrt(num.trial(exp));
        catch
            y_err = baseY.increasing.err./sqrt(num.trial(exp));
        end
        plot_error_fills(true, x, y, y_err, grouped(exp).color,fig_type,FA);
        plot(x, y, 'color', grouped(exp).color,'linewidth', LW)
    end
    % formatting
    xlabel('temperature')
    xlim(xlims)
    ylabel([param ' %'])
    title('warming')       
    if autoLims
        ylims = [ylims; ylim];
    end
formatFig(fig, blkbgd,[r c]);
for i = 1:2
    subplot(r,c,i); hold on
    if autoLims
        ylim([min(ylims(:,1)),max(ylims(:,2))]);
    else 
        ylim(ylims)
    end
end
set(gca,'ycolor', 'none')

save_figure(fig,[figDir, param ' occ temp tuning curve separated heating and cooling'],fig_type);


% Temp-occupancy correlation?













    













































































