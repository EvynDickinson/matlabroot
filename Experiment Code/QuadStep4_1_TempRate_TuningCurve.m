
% Temp-rate distance tuning curve for a single genotype or across genotypes
% Show how a single genotypes compares in food attraction:temperature for
% different rates of temperature change

% THIS ASSUMES ALL LOADED GROUPS HAVE THE SAME WITHIN-GROUPING TEMPERATURE
% PROTOCOLS
% DOESN'T WORK FOR TEMP PROTOCOLS WITH MORE THAN 1 HEATING AND COOLING TEMP
% RATE

%% Select data groups to compare

clear; close all; clc

% Select processed data structures to compare:
baseFolder = getCloudPath;  
structFolder = [baseFolder 'Data structures\'];
list_dirs = dir(structFolder);
list_dirs = {list_dirs(:).name};
list_dirs(1:2) = [];
expIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'multiple','ListSize',[300,450]);
expNames = list_dirs(expIdx); %name of experiment groups selected
num.exp = length(expIdx);  %number of groups selected

% Load selected experiment data groups:plotData
for i = 1:num.exp
    data(i) = load([structFolder expNames{i} '\' expNames{i} ' post 3.1 data.mat']);
end

clear list_dirs expIdx dirIdx
% Set up base variables
initial_vars = who;
initial_vars = [initial_vars(:); 'initial_vars'; 'grouped'; 'expGroup'; 'saveDir'; 'mat';'expOrder'];
initial_vars = unique(initial_vars);

% Save data / make new grouped data folder
switch questdlg('Select data saving format:','','new structure','existing structure', 'cancel','new structure');
    case 'new structure'
        expGroup = char(inputdlg('Structure name:'));
        saveDir = [baseFolder 'Grouped Data Structures\' expGroup '\'];
        if ~exist(saveDir,'dir')
            mkdir(saveDir);
        end
        save([saveDir expGroup ' data.mat'],'-v7.3');
        disp([expGroup ' saved'])
    case 'existing structure'
        list_dirs = dir([baseFolder 'Grouped Data Structures\']);
        list_dirs = {list_dirs(:).name};
        list_dirs(1:2) = [];
        dirIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[300,450]);
        expGroup = list_dirs{dirIdx}; %name of experiment groups selected
        saveDir = [baseFolder 'Grouped Data Structures\' expGroup '\'];
        save([saveDir expGroup ' data.mat'],'-v7.3');
        disp([expGroup ' saved'])
    case 'cancel'
        return
end


%% ANALYSIS: get average timecourse traces for each group

clearvars('-except',initial_vars{:})
grouped = struct;

if strcmp(expGroup,'WT linear recovery caviar')
    expOrder = [];
    expList = {'Berlin WT','CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'};
    colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
    for ii = 1:num.exp
        expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar']));
%         grouped(expOrder(ii)).color = Color(colors{(ii)});
    end
else
    expOrder = 1:num.exp;
    colors = {'BlueViolet','red','white','Turquoise','Gold','pink','Orange'};
    % colors = {'red','yellow','dodgerblue','Gold','pink','Orange'};
    % colors = {'BlueViolet','white','turquoise','Gold','pink','Orange'};
    % colors = {'BlueViolet','gold','white','turquoise','pink','Orange'};
    % colors = {'Teal','gold','white', 'magenta','dodgerblue','Orange'};
    % colors = {'White','magenta','dodgerblue','Orange'};
end


for i = 1:num.exp % FOR EACH DATA GROUP
    % GENERAL
    grouped(i).name = data(i).ExpGroup;
    grouped(expOrder(i)).color = Color(colors{(i)});

    % TIME COURSE DATA
    num.trial(i) = data(i).ntrials;
    [time,temp,speed,distance] = deal([]);
    for trial = 1:num.trial(i)
        time = autoCat(time,data(i).data(trial).occupancy.time,false,true);
        temp = autoCat(temp,data(i).data(trial).occupancy.temp,false,true);
        speed = autoCat(speed,data(i).data(trial).speed.avg,false,true);
    % movement = autoCat(speed,data(i).data(trial).speed.avg,false,true);
        distance = autoCat(distance,data(i).data(trial).dist2wells(:,data(i).T.foodLoc(trial)),...
                   false,true);
    end
    grouped(i).time = mean(time,2,'omitnan');
    grouped(i).temp = mean(temp,2,'omitnan');
    grouped(i).speed.all = speed;
    grouped(i).speed.avg = mean(speed,2,'omitnan');
    grouped(i).speed.err = std(speed,0,2,'omitnan')/sqrt(num.trial(i));
    grouped(i).speed.zscore.all = zscore(grouped(i).speed.all);
    grouped(i).speed.zscore.avg = mean(grouped(i).speed.zscore.all,2,'omitnan');
    grouped(i).dist.all = distance;
    grouped(i).dist.avg = mean(distance,2,'omitnan');
    grouped(i).dist.err = std(distance,0,2,'omitnan')/sqrt(num.trial(i));
    grouped(i).dist.zscore.all = zscore(grouped(i).dist.all);
    grouped(i).dist.zscore.avg = mean(grouped(i).dist.zscore.all,2,'omitnan');

    % AVG POSITION BINNED BY TEMP
    for trial = 1:num.trial(i)
        tempList = data(i).G(trial).TR.temps;
        y = data(i).G(trial).TR.data(:,3);
        for tt = 1:length(tempList)
            grouped(i).dist.tempBinned(trial,tt) = mean(y(data(i).G(trial).TR.tempIdx==tt),'omitnan');
        end
        grouped(i).dist.tempList(trial,:) = tempList;
    end
    grouped(i).dist.distavgbytemp = [grouped(i).dist.tempList(1,:)',...
                                     mean(grouped(i).dist.tempBinned,1,'omitnan')'];
    grouped(i).dist.distavgbytemp_err = [grouped(i).dist.tempList(1,:)',...
                                         (std(grouped(i).dist.tempBinned,0,'omitnan')/sqrt(num.trial(i)))'];

    % BINNED
    [tempRates,decreasing,increasing,temperatures] = deal([]);
    for trial = 1:num.trial(i)
        % Account for multiple numbers of rate trials
        rates = data(i).G(1).TR.rates;
        nRates(i) = size(rates,2);
        if nRates(i)==2 
            blankdata = data(i).G(trial).TR.dist_mat.avg;
        elseif nRates == 3
            blankdata = data(i).G(trial).TR.dist_mat.avg;
            blankdata(rates==0,:) = [];
            rates(rates==0) = [];
        else 
            warndlg('Temp protocol has too many rates for this analysis')
            return
        end
        downIdx = find(rates<0);
        upIdx = find(rates>0);
        decreasing(:,trial) = blankdata(downIdx,:);
        increasing(:,trial) = blankdata(upIdx,:);
        tempRates = autoCat(tempRates,rates,true,true);
        temperatures(:,trial) = data(i).G(trial).TR.temps;
    end
    grouped(i).increasing.temps = median(temperatures,2);
    grouped(i).increasing.all = increasing;
    grouped(i).increasing.avg = mean(increasing,2,'omitnan');
    grouped(i).increasing.zscore.all = zscore(increasing,0,'omitnan');
    grouped(i).increasing.zscore.avg = mean(grouped(i).increasing.zscore.all,2,'omitnan');
    grouped(i).increasing.err = std(increasing,0,2,'omitnan')/sqrt(num.trial(i)); 
    grouped(i).increasing.rate = median(tempRates(:,upIdx));
    grouped(i).decreasing.temps = median(temperatures,2);
    grouped(i).decreasing.all = decreasing;
    grouped(i).decreasing.avg = mean(decreasing,2,'omitnan');
    grouped(i).decreasing.zscore.all = zscore(decreasing,0,'omitnan');
    grouped(i).decreasing.zscore.avg = mean(grouped(i).decreasing.zscore.all,2,'omitnan');
    grouped(i).decreasing.err = std(decreasing,0,2,'omitnan')/sqrt(num.trial(i)); 
    grouped(i).decreasing.rate = median(tempRates(:,downIdx));

end

% RAMP-TO-RAMP COMPARISONS
binWidth = 0.5; %temp bin increment
mat = struct;
% Extract ramp by ramp information from each trial
for i = 1:num.exp
    for trial = 1:num.trial(i)
        tempPoints = getTempTurnPoints(data(i).T.TempProtocol{trial});
        foodLoc = data(i).T.foodLoc(trial);
        threshLow = tempPoints.threshLow;
        threshHigh = tempPoints.threshHigh;
        dist = data(i).data(trial).occupancy.dist2wells(:,foodLoc);
        temp = data(i).data(trial).occupancy.temp;
        tempBins = floor(threshLow):binWidth:ceil(threshHigh); 
        if tempBins(end)<ceil(threshHigh)
            tempBins(end+1) = ceil(threshHigh) + binSpace;
        end

        for idx = 1:tempPoints.nDown
            % COOLING
            downROI = tempPoints.down(idx,1):tempPoints.down(idx,2);
            mat(i).cooling(idx).dist(:,trial) = dist(downROI);
            mat(i).cooling(idx).temp(:,trial) = temp(downROI);
            mat(i).cooling(idx).tempIdx(:,trial) = discretize(temp(downROI),tempBins);

            % HEATING
            upROI = tempPoints.up(idx,1):tempPoints.up(idx,2);
            mat(i).heating(idx).dist(:,trial) = dist(upROI);
            mat(i).heating(idx).temp(:,trial) = temp(upROI);
            mat(i).heating(idx).tempIdx(:,trial) = discretize(temp(upROI),tempBins);

            % HYSTERESIS MEASURE
            for n = 1:length(tempBins)-1 %pull avg temp for this particular cooling ramp
                %cooling
                C_loc = mat(i).cooling(idx).tempIdx(:,trial)==n;
                mat(i).cooling(idx).bintemp(n,trial) = mean(mat(i).cooling(idx).dist(C_loc,trial),'omitnan');
                %heating
                H_loc = mat(i).heating(idx).tempIdx(:,trial)==n;
                mat(i).heating(idx).bintemp(n,trial) = mean(mat(i).heating(idx).dist(H_loc,trial),'omitnan');
            end
            mat(i).hysteresis(idx).all(:,trial) = mat(i).cooling(idx).bintemp(:,trial)-mat(i).heating(idx).bintemp(:,trial);
            mat(i).hysteresis(idx).sum(trial) = sum(mat(i).hysteresis(idx).all(:,trial),'omitnan');
            mat(i).hysteresis(idx).tempbins = tempBins;
            mat(i).distHist(:,trial,idx) = mat(i).hysteresis(idx).all(:,trial); %hyst by temp bin
            mat(i).cumHist(trial,idx) = mat(i).hysteresis(idx).sum(trial); %cumulative hysteresis
        end
    end
end


disp('Next')



%% FIGURE: Basic over-lap of time-trials and temperature protocols
clearvars('-except',initial_vars{:})
plot_err = false;

blkbgd = true;



% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8]; %distance from food timecourse %TODO: normalize this to something more intuitive? 
sb(3).idx = [10,11,13,14]; %speed timecourse
sb(4).idx = 3:3:15; %binned distance alignment

LW = 0.75;
sSpan = 180;

if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end

% FIGURE:
fig = figure; set(fig,'color','w',"Position",[1934 468 1061 590])
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',LW,'Color',kolor)
    
    %distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).dist.avg,'moving',sSpan);
        y_err = smooth(grouped(i).dist.err,'moving',sSpan);
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        end
        plot(x,y,'LineWidth',LW,'Color',kolor)

%     %speed
    subplot(r,c,sb(3).idx); hold on
        y = smooth(grouped(i).speed.avg,'moving',sSpan);
        y_err = smooth(grouped(i).speed.err,'moving',sSpan);
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        end
        plot(x,y,'LineWidth',LW,'Color',kolor)

    %temp rate depended distance
    subplot(r,c,sb(4).idx); hold on
        x = grouped(i).dist.distavgbytemp(:,1);
        y = grouped(i).dist.distavgbytemp(:,2);
        y_err = grouped(i).dist.distavgbytemp_err(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];

        plot(x,y,'color',kolor,'linewidth',LW+1)
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
        dataString{i} = grouped(i).name;

%         %increasing 
%         x = grouped(i).increasing.temps;
%         y = grouped(i).increasing.avg;
%         y_err = grouped(i).increasing.err;
%         loc = isnan(y) | isnan(y_err);% remove nans 
%         y(loc) = []; x(loc) = []; y_err(loc) = [];
% 
%         if plot_err
%             fill_data = error_fill(x, y, y_err);
%             h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
%             set(h, 'facealpha', 0.2)
%         end
%         plot(x,y,'LineWidth',LW+0.5,'Color',kolor,'linestyle','-')
%         %decreasing 
%         x = grouped(i).decreasing.temps;
%         y = grouped(i).decreasing.avg;
%         y_err = grouped(i).decreasing.err;
%         loc = isnan(y) | isnan(y_err);% remove nans 
%         y(loc) = []; x(loc) = []; y_err(loc) = [];
% 
%         if plot_err
%             fill_data = error_fill(x, y, y_err);
%             h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
%             set(h, 'facealpha', 0.2)
%         end
%         plot(x,y,'LineWidth',LW+.5,'Color',kolor,'linestyle','--','HandleVisibility','off');
% 
%         % Names and Colors of included data
%         dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx) 
ylabel('\circC')
set(gca,"XColor",backColor)
% distance
subplot(r,c,sb(2).idx) 
ylabel('distance (mm)')
set(gca,"XColor",backColor)
% speed
subplot(r,c,sb(3).idx) 
ylabel('speed (mm/s)')
xlabel('time (min)')
% temp rate 
subplot(r,c,sb(4).idx) 
ylabel('distance (mm)')
xlabel('temp (\circC)')
% 
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' timecourse summary'],'-png');

%% FIGURE: cumulative hysteresis for each genotype / trial

clearvars('-except',initial_vars{:})
LW = 0.75;
buff = 0.2;
SZ = 50;
r = 1; %rows
c = 3; %columns
plot_err = false;


% FIGURE:
fig = figure; set(fig,'color','w',"Position",[1932 690 1050 438])
% Hystersis
subplot(r,c,1)
hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    %increasing 
    x = grouped(i).increasing.temps;
    y = grouped(i).increasing.avg;
    y_err = grouped(i).increasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans 
    y(loc) = []; x(loc) = []; y_err(loc) = [];

    if plot_err
        fill_data = error_fill(x, y, y_err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
        set(h, 'facealpha', 0.2)
    end
    plot(x,y,'LineWidth',LW+0.5,'Color',kolor,'linestyle','-')
    %decreasing 
    x = grouped(i).decreasing.temps;
    y = grouped(i).decreasing.avg;
    y_err = grouped(i).decreasing.err;
    loc = isnan(y) | isnan(y_err);% remove nans 
    y(loc) = []; x(loc) = []; y_err(loc) = [];

    if plot_err
        fill_data = error_fill(x, y, y_err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
        set(h, 'facealpha', 0.2)
    end
    plot(x,y,'LineWidth',LW+.5,'Color',kolor,'linestyle','--','HandleVisibility','off');

    % Names and Colors of included data
    dataString{i} = grouped(i).name;
end
subplot(r,c,1) 
ylabel('distance (mm)')
xlabel('temp (\circC)')

% Pull difference in distance heating-cooling
subplot(r,c,2)
hold on
for i = 1:num.exp
    x = repmat(grouped(i).decreasing.temps,[1,num.trial(i)]);
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    kolor = grouped(i).color;
%     plot(x,y,'color',kolor,'LineWidth',LW); 
    plot(mean(x,2),mean(y,2),'color',kolor,'LineWidth',2)
end
h_line(0,'w',':',1)
xlabel('temp (\circC)')
ylabel('distance difference (mm)')

% Cumulative difference in proximity
subplot(r,c,3)
hold on
for ii = 1:num.exp
    i = expOrder(ii);

    kolor = grouped(i).color;
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    plotY = sum(y,1,'omitnan');
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i))); 
    scatter(x,plotY,SZ,kolor,"filled","o")
    plot([ii-buff,ii+buff],[mean(plotY),mean(plotY)],'color','w','LineWidth',2)
end
xlim([0.5,num.exp+0.5])
h_line(0,'w',':',1)
ylabel('cumulatice difference (mm)')

formatFig(fig,true,[r,c]);
set(gca,'XTick',[],'xcolor','k')
xlabel('Group','color','w')

% save figure
save_figure(fig,[saveDir expGroup ' hysteresis summary'],'-png');

%% FIGURE: ramp by ramp hysteresis comparision

clearvars('-except',initial_vars{:})


% FIGURE: plot the heating vs cooling plots for each of the four ramps
% acros each of the different trials
buff = 0.3;
LW = 0.75;

fig = figure; set(fig,'pos',[1998 582 431 497]); hold on
    for ii = 1:num.exp
        i = expOrder(ii);
        kolor = grouped(i).color;
        nRamps = size(mat(i).hysteresis,2);
        x = repmat(linspace(ii-buff,ii+buff,nRamps),[num.trial(i),1]);
        plot(x',mat(i).cumHist','color',kolor,'linewidth',LW,...
             'Marker','o', 'MarkerFaceColor',kolor)
    end

h_line(0,'w',':',1)
formatFig(fig,true);
set(gca,'XTick',1:num.exp)
xlabel('Group')
ylabel('Cumulative dist difference (cool-heat)')

% save figure
save_figure(fig,[saveDir expGroup ' ramp by ramp cumulative hysteresis'],'-png');


%% ANALYSIS AND FIGURES: Event-aligned comparisons

clearvars('-except',initial_vars{:})
% ylimits = [-8,3];
ylimits = [-15,15];
timeROI = 60; % how many minutes to look at behavior after each event
duration = ceil(timeROI*3*60);

sections = {'increasing','decreasing','holding'};
s_color = {'red','dodgerblue','white'};

for i = 1:num.exp
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    for ss = 1:length(sections)
        switch sections{ss}
            case 'increasing'
                tpBin = 'up';
                nrr = tp.nUp;
            case 'decreasing'
                tpBin = 'down';
                nrr = tp.nDown;
            case 'holding'
                tpBin = 'hold';
                nrr = tp.nHold;
        end
        temp = [];
        for rr = 1:nrr
            ROI = tp.(tpBin)(rr,1):tp.(tpBin)(rr,1)+duration;
            temp(:,:,rr) = grouped(i).dist.all(ROI,:);
        end
    
        temp_norm = temp-mean(temp(1:10,:,:),'omitnan');
        temp_avg = mean(temp_norm,3);
        % add to the grouped data
        grouped(i).aligned.([sections{ss} '_avg']) = temp_avg;
        grouped(i).aligned.([sections{ss} '_norm']) = temp_norm;
        grouped(i).aligned.([sections{ss} '_all']) = temp;
        grouped(i).aligned.([sections{ss} '_SEM']) = std(temp_avg,0,2,'omitnan')/sqrt(num.trial(i));
        grouped(i).aligned.([sections{ss} '_MEAN']) = mean(temp_avg,2,'omitnan');
    end
end


% % FIGURES: SINGLE EXPERIMENT COMPARISON
% dispROI = 15;
% duration = ceil(dispROI*3*60);
% x = linspace(0,dispROI,duration+1);
% LW = 1.5;
% 
% for i = 1:num.exp
% 
%     fig = figure; set(fig,'pos',[2130 275 428 534])
%     hold on
%     for ss = 1:length(sections)
%         y = grouped(i).aligned.([sections{ss} '_MEAN'])(1:duration+1);
%         y_err = grouped(i).aligned.([sections{ss} '_SEM'])(1:duration+1);
%     
%         fill_data = error_fill(x, y, y_err);
%         h = fill(fill_data.X, fill_data.Y, Color(s_color{ss}), 'EdgeColor','none','HandleVisibility','off');
%         set(h, 'facealpha', 0.4)
%         plot(x, y,'color',Color(s_color{ss}),'LineWidth',LW)
%     end
%     
%     xlabel('time (min)')
%     ylabel('distance from food (mm)')
%     formatFig(fig,true);
%     title([grouped(i).name],'Color','w','FontSize',12,'FontName','times')
%     
%     save_figure(fig,[saveDir grouped(i).name...
%                 ' event aligned distance -duration ' num2str(dispROI) ' min'],...
%                 '-png',true);
% end

% FIGURE: CROSS EXPERIMENT COMPARISION (WITHIN HEATING,COOLING,HOLDING)
r = 1;
c = 3;

dispROI = 50;
duration = ceil(dispROI*3*60);
x = linspace(0,dispROI,duration+1);
LW = 1.5;
SEM_shading = false;
sSpan = 1;


fig = figure; set(fig,'pos',[1932 690 1050 438])
for ss = 1:length(sections) %
    subplot(r,c,ss); hold on
    for i = 1:num.exp
        y = grouped(i).aligned.([sections{ss} '_MEAN'])(1:duration+1);
        if SEM_shading
            y_err = grouped(i).aligned.([sections{ss} '_SEM'])(1:duration+1);
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y,grouped(i).color, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.4)
        end
        plot(x, smooth(y,'moving',sSpan),'color',grouped(i).color,'LineWidth',LW)
    end
    xlabel('time (min)')
    ylabel('distance from food (mm)')
    title(sections{ss})
    ylim(ylimits)
end
formatFig(fig,true,[r,c]);

save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' ...
            num2str(sSpan) ' duration ' num2str(dispROI) ' min'],'-png',true);
    


%% ANALYSIS AND FIGURES: COM of postions for each trial
% vectors to fly 'mass' in arena at different temperatures from food

clearvars('-except',initial_vars{:})


% 1) for a given experiment, find the avg fly position for each frame. 
% 2) plot the avg position on a 'blank' arena with the food loc indicated
% color coded by temperature at the time
% 3) add the avg vector line for binned temperatures
% 4) figure out how to rotate the arena or overlay vectors to collapse
% across trials within an experimental group e.g. berlin 

 
switch questdlg('Save trial by trial figures?','','yes','no','cancel','no')
    case 'yes' 
        saveFig = true;
    case 'no'
        saveFig = false;
    case 'cancel'
        return
end

for i = 1:num.exp 
  for trial = 1:num.trial(i)

    % get arena information
    well_loc = data(i).T.foodLoc(trial);
    C = data(i).data(trial).data.centre;
    r = data(i).data(trial).data.r;
    well_C = data(i).data(trial).data.wellcenters(:,well_loc);
    arena_C = data(i).data(trial).data.wellcenters(:,5);   
    
    % positions for all 'flies' over time
    x = data(i).data(trial).data.x_loc; 
    y = data(i).data(trial).data.y_loc;
    x_avg = mean(x,2,'omitnan');
    y_avg = mean(y,2,'omitnan');
    % temperature for each frame
    temperature = data(i).data(trial).occupancy.temp;
    temp = [temperature,x_avg,y_avg]; 
    
    % smooth into 1 minute bins
    smoothData = [];
    bin = 1*60*3;
    for c = 1:size(temp,2)
        smoothData(:,c) = smooth(temp(:,c),'moving',bin);
    end
    dsData = smoothData(1:bin:end,:); % select one point for each minute
    
    % get temp color distribution (1/4 degree incrememts)
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    cMapRange = tp.threshLow:0.25:tp.threshHigh;
    dsData(:,4) = discretize(dsData(:,1),cMapRange); %binned temp category
    dsData(isnan(dsData(:,4)),:) = []; %remove data points outside allowed range
    %avg position for each binned temp
    position = [];
    ntemps = length(cMapRange)-1;
    for g = 1:ntemps
        loc = dsData(:,4)==g;
        position(g,:) = mean(dsData(loc,2:3),1,'omitnan');
    end
    
    g1 = floor(ntemps/2);
    g2 = ntemps-g1;
    g1_cMap = Color('darkblue','grey',g1); %deepskyblue
    g2_cMap = Color('grey','red',g2);
    cMap = [g1_cMap;g2_cMap];

    kolors = cMap(dsData(:,4),:); 
    
    % save data into structure:
    mat(i).position(trial).data = position;
    mat(i).position(trial).tempBins = cMapRange;
    mat(i).position(trial).cMap = cMap;
    mat(i).position(trial).wells = data(i).data(trial).data.wellcenters;
    mat(i).position(trial).wellLoc = well_loc;
    

    % Plot the avg position within in the arena
    if saveFig
        SZ = 40;
        fig = figure; set(fig, 'pos',[1970 552 484 384]); 
        hold on
            scatter(dsData(:,2), dsData(:,3), 5, kolors,'filled')
            for well = 1:4
                C = data(i).data(trial).data.wellcenters(:,well);
                scatter(C(1),C(2),SZ,Color('grey'),'filled')
            end
            scatter(well_C(1),well_C(2),SZ,'y','filled')
            viscircles([arena_C(1),arena_C(2)],r,'Color','w');
            scatter(position(:,1),position(:,2),15,cMap,'filled','o','MarkerEdgeColor','w')
            axis equal square;
            formatFig(fig,true);
            set(gca,'XColor','k','YColor','k');
            % set color bar information
            colorData = uint8(round(cMap.*255)); % convert color map to uint8
            colormap(colorData);
            c = colorbar('color','w');
            set(c,'Ticks',[0,1],'TickLabels',[tp.threshLow,tp.threshHigh]);
            c.Label.String = 'Temperature (\circC)';
            c.Label.VerticalAlignment = "bottom";
        title([data(i).T.Genotype{trial} ' ' data(i).T.Date{trial}],'color','w')
        
        figFolder = [data(i).figDir 'avg position by trial/'];
        if ~exist(figFolder,'dir')
            mkdir(figFolder);
        end
        save_figure(fig,[figFolder data(i).T.Date{trial} ' ' data(i).T.ExperimentID{trial} ' ' data(i).T.Arena{trial}],'-png',true);
    end
  end
end




%% FIGURE: Register position COM to common frame

clearvars('-except',initial_vars{:})
blkbgk = true;
if blkbgk
    backColor = 'k';
    foreColor = 'w';
else
    backColor = 'w';
    foreColor = 'k';
end

SZ = 40;


for i = 1:num.exp
    fig = figure; hold on
    for trial = 1:num.trial(i)
        
        wells = mat(i).position(trial).wells;
        wellLoc = mat(i).position(trial).wellLoc;
        x = mat(i).position(trial).data(:,1);
        y = mat(i).position(trial).data(:,2);
        kolor = mat(i).position(trial).cMap;
        
        % Make food well the origin
        x_offset = wells(1,wellLoc);
        y_offset = wells(2,wellLoc);
        wells_x = wells(1,:)-x_offset;
        wells_y = wells(2,:)-y_offset;
        X = x-x_offset;
        Y = y-y_offset;
        
        % Rotate to correct orientation
        switch wellLoc
            case 1
                plotData(:,1) = Y;
                plotData(:,2) = -X;
                WELLS(:,1) = wells_y;
                WELLS(:,2) = -wells_x;
            case 2 
                plotData(:,1) = X;
                plotData(:,2) = -Y;
                WELLS(:,1) = wells_x;
                WELLS(:,2) = -wells_y;
            case 3
                plotData(:,1) = -Y;
                plotData(:,2) = X;
                WELLS(:,1) = -wells_y;
                WELLS(:,2) = wells_x;
            case 4 
                plotData(:,1) = X;
                plotData(:,2) = Y;
                WELLS(:,1) = wells_x;
                WELLS(:,2) = wells_y;
        end
        
        % PLOT
        scatter(WELLS(1:4,1),WELLS(1:4,2),SZ,foreColor,'filled')
        scatter(WELLS(wellLoc,1),WELLS(wellLoc,2),SZ,'green','filled')
        scatter(plotData(:,1),plotData(:,2),15,kolor,'filled')
        
    end
    
%     viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color',grouped(i).color);
        viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color','k');

    
    axis square; 
    axis equal;
    formatFig(fig,blkbgk);
    
    set(gca,'XColor',backColor,'YColor',backColor);
    % set color bar information
    colorData = uint8(round(kolor.*255)); % convert color map to uint8
    colormap(colorData);
    c = colorbar('color',foreColor);
    set(c,'Ticks',[0,1],'TickLabels',[mat(i).position(trial).tempBins(1),mat(i).position(trial).tempBins(end)]);
    c.Label.String = 'Temperature (\circC)';
    c.Label.VerticalAlignment = "bottom";
    title([data(i).ExpGroup],'color',foreColor)
    save_figure(fig,[saveDir expGroup grouped(i).name ' COM position'],'-png',true);

end


%% FIGURE: Normalized increasing/decreasing temp responses

clearvars('-except',initial_vars{:})

plot_err = false;
blkbgd = true;

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse
sb(3).idx = 3:3:15; %binned distance alignment

LW = 1.5;
sSpan = 180;

if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end

% FIGURE:
fig = figure; set(fig,'color','w',"Position",[1934 468 1061 590])
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',LW,'Color',kolor)
    
%     %distance
%     subplot(r,c,sb(2).idx); hold on
%         y = smooth(grouped(i).dist.avg,'moving',sSpan);
%         y_err = smooth(grouped(i).dist.err,'moving',sSpan);
%         if plot_err
%             fill_data = error_fill(x, y, y_err);
%             h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%             set(h, 'facealpha', 0.2)
%         end
%         plot(x,y,'LineWidth',LW,'Color',kolor)

    %zscore distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).dist.zscore.avg,'moving',sSpan);
%         y_err = smooth(grouped(i).speed.err,'moving',sSpan);
%         if plot_err
%             fill_data = error_fill(x, y, y_err);
%             h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%             set(h, 'facealpha', 0.2)
%         end
        plot(x,y,'LineWidth',LW,'Color',kolor)

    %temp rate depended distance
    subplot(r,c,sb(3).idx); hold on
        %increasing 
        x = grouped(i).increasing.temps;
        y = grouped(i).increasing.zscore.avg;
        y_err = grouped(i).increasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans 
        y(loc) = []; x(loc) = []; y_err(loc) = [];

        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
        plot(x,y,'LineWidth',LW,'Color',kolor,'linestyle','-')
        %decreasing 
        x = grouped(i).decreasing.temps;
        y = grouped(i).decreasing.zscore.avg;
        y_err = grouped(i).decreasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans 
        y(loc) = []; x(loc) = []; y_err(loc) = [];

        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
        plot(x,y,'LineWidth',LW,'Color',kolor,'linestyle','--','HandleVisibility','off');

        % Names and Colors of included data
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx) 
ylabel('\circC')
set(gca,"XColor",backColor)
% distance
subplot(r,c,sb(2).idx) 
ylabel('distance zscore')
xlabel('time (min)')

% temp rate 
subplot(r,c,sb(3).idx) 
ylabel('distance (zscore)')
xlabel('temp (\circC)')
% 
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' zscore timecourse summary'],'-png');

%% FIGURE: Position averages across temp (regardless heat/cool)


clearvars('-except',initial_vars{:})

plot_err = false;
blkbgd = true;

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse
sb(3).idx = 3:3:15; %binned distance alignment

LW = 1.5;
SZ = 50;
sSpan = 180;

if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end
buff = 0.2;

% FIGURE:
fig = figure; set(fig,'color','w',"Position",[2108 475 453 590]); hold on
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
    y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
%     datastats.all = autoCat(datastats.all,y,false);
%     datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
    scatter(x,y,SZ,kolor,'filled')
    plot([ii-buff,ii+buff],[mean(y),mean(y)],'color', kolor,'linewidth',LW)
end
set(gca,'xtick',0:num.exp+1,'XTickLabel',[],'XColor',backColor)
xlim([0,num.exp+1])
ylim([0,22])
% xlabel('Fly Strain')
ylabel('temp-driven spatial range (mm)')
formatFig(fig,true);
save_figure(fig,[saveDir expGroup ' distance modulation by temperature'],'-png');

% STATS:
[datastats.all, datastats.id] = deal([]);
for i = 1:num.exp
    y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
    datastats.all = autoCat(datastats.all,y,false);
    datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
end

% determine which groups differ from each other
[~,~,stats] = anova1(datastats.all,datastats.id);
[c,~,~,~] = multcompare(stats);

% bonferonni multiple comparisons correction
alpha = 0.05; %significance level
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nTemp-driven position range statistics\n\n')
[Group1,Group2,P_Value] = deal([]);
idx = 0;
for i = 1:length(significantHypotheses)
    if significantHypotheses(i)
        idx = idx+1;
        Group1{idx,1} = expNames{c(i,1)};
        Group2{idx,1} = expNames{c(i,2)};
        P_Value(idx,1) = c(i,6);
    end
end
sig_comp = table(Group1,Group2,P_Value);
disp(sig_comp)

%% FIGURE & STATS: distance traveled over temperature variation

clearvars('-except',initial_vars{:})
y_limits = [1,25];
plot_err = false;
blkbgd = true;

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse
sb(3).idx = 3:3:15; %binned distance alignment

LW = 1.5;
SZ = 50;
sSpan = 180;

if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end
buff = 0.2;

% FIGURE:
fig = figure; set(fig,'color','w',"Position",[2289 468 706 590]); hold on %[1934 468 1061 590]
for i = 1:num.exp
    kolor = grouped(i).color;
    x = shuffle_data(linspace(i-buff,i+buff,num.trial(i)));
    y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
%     datastats.all = autoCat(datastats.all,y,false);
%     datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
    scatter(x,y,SZ,kolor,'filled')
    plot([i-buff,i+buff],[mean(y),mean(y)],'color', kolor,'linewidth',LW)
end
set(gca,'xtick',0:num.exp+1,'XTickLabel',[])
xlim([0,num.exp+1])
ylim(y_limits)
xlabel('Fly Strain')
ylabel('distance traveled (mm)')
formatFig(fig,true);
save_figure(fig,[saveDir expGroup ' distance modulation by temperature'],'-png');

% STATS:
[datastats.all, datastats.id] = deal([]);
for i = 1:num.exp
    y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
    datastats.all = autoCat(datastats.all,y,false);
    datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
end

[~,~,stats] = anova1(datastats.all,datastats.id);
[c,~,~,~] = multcompare(stats);
stats_tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
disp(table(expNames',(1:num.exp)', 'VariableNames',{'Name', 'ID number'}))
disp(stats_tbl)

% determine which groups differ from each other

%% TODO FIGURE & STATS: Average distance from the food and closest sustained distance from the food (for a temperature bin)
% TODO: combine this with the above section and automate the stats with
% significance lines
clearvars('-except',initial_vars{:})
blkbgd = true;
y_limits = [0,30];




 

plotdata = [];
for i = 1:num.exp
    [M,avg] = deal([]);
    x = grouped(i).increasing.all;
    for trial = 1:num.trial(i)
        temp = smooth(x(:,trial),2,'moving');
        M(trial) = min(temp);
        avg(trial) = mean(temp,'omitnan');
    end
    plotdata(i).min = M;
    plotdata(i).avg_dist = avg;
end

% FIGURE:
r = 1;
c = 2;
LW = 1.5;
SZ = 50;

if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end
buff = 0.2;


fig = figure; set(fig,'color','w',"Position",[2126 458 660 590])
% minimum distance to food (for a 1deg C period)
subplot(r,c,1); hold on
    for i = 1:num.exp
        kolor = grouped(i).color;
        x = shuffle_data(linspace(i-buff,i+buff,num.trial(i)));
        y = plotdata(i).min;
        scatter(x,y,SZ,kolor,'filled')
        plot([i-buff,i+buff],[mean(y),mean(y)],'color', kolor,'linewidth',LW)
    end
    set(gca,'xtick',0:num.exp+1,'XTickLabel',[])
    xlim([0,num.exp+1])
    ylim(y_limits)
    xlabel('Fly Strain')
    ylabel('peak food proximity (mm)')

% minimum distance to food (for a 1deg C period)
subplot(r,c,2); hold on
    for i = 1:num.exp
        kolor = grouped(i).color;
        x = shuffle_data(linspace(i-buff,i+buff,num.trial(i)));
        y = plotdata(i).avg_dist;
        scatter(x,y,SZ,kolor,'filled')
        plot([i-buff,i+buff],[mean(y),mean(y)],'color', kolor,'linewidth',LW)
    end
    set(gca,'xtick',0:num.exp+1,'XTickLabel',[])
    xlim([0,num.exp+1])
    ylim(y_limits)
    xlabel('Fly Strain')
    ylabel('avg distance to food (mm)')


formatFig(fig,true,[r,c]);
save_figure(fig,[saveDir expGroup ' food proximity modulation by temperature'],'-png');


% % STATS:
% [datastats.all, datastats.id] = deal([]);
% for i = 1:num.exp
%     y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
%     datastats.all = autoCat(datastats.all,y,false);
%     datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
% end
% 
% [~,~,stats] = anova1(datastats.all,datastats.id);
% [c,~,~,gnames] = multcompare(stats);
% stats_tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% disp(table(expNames',(1:num.exp)', 'VariableNames',{'Name', 'ID number'}))
% disp(stats_tbl)


%% FIGURE AND ANALYSIS: temp range at minimum and max distances
clearvars('-except',initial_vars{:})

for tt = 1:2
  switch tt
     case 1
        type = 'increasing';
     case 2
        type = 'decreasing';
  end

    % Pull the avg min|max distance and temperature index for each trial
    [min_dist,max_dist,I_min,I_max]  = deal(nan([max(num.trial),num.exp]));
    for i = 1:num.exp
        for trial = 1:num.trial(i)
            y = grouped(i).(type).all(:,trial);
            % Minimum distance
            [min_dist(trial,i), I_min(trial,i)] = min(y); %I is the temp index 
            % Maximum distance
            [max_dist(trial,i), I_max(trial,i)] = max(y);
        end
    end 
    
    % For each trial -- find which temperatures have a mean distance that falls
    % within the SE of the distance for the 'min|max' distance found above
    [minimumTempRangeStart,minimumTempRangeEnd,mimimumTemp] = deal(nan([max(num.trial),num.exp]));
    [maxTempRangeStart,maxTempRangeEnd,maxTemp] = deal(nan([max(num.trial),num.exp]));
    for i = 1:num.exp
        for trial = 1:num.trial(i)
            rateIDX = 3; %3=increasing
            tempList = data(i).G(trial).TR.temps;
            y = data(i).G(trial).TR.data(:,3); % distance from food for this whole trial
            positionAvg = grouped(i).(type).all(:,trial); % avg distance at temps for trial
            
            % ------ minimum temp -------
            % find the SE for the distance at the minumum temperature:
            loc = (data(i).G(trial).TR.tempIdx==I_min(trial,i)) &...
                  (data(i).G(trial).TR.rateIdx==rateIDX);
            testData = y(loc);
            SE = std(testData,'omitnan');
            boundLimit = min_dist(trial,i)+SE; %since this is minumum value, the limit would only exist above
            tempLocWithinBounds = positionAvg<=boundLimit;
            
            % find other temps that have their distance within 1 STD of the min food distance
            upTemp = tempLocWithinBounds(I_min(trial,i):end);
            upLoc = find(diff(upTemp)<0);
            if ~isempty(upLoc)
                tempLocWithinBounds(I_min(trial,i)+upLoc(1):end) = false;
            end
            downTemp = tempLocWithinBounds(1:I_min(trial,i));
            downLoc = find(diff(downTemp)>0);
            if ~isempty(downLoc)
                tempLocWithinBounds(1:downLoc(end)) = false;
            end
            tempswithinBounds = tempList(tempLocWithinBounds);
    
            % Save the temp ranges for minumum distance:
            minimumTempRangeStart(trial,i) = tempswithinBounds(1);
            minimumTempRangeEnd(trial,i) = tempswithinBounds(end);
            mimimumTemp(trial,i) = tempList(I_min(trial,i));
        
             % ------ maximum temp -------
            % find the SE for the distance at the minumum temperature:
            loc = (data(i).G(trial).TR.tempIdx==I_max(trial,i)) & ...
                  (data(i).G(trial).TR.rateIdx==rateIDX);
            testData = y(loc);
            SE = std(testData,'omitnan');
            boundLimit = max_dist(trial,i)-SE; %since this is maximum value, the limit would only exist below
            tempLocWithinBounds = positionAvg>=boundLimit;
            
            % find other temps that have their distance within 1 STD of the max food distance
            upTemp = tempLocWithinBounds(I_max(trial,i):end);
            upLoc = find(diff(upTemp)<0);
            if ~isempty(upLoc)
                tempLocWithinBounds(I_max(trial,i)+upLoc(1):end) = false;
            end
            downTemp = tempLocWithinBounds(1:I_max(trial,i));
            downLoc = find(diff(downTemp)>0);
            if ~isempty(downLoc)
                tempLocWithinBounds(1:downLoc(end)) = false;
            end
            tempswithinBounds = tempList(tempLocWithinBounds);
    
            % Save the temp ranges for minumum distance:
            maxTempRangeStart(trial,i) = tempswithinBounds(1);
            maxTempRangeEnd(trial,i) = tempswithinBounds(end);
            maxTemp(trial,i) = tempList(I_max(trial,i));
        end
        grouped(i).(type).minDist.TempRange = [minimumTempRangeStart(:,i),minimumTempRangeEnd(:,i)];
        grouped(i).(type).minDist.Temp = mimimumTemp(:,i);
        grouped(i).(type).maxDist.TempRange = [maxTempRangeStart(:,i),maxTempRangeEnd(:,i)];
        grouped(i).(type).maxDist.Temp = maxTemp(:,i);
    end
end

% FIGURE: 
buff = 0.15;
r = 1;
c = 2;
yLimits = [14,24];
fig = figure; 
for i = 1:num.exp
    kolor = grouped(i).color;
    % MAX DISTANCE TEMP RANGE
    subplot(r,c,1); hold on
    % -- increasing ---
    y1 = mean(grouped(i).increasing.maxDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).increasing.maxDist.TempRange(:,2),'omitnan');
    plot([i+buff,i+buff],[y1,y2],'color',kolor,'LineWidth',1.5) 
    scatter(i+buff,mean(grouped(i).increasing.maxDist.Temp,'omitnan'),50,kolor,'filled')
    % -- decreasing ---
    y1 = mean(grouped(i).decreasing.maxDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).decreasing.maxDist.TempRange(:,2),'omitnan');
    plot([i-buff,i-buff],[y1,y2],'color',kolor,'LineWidth',1.5,'LineStyle',':') 
    scatter(i-buff,mean(grouped(i).decreasing.maxDist.Temp,'omitnan'),50,kolor,'filled')

    % MIN DISTANCE TEMP RANGE
    subplot(r,c,2); hold on
    % -- increasing ---
    y1 = mean(grouped(i).increasing.minDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).increasing.minDist.TempRange(:,2),'omitnan');
    plot([i+buff,i+buff],[y1,y2],'color',kolor,'LineWidth',1.5) 
    scatter(i+buff,mean(grouped(i).increasing.minDist.Temp,'omitnan'),50,kolor,'filled')
    % -- decreasing ---
    y1 = mean(grouped(i).decreasing.minDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).decreasing.minDist.TempRange(:,2),'omitnan');
    plot([i-buff,i-buff],[y1,y2],'color',kolor,'LineWidth',1.5,'linestyle',':') 
    scatter(i-buff,mean(grouped(i).decreasing.minDist.Temp,'omitnan'),50,kolor,'filled')
end
subplot(r,c,2)
xlim([0,num.exp+1]); ylim(yLimits)
set(gca,'xtick',0:num.exp+1,'XTickLabel',[])
xlabel(' ')
ylabel('temp (\circC) when closest to food')
subplot(r,c,1)
xlim([0,num.exp+1]); ylim(yLimits)
set(gca,'xtick',0:num.exp+1,'XTickLabel',[])
xlabel(' ')
ylabel('temp (\circC) when furthest from food')
formatFig(fig,true,[r,c]);
    
save_figure(fig,[saveDir expGroup ' min and max ranges'],'-png');


%% TODO ANALYSIS AND FIGURE: Closest neighbor analysis
clearvars('-except',initial_vars{:})

data(1).data(1).data.x_loc



%% TODO ANALYSIS AND FIGURE: ramp to ramp comparisons of movement
clearvars('-except',initial_vars{:})


%% TODO ANALYSIS AND FIGURE: 3D temperature modulation of behavior

clearvars('-except',initial_vars{:})

% hysteresis 

% distance range

% 

















%% FIGURE: TODO Comparison of temperature driven movement

% temp migration:  avg distance at hottest - avg distance at coolest
% distance at warmest
% distance at coolest
% dashed line between them?


fig = figure;





% FIGURE:
fig = figure; set(fig,'color','w',"Position",[1934 468 1061 590])

% Pull difference in distance heating-cooling
subplot(r,c,1)
hold on
for i = 1:num.exp
    x = repmat(grouped(i).decreasing.temps,[1,num.trial(i)]);
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    kolor = grouped(i).color;
    plot(x,y,'color',kolor,'LineWidth',LW); 
    plot(mean(x,2),mean(y,2),'color',kolor,'LineWidth',3.5)
end
h_line(0,'w',':',1)
xlabel('temp (\circC)')
ylabel('distance difference (mm)')


%% FIGURE: TODO COM for ramp up and ramp down, without holding period






































