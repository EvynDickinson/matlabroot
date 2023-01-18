
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

%% ANALYSIS: organize data for each group
disp(expNames')
clearvars('-except',initial_vars{:})
initial_vars = [initial_vars(:); 'initial_vars'; 'grouped'; 'expGroup'; 'saveDir'; 'mat';'expOrder'];
initial_vars = unique(initial_vars);

grouped = struct;

if strcmp(expGroup,'WT linear recovery caviar')
    expOrder = [];
    expList = {'Berlin WT','CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe'};
%     expList = {'Swedish', 'Berlin WT', 'OregonR','CantonS', 'Malawi','Zimbabwe'};
    colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White'};
    for ii = 1:num.exp
        expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar']));
%         grouped(expOrder(ii)).color = Color(colors{(ii)});
    end
elseif strcmp(expGroup,'All linear recovery ramp')
    expOrder = [];
    expList = {'Berlin WT','CantonS', 'OregonR', 'Swedish', 'Malawi', 'Zimbabwe','IsoD1','W1118'};
    colors = {'DarkOrchid','DeepSkyBlue','LimeGreen','Red','Gold','White','LightSalmon','Blue','DeepPink'};
    for ii = 1:num.exp
        expOrder(ii) = find(strcmp(expNames,[expList{ii} ' linear recovery ramp caviar']));
%         grouped(expOrder(ii)).color = Color(colors{(ii)});
    end
else
    expOrder = 1:num.exp;
    
%     colors = {'BlueViolet','red','white','Turquoise','Gold','pink','Orange'};

%     % Temperature shift experiments
%     expOrder = [3, 1, 2]; %lowest to highest
%     colors = {'dodgerblue','red','yellow','Gold','pink','Orange'};
    
%     colors = {'BlueViolet','white','turquoise','Gold','pink','Orange'}; % ===food no-food comp===
%     colors = {'white','BlueViolet','turquoise','Gold','pink','Orange'};
%     colors = {'white','LightSalmon','DeepPink','turquoise','Gold','pink','Orange'};
%     colors = {'Blue','lightskyblue','white','turquoise','pink','Orange'};
%     colors = {'red','yellow','dodgerblue','Gold','pink','Orange'};
%     colors = {'Teal','white','gold', 'magenta','dodgerblue','Orange'};
%===UAS linecomparisons===
%     colors = {'Grey','DeepSkyBlue','Blue', 'White','DeepPink'}; 
%     colors = {'BlueViolet','deeppink','orange', 'magenta','dodgerblue','Orange'};
% ===sex comparison colors===
    colors = {'White','magenta','dodgerblue','Orange'}; 
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
plot_err = true;
blkbgd = true;
% Y limit ranges
speed_lim = [0,6.2]; %speed
dist_lim = [10, 31]; %distance
% dt_lim = [14, 28];        %distance-temp
dt_lim = [10, 28];        %distance-temp
nMax = num.exp;%

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
for i = 1:nMax
%     i = expOrder(ii);
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
        plot(x,y,'LineWidth',LW,'Color',kolor)
        ylim(dist_lim)

%     %speed
    subplot(r,c,sb(3).idx); hold on
        y = smooth(grouped(i).speed.avg,'moving',sSpan);
        y_err = smooth(grouped(i).speed.err,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        ylim(speed_lim)

    %temp dependent distance
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
% temp-distance relationship 
subplot(r,c,sb(4).idx) 
ylabel('distance (mm)')
xlabel('temp (\circC)')
ylim(dt_lim)
h_line([18.1,36.2],'grey',':',1)
% 
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' timecourse summary'],'-png');

%% FIGURE & STATS: cumulative hysteresis for each genotype / trial

clearvars('-except',initial_vars{:})
LW = 0.75;
buff = 0.2;
SZ = 50;
r = 1; %rows
c = 3; %columns
plot_err = false;
plotSig = true; %plot significance stars


% FIGURE:
fig = figure; set(fig,'color','w',"Position",[1932 517 1050 611])
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

% STATS: are the means of any groups different from zero?
[p, mlt, id] = deal([]); 
for ii = 1:num.exp
    i = expOrder(ii);
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    plotY = sum(y,1,'omitnan');
    [~,p(ii)] = ttest(plotY); 
    group_name{ii} = expNames{i};
    %multicompare
    mlt = autoCat(mlt, plotY',false);
    id = autoCat(id,i*ones(length(plotY),1),false);
end
%Bonferonni correction: 
alpha = 0.05;
m = num.exp;
p_limit = alpha/m;
h = p<=p_limit;
stats_tbl = table(group_name',h',p','VariableNames',{'group','significant','p value'});
disp(stats_tbl)

% add significance stars to the figure:
if plotSig
    y_pos = rangeLine(fig,1);
    subplot(r,c,3); hold on
    for ii = 1:num.exp
        if h(ii)
            scatter(ii,y_pos,100,'w','*')
        end
    end
end

% Multicompare across the groups for significance
% STATS:

% determine which groups differ from each other
[~,~,stats] = anova1(mlt(:),id(:),'off');
if strcmp(getenv('COMPUTERNAME'),'ACADIA')
    [c,~,~,~] = multcompare(stats,.05,'off');
else
    [c,~,~,~] = multcompare(stats,[],'off');
end
% bonferonni multiple comparisons correction
alpha = 0.05; %significance level
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nPosition hysteresis cross group comparison statistics\n\n')
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
switch questdlg('Short or long display range?','','Short (15 min)','Long (50 min)','Cancel','Short (15 min)')
    case 'Short (15 min)'
        ylimits = [-8,4];
        dispROI = 15;
    case 'Long (50 min)'   
        ylimits = [-15,15];
        dispROI = 50;
end
      

% 
timeROI = 60; % how many minutes to look at behavior after each event
duration = ceil(timeROI*3*60);

sections = {'increasing','decreasing','holding'};
s_color = {'red','dodgerblue','white'};
temp_limits = [nan,nan];
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
        [temp,temperature] = deal([]);
        for rr = 1:nrr
            ROI = tp.(tpBin)(rr,1):tp.(tpBin)(rr,1)+duration;
            temp(:,:,rr) = grouped(i).dist.all(ROI,:);
            temperature(:,rr) = grouped(i).temp(ROI);
        end
        
        temp_norm = temp-mean(temp(1:10,:,:),'omitnan'); %normalize to zero distance
        temp_avg = mean(temp_norm,3);
        
        % add to the grouped data
        grouped(i).aligned.([sections{ss} '_avg']) = temp_avg;
        grouped(i).aligned.([sections{ss} '_norm']) = temp_norm;
        grouped(i).aligned.([sections{ss} '_all']) = temp;
        grouped(i).aligned.([sections{ss} '_SEM']) = std(temp_avg,0,2,'omitnan')/sqrt(num.trial(i));
        grouped(i).aligned.([sections{ss} '_MEAN']) = mean(temp_avg,2,'omitnan');
        grouped(i).aligned.([sections{ss} '_temperature']) = mean(temperature,2,'omitnan');
    end
    temp_limits(1) = min([temp_limits(1),tp.threshLow]);
    temp_limits(2) = max([temp_limits(2),tp.threshHigh]);
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
r = 5;
c = 3;
sb(1).idx = 1; sb(2).idx = 2; sb(3).idx = 3; % temperature ramps
sb(4).idx = 4:3:(r*c); sb(5).idx = 5:3:(r*c);  sb(6).idx = 6:3:(r*c);

% dispROI = 50;
duration = ceil(dispROI*3*60);
x = linspace(0,dispROI,duration+1);
LW = 1.5;
SEM_shading = false;
sSpan = 1;


fig = figure; set(fig,'pos',[1932 586 1050 542])
for ss = 1:length(sections) %
    % temp ramp
    subplot(r,c,sb(ss).idx); hold on
    for i = 1:num.exp
        y = grouped(i).aligned.([sections{ss} '_temperature'])(1:duration+1);
        plot(x,y,'color', grouped(i).color,'linewidth',LW)
    end
    ylabel('\circC')
    ylim(temp_limits) % TODO 1/4[tp.threshLow,tp.threshHigh]
    % event-aligned plot
    subplot(r,c,sb(ss+3).idx); hold on
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
%     title(sections{ss})
    ylim(ylimits)
end
formatFig(fig,true,[r,c],sb);
for ss = 1:length(sections) %
    subplot(r,c,sb(ss).idx)
    set(gca,'xcolor','k')
end

save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' ...
            num2str(sSpan) ' duration ' num2str(dispROI) ' min'],'-png',true);



% % FIGURE: CROSS EXPERIMENT COMPARISION (WITHIN HEATING,COOLING,HOLDING -- NO TEMP PLOTS)
% r = 1;
% c = 3;
% 
% dispROI = 50;
% duration = ceil(dispROI*3*60);
% x = linspace(0,dispROI,duration+1);
% LW = 1.5;
% SEM_shading = false;
% sSpan = 1;
% 
% 
% fig = figure; set(fig,'pos',[1932 690 1050 438])
% for ss = 1:length(sections) %
%     subplot(r,c,ss); hold on
%     for i = 1:num.exp
%         y = grouped(i).aligned.([sections{ss} '_MEAN'])(1:duration+1); 
%         if SEM_shading
%             y_err = grouped(i).aligned.([sections{ss} '_SEM'])(1:duration+1);
%             fill_data = error_fill(x, y, y_err);
%             h = fill(fill_data.X, fill_data.Y,grouped(i).color, 'EdgeColor','none','HandleVisibility','off');
%             set(h, 'facealpha', 0.4)
%         end
%         plot(x, smooth(y,'moving',sSpan),'color',grouped(i).color,'LineWidth',LW)
%     end
%     xlabel('time (min)')
%     ylabel('distance from food (mm)')
%     title(sections{ss})
%     ylim(ylimits)
% end
% formatFig(fig,true,[r,c]);

% save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' ...
%             num2str(sSpan) ' duration ' num2str(dispROI) ' min'],'-png',true);
    
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
    if ~strcmpi(getenv('COMPUTERNAME'),'EVYNPC')
    
        % set color bar information
        colorData = uint8(round(kolor.*255)); % convert color map to uint8
        colormap(colorData);
        c = colorbar('color',foreColor);
        set(c,'Ticks',[0,1],'TickLabels',[mat(i).position(trial).tempBins(1),mat(i).position(trial).tempBins(end)]);
        c.Label.String = 'Temperature (\circC)';
        c.Label.VerticalAlignment = "bottom";
    end
    title([data(i).ExpGroup],'color',foreColor)
    save_figure(fig,[saveDir expGroup grouped(i).name ' COM position'],'-png',true);

end

%% ANALYSIS AND FIGURES: COM of postions for each trial divided by heating / cooling
% periods

% vectors to fly 'mass' in arena at different temperatures from food

clearvars('-except',initial_vars{:})

types = {'hold', 'up', 'down'};

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
    temporary = [temperature,x_avg,y_avg]; 
    
    % get temp color distribution (1/4 degree incrememts)
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    cMapRange = tp.threshLow:0.25:tp.threshHigh;
    
    bin = 1*60*3;
    ntemps = length(cMapRange)-1;
    % color map information
    g1 = floor(ntemps/2);
    g2 = ntemps-g1;
    g1_cMap = Color('darkblue','grey',g1); %deepskyblue
    g2_cMap = Color('grey','red',g2);
    cMap = [g1_cMap;g2_cMap];
   
    for type = 1:length(types)
        [dsData,position] = deal([]);
        for ramp = 1:size(tp.(types{type}),1)
            smoothData = [];
            ROI = tp.(types{type})(ramp,1):tp.(types{type})(ramp,2); %time points for first type period (e.g. first hold etc)
            % smooth into 1 minute bins to select one point for each minute
            for c = 1:size(temporary,2)
                smoothData(:,c) = smooth(temporary(ROI,c),'moving',bin);
            end
            dsData = autoCat(dsData, smoothData(1:bin:end,:)); % select one point for each minute
        end
        dsData(:,4) = discretize(dsData(:,1),cMapRange); %binned temp category
        dsData(isnan(dsData(:,4)),:) = []; %remove data points outside allowed range
        
        %avg position for each binned temp
        for g = 1:ntemps
            loc = dsData(:,4)==g;
            position(g,:) = mean(dsData(loc,2:3),1,'omitnan');
        end
        % save data into structure:
        mat(i).position(trial).(types{type}).data = position;
        mat(i).position(trial).(types{type}).tempBins = cMapRange;
        mat(i).position(trial).(types{type}).cMap = cMap;
        mat(i).position(trial).(types{type}).wells = data(i).data(trial).data.wellcenters;
        mat(i).position(trial).(types{type}).wellLoc = well_loc;
    end
  end
end


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
types = {'hold','down','up'};

for i = 1:num.exp
  fig = figure; set(fig, 'pos',[41 282 1368 557])%[26 722 1077 568])
     for trial = 1:num.trial(i)
        wells = mat(i).position(trial).wells;
        wellLoc = mat(i).position(trial).wellLoc;
        kolor = mat(i).position(trial).cMap;
        
        for tt = 1:3 
            subplot(1,3,tt); hold on
            x = mat(i).position(trial).(types{tt}).data(:,1);
            y = mat(i).position(trial).(types{tt}).data(:,2);

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
     end
     % Formatting
     formatFig(fig,blkbgk,[1,3]);
     for tt = 1:3   
        subplot(1,3,tt)
        viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color','k');
        axis square; 
        axis equal;
        title(types{tt},'color','w')
        set(gca,'XColor',backColor,'YColor',backColor);
        if ~strcmpi(getenv('COMPUTERNAME'),'EVYNPC')
            % set color bar information
            colorData = uint8(round(kolor.*255)); % convert color map to uint8
            colormap(colorData);
            c = colorbar('color',foreColor);
            set(c,'Ticks',[0,1],'TickLabels',[mat(i).position(trial).tempBins(1),mat(i).position(trial).tempBins(end)]);
            c.Label.String = 'Temperature (\circC)';
            c.Label.VerticalAlignment = "bottom";
        end
     end
    
%     title([data(i).ExpGroup],'color',foreColor)
    save_figure(fig,[saveDir expGroup grouped(i).name ' rate divided COM position'],'-png',true);

end

%% FIGURE: Normalized increasing/decreasing temp responses
% TODO: add option to plot the heating | cooling on different subplots
clearvars('-except',initial_vars{:})

sep_h_c = true; %separate heating and cooling ramps
plot_err = false;
blkbgd = true;

% set up figure aligments
if sep_h_c 
    r = 5; %rows
    c = 11; %columns
    sb(1).idx = 1:6; %temp timecourse
    sb(2).idx = [12:(12+5),23:(23+5),34:(34+5),45:(45+5)]; %distance from food timecourse
    sb(3).idx = [8:c:r*c,9:c:r*c]; %binned distance alignment cooling
    sb(4).idx = [10:c:r*c,11:c:r*c]; %binned distance alignment heating
    sb(5).idx = 7:c:r*c;%buffer space
else
    r = 5; %rows
    c = 3; %columns
    sb(1).idx = 1:2; %temp timecourse
    sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse
    sb(3).idx = 3:3:15; %binned distance alignment
end

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
fig = figure; set(fig,'color','w',"Position",[1725 615 1140 590])
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',LW,'Color',kolor)

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

    %cooling dependent distance
    subplot(r,c,sb(3).idx); hold on
        
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
    
    if sep_h_c
%         subplot(r,c,sb(5).idx);
        %heating dependent distance
        subplot(r,c,sb(4).idx); hold on
    end
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
ylimits = ylim;
ylabel('distance (zscore)')
xlabel('temp (\circC)')
if sep_h_c
    subplot(r,c,sb(4).idx) 
    xlabel('temp (\circC)')
    set(gca,'ycolor',backColor)
    
    %set same y scale range
    ylimits_2 = ylim;
    ymax = max([ylimits(2),ylimits_2(2)]);
    ymin = min([ylimits(1),ylimits_2(1)]);
    ylim([ymin, ymax])
    subplot(r,c,sb(3).idx) 
    ylim([ymin, ymax])
end
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' zscore timecourse summary'],'-png');


% FIGURE: Separated heating & cooling 
fig = figure; set(fig,'color','w',"Position",[1725 501 1032 704])
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;
    %cooling dependent distance
    subplot(1,2,1); hold on
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
        plot(x,y,'LineWidth',LW,'Color',kolor,'HandleVisibility','off');
        title('Heating')
        xlabel('Temperature (\circC)')
        ylabel('Distance z-score')
    %heating dependent distance
    subplot(1,2,2); hold on
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
        plot(x,y,'LineWidth',LW,'Color',kolor)       
        title('Cooling')
        xlabel('Temperature (\circC)')
        ylabel('Distance z-score')
        % Names and Colors of included data
        dataString{i} = grouped(i).name;
end
% FORMATING AND LABELS
formatFig(fig,blkbgd,[1,2]);
fig = matchAxis(fig,true);
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 6)

save_figure(fig,[saveDir expGroup ' zscore heating and cooling overlay'],'-png');

%% FIGURE & STATS: Distance traveled over temperature variation
% TODO: fix error message
% Error using multcompare
% ALPHA must be a scalar between 0 and 1, exclusive.
% Error in QuadStep4_1_TempRate_TuningCurve (line 1273)

clearvars('-except',initial_vars{:})
autoAxis = true;
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
if ~autoAxis
    ylim([0,22])
end
% xlabel('Fly Strain')
ylabel('temp-driven spatial range (mm)')
formatFig(fig,true);
set(gca,'XColor',backColor)
save_figure(fig,[saveDir expGroup ' distance modulation by temperature'],'-png');

% STATS:
[datastats.all, datastats.id] = deal([]);
for i = 1:num.exp
    y = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
    datastats.all = autoCat(datastats.all,y,false);
    datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
end

% determine which groups differ from each other
[~,~,stats] = anova1(datastats.all,datastats.id,'off');
if strcmp(getenv('COMPUTERNAME'),'ACADIA')
    [c,~,~,~] = multcompare(stats);
else
    [c,~,~,~] = multcompare(stats,[],'off');
end
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

%% FIGURE & STATS: Min and max distance from the food (for a temperature bin)
clearvars('-except',initial_vars{:})
blkbgd = true;

y_limits = [5,35];
if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end
buff = 0.2;
r = 2; 
c = 2;
SZ = 30;
LW = 1.5;

% GATHER DATA
for tt = 1:2
  switch tt
     case 1
        type = 'increasing';
     case 2
        type = 'decreasing';
  end

    % Pull the avg min|max distance
    [min_dist,max_dist]  = deal(nan([max(num.trial),num.exp]));
    for i = 1:num.exp
        for trial = 1:num.trial(i)
            y = grouped(i).(type).all(:,trial);
            % Minimum distance
            min_dist(trial,i) = min(y); %I is the temp index 
            % Maximum distance
            max_dist(trial,i) = max(y);
        end
        grouped(i).(type).minDist.dist = min_dist(:,i);
        grouped(i).(type).maxDist.dist = max_dist(:,i);
    end 
end


fig = figure; set(fig,'pos',[103 558 631 718])
for ii = 1:num.exp
    i = expOrder(ii); %plot in correct order
    kolor = grouped(i).color;
    
    % Increasing temp rate -> minimum distance to food
    subplot(r,c,1); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).increasing.minDist.dist;
        y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)
        
    % Increasing temp rate -> maximum distance to food
    subplot(r,c,2); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).increasing.maxDist.dist;
         y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)
    % Decreasing temp rate -> minimum distance to food
    subplot(r,c,3); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).decreasing.minDist.dist;
         y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)
    % Decreasing temp rate -> maximum distance to food
    subplot(r,c,4); hold on
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        y = grouped(i).decreasing.maxDist.dist;
         y(isnan(y)) = [];
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[mean(y),mean(y)],'color',foreColor,'linewidth',LW)
end
%labels and formatting
formatFig(fig,true,[r,c]);
for ii = 1:4
    subplot(r,c,ii)
    set(gca,'xcolor', backColor)
    xlim([0,num.exp+1])
    ylim(y_limits)
    ylabel('food distance (mm)')
    switch ii
        case 1 
            title('Min | Heating','color',foreColor)
        case 2
            title('Max | Heating','color',foreColor)
        case 3
            title('Min | Cooling','color',foreColor)
        case 4
            title('Max | Cooling','color',foreColor)
    end
end
save_figure(fig,[saveDir expGroup ' food proximity modulation by temperature'],'-png');


% STATS: are the distances different from each other for each condition (min/max | heat/cool)?
for ll = 1:2 
    switch ll
        case 1 
            distGroup = 'minDist'; 
        case 2 
            distGroup = 'maxDist'; 
    end
    for tt = 1:2
        switch tt
         case 1
            type = 'increasing';
         case 2
            type = 'decreasing';
        end

        [datastats.all, datastats.id] = deal([]);
        for i = 1:num.exp
            y = grouped(i).(type).(distGroup).dist; y(isnan(y)) = [];
            datastats.all = autoCat(datastats.all,y',false);
            datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
        end
        % % determine which groups differ from each other
        [~,~,stats] = anova1(datastats.all,datastats.id,'off'); close all
        [c,~,~,~] = multcompare(stats,[],'off');

        % bonferonni multiple comparisons correction
        alpha = 0.05; %significance level
        m = size(c,1); %number of hypotheses
        sigThreshold = alpha/m;
        %find p-values that fall under the threshold
        significantHypotheses = c(:,6)<=sigThreshold;
        fprintf(['\n\n' type ' temp rate | ' distGroup ' comparison\n\n'])
        [Group1,Group2,P_Value] = deal([]);
        idx = 0;
        if any(significantHypotheses)
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
        else
            disp('no significance')
        end

    end
end

%% ANALYSIS AND FIGURE: temp range at minimum and max distances
% TODO: statistical comparisions within groups 1/2/23
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
% yLimits = [14,24];
yLimits = [14,28];
fig = figure; 
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    % MAX DISTANCE TEMP RANGE
    subplot(r,c,1); hold on
    % -- increasing ---
    y1 = mean(grouped(i).increasing.maxDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).increasing.maxDist.TempRange(:,2),'omitnan');
    plot([ii+buff,ii+buff],[y1,y2],'color',kolor,'LineWidth',1.5) 
    scatter(ii+buff,mean(grouped(i).increasing.maxDist.Temp,'omitnan'),50,kolor,'filled')
    % -- decreasing ---
    y1 = mean(grouped(i).decreasing.maxDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).decreasing.maxDist.TempRange(:,2),'omitnan');
    plot([ii-buff,ii-buff],[y1,y2],'color',kolor,'LineWidth',1.5,'LineStyle',':') 
    scatter(ii-buff,mean(grouped(i).decreasing.maxDist.Temp,'omitnan'),50,kolor,'filled')

    % MIN DISTANCE TEMP RANGE
    subplot(r,c,2); hold on
    % -- increasing ---
    y1 = mean(grouped(i).increasing.minDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).increasing.minDist.TempRange(:,2),'omitnan');
    plot([ii+buff,ii+buff],[y1,y2],'color',kolor,'LineWidth',1.5) 
    scatter(ii+buff,mean(grouped(i).increasing.minDist.Temp,'omitnan'),50,kolor,'filled')
    % -- decreasing ---
    y1 = mean(grouped(i).decreasing.minDist.TempRange(:,1),'omitnan');
    y2 = mean(grouped(i).decreasing.minDist.TempRange(:,2),'omitnan');
    plot([ii-buff,ii-buff],[y1,y2],'color',kolor,'LineWidth',1.5,'linestyle',':') 
    scatter(ii-buff,mean(grouped(i).decreasing.minDist.Temp,'omitnan'),50,kolor,'filled')
end
formatFig(fig,true,[r,c]);
subplot(r,c,2)
xlim([0,num.exp+1]); ylim(yLimits)
set(gca,'xcolor', 'k') % 'xtick',0:num.exp+1,'XTickLabel',[])
xlabel(' ')
ylabel('temp (\circC) when closest to food')
subplot(r,c,1)
xlim([0,num.exp+1]); ylim(yLimits)
set(gca,'xcolor', 'k')%  set(gca,'xtick',0:num.exp+1,'XTickLabel',[])
xlabel(' ')
ylabel('temp (\circC) when furthest from food')


% TODO STATS:  is there overlap between range? that should tell something??  

    
save_figure(fig,[saveDir expGroup ' min and max ranges'],'-png');

%% !SLOW! ANALYSIS AND FIGURE: Nearest neighbor analysis
% WARNING: TAKES A MINUTE OR TWO (COMPUTATION HEAVY)
clearvars('-except',initial_vars{:})
pix2mm = 12.8;

% Calculate nearest neighbor distance for each frame
if ~isfield(grouped,'NN')
    for i = 1:num.exp
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
            grouped(i).NN.data(trial).all = NN./pix2mm; %get mm from pixel distance
            disp(['Done exp ' num2str(i) ' trial ' num2str(trial)])
        end  
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
end

% ==== TIME COURSE FIGURE ====
blkbgd = true;
% set up figure aligments
r = 4; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8,10,11]; % nearest neighbor 
sb(3).idx = [3,6,9,12]; % nearest neighbor distance by temperature
LW = 0.75;
sSpan = 180;
if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end

% ======== TIMECOURSE FIGURE =========
fig = figure; set(fig,'pos',[2111 307 1258 732]); hold on
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

formatFig(fig, true,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca,'xcolor', backColor)

% save figure
save_figure(fig,[saveDir expGroup ' nearest neighbor timecourse'],'-png');


% ======== DISTANCE TO FOOD VS N.N. DISTANCE FIGURE =========
LW = 1;
fig = figure; set(fig,'pos',[2111 307 497 732]); hold on

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
formatFig(fig, true);

save_figure(fig,[saveDir expGroup ' nearest neighbor vs food distance'],'-png');

%% FIGURE: 3D temperature modulation of behavior
clearvars('-except',initial_vars{:})

%for each trial -- plot it in 3D space based on: 
% X) distance range
% Y) minimum distance to food (either during heating/cooling)
% Z) hysteresis 

SZ = 75;

fig = figure; set(fig,'pos',[-1054 648 1016 592]);
for i = 1:num.exp
    kolor = grouped(i).color;
    x = range([grouped(i).increasing.all;grouped(i).decreasing.all],1);
    y = min([grouped(i).increasing.minDist.dist';grouped(i).decreasing.minDist.dist']);
    y(isnan(y)) = [];
    z = sum(grouped(i).decreasing.all-grouped(i).increasing.all,'omitnan');
    
    scatter3(x,y,z,SZ,kolor,'filled') 
    hold on
end
% Labels and formatting
formatFig(fig,true);    
xlabel('range (mm)')
ylabel('min distance (mm)')
zlabel('hysteresis (mm)')
set(gca,'zcolor','w')

save_figure(fig,[saveDir expGroup ' 3D space modulation by temperature'],'-png');

%% FIGURE: speed histogram
clearvars('-except',initial_vars{:})

% Simple histogram of speed across full experiment
bins = 0:0.5:15;
% xlimit = [-0.5,7];
xlimit = [-.5,20];

[nrows, ncols] = subplot_numbers(num.exp, 3);
if num.exp==3
    nrows = 1;
    ncols = 3;
end
if nrows==1
    figPos = [-979 638 953 439];
else
    figPos = [132 133 953 1021];
end


fig = figure; set(fig,'pos',figPos)% );
for ii = 1:num.exp
    i = expOrder(ii);
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,true);
    end
    % histogram for each group
    subplot(nrows,ncols,ii)
    histogram(temp,bins,'FaceColor',grouped(i).color,'facealpha',1)
    v_line(mean(temp,'omitnan'),'w','--',1.5)
%     title(grouped(i).name)
    xlabel('speed (mm/s)')
    ylabel('count')
    xlim(xlimit)
end
formatFig(fig,true,[nrows,ncols]);
save_figure(fig,[saveDir expGroup ' speed histograms'],'-png',true);

%% FIGURE AND STATS: average experiment speed
% Comparision of avg speed and range across full experiments
clearvars('-except',initial_vars{:})
plot_err = true;
buff = 0.2;
SZ = 50;
LW = 1.5;
r = 1;
c = 3;
sb(1).idx = 1:2; %avg speed per temp
sb(2).idx = 3;   %avg speed


fig = figure; set(fig,'pos',[2108 475 1000 590]);
subplot(r,c,sb(1).idx); hold on
    %pull temp avg data
    [speed,speed_err,temp] = deal([]);
    for i = 1:num.exp
        temp(:,i) = data(i).G(1).TR.temps;
        tempIdx = data(i).G(1).TR.tempIdx;
        for t = 1:size(temp,1)-1
            loc = (tempIdx==t);
            raw = mean(grouped(i).speed.all(loc,:),'omitnan');
            y = mean(raw);
            y_err = std(raw)./sqrt(num.trial(i));
            speed(t,i) = y;
            speed_err(t,i) = y_err;
        end
    end
    for i = 1:num.exp
        kolor = grouped(i).color;
        loc = ~isnan(speed(:,i));
        x = temp(loc,i);
        y = speed(loc,i);
        y_err = speed_err(loc,i);
        if plot_err
           fill_data = error_fill(x, y, y_err);
          h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
          set(h, 'facealpha', 0.3)
        end
        plot(x,y,'color',kolor,'linewidth', LW)
    end
    xlabel('Temperature (\circC)')
    ylabel('Speed (mm/s)')

subplot(r,c,sb(2).idx);hold on
for ii = 1:num.exp
    i = expOrder(ii);
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,false);
    end
    y = mean(temp,1,'omitnan');
    y_avg = mean(mean(temp,'omitnan'));
    x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
    scatter(x,y,SZ,grouped(i).color,'filled')
    plot([ii-buff,ii+buff],[y_avg,y_avg],'color', grouped(i).color,'linewidth',LW)
end
xlim([0,num.exp+(buff*2)])
ylabel('avg speed (mm/s)')
formatFig(fig,true,[r,c],sb);
set(gca,'xcolor','k')
save_figure(fig,[saveDir expGroup ' avg speed'],'-png',true);

% AVG SPEED STATS:
[datastats.all, datastats.id] = deal([]);
for i = 1:num.exp
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,false);
    end
    y = mean(temp,1,'omitnan');

    datastats.all = autoCat(datastats.all,y,false);
    datastats.id = autoCat(datastats.id,repmat(i,[1,num.trial(i)]),false);
end

% determine which groups differ from each other
[~,~,stats] = anova1(datastats.all,datastats.id,'off');
[c,~,~,~] = multcompare(stats,[],'off');

% bonferonni multiple comparisons correction
alpha = 0.05; %significance level
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nAverage speed statistics\n\n')
[Group1,Group2,P_Value] = deal([]);
idx = 0;
if ~any(significantHypotheses)
    disp('No statistical differences in avg speed between groups')
else
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
end

%% FIGURE: cumulative distribution function of speed
clearvars('-except',initial_vars{:})
LW = 1;
xlimit = [0,10];
% xlimit = [0,20];

fig = figure; hold on
for i = 1:num.exp
    temp = [];
    for trial = 1:num.trial(i)
        temp = autoCat(temp,data(i).data(trial).speed.avg,true);
    end
    h = cdfplot(temp);
    set(h,'color',grouped(i).color,'linewidth',LW)
end
xlim(xlimit)
xlabel('Speed (mm/s)')
ylabel('Empirical CDF')
title('')
formatFig(fig, true);   
set(gca,'xgrid','off','ygrid','off')
set(gca,'ytick',0:0.2:1)
save_figure(fig,[saveDir expGroup ' speed CDF'],'-png',true);
    
%% FIGURE & STATS: speed hysteresis across groups
clearvars('-except',initial_vars{:})
plotSig = true;
LW = 0.75;
buff = 0.2;
SZ = 50;
r = 1; %rows
c = 3; %columns

[groupName, stats] = deal([]);

fig = figure; set(fig,'color','w','pos',[1932 481 1268 647])%[1932 517 1050 611]
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    x = data(i).G(1).TR.temps(1:end-1);
    [h_speed,c_speed] = deal([]);
    rateIdx = data(i).G(1).TR.rateIdx;
    
    % heating locations 
    idx = find(data(i).G(1).TR.rates>0);
    heatloc = rateIdx==idx;
    h_style = '-'; 
    % cooling locations
    idx = find(data(i).G(1).TR.rates<0);
    coolloc = rateIdx==idx;
    c_style = '--';

    %pull temp avg data
    tempIdx = data(i).G(1).TR.tempIdx;
    speed_data = grouped(i).speed.all;
    for t = 1:length(x)
        temp_loc = (tempIdx==t);
        % heating
        loc = temp_loc & heatloc;
        raw = mean(speed_data(loc,:),'omitnan');
        h_speed(t,:) =  raw;
        % cooling
        loc = temp_loc & coolloc;
        raw = mean(speed_data(loc,:),'omitnan');
        c_speed(t,:) =  raw;
    end
    
    % plot heat and cooling lines
    subplot(r,c,1); hold on
        heat_y = mean(h_speed,2,'omitnan'); 
        cool_y = mean(c_speed,2,'omitnan'); 
        plot(x,heat_y,'color',kolor,'linewidth',LW,'linestyle',h_style)
        plot(x,cool_y,'color',kolor,'linewidth',LW,'linestyle',c_style)

    % plot speed hysteresis for each temp bin
    subplot(r,c,2); hold on
        y = cool_y-heat_y;
        plot(x,y,'color',kolor,'linewidth',LW)
    
    % plot cumulative hysteresis for each trial
    subplot(r,c,3); hold on
        y = sum((c_speed-h_speed),1,'omitnan');
        y_avg = mean(y);
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i))); 
        scatter(x,y,SZ,kolor,'filled')
        plot([ii-buff,ii+buff],[y_avg,y_avg],'color',kolor,'linewidth',1.3)
        % save cumulative difference for stats
        stats = autoCat(stats,y,true);
        groupName{ii} = grouped(i).name;
end
    
%formatting and labels
formatFig(fig,true,[r,c]);    
subplot(r,c,1)
    xlabel('Temp (\circC)')
    ylabel('Speed (mm/s)')
subplot(r,c,2)
    xlabel('Temp (\circC)')
    ylabel('Speed diff (cool-heat)')
    h_line(0,'w',':',1)
subplot(r,c,3)
    xlim([0.5,num.exp+.5])
    h_line(0,'w',':',1)
    ylabel('Cumulative speed difference')
    set(gca,'xcolor','k')
    
   
% SIG HYSTERESIS STATS: are the means of any groups different from zero?
p = [];
for ii = 1:num.exp
    [~,p(ii)] = ttest(stats(ii,:)); 
end
%Bonferonni correction: 
alpha = 0.05;
m = num.exp;
p_limit = alpha/m;
h = p<=p_limit;
stats_tbl = table(groupName',h',p','VariableNames',{'group','significant','p value'});
disp(stats_tbl)

% add significance stars to the figure:
if plotSig
    y_pos = rangeLine(fig,10,false);
    subplot(r,c,3); hold on
    for ii = 1:num.exp
        if h(ii)
            scatter(ii,y_pos,100,'w','*')
        end
    end
end   
    
% SIG DIFF BETWEEN GROUPS STATS:
%format for anova
stats_id = repmat((1:num.exp)',[1,size(stats,2)]);
stats_id = reshape(stats_id,[numel(stats),1]);
x = reshape(stats,[numel(stats),1]);
loc = isnan(x);
stats_id(loc) = [];
x(loc) = [];

% determine which groups differ from each other
[~,~,data_stats] = anova1(x,stats_id,'off');
[c,~,~,~] = multcompare(data_stats,[],'off');

% bonferonni multiple comparisons correction
alpha = 0.05; %significance level
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nSpeed hysteresis comparison statistics:\n\n')
[Group1,Group2,P_Value] = deal([]);
idx = 0;
if any(significantHypotheses)
    for i = 1:length(significantHypotheses)
        if significantHypotheses(i)
            idx = idx+1;
            Group1{idx,1} = groupName{c(i,1)};
            Group2{idx,1} = groupName{c(i,2)};
            P_Value(idx,1) = c(i,6);
        end
    end
    sig_comp = table(Group1,Group2,P_Value);
    disp(sig_comp)
else
    disp('no significant differences between groups')
end


% save figure
save_figure(fig,[saveDir expGroup ' speed hysteresis summary'],'-png');  

%% FIGURE: difference between a 'control' and 'test' group (i.e. food | no food)

clearvars('-except',initial_vars{:})

if ~(num.exp==2)
    warndlg('Need to input paired exp numbers')
    return
end
plot_err = true;
blkbgd = true;


LW = 0.75;
sSpan = 180;

if blkbgd
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end
compareY = [];
% FIGURE:
fig = figure; set(fig,'color','w',"Position",[2414 500 397 590])
hold on
for i = 1:num.exp
    x = grouped(i).time;
    kolor = grouped(i).color;
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
        compareY(:,i) = y;
end
diffY = compareY(:,2)-compareY(:,1);
plot(x, diffY,'color',Color('lightpink'),'linewidth',2);

% FORMATING AND LABELS
formatFig(fig,blkbgd);

% temp-distance relationship 
ylabel('distance (mm)')
xlabel('temp (\circC)')
ylim([0,35])
% 
legend(dataString,'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' distance difference between groups'],'-png');

%% FIGURE: Speed - distance correlation
clearvars('-except',initial_vars{:})
plot_err = true;
corr_coef = [];
num_points = 300;
buff = 0.2;
r = 1;
c = 3;
sb(1).idx = 1:2;
sb(2).idx = 3;

% get correlation data
for i = 1:num.exp
    % get speed / distance information
    speed = grouped(i).speed.all;   %speed for each trial within the exp.
    fly_dist = grouped(i).dist.all; %dist for each trial within the exp

    for trial = 1:num.trial(i)
        temp = [];
        temp(:,1) = speed(:,trial);
        temp(:,2) = fly_dist(:,trial);
        loc = any(isnan(temp),2);
        temp(loc,:) = [];
        % speed-distance correlation
        rho = corr(temp); 
        corr_coef(i).all(trial) = rho(1,2);
    end
 end
     
% PLOT FIGURE
fig = figure; set(fig,'pos',[201 528 1070 618])
% linear fit of speed-distance to food|source
subplot(r,c,sb(1).idx); hold on
    for i = 1:num.exp

        % get speed / distance information
        speed = grouped(i).speed.all;   %speed for each trial within the exp.
        fly_dist = grouped(i).dist.all; %dist for each trial within the exp
        % reshape & remove NaNs
        y = [];
        y(:,1) = speed(:);
        y(:,2) = fly_dist(:);
        loc = any(isnan(y),2);
        y(loc,:) = [];

        % speed-distance correlation
        rho = corr(y); 
        corr_coef(i).group = rho(1,2);

        % sort by speed
        sorted_data = sortrows(y,1);
        idx = randi(length(y),[num_points,1]);
        X = sorted_data(idx,1);
        Y = sorted_data(idx,2);

        % Plot random selection of points & linear regression line of best fit
        kolor = grouped(i).color;
        scatter(X,Y,10,kolor)
        fitX = sorted_data(:,1);
        [p,S] = polyfit(fitX,sorted_data(:,2),1); 
        [y_fit,delta] = polyval(p,fitX,S);
        plot(fitX,y_fit,'color',kolor,'linestyle','-','linewidth',1)
        if plot_err
            plot(fitX,y_fit+2*delta,'color',kolor,'linestyle','--','linewidth',1)
            plot(fitX,y_fit-2*delta,'color',kolor,'linestyle','--','linewidth',1)
        end
    end
    xlim([0,15])
    xlabel('Speed (mm/s)')
    ylabel('Distance from source (mm)')
    
% correlation coefficients
subplot(r,c,sb(2).idx); hold on
     for ii = 1:num.exp
       i = expOrder(ii);
       kolor = grouped(i).color;
       xlow = ii-buff-0.1;
       xhigh = ii+buff+0.1;
       x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
       y = corr_coef(i).all;
       y_avg = mean(corr_coef(i).all);
       scatter(x,y,50,kolor,'filled')
       plot([xlow,xhigh],[corr_coef(i).group,corr_coef(i).group],'color',kolor,'linestyle',':','linewidth',1.5)
       plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',1.5)
     end
     xlim([0.5,num.exp+.5])       
     ylabel('correlation coefficient')
     h_line(0,'w',':',1)
formatFig(fig,true,[r,c],sb);
subplot(r,c,sb(2).idx)
set(gca,'xcolor','k')

% save figure
save_figure(fig,[saveDir expGroup ' speed distance correlation'],'-png');  







% 
% % Save correlation coefficient?
% switch questdlg('Save correlation coefficient?')
%     case 'Yes'
% %         save()
% end


%    
% % SIG HYSTERESIS STATS: are the means of any groups different from zero?
% p = [];
% for ii = 1:num.exp
%     [~,p(ii)] = ttest(stats(ii,:)); 
% end
% %Bonferonni correction: 
% alpha = 0.05;
% m = num.exp;
% p_limit = alpha/m;
% h = p<=p_limit;
% stats_tbl = table(groupName',h',p','VariableNames',{'group','significant','p value'});
% disp(stats_tbl)
% 
% % add significance stars to the figure:
% if plotSig
%     y_pos = rangeLine(fig,10,false);
%     subplot(r,c,3); hold on
%     for ii = 1:num.exp
%         if h(ii)
%             scatter(ii,y_pos,100,'w','*')
%         end
%     end
% end   
%     
% % SIG DIFF BETWEEN GROUPS STATS:
% %format for anova
% stats_id = repmat((1:num.exp)',[1,size(stats,2)]);
% stats_id = reshape(stats_id,[numel(stats),1]);
% x = reshape(stats,[numel(stats),1]);
% loc = isnan(x);
% stats_id(loc) = [];
% x(loc) = [];
% 
% % determine which groups differ from each other
% [~,~,data_stats] = anova1(x,stats_id,'off');
% [c,~,~,~] = multcompare(data_stats,[],'off');
% 
% % bonferonni multiple comparisons correction
% alpha = 0.05; %significance level
% m = size(c,1); %number of hypotheses
% sigThreshold = alpha/m;
% %find p-values that fall under the threshold
% significantHypotheses = c(:,6)<=sigThreshold;
% fprintf('\n\nSpeed hysteresis comparison statistics:\n\n')
% [Group1,Group2,P_Value] = deal([]);
% idx = 0;
% if any(significantHypotheses)
%     for i = 1:length(significantHypotheses)
%         if significantHypotheses(i)
%             idx = idx+1;
%             Group1{idx,1} = groupName{c(i,1)};
%             Group2{idx,1} = groupName{c(i,2)};
%             P_Value(idx,1) = c(i,6);
%         end
%     end
%     sig_comp = table(Group1,Group2,P_Value);
%     disp(sig_comp)
% else
%     disp('no significant differences between groups')
% end
% 

%% FIGURE: Temperature - distance correlation
clearvars('-except',initial_vars{:})
corr_coef = [];
buff = 0.2;

% get correlation data
for i = 1:num.exp
    pooledData = [];
    % get speed / distance information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp;
        y = grouped(i).dist.all(:,trial);
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
fig = figure; set(fig,'color','w',"Position",[2108 475 453 590]); hold on
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
   plot([xlow,xhigh],[corr_coef(i).group,corr_coef(i).group],'color',kolor,'linestyle',':','linewidth',1.5)
   plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',1.5)
 end
 xlim([0.5,num.exp+.5])       
 ylabel('temp-distance corr. coef.')
 h_line(0,'w',':',1)    
 formatFig(fig,true);    
 set(gca,'xcolor','k')
 
% save figure
save_figure(fig,[saveDir expGroup ' temp distance correlation'],'-png');  

%% TODO ANALYSIS AND FIGURE: ramp to ramp comparisons of movement
% do we need a fourth ramp??
clearvars('-except',initial_vars{:})


% range traveled over ramp

% hysteresis over ramp

% min distance over ramp 

%% FIGURE: avg distance between wells using the pix2mm conversion
pix2mm = 12.8;

d = [];
for i = 1:num.exp 
  for trial = 1:num.trial(i)
    % distance between well 1 and 3
    well_1 = data(i).data(trial).data.wellcenters(:,1);
    well_3 = data(i).data(trial).data.wellcenters(:,3);
    d = [d;sum((well_1-well_3).^2).^0.5];
    % distance between well 2 and 4
    well_1 = data(i).data(trial).data.wellcenters(:,2);
    well_3 = data(i).data(trial).data.wellcenters(:,4);
    d = [d;sum((well_1-well_3).^2).^0.5];
  end
end
    
well_distance = d./pix2mm;
% mean(well_distance)

fig = figure; set(fig, 'pos',[2319 619 390 573]);
violin(well_distance,'xlabel',{''},'facecolor',...
       Color('teal'),'edgecolor','w','alp',1,...
       'bw',0.35,'mc','w','medc','r--')
ylabel('well distance (mm)','FontSize',14)
formatFig(fig,true);
h_line(36.2,'gold','--',1.5)
ylim([35,40])
set(gca,'xcolor','k','ytick',35:1:40)
save_figure(fig,[saveDir expGroup ' physcial well distance selection'],'-png'); 

%% FIGURE & STATS: Hysteresis within a specific time range
clearvars('-except',initial_vars{:})
LW = 0.75;
buff = 0.2;
SZ = 50;
r = 1; %rows
c = 3; %columns
plot_err = false;
plotSig = true; %plot significance stars

ROI = [14,18];

% FIGURE:
fig = figure; set(fig,'color','w',"Position",[1932 517 1050 611])
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
xlim([ROI(1),ROI(2)])

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
xlim([ROI(1),ROI(2)])
xlabel('temp (\circC)')
ylabel('distance difference (mm)')

% Cumulative difference in proximity
subplot(r,c,3)
hold on

%TODO: update for smaller range or range that doesn't match exactly the
%temperature bins
for ii = 1:num.exp
    i = expOrder(ii);
    kolor = grouped(i).color;
    %find temp ROI region: 
    loc = [];
    loc(1) = find(grouped(i).decreasing.temps==ROI(1));
    loc(2) = find(grouped(i).decreasing.temps==ROI(2));
    y1 = grouped(i).decreasing.all(loc(1):loc(2),:);
    y2 = grouped(i).increasing.all(loc(1):loc(2),:);
    y = y1-y2;
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

% STATS: are the means of any groups different from zero?
[p, mlt, id] = deal([]); 
for ii = 1:num.exp
    i = expOrder(ii);
    loc = [];
    loc(1) = find(grouped(i).decreasing.temps==ROI(1));
    loc(2) = find(grouped(i).decreasing.temps==ROI(2));
    y1 = grouped(i).decreasing.all(loc(1):loc(2),:);
    y2 = grouped(i).increasing.all(loc(1):loc(2),:);
    y = y1-y2;
    plotY = sum(y,1,'omitnan');
    [~,p(ii)] = ttest(plotY); 
    group_name{ii} = expNames{i};
    %multicompare
    mlt = autoCat(mlt, plotY',false);
    id = autoCat(id,i*ones(length(plotY),1),false);
end
%Bonferonni correction: 
alpha = 0.05;
m = num.exp;
p_limit = alpha/m;
h = p<=p_limit;
stats_tbl = table(group_name',h',p','VariableNames',{'group','significant','p value'});
disp(stats_tbl)

% add significance stars to the figure:
if plotSig
    y_pos = rangeLine(fig,1);
    subplot(r,c,3); hold on
    for ii = 1:num.exp
        if h(ii)
            scatter(ii,y_pos,100,'w','*')
        end
    end
end


% STATS:
% determine which groups differ from each other
[~,~,stats] = anova1(mlt(:),id(:),'off');
if strcmp(getenv('COMPUTERNAME'),'ACADIA')
    [c,~,~,~] = multcompare(stats,.05,'off');
else
    [c,~,~,~] = multcompare(stats,[],'off');
end
% bonferonni multiple comparisons correction
alpha = 0.05; %significance level
m = size(c,1); %number of hypotheses
sigThreshold = alpha/m;
%find p-values that fall under the threshold
significantHypotheses = c(:,6)<=sigThreshold;
fprintf('\n\nPosition hysteresis cross group comparison statistics\n\n')
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

% save figure
save_figure(fig,[saveDir 'Hysteresis summary ' num2str(ROI(1)) ' to ' num2str(ROI(2)) ' deg'],'-png');

%% FIGURE: Event aligned with preceding data before event

clearvars('-except',initial_vars{:})
switch questdlg('Short or long display range?','','Short (15 min)','Long (50 min)','Cancel','Short (15 min)')
    case 'Short (15 min)'
%         ylimits = [-8,4];
        ylimits = [-5,2];
        dispROI = 10;
        pre_disp = 10;
    case 'Long (50 min)'   
        ylimits = [-15,15];
        dispROI = 50;
        pre_disp = 20;
end
        

% SELECT THE TIME ROIs
timeROI = 60; % how many minutes to look at behavior after each event
preROI = pre_disp;  % how many minutes to look at behavior before each event
post_duration = ceil(timeROI*3*60);
pre_duration = ceil(preROI*3*60);
total_duration = pre_duration+post_duration;

sections = {'increasing','decreasing','holding'};
s_color = {'red','dodgerblue','white'};
temp_limits = [nan,nan];
for i = 1:num.exp
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    for ss = 1:length(sections)
        switch sections{ss}
            case 'increasing'
                tpBin = 'up';
                nrr = tp.nUp; %number of events for this category
            case 'decreasing'
                tpBin = 'down';
                nrr = tp.nDown;
            case 'holding'
                tpBin = 'hold';
                nrr = tp.nHold;
        end
        [temp,temperature] = deal([]);
        for rr = 1:nrr
            ROI = tp.(tpBin)(rr,1)-pre_duration:tp.(tpBin)(rr,1)+post_duration;
            if any(ROI<=0)
                loc = (ROI<=0);
                ROI(loc) = [];
                xx = find(loc); 
                MT = nan(xx(end),num.trial(i));
                temp(:,:,rr) = [MT ; grouped(i).dist.all(ROI,:)];
                temperature(:,rr) = [MT(:,1); grouped(i).temp(ROI)];
            else
                temp(:,:,rr) = grouped(i).dist.all(ROI,:);
                temperature(:,rr) = grouped(i).temp(ROI);
            end
            
        end
        
        temp_norm = temp-mean(temp(pre_duration-5:pre_duration+5,:,:),'omitnan'); %normalize to zero distance
        temp_avg = mean(temp_norm,3,'omitnan');
        
        % add to the grouped data
        grouped(i).ext_aligned.([sections{ss} '_avg']) = temp_avg;
        grouped(i).ext_aligned.([sections{ss} '_norm']) = temp_norm;
        grouped(i).ext_aligned.([sections{ss} '_all']) = temp;
        grouped(i).ext_aligned.([sections{ss} '_SEM']) = std(temp_avg,0,2,'omitnan')/sqrt(num.trial(i));
        grouped(i).ext_aligned.([sections{ss} '_MEAN']) = mean(temp_avg,2,'omitnan');
        grouped(i).ext_aligned.([sections{ss} '_temperature']) = mean(temperature,2,'omitnan');
    end
    temp_limits(1) = min([temp_limits(1),tp.threshLow]);
    temp_limits(2) = max([temp_limits(2),tp.threshHigh]);
    grouped(i).ext_aligned.postTime = timeROI;
    grouped(i).ext_aligned.preTime = preROI ;
    grouped(i).ext_aligned.post_duration = post_duration;
    grouped(i).ext_aligned.pre_duration = pre_duration;  
end


% FIGURE: CROSS EXPERIMENT COMPARISION (WITHIN HEATING,COOLING,HOLDING)
r = 5;
c = 3;
sb(1).idx = 1; sb(2).idx = 2; sb(3).idx = 3; % temperature ramps
sb(4).idx = 4:3:(r*c); sb(5).idx = 5:3:(r*c);  sb(6).idx = 6:3:(r*c);


post_dur = ceil(dispROI*3*60);
pre_dur = ceil(pre_disp*3*60);
x = [linspace(-pre_disp,0,pre_dur), linspace(0,dispROI,post_dur+1)];
LW = 1.5;
SEM_shading = false;
sSpan = 30;
xlimit = [-pre_disp,dispROI];

fig = figure; set(fig,'pos',[1932 586 1050 542])
for ss = 1:length(sections) %
    % temp ramp
    subplot(r,c,sb(ss).idx); hold on
    for i = 1:num.exp
        y = grouped(i).ext_aligned.([sections{ss} '_temperature'])(1:length(x));
        plot(x,y,'color', grouped(i).color,'linewidth',LW)
    end
    ylabel('\circC')
    ylim(temp_limits) % TODO 1/4[tp.threshLow,tp.threshHigh]
    v_line(0,'w',':',1)
    xlim(xlimit)
    % event-aligned plot
    subplot(r,c,sb(ss+3).idx); hold on
    for i = 1:num.exp
        y = grouped(i).ext_aligned.([sections{ss} '_MEAN'])(1:length(x)); 
        if SEM_shading
            y_err = grouped(i).ext_aligned.([sections{ss} '_SEM'])(1:length(x));
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y,grouped(i).color, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.4)
        end
        plot(x, smooth(y,sSpan,'moving'),'color',grouped(i).color,'LineWidth',LW)
    end
    xlabel('time (min)')
    ylabel('distance (mm)')
%     title(sections{ss})
    ylim(ylimits)
    v_line(0,'w',':',1)
    xlim(xlimit)
end
formatFig(fig,true,[r,c],sb);
for ss = 1:length(sections) %
    subplot(r,c,sb(ss).idx)
    set(gca,'xcolor','k')
end

save_figure(fig,[saveDir expGroup ' event aligned distance - smoothed ' ...
            num2str(sSpan) ' duration ' num2str(dispROI) ' min with pretime'],'-png',true);



% % 
    
%% FIGURE: (WT) latutide organized 
clearvars('-except',initial_vars{:})
corr_coef = [];
buff = 0.2;

lat_list = {'Swedish', 'Berlin', 'Oregon','Canton','Malawi', 'Zimbabwe'};
latitudes = [60.1282,  52.5200,    43.8041,   40.7989, -13.2543, -19.0154];

% get correlation data
for i = 1:num.exp
    pooledData = [];
    % get speed / distance information
    for trial = 1:num.trial(i)
        x = data(i).data(trial).occupancy.temp;
        y = grouped(i).dist.all(:,trial);
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
fig = figure; set(fig,'color','w',"Position",[2108 475 453 590]); hold on
hold on
 for ii = 1:num.exp
   i = find(contains(expNames,lat_list(ii)));
     
   kolor = grouped(i).color;
   lat = latitudes(ii);
   xlow = lat-2;
   xhigh = lat+2;
%    x = shuffle_data(linspace(lat-buff,lat+buff,num.trial(i)));
   x = lat*ones(1,num.trial(i));
   y = corr_coef(i).all;
   y_avg = mean(corr_coef(i).all);
   scatter(x,y,50,kolor,'filled')
%    plot([xlow,xhigh],[corr_coef(i).group,corr_coef(i).group],'color',kolor,'linestyle',':','linewidth',1.5)
%    plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',1.5)
 end
 xlim([-30,70])       
 ylabel('temperature-distance correlation')
 h_line(0,'w',':',1)    
 formatFig(fig,true);    
 set(gca,'xcolor','w')
 xlabel('latitude (\circ)')
% save figure
save_figure(fig,[saveDir expGroup ' temp distance correlation by latitude 2'],'-png');  


    





























