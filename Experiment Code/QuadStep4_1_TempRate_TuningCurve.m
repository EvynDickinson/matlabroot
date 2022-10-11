
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

% Load selected experiment data groups:
for i = 1:num.exp
    data(i) = load([structFolder expNames{i} '\' expNames{i} ' post 3.1 data.mat']);
end

clear list_dirs expIdx dirIdx
% Set up base variables
initial_vars = who;
initial_vars = [initial_vars(:); 'initial_vars'; 'grouped'; 'expGroup'; 'saveDir'; 'mat'];
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

% colors for each group
colors = {'BlueViolet',...         'Orange',...
         'red',...
         'white',...
         'Turquoise',...
         'Gold',...
         'pink'};

grouped = struct;
for i = 1:num.exp % FOR EACH DATA GROUP
    % GENERAL
    grouped(i).name = data(i).ExpGroup;
    grouped(i).color = Color(colors{i});

    % TIME COURSE DATA
    num.trial(i) = data(i).ntrials;
    [time,temp,speed,distance] = deal([]);
    for trial = 1:num.trial(i)
        time = autoCat(time,data(i).data(trial).occupancy.time,false,true);
        temp = autoCat(temp,data(i).data(trial).occupancy.temp,false,true);
        speed = autoCat(speed,data(i).data(trial).speed.avg,false,true);
        distance = autoCat(distance,data(i).data(trial).dist2wells(:,data(i).T.foodLoc(trial)),...
                   false,true);
    end
    grouped(i).time = mean(time,2,'omitnan');
    grouped(i).temp = mean(temp,2,'omitnan');
    grouped(i).speed.all = speed;
    grouped(i).speed.avg = mean(speed,2,'omitnan');
    grouped(i).speed.err = std(speed,0,2,'omitnan')/sqrt(num.trial(i));
    grouped(i).dist.all = distance;
    grouped(i).dist.avg = mean(distance,2,'omitnan');
    grouped(i).dist.err = std(distance,0,2,'omitnan')/sqrt(num.trial(i));
    
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
    grouped(i).increasing.err = std(increasing,0,2,'omitnan')/sqrt(num.trial(i)); 
    grouped(i).increasing.rate = median(tempRates(:,upIdx));
    grouped(i).decreasing.temps = median(temperatures,2);
    grouped(i).decreasing.all = decreasing;
    grouped(i).decreasing.avg = mean(decreasing,2,'omitnan');
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

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8]; %distance from food timecourse %TODO: normalize this to something more intuitive? 
sb(3).idx = [10,11,13,14]; %speed timecourse
sb(4).idx = 3:3:15; %binned distance alignment

LW = 1.5;
sSpan = 180;

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

    %speed
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
        plot(x,y,'LineWidth',LW,'Color',kolor,'linestyle','-')
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
        plot(x,y,'LineWidth',LW,'Color',kolor,'linestyle','--','HandleVisibility','off');

        % Names and Colors of included data
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,true,[r,c],sb);
% temp
subplot(r,c,sb(1).idx) 
ylabel('\circC')
set(gca,"XColor",'k')
% distance
subplot(r,c,sb(2).idx) 
ylabel('distance (mm)')
set(gca,"XColor",'k')
% speed
subplot(r,c,sb(3).idx) 
ylabel('speed (mm/s)')
xlabel('time (min)')
% temp rate 
subplot(r,c,sb(4).idx) 
ylabel('distance (mm)')
xlabel('temp (\circC)')
% 
legend(dataString,'textcolor', 'w', 'location', 'northeast', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' timecourse summary'],'-png');

%% FIGURE: cumulative hysteresis for each genotype / trial

clearvars('-except',initial_vars{:})
LW = 0.75;
buff = 0.2;
SZ = 80;
r = 1; %rows
c = 2; %columns


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

% Cumulative difference in proximity
subplot(r,c,2)
hold on
for i = 1:num.exp
    kolor = grouped(i).color;
    y = grouped(i).decreasing.all-grouped(i).increasing.all;
    plotY = sum(y,1,'omitnan');
    x = shuffle(linspace(i-buff,i+buff,num.trial(i)));
    scatter(x,plotY,SZ,kolor,"filled","o")
    plot([i-buff,i+buff],[mean(plotY),mean(plotY)],'color','w','LineWidth',2)
end
xlim([0.5,num.exp+0.5])
h_line(0,'w',':',1)
xlabel('Group')
set(gca,'XTick',1:num.exp)

ylabel('cumulatice distance difference (mm)')

formatFig(fig,true,[r,c]);

% save figure
save_figure(fig,[saveDir expGroup ' hysteresis summary'],'-png');

%%  FIGURE: ramp by ramp hysteresis comparision

clearvars('-except',initial_vars{:})


% FIGURE: plot the heating vs cooling plots for each of the four ramps
% acros each of the different trials
buff = 0.3;
LW = 0.75;

fig = figure; set(fig,'pos',[1998 582 431 497]); hold on
    for i = 1:num.exp
        kolor = grouped(i).color;
        nRamps = size(mat(i).hysteresis,2);
        x = repmat(linspace(i-buff,i+buff,nRamps),[num.trial(i),1]);
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


%% Event-aligned comparisons
%TODO
clearvars('-except',initial_vars{:})

timeROI = 25; % how many minutes to look at behavior after each event
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




% SINGLE EXPERIMENT COMPARISON
x = linspace(0,timeROI,duration+1);
LW = 1.5;

for i = 1:num.exp

    fig = figure; set(fig,'pos',[2130 275 428 534])
    hold on
    for ss = 1:length(sections)
        y = grouped(i).aligned.([sections{ss} '_MEAN']);
        y_err = grouped(i).aligned.([sections{ss} '_SEM']);
    
        fill_data = error_fill(x, y, y_err);
        h = fill(fill_data.X, fill_data.Y, Color(s_color{ss}), 'EdgeColor','none','HandleVisibility','off');
        set(h, 'facealpha', 0.4)
        plot(x, y,'color',Color(s_color{ss}),'LineWidth',LW)
    end
    
    xlabel('time (min)')
    ylabel('distance from food (mm)')
    formatFig(fig,true);
    title([grouped(i).name],'Color','w','FontSize',12,'FontName','times')
    
    save_figure(fig,[saveDir grouped(i).name ' event aligned distance'],'-png',true);
end

% CROSS EXPERIMENT COMPARISION (WITHIN HEATING,COOLING,HOLDING)












































