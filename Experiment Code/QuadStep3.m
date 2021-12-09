    

  
%% Load in multiple trials that are grouped into a structure
clear

% baseFolder = getCloudPath;
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
baseFolder = getCloudPath;

% Select structure to load:
[~,~,structInfo] = getExcelStructureNames(true);
ExpGroup = structInfo.StructName;
ntrials = structInfo.numTrials;
      
% Make a data folder for structure images
figDir = [baseFolder 'Data structures\' ExpGroup '\'];
if ~exist(figDir, 'dir'); mkdir(figDir); end

% Load data from each trial in the structure
data = [];
fprintf('\nLoading trials: \n')
for trial = 1:ntrials
    % print the experiments as they are loaded
    trialExpID = excelfile{structInfo.rowNum(trial), Excel.expID};
    trialDate = excelfile{structInfo.rowNum(trial), Excel.date};
    trialArena = excelfile{structInfo.rowNum(trial), Excel.arena};
    trialName = [trialExpID ' Arena ' trialArena ' ' trialDate];
    disp(trialName)
    
    % build the path for the trial data
    dirc = [baseFolder, trialDate, '\Arena ' trialArena '\analysis\' trialExpID trialArena  ' timecourse data.mat'];
    
    % load data
    todel = load(dirc);
    variList = fieldnames(todel);
%     todel = load(dirc, varList{:});
    data(trial).trialName = trialName;
    for ii = 1:length(variList)
        data(trial).(variList{ii}) = todel.(variList{ii});
    end
end; clear varList

% parameters:
pix2mm = 12.8; %conversion from pixels to mm for these videos

initial_vars = {'baseFolder', 'data', 'ExpGroup', 'ntrials', 'initial_vars', 'figDir', 'pix2mm'};
clearvars('-except',initial_vars{:})
fprintf('Data loaded\n')

%% visual check of temperature alignment across the experiments:
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.time;
    Y = data(trial).occupancy.temp;
    plot(X, Y, 'linewidth', 1)
end
xlabel('Time (min)')
ylabel('Temp (\circ)')
title({'temperature alignment across experiments';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true)

save_figure(fig, [figDir ExpGroup ' temperature alignment'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Distance to food vs temperature (no time component --  quick and dirty)
sSpan = 240;
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.temp;
    for well = 1:4
        try Y = (data(trial).dist2wells(well).N(:,1));
        catch 
            Y = data(trial).occupancy.dist2wells(well).N(:,1);
            data(trial).dist2wells = data(trial).occupancy.dist2wells;
        end
        Y = Y./pix2mm; % convert the pixel values to mm
        % reorder the data by temperature:
        [plotX, idx] = sort(X);
        plot(plotX, smooth(Y(idx),sSpan), 'linewidth', 1,'color', pullFoodColor(data(trial).wellLabels{well}))
    end
end
ylabel('Distance to well (mm)')
xlabel('Temp (\circC)')
title({'Location from food sources by temperature';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true)

save_figure(fig, [figDir ExpGroup ' all temp vs distance'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Average across trials:

sSpan = 240;
[plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        x = data(trial).occupancy.temp';
        y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        % plant food well:
        if strfind(data(trial).wellLabels{well},'Plant')
            plant = [plant; x, y];
        end
        % yeast food well:
        if strfind(data(trial).wellLabels{well},'Yeast')
            yeast = [yeast; x, y];
        end
        % empty food well:
        if strfind(data(trial).wellLabels{well},'Empty')
            empty = [empty; x, y];
        end 
    end
end 

numFitPoints = 1000;
threshHigh = 19.88;
threshLow = 8.02;
plotData = [];
% Plant data:
for ii = 1:3
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); %#ok<*ASGLU>
    switch ii
        case 1 % plant food
            inputData = plant;
            kolor = pullFoodColor('Plant');
        case 2 % yeast food
            inputData = yeast;
            kolor = pullFoodColor('Yeast');
        case 3 % empty
            inputData = empty;
            kolor = pullFoodColor('Empty');
    end
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    coefficients = polyfit(smoothData(:,1), smoothData(:,2),1);
    xFit = linspace(smoothData(1,1), smoothData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, smoothData(:,1));
    yResid = smoothData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(smoothData(:,2))-1) * var(smoothData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(smoothData(:,1),smoothData(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(ii).smoothed = smoothData;
    plotData(ii).bestfit = [smoothData(:,1),yfit];
    plotData(ii).color = kolor;
    plotData(ii).rsq = round(rsq,2);
    plotData(ii).R = round(R,2);
end

% FIGURE: 
LW = 2;
fig = getfig; hold on
for ii = 1:3
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color)
end
for ii = 1:3
%     plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'r', 'linewidth', LW+1)
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'w', 'linewidth', LW)
end
ylim([10,40])

% Labels:
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);

% label key:
str = ['Plant :   R = ' num2str(plotData(1).R) '  R^2 = ' num2str(plotData(1).rsq)];
text(8.5,200, str, 'Color', plotData(1).color, 'FontSize', 12);
str = ['Yeast :  R = ' num2str(plotData(2).R) '  R^2 = ' num2str(plotData(2).rsq)];
text(8.5,180, str, 'Color', plotData(2).color, 'FontSize', 12);
str = ['Empty : R = ' num2str(plotData(3).R) '  R^2 = ' num2str(plotData(3).rsq)];
text(8.5,160, str, 'Color', plotData(3).color, 'FontSize', 12);

save_figure(fig, [figDir ExpGroup ' avg temp vs distance'], '-png');
clearvars('-except',initial_vars{:})

%% Plot all trials: distance vs time % WORKING HERE TODO
% ADD IN THE TEMPERATURE LINE AS WELL
sSpan = 240;
nrows = 3;
ncols = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;
frame_buffer = 20;

fig = figure; set(fig, 'pos', [96 444 660 482]);
subplot(nrows, ncols, sb(1).idx)
    hold on
    for trial = 1:ntrials
        X = data(trial).occupancy.time;
        Y = data(trial).occupancy.temp;
        plot(X, Y, 'linewidth', 1)
    end
    xlabel('Time (min)')
    ylabel('Temp (\circ)')
    title({'Distance from food sources';...
          ['N = ' num2str(ntrials)]})

subplot(nrows, ncols, sb(2).idx)
    hold on
    [yAvg(1).N,yAvg(2).N,yAvg(3).N,yAvg(4).N] = deal(nan([length(X)+frame_buffer,ntrials]));
    for trial = 1:ntrials
        X = data(trial).occupancy.time;
        for well = 1:4
            Y = (data(trial).dist2wells(well).N(:,1));
            Y = Y./pix2mm; % convert the pixel values to mm
            yAvg(well).N(1:length(Y),trial) = Y;
            plot(X, smooth(Y,sSpan), 'linewidth', 0.5,'color', pullFoodColor(data(trial).wellLabels{well}))
        end
    end
    % plot the avg:
    for well = 1:4
        Y = nanmean(yAvg(well).N,2);
        Y(isnan(Y)) = [];
        % TODO update the time component here to match the longest!! 
        plot(X, Y, 'linewidth', 2,'color', pullFoodColor(data(trial).wellLabels{well}))
    end
    ylabel('Distance (mm)')
    xlabel('Time (min)')
    
formatFig(fig, true,[nrows,ncols],sb)

save_figure(fig, [figDir ExpGroup ' all distance vs time'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')






%% Plot female vs. grouped|male trials:

if strcmpi(ExpGroup,'cooling_warming_ramp')
    female_exps = [7,8];
elseif strcmpi(ExpGroup,'warming_cooling_ramp')
    female_exps = [2,4];
end

sSpan = 240;
[F_plant, F_yeast, F_empty, plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    if ismember(trial,female_exps)
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                F_plant = [F_plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                F_yeast = [F_yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                F_empty = [F_empty; x, y];
            end 
        end
    else
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                plant = [plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                yeast = [yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                empty = [empty; x, y];
            end 
        end
    end
end 

numFitPoints = 1000;
threshHigh = 19.88;
threshLow = 8.02;
plotData = [];
% sort data
for ii = 1:6 %for each well type
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch ii
        case 1 % plant food
            inputData = plant;
            kolor = pullFoodColor('Plant');
        case 2 % yeast food
            inputData = yeast;
            kolor = pullFoodColor('Yeast');
        case 3 % empty
            inputData = empty;
            kolor = pullFoodColor('Empty');
        case 4 % FM plant
            inputData = F_plant;
            kolor = Color('PaleGreen');
        case 5 % FM yeast
            inputData = F_yeast;
            kolor = Color('LemonChiffon');
        case 6 % FM empty
            inputData = F_empty;
            kolor = Color('lightgrey');
    end
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    coefficients = polyfit(smoothData(:,1), smoothData(:,2),1);
    xFit = linspace(smoothData(1,1), smoothData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, smoothData(:,1));
    yResid = smoothData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(smoothData(:,2))-1) * var(smoothData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(smoothData(:,1),smoothData(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(ii).smoothed = smoothData;
    plotData(ii).bestfit = [smoothData(:,1),yfit];
    plotData(ii).color = kolor;
    plotData(ii).rsq = round(rsq,2);
    plotData(ii).R = round(R,2);
end

% FIGURE: 
LW = 2;
fig = getfig; hold on
for ii = 1:6
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color, 'filled')
end
for ii = 1:6
%     plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'r', 'linewidth', LW+1)
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'w', 'linewidth', LW)
end
ylim([10,40])

% Labels:
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);

% label key:
str = {'Plant mixed', 'Yeast Mixed', 'Empty Mixed',...
       'Plant female', 'Yeast female', 'Empty female'};
for ii = 1:6
    legend_str{ii} = [str{ii}  ' R = ' num2str(plotData(ii).R) '  R^2 = ' num2str(plotData(1).rsq)];
end
l = legend(legend_str); 
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.7081 0.8031 0.1952 0.1606]);

save_figure(fig, [figDir ExpGroup ' female-male temp vs distance'], '-png');
clearvars('-except',initial_vars{:})

%% Plot male and female and mixed arenas for distance vs. temp

if strcmpi(ExpGroup,'cooling_warming_ramp')
    female_exps = [7,8];
    male_exps = [5,6];
elseif strcmpi(ExpGroup,'warming_cooling_ramp')
    female_exps = [2,4];
    male_exps = [1,3];
end

sSpan = 360;
[M_plant, M_yeast, M_empty, F_plant, F_yeast, F_empty, plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    if ismember(trial,female_exps)
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                F_plant = [F_plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                F_yeast = [F_yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                F_empty = [F_empty; x, y];
            end 
        end
    elseif ismember(trial,male_exps)
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                M_plant = [M_plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                M_yeast = [M_yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                M_empty = [M_empty; x, y];
            end 
        end
    else
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                plant = [plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                yeast = [yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                empty = [empty; x, y];
            end 
        end
    end
end 

numFitPoints = 1000;
threshHigh = 19.88;
threshLow = 8.02;
plotData = [];
% sort data
for ii = 1:9 %for each well type
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch ii
        case 1 % plant food
            inputData = plant;
            kolor = Color('MediumSeaGreen');
        case 2 % yeast food
            inputData = yeast;
            kolor = pullFoodColor('Yeast');
        case 3 % empty
            inputData = empty;
            kolor = Color('white');
        case 4 % F plant
            inputData = F_plant;
            kolor = Color('PaleGreen');
        case 5 % F yeast
            inputData = F_yeast;
            kolor = Color('LemonChiffon');
        case 6 % F empty
            inputData = F_empty;
            kolor = Color('lightgrey');
        case 7 % M plant
            inputData = M_plant;
            kolor = Color('DarkGreen');
        case 8 % M yeast
            inputData = M_yeast;
            kolor = Color('DarkOrange');
        case 9 % M empty
            inputData = M_empty;
            kolor = Color('Gray');
            
    end
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    coefficients = polyfit(smoothData(:,1), smoothData(:,2),1);
    xFit = linspace(smoothData(1,1), smoothData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, smoothData(:,1));
    yResid = smoothData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(smoothData(:,2))-1) * var(smoothData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(smoothData(:,1),smoothData(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(ii).smoothed = smoothData;
    plotData(ii).bestfit = [smoothData(:,1),yfit];
    plotData(ii).color = kolor;
    plotData(ii).rsq = round(rsq,2);
    plotData(ii).R = round(R,2);
end


% FIGURE: 
nrows = 1;
ncol = 6;
sb(1).idx = 1:3;
sb(2).idx = 4:5;
sb(3).idx = 6;
SZ = 120;
LW = 4; 
fig = getfig; 
% avg data
subplot(nrows, ncol, sb(1).idx)
hold on
for ii = 1:9
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color, 'filled')
end
ylim([10,40])
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
% best fit line
subplot(nrows, ncol, sb(2).idx)
hold on
for ii = 1:9
    slope(ii).x = plotData(ii).bestfit([1,end],1);
    slope(ii).y = plotData(ii).bestfit([1,end],2);
    slope(ii).m = (slope(ii).y(2)-slope(ii).y(1))/(slope(ii).x(2)-slope(ii).x(1));
    plot(slope(ii).x, slope(ii).y, 'color', plotData(ii).color, 'linewidth', LW)
end
ylim([10,40])
xlabel('temperature (\circC)')
% ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
% slope points
subplot(nrows, ncol, sb(3).idx)
hold on
for ii = 1:9
    slope(ii).m = (slope(ii).y(2)-slope(ii).y(1))/(slope(ii).x(2)-slope(ii).x(1));
    scatter(1, slope(ii).m, SZ, plotData(ii).color, 'filled')
end
xlim([0,2])
hline(0,'w:')
% xlabel('temperature (\circC)')
ylabel('slope of best fit line')
% ylim([-0.8,0.2])

formatFig(fig,true, [nrows,ncol], sb);
set(gca, 'XColor', 'k')

subplot(nrows, ncol, sb(1).idx)
%legend color code:
str = {'Plant', 'Yeast', 'Empty'};
% for ii = 1:9
%     legend_str{ii} = [str{ii}  ': R = ' num2str(plotData(ii).R) '  R^2 = ' num2str(plotData(ii).rsq)];
% end
l = legend(str); 
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.4278 0.8448 0.0662 0.0689]);

save_figure(fig, [figDir ExpGroup ' temp vs distance by sex'], '-png');
clearvars('-except',initial_vars{:})





%% Group occupancy across the trials % still in progress!!!

% for now: average across the experiments (assuming they are TEMP locked)
[A,B,C,D] = deal([]);
for trial = 1:ntrials
    A = [A,data(trial).occupancy.occ(:,1)]; %(:,well) % well occupancy probability for each well
    B = [B, data(trial).occupancy.occ(:,2)];
    C = [C, data(trial).occupancy.occ(:,3)];
    D = [D, data(trial).occupancy.occ(:,4)];
end
% average across the trials TODO update the timing according to temp alignment
for trial = 1:ntrials
    temperature(trial).log = data(trial).tempLog(:,2);
end
% TEMP ALIGNED TIME (ONLY WORKS FOR TEMP HOLD OR SAME GROUP TRIALS
TAT = data(1).time();
temperature = data(1).occupancy.temp;
wellLabels = data(1).wellLabels;

nrows = 4;
ncols = 1;
subs(1).idx = 1;
subs(2).idx = 2:4;

fig = getfig;
subplot(nrows, ncols, subs(1).idx)
    plot(TAT,temperature, 'linewidth', 2, 'color', 'w')
    ylim([5,30])
    ylabel('Temp (\circC)')
subplot(nrows, ncols, subs(2).idx)
    sSpan = 240;
    hold on
    Y = [smooth(mean(A,2),sSpan), smooth(mean(B,2),sSpan),...
         smooth(mean(C,2),sSpan), smooth(mean(D,2),sSpan)];
    h = area(TAT,Y);
    for well = 1:4
        h(well).FaceColor = pullFoodColor(data(1).wellLabels{well});
    end
    set(gca, 'tickdir', 'out')
    l = legend(strrep(wellLabels,'_','-'));
    set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k',...
    'position', [0.7457 0.6520 0.0883 0.0772]);% [0.8780 0.8119 0.0963 0.1126])%

ylabel('Occupancy')
xlabel('Time (min)')
formatFig(fig, true, [nrows, ncols], subs)
subplot(nrows, ncols, subs(1).idx)
set(gca, 'XColor', 'k')
title(strrep(ExpGroup,'_','-'),'color', 'w')

save_figure(fig, [figDir 'Quadrant Occupation'], '-png');

clearvars('-except',initial_vars{:})

%% Temperature and occupancy relationship: [takes a few seconds to run]
wellLabels = data(1).wellLabels;
% organize the occupation by fly count in each temperature
% --> temperature for each fly at each timepoint
tic
[sumA,sumB,sumC,sumD,tempAll] = deal(zeros(30,1));
for trial = 1:ntrials
    [A,B,C,D] = deal([]);
    temperature = round(data(trial).occupancy.temp);
    raw = data(trial).occupancy.count;
    for ii = 1:size(raw,1) %for each data point
        A = [A; temperature(ii)*ones(raw(ii,1),1)];
        B = [B; temperature(ii)*ones(raw(ii,2),1)];
        C = [C; temperature(ii)*ones(raw(ii,3),1)];
        D = [D; temperature(ii)*ones(raw(ii,4),1)];
    end
    for temp = 5:30
        sumA(temp) = sumA(temp) + sum(A==temp);
        sumB(temp) = sumB(temp) + sum(B==temp);
        sumC(temp) = sumC(temp) + sum(C==temp);
        sumD(temp) = sumD(temp) + sum(D==temp);
        tempAll(temp) = tempAll(temp) + sum(temperature==temp);
    end
end
toc

% normalize the number of flies for each temperature but the total time
% spent at a given temperature:
tempFrac = tempAll./(sum(tempAll));
ratio = 1; %each point will represent 20 flies normalized for total time spent at each temperature
arenaA = round(sumA.*tempFrac./ratio);
arenaB = round(sumB.*tempFrac./ratio);
arenaC = round(sumC.*tempFrac./ratio);
arenaD = round(sumD.*tempFrac./ratio);


% increase points for each temperature

[yA,yB,yC,yD] = deal([]);
for temp = 1:length(arenaA)
    yA = [yA; temp*ones(arenaA(temp),1)];
    yB = [yB; temp*ones(arenaB(temp),1)];
    yC = [yC; temp*ones(arenaC(temp),1)];
    yD = [yD; temp*ones(arenaD(temp),1)];
end
xA = 1*ones(length(yA),1);
xB = 2*ones(length(yB),1);
xC = 3*ones(length(yC),1);
xD = 4*ones(length(yD),1);

% Plot the occupation by temperature plots
sz = 30;
fig = getfig;
    hold on
    swarmchart(xA,yA,sz, pullFoodColor(wellLabels{1}), 'filled');
    swarmchart(xB,yB,sz, pullFoodColor(wellLabels{2}), 'filled');
    swarmchart(xC,yC,sz, pullFoodColor(wellLabels{3}), 'filled');
    swarmchart(xD,yD,sz, pullFoodColor(wellLabels{4}), 'filled');

% swarmchart(x,y,sz,c)














%% Collapse Well ROIs across all trials for each temperature to create 'occupancyFrames'
gridSize = [50,50];
nvids = length(data(1).temp);

% loop through all trials
for trial = 1:ntrials
  temp = data(trial).temp; 
  for ii = 1:4 % well#
    % crop each of the frames from the video into well rois
    frames = data(trial).frame;
    mask = data(trial).quadMask(ii).mask;
    props = regionprops(~mask, 'BoundingBox');
    mask = repmat(mask,[1,1,size(frames,3)]); %mask for all frames
    frames(mask) = 0; %mask out other data outside ROI
    for vid = 1:nvids
        img = imcrop(frames(:,:,vid), props.BoundingBox); %crop to the circle
        maskedImage(:,:,vid) = img;
        roiImage(:,:,vid) = imresize(img, gridSize); %bin across pixels
    end
    well(ii).img = roiImage;
    well(ii).maskedImg = maskedImage;
    well(ii).temp = temp;
  end
  data(trial).well = well;
end

% align across videos roughly (TODO: refine this later)
Wells = [];
for roi = 1:4
    Wells(roi).name = data(1).trialData(1).params.(['well_' num2str(roi)]);
    Wells(roi).img = zeros(size(data(1).well(roi).img));
    for trial = 1:ntrials
        Wells(roi).img = Wells(roi).img + data(trial).well(roi).img;
        Wells(roi).occ(:,trial) = data(trial).y(:,roi); % determine this...
    end
    % find density range for the rois:
    M(roi) = max(Wells(roi).img,[],'all');
end

% Normalize the well density figure:
gridS = 10;
n = gridS+2;
well = [];
for roi = 1:4
    well(roi).img = Wells(roi).img / max(M);
    dummy = nan(n,n,nvids);
    img = imresize(well(roi).img, [gridS,gridS]);
    dummy(2:gridS+1,2:gridS+1,:) = img;
    well(roi).bin = dummy;
end
% grouped image:
occupancyFrames = [well(1).bin, well(2).bin;...
                   well(4).bin, well(3).bin];
               
occupancyNames = {Wells(1).name, Wells(2).name;...
                   Wells(4).name, Wells(3).name};              

% quick vid preview of any structure over time
fig = figure; set(fig, 'color', 'k');
disp(occupancyNames)
for vid = 1:nvids
    img = occupancyFrames(:,:,vid);
    imagesc(img)
    pause(0.1)
end
     






%% IDEAS for later:




























