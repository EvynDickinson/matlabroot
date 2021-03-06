
% This script loads in data from the GroupDataGUI to compare across the conditions to
% find any other trends...

%% LOAD: ACROSS DATA STRUCUTRES load data:
clear

%select data structure(s) from folder names:
[baseFolder, folder] = getCloudPath(3);
list_dirs = dir(folder); 
list_dirs = {list_dirs(:).name};
list_dirs(1:2) = [];
idx = listdlg('ListString', list_dirs);

% loop through for each specific data structure
var2load = {'dist2wells', 'wellLabels', 'occupancy'}; %variables to pull: 'data', 
nstruct = length(idx);
group = struct;
for n = 1:nstruct
    group(n).dataDir = [folder list_dirs{idx(n)} '/'];
    group(n).ExpGroup = list_dirs{idx(n)};

    % load the directory list:
    load([group(n).dataDir 'fileList.mat'])
    group(n).T = T;

    % extract information
    group(n).ntrials = size(T,1);
    dates = T.Date;
    arenas = T.Arena;
    expID = T.ExperimentID;
    
    % load data
    for trial = 1:group(n).ntrials
        filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/analysis/'];
        temp = load([filePath expID{trial} arenas{trial} ' timecourse data.mat'], var2load{:});
    %     variList = fieldnames(todel);
        for ii = 1:length(var2load)
            try data(trial).(var2load{ii}) = temp.(var2load{ii});
            catch data(trial).(var2load{ii}) = [];
            end
        end
        disp(trial)
        temp = [];
    end
    group(n).data = data;
    disp(['Group ' num2str(n) ' done'])
end

pix2mm = 12.8; %conversion from pixels to mm for these videos

initial_vars = {'group','nstruct', 'baseFolder','folder','initial_vars', 'pix2mm', 'structName','structDir'};
                
% optional save loaded data:
switch questdlg('Save loaded data?')
    case 'Yes'
        structName = inputdlg('Structure name?');
        structName = structName{:};
        structDir = [folder structName];
        if ~exist(structDir,'dir')
            mkdir(structDir)
        end
        save([structDir '/' structName])
    case 'No'
end

clearvars('-except',initial_vars{:})
fprintf('Data loaded\n')



%% Proximity index vs temp between arenas and food type:
sSpan = 240;

for gg = 1:nstruct
    [plant, yeast, empty] = deal([]);
    for trial = 1:group(gg).ntrials
        for well = 1:4
            x = group(gg).data(trial).occupancy.temp';
            try y = group(gg).data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm; %convert the data to mm
            catch y = (group(gg).data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            end
            % plant food well:
            if strfind(group(gg).data(trial).wellLabels{well},'Plant')
                plant = [plant; x, y];
            end
            % yeast food well:
            if strfind(group(gg).data(trial).wellLabels{well},'Yeast')
                yeast = [yeast; x, y];
            end
            % empty food well:
            if strfind(group(gg).data(trial).wellLabels{well},'Empty')
                empty = [empty; x, y];
            end 
        end
    end 
    group(gg).plant = plant;
    group(gg).yeast = yeast;
    group(gg).empty = empty;
end

numFitPoints = 1000;
threshHigh = 19.88;
threshLow = 8.02;
plotData = [];

% PLOT DATA:
nrows = 2;
ncols = 2;
sb(1).idx = 1;
sb(2).idx = 2;
sb(3).idx = 3;
sb(4).idx = 4;
LW = 2;
fig = getfig;  
for tt = 1:3 % subplot for each food type:   
    foodOpt = tt; % plant | yeast | empty
    kolors = {'purple', 'cyan', 'orange', 'green'}; % arenas A, B, C, D

    for gg = 1:nstruct %arenas A-D
        [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
        switch foodOpt
            case 1 % plant food
                inputData = group(gg).plant;
                plotname = 'Plant';
                kolor = pullFoodColor('Plant');
            case 2 % yeast food
                inputData = group(gg).yeast;
                plotname = 'Yeast';
                kolor = pullFoodColor('Yeast');
            case 3 % empty
                inputData = group(gg).empty;
                plotname = 'Empty';
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
        plotData(gg).smoothed = smoothData;
        plotData(gg).bestfit = [smoothData(:,1),yfit];
        plotData(gg).color = Color(kolors{gg});
        plotData(gg).rsq = round(rsq,2);
        plotData(gg).R = round(R,2);
    end    
    % distance vs temp lines of best fit
    subplot(nrows, ncols, sb(tt+1).idx); hold on
    for ii = 1:nstruct
        scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color)
    end
    for ii = 1:nstruct
        plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'w', 'linewidth', LW+1)
        plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', plotData(ii).color, 'linewidth', LW)
    end
    ylim([10,40])
    % Labels:
    xlabel('temperature (\circC)')
    ylabel('distance from well (mm)')
    title(plotname)

    % Plot the R2 value in a diff graph...
    subplot(nrows, ncols, sb(1).idx); hold on
    for ii = 1:size(plotData,2)
        Z(ii) = plotData(ii).rsq;
    end
    plot(1:4,Z,'color', kolor, 'LineStyle',':','LineWidth',LW,'Marker','*') 
end


% slope of line of best for across arenas
subplot(nrows, ncols, sb(1).idx)
ylim([-0.1,1.1])
xlim([0,5])
ax = gca;
set(ax,'YTick', 0:0.2:1,'XTick',1:4,'XTickLabel',{'A','B','C','D'})
ylabel('R^2')
xlabel('Arena')

formatFig(fig,true,[nrows,ncols],sb);


save_figure(fig, [structDir '\Dist vs temp across arenas ' plotname], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: trace every distance to each well by well location NOT well contents 
kolors = {'DeepPink', 'Orange', 'Lime', 'DodgerBlue'};

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
        plot(plotX, smooth(Y(idx),sSpan), 'linewidth', 1,'color', Color(kolors{well}))
    end
end
ylabel('Distance to well (mm)')
xlabel('Temp (\circC)')
title({'Location from well locations by temperature';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' all temp vs well location'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Average Distance vs temp across trials:
kolors = {'DeepPink', 'Orange', 'Lime', 'DodgerBlue'};
sSpan = 900;
[well_1, well_2, well_3, well_4] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        end
        % well 1:
        switch well
            case 1
                well_1 = [well_1; x, y];
            case 2
                well_2 = [well_2; x, y];
            case 3
                well_3 = [well_3; x, y];
            case 4
                well_4 = [well_4; x, y];
        end 
    end
end 

numFitPoints = 1000;
threshHigh = 19.88;
threshLow = 8.02;
plotData = [];
for ii = 1:4
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); %#ok<*ASGLU>
    switch ii
        case 1 % well 1
            inputData = well_1;
        case 2 % well 2
            inputData = well_2;
        case 3 % well 3
            inputData = well_3;
        case 4
            inputData = well_4;
    end
    kolor = Color(kolors{ii});
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
for ii = 1:4
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color,'filled')
end
for ii = 1:4
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'w', 'linewidth', LW+1)
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', plotData(ii).color,...
         'linewidth', LW, 'LineStyle', '-.')
end
ylim([10,40])

% Labels:
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 25)
l = legend({'Well 1', 'Well 2', 'Well 3', 'Well 4'});
set(l, 'TextColor', 'w')

% TODO: include histograms of the food:well location to show that it's a
% random spread?


% % label key:
% str = ['Plant :   R = ' num2str(plotData(1).R) '  R^2 = ' num2str(plotData(1).rsq)];
% text(8.5,16, str, 'Color', plotData(1).color, 'FontSize', 12);
% str = ['Yeast :  R = ' num2str(plotData(2).R) '  R^2 = ' num2str(plotData(2).rsq)];
% text(8.5,15, str, 'Color', plotData(2).color, 'FontSize', 12);
% str = ['Empty : R = ' num2str(plotData(3).R) '  R^2 = ' num2str(plotData(3).rsq)];
% text(8.5,14, str, 'Color', plotData(3).color, 'FontSize', 12);

save_figure(fig, [figDir ExpGroup ' avg temp vs well location'], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: histogram of well + food pairings

kolors = {'DeepPink', 'Orange', 'Lime', 'DodgerBlue'};
sSpan = 240;
[plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        labs = data(trial).wellLabels{well};
        if strfind(labs,'Plant')
            plant = [plant,well];
        end
        if strfind(labs, 'Yeast')
            yeast = [yeast,well];
        end
        if strfind(labs, 'Empty')
            empty = [empty,well];
        end
    end
end 


for well = 1:4
    well_ID(1,well) = sum(plant==well);
    well_ID(2,well) = sum(yeast==well);
    well_ID(3,well) = sum(empty==well);
end
%normalize
kolors = {'Plant', 'Yeast', 'Empty'};
wellID_per = round((well_ID./sum(well_ID,2))*100);
% PLOT
fig = figure;
hold on
for ii = 1:3
%     plot(1:4,wellID_per(ii,:),'Marker', '*', 'Color', pullFoodColor(kolors{ii}), 'LineWidth', 2)
    plot(1:4,well_ID(ii,:),'Marker', '*', 'Color', pullFoodColor(kolors{ii}), 'LineWidth', 2)
end
xlim([0,5])
ylim([0,30])
set(gca, 'XTick', 1:4)
xlabel('Well')
ylabel('Instances')
title({'Well location compared to food identity'; '  '})
fig = formatFig(fig, true);
l = legend(kolors);
set(l, 'TextColor', 'w')


% save_figure(fig, [figDir ExpGroup ' well to food identity spread'], '-png');
save_figure(fig, [figDir ExpGroup ' well to food identity spread absolute'], '-png');

% ---------------------------------
% Sum of each within wells 1+4 vs 2+3:
pair_1 = [1,4];
pair_2 = [2,3];
P_1 = sum(well_ID(:,pair_1),2);
P_2 = sum(well_ID(:,pair_2),2);

k_list = {'green', 'gold', 'grey'};
fig = figure;
h = bar([P_1';P_2']);
for ii = 1:3
    h(ii).FaceColor = Color(k_list{ii});
end
ylabel('Instances')
set(gca, 'XTickLabels', {'Wells 1 & 4', 'Wells 2 & 3'},'FontSize', 20)
fig = formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' well counts paired'], '-png');

clearvars('-except',initial_vars{:})

%%  FIGURE: bin preference by temperature... might to be cleaner to look at...















