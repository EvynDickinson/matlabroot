
%    GroupDataGUI

%% LOAD: multiple trials that are grouped into a structure
clear; warning off

if strcmpi('Yes', questdlg('Use excel named structure?','','Yes','No', 'Cancel', 'No'))
    % baseFolder = getCloudPath;
    [excelfile, Excel, xlFile] = load_QuadBowlExperiments;
    baseFolder = getCloudPath;
    
    % Select structure to load:
    [~,~,structInfo] = getExcelStructureNames(true);
    ExpGroup = structInfo.StructName;
    ntrials = structInfo.numTrials;
    
    % Make a data folder for structure images
    figDir = [baseFolder 'Data structures/' ExpGroup '/'];
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
        dirc = [baseFolder, trialDate, '/Arena ' trialArena '/analysis/' trialExpID trialArena  ' timecourse data.mat'];
        
        % load data
        todel = load(dirc);
        variList = fieldnames(todel);
    %     todel = load(dirc, varList{:});
        data(trial).trialName = trialName;
        for ii = 1:length(variList)
            data(trial).(variList{ii}) = todel.(variList{ii});
        end
    end

else %select data structure from folder names out of GroupDataGUI:
    [baseFolder, folder] = getCloudPath(3);
    list_dirs = dir(folder); 
    list_dirs = {list_dirs(:).name};
    list_dirs(1:2) = [];
    idx = listdlg('ListString', list_dirs,'ListSize', [250, 400]);
    
    figDir = [folder list_dirs{idx} '/'];
    ExpGroup = list_dirs{idx};
    
    % load the directory list:
    load([figDir 'fileList.mat'])
    
    % extract information
    ntrials = size(T,1);
    dates = T.Date;
    arenas = T.Arena;
    expID = T.ExperimentID;
    
    % load data
    n = cell(1,ntrials);
    data = struct('dist2wells', n, 'wellLabels', n, 'occupancy', n,'nflies',n,'data', n);  
    var2load = fieldnames(data);
    
    for trial = 1:ntrials
        filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/analysis/'];
        if ~isfolder(filePath)
            filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/'];
        end
        
        temp = load([filePath expID{trial} arenas{trial} ' timecourse data.mat'], var2load{:});
        for ii = 1:length(var2load)
            try data(trial).(var2load{ii}) = temp.(var2load{ii});
            catch
                try data(trial).(var2load{ii}) = temp.data.(var2load{ii});
                catch
                    if strcmp(var2load{ii},'dist2wells')
                        try data(trial).(var2load{ii}) = temp.data.dist2well;
                        catch; data(trial).(var2load{ii}) = [];
                        end
                    else
                        data(trial).(var2load{ii}) = [];
                    end
                end
            end
        end

        % load speed data
        temp = load([filePath expID{trial} ' speed data.mat']);
        data(trial).speed = temp.speed;
        disp([expID{trial} arenas{trial}])
    end
end

%Pull and reorganize data within the structure:

% pix2mm = 12.8; %conversion from pixels to mm for these videos
initial_vars = {'ExpGroup','baseFolder', 'T', 'data', 'figDir', 'filePath',...
                'initial_vars', 'folder', 'ntrials', 'pix2mm'};
clearvars('-except',initial_vars{:})
if questdlg('Save loaded data?')
    save([figDir ExpGroup ' raw'],'-v7.3')
end
fprintf('Data loaded\n')

%% ANALYSIS: General data organization (forward and backwards compatability
if ~exist('pix2mm','var')
    pix2mm = 12.8;
end
% fill out the data structure
for trial = 1:ntrials
    if isempty(data(trial).occupancy)
        data(trial).occupancy = data(trial).data.occupancy;
    end
    if isempty(data(trial).nflies)
        data(trial).nflies = data(trial).data.nflies;
    end
    if isempty(data(trial).wellLabels)
        data(trial).wellLabels = data(trial).data.wellLabels;
    end
    if isempty(data(trial).dist2wells)
        try
            data(trial).dist2wells = data(trial).data.occupancy.dist2wells;
        catch
            for well = 1:4
                data(trial).dist2wells(:,well) = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm;
            end
        end
    elseif isfield(data(trial).dist2wells, 'N')
        temp = data(trial).dist2wells;
        data(trial).dist2wells = [];
        for well = 1:4
            data(trial).dist2wells(:,well) = temp(well).N(:,1)./pix2mm;
        end
    end
end

% group names into 'classes' and extract the food group identities
[foodName,foodLoc,foodCat] = deal([]);
bin_foods = {'Plant', 'Molasses', 'Glucose', 'Water', 'ACV', 'Sugar','Merlot','glucose','sucrose','Caviar'};
for trial = 1:ntrials
    loc = ~strcmp(data(trial).wellLabels,'Empty');
    foodName{trial,1} = data(trial).wellLabels{loc};
    foodLoc(trial,1) = find(loc);
    % find the 'category' of food:
    for ii = 1:length(bin_foods)
        if strfind(foodName{trial}, bin_foods{ii})
            foodCat{trial,1} = bin_foods{ii};
        end
    end
end

disp(table(foodName, foodLoc, foodCat))
% TODO: add a check here to not double add the columns of data
T = addvars(T, foodName, foodLoc, foodCat);
foodNames = unique(foodName);
nfoods = length(foodNames);

initial_vars{end+1} = 'nfoods';
initial_vars{end+1} = 'foodNames';
clearvars('-except',initial_vars{:})

%% FIGURE: visual check of temperature alignment across the experiments:
fig = figure; set(fig, 'position', [107 595 1150 235]); hold on
for trial = 1:ntrials
    X = data(trial).occupancy.time;
    Y = data(trial).occupancy.temp;
    plot(X, Y, 'linewidth', 1) 
end
xlabel('Time (min)')
ylabel('Temp (\circ)')
title({'temperature alignment across experiments';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true);

save_figure(fig, [figDir 'temperature alignment'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

% % POSTER CODE:
% CList = {'DodgerBlue', 'Magenta'};
% ii = 0;
% fig = figure; set(fig, 'position', [107 595 1150 235]); hold on
% for trial = [1,ntrials]
%     ii = ii+1;
%     X = data(trial).occupancy.time;
%     Y = data(trial).occupancy.temp;
%     plot(X, Y, 'linewidth', 2, 'color', Color(CList{ii}))
% end
% xlabel('Time (min)')
% ylabel('Temp (\circ)')
% formatFig(fig, false);
% save_figure(fig, [figDir 'temperature alignment'], '-pdf');

%% FIGURE: Distance to food vs temperature -- ALL trial lines
sSpan = 360;
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.temp;
    for well = 1:4
        Y = data(trial).dist2wells(:,well);
        % reorder the data by temperature:
        [plotX, idx] = sort(X);
        plot(plotX, smoothdata(Y(idx),'movmean', sSpan), 'linewidth', 1,'color', pullFoodColor(data(trial).wellLabels{well}))
    end
end
ylabel('Distance to well (mm)')
xlabel('Temp (\circC)')
title({'Location from food sources by temperature';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true);
save_figure(fig, [figDir ExpGroup ' all temp vs distance'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% FIGURE: Average distance | occupancy vs temp across trials:

inputVar =  questdlg('Which data type to compare?','','distance','occupation probability','Cancel','distance');

food = struct;
switch inputVar
    case 'distance'
        ylab = 'Distance from well (mm)';
        L_loc = 'southwest';
        yLimit = [10,35];
    case 'occupation probability'
        ylab = inputVar;
        L_loc = 'northwest';
        yLimit = [0,1];
    case 'Cancel'
        return
end

%allocate empty data structure & set params
for ii = 1:nfoods+1
    food(ii).N = [];
    if ii <= nfoods
        food(ii).name = foodNames{ii};
    else 
        food(ii).name = 'Empty';
    end
end
for trial = 1:ntrials
    x = data(trial).occupancy.temp; % temperature data
    % Screen out data pre/post data 
    tempPoints = getTempTurnPoints(T.TempProtocol{trial});
    roi = sort([tempPoints.DownROI,tempPoints.UpROI]);
    for well = 1:4
        switch inputVar
            case 'distance'
                y = data(trial).dist2wells(:,well);
            case 'occupation probability'
                y = data(trial).occupancy.occ(:,well);
        end
        if well==T.foodLoc(trial)
            foodType = find(strcmp(T.foodName{trial},foodNames));
            food(foodType).N = [food(foodType).N; x(roi), y(roi)];
        else % empty well
            food(nfoods+1).N = [food(nfoods+1).N; x(roi), y(roi)];
        end
    end
end 

% Temperature range and bin size formatting
[threshHigh, threshLow] = getTempThresholds(T.TempProtocol{1});
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
if t_roi(end)<ceil(threshHigh)
    t_roi(end+1) = ceil(threshHigh) + binSpace;
end
    
% Pull plotting data:
for ii = 1:nfoods+1 %each food type
    % cut off the high and low ends of data (to clean):
    loc = food(ii).N(:,1)>threshHigh | food(ii).N(:,1)<threshLow;
    food(ii).N(loc,:) = [];
    % sort all the data by temperature:
    food(ii).N(:,3) = discretize(food(ii).N(:,1),t_roi);
    for tt = 1:length(t_roi)
        loc = food(ii).N(:,3)==tt;
        food(ii).avg(tt) = mean(food(ii).N(loc,2),1,'omitnan');
        food(ii).err(tt) = std(food(ii).N(loc,2),0,1,'omitnan')/sqrt(ntrials);
    end
end
    
% FIGURE: plot the avg. temp vs. distance data
CList = Color('deepskyblue', 'magenta', nfoods); %steelblue
CList(nfoods+1,:) = Color('White');

fig = figure; set(fig,'pos', [67 82 675 692]);
hold on
for ii = 1:nfoods+1
    kolor = CList(ii,:);
    x = t_roi(1:end-1);
    y = food(ii).avg(1:end-1);
    y_err = food(ii).err(1:end-1);
    fill_data = error_fill(x,y, y_err);
    h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
    set(h, 'facealpha', 0.35)
    plot(x, y, 'Color', kolor, 'LineWidth', 2)
%     plot(x, y+y_err, 'Color', kolor, 'LineWidth', 0.25)
%     plot(x, y-y_err, 'Color', kolor, 'LineWidth', 0.25)
end
%Labels:
% xlim([7,20])
ylim(yLimit)
xlabel('temperature (\circC)')
ylabel(ylab)
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 18)
str = strrep([foodNames; 'Empty'],'_',' ');
for i = 1:length(str)
    leg_str{(i*2)-1} = '';
    leg_str{i*2} = str{i};
end
legend(leg_str,'textcolor', 'w', 'location', L_loc, 'box', 'off','fontsize',12)

save_figure(fig, [figDir 'Temp vs ' inputVar ' bin size ' num2str(binSpace)], '-png');
clearvars('-except',initial_vars{:})

% % POSTER CODE:
% % FIGURE: plot the avg. temp vs. distance data
% CList = Color('SteelBlue', 'magenta', nfoods);
% CList(nfoods+1,:) = Color('black');
% 
% fig = figure; set(fig,'pos', [-789 398 480 717]);
% hold on
% for ii = 1:nfoods+1
%     kolor = CList(ii,:);
%     x = t_roi(1:end-1);
%     y = food(ii).avg(1:end-1);
%     y_err = food(ii).err(1:end-1);
% %     fill_data = error_fill(x,y, y_err);
% %     h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
% %     set(h, 'facealpha', 0.35)
%     plot(x, y, 'Color', kolor, 'LineWidth', 2)
%     plot(x, y+y_err, 'Color', kolor, 'LineWidth', 0.25,'HandleVisibility','off')
%     plot(x, y-y_err, 'Color', kolor, 'LineWidth', 0.25,'HandleVisibility','off')
% end
% %Labels:
% % xlim([7,20])
% ylim([10,35])
% xlabel('temperature (\circC)')
% ylabel('distance from well (mm)')
% title('')
% formatFig(fig,false);
% ax = gca;
% set(ax, 'FontSize', 18)
% str = strrep([foodNames; 'Empty'],'_',' ');
% str = {'Glucose','Molasses','Plant','Empty'};
% legend(str,'textcolor', 'k', 'location', L_loc, 'box', 'off','fontsize',12)
% save_figure(fig, ['G:\My Drive\Presentations\SRT May 2022\Temp vs ' inputVar ' bin size ' num2str(binSpace)], '-pdf');
% clearvars('-except',initial_vars{:})

%% FIGURE: Average movement | cluserting | eccentricity vs temp across trials:

inputVar =  questdlg('Which data type to compare?','','movement','clustering','eccentricity','movement');
switch inputVar
    case 'movement'
        ylab = 'movement (~mm/s)';
        L_loc = 'southeast';
    case 'clustering'
        ylab = 'Inter-fly-distance (mm)';
        L_loc = 'northwest';
    case 'eccentricity'
        ylab = 'eccentricity (mm)';
        L_loc = 'northwest';
    case ''
        return
end

%allocate empty data structure & set params
food = struct;
food.N = [];
for trial = 1:ntrials
    % Screen out data pre/post data 
    tempPoints = getTempTurnPoints(T.TempProtocol{trial});
    roi = sort([tempPoints.DownROI,tempPoints.UpROI,tempPoints.HoldROI]);
    switch inputVar
        case 'movement'
            x = data(trial).occupancy.temp(1:end-1);
            y = data(trial).occupancy.movement;
        case 'clustering'
            y = data(trial).occupancy.IFD';
            x = data(trial).occupancy.temp;
        case 'eccentricity'
            y = data(trial).occupancy.eccentricity;
            x = data(trial).occupancy.temp;
    end 
    food.N = [food.N; x(roi),y(roi)];
end 

% Temperature range and bin size formatting
[threshHigh, threshLow] = getTempThresholds(T.TempProtocol{1});
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
if t_roi(end)<ceil(threshHigh)
    t_roi(end+1) = ceil(threshHigh) + binSpace;
end
    
% cut off the high and low ends of data (to clean):
loc = food.N(:,1)>threshHigh | food.N(:,1)<threshLow;
food.N(loc,:) = [];
% sort all the data by temperature:
food.N(:,3) = discretize(food.N(:,1),t_roi);
for tt = 1:length(t_roi)
    loc = food.N(:,3)==tt;
    food.avg(tt) = mean(food.N(loc,2),1,'omitnan');
    food.err(tt) = std(food.N(loc,2),0,1,'omitnan')/sqrt(ntrials);
end

    
% FIGURE: plot the avg. temp vs. distance data
kolor = Color('Teal');

fig = figure; set(fig,'pos', [67 82 675 692]);
hold on
x = t_roi(1:end-1);
y = food.avg(1:end-1);
y_err = food.err(1:end-1);


fill_data = error_fill(x,y, y_err);
h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
set(h, 'facealpha', 0.3)
plot(x, y, 'Color', kolor, 'LineWidth', 2)
plot(x, y+y_err, 'Color', kolor, 'LineWidth', 0.25)
plot(x, y-y_err, 'Color', kolor, 'LineWidth', 0.25)


%Labels:
% xlim([7,20])
xlabel('temperature (\circC)')
ylabel(ylab)
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 18)

save_figure(fig, [figDir ExpGroup ' temp vs ' inputVar ' bin size ' num2str(binSpace)], '-png');
clearvars('-except',initial_vars{:})

% 
% % POSTER CODE:   
% % FIGURE: plot the avg. temp vs. distance data
% kolor = Color('black');
% 
% fig = figure; set(fig,'pos', [67 82 675 692]);
% hold on
% x = t_roi;
% y = food.avg./pix2mm;
% y_err = food.err./pix2mm;
% 
% plot(x, y, 'Color', kolor, 'LineWidth', 2)
% plot(x, y+y_err, 'Color', kolor, 'LineWidth', 0.25)
% plot(x, y-y_err, 'Color', kolor, 'LineWidth', 0.25)
% 
% %Labels:
% % xlim([7,20])
% xlabel('temperature (\circC)')
% ylabel(ylab)
% % ylabel(ylab)
% formatFig(fig);
% ax = gca;
% set(ax, 'FontSize', 18)
% 
% save_figure(fig, [figDir ExpGroup ' temp vs ' inputVar ' bin size ' num2str(binSpace)], '-pdf');
% clearvars('-except',initial_vars{:})

%% FIGURE: group movement histogram

mv = [];
for trial = 1:ntrials
    y = data(trial).occupancy.movement;
    mv = [mv; y];
end

fig = figure; hold on
h = histogram(mv,30);
h.FaceColor = Color('teal');
h.FaceAlpha = 1;
h.EdgeColor = 'w';
xlabel('Movement speed (~mm/s)')
ylabel('Count')
formatFig(fig,true);
set(gca, 'fontsize', 20)

save_figure(fig, [figDir 'Movement histogram'], '-png');
clearvars('-except',initial_vars{:})

%% ANALYSIS: Get temp rate information from all trials
G = struct;
[threshHigh, threshLow] = getTempThresholds(T.TempProtocol);
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
vars = [initial_vars(:)', 'trial', 'threshHigh', 'threshLow', 'binSpace','G', 'vars'];

% Get the temp rate, temp, and distance from food
for trial = 1:ntrials
%     if isfield(data(trial).data, 'TR')
%         G(trial).TR = data(trial).data.TR;
%         continue
%     end
    [TR,T_rates,plotData] = deal([]); % temp rate data will go here then get saved into Arena Data

    % temperature
    temp = data(trial).occupancy.temp; temp = reshape(temp,numel(temp),1);
    time = data(trial).occupancy.time; time = reshape(time,numel(time),1);
    plotData(:,1) = temp(1:end-1);
    dT = diff(smooth(temp,180)); %180 = 1 minute smoothing kernal (for 3fps data)
    dt = diff(time); 
    plotData(:,2) = dT./dt; 

    % distance from food
    wellLabels = data(trial).wellLabels;
    for well = 1:4
        [kolor,num] = pullFoodColor(wellLabels{well});
        if ~strcmp(wellLabels{well},'Empty')
            y = data(trial).dist2wells(:,well);
            plotData(:,3) = y(1:end-1);
            foodWell = well;
            break
        end
    end

    % movement
    plotData(:,6) = data(trial).occupancy.movement;

    % find the mean temp rate during each ramp period:
    tPoints = getTempTurnPoints(T.TempProtocol{trial}); %accomodates multiple temp protocols within the data group
    if strcmp(T.TempProtocol{trial},'linear_ramp_with_recovery_23-15')
        tPoints.rates = [tPoints.rates(1), 0, tPoints.rates(2)];
        tPoints.nRates = 3;
    end
%     threshLow = tPoints.threshLow;
%     threshHigh = tPoints.threshHigh;
    keepLoc = false(1,size(plotData,1));
%     roi = sort([tPoints.DownROI,tPoints.UpROI,tPoints.HoldROI]);
%     keepLoc(roi) = true;
    for ii = 1:length(tPoints.transitions)-1
        roi = tPoints.transitions(ii):tPoints.transitions(ii+1);
%         distD = median(plotData(roi,2));
        rateAvg = mean(plotData(roi,2));
        plotData(roi,4) = round(rateAvg*ones(range(roi)+1,1),2);
        keepLoc(roi) = true;
    end
    plotData(~keepLoc,:) = nan; % exclude data outside prescribed ROI regions
%     figure; plot(time(1:end-1),plotData(:,2)); hold on; scatter(time(keepLoc),plotData(keepLoc,2))

    % Temp-rate identification and sorting: 
    buffSize = 0.05;
    if any(ismember(strsplit(T.TempProtocol{trial},'_'),'sweeps')) %account for high freq jitter
        buffSize = 0.1;
    end
    for ii = 1:tPoints.nRates
        edges(ii,:) = [tPoints.rates(ii)-buffSize, tPoints.rates(ii)+buffSize];
    end  
    rateData = plotData(:,4);
    nRates = tPoints.nRates;
    if ~isempty(tPoints.hold) && ~any(tPoints.rates==0) %if there is hold data but no 0 rate included...
        rateIdx = discretize(rateData,nRates+1);
        nRates = nRates+1;
        % figure; plot(rateIdx)
    else
        rateIdx = discretize(rateData,nRates);
    end
    for ii = 1:nRates %TODO: adjust this to account for the zero rate of change?
        TD = rateData(rateIdx==ii);
        T_rates(ii) = round(mean(TD,'omitnan'),2);
        % do these match the assumed rates?
        idx = find(T_rates(ii) > edges(:,1) & T_rates(ii) < edges(:,2));
        if isempty(idx)

            % plot the raw vs smoothed temp rate data
            sb(1).idx = 1; sb(2).idx = 2:4; r = 4; c = 1;
            fig = figure; set(fig,'pos',[856 299 678 598]); 
                subplot(r,c,sb(1).idx); plot(time,temp,'color','w'); ylabel('temp (\circC)')
                hold on; yyaxis right; plot(time(1:end-1),rateIdx,'color','m')
                subplot(r,c,sb(2).idx); plot(time(1:end-1),diff(temp)./diff(time),'color',Color('grey'))
                hold on; plot(time(1:end-1),plotData(:,2),'color', 'yellow','linewidth',0.5)
                ylim([-2,2]); y1 = rangeLine(fig,0.15);% plot event regions
                for i = 1:size(tPoints.up,1)
                    plot(time(tPoints.up(i,:)),[y1,y1],'color','red','linewidth',2.5)
                end
                for i = 1:size(tPoints.down,1)
                    plot(time(tPoints.down(i,:)),[y1,y1],'color','blue','linewidth',2.5)
                end
%                 v_line(time(tPoints.transitions),'w','-',1)
                xlabel('time (min)'); ylabel('Temp rate of change')
            formatFig(fig,true,[r,c],sb)

            warndlg('Temp rate not within expected range')
            return
        end
        % set the rate value to the uniform cross-group value:
        rateData(rateIdx==ii) = tPoints.rates(ii);
        T_rates(ii) = tPoints.rates(ii);
%        figure; hold on; plot(rateData); plot(plotData(:,4))     
    end
    plotData(:,4) = rateData;
    plotData(:,5) = rateIdx; %rate bin
    TR.rateIdx = rateIdx;
    TR.rates = T_rates;
    TR.data = plotData;

    % Temperature range and bin size formatting
    t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
    if t_roi(end)<ceil(threshHigh)
        t_roi(end+1) = ceil(threshHigh) + binSpace;
    end
    nTemps = length(t_roi);
    tempIdx = discretize(plotData(:,1), t_roi);
    TR.nTemps = nTemps;
    TR.tempIdx = tempIdx;
    TR.temps = t_roi;
    
    
    % turn coordinates of heatmap into a vector to fill in 2D
    [HM,FC,FD,dist_mat,FC_mat,mov_mat,MM] = deal([]);
    flyCount = data(trial).occupancy.flyCount;
    if size(flyCount,2)>1 %aka new data with white arena
        flyCount = flyCount(:,Alphabet(T.Arena{trial}));
    end
    for col = 1:nTemps
      for row = 1:nRates
        % Pull the data that matches this category
        loc = (rateIdx==row) & (tempIdx==col);

        % fly count data: (for tracking assessment)
        fc = flyCount(loc);
        FC(row,col).data = fc;
        FC(row,col).mean = mean(fc,'omitnan');
        FC(row,col).err = std(fc,0,'omitnan');
        FC_mat.avg(row,col) = FC(row,col).mean;
        FC_mat.err(row,col) = FC(row,col).err;

        % distance data:
        dist = plotData(loc,3);
        HM(row,col) = mean(dist,'omitnan'); %heat map
        FD(row,col).data = dist;
        FD(row,col).mean = mean(dist,'omitnan');
        FD(row,col).err = std(dist,0,'omitnan');
        dist_mat.avg(row,col) = FD(row,col).mean;
        dist_mat.err(row,col) = FD(row,col).err;

        % movement data:
        movem = plotData(loc,6);
        MM(row,col).data = movem;
        MM(row,col).mean = mean(movem,'omitnan');
        MM(row,col).err = std(movem,0,'omitnan');
        mov_mat.avg(row,col) = MM(row,col).mean;
        mov_mat.err(row,col) = MM(row,col).err;
      end
    end
    TR.heatmap = HM;
    TR.dist_mat = dist_mat;
    TR.movement = mov_mat;
    TR.FC_mat = FC_mat;
    TR.nRates = length(TR.rates);
    G(trial).TR = TR;

clearvars('-except',vars{:}) 
end

%% FIGURE: Temp hysteresis - distance to food | movement 
clearvars('-except',vars{:}) 
dataType =  questdlg('Which data type to compare?','','distance','movement','Cancel','distance');

switch dataType
    case 'distance'
        ylab = 'Distance from well (mm)';
        L_loc = 'southwest';
    case 'movement'
        ylab = 'Movement (au)';
        L_loc = 'northwest';
    case 'Cancel'
        return
end

% Find the total number and id of temp rates:
allRates=[];
for trial = 1:ntrials
    allRates = [allRates,G(trial).TR.rates];
end
tRates = sort(unique(allRates));
nRates = length(tRates);
nTemps = G(1).TR.nTemps; %these should all be the same since they're held constant above

% Group the data for each temp rate
tempData = nan(nRates,nTemps,ntrials);
for trial = 1:ntrials
    for rr = 1:nRates
        idx = find(G(trial).TR.rates==tRates(rr));
        if isempty(idx)
            continue
        end
        if strcmp(dataType,'distance')
            tempData(rr,:,trial) = G(trial).TR.heatmap(idx,:);
        elseif strcmp(dataType,'movement')
            tempData(rr,:,trial) = G(trial).TR.movement.avg(idx,:);
        end
    end
end

% Find the mean and err of each temp bin:
plotData.avg = mean(tempData,3,'omitnan');
plotData.err = std(tempData,0,3,'omitnan')./sqrt(ntrials);

heatMapData = plotData.avg;
t_roi = G(trial).TR.temps;
title_str = [ExpGroup ' (n = ' num2str(ntrials) ')'];

% ========= HeatMap of dT/dt vs T =============
fig = figure; set(fig, 'pos', [560 127 983 417]);
    hold on
    imAlpha=ones(size(heatMapData));
    imAlpha(isnan(heatMapData))=0;
    imagesc(heatMapData,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    axis tight
    title(title_str)
    % Axes formatting
    ax = gca;
    fig = formatFig(fig, true);
    XtickNum = ax.XTick;
    ax.XTickLabel = t_roi(XtickNum);
    YtickNum = ax.YTick;
    try ax.YTickLabel = tRates(YtickNum);
        
    catch
        ax.YTick = 1:length(tRates);
        ax.YTickLabel = tRates;
    end
    ylabel('\DeltaT/dt (\circC/min)')
    xlabel('Temp (\circC)')
    % Colorbar formatting
    cbh = colorbar(); 
    cbh.Label.String = ylab;
    cbh.Color = Color('white');
    % flip colormap around to make yellow closer to food
    if strcmp(dataType,'distance')
        cmap = colormap;
        set(gca, 'colormap', flip(cmap))
    end

save_figure(fig, [figDir 'Temp_rate ' dataType ' heatmap'], '-png');

% % POSTER CODE
% % ========= HeatMap of dT/dt vs T =============
% fig = figure; set(fig, 'pos', [560 127 983 417]);
%     hold on
%     imAlpha=ones(size(heatMapData));
%     imAlpha(isnan(heatMapData))=0;
%     imagesc(heatMapData,'AlphaData',imAlpha);
%     set(gca,'color',[1 1 1]);
%     axis tight
%    
%     % Axes formatting
%     ax = gca;
%     fig = formatFig(fig);
%     XtickNum = ax.XTick;
%     ax.XTickLabel = t_roi(XtickNum);
%     YtickNum = ax.YTick;
%     try ax.YTickLabel = tRates(YtickNum);
%         
%     catch
%         ax.YTick = 1:length(tRates);
%         ax.YTickLabel = tRates;
%     end
%     ylabel('\DeltaT/t (\circC/min)')
%     xlabel('temperature (\circC)')
%     % Colorbar formatting
%     cbh = colorbar(); 
%     cbh.Label.String = ylab;
%     cbh.Color = Color('black');
%     % flip colormap around to make yellow closer to food
%     if strcmp(dataType,'distance')
%         cmap = colormap;
%         set(gca, 'colormap', flip(cmap))
%     end
% save_figure(fig, [figDir 'Temp_rate ' dataType ' heatmap'], '-pdf');
% 
% ========== Line plots of each rate comparison ==============
LS = {'--','-.','-'}; %cooling|stationary|heating

fig = figure;
hold on
idx = 0; str = [];
for rr = 1:nRates
    if tRates(rr)>0
        lstyle = LS{3};
        kolor = Color('red');
    elseif tRates(rr)<0
        lstyle = LS{1};
        kolor = Color('dodgerblue');
    elseif tRates(rr)==0
        continue
        lstyle = LS{2};
        kolor = Color('white');
    end
%     kolor = pullFoodColor(tRates(rr));
    x = t_roi;
    y = plotData.avg(rr,:);
    y_err = plotData.err(rr,:);
    %sort nans to allow for error fill plot
    keep_loc = ~isnan(y) | ~isnan(y_err);
    x = x(keep_loc); y = y(keep_loc); y_err = y_err(keep_loc);
    if nRates<4
        fill_data = error_fill(x, y, y_err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor', 'none');
        set(h, 'facealpha', 0.3)
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        idx = idx + 1;
        str{idx} = [num2str(tRates(rr)) '\circC/min'];
    end
    plot(x,y,'color', kolor, 'linewidth', 2.5, 'linestyle', lstyle)%Color(cList{rr})
%     if nRates<4
%         plot(x,y+y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%         plot(x,y-y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%     end
end
% title(ExpGroup)
ylabel(ylab)
xlabel('Temperature (\circC)')
formatFig(fig, true);
% legend:
% idx = 0; str = [];
% for ii = 1:nRates
%     idx = idx+1;
%     str{idx} = [num2str(tRates(ii)) '\circC/min'];
%     if nRates<4
%         str{idx+1} = ''; str{idx+2} = ''; idx = idx+2;
%     end
% end
legend(str,'textcolor', 'w', 'location',L_loc, 'box', 'off')
set(gca,'fontsize', 14)
%Save figure
save_figure(fig, [figDir 'temp_rate ' dataType ' all rates demo'], '-png');
% 

% 
% % ---- POSTER CODE -----
% LS = {'--','-.','-'}; %cooling|stationary|heating
% fig = figure; set(fig, 'pos',[-768 538 341 645])
% hold on
% for rr = 1:nRates
%     if tRates(rr)>0
%         lstyle = LS{3};
%         kolor = Color('darkred');
%     elseif tRates(rr)<0
%         lstyle = LS{1};
%         kolor = Color('navy');
%     elseif tRates(rr)==0
%         continue
%         lstyle = LS{2};
%         kolor = Color('white');
%     end
% %     kolor = pullFoodColor(tRates(rr));
%     x = t_roi;
%     y = plotData.avg(rr,:);
%     y_err = plotData.err(rr,:);
% %     plot(x,y,'color', kolor, 'linewidth', 2.5)%Color(cList{rr})
%     plot(x,y,'color', Color('darkgreen'), 'linewidth', 2.5,'linestyle', lstyle)
%     if nRates<4
% %         plot(x,y+y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
% %         plot(x,y-y_err,'color', kolor, 'linewidth', 0.5,'linestyle', lstyle)
%         
%         plot(x,y+y_err,'color', Color('darkgreen'), 'linewidth', 0.5, 'linestyle', lstyle)
%         plot(x,y-y_err,'color', Color('darkgreen'), 'linewidth', 0.5,'linestyle', lstyle)
%     end
% end
% % title(ExpGroup)
% ylabel('movement (~mm/s)')
% % ylabel('distance from food (mm)')
% xlabel('temperature (\circC)')
% formatFig(fig, false);
% axis tight
% xlim([threshLow-2, threshHigh+2])
% set(gca,'fontsize', 15)
% ylim([10,28])
% %Save figure
% save_figure(fig, [figDir 'temp_rate ' dataType ' Berlin vs Cantons all rates demo'], '-pdf');

% save_figure(fig, 'G:\My Drive\Jeanne Lab\Presentations\Neuroscience retreat poster 5.25.2022\Berlin vs Cantons movement all rates demo', '-pdf');

% ========== Line plots of separated by rate MANUAL ADJUST FOR MORE RATES ==============
% TODO: adjust this to account for actual data...
% if nRates > 3

%     LS = {'--','-.','-'}; %cooling|stationary|heating
%     % legStr = {'SEM','Cooling','','','SEM','Heating','',''};
%     legStr = {'Cooling','','','Heating','',''};
%     ncol = 2;
%     nrow = 2;
% 
%     fig = figure; set(fig, 'pos',[296 35 1224 956])
%     % All trials
%         subplot(nrow, ncol, 1); hold on
%         for rr = 1:nRates
%             kolor = pullFoodColor(tRates(rr));
%             if tRates(rr)>0
%                 lstyle = LS{3};
%             elseif tRates(rr)<0
%                 lstyle = LS{1};
%             elseif tRates(rr)==0
%                 lstyle = LS{2};
%                 continue
%             end
%             x = t_roi;
%             y = plotData.avg(rr,:);
%             plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
%         end
%         title([ExpGroup ' (n = ' num2str(ntrials) ')'])
%         ylabel('Distance to food (mm)')
%         xlabel('Temperature (\circC)')
%     % Fast
%         subplot(nrow, ncol, 2); hold on
%         for rr = [1,7]
%             kolor = pullFoodColor(tRates(rr));
%             if tRates(rr)>0
%                 lstyle = LS{3};
%             elseif tRates(rr)<0
%                 lstyle = LS{1};
%             elseif tRates(rr)==0
%                 lstyle = LS{2};
%                 continue
%             end
%             x = t_roi(1:end-1);
%             y = plotData.avg(rr,1:end-1);
%     %         kolor = Color(cList{rr});
%             err = plotData.err(rr,1:end-1);
%     %         fill_data = error_fill(x, y, err);
%     %         h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%     %         set(h, 'facealpha', 0.2)
%             plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
%             plot(x,y+err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%             plot(x,y-err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%         end
%         title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
%         ylabel('Distance to food (mm)')
%         xlabel('Temperature (\circC)')
%         l = legend(legStr);
%         set(l,'textcolor', 'w','box', 'off')
%     % Medium
%         subplot(nrow, ncol, 3); hold on
%         for rr = [2,6]
%             if tRates(rr)>0
%                 lstyle = LS{3};
%             elseif tRates(rr)<0
%                 lstyle = LS{1};
%             elseif tRates(rr)==0
%                 lstyle = LS{2};
%                 continue
%             end
%             x = t_roi(1:end-1);
%             y = plotData.avg(rr,1:end-1);
%     %         kolor = Color(cList{rr});
%             kolor = pullFoodColor(tRates(rr));
%             err = plotData.err(rr,1:end-1);
%     %         fill_data = error_fill(x, y, err);
%     %         h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%     %         set(h, 'facealpha', 0.2)
%             plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
%             plot(x,y+err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%             plot(x,y-err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%         end
%         title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
%         ylabel('Distance to food (mm)')
%         xlabel('Temperature (\circC)')
%         l = legend(legStr);
%         set(l,'textcolor', 'w','box', 'off')
%     % Slow
%         subplot(nrow, ncol, 4); hold on
%         for rr = [3,5]
%             if tRates(rr)>0
%                 lstyle = LS{3};
%             elseif tRates(rr)<0
%                 lstyle = LS{1};
%             elseif tRates(rr)==0
%                 lstyle = LS{2};
%                 continue
%             end
%             x = t_roi(1:end-1);
%             y = plotData.avg(rr,1:end-1);
%     %         kolor = Color(cList{rr});
%             kolor = pullFoodColor(tRates(rr));
%             err = plotData.err(rr,1:end-1);
%     %         fill_data = error_fill(x, y, err);
%     %         h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%     %         set(h, 'facealpha', 0.2)
%             plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
%             plot(x,y+err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%             plot(x,y-err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
%         end
%         title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
%         ylabel('Distance to food (mm)')
%         xlabel('Temperature (\circC)')
%         l = legend(legStr);
%         set(l,'textcolor', 'w','box', 'off')
% 
%     fig = formatFig(fig, true,[nrow,ncol]);
%     save_figure(fig, [figDir ExpGroup ' temp rate food preference dependence ' ExpGroup], '-png');
% end
clearvars('-except',vars{:})

%% FIGURE: Temp hysteresis across different foods and summed for all rates

foodCat = unique(T.foodCat);
n_food = length(foodCat);

% Find the total number and id of temp rates:
allRates=[];
for trial = 1:ntrials
    allRates = [allRates,G(trial).TR.rates];
end
tRates = sort(unique(allRates));
nRates = length(tRates);
nTemps = G(1).TR.nTemps; %these should all be the same since they're held constant above
abs_rates = unique(abs(tRates));

%allocate empty data structure & set params
HM = struct;
for ii = 1:n_food % rotate through food types
    % Group the data for each temp rate & food type:
    tempData = nan(nRates,nTemps,1);
    loc = 0;
    for trial = 1:ntrials
        % is this trial within the food bin?
        if strcmp(foodCat{ii},T.foodCat{trial})
            loc = loc + 1;
            for rr = 1:nRates
                idx = find(G(trial).TR.rates==tRates(rr));
                if isempty(idx)
                    continue
                end
                tempData(rr,:,loc) = G(trial).TR.heatmap(idx,:);
            end
        end
    end
    % Find the mean and err of each temp bin for this food:
    HM(ii).all = tempData;
    HM(ii).avg = mean(tempData,3,'omitnan');
    HM(ii).err = std(tempData,0,3,'omitnan')./sqrt(loc);
    HM(ii).name = foodCat{ii};
end

% Find the hysteresis index for each food type:



% FIGURE: same temp rates compared across foods:
LS = {':','-.','-'}; %cooling|stationary|heating
x = G(trial).TR.temps;
LW = 2; 
buff = 0.3;
for rr = 1:length(abs_rates) % Absolute rate of change (to cmp heat v cool)
  % FIGURE FOR EACH TEMPERATURE RATE    
  rateList = find(abs(tRates)==abs_rates(rr));
  row = 1; col = 3; 
  sb(1).idx = 1:2;
  sb(2).idx = 3;
  fig = figure; 
  set(fig, 'pos', [57 55 1120 642]);
  % line graph:
  subplot(row,col,sb(1).idx)
  hold on; 
    idx = 1; str = [];
    for ii = 1:n_food % cycle through each food group
        for jj = rateList
            plotRate = tRates(jj);
            if plotRate>0
                lstyle = LS{3};
                str_tag = 'heating';
            elseif plotRate<0
                lstyle = LS{1};
                str_tag = 'cooling';
            elseif plotRate==0
                lstyle = LS{2};
                str_tag = 'holding';
            end
            kolor = pullFoodColor(foodCat{ii});
            y = HM(ii).avg(jj,:);
            y_err = HM(ii).err(jj,:);
            plot(x,y,'color', kolor, 'linewidth', LW, 'linestyle', lstyle);
            str{idx} = [foodCat{ii} ' ' str_tag]; idx = idx+1;
        end
    end
    % legends, labels, keys
    title(['\pm' num2str(abs_rates(rr)) '\circC/min temp ramps'])
    ylabel('Distance to food (mm)')
    xlabel('Temperature (\circC)')
    legend(str,'textcolor', 'w', 'location', 'northeast', 'box', 'off','fontsize', 10)
    
    %  Hysteresis index numbers
    subplot(row,col,sb(2).idx)
    hold on; 
    for ii = 1:n_food
        kolor = pullFoodColor(foodCat{ii});
        hys_idx = [];
        for trial = 1:size(HM(ii).all,3)
            % calculate hysteresis index
            up = HM(ii).all(1,:,trial);
            down = HM(ii).all(2,:,trial);
            hys_idx(trial) = sum(up-down,'omitnan');
        end
        % plot
   
        x = shuffle_data(linspace(ii-buff,ii+buff,length(hys_idx)));
        scatter(x,hys_idx,75,kolor,'filled')
        plot([ii-(buff*1.25),ii+(buff*1.25)],[mean(hys_idx),mean(hys_idx)],'color',kolor,'linewidth',LW)
    end
    % labels, legends, etc
    ax = gca;
    set(ax,'XTick',1:n_food,'XTickLabels',foodCat,'XTickLabelRotation',20)
    ylabel('distance hysteresis (mm)')
    formatFig(fig, true,[row,col],sb);
    
    subplot(row,col,sb(1).idx)
    set(gca, 'fontsize', 18)
    subplot(row,col,sb(2).idx)
    set(gca, 'fontsize', 18)
%Save figure
% save_figure(fig, ['temp hysteresis ' foodCat{ii} ' food ' num2str(abs_rates(rr))], '-png');
save_figure(fig, [figDir 'temp hysteresis ' foodCat{ii} ' food ' num2str(abs_rates(rr))], '-png');
end

clearvars('-except',vars{:})

%% FIGURE: Heat map of location within the arena at key points during the temp ramp
% use the temperatures : 8:2:22 for key points to check location of flies
clearvars('-except',vars{:})
tempList = [15,18.5,23];
% tempList = [15,18.5,23];
rateList = [0.125,-0.125];
% tempList = [8,12,17,22];
% rateList = [0.16,-0.16];
% tempList = [15,20,25,30];
% rateList = [0.15,-0.15];

% tempList = [8,12,17,22];
% rateList = [0.16,-0.16];

n = 10; % number of spatial bins
tt = 1; 
buffer = 0.5; % temperature buffer around target temperature
idx = 1;
for rr = 1:length(rateList)
  for tt = 1:length(tempList)   
    temp = tempList(tt); % temp in C
    rate = find(rateList(rr)==G(trial).TR.rates);
    HM = zeros(n);
    for trial = 1:ntrials
        % frame location selection
        tLoc = G(trial).TR.data(:,1)>=temp-buffer & G(trial).TR.data(:,1)<=temp+buffer;
        rLoc = G(trial).TR.rateIdx==rate;
        frames = tLoc & rLoc;

        % pull the position data for these frames and concatenate into a large
        % structure for this temp rate and temp location
        x = data(trial).data.x_loc(frames,:);
        y = data(trial).data.y_loc(frames,:);
        X = reshape(x,numel(x),1);
        Y = reshape(y,numel(y),1);
        pos = [X,Y];

        % get the 'square' units for partitioning space
        C = data(trial).data.centre;
        r = data(trial).data.r;
        x_edge = linspace(C(1)-r,C(1)+r,n);
        y_edge = linspace(C(2)-r,C(2)+r,n);

        % find x and y that are within each 'box'
        xInd = discretize(pos(:,1),x_edge);
        yInd = discretize(pos(:,2),y_edge);

        % find the number of flies within each spatial bin:
        for row = 1:n
            for col = 1:n
                nflies(row,col) = sum(yInd==row & xInd==col);
            end
        end

        % Rotate the matrix if needed to align to a well position of '2'
        k = T.foodLoc(trial)-2;
        B = rot90(nflies,k);

%         fig = figure; imagesc(B);
%         uiwait(fig)

        % Save matrix to the structure:
        HM = HM + B;
    end
    plotData(rr,tt).HM = HM;
    cmaps(idx,1) = min(min(HM));
    cmaps(idx,2) = max(max(HM));
    idx = idx + 1;
  end
end

% figures   
idx = 0;
fig = figure; set(fig, 'color', 'k', 'position',  [547 304 992 555]); %[32 70 1318 435]
for rr = 1:length(rateList)
  for tt = 1:length(tempList)
    idx = idx+1;
    subplot(length(rateList),length(tempList),idx)
        imagesc(plotData(rr,tt).HM)
        ax = gca;
        set(ax, 'xscale', 'log', 'yscale', 'log')
        axis tight
        set(ax, 'XColor', 'k','YColor', 'k', 'XTick', [],'YTick', []);
        c = colorbar;
        c.Label.String = 'Number of flies';
        c.Label.Color = 'w';
        c.Color = 'w';
        title([num2str(tempList(tt)) ' \circC at ' num2str(rateList(rr)) ' \circC/min'],'color', 'w')
        caxis([min(cmaps(:,1)) max(cmaps(:,2))]) 
        axis square
  end
end

save_figure(fig, [figDir 'Four temp hysteresis position heatmaps'], '-png');
clearvars('-except',vars{:})

%% FIGURE: Heat map of location within the arena at key points during the temp ramp for WHITE BACKGROUND
% save images as individual figures, not a group figure
% use the temperatures : 8:2:22 for key points to check location of flies
clearvars('-except',vars{:})

% tempList = [8,12,17,22];
% rateList = [0.16,-0.16];
% tempList = [15,20,25,30];
% rateList = [0.15,-0.15];

tempList = [8,12,17,22];
rateList = [0.16,-0.16];

n = 10; % number of spatial bins
tt = 1; 
buffer = 0.5; % temperature buffer around target temperature
idx = 1;


for rr = 1:length(rateList)
  for tt = 1:length(tempList) 
    frameCount = 0;  
    temp = tempList(tt); % temp in C
    rate = find(rateList(rr)==G(trial).TR.rates);
    HM = zeros(n);
    for trial = 1:ntrials
        % frame location selection
        tLoc = G(trial).TR.data(:,1)>=temp-buffer & G(trial).TR.data(:,1)<=temp+buffer;
        rLoc = G(trial).TR.rateIdx==rate;
        frames = tLoc & rLoc;
        frameCount = frameCount + sum(frames); % to normalize to the number of flies, not frames

        % pull the position data for these frames and concatenate into a large
        % structure for this temp rate and temp location
        x = data(trial).data.x_loc(frames,:);
        y = data(trial).data.y_loc(frames,:);
        X = reshape(x,numel(x),1);
        Y = reshape(y,numel(y),1);
        pos = [X,Y];

        % get the 'square' units for partitioning space
        C = data(trial).data.centre;
        r = data(trial).data.r;
        x_edge = linspace(C(1)-r,C(1)+r,n);
        y_edge = linspace(C(2)-r,C(2)+r,n);

        % find x and y that are within each 'box'
        xInd = discretize(pos(:,1),x_edge);
        yInd = discretize(pos(:,2),y_edge);

        % find the number of flies within each spatial bin:
        for row = 1:n
            for col = 1:n
                nflies(row,col) = sum(yInd==row & xInd==col);
            end
        end

        % Rotate the matrix if needed to align to a well position of '2'
        k = T.foodLoc(trial)-2;
        B = rot90(nflies,k);

%         fig = figure; imagesc(B);
%         uiwait(fig)

        % Save matrix to the structure:
        HM = HM + B;
    end
    plotData(rr,tt).HM = HM./frameCount;
    cmaps(idx,1) = min(min(HM));
    cmaps(idx,2) = max(max(HM));
    idx = idx + 1;
  end
end


% Find the max and min 'avg' number of flies...
ii = 0;
for rr = 1:length(rateList)
  for tt = 1:length(tempList)
      ii = ii + 1;
      flyMax(ii) = max(max((plotData(rr,tt).HM)));
  end
end



% figures   
for rr = 1:length(rateList)
  for tt = 1:length(tempList)
      
    fig = figure; set(fig, 'color', 'w', 'position', [-931 697 534 421])  

        imagesc(plotData(rr,tt).HM)
        
        axis square
        ax = gca;
        set(ax, 'xscale', 'log', 'yscale', 'log')
        axis tight
        set(ax, 'XColor', 'k','YColor', 'k', 'XTick', [],'YTick', []);      
        c = colorbar;
        c.Label.String = 'Number of flies';
        c.Label.FontSize = 15;
        c.FontSize = 13;
        c.Label.Color = 'k';
        c.Color = 'k';
        title([num2str(tempList(tt)) ' \circC at ' num2str(rateList(rr)) ' \circC/min'],'color', 'k')
        caxis([0 max(flyMax)]) 
        save_figure(fig, [figDir 'Place Diagrams/Position heatmaps rate ' num2str(rr) ' temp ' num2str(tempList(tt))], '-pdf');
  end
end

clearvars('-except',vars{:})

%% FIGURE: Plot the time course for all trials with the same temperature protocol
clearvars('-except',vars{:})

shaded_Err = true;

tempRegimes = unique(T.TempProtocol);
n = size(tempRegimes,1);
% CList = Color('black','green', n);
CList = Color('cyan','purple', n);
LW = 1;
sSpan = 180;
row = 5;
col = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;
sb(3).idx = 4:5;

fig = figure; set(fig, 'pos',[36 77 1082 625]); 
    for type = 1:n
        [time,temp,movement,distance] = deal([]);
        % pull the data for each plotting category
        for trial = 1:ntrials
            if strcmp(T.TempProtocol{trial},tempRegimes{type})
                % sort time
                curr_time = data(trial).occupancy.time;
                if isempty(time) % first time for this temp protocol
                    time = curr_time(1:end-1);
                    roi = 1:length(time);
                    temp(:,1) = data(trial).occupancy.temp(roi);
%                     distance = data(trial).occupancy.dist2wells(roi,T.foodLoc(trial));
                    distance = data(trial).dist2wells(roi,T.foodLoc(trial));
                    movement = data(trial).occupancy.movement(roi);
                else
                    if length(curr_time)<length(time) % shorter data 
                        time = curr_time;
                    end
                    roi = 1:length(time)-1;
                    % other parameters
                    try temp = [temp(roi,:),data(trial).occupancy.temp(roi)];
                    catch
                        temp = [temp(roi,:),data(trial).occupancy.temp(roi)'];
                    end
%                     distance = [distance(roi,:),data(trial).occupancy.dist2wells(roi,T.foodLoc(trial))];
                    distance = [distance(roi,:),data(trial).dist2wells(roi,T.foodLoc(trial))];
                    movement = [movement(roi,:),data(trial).occupancy.movement(roi)];
                end
                
            end
        end
        time = time(roi);
        kolor = CList(type,:);
        % Plot the temperature protocol averages 
        % TEMPERATURE
        subplot(row,col,sb(1).idx); hold on
        plot(time,smooth(mean(temp,2),sSpan),'color',kolor,'linewidth',LW)
        % DISTANCE
        subplot(row,col,sb(2).idx); hold on
        plot(time,smooth(mean(distance,2),sSpan),'color',kolor,'linewidth',LW)
        % MOVEMENT
        subplot(row,col,sb(3).idx); hold on
        plot(time,smooth(mean(movement,2),sSpan),'color',kolor,'linewidth',LW)
    end
    % labels and formatting
    formatFig(fig,true, [row,col],sb);
    subplot(row,col,sb(1).idx);
    ylabel('\circC')
    set(gca,'XColor', 'k')
%     title(ExpGroup,'color', 'w')
    subplot(row,col,sb(2).idx);
    ylabel('Distance to food (mm)')
    set(gca,'XColor', 'k')
    subplot(row,col,sb(3).idx);
    ylabel('Movement (~mm/s)')
    xlabel('Time (min)')

save_figure(fig, [figDir 'temp protocols time course single trial'], '-png');

% TODO: plot the temp distance relationship across multiple genotypes...?

 
%% FIGURE: time course comparison across genotypes
clearvars('-except',vars{:})

genotypes = unique(T.Genotype);
n = size(genotypes,1);
CList = {'BlueViolet', 'white', 'DeepPink', 'Orange', 'Lime', 'DodgerBlue', 'Teal', 'Red'};

LW = 1.5;
sSpan = 360;
row = 5;
col = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;
sb(3).idx = 4:5;

fig = figure; set(fig, 'pos',[36 77 1082 625]); 
    for type = 1:n
        [time,temp,movement,distance] = deal([]);
        % pull the data for each plotting category
        for trial = 1:ntrials
            if strcmp(T.Genotype{trial},genotypes{type})
                % sort time
                curr_time = data(trial).occupancy.time;
                if isempty(time) % first time for this temp protocol
                    time = curr_time(1:end-1);
                    roi = 1:length(time);
                    temp = data(trial).occupancy.temp(roi);
                    distance = data(trial).occupancy.dist2wells(roi,T.foodLoc(trial));
                    movement = data(trial).occupancy.movement(roi);
                else
                    if length(curr_time)<length(time) % shorter data 
                        time = curr_time;
                    end
                    roi = 1:length(time)-1;
                    % other parameters
                    temp = [temp(roi,:),data(trial).occupancy.temp(roi)];
                    distance = [distance(roi,:),data(trial).occupancy.dist2wells(roi,T.foodLoc(trial))];
                    movement = [movement(roi,:),data(trial).occupancy.movement(roi)];
                end
                
            end
        end
        time = time(roi);
        kolor = Color(CList{type});
        % Plot the temperature protocol averages 
        % TEMPERATURE
        subplot(row,col,sb(1).idx); hold on
        plot(time,smooth(mean(temp,2),sSpan),'color',kolor,'linewidth',LW)
        % DISTANCE
        subplot(row,col,sb(2).idx); hold on
        plot(time,smooth(mean(distance,2),sSpan),'color',kolor,'linewidth',LW)
        % MOVEMENT
        subplot(row,col,sb(3).idx); hold on
        plot(time,smooth(mean(movement,2),sSpan),'color',kolor,'linewidth',LW)
    end
    % labels and formatting
    formatFig(fig,true, [row,col],sb);
    subplot(row,col,sb(1).idx);
    ylabel('\circC')
    set(gca,'XColor', 'k')
    title(ExpGroup,'color', 'w')
    subplot(row,col,sb(2).idx);
    ylabel('Distance to food (mm)')
    set(gca,'XColor', 'k')
    subplot(row,col,sb(3).idx);
    ylabel('Movement (au)')
    xlabel('Time (min)')
legend(strrep(genotypes,'_',' '),'textcolor', 'w', 'box', 'off')

save_figure(fig, [figDir 'time course comparison'], '-png');

tbl = table(T.Genotype,T.TempProtocol,T.foodCat,T.NumFlies,'VariableNames',{'Genotype','Protocol','Food','N flies'});
disp(tbl)

% % POSTER CODE
% CList = {'Orange', 'mediumvioletred', 'darkgreen', 'Orange', 'Lime', 'DodgerBlue', 'Teal', 'Red'};
% 
% fig = figure; set(fig, 'pos',[-1020 575 721 614]); 
%     for type = 2:n
%         [time,temp,movement,distance] = deal([]);
%         % pull the data for each plotting category
%         for trial = 1:ntrials
%             if strcmp(T.Genotype{trial},genotypes{type})
%                 % sort time
%                 curr_time = data(trial).occupancy.time;
%                 if isempty(time) % first time for this temp protocol
%                     time = curr_time(1:end-1);
%                     roi = 1:length(time);
%                     temp = data(trial).occupancy.temp(roi);
%                     distance = data(trial).occupancy.dist2wells(roi,T.foodLoc(trial));
%                     movement = data(trial).occupancy.movement(roi);
%                 else
%                     if length(curr_time)<length(time) % shorter data 
%                         time = curr_time;
%                     end
%                     roi = 1:length(time)-1;
%                     % other parameters
%                     temp = [temp(roi,:),data(trial).occupancy.temp(roi)];
%                     distance = [distance(roi,:),data(trial).occupancy.dist2wells(roi,T.foodLoc(trial))];
%                     movement = [movement(roi,:),data(trial).occupancy.movement(roi)];
%                 end
%                 
%             end
%         end
%         time = time(roi);
%         kolor = Color(CList{type});
%         % Plot the temperature protocol averages 
%         % TEMPERATURE
%         subplot(row,col,sb(1).idx); hold on
%         plot(time,smooth(mean(temp,2),sSpan),'color','k','linewidth',LW)
%         % DISTANCE
%         subplot(row,col,sb(2).idx); hold on
%         plot(time,smooth(mean(distance,2),sSpan),'color',kolor,'linewidth',LW)
%         % MOVEMENT
%         subplot(row,col,sb(3).idx); hold on
%         plot(time,smooth(mean(movement,2),sSpan),'color',kolor,'linewidth',LW)
%     end
%     % labels and formatting
%     formatFig(fig,false, [row,col],sb);
%     subplot(row,col,sb(1).idx);
%     ylabel('\circC')
%     set(gca,'XColor', 'k')
%     title(ExpGroup,'color', 'w')
%     subplot(row,col,sb(2).idx);
%     ylabel('Distance to food (mm)')
%     set(gca,'XColor', 'k')
%     subplot(row,col,sb(3).idx);
%     ylabel('Movement (au)')
%     xlabel('Time (min)')
% legend(strrep(genotypes,'_',' '),'textcolor', 'w', 'box', 'off')
% 
% save_figure(fig, [figDir 'time course comparison'], '-pdf');

% 








%% Compare MOVEMENT plots across genotypes...
clearvars('-except',vars{:})

inputVar =  questdlg('Which data type to compare?','','distance','movement','Cancel','distance');

switch inputVar
    case 'distance'
        ylab = 'Distance from well (mm)';
%         L_loc = 'southwest';
    case 'movement'
        ylab = 'Movement (au)';
%         L_loc = 'northwest';
    case 'Cancel'
        return
end

ID_mat = nan(ntrials,3);

% what are the temp protocols: 
tempList = unique(T.TempProtocol);
nTempLists = size(tempList,1);
for ii = 1:nTempLists
    ID_mat(strcmp(T.TempProtocol,tempList(ii)),1) = ii;
end

% what are the genotypes:
genoList = unique(T.Genotype);
ngenoList = size(genoList,1);
for ii = 1:ngenoList
    ID_mat(strcmp(T.Genotype,genoList(ii)),2) = ii;
end

% Genotype Counts
for ii = 1:ngenoList
    disp([genoList{ii} ': n = ' num2str(sum(strcmp(genoList(ii),T.Genotype)))])
end

% what are the foods/well options:
foodList = unique(T.foodCat);
nfoodList = size(foodList,1);
for ii = 1:nfoodList
    ID_mat(strcmp(T.foodCat,foodList(ii)),3) = ii;
end

% Find unique instances across all trials:
ID_List = string(ID_mat(:,1));
for ii = 2:3
    ID_List = [ID_List+string(ID_mat(:,ii))];
end
IDs = unique(ID_List);
nIDs = size(IDs,1);

% Make an empty structure in which to group the data
plotData = struct;
for ii = 1:nIDs
    plotData(ii).h = []; %heating
    plotData(ii).c = []; %cooling
    plotData(ii).r = []; %temp rate list
end

% Compare across all conditions (overlays all temp rates...)
% CList = {'BlueViolet', 'DeepPink','Orange','Lime','DodgerBlue','Teal','Red',...
%          'Turquoise', 'DarkRed', 'Indigo', 'Plum'};
colors = {'BlueViolet',...
         'white',...
         'Indigo',...
         'Plum',...
         'Thistle',...
         'Teal',...
         'Turquoise',...
         'Aquamarine',...
         'Red',...
         'DarkRed',...
         'Orange',...
         'Gold'};
% CList = colors([1,5,10,8,2,3,6,11,9,4,7]);
CList = colors([1,2,5,10,8,2,3,6,11,9,4,7]);

     
% group together trials with the same identity
for trial = 1:ntrials
    id = find(strcmp(ID_List(trial),IDs));
    kolor = Color(CList{id});
    plotData(id).color = kolor;
    
    % collapse all temp rates into one:
    for ii = 1:G(trial).TR.nRates
        plotData(id).r = [plotData(id).r, G(trial).TR.rates(ii)];
        % Account for appropriate data type
        switch inputVar
            case 'distance'
                temp = G(trial).TR.heatmap(ii,:);
            case 'movement'
                temp = G(trial).TR.movement.avg(ii,:);
        end
        if G(trial).TR.rates(ii)<0 % Cooling data
            plotData(id).c = [plotData(id).c; temp];
        elseif G(trial).TR.rates(ii)>0 % Heating data
            plotData(id).h = [plotData(id).h; temp];
        end
    end
end
    
% Average across trials and 'uncode' the parameter names
for id = 1:nIDs
    plotData(id).c_avg = mean(plotData(id).c,1,'omitnan');
    plotData(id).h_avg = mean(plotData(id).h,1,'omitnan');
    % decode parameter names
    a = char(IDs(id));
    plotData(id).name = [tempList{str2double(a(1))} ' ' genoList{str2double(a(2))} ' ' foodList{str2double(a(3))}];
end
    
% FIGURE: Plot the temp-rate food proximity tuning curves
temp = G(1).TR.temps; %all trials should have same temp ids
LW = 2;
lStr = [];
fig = figure; set(fig, 'pos', [210 121 977 660])
    hold on
    ii = 1;
    for id = 1:nIDs
        plot(temp, plotData(id).c_avg,'color', plotData(id).color,'linewidth',LW,'linestyle', '--')
        plot(temp, plotData(id).h_avg,'color', plotData(id).color,'linewidth',LW,'linestyle', '-')
        % legend string
        if nIDs<5
            lStr{ii} = [plotData(id).name ' cooling'];
            lStr{ii+1} = [plotData(id).name ' heating'];
            ii = ii+2;
        else
            lStr{ii} = plotData(id).name;
            lStr{ii+1} = '';
            ii = ii+2;
        end
    end
    %labels and formatting
    xlabel('Temperature (\circC)')
    ylabel(ylab)
    formatFig(fig, true);
    set(gca,'fontsize', 18)
    xlim([floor(temp(1)-1),temp(end)+ceil(range(temp)*0.33)])
    %legend
    lStr = strrep(lStr,'_',' ');
    legend(lStr,'textcolor', 'w', 'box', 'off','fontsize', 8)
    
save_figure(fig, [figDir inputVar ' tuning overlay'], '-png');



%% Compare hysteresis across different genotypes -- UNFINISHED


clearvars('-except',vars{:}) 
dataType =  questdlg('Which data type to compare?','','distance','movement','Cancel','distance');

switch dataType
    case 'distance'
        ylab = 'Distance from well (mm)';
        L_loc = 'southwest';
    case 'movement'
        ylab = 'Movement (au)';
        L_loc = 'northwest';
    case 'Cancel'
        return
end

% Find the total number and id of temp rates:
allRates=[];
for trial = 1:ntrials
    allRates = [allRates,G(trial).TR.rates];
end
tRates = sort(unique(allRates));
nRates = length(tRates);
nTemps = G(1).TR.nTemps; %these should all be the same since they're held constant above

% Group the data for each temp rate
tempData = nan(nRates,nTemps,ntrials);
for trial = 1:ntrials
    for rr = 1:nRates
        idx = find(G(trial).TR.rates==tRates(rr));
        if isempty(idx)
            continue
        end
        if strcmp(dataType,'distance')
            tempData(rr,:,trial) = G(trial).TR.heatmap(idx,:);
        elseif strcmp(dataType,'movement')
            tempData(rr,:,trial) = G(trial).TR.movement.avg(idx,:);
        end
    end
end

% Find the mean and err of each temp bin:
plotData.avg = mean(tempData,3,'omitnan');
plotData.err = std(tempData,0,3,'omitnan')./sqrt(ntrials);

heatMapData = plotData.avg;
t_roi = G(trial).TR.temps;
title_str = [ExpGroup ' (n = ' num2str(ntrials) ')'];

% ========= HeatMap of dT/dt vs T =============
fig = figure; set(fig, 'pos', [560 127 983 417]);
    hold on
    imAlpha=ones(size(heatMapData));
    imAlpha(isnan(heatMapData))=0;
    imagesc(heatMapData,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    axis tight
    title(title_str)
    % Axes formatting
    ax = gca;
    fig = formatFig(fig, true);
    XtickNum = ax.XTick;
    ax.XTickLabel = t_roi(XtickNum);
    YtickNum = ax.YTick;
    try ax.YTickLabel = tRates(YtickNum);
        
    catch
        ax.YTick = 1:length(tRates);
        ax.YTickLabel = tRates;
    end
    ylabel('\DeltaT/dt (\circC/min)')
    xlabel('Temp (\circC)')
    % Colorbar formatting
    cbh = colorbar(); 
    cbh.Label.String = ylab;
    cbh.Color = Color('white');
    % flip colormap around to make yellow closer to food
    if strcmp(dataType,'distance')
        cmap = colormap;
        set(gca, 'colormap', flip(cmap))
    end

save_figure(fig, [figDir 'Temp_rate ' dataType ' heatmap'], '-png');

% % POSTER CODE
% % ========= HeatMap of dT/dt vs T =============
% fig = figure; set(fig, 'pos', [560 127 983 417]);
%     hold on
%     imAlpha=ones(size(heatMapData));
%     imAlpha(isnan(heatMapData))=0;
%     imagesc(heatMapData,'AlphaData',imAlpha);
%     set(gca,'color',[1 1 1]);
%     axis tight
%    
%     % Axes formatting
%     ax = gca;
%     fig = formatFig(fig);
%     XtickNum = ax.XTick;
%     ax.XTickLabel = t_roi(XtickNum);
%     YtickNum = ax.YTick;
%     try ax.YTickLabel = tRates(YtickNum);
%         
%     catch
%         ax.YTick = 1:length(tRates);
%         ax.YTickLabel = tRates;
%     end
%     ylabel('\DeltaT/t (\circC/min)')
%     xlabel('temperature (\circC)')
%     % Colorbar formatting
%     cbh = colorbar(); 
%     cbh.Label.String = ylab;
%     cbh.Color = Color('black');
%     % flip colormap around to make yellow closer to food
%     if strcmp(dataType,'distance')
%         cmap = colormap;
%         set(gca, 'colormap', flip(cmap))
%     end
% save_figure(fig, [figDir 'Temp_rate ' dataType ' heatmap'], '-pdf');
% 
% ========== Line plots of each rate comparison ==============
LS = {'--','-.','-'}; %cooling|stationary|heating

fig = figure;
hold on
for rr = 1:nRates
    if tRates(rr)>0
        lstyle = LS{3};
        kolor = Color('red');
    elseif tRates(rr)<0
        lstyle = LS{1};
        kolor = Color('dodgerblue');
    elseif tRates(rr)==0
        continue
        lstyle = LS{2};
        kolor = Color('white');
    end
%     kolor = pullFoodColor(tRates(rr));
    x = t_roi;
    y = plotData.avg(rr,:);
    y_err = plotData.err(rr,:);
    plot(x,y,'color', kolor, 'linewidth', 2.5, 'linestyle', lstyle)%Color(cList{rr})
    if nRates<4
        plot(x,y+y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
        plot(x,y-y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
    end
end
% title(ExpGroup)
ylabel(ylab)
xlabel('Temperature (\circC)')
formatFig(fig, true);
% legend:
idx = 0; str = [];
for ii = 1:nRates
    idx = idx+1;
    str{idx} = [num2str(tRates(ii)) '\circC/min'];
    if nRates<4
        str{idx+1} = ''; 
        str{idx+2} = ''; 
        idx = idx+2;
    end
end
legend(str,'textcolor', 'w', 'location',L_loc, 'box', 'off')
set(gca,'fontsize', 18)
%Save figure
save_figure(fig, [figDir 'temp_rate ' dataType ' all rates demo'], '-png');
% 
