
% Average Center of Mass of the flies in the arena

%% ANALYSIS AND FIGURES: COM of postions for each trial
% vectors to fly 'mass' in arena at different temperatures from food
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
time_bins = 1;      % 1 minute time bins for position
color_bin = 0.25;   % how many degrees (C) for a color bin

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
    bin = time_bins*f2m;
    for c = 1:size(temp,2)
        smoothData(:,c) = smooth(temp(:,c),'moving',bin);
    end
    dsData = smoothData(1:bin:end,:); % select one point for each minute

    % get temp color distribution (1/4 degree incrememts)
    tp = getTempTurnPoints(data(i).temp_protocol);
    cMapRange = tp.threshLow:color_bin:tp.threshHigh;
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
        fig = getfig('',true,[484 384]);
        hold on
            scatter(dsData(:,2), dsData(:,3), 5, kolors,'filled')
            for well = 1:4
                C = data(i).data(trial).data.wellcenters(:,well);
                scatter(C(1),C(2),SZ,Color('grey'),'filled')
            end
            scatter(well_C(1),well_C(2),SZ,'y','filled')
            viscircles([arena_C(1),arena_C(2)],r,'Color',foreColor);
            scatter(position(:,1),position(:,2),15,cMap,'filled','o','MarkerEdgeColor',foreColor)
            axis equal square;
            formatFig(fig,true);
            set(gca,'XColor',backColor,'YColor',backColor);
            % set color bar information
            colorData = uint8(round(cMap.*255)); % convert color map to uint8
            colormap(colorData);
            c = colorbar('color',foreColor);
            set(c,'Ticks',[0,1],'TickLabels',[tp.threshLow,tp.threshHigh]);
            c.Label.String = 'Temperature (\circC)';
            c.Label.VerticalAlignment = "bottom";
        title([data(i).T.Genotype{trial} ' ' data(i).T.Date{trial}],'color',foreColor)

        figFolder = [data(i).figDir 'avg position by trial/'];
        if ~exist(figFolder,'dir')
            mkdir(figFolder);
        end
        save_figure(fig,[figFolder data(i).T.Date{trial} ' ' data(i).T.ExperimentID{trial} ' ' data(i).T.Arena{trial}],fig_type,true);
    end
  end
end

%% FIGURE: Register position COM to common frame
% TODO -- account for the new camera orientation here!!! 1.24.24
clearvars('-except',initial_vars{:})
OG_Orientation = datetime('10.20.2023','InputFormat','MM.dd.yyyy');

autoSave = true; %autosave the resulting figures
[foreColor,backColor] = formattingColors(blkbgd);

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
        
        trial_date = datetime(data(i).T.Date{trial},'InputFormat','MM.dd.yyyy');

        if trial_date > OG_Orientation %new camera orientation
            % Rotate to correct orientation
            switch wellLoc
                case 1
                    plotData(:,1) = X;
                    plotData(:,2) = Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
                case 2
                    plotData(:,1) = Y;
                    plotData(:,2) = -X;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 3
                    plotData(:,1) = X;
                    plotData(:,2) = -Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 4
                    plotData(:,1) = -Y;
                    plotData(:,2) = X;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
            end
        else % Rotate to correct orientation with older camera arrangement
            
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

        end

        % PLOT
        scatter(WELLS(1:4,1),WELLS(1:4,2),SZ,foreColor,'filled')
        scatter(WELLS(wellLoc,1),WELLS(wellLoc,2),SZ,'green','filled')
        scatter(plotData(:,1),plotData(:,2),15,kolor,'filled')

    end

    viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color',grouped(i).color);
%         viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color',foreColor);
    axis square; axis equal;
    formatFig(fig,blkbgd);
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
    save_figure(fig,[saveDir expGroup grouped(i).name ' COM position'],fig_type,autoSave);

end

%% ANALYSIS: COM of postions for each trial divided by heating / cooling periods
% vectors to fly 'mass' in arena at different temperatures from food
clearvars('-except',initial_vars{:})
time_bins = 1;      % 1 minute time bins for position
color_bin = 0.25;   % how many degrees (C) for a color bin
types = {'hold','down', 'up'};

for i = 1:num.exp
  for trial = 1:num.trial(i)

    % get arena information
    well_loc = data(i).T.foodLoc(trial);
    C = data(i).data(trial).data.centre;
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

    % get temp color distribution
    tp = getTempTurnPoints(data(i).T.TempProtocol{1});
    cMapRange = tp.threshLow:color_bin:tp.threshHigh;

    bin = time_bins*f2m;
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

%% FIGURES: COM of postions for each trial divided by heating / cooling periods


%% ANALYSIS & FIGURES: heatmap of fly position within arena at various points in time


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

    



%    % BINNED
%     [tempRates,decreasing,increasing,temperatures] = deal([]);
%     for trial = 1:num.trial(i)
%         % Account for multiple numbers of rate trials
%         rates = data(i).G(1).TR.rates;
%         nRates(i) = size(rates,2);
%         if nRates(i)==2
%             blankdata = data(i).G(trial).TR.dist_mat.avg;
%         elseif nRates == 3
%             blankdata = data(i).G(trial).TR.dist_mat.avg;
%             blankdata(rates==0,:) = [];
%             rates(rates==0) = [];
%         else
%             warndlg('Temp protocol has too many rates for this analysis')
%             return
%         end
%         downIdx = find(rates<0);
%         upIdx = find(rates>0);
%         decreasing(:,trial) = blankdata(downIdx,:);
%         increasing(:,trial) = blankdata(upIdx,:);
%         tempRates = autoCat(tempRates,rates,true,true);
%         temperatures(:,trial) = data(i).G(trial).TR.temps;
% 




    x_avg = mean(x,2,'omitnan');
    y_avg = mean(y,2,'omitnan');
    % temperature for each frame
    temperature = data(i).data(trial).occupancy.temp;
    temp = [temperature,x_avg,y_avg];

    % smooth into 1 minute bins
    smoothData = [];
    bin = time_bins*f2m;
    for c = 1:size(temp,2)
        smoothData(:,c) = smooth(temp(:,c),'moving',bin);
    end
    dsData = smoothData(1:bin:end,:); % select one point for each minute

    % get temp color distribution (1/4 degree incrememts)
    tp = getTempTurnPoints(data(i).temp_protocol);
    cMapRange = tp.threshLow:color_bin:tp.threshHigh;
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
        fig = getfig('',true,[484 384]);
        hold on
            scatter(dsData(:,2), dsData(:,3), 5, kolors,'filled')
            for well = 1:4
                C = data(i).data(trial).data.wellcenters(:,well);
                scatter(C(1),C(2),SZ,Color('grey'),'filled')
            end
            scatter(well_C(1),well_C(2),SZ,'y','filled')
            viscircles([arena_C(1),arena_C(2)],r,'Color',foreColor);
            scatter(position(:,1),position(:,2),15,cMap,'filled','o','MarkerEdgeColor',foreColor)
            axis equal square;
            formatFig(fig,true);
            set(gca,'XColor',backColor,'YColor',backColor);
            % set color bar information
            colorData = uint8(round(cMap.*255)); % convert color map to uint8
            colormap(colorData);
            c = colorbar('color',foreColor);
            set(c,'Ticks',[0,1],'TickLabels',[tp.threshLow,tp.threshHigh]);
            c.Label.String = 'Temperature (\circC)';
            c.Label.VerticalAlignment = "bottom";
        title([data(i).T.Genotype{trial} ' ' data(i).T.Date{trial}],'color',foreColor)

        figFolder = [data(i).figDir 'avg position by trial/'];
        if ~exist(figFolder,'dir')
            mkdir(figFolder);
        end
        save_figure(fig,[figFolder data(i).T.Date{trial} ' ' data(i).T.ExperimentID{trial} ' ' data(i).T.Arena{trial}],fig_type,true);
    end
  end
end


% use the temperatures : 8:2:22 for key points to check location of flies
clearvars('-except',initial_vars{:})

n = 10; % number of spatial bins
tt = 1; 
buffer = 0.5; % temperature buffer around target temperature

% 1) select all frames within a given temperature bin (from temp rate
% analysis in step 3.1) 
% 2) align the grid to a uniform heading
% 3) normalize probabiltiy for each of the 2d space bins

tt = 24.5; 

clearvars('-except',initial_vars{:})
OG_Orientation = datetime('10.20.2023','InputFormat','MM.dd.yyyy');

autoSave = true; %autosave the resulting figures
[foreColor,backColor] = formattingColors(blkbgd);


for i = 1:num.exp
    fig = figure; hold on
    for trial = 1:num.trial(i)
        
        % Find the position data for all flies for the selected timepoint


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
        
        trial_date = datetime(data(i).T.Date{trial},'InputFormat','MM.dd.yyyy');

        if trial_date > OG_Orientation %new camera orientation
            % Rotate to correct orientation
            switch wellLoc
                case 1
                    plotData(:,1) = X;
                    plotData(:,2) = Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
                case 2
                    plotData(:,1) = Y;
                    plotData(:,2) = -X;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 3
                    plotData(:,1) = X;
                    plotData(:,2) = -Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 4
                    plotData(:,1) = -Y;
                    plotData(:,2) = X;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
            end
        else % Rotate to correct orientation with older camera arrangement
            
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

        end

        % PLOT
        scatter(WELLS(1:4,1),WELLS(1:4,2),SZ,foreColor,'filled')
        scatter(WELLS(wellLoc,1),WELLS(wellLoc,2),SZ,'green','filled')
        scatter(plotData(:,1),plotData(:,2),15,kolor,'filled')

    end

    viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color',grouped(i).color);
%         viscircles([WELLS(5,1),WELLS(5,2)],data(i).data(trial).data.r,'Color',foreColor);
    axis square; axis equal;
    formatFig(fig,blkbgd);
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
    save_figure(fig,[saveDir expGroup grouped(i).name ' COM position'],fig_type,autoSave);

end

% % TIDY WAY TO MAKE HEATMAPS 
% idx = 1;
% for rr = 1:length(rateList)
%   for tt = 1:length(tempList)   
%     temp = tempList(tt); % temp in C
%     rate = find(rateList(rr)==G(trial).TR.rates);
%     HM = zeros(n);
%     for trial = 1:ntrials
%         % frame location selection
%         tLoc = G(trial).TR.data(:,1)>=temp-buffer & G(trial).TR.data(:,1)<=temp+buffer;
%         rLoc = G(trial).TR.rateIdx==rate;
%         frames = tLoc & rLoc;
% 
%         % pull the position data for these frames and concatenate into a large
%         % structure for this temp rate and temp location
%         x = data(trial).data.x_loc(frames,:);
%         y = data(trial).data.y_loc(frames,:);
%         X = reshape(x,numel(x),1);
%         Y = reshape(y,numel(y),1);
%         pos = [X,Y];
% 
%         % get the 'square' units for partitioning space
%         C = data(trial).data.centre;
%         r = data(trial).data.r;
%         x_edge = linspace(C(1)-r,C(1)+r,n);
%         y_edge = linspace(C(2)-r,C(2)+r,n);
% 
%         % find x and y that are within each 'box'
%         xInd = discretize(pos(:,1),x_edge);
%         yInd = discretize(pos(:,2),y_edge);
% 
%         % find the number of flies within each spatial bin:
%         for row = 1:n
%             for col = 1:n
%                 nflies(row,col) = sum(yInd==row & xInd==col);
%             end
%         end
% 
%         % Rotate the matrix if needed to align to a well position of '2'
%         k = T.foodLoc(trial)-2;
%         B = rot90(nflies,k);
% 
% %         fig = figure; imagesc(B);
% %         uiwait(fig)
% 
%         % Save matrix to the structure:
%         HM = HM + B;
%     end
%     plotData(rr,tt).HM = HM;
%     cmaps(idx,1) = min(min(HM));
%     cmaps(idx,2) = max(max(HM));
%     idx = idx + 1;
%   end
% end
% 
% % figures   
% idx = 0;
% fig = figure; set(fig, 'color', 'k', 'position',  [547 304 992 555]); %[32 70 1318 435]
% for rr = 1:length(rateList)
%   for tt = 1:length(tempList)
%     idx = idx+1;
%     subplot(length(rateList),length(tempList),idx)
%         imagesc(plotData(rr,tt).HM)
%         ax = gca;
%         set(ax, 'xscale', 'log', 'yscale', 'log')
%         axis tight
%         set(ax, 'XColor', 'k','YColor', 'k', 'XTick', [],'YTick', []);
%         c = colorbar;
%         c.Label.String = 'Number of flies';
%         c.Label.Color = 'w';
%         c.Color = 'w';
%         title([num2str(tempList(tt)) ' \circC at ' num2str(rateList(rr)) ' \circC/min'],'color', 'w')
%         caxis([min(cmaps(:,1)) max(cmaps(:,2))]) 
%         axis square
%   end
% end
% 
% save_figure(fig, [figDir 'Four temp hysteresis position heatmaps'], '-png');
% clearvars('-except',vars{:})