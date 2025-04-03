
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
    tp = getTempTurnPoints(data(i).temp_protocol);
    fps = tp.fps;

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
    bin = time_bins*(fps*60); % convert by frames per minute
    for c = 1:size(temp,2)
        smoothData(:,c) = smooth(temp(:,c),'moving',bin);
    end
    dsData = smoothData(1:bin:end,:); % select one point for each minute

    % get temp color distribution (1/4 degree incrememts)
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

    bin = time_bins*(tp.fps*60);
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


% %% ANALYSIS: normalize fly position within arena 
% clearvars('-except',initial_vars{:})
% OG_Orientation = datetime('10.20.2023','InputFormat','MM.dd.yyyy');
% 
% for i = 1:num.exp
%   % BINNED
%   [tempRates,decreasing,increasing,temperatures] = deal([]);
%   % Pull temperature and rate information for the temp protocol
%   temp_rates = data(i).G(1).TR.rates;
%   temp_list = data(i).G(1).TR.temps;
%   nrates = data(i).G(1).TR.nRates;
%   ntemps = data(i).G(1).TR.nTemps;
%   rate_idx = data(i).G(1).TR.rateIdx;
%   temp_idx = data(i).G(1).TR.tempIdx;
%   loc_mat = struct;
% 
%   % find frame index for each temp bin & rate of change
%   for tt = 1:ntemps
%     for rr = 1:nrates
%         rateAligned = rate_idx==rr;
%         tempAligned = temp_idx==tt;
%         loc = find(rateAligned & tempAligned);
%         if isempty(loc)
%             loc = nan;
%         end
%         loc_mat(rr,tt).frames = loc;
%         loc_mat(rr,tt).x = []; %set empty space for appending
%         loc_mat(rr,tt).y = [];
%     end
%   end
%   wellXY = [];
%   for trial = 1:num.trial(i)
%     trial_date = datetime(data(i).T.Date{trial},'InputFormat','MM.dd.yyyy');
%     % get arena information
%     well_loc = data(i).T.foodLoc(trial);
%     wells = data(i).data(trial).data.wellcenters;
%     % find offset to make the food well the origin
%     x_offset = wells(1,well_loc);
%     y_offset = wells(2,well_loc);
%     wells_x = wells(1,:)-x_offset;
%     wells_y = wells(2,:)-y_offset;
% 
%     X = data(i).data(trial).data.x_loc;
%     Y = data(i).data(trial).data.y_loc;
%     X = X-x_offset;
%     Y = Y-y_offset;
%     % save the position normalized data into the grouped structure or
%     % something
% 
%     [WELLS,x_data,y_data] = deal([]);
%     if trial_date > OG_Orientation %new camera orientation
%         % Rotate to correct orientation
%         switch well_loc
%             case 1
%                 x_data = X;
%                 y_data = Y;
%                 WELLS(:,1) = wells_x;
%                 WELLS(:,2) = wells_y;
%             case 2
%                 x_data = Y;
%                 y_data = -X;
%                 WELLS(:,1) = wells_y;
%                 WELLS(:,2) = -wells_x;
%             case 3
%                 x_data = X;
%                 y_data = -Y;
%                 WELLS(:,1) = wells_x;
%                 WELLS(:,2) = -wells_y;
%             case 4
%                 x_data = -Y;
%                 y_data = X;
%                 WELLS(:,1) = -wells_y;
%                 WELLS(:,2) = wells_x;
%         end
%     else % Rotate to correct orientation with older camera arrangement            
%         switch well_loc
%             case 1
%                 x_data = Y;
%                 y_data = -X;
%                 WELLS(:,1) = wells_y;
%                 WELLS(:,2) = -wells_x;
%             case 2
%                 x_data = X;
%                 y_data = -Y;
%                 WELLS(:,1) = wells_x;
%                 WELLS(:,2) = -wells_y;
%             case 3
%                 x_data = -Y;
%                 y_data = X;
%                 WELLS(:,1) = -wells_y;
%                 WELLS(:,2) = wells_x;
%             case 4
%                 x_data = X;
%                 y_data = Y;
%                 WELLS(:,1) = wells_x;
%                 WELLS(:,2) = wells_y;
%         end
% 
%     end
% 
%     wellXY.x(:,trial) = WELLS(:,1);
%     wellXY.y(:,trial) = WELLS(:,2);
% 
%     % For each temp and rate, pool all (now normalized) fly positions
%     for tt = 1:ntemps
%       for rr = 1:nrates
%         frame_idx = loc_mat(rr,tt).frames;
%         % shift positions for food well as origin
%         if ~isnan(frame_idx)
%             x = x_data(frame_idx,:);
%             y = y_data(frame_idx,:);
% 
%             % save data into structure:
%             loc_mat(rr,tt).data(trial).pos = [x(:),y(:)];
%             loc_mat(rr,tt).x = [loc_mat(rr,tt).x; x(:)];
%             loc_mat(rr,tt).y = [loc_mat(rr,tt).y; y(:)];
%         else
%             loc_mat(rr,tt).data(trial).pos = [nan nan];
%             loc_mat(rr,tt).x = nan;
%             loc_mat(rr,tt).y = nan;
%         end
%       end
%     end
%   end
% 
%   % save trial data into broad structure
%   grouped(i).position.loc = loc_mat;
%   grouped(i).position.temp_rates = temp_rates;
%   grouped(i).position.temp_list = temp_list;
%   grouped(i).position.well_pos = wellXY;
% end

%% FIGURE: plot heatmap of fly position within arena
clearvars('-except',initial_vars{:})
save_path = createFolder([saveDir 'COM/']);
autoSave = true;

fig_type = '-pdf';
blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);

% Find the occupancy for each bin:
r = data(1).data(1).data.r; %pixel radius of the arena
n = 26; % number of spatial bins
autoLim = false;
axis_limits = [0, 0.01];

expList = 1:num.exp;

% Set Temperature
for temp = [15:4:35] %[17,18, 20,30,32]%16:2:35 %(17:2:25)

plotData = [];
max_occ = [];
% GROUP DATA
for i = expList
    % get the 'square' units for partitioning space
    Cx = mean(grouped(i).position.well_pos.x(5,:)); %center X
    Cy = mean(grouped(i).position.well_pos.y(5,:)); %center Y
    x_edge = linspace(Cx-r,Cx+r,n);
    y_edge = linspace(Cy-r,Cy+r,n);

    % determine what a circle around the arena would look like:
    % r needs to be transformed into unit square space...
    square_unit = mean(diff(x_edge)); % pixel size for one bin
    circ_r = r/square_unit; % arena radius in bin size
    circ_X = discretize(Cx, x_edge);
    circ_Y = discretize(Cy, y_edge);

    % find the temp bin that corresponds to the selected temp:
    [~,idx] = min(abs(grouped(i).position.temp_list-temp));
    nRates = length(grouped(i).position.temp_rates);
    nflies = [];
    for rr = 1:nRates
        % find x and y that are within each 'box'
        x = grouped(i).position.loc(rr,idx).x;
        y = grouped(i).position.loc(rr,idx).y;
        nanLoc = isnan(x)| isnan(y);
        x(nanLoc) = [];
        y(nanLoc) = [];

        xInd = discretize(x,x_edge);
        yInd = discretize(y,y_edge);
    
        % find the number of flies within each spatial bin:
        for row = 1:n
            for col = 1:n
                nflies(row,col) = sum(yInd==row & xInd==col);
            end
        end
        % turn to prob and not direct occupancy
        plotData(i,rr).data = nflies./sum(sum(nflies));
        
        max_occ = max([max_occ,max(max(plotData(i,rr).data))]);

        % Find the wells within the binned space
        % wellX = (grouped(i).position.well_pos.x(1:4,:)); 
        % wellY = (grouped(i).position.well_pos.y(1:4,:));
        % wellX = wellX(:);
        % wellY = wellY(:);
        xInd = discretize(0,x_edge);
        yInd = discretize(0,y_edge);

        plotData(i,rr).wells = [xInd,yInd];
    end
end

disp(['Max occupancy: ' num2str(max_occ)])

% PLOT 
fig_W = 20 + (400*nRates);

for i = expList
    fig = getfig('',false,[fig_W, 340]); 
    for rr = 1:nRates
        subplot(1,nRates,rr)
        hold on
        imagesc(plotData(i,rr).data); hold on
        scatter(plotData(i,rr).wells(:,1),plotData(i,rr).wells(:,2),10,'r','filled')
        axis tight;
        axis square;
        % h = drawcircle('Center',[circ_X,circ_Y],'Radius',circ_r,'StripeColor',foreColor);
        v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
    end
    formatFig(fig, blkbgd,[1,nRates]);
    for rr = 1:nRates
        subplot(1,nRates,rr)
        set(gca,'XColor',backColor,'Ycolor',backColor,'XTick', [],'YTick', [])
        t_str = [num2str(grouped(i).position.temp_rates(rr)) '\circC/min | ' num2str(temp) '\circC'];
        title({grouped(i).name; t_str},'color',foreColor,'fontsize', 12)
        
        % set(gca,'ColorScale','log')

        c = colorbar;
        c.Label.String = 'Occupancy Probability';
        c.Label.Color = foreColor;
        c.Color = foreColor;
        if autoLim
            clim([0,max_occ]) 
        else
            clim(axis_limits)
        end
    end
    colormap(flipud(gray))
    % save the figure to a folder specific to that cohort?
    save_figure(fig,[save_path grouped(i).name ' ' num2str(temp) ' deg'], fig_type,autoSave,true);
end

end



% make mask for arena sizing...

%% FIGURE (TEMP HOLDS ONLY): 2D spatial occupancy in arena histogram
% TODO 4/3 get this to work -- just saved all aligned position data to
% grouped.position.trial.x/y so use that
clearvars('-except',initial_vars{:})
save_path = createFolder([saveDir 'COM/']);
autoSave = true;

fig_type = '-pdf';
blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);

% Find the occupancy for each bin:
r = data(1).data(1).data.r; %pixel radius of the arena
n = 26; % number of spatial bins
autoLim = false;
axis_limits = [0, 0.01];

expList = 1:2:num.exp; % take the caviar trials first

% % Set Temperature
% for temp = [15:4:35] %[17,18, 20,30,32]%16:2:35 %(17:2:25)

plotData = [];
max_occ = [];
% GROUP DATA
for i = expList
    % get the 'square' units for partitioning space  (aka grid lines for subdividing space) 
    Cx = mean(grouped(i).position.well_pos.x(5,:)); %center X
    Cy = mean(grouped(i).position.well_pos.y(5,:)); %center Y
    x_edge = linspace(Cx-r,Cx+r,n);
    y_edge = linspace(Cy-r,Cy+r,n);

    % determine what a circle around the arena would look like:
    % r needs to be transformed into unit square space...
    square_unit = mean(diff(x_edge)); % pixel size for one bin
    circ_r = r/square_unit; % arena radius in bin size
    circ_X = discretize(Cx, x_edge);
    circ_Y = discretize(Cy, y_edge);

    % find the temp bin that corresponds to the selected temp:
    % [~,idx] = min(abs(grouped(i).position.temp_list-temp));
    % nRates = length(grouped(i).position.temp_rates);
    % nflies = [];
    % for rr = 1:nRates
        % find x and y that are within each 'box'
        x = data(i).data(trial).data.x_loc; % x locations for the entire experiment
        y = data(i).data(trial).data.y_loc; % x locations for the entire experiment

        % x = grouped(i).position.loc(rr,idx).x;
        % y = grouped(i).position.loc(rr,idx).y;
        nanLoc = isnan(x)| isnan(y);
        x(nanLoc) = [];
        y(nanLoc) = [];

        xInd = discretize(x,x_edge);
        yInd = discretize(y,y_edge);
    
        % find the number of flies within each spatial bin:
        for row = 1:n
            for col = 1:n
                nflies(row,col) = sum(yInd==row & xInd==col);
            end
        end
        % turn to prob and not direct occupancy
        plotData(i,rr).data = nflies./sum(sum(nflies));
        
        max_occ = max([max_occ,max(max(plotData(i,rr).data))]);

        % Find the wells within the binned space
        % wellX = (grouped(i).position.well_pos.x(1:4,:)); 
        % wellY = (grouped(i).position.well_pos.y(1:4,:));
        % wellX = wellX(:);
        % wellY = wellY(:);
        xInd = discretize(0,x_edge);
        yInd = discretize(0,y_edge);

        plotData(i,rr).wells = [xInd,yInd];
    end
end

disp(['Max occupancy: ' num2str(max_occ)])

% PLOT 
fig_W = 20 + (400*nRates);

for i = expList
    fig = getfig('',false,[fig_W, 340]); 
    for rr = 1:nRates
        subplot(1,nRates,rr)
        hold on
        imagesc(plotData(i,rr).data); hold on
        scatter(plotData(i,rr).wells(:,1),plotData(i,rr).wells(:,2),10,'r','filled')
        axis tight;
        axis square;
        % h = drawcircle('Center',[circ_X,circ_Y],'Radius',circ_r,'StripeColor',foreColor);
        v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
    end
    formatFig(fig, blkbgd,[1,nRates]);
    for rr = 1:nRates
        subplot(1,nRates,rr)
        set(gca,'XColor',backColor,'Ycolor',backColor,'XTick', [],'YTick', [])
        t_str = [num2str(grouped(i).position.temp_rates(rr)) '\circC/min | ' num2str(temp) '\circC'];
        title({grouped(i).name; t_str},'color',foreColor,'fontsize', 12)
        
        % set(gca,'ColorScale','log')

        c = colorbar;
        c.Label.String = 'Occupancy Probability';
        c.Label.Color = foreColor;
        c.Color = foreColor;
        if autoLim
            clim([0,max_occ]) 
        else
            clim(axis_limits)
        end
    end
    colormap(flipud(gray))
    % save the figure to a folder specific to that cohort?
    save_figure(fig,[save_path grouped(i).name ' ' num2str(temp) ' deg'], fig_type,autoSave,true);
end

end






 