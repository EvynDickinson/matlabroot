


%% plot the avg and err of a variable for both flies along with the temperature 
clearvars('-except',initial_var{:})
figFolder = createFolder([figDir 'timecourse/']);
[foreColor, ~] = formattingColors(blkbgd); %get background colors
LW = 1.5; %linewidth for plotting
plot_err = true; % plot the error on the graph?
fig_type = '-png';
sSpan = 15*fly(1).fps;

opt_list = {'distance to food', 'inter-fly-distance', 'eccentricity', 'food quadrant', 'food circle', 'turning','sleep','outer ring','courtship index', 'flies on food'}; 
list_idx = listdlg('promptstring','Select the variable to plot', 'ListString',opt_list,'ListSize',[200,180]);
if isempty(list_idx)
        disp('No variable selected')
        return
end  

% defaults 
multiplier = 1;
sem_scale = 1; %(1/num.trials);
axis_dir = 'normal';

% pre-load variable specific data format
switch opt_list{list_idx}
    case 'distance to food'
        y_label = [opt_list{list_idx} ' (mm)'];
        axis_dir = 'reverse';
        d_type = 2; % M F separate
        d_name = 'dist2food';
    case 'inter-fly-distance'
        y_label = [opt_list{list_idx} ' (mm)'];
        d_type = 1; % combined MF
        d_name = 'IFD';
    case 'courtship index'
        y_label = opt_list{list_idx};
        d_type = 1; % combined MF
        d_name = 'CI';
    case 'flies on food'
        y_label = opt_list{list_idx};
        d_type = 2; % M F separate
        d_name = 'FlyOnFood';
        multiplier = 100;
    case 'eccentricity'
        y_label = [opt_list{list_idx} ' (mm)'];
        d_type = 2; % M F separate
        d_name = opt_list{list_idx};
    case 'food quadrant'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = 'foodQuad';
        multiplier = 100;
    case 'food circle'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = 'foodcircle';
        multiplier = 100;
    case 'outer ring'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = 'OutterRing';
        multiplier = 100;
    case 'turning'
        y_label = [opt_list{list_idx} ' (mm/s)'];
        d_type = 2; % M F separate
        d_name = opt_list{list_idx};
    case 'sleep'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = opt_list{list_idx};
        multiplier = 100;
end

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; % variable to examine timecourse
sb(3).idx = 3:c:r*c; % temp binned variable alignment

fig = getfig('',1);
% temperature time course
subplot(r,c,sb(1).idx); hold on
x = data.time;
y = data.temp;
plot(x,y,'color', foreColor, 'linewidth',LW)
% TODO: update this code to load and parse the data correctly for diff
% sized data structures (e.g., 2 flies vs single column)  
raw_data = data.(d_name);

% time course of variable
subplot(r,c,sb(2).idx); hold on
x = data.time;
switch d_type
    case 1 % single value for M & F
            y_all = raw_data;
            y_avg = smooth(mean(y_all,2,'omitnan'),sSpan,'moving').*multiplier;
            y_err = smooth(std(y_all, 0,2,'omitnan'),sSpan,'moving').*multiplier.*sem_scale;
            plot_error_fills(plot_err, x, y_avg, y_err, foreColor,fig_type, 0.35);
            plot(x,y_avg,'color',foreColor,'linewidth',LW+1)
    case 2 % separate values for M & F 
        for sex = 1:2
            y_all = squeeze(raw_data(:,sex,:));
            y_avg = smooth(mean(y_all,2,'omitnan'),sSpan,'moving').*multiplier;
            y_err = smooth(std(y_all, 0,2,'omitnan'),sSpan,'moving').*multiplier.*sem_scale;
            plot_error_fills(plot_err, x, y_avg, y_err, data.color(sex,:),fig_type, 0.35);
            plot(x,y_avg,'color',data.color(sex,:),'linewidth',LW+1)
        end
end

%temp dependent variable
subplot(r,c,sb(3).idx); hold on
x = data.tempbin.temps;
    for sex = 1:2
        % setup the data for extraction
        if d_type==1
            y_all = raw_data;
            kolor = foreColor;
        else
            y_all = squeeze(raw_data(:,sex,:));
            kolor = data.color(sex,:);
        end
        if d_type==1 && sex==2
            continue
        end
        %extract data
        [c_avg,w_avg,c_err,w_err] = deal([]);
        for t = 1:length(x)
            % cooling
            loc = data.tempbin.cooling(:,t);
            y_sel = mean(y_all(loc,:),1,'omitnan').*multiplier;
            c_avg(t) = mean(y_sel,'omitnan');
            c_err(t) = std(y_sel, 0,2,'omitnan').*sem_scale;
            % warming
            loc = data.tempbin.warming(:,t);
            y_sel = mean(y_all(loc,:),1,'omitnan').*multiplier;
            w_avg(t) = mean(y_sel,'omitnan');
            w_err(t) = std(y_sel, 0,2,'omitnan').*sem_scale;
            % % cooling
            % loc = data.tempbin.cooling(:,t);
            % y_sel = mean(squeeze(raw_data(loc,sex,:)),1,'omitnan');
            % c_avg(t) = mean(y_sel,'omitnan');
            % c_err(t) = std(y_sel, 0,2,'omitnan');
            % % warming
            % loc = data.tempbin.warming(:,t);
            % y_sel = mean(squeeze(raw_data(loc,sex,:)),1,'omitnan');
            % w_avg(t) = mean(y_sel,'omitnan');
            % w_err(t) = std(y_sel, 0,2,'omitnan');
        end
        % plot data
        plot_error_fills(plot_err, x, c_avg, c_err, kolor,fig_type);
        plot_error_fills(plot_err, x, w_avg, c_err, kolor,fig_type);
        plot(x,c_avg,'color',kolor,'linewidth',LW+1,'LineStyle','--')
        plot(x,w_avg,'color',kolor,'linewidth',LW+1,'LineStyle','-')
    end

% format figure: 
formatFig(fig, blkbgd, [r,c],sb);
subplot(r,c,sb(1).idx);
    set(gca, 'xcolor', 'none')
    ylabel('\circC')
subplot(r,c,sb(2).idx);    
    ylabel(y_label)
    xlabel('time (min)')
    set(gca, 'YDir',axis_dir)
     ylims = ylim;
    if ylims(1)<0
        ylims(1) = 0;
    end
    if ylims(2) >100
        ylims(2) = 100;
    end
    ylim(ylims)
subplot(r,c,sb(3).idx);
    ylabel(y_label)
    xlabel('temp (\circC)')
    set(gca, 'YDir',axis_dir)
     ylims = ylim;
    if ylims(1)<0
        ylims(1) = 0;
    end
    if ylims(2) >100
        ylims(2) = 100;
    end
    ylim(ylims)

save_figure(fig, [figFolder opt_list{list_idx}],fig_type);


%%
%% Total count variables: plot the avg and err of a variable for both flies along with the temperature 
clearvars('-except',initial_var{:})
figFolder = createFolder([figDir 'timecourse/']);
[foreColor, ~] = formattingColors(blkbgd); %get background colors
LW = 1.5; %linewidth for plotting
plot_err = true; % plot the error on the graph?
fig_type = '-png';
sSpan = 1;

opt_list = {'courtship index', 'flies on food','sleep'}; 
list_idx = listdlg('promptstring','Select the variable to plot', 'ListString',opt_list,'ListSize',[200,180]);
if isempty(list_idx)
        disp('No variable selected')
        return
end  

% defaults 
multiplier = 1;
axis_dir = 'normal';

% pre-load variable specific data format
switch opt_list{list_idx}
    case 'distance to food'
        y_label = [opt_list{list_idx} ' (mm)'];
        axis_dir = 'reverse';
        d_type = 2; % M F separate
        d_name = 'dist2food';
    case 'inter-fly-distance'
        y_label = [opt_list{list_idx} ' (mm)'];
        d_type = 1; % combined MF
        d_name = 'IFD';
    case 'courtship index'
        y_label = opt_list{list_idx};
        d_type = 1; % combined MF
        d_name = 'CI';
    case 'flies on food'
        y_label = opt_list{list_idx};
        d_type = 2; % M F separate
        d_name = 'FlyOnFood';
    case 'eccentricity'
        y_label = [opt_list{list_idx} ' (mm)'];
        d_type = 2; % M F separate
        d_name = opt_list{list_idx};
    case 'food quadrant'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = 'foodQuad';
        multiplier = 100;
    case 'food circle'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = 'foodcircle';
        multiplier = 100;
    case 'outer ring'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = 'OutterRing';
        multiplier = 100;
    case 'turning'
        y_label = [opt_list{list_idx} ' (mm/s)'];
        d_type = 2; % M F separate
        d_name = opt_list{list_idx};
    case 'sleep'
        y_label = [opt_list{list_idx} ' (%)'];
        d_type = 2; % M F separate
        d_name = opt_list{list_idx};
end

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; % variable to examine timecourse
sb(3).idx = 3:c:r*c; % temp binned variable alignment

fig = getfig('',1);
% temperature time course
subplot(r,c,sb(1).idx); hold on
x = data.time;
y = data.temp;
plot(x,y,'color', foreColor, 'linewidth',LW)
% TODO: update this code to load and parse the data correctly for diff
% sized data structures (e.g., 2 flies vs. single column)  
raw_data = data.(d_name);

% time course of variable
subplot(r,c,sb(2).idx); hold on
x = data.time;
switch d_type
    case 1 % single value for M & F
            y_all = raw_data;
            y_avg = smooth(sum(y_all,2,'omitnan'),sSpan,'moving').*multiplier;
            plot(x,y_avg,'color',foreColor,'linewidth',LW+1)
    case 2 % separate values for M & F 
        for sex = 1:2
            y_all = squeeze(raw_data(:,sex,:));
            y_avg = smooth(sum(y_all,2,'omitnan'),sSpan,'moving').*multiplier;
            plot(x,y_avg,'color',data.color(sex,:),'linewidth',LW+1)
        end
end

%temp dependent variable
subplot(r,c,sb(3).idx); hold on
x = data.tempbin.temps;
    for sex = 1:2
        % setup the data for extraction
        if d_type==1
            y_all = raw_data;
            kolor = foreColor;
        else
            y_all = squeeze(raw_data(:,sex,:));
            kolor = data.color(sex,:);
        end
        if d_type==1 && sex==2
            continue
        end
        %extract data
        [c_avg,w_avg,c_err,w_err] = deal([]);
        for t = 1:length(x)
            % cooling
            loc = data.tempbin.cooling(:,t);
            y_sel = sum(y_all(loc,:),1,'omitnan').*multiplier;
            c_avg(t) = sum(y_sel,'omitnan');
            % c_err(t) = std(y_sel, 0,2,'omitnan');
            % warming
            loc = data.tempbin.warming(:,t);
            y_sel = sum(y_all(loc,:),1,'omitnan').*multiplier;
            w_avg(t) = sum(y_sel,'omitnan');
            % w_err(t) = std(y_sel, 0,2,'omitnan');
            % % cooling
            % loc = data.tempbin.cooling(:,t);
            % y_sel = mean(squeeze(raw_data(loc,sex,:)),1,'omitnan');
            % c_avg(t) = mean(y_sel,'omitnan');
            % c_err(t) = std(y_sel, 0,2,'omitnan');
            % % warming
            % loc = data.tempbin.warming(:,t);
            % y_sel = mean(squeeze(raw_data(loc,sex,:)),1,'omitnan');
            % w_avg(t) = mean(y_sel,'omitnan');
            % w_err(t) = std(y_sel, 0,2,'omitnan');
        end
        % plot data
        % plot_error_fills(plot_err, x, c_avg, c_err, kolor,fig_type);
        % plot_error_fills(plot_err, x, w_avg, c_err, kolor,fig_type);
        plot(x,c_avg,'color',kolor,'linewidth',LW+1,'LineStyle','--')
        plot(x,w_avg,'color',kolor,'linewidth',LW+1,'LineStyle','-')
    end

% format figure: 
formatFig(fig, blkbgd, [r,c],sb);
subplot(r,c,sb(1).idx);
    set(gca, 'xcolor', 'none')
    ylabel('\circC')
subplot(r,c,sb(2).idx);    
    ylabel(y_label)
    xlabel('time (min)')
    set(gca, 'YDir',axis_dir)
subplot(r,c,sb(3).idx);
    ylabel(y_label)
    xlabel('temp (\circC)')
    set(gca, 'YDir',axis_dir)

save_figure(fig, [figFolder opt_list{list_idx}],fig_type);

%% 


%% ANALYSIS: Find the quadrants with the highest and lowest occupancy over the exp for null comparison
% as defined by the 7% space around the well (aka using the 'occupancy' measure
clearvars('-except',initial_var{:})

occ_idx = (nan([num.trials,2])); % initialize an empty occupancy index

for trial = 1:num.trials

    
    quad_avg = mean(data(trial).data(trial).data.occ_P,1,'omitnan'); % pull the percent occupancy for each food well
    [~, low_idx] = min(quad_avg); % which quad was lowest occ on avg
    [~, high_idx] = max(quad_avg); % which quad was highest occ on avg
    grouped(trial).occ_idx(trial,1) = low_idx; % register into the index the quad identity for lowest occupancy 
    grouped(trial).occ_idx(trial,2) = high_idx; % register into the index the quad identity for highest occupancy 
end
% add a logical for if this is an empty trial (make it easier for switching code sections later
if strcmp(data(trial).foodNames, 'Movement')
    data(trial).emptytrial = true;
else
    data(trial).emptytrial = false;
end



%% ANALYSIS: normalize fly position within arena to common orientation (food location) 
clearvars('-except',initial_vars{:})
OG_Orientation = datetime('10.20.2023','InputFormat','MM.dd.yyyy'); % camera & lens change accounting

for type = 1:3 % assigned food well, lowest occupied, highest occupied
    switch type
        case 1 % food assigned position 
            field_label = 'position';
        case 2 % lowest occupancy position aligned
            field_label = 'position_low';
        case 3 % highest occupancy position aligned
            field_label = 'position_high';
    end

    for i = 1:num.exp
      % BINNED
      [tempRates,decreasing,increasing,temperatures] = deal([]);
      % Pull temperature and rate information for the temp protocol
      temp_rates = data(i).G(1).TR.rates;
      temp_list = data(i).G(1).TR.temps;
      nrates = data(i).G(1).TR.nRates;
      ntemps = data(i).G(1).TR.nTemps;
      rate_idx = data(i).G(1).TR.rateIdx;
      temp_idx = data(i).G(1).TR.tempIdx;
      loc_mat = struct;
    
      % find frame index for each temp bin & rate of change
      for tt = 1:ntemps
        for rr = 1:nrates
            rateAligned = rate_idx==rr;
            tempAligned = temp_idx==tt;
            loc = find(rateAligned & tempAligned);
            if isempty(loc)
                loc = nan;
            end
            loc_mat(rr,tt).frames = loc;
            loc_mat(rr,tt).x = []; % set empty space for appending
            loc_mat(rr,tt).y = [];
        end
      end
      wellXY = [];
      for trial = 1:num.trial(i)
        trial_date = datetime(data(i).T.Date{trial},'InputFormat','MM.dd.yyyy');
        
        % get arena information
         switch type
             case 1 % food assigned position 
                well_loc = data(i).T.foodLoc(trial);
             case 2 % lowest occupancy position aligned
                well_loc = grouped(i).occ_idx(trial,1);
             case 3 % highest occupancy position aligned
                well_loc = grouped(i).occ_idx(trial,2);
        end
        wells = data(i).data(trial).data.wellcenters;
    
        % % find offset to make the food well the origin
        % x_offset = wells(1,well_loc);
        % y_offset = wells(2,well_loc);
        % wells_x = wells(1,:)-x_offset;
        % wells_y = wells(2,:)-y_offset;
        % 
        % X = data(i).data(trial).data.x_loc;
        % Y = data(i).data(trial).data.y_loc;
        % X = X-x_offset;
        % Y = Y-y_offset;
        % % save the position normalized data into the grouped structure or something
    
        % use the arena center as rotation center rather than the wells
        arena_center_x = mean(wells(1,:));
        arena_center_y = mean(wells(2,:));
        
        % offset all data by arena center
        X = data(i).data(trial).data.x_loc - arena_center_x;
        Y = data(i).data(trial).data.y_loc - arena_center_y;
        
        wells_x = wells(1,:) - arena_center_x;
        wells_y = wells(2,:) - arena_center_y;
    
        [WELLS,x_data,y_data] = deal([]);
        if trial_date > OG_Orientation %new camera orientation
            % Rotate to correct orientation
            switch well_loc
                case 1
                    x_data = X;
                    y_data = Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
                case 2
                    x_data = Y;
                    y_data = -X;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 3
                    x_data = X;
                    y_data = -Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 4
                    x_data = -Y;
                    y_data = X;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
            end
        else % Rotate to correct orientation with older camera arrangement            
            switch well_loc
                case 1
                    x_data = Y;
                    y_data = -X;
                    WELLS(:,1) = wells_y;
                    WELLS(:,2) = -wells_x;
                case 2
                    x_data = X;
                    y_data = -Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = -wells_y;
                case 3
                    x_data = -Y;
                    y_data = X;
                    WELLS(:,1) = -wells_y;
                    WELLS(:,2) = wells_x;
                case 4
                    x_data = X;
                    y_data = Y;
                    WELLS(:,1) = wells_x;
                    WELLS(:,2) = wells_y;
            end
    
        end
        
        wellXY.x(:,trial) = WELLS(:,1);
        wellXY.y(:,trial) = WELLS(:,2);
    
        % For each temp and rate, pool all (now normalized) fly positions
        for tt = 1:ntemps
          for rr = 1:nrates
            frame_idx = loc_mat(rr,tt).frames;
            % shift positions for food well as origin
            if ~isnan(frame_idx)
                x = x_data(frame_idx,:);
                y = y_data(frame_idx,:);
                
                % save data into structure:
                loc_mat(rr,tt).data(trial).pos = [x(:),y(:)];
                loc_mat(rr,tt).x = [loc_mat(rr,tt).x; x(:)];
                loc_mat(rr,tt).y = [loc_mat(rr,tt).y; y(:)];
            else
                loc_mat(rr,tt).data(trial).pos = [nan nan];
                loc_mat(rr,tt).x = nan;
                loc_mat(rr,tt).y = nan;
            end
          end
        end
        grouped(i).(field_label).trial(trial).x = x_data; % this is the full x data for each trial reoriented! 
        grouped(i).(field_label).trial(trial).y = y_data; % this is the full y data for each trial reoriented! 
      end
      
      % save trial data into broad structure
      grouped(i).(field_label).loc = loc_mat;
      grouped(i).(field_label).temp_rates = temp_rates;
      grouped(i).(field_label).temp_list = temp_list; %
      grouped(i).(field_label).well_pos = wellXY; %position of the wells in the newly oriented arena
    end

end












%% FIGURE: Heatmap of the fly positions within the arena across temperature bins
clearvars('-except',initial_var{:})

% Need to rotate and align the fly arena data first...


save_path = createFolder([figDir 'Position Heat Maps/']);
autoSave = true;

% fig_type = '-pdf';
% blkbgd = false;
% [foreColor,backColor] = formattingColors(blkbgd);

% Find the occupancy for each bin:
% r = data(1).data(1).data.r; %pixel radius of the arena

n = 26; % number of spatial bins
autoLim = false;
axis_limits = [0, 0.01];

% Set Temperature
for temp = 15:4:35  % [17,18, 20,30,32] % 16:2:35  % (17:2:25)

    plotData = [];
    max_occ = [];
    % GROUP DATA
    
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

%% DETERMINE WHAT IS UP WITH SPEED DATA: 

figure; 
for ii = 1:num.trials
    subplot(7,2,ii);
    hold on
    histogram(fly(ii).m.speed,'FaceColor','b')
    histogram(fly(ii).f.speed, 'FaceColor', 'm')
end



figure; 
for ii = 1:num.trials
    subplot(7,2,ii);
    hold on
    histogram(data.speed(:,M,ii),'FaceColor','b')
    histogram(data.speed(:,F,ii), 'FaceColor', 'm')
end

figure; 
histogram(data.speed(:,:,:))
xlim()


%% Sort out speed and jump locations???
 
% Inter-fly-distance from the fly's center point
trial = 13;
maxtrial = 3 ; % what peak speed do you want to look at? max = 1, second fastest = 2 etc....

speed_thresh = 100;
x1 = fly(trial).m.pos(:,body.center,1); % x location for male center
y1 = fly(trial).m.pos(:,body.center,2);
x2 = fly(trial).f.pos(:,body.center,1); % x location for female center
y2 = fly(trial).f.pos(:,body.center,2);

D = hypot(diff(x1), diff(y1)); 
test = smooth(D,6,'moving');
a = find(D>=speed_thresh); % where does speed exceed 50mm/s
[maxSpeed, b] = maxk(D,10); % custom selected range that has two high speed frames
if isempty(b)
    disp('no examples of excess speed')
    return
end
B = b(maxtrial)-3:b(maxtrial)+3;
printstring = '\n\t trial had %i frames above %d mm/s\n\t max speed: %i mm/s\n';
titleString = sprintf(printstring, length(a), speed_thresh, round(maxSpeed(maxtrial)));
% fprintf(printstring, length(a), speed_thresh, maxSpeed)
disp(titleString)

fig = figure; hold on
    plot(D,'color', Color('vaporwavepurple'))
    plot(test, 'color', Color('vaporwaveyellow'))
    h_line(50,'vaporwaveblue', '-',2)
    v_line([B(1),B(end)], 'metrored','-',2)
    xlabel('frame number')
    ylabel('fly speed (mm/s)')
    formatFig(fig, blkbgd);
    title(titleString)

% plot out the fly locations for these frames

fig = figure; hold on
for ii = 1:length(B)
    mColor = Color('vaporwaveblue');
    fColor = (Color('vaporwavepink'));
    if B(ii)==b(maxtrial) || B(ii)==b(maxtrial)+1
        mColor = Color('floralblue');
        fColor = Color('floraldarkpink');
    end
    sex = 'm';
    x =  fly(trial).(sex).pos(B(ii),:,1);
    y =  fly(trial).(sex).pos(B(ii),:,2);
    plotFlySkeleton(fig, x,y,mColor,true);
    sex = 'f';
    x =  fly(trial).(sex).pos(B(ii),:,1);
    y =  fly(trial).(sex).pos(B(ii),:,2);
    plotFlySkeleton(fig, x,y,fColor,true);
end
% format
title(titleString)
axis equal
formatFig(fig, blkbgd);
disp('Speeds')
disp(D(B))
disp('Frame Numbers')
disp(B')