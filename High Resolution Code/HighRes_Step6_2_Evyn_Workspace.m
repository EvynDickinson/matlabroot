


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






%% TODO : suppose also need to look at the total number of encounters, 
% since flies are also chasing and seeking (but that should show up as an
% increase in the proportion of encounters to courtship rate as well





%%



fig = getfig('',1); 
[p, yAvg, nAvg] = deal([]);
for trial = 1:num.trials
  subplot(r,c,trial); hold on
    enc_end = encounters(trial).locs(:,2); % end frame of an encounter bout
    enc_start = encounters(trial).locs(:,1); % start frame of an encounter bout
    % time between end of last encounter and start of this encounter
    frames_since_enc = enc_start(2:end) - enc_end(1:end-1); % in frames
    time_since_enc = frames_since_enc./fps;  % converted to seconds
    
    % quick scatter plot of 'yes' vs 'no' 
    loc_yes = [logical(encounters(trial).locs(2:end,3))]; % logical of courtship attempts for a given encounter
    loc_no = ~loc_yes;
    plot_yes = time_since_enc(loc_yes);
    plot_no = time_since_enc(loc_no);

    % attempted courtship
    scatter(ones(size(plot_yes)), plot_yes, SZ, foreColor, 'filled',...
        'MarkerFaceAlpha', FA, 'xjitter', 'density')
    avg = mean(plot_yes, 'omitnan');
    yAvg(trial) = avg;
    plot([1-buff, 1+buff], [avg,avg],'color', yColor,...
        'linewidth', LW)
    % did not attempt courtship
    scatter(2*ones(size(plot_no)), plot_no, SZ, foreColor, 'filled',...
        'MarkerFaceAlpha', FA, 'xjitter', 'density')
    avg = mean(plot_no, 'omitnan');
    nAvg(trial) = avg;
    plot([2-buff, 2+buff], [avg,avg],'color', nColor,...
        'linewidth', LW)

    set(gca, 'yscale', 'log')
    set(gca, 'XTick',1:2, 'xticklabel', {'Y','N'})

    % t-test of duration difference 
    [~, p(trial)] = ttest2(plot_yes, plot_no); % welchs t-test (unpaired)

end
formatFig(fig, blkbgd, [r,c,]);
matchAxis(fig,true);

edge_idx = 1:c:r*c;
bottom_idx = ((r-1)*c)+1:r*c;

% correct for MC with bonferonnis :
hC = p<=sig_alpha; % corrected significance
h = p<0.05; % uncorrected significance
for trial = 1:num.trials
    subplot(r,c,trial); 
    if hC(trial)
        title('*', 'color', foreColor,'fontsize', 25)
    elseif h(trial)
        title('* nc', 'color', Color('grey'),'fontsize', 20)
    end
    if ~any(trial==edge_idx)
        set(gca, 'ycolor', 'none')
    end
end





%%
    

types = {'WT', 'WS', 'CT', 'CS'};
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
nTypes = length(types);

tb = data.tempbin;
[tot_encounters,tot_courtship,court_rate, encounter_rate] = deal(nan([num.trials,4]));

for i = 1:nTypes % for each of the temp regime types
    roi = find(tb.(types{i})); % find the frames within the region
    time_period = length(roi); % how many frames possible could there be for encounters (to normalize in case there are difference in the time periods of the diff temp regimess
    % disp(time_period)
    % which encounter locations are within the temp regime
    for trial = 1:num.trials
        % encounters within this temperature regime
        encounter_locs = ismember(encounters(trial).locs(:,1),roi);
        tot_encounters(trial,i) = sum(encounter_locs);
        tot_courtship(trial,i) = sum(encounters(trial).locs(encounter_locs,3));
        court_rate(trial,i) = (tot_courtship(trial,i)/tot_encounters(trial,i))*100;
    end
    encounter_rate(:,i) = (tot_encounters(:,i)/time_period).*100;
end

% PLOT THE PERCENTAGE OF FRAMES WITH AN ENCOUNTER PER TEMP REGIME
y_avg = median(encounter_rate,1,'omitnan');
y_sem = std(encounter_rate,0,1,'omitnan')./sqrt(num.trials);
% plot the encounter rate data (bar graph)
% Create a bar plot for the average CI data across temperature regimes
sz = 75;
fig = getfig('', 1, [436 620]);
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:nTypes, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:nTypes,[num.trials, 1]);
    scatter(x(:), encounter_rate(:), sz, Color('grey'), 'filled', ...
        'xjitter','density','xjitterwidth', 0.3, 'MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels,'XTickLabelRotation',30);
    ylabel('Rate of encounters (% of total time)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
    % ylim([-3,70])
save_figure(gcf, [figDir 'rate of encounters across temp regimes'], fig_type);

% run simple stats:
% Warm threat vs warm safe
[h,p,~,stats] = ttest(encounter_rate(:,1), encounter_rate(:,2)); 
fprintf('encounter rate warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% cool threat vs cool safe
[h,p,~,stats] = ttest(encounter_rate(:,3), encounter_rate(:,4)); 
fprintf('encounter rate cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

%%

%% Where are flies sleeping in the arena? Does this change over temperature?
clearvars('-except',initial_var{:})

[title_str,pName,y_dir,y_lab,nullD,scaler,dType,dir_end,sexSep,ylimits] = PlotParamSelectionHR(true);


loc = data.sleep;
asleep = data.dist2food(loc);
awake = data.dist2food(~loc);
binedges = 0:ceil(max(max(max(data.dist2food))));
sleep_color = Color('vaporwaveblue');
awake_color = Color('grey');

fig = getfig('',1); hold on
    yyaxis left 
    histogram(awake, binedges, 'FaceColor', awake_color)
    yyaxis right 
    histogram(asleep, binedges, 'FaceColor', sleep_color,...
        'FaceAlpha',0.8)
    formatFig(fig, blkbgd);
    xlabel('distance to food (mm)')
    yyaxis left 
    set(gca, 'ycolor', awake_color)
    ylabel('awake fly count (#)')
    yyaxis right
    set(gca, 'ycolor', sleep_color)
    ylabel('asleep fly count (#)')
save_figure(fig, [figDir 'sleeping distance to food histogram']);


%% 
% fly sleep based on region location not distance
loc = data.sleep;

% compare the norm distibution to % of flies asleep for each region
% food quad


foodQuad = data.foodQuad(loc);
awake = data.dist2food(~loc);
binedges = 0:ceil(max(max(max(data.dist2food))));
sleep_color = Color('vaporwaveblue');
awake_color = Color('grey');

fig = getfig('',1); hold on
    yyaxis left 
    histogram(awake, binedges, 'FaceColor', awake_color)
    yyaxis right 
    histogram(asleep, binedges, 'FaceColor', sleep_color,...
        'FaceAlpha',0.8)
    formatFig(fig, blkbgd);
    xlabel('distance to food (mm)')
    yyaxis left 
    set(gca, 'ycolor', awake_color)
    ylabel('awake fly count (#)')
    yyaxis right
    set(gca, 'ycolor', sleep_color)
    ylabel('asleep fly count (#)')
save_figure(fig, [figDir 'sleeping distance to food histogram']);

%% Figure out why the speedTest frames are not aligning across computers

clearvars('-except',initial_var{:})

keyFramesMAC = keyFrames;
load('keyframes.mat')

% compare the two keyframes data (these should be identical) 
errList = [];
for trial = 1:num.trials
    % key frame comparisons
    a = keyFramesMAC(trial).frame_pairs(:)';
    b = keyFrames(trial).frame_pairs(:)';
    errList(trial,1) = sum(~ismember(a,b));
    errList(trial,2) = sum(~ismember(b,a));
end


%% 

%% ANALYSIS & FIGURES: identify likely frame swap locations
% TODO: 2.25.26 convert this to something that can be used in 
% the main data analysis pipeline
initial_var = add_var(initial_var, 'keyFrames');
clearvars('-except',initial_var{:})

speed_thresh = 35; %mm/s speed threshold
skip_threshold = 3; % how many frames for a confident swap pair in time

% data structure holding the frame data
keyFrames = [];  

for trial = 1:num.trials
    % MALE speed data for current trial 
    mspeed = squeeze(data.speed(:,M,trial));
    mspeed_loc = mspeed>=speed_thresh; % frames with above threshold speed
    mspeed_frames = find(mspeed_loc); % pull frames above the limit speed
    keyFrames(trial).mspeed_frames = mspeed_frames;
    
    % FEMALE speed data for current trial 
    fspeed = squeeze(data.speed(:,F,trial));
    fspeed_loc = fspeed>=speed_thresh; % frames with above threshold speed
    fspeed_frames = find(fspeed_loc); % pull frames above the limit speed
    keyFrames(trial).fspeed_frames = fspeed_frames;
    
    % Possible swap locations: (based on high speed for both flies)
    double_speed = squeeze(data.speed(:,:,trial)); % speed in both sexes
    double_speed_loc = double_speed>=speed_thresh;
    swap_ConfidentFrames = find(sum(double_speed_loc,2)==2); % where are both M&F speeding
    
    % Possible swap locations if one of the sexes does not have a labeled skeleton
    a = fspeed_loc & isnan(mspeed);
    b = mspeed_loc & isnan(fspeed);
    possible_swapFrames = find(a | b);
    if ~isempty(possible_swapFrames)
        allSwaps = sort([swap_ConfidentFrames; possible_swapFrames]); % all frames with double speed MF or single+nan
        keyFrames(trial).swap_LikelyFrames = allSwaps; % all frames with double speed M F
    else % no changes due to new pairings, so the same frames as confident
        keyFrames(trial).swap_LikelyFrames = swap_ConfidentFrames;
    end
end

% PAIR LIKELY SWAP FRAMES FOR EACH OF THE TRIALS
nNeighbors = 5; % how many preceeding frames to look at for close distance?
for trial = 1:num.trials
    a = keyFrames(trial).swap_LikelyFrames;% current list of possible frames for pairs
    % Find differences between consecutive frames
    frame_diffs = diff(a);
    % Find where frames are close together (potential pairs)
    is_close = frame_diffs <= skip_threshold;
    % Extract pairs
    pairs = [];
    ii = 1; % index location within frame list
    while ii <= length(is_close)
        if is_close(ii)
            % Found a pair: frames a(i) and a(i+1)
            pairs = [pairs; a(ii), a(ii+1)];
            
            % Skip the next position to avoid overlapping pairs
            % (e.g., if frames 1,2,3 are all close, we want pairs [1,2] and not [2,3])
            ii = ii + 2;
        else
            ii = ii + 1;
        end
    end
    
    % Second filter: body size and position relative to past trajectory
    mPos = squeeze(fly(trial).m.pos(:,body.center,:));
    fPos = squeeze(fly(trial).f.pos(:,body.center,:));

    % for each of the likely speed frame pairs
    likely_pairs = pairs;
    rois = likely_pairs(:,1) + ((-nNeighbors-1): -1);
   
    % check location alignment for each swap pair
    likelyswitch = false([size(likely_pairs, 1), 1]);
    for ii = 1:size(likely_pairs,1)
        
        % frame locations for pre and during swaps
        pre_roi = rois(ii,:);
        dur_roi = likely_pairs(ii,1);
        
        % distance of each during swap roi to pre-swap roi relative
        x = mPos(dur_roi,1) - mPos(pre_roi,1); 
        y = mPos(dur_roi,2) - mPos(pre_roi,2);
        dM = mean(hypot(x,y),'omitnan');  % male track distance to assigned male (distance to male)
        x = mPos(dur_roi,1) - fPos(pre_roi,1); 
        y = mPos(dur_roi,2) - fPos(pre_roi,2);
        dF = mean(hypot(x,y),'omitnan');  % male track distance to assigned female (distance to female)

        % determine if likely switch
        likelyswitch(ii) = dM>=dF; % female would be closer than correct male
    end
    % extract the likely swap locations
    if any(likelyswitch)
        likely_pairs = pairs(likelyswitch,:); % pull only the likely pairs
        keyFrames(trial).swap_pairs = likely_pairs;
    end
end
keyFrames = rmfield(keyFrames, 'swap_LikelyFrames'); % remove excess field

% Quick look at the jump vs swap stats
m_frames = arrayfun(@(x) size(x.mspeed_frames, 1), keyFrames);
f_frames = arrayfun(@(x) size(x.fspeed_frames, 1), keyFrames);
pair_frames = arrayfun(@(x) numel(x.swap_pairs), keyFrames);

% relative percentages
tot_frames = size(fly(1).T,1);
m_framesP = mean((m_frames/tot_frames)*100);
f_framesP = mean((f_frames/tot_frames)*100);
swap_percentage =  mean((pair_frames ./ tot_frames)*100,'omitnan');
fprintf('\n%2.3g  percent of total frames are high speed \n', (m_framesP + f_framesP))
fprintf('\n%2.5g  percent of total frames are confident swaps',swap_percentage);

% FIX THE SWAPS:


