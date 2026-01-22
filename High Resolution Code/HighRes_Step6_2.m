
% Generate figures to compare different behaviors and frequencies of
% behaviors over time:

%% RUN THIS SECTION BEFORE PROCEEDING
clearvars('-except',initial_var{:})

% ADD SPEED TO THE DATA STRUCTURE
data.speed = nan(size(data.sleep));
for trial = 1:num.trials
    data.speed(:,F,trial) = fly(trial).f.speed;
    data.speed(:,M,trial) = fly(trial).m.speed;
end

% make an inner food quad region ROI
data.innerFoodQuad = data.foodQuad;
loc = logical(replaceNaN(data.OutterRing,false));
data.innerFoodQuad(loc) = false;

% make a courtship all (un-time-restricted) structure 
data.CI_all = replaceNaN(data.wing_ext_all,false) |...
                        replaceNaN(data.chase_all,false) |...
                        replaceNaN(data.circling_all,false);

% add new variables to the saved variable list
initial_var{end+1} = 'encounters';
initial_var{end+1} = 'foreColor';
foreColor = formattingColors(blkbgd); % get background colors

% build indexes for warm threat, cool threat, warm safe, cool safe
threshTemp = 25; % 'neutral temp'
data.tempbin.WT = data.tempbin.h_idx & data.temp>threshTemp; % warm threat
data.tempbin.WS = data.tempbin.c_idx & data.temp>threshTemp; % warm safe
data.tempbin.CT = data.tempbin.c_idx & data.temp<threshTemp; % cool threat
data.tempbin.CS = data.tempbin.h_idx & data.temp<threshTemp; % cool safe
data.tempbin.SS = ~data.tempbin.h_idx & ~data.tempbin.c_idx; % static safe

% adjust the overshoot temps in F LRR to not count towards a temp region
if contains(groupName, 'F LRR 25-17')
    data.tempbin.WT(data.tempbin.WT) = false;
    data.tempbin.WS(data.tempbin.WS) = false;
end

%% Simple comparison across flies: distance to food over time
clearvars('-except',initial_var{:})

% plot the data:
lw = 1;
plotSkip = 120;
sSpan = plotSkip*fly(1).fps; %  60 second smoothing


% fig = getfig('',1); hold on
% for sex = 1:2
%     % plot all the trials
%     x = repmat(data.time, [1,num.trials]);
%     y = squeeze(data.dist2food(:,sex,:));
%     plot(x,y,'color', data.color(sex,:),'LineWidth',lw)
% end
% xlabel('time (min)')
% ylabel('distance to food (mm)')
% formatFig(fig);
% save_figure(fig, [figDir 'test fly distance to food all trials'],fig_type)


% Plot the average distance to food across the trials: 
r = 5;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2:r;
fig = getfig('',1); 
    subplot(r,c,sb(1).idx);
        plot(data.time,data.temp,'color', foreColor, 'linewidth', 1)
        ylabel('temp (\circC)')
    subplot(r,c,sb(2).idx);
        hold on
        for sex = 1:2
            % plot all the trials
            x = data.time;
            y = squeeze(data.dist2food(:,sex,:));
            y_avg = mean(y,2,'omitnan');
            y_avg = smooth(y_avg, sSpan, 'moving');
            y_err = std(y,0,2,'omitnan')/sqrt(num.trials);
            y_err = smooth(y_err, sSpan, 'moving');
            h = plot_error_fills(true, x(1:plotSkip:end), y_avg(1:plotSkip:end), y_err(1:plotSkip:end), data.color(sex,:),fig_type);
            plot(x(1:plotSkip:end),y_avg(1:plotSkip:end),'color', data.color(sex,:),'LineWidth',lw)
        end
        xlabel('time (min)')
        ylabel('distance to food (mm)')
    formatFig(fig,blkbgd,[r c],sb);
    subplot(r,c,sb(1).idx)
    set(gca, 'xcolor', 'none')
    if strcmp(groupName, 'Berlin LTS caviar')
        ylim([15,35])
        set(gca, 'ytick', 15:10:35)
    end
save_figure(fig, [figDir 'avg distance to food M and F'],fig_type);

% % TODO: 5.13 MAKE A THING TO PLOT TEMP, FOOD METRIC, AND THEN THE
% TEMP-TUNING CURVE

% clearvars('-except',initial_vars{:})
% plot_err = true;
% autoLim = true;
% % Y limit ranges
% dist_lim = [5,35];       %distance
% dt_lim = [10,32];        %distance-temp
% auto_time = true;      % automated time axis limits
% time_lim = [0,400];     %time limit (x-axis)
% nMax =  num.exp;%
% [~,backColor] = formattingColors(blkbgd); %get background colors
% 
% % set up figure aligments
% r = 5; %rows
% c = 3; %columns
% sb(1).idx = [1,2]; %temp timecourse
% sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive?
% sb(3).idx = 3:c:r*c; %binned distance alignment
% 
% LW = 0.75;
% sSpan = 360;
% dataString = cell([1,num.exp]);
% 
% % FIGURE:
% fig = getfig('',true);
% for i = 1:nMax
% %     i = expOrder(ii);
%     x = grouped(i).time;
%     kolor = grouped(i).color;
% 
%     %temp
%     subplot(r,c,sb(1).idx); hold on
%         y = grouped(i).temp;
%         plot(x,y,'LineWidth',2,'Color',kolor)
% 
%     %distance
%     subplot(r,c,sb(2).idx); hold on
%         y = smooth(grouped(i).dist.avg,'moving',sSpan);
% %         y_err = smooth(grouped(i).dist.err,'moving',sSpan);
%         plot(x,y,'LineWidth',LW,'Color',kolor)
%         if ~autoLim
%             ylim(dist_lim)
%         end
% 
%     %temp dependent distance
%     subplot(r,c,sb(3).idx); hold on
%         x = grouped(i).dist.distavgbytemp(:,1);
%         y = grouped(i).dist.distavgbytemp(:,2);
%         y_err = grouped(i).dist.distavgbytemp_err(:,2);
% 
%         plot(x,y,'color',kolor,'linewidth',LW+1)
%         plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
%         dataString{i} = grouped(i).name;
% end
% 
% % FORMATING AND LABELS
% formatFig(fig,blkbgd,[r,c],sb);
% % temp
% subplot(r,c,sb(1).idx)
% ylabel('\circC')
% set(gca,"XColor",backColor)
% if ~auto_time
%     xlim(time_lim)
% end
% % distance
% subplot(r,c,sb(2).idx)
% ylabel('proximity to food (mm)')
% xlabel('time (min)')
% set(gca,'ydir','reverse')
% if ~auto_time
%     xlim(time_lim)
% end
% % temp-distance relationship
% subplot(r,c,sb(3).idx)
% ylabel('proximity to food (mm)')
% xlabel('temp (\circC)')
% if ~autoLim
%     ylim(dt_lim)
% end
% h_line(18.1,'grey',':',1) %36.2
% set(gca,'ydir','reverse')
% %
% % legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)
% 
% % % get the avg distance over the whole experiment
% % roi = [1000 159935];
% % [dist_mean, dist_err] = deal([]);
% % for i = 1:num.exp
% %     dist_all = mean(grouped(i).dist.all(roi(1):roi(2),:),1,'omitnan');
% %     dist_err(i) = std(dist_all);
% %     dist_mean(i) = mean(dist_all);
% % end
% %
% % e = errorbar([17,20,23,25],dist_mean,dist_err);
% % e.Marker  = 'o';
% % e.Color = 'r';
% % e.MarkerSize = 10;
% 
% % save figure
% save_figure(fig,[saveDir expGroup ' timecourse summary no speed food only'],fig_type);

%% FIGURE: generate temperature-tuning curves for metrics of your selection
clearvars('-except',initial_var{:})

% select data metric to plot 
% TODO: make this it's own function that has the appropriate information
% tied to each of the different data types (e.g. sex split, percent or
% distance based etc.)
fields = {'OutterRing', 'sleep',  'foodQuad','speed','innerFoodQuad','dist2food', 'FlyOnFood', 'foodcircle',...
              'eccentricity', 'CI',  'circling_1sec', 'circling_all', 'court_chase', 'chase_all', 'wing_ext', 'wing_ext_all'};
idx = listdlg("PromptString",'Select data type', 'SelectionMode','single', 'ListSize',[200, 300], 'ListString',fields);
sel_field = fields{idx};
double_field = size(data.(sel_field),2)==2; % if there are two fields for both sexes
maxRoi = 730000; % cutoff for food quality loss in the smaller plate

% find averages for each fly for each temperature bin
temps = data.tempbin.temps;
nTemps = length(temps);
[pData.cooling.avg,pData.cooling.sem, pData.warming.avg,pData.warming.sem]  = deal(nan([nTemps,1]));
typeNames = {'cooling', 'warming'};

for type = 1:2 % heating and cooling
    type_name = typeNames{type};
    for t = 1:nTemps
        roi = data.tempbin.(type_name)(:,t);
        roi(maxRoi:end) = false;
        if double_field
            raw = [squeeze(data.(sel_field)(roi,1,:)), squeeze(data.(sel_field)(roi,2,:))];
        else 
            raw = data.(sel_field)(roi,:);
        end
        processed =  mean(raw,1,'omitnan');
        pData.(type_name).avg(t) = mean(processed,'omitnan');
         pData.(type_name).sem(t) = std(processed,'omitnan')/sqrt(length(processed));
    end
end

% plot the data
if strcmp(groupName,'Berlin LTS caviar')
    xlimits = [13, 37];
end

r = 1; 
c = 2; 
LW = 1;
kolor = Color('dodgerBlue'); %foreColor;
plot_err = true;

switch sel_field
    case 'OutterRing'
        scaler = 100;
        ylabel_str = 'edge occupancy (% flies)';
        ylimits = [0, 50];
    case 'speed'
        scaler = 1;
        ylabel_str = 'speed (cm/s)';
        ylimits = [0 18];
    case 'foodQuad'
        scaler = 100;
        ylabel_str = 'food quadrant (% flies)';
        ylimits = [0 90];
    case 'innerFoodQuad'
        scaler = 100;
        ylabel_str = 'inner food quadrant (% flies)';
        ylimits = [0 90];
    case 'foodcircle'
        scaler = 100;
        ylabel_str = 'food circle (% flies)';
        ylimits = [0 65];
    case 'sleep'
        scaler = 100;
        ylabel_str = 'sleeping (% flies)';
        ylimits = [0 60];
    case 'CI'
        scaler = 100;
        ylabel_str = 'sleeping (% flies)';
        ylimits = [0 60];
end

fig = getfig([sel_field ' Tuning Curve'], 1);
    for type = 1:2 % cooling and warming
        subplot(r, c, type)
        hold on
        x = temps;
        y = pData.(typeNames{type}).avg .* scaler;
        y_err = pData.(typeNames{type}).sem .* scaler;
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.4);
        plot(x,y,'color',kolor,'linewidth',LW+1)
        
        % formatting
        if type == 1
            set(gca, 'xdir', 'reverse')
            ylabel(ylabel_str)
        else
            set(gca, 'ycolor', 'none')
        end
        xlabel('temp (\circC)')
        xlim(xlimits)
        ylim(ylimits)
        title(typeNames{type})
    end
formatFig(fig, blkbgd,[r,c]);
            
save_figure(fig, [figDir sel_field ' tuning curve'])



%% plot the female and male positions within the arena ...
% TODO: need to rotate the arena to match the food alignment across trials ...
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); % get background colors

a = inputdlg(['There are ' num2str(num.trials) ' figures, how many columns?']);
c = str2double(a{:});
r = ceil(num.trials/c);

% female position figure: 
fig = getfig('female fly positions in the arena', 1,[1064 774]);
for i = 1:num.trials
    subplot(r,c,i); % male position in arena
    x = fly(i).f.pos(:,body.center,1);
    y = fly(i).f.pos(:,body.center,2);
    scatter(x,y,3, data.color(F,:))
    hold on
    scatter(fly(i).well.food(1),fly(i).well.food(2),10, foreColor)
    viscircles(fly(i).well.center(5,:),30/fly(i).pix2mm,'Color',foreColor)
    axis square equal
    set(gca, 'XColor','none', 'ycolor', 'none')
end
save_figure(fig, [figDir 'female fly positions in arena'],fig_type)

% male position figure: 
fig = getfig('male fly positions in the arena', 1,[1064 774]);
for i = 1:num.trials
    subplot(r,c,i); % male position in arena
    x = fly(i).m.pos(:,body.center,1);
    y = fly(i).m.pos(:,body.center,2);
    scatter(x,y,3, data.color(M,:))
    hold on
    scatter(fly(i).well.food(1),fly(i).well.food(2),10, foreColor)
    viscircles(fly(i).well.center(5,:),30/fly(i).pix2mm,'Color',foreColor);
    axis square equal
    set(gca, 'XColor','none', 'ycolor', 'none')
end
save_figure(fig, [figDir 'male fly positions in arena'],fig_type)


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


%% Courtship frequency figure: TODO-check that this is aligned time! 
clearvars('-except',initial_var{:})
% when and where is the concentration of courtship over experimental time?
foreColor = formattingColors(blkbgd); % get background colors
sex = 1;
% compile the data: 
[pD(1).x, pD(2).x, pD(1).y,pD(2).y] = deal([]);
for i = 1:num.trials

        x = fly(i).time;
        y = fly(i).T.CI;
        pD(sex).x = [pD(sex).x, x];
        pD(sex).y = [pD(sex).y, y];

end

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14];
sb(3).idx = 3:c:r*c; %binned distance alignment
LW = 1;

fig = getfig('',0); 

subplot(r,c,sb(1).idx);
hold on
for i = 1:num.trials
    x = fly(i).time;
    y = fly(i).T.temperature;
    plot(x,y,'color', foreColor,'LineWidth', LW)
end
ylabel('temp (\circC)')

subplot(r,c,sb(2).idx);
hold on
for sex = 1:2
    x = mean(pD(sex).x,2);
    y = sum(pD(sex).y,2);
    plot(x,smooth(y,5*60,'moving'),'color', foreColor,'LineWidth',LW)
end
xlabel('time (min)')
ylabel('courtship index sum')

formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca, 'xcolor', 'none')

% TODO: 
% plot out by temp region (improving vs worsening) the avg 
% plot out the avg over time -- running avg. 

%% Courtship Index: 
% TODO: update this to work with the current data...
clearvars('-except',initial_var{:})
% when and where is the concentration of courtship over experimental time?
[foreColor, ~] = formattingColors(blkbgd); %get background colors

% compile the data: 
[pD(1).x, pD(2).x, pD(1).y,pD(2).y] = deal([]);
for i = 1:num.trials
    sex = 1;
    % for sex = 1:2
        x = fly(i).time;
        y = fly(i).T.CI;
        pD(sex).x = [pD(sex).x, x];
        pD(sex).y = [pD(sex).y, y];
    % end
end

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14];
sb(3).idx = 3:c:r*c; %binned distance alignment
LW = 1;

fig = getfig('',1); 

subplot(r,c,sb(1).idx);
hold on
for i = 1:num.trials
    x = fly(i).time;
    y = fly(i).T.temperature;
    plot(x,y,'color', foreColor,'LineWidth', LW)
end
ylabel('temp (\circC)')

subplot(r,c,sb(2).idx);
hold on
for sex = 1:2
    x = mean(pD(sex).x,2);
    y = sum(pD(sex).y,2);
    plot(x,smooth(y,5*60,'moving'),'color', foreColor,'LineWidth',lw)
end
xlabel('time (min)')
ylabel('courtship index sum')

formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca, 'xcolor', 'none')

%% TODO: group behavior state transition map

%% TODO:  grouped Courtship behavior frequency time course
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); %get background colors
kolor = Color('gold');
% simple plot of when the different courtship behaviors are happening
r = 5;
c = 1;
lw = 2;

time = mean(data.time,2,'omitnan');
temp = mean(data.temperature,2,'omitnan');

fig = getfig('',0);
% time
subplot(r,c,1); hold on 
% x = data.time;
% y = data.temp;
% plot(x,y,'color', foreColor, 'linewidth', lw)
% ylabel('\circC')
plot(time, temp,'color', foreColor, 'linewidth', lw)
ylabel('\circC')

% full courtship index
subplot(r,c,2); hold on 
    x = time;
    y = sum(data.CI,2);
    plot(x,y,'color', foreColor, 'linewidth', lw)
    ylabel('CI')

% wing extension
subplot(r,c,3); hold on 
    x = time;
    y = sum(data.wing_ext_all,2);
    plot(x,y,'color', foreColor, 'linewidth', lw)
    y = sum(data.wing_ext,2);
    plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('wing ext')

% chasing
subplot(r,c,4); hold on 
    x = time;
    y = sum(data.chase_all,2);
    plot(x,y,'color', foreColor, 'linewidth', lw)
    y = sum(data.court_chase,2);
    plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('chase')

% circling
subplot(r,c,5); hold on 
    x = time;
    y = sum(data.circling_all,2);
    plot(x,y,'color', foreColor, 'linewidth', lw)
    y = sum(data.circling_1sec,2);
    plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('circling')
    xlabel('time (min)')

% formating: 
formatFig(fig, blkbgd,[r,c]);
trans = [data.cooling_idx, data.warming_idx(2)];
for i = 1:r
    subplot(r,c,i)
    v_line(data.time(trans),'r', '-',1)
    if i<r
        set(gca, 'xcolor', 'none')
    end
end

save_figure(fig, [figDir 'courtship behaviors over time'],fig_type);


%% FIGURE: courtship temp tuning curves
% TODO: turn this into a temperature tuning curve for each of these metrics
ntemps = length(data.tempbin.temps);
[TC.CI.avg, TC.CI.std] = deal(nan(ntemps,2)); % warming then cooling for columns
for t = 1:ntemps
    for type = 1:2 %warming then cooling
        switch type 
            case 1 
                ROI = data.tempbin.warming(:,t);
            case 2
                ROI = data.tempbin.cooling(:,t);
        end
        y = (sum(data.CI(ROI,:),2)./num.trials)*100; % percent of flies doing 'official' courtship
        y_avg = mean(y,'omitnan');
        y_err = std(y, 0,1,'omitnan');
        TC.CI.avg(t,type) = y_avg;
        TC.CI.std(t,type) = y_err;
    end
end

% Plot
sSpan = 8;

fig = getfig('', 1,[486 680]);
hold on
    x = data.tempbin.temps';
    for type = 1:2
        switch type
            case 1
                kolor = Color('red');
            case 2
                kolor = Color('dodgerblue');
        end
        % smooth / format
        y = smooth(TC.CI.avg(:,type),sSpan,'moving');
        y_err = smooth(TC.CI.std(:,type),sSpan,'moving')./sqrt(num.trials); % SEM     
        y_low = y-y_err;
        y_high = y+y_err;
        y_low(y_low<0) = 0; % threshold error to zero
        % plot 
        plot_error_fills(blkbgd, x, y, y_err, kolor,fig_type, 0.35);
        if ~blkbgd
            plot(x, y_low, 'color', kolor, 'linewidth',0.25)
            plot(x, y_high, 'color', kolor, 'linewidth',0.25)
        end
        plot(x, y, 'color', kolor, 'linewidth',1)
    end

 %   formatting
 formatFig(fig, blkbgd);
 xlabel('temperature (\circC)')
 ylabel('flies courting (%)')

 save_figure(fig, [figDir 'courtship CI index temp tuning curve'],fig_type);

 


%% TODO: behavior probabilities for states that happen before sleep & how long it took between them...
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); %get background colors

% what is the behavior that happened just prior to sleep and how long ago
% did it happen?
[behavior_next, behavior_last] = deal([]);
index = 1;
for i = 1:num.trials
    for sex = 1:2
        loc = find(diff(data.sleep(:,sex,i))==1)+1;  % frame number of when sleep starts
        sleep_off = find(diff(data.sleep(:,sex,i))==-1);
        % behavior comparisions
        if ~isempty(loc)
            for sleep = 1:length(loc)
                % last behavior that happened before sleep
                outRing_last = find(data.OutterRing(1:loc(sleep),sex,i), 1, 'last' ); % outer ring 
                outRing_last = empty2nan(outRing_last);
                onFood_last = find(data.FlyOnFood(1:loc(sleep),sex,i), 1, 'last' ); % on food
                onFood_last = empty2nan(onFood_last);
                if sex==1
                    courtship_last = find(data.FlyOnFood(1:loc(sleep),sex,i), 1, 'last' ); % on food
                    courtship_last = empty2nan(courtship_last);
                else
                    courtship_last = nan;
                end
                [frame_last, idx_last] = max([outRing_last, onFood_last, courtship_last]);
                %save into large matrix
                behavior_last(index,:) = [i, loc(sleep), idx_last, frame_last, outRing_last, onFood_last, courtship_last];

                % next behavior that happened after sleep
                outRing_next = find(data.OutterRing(sleep_off(sleep)+1:end,sex,i), 1, 'first' )+sleep_off(sleep)+1; % outer ring 
                outRing_next = empty2nan(outRing_next);
                onFood_next = find(data.FlyOnFood(sleep_off(sleep)+1:end,sex,i), 1, 'first' )+sleep_off(sleep)+1; % on food
                onFood_next = empty2nan(onFood_next);
                if sex==1
                    courtship_next = find(data.FlyOnFood(sleep_off(sleep)+1:end,sex,i), 1, 'first' )+sleep_off(sleep)+1; % on food
                    courtship_next = empty2nan(courtship_next);
                else
                    courtship_next = nan;
                end
                [frame_next, idx_next] = min([outRing_next, onFood_next, courtship_next]);
                %save into large matrix
                behavior_next(index,:) = [i, sleep_off(sleep), idx_next, frame_next, outRing_next, onFood_next, courtship_next];
                
                % increase the matrix count: 
                index = index+1;
            end
        end
    end
    disp(i)
end

% figure -- what were the last and next behaviors?
time_to_sleep = ((behavior_last(:,2)-behavior_last(:,5:7))./fly(1).fps)./60; % in minutes
time_since_sleep = ((behavior_next(:,5:7) - behavior_next(:,2))./fly(1).fps)./60; % in minutes

beh_list = {'outer ring', 'food', 'courtship'};

fig = getfig('',0);
subplot(1,2,1); hold on
    x = repmat([1,2,3],[size(time_to_sleep,1),1]);
    plot(x', time_to_sleep')
    scatter(x, time_to_sleep,35,foreColor,'filled')
    xlim([0,4])
    set(gca, 'xtick', 1:3,'xticklabel', beh_list)
    ylabel('time from last behavior to sleep (min)')
    title('Previous Behavior')
subplot(1,2,2); hold on
    x = repmat([1,2,3],[size(time_since_sleep,1),1]);
    plot(x', time_since_sleep')
    scatter(x, time_since_sleep,35,foreColor,'filled')
    xlim([0,4])
    set(gca, 'xtick', 1:3,'xticklabel', beh_list)
    ylabel('time to next behavior after sleep (min)')
    title('Subsequent Behavior')
formatFig(fig, blkbgd,[1,2]);

save_figure(fig, [figDir 'behaviors before and after sleep'],fig_type);

%% TODO: add a frequency scatter plot that shows the number of instances per temperature region
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); %get background colors

temp_regimes = {'hold', 'cooling', 'warming', 'hold'};
ntypes = length(temp_regimes);
data_type = 'CI';

% data_type = 'wing_ext_all';
% data_type = 'chase_all';

% Extract data for plotting
plotData = [];
for t = 1:ntypes
    switch t
        case 1 % hold
            t_name = 'pre hold';
            idx = 1:data.cooling_idx(1)-1;
        case 2 % cooling
            t_name = 'cooling';
            idx = data.cooling_idx(1):data.cooling_idx(2);
        case 3 % warming
            t_name = 'warming';
            idx = data.warming_idx(1):data.warming_idx(2);
        case 4 % post hold
            t_name = 'post hold';
            idx = data.warming_idx(2)+1:length(data.temp);
    end
    
    % have some switch mechanism here to look at different types of parameters
    raw_data = data.(data_type);
    
    y_raw = sum(raw_data(idx,:),'omitnan');
    y = y_raw./(length(idx)/(fly(1).fps*60)); % instances per minute
    
    plotData(t,:) = y;
end

% Visualize Data:
buff = 0.2;
avg_buff = 0.3;
sz = 60;
lw = 2;
cList = {'grey', 'dodgerblue', 'red', 'grey'};

fig = getfig('',1,[547 526]);
hold on
for t = 1:ntypes
    x = shuffle_data(linspace(t-buff,t+buff,num.trials));
    y = plotData(t,:);
    scatter(x,y,sz, Color(cList{t}), "filled")
    x = [t-avg_buff, t+avg_buff];
    y_mean = mean(y, 'omitnan');
    plot(x,[y_mean, y_mean],"Color",foreColor, 'linewidth', lw)
    y_err = std(y, 'omitnan')/num.trials;
    y_sem = [y_mean-y_err, y_mean+y_err];
    plot(x,[y_sem(1), y_sem(1)], 'Color',foreColor, 'linewidth', 0.5, 'linestyle', '--')
    plot(x,[y_sem(2), y_sem(2)], 'Color',foreColor, 'linewidth', 0.5, 'linestyle', '--')
end
set(gca, 'xtick', 1:ntypes,'xticklabel', temp_regimes)
ylabel([strrep(data_type,'_', '-') ' frequency (#/min)'])
formatFig(fig, blkbgd);
save_figure(fig, [figDir 'temp regime binned frequency of ' data_type],fig_type);

% run statistical tests: are they different? 
statData = plotData';
[p,tbl,stats] = kruskalwallis(statData,[],'off');
fig2 = getfig('', 1, [633 580]);
c = multcompare(stats);
tble = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
formatFig(fig2, blkbgd);

%% TODO: duration on food across different temperatures
% Question: when are the flies going onto the food and how long are they staying
% there? How does this change for different temperature regimes?
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); %get background colors

temp_regimes = {'hold', 'cooling', 'warming', 'hold'};
ntypes = length(temp_regimes);
% data_type = 'CI';
% data_type = 'wing_ext';
% data_type = 'chase_all';
data_type = 'FlyOnFood';

% Extract data for plotting
plotData = [];
for t = 1:ntypes
    switch t
        case 1 % hold
            t_name = 'pre hold';
            idx = 1:data.cooling_idx(1)-1;
        case 2 % cooling
            t_name = 'cooling';
            idx = data.cooling_idx(1):data.cooling_idx(2);
        case 3 % warming
            t_name = 'warming';
            idx = data.warming_idx(1):data.warming_idx(2);
        case 4 % post hold
            t_name = 'post hold';
            idx = data.warming_idx(2)+1:length(data.temp);
    end

    % have some switch mechanism here to look at different types of parameters
    raw_data = data.(data_type);
    
    
    y_raw = sum(raw_data(idx,:),'omitnan');
    y = y_raw./(length(idx)/(fly(1).fps*60)); % instances per minute
    
    plotData(t,:) = y;
end

% Visualize Data:
buff = 0.2;
avg_buff = 0.3;
sz = 60;
lw = 2;
cList = {'grey', 'dodgerblue', 'red', 'grey'};

fig = getfig('',1,[547 526]);
hold on
for t = 1:ntypes
    x = shuffle_data(linspace(t-buff,t+buff,num.trials));
    y = plotData(t,:);
    scatter(x,y,sz, Color(cList{t}), "filled")
    x = [t-avg_buff, t+avg_buff];
    y_mean = mean(y, 'omitnan');
    plot(x,[y_mean, y_mean],"Color",foreColor, 'linewidth', lw)
    y_err = std(y, 'omitnan')/num.trials;
    y_sem = [y_mean-y_err, y_mean+y_err];
    plot(x,[y_sem(1), y_sem(1)], 'Color',foreColor, 'linewidth', 0.5, 'linestyle', '--')
    plot(x,[y_sem(2), y_sem(2)], 'Color',foreColor, 'linewidth', 0.5, 'linestyle', '--')
end
set(gca, 'xtick', 1:ntypes,'xticklabel', temp_regimes)
ylabel([strrep(data_type,'_', '-') ' frequency (#/min)'])
formatFig(fig, blkbgd);
save_figure(fig, [figDir 'temp regime binned frequency of ' data_type],fig_type);

%% FIGURE: rasterplot of frequency of courtship-like events over the timecourse

% bin by time?

% look at the total number for each group...

% bin by 1 minute increments

fps = fly(1).fps;
binWidth = 1; % seconds

binedges = 1:binWidth*fps:size(data.time,1);
pD = [];
[pD.wingext, pD.chase, pD.circling, pD.CI, pD.temp,pD.time] = deal(nan(length(binedges)-1,1));
for i = 1:length(binedges)-1
    roi = binedges(i):binedges(i+1);
    pD.wingext(i) = sum(sum(data.wing_ext_all(roi,:),'omitnan'),'omitnan');
    pD.chase(i) = sum(sum(data.chase_all(roi,:),'omitnan'),'omitnan');
    pD.circling(i) = sum(sum(data.circling_all(roi,:),'omitnan'),'omitnan');
    pD.CI(i) = sum(sum(data.CI(roi,:),'omitnan'),'omitnan');
    pD.temp(i) = mean(data.temp(roi),'omitnan');
    pD.time(i) = mean(data.time(roi),'omitnan');
end
% turn these into rates
params = {'chase', 'circling', 'wingext'};
sSpan = 10;
for i = 1:length(params)
    pD.(params{i}) = (pD.(params{i})./binWidth)./num.trials;
    pD.(params{i}) = smooth(pD.(params{i}),sSpan,'moving');
end


r = 7;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;
sb(3).idx = 4:5;
sb(4).idx = 6:7;

LW = 1;
params = {'temp', 'chase', 'circling', 'wingext'};

fig = getfig('',1,[665 680]);
for i = 1:length(params)
    subplot(r,c,sb(i).idx)
    plot(pD.time, pD.(params{i}),'color', Color('black'),'linewidth', LW)
    ylabel([params{i} ' (hz)'])
end
formatFig(fig, blkbgd,[r c],sb);
for i = 1:length(params)-1
    subplot(r,c,sb(i).idx)
    set(gca,'xcolor','none')
end
save_figure(fig, [figDir 'courstship measures over time'],'-pdf');


% Temp region summaries:
tROI = [1,989; 990 1959; 1960 2934; 2935 3839]; % hold, cooling, warming, hold
figure; 
plot(pD.temp)
v_line(tROI)

PD = [];
for t = 1:4 % hold, cooling, warming, hold
    for i = 2:4 % params other than temp
        PD(t,i-1) = sum(pD.(params{i})(tROI(t,1):tROI(t,2)),'omitnan')./(tROI(t,2)-tROI(t,1));
    end
end


cData = [Color('grey'); Color('dodgerblue'); Color('red'); Color('grey')];
n = 3;
fig = getfig('',0); 
for i = 1:n
    subplot(1,n,i)
    h = bar(PD(:,i));
    h.FaceColor = 'flat';
    h.CData = cData;

    ylabel([params{i+1} ' (hz)'])
end
formatFig(fig,blkbgd,[1,n])
for i = 1:n
    subplot(1,n,i)
    set(gca, 'xcolor', 'none')
end
save_figure(fig, [figDir 'courstship measures avg over temp region'],'-pdf');


%% Trial by trial comparision of the different courtship measures 
clearvars('-except',initial_var{:})

% UPDATE THIS FOR NOT LTS: 
data.temp = mean(data.temperature,2,'omitnan');

% bin by 1 minute increments

fps = fly(1).fps;
binWidth = 1; % seconds

binedges = 1:binWidth*fps:size(data.time,1);
pD = [];
[pD.wingext, pD.chase, pD.circling, pD.CI] = deal(nan(length(binedges)-1,num.trials));
[pD.temp,pD.time] = deal(nan(length(binedges)-1,1));

for i = 1:length(binedges)-1
    roi = binedges(i):binedges(i+1);
    pD.wingext(i,:) = (sum(data.wing_ext_all(roi,:),1,'omitnan'));
    pD.chase(i,:) = (sum(data.chase_all(roi,:),1,'omitnan'));
    pD.circling(i,:) = (sum(data.circling_all(roi,:),1,'omitnan'));
    pD.CI(i,:) = (sum(data.CI(roi,:),1,'omitnan'));
    pD.temp(i) = mean(mean(data.temp(roi),'omitnan'));
    pD.time(i) = mean(mean(data.time(roi),'omitnan'));
end
% turn these into rates
params = {'chase', 'circling', 'wingext'};
sSpan = 10;
for i = 1:length(params)
    pD.(params{i}) = (pD.(params{i})./binWidth)./num.trials;
    % pD.(params{i}) = smooth(pD.(params{i}),sSpan,'moving');
end

% Temp region summaries:
tROI = [1,989; 990 1959; 1960 2934; 2935 3839]; % hold, cooling, warming, hold
figure; 
plot(pD.temp)
v_line(tROI)

% Bin the data for each of the temp epochs
PD = [];
for t = 1:4 % hold, cooling, warming, hold
    for i = 1:3 % params other than temp
        PD(t,1:num.trials,i) = sum(pD.(params{i})(tROI(t,1):tROI(t,2),:),'omitnan')./(tROI(t,2)-tROI(t,1));
    end
end

% FIGURE: visualize courtship measures for each of the trials
r = 1;
c = 3;
sz = 30;
kolor = {'grey', 'dodgerblue', 'red', 'grey'};

fig = getfig('', 0);
    for i = 1:c
        subplot(r,c,i); hold on
        for t = 1:4
            x = t*ones(1,num.trials);
            y = PD(t,:,i);
            scatter(x,y,sz, Color(kolor{t}),'filled','XJitter','density')
            y_avg = mean(y,'omitnan');
            plot([t-0.3,t+0.3],[y_avg, y_avg],'color', Color(kolor{t}),'linewidth', 2)
        end
        title(params{i})
    end
    formatFig(fig, blkbgd,[r,c]);
    for i = 1:c
        subplot(r,c,i); hold on
        set(gca, 'xcolor', 'none')
    end

gNames = [ones(num.trials,1),2*ones(num.trials,1),3*ones(num.trials,1),4*ones(num.trials,1)];
% STATS:
stats = [];
for i = 1:c
    inputData = squeeze(PD(:,:,i)');
    [~,~,stats(i).s] = anova1(inputData);
    alpha = 0.05; %significance level
    [stats(i).c,~,~,~] = multcompare(stats(i).s,alpha,'off');
end

% so far not 

%% FIGURE:  plot major components of courtship behavior over time just raw numbers
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); %get background colors
kolor = Color('gold');

% simple plot of when the different courtship behaviors are happening
r = 5;
c = 1;
lw = 2;
sSpan = 60*30; % 1 min smoothing bin for the Courtship Index in time

time = mean(data.time,2,'omitnan');
temp = mean(data.temperature,2,'omitnan');

fig = getfig('',0);
% time
subplot(r,c,1); hold on 
% x = data.time;
% y = data.temp;
% plot(x,y,'color', foreColor, 'linewidth', lw)
% ylabel('\circC')
plot(time, temp,'color', foreColor, 'linewidth', lw)
ylabel('\circC')

% full courtship index
subplot(r,c,2); hold on 
    x = time;
    y = sum(data.CI,2);
    y1 = smooth(y,sSpan, 'moving');
    plot(x,y,'color', foreColor, 'linewidth', lw)
    ylabel('CI')
    yyaxis right
    plot(x,y1,'color', Color('dodgerblue'), 'linewidth', lw)
    ylabel('CI minute avg')

% wing extension
subplot(r,c,3); hold on 
    x = time;
    y = sum(data.wing_ext_all,2);
    y1 = smooth(y,sSpan, 'moving');
    plot(x,y,'color', foreColor, 'linewidth', lw)
    y = sum(data.wing_ext,2);
    plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('wing ext')
    yyaxis right
    plot(x,y1,'color', Color('dodgerblue'), 'linewidth', lw)

% chasing
subplot(r,c,4); hold on 
    x = time;
    y = sum(data.chase_all,2);
    y1 = smooth(y,sSpan, 'moving');
    plot(x,y,'color', foreColor, 'linewidth', lw)
    y = sum(data.court_chase,2);
    plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('chase')
    yyaxis right
    plot(x,y1,'color', Color('dodgerblue'), 'linewidth', lw)

% circling
subplot(r,c,5); hold on 
    x = time;
    y = sum(data.circling_all,2);
    y1 = smooth(y,sSpan, 'moving');
    plot(x,y,'color', foreColor, 'linewidth', lw)
    y = sum(data.circling_1sec,2);
    plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('circling')
    xlabel('time (min)')
    yyaxis right
    plot(x,y1,'color', Color('dodgerblue'), 'linewidth', lw)

% formating: 
formatFig(fig, blkbgd,[r,c]);
% plot the transition points between heating/cooling regimes...
% trans = [data.cooling_idx, data.warming_idx(2)];
for i = 1:r
    subplot(r,c,i)
    % v_line(data.time(trans),'r', '-',1)
    if i<r
        set(gca, 'xcolor', 'none')
    end
end

save_figure(fig, [figDir 'courtship behaviors over time'],fig_type);

%% FIGURE:  plot SMOOTHED avg major components of courtship behavior over time
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); %get background colors
kolor = Color('dodgerblue');

% simple plot of when the different courtship behaviors are happening
r = 5;
c = 1;
lw = 2;
sSpan = 60*30; % 1 min smoothing bin for the Courtship Index in time

time = mean(data.time,2,'omitnan');
temp = mean(data.temperature,2,'omitnan');

fig = getfig('',0);
% time
subplot(r,c,1); hold on 
% x = data.time;
% y = data.temp;
% plot(x,y,'color', foreColor, 'linewidth', lw)
% ylabel('\circC')
plot(time, temp,'color', foreColor, 'linewidth', lw)
ylabel('\circC')

% full courtship index
subplot(r,c,2); hold on 
    x = time;
    y = mean(data.CI,2,'omitnan');
    y1 = smooth(y,sSpan, 'moving');
    % plot(x,y,'color', foreColor, 'linewidth', lw)
    % ylabel('CI')
    plot(x,y1,'color', kolor, 'linewidth', lw)
    ylabel('CI min avg')

% wing extension
subplot(r,c,3); hold on 
    x = time;
    y = mean(data.wing_ext_all,2,'omitnan');
    y1 = smooth(y,sSpan, 'moving');
    % plot(x,y,'color', foreColor, 'linewidth', lw)
    % y = sum(data.wing_ext,2);
    % plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('wing ext')
    % yyaxis right
    plot(x,y1,'color', kolor, 'linewidth', lw)

% chasing
subplot(r,c,4); hold on 
    x = time;
    y = mean(data.chase_all,2,'omitnan');
    y1 = smooth(y,sSpan, 'moving');
    % plot(x,y,'color', foreColor, 'linewidth', lw)
    % y = sum(data.court_chase,2);
    % plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('chase')
    % yyaxis right
    plot(x,y1,'color', kolor, 'linewidth', lw)

% circling
subplot(r,c,5); hold on 
    x = time;
    y = mean(data.circling_all,2,'omitnan');
    y1 = smooth(y,sSpan, 'moving');
    % plot(x,y,'color', foreColor, 'linewidth', lw)
    % y = sum(data.circling_1sec,2);
    % plot(x,y,'color', kolor, 'linewidth', lw)
    ylabel('circling')
    xlabel('time (min)')
    % yyaxis right
    plot(x,y1,'color', kolor, 'linewidth', lw)

% formating: 
formatFig(fig, blkbgd,[r,c]);
% plot the transition points between heating/cooling regimes...
% trans = [data.cooling_idx, data.warming_idx(2)];
for i = 1:r
    subplot(r,c,i)
    % v_line(data.time(trans),'r', '-',1)
    if i<r
        set(gca, 'xcolor', 'none')
    end
end

save_figure(fig, [figDir 'courtship behaviors smoothed over time'],fig_type);

%% Figure: scatter plot the total quantity of courtship for each pair of flies
clearvars('-except',initial_var{:})

SZ = 40;
kolor = Color('dodgerblue');
startX = 1;
stopX = 2;
buff = 0.5;
xlimits = [startX-1, stopX+1];
LW = 2;

x = shuffle_data(linspace(startX,stopX,num.trials));

fieldList = {'CI', 'wing_ext_all', 'chase_all', 'circling_all'};
ylabList = {'CI', 'wing extension', 'chase', 'circling'};
ylimits = [7 4 25 3.5];
r = 1;
c = length(fieldList);

fig = getfig('',0);

for i = 1:c
subplot(r,c,i); hold on 
    y = sum(data.(fieldList{i}),1,'omitnan');
    y = (y./fps)./60; % how many minutes of behavior was there
    scatter(x, y, SZ, kolor,'filled')
    plot([startX-buff, stopX+buff], [mean(y), mean(y)],'color', kolor,'linewidth', LW)
    xlim(xlimits)
    ylabel(['total ' ylabList{i} ' (min)'])
end
formatFig(fig, blkbgd,[r c]);
for i = 1:c
    subplot(r,c,i)
    set(gca, 'xcolor', 'none')
    ylim([0, ylimits(i)])
end

save_figure(fig, [figDir 'total time for each behavior'],fig_type);


%% TODO: What are the speeds within each region of the arena? on avg and by temperature

%% FIGURE: overlay of courtship, sleep, escape, and food attraction
clearvars('-except',initial_var{:})
[foreColor, ~] = formattingColors(blkbgd); %get background colors

ntemps = length(data.tempbin.temps);
[TC.CI.avg, TC.CI.std] = deal(nan(ntemps,2)); % cooling then warming for columns
[TC.sleep.avg, TC.sleep.std] = deal(nan(ntemps,2)); 
[TC.food.avg, TC.food.std] = deal(nan(ntemps,2)); 
[TC.escape.avg, TC.escape.std] = deal(nan(ntemps,2)); 

paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing'};
colors = {'purple', 'dodgerblue', 'gold', 'red'}; % CI, sleep, food, escape
roi = 1:640000;

for t = 1:ntemps
    for type = 1:2 % cooling then warming
        switch type 
            case 1 
                ROI = data.tempbin.cooling(roi,t);
            case 2
                ROI = data.tempbin.warming(roi,t);
        end
        % Courtship Index
        y = (sum(data.CI(ROI,:),2)./num.trials)*100; % percent of flies doing 'official' courtship
        y_avg = mean(y,'omitnan');
        y_err = std(y, 0,1,'omitnan');
        TC.CI.avg(t,type) = y_avg;
        TC.CI.std(t,type) = y_err;
        % sleep, food, escape
        for i = 2:4 % for the other three parameter metrics
            y = squeeze(mean(data.(paramList{i})(ROI,:,:),2,'omitnan'));
            y = mean(y*100,1); % convert to percent of flies per trial that are sleeping
            y_avg = mean(y,'omitnan');
            y_err = std(y, 0,2,'omitnan');
            TC.(paramList{i}).avg(t,type) = y_avg;
            TC.(paramList{i}).std(t,type) = y_err;
        end
    end
end

% Plot
plotErr = false;
sSpan = 8;
r = 1;
c = 2; 

x = data.tempbin.temps';

fig = getfig('', 1,[774 680]);
% cooling
for type = 1:2
    subplot(r,c,type); hold on
    for param = length(paramList):-1:1
        % if param == 1
        %     yyaxis right
        % else 
        %     yyaxis left
        % end
        kolor = Color(colors{param});
        % smooth / format
        y = smooth(TC.(paramList{param}).avg(:,type),sSpan,'moving');
        y_err = smooth(TC.(paramList{param}).std(:,type),sSpan,'moving')./sqrt(num.trials); % SEM     
        y_low = y-y_err;
        y_high = y+y_err;
        y_low(y_low<0) = 0; % threshold error to zero
        % plot 
        plot_error_fills(blkbgd, x, y, y_err, kolor,fig_type, 0.35);
        if ~blkbgd
            plot(x, y_low, 'color', kolor, 'linewidth',0.25)
            plot(x, y_high, 'color', kolor, 'linewidth',0.25)
        end
        plot(x, y, 'color', kolor, 'linewidth',1)
    end
end
% formatting
formatFig(fig, blkbgd, [r,c]);
matchAxis(fig, true);
subplot(r,c,2)
set(gca, 'ycolor', 'none')
xlabel('temperature (\circC)')
title('warming','color', foreColor,'fontname', 'Arial','FontAngle','italic')
subplot(r,c,1)
ylabel('flies (%)')
set(gca, 'xdir', 'reverse')
xlabel('temperature (\circC)')
title('cooling','color', foreColor,'fontname','Arial','FontAngle','italic')
if strcmp(groupName, 'Berlin LTS caviar')
    for i = 1:2
        subplot(r,c,i)
        set(gca, 'xtick', 15:5:35)
        xlim([13, 37])
    end
end


 % save_figure(fig, [figDir 'courtship CI index temp tuning curve'],fig_type);


%% FIGURE: z-score overlay of courtship, sleep, escape, and food attraction
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors
paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing'};
colors = {'purple', 'dodgerblue', 'gold', 'red'}; % CI, sleep, food, escape
ntemps = length(data.tempbin.temps);
nParams = length(paramList);
TC = struct;
for i = 1:nParams
    [TC.(paramList{i}).avg, TC.(paramList{i}).std] = deal(nan(ntemps,2)); % cooling then warming for columns
end

roi = 1:640000; % eliminate the last bit of the experiment since they fall off the food...for the LTS trials

% preformat the z-score for the data that will be extracted and used: 
Zdata = [];
Zdata(:,1) = (sum(data.CI,2)./num.trials)*100; % courtship index
for i = 2:4
    Zdata(:,i) = sum(sum(data.(paramList{i}),2),3)./(num.trials*2)*100;
end
Zdata(isnan(Zdata)) = 0;
Zdata = zscore(Zdata);

% % dirty time course plot of the zscores
% figure;
% hold on
% for i = 1:4
%     y = smooth(Zdata(:,i),5000,'moving');
%     plot(y(1:100:end))
%     disp(i)
% end

for t = 1:ntemps
    for type = 1:2 % cooling then warming
        switch type 
            case 1 
                ROI = data.tempbin.cooling(roi,t);
            case 2
                ROI = data.tempbin.warming(roi,t);
        end
        
        for i = 1:nParams % courtship, sleep, food, escape
            y = Zdata(ROI,i);
            y_avg = mean(y,'omitnan');
            TC.(paramList{i}).avg(t,type) = y_avg;
            if ~isempty(y)
                y_err = std(y, 0,1,'omitnan');
                TC.(paramList{i}).std(t,type) = y_err;
            end
        end
    end
end

% Plot
plotErr = false;
sSpan = 8;
r = 1;
c = 2; 

x = data.tempbin.temps';

fig = getfig('', 1,[774 680]);
% cooling
for type = 1:2
    subplot(r,c,type); hold on
    for param = 1: nParams
        % if param == 1
        %     yyaxis right
        % else 
        %     yyaxis left
        % end
        kolor = Color(colors{param});
        % smooth / format
        y = smooth(TC.(paramList{param}).avg(:,type),sSpan,'moving');
        y_err = smooth(TC.(paramList{param}).std(:,type),sSpan,'moving')./sqrt(num.trials); % SEM     
        y_low = y-y_err;
        y_high = y+y_err;
        % plot 
        plot_error_fills(blkbgd, x, y, y_err, kolor,fig_type, 0.35);
        if ~blkbgd
            plot(x, y_low, 'color', kolor, 'linewidth',0.25)
            plot(x, y_high, 'color', kolor, 'linewidth',0.25)
        end
        plot(x, y, 'color', kolor, 'linewidth',1)
    end
end
% formatting
formatFig(fig, blkbgd, [r,c]);
matchAxis(fig, true);
subplot(r,c,2)
set(gca, 'ycolor', 'none')
xlabel('temperature (\circC)')
title('warming','color', foreColor,'fontname', 'Arial','FontAngle','italic')
subplot(r,c,1)
ylabel('zscore')
set(gca, 'xdir', 'reverse')
xlabel('temperature (\circC)')
title('cooling','color', foreColor,'fontname','Arial','FontAngle','italic')
if strcmp(groupName, 'Berlin LTS caviar')
    for i = 1:2
        subplot(r,c,i)
        set(gca, 'xtick', 15:5:35)
        xlim([13, 37])
    end
end

%% FIGURE: normalized 0 to max temp tuning curves for the four essential behaviors
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors
paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing'};
colors = {'purple', 'dodgerblue', 'gold', 'red'}; % CI, sleep, food, escape
ntemps = length(data.tempbin.temps);
nParams = length(paramList);
TC = struct;
for i = 1:nParams
    [TC.(paramList{i}).avg, TC.(paramList{i}).std] = deal(nan(ntemps,2)); % cooling then warming for columns
end

% roi = 1:640000; % eliminate the last bit of the experiment since they fall off the food...for the LTS trials

% preformat the z-score for the data that will be extracted and used: 
Zdata = [];
Zdata(:,1) = (sum(data.CI,2)./num.trials); % courtship index
for i = 2:nParams
    Zdata(:,i) = sum(sum(data.(paramList{i}),2),3)./(num.trials*2);
end
Zdata(isnan(Zdata)) = 0;
% for i = 1:nParams
%     Zdata(:,i) = smooth(Zdata(:,i),200,'moving');
%     Zdata(:,i) = rescale(Zdata(:,i));
% end
% mVal = max(Zdata);
% Zdata = Zdata./mVal; % normalize on a [0-1] scale

% % dirty time course plot of the zscores
% figure;
% hold on
% for i = 1:4
%     % y = smooth(Zdata(:,i),5000,'moving');
%     y = Zdata(:,i);
%     plot(y(1:100:end))
% end

for t = 1:ntemps
    for type = 1:2 % cooling then warming
        switch type 
            case 1 
                ROI = data.tempbin.cooling(:,t);
            case 2
                ROI = data.tempbin.warming(:,t);
        end
        
        for i = 1:nParams % courtship, sleep, food, escape
            y = Zdata(ROI,i);
            y_avg = mean(y,'omitnan');
            TC.(paramList{i}).avg(t,type) = y_avg;
            if ~isempty(y)
                y_err = std(y, 0,1,'omitnan');
                TC.(paramList{i}).std(t,type) = y_err;
            end
        end
    end
end

% Plot
plotErr = false;
sSpan = 8;
r = 1;
c = 2; 

x = data.tempbin.temps';

% fig = getfig('', 1,[774 680]);
% for param = 1: nParams
%     raw_y = TC.(paramList{param}).avg;
%     for i = 1:size(raw_y,2)
%         raw_y(:,i) = smooth(raw_y(:,i), sSpan, 'moving');
%     end
%     scaleY = rescale(raw_y);
%     for type = 1:2
%         subplot(r,c,type); hold on
%         kolor = Color(colors{param});
%         y = scaleY(:,type);
%         plot(x, y, 'color', kolor, 'linewidth',2)
%     end
% end

fig = getfig('', 0,[1344 680]);
for type = 1:2
    subplot(r,c,type); hold on
    for param = 1: nParams
        raw_y = smooth(TC.(paramList{param}).avg(:,type), sSpan, 'moving');
        scaleY = rescale(raw_y);
        kolor = Color(colors{param});
        plot(x, scaleY, 'color', kolor, 'linewidth',2)
    end
end

% formatting
formatFig(fig, blkbgd, [r,c]);
matchAxis(fig, true);
subplot(r,c,2)
set(gca, 'ycolor', 'none')
xlabel('temperature (\circC)')
title('warming','color', foreColor,'fontname', 'Arial','FontAngle','italic')
subplot(r,c,1)
ylabel('flies (norm %)')
set(gca, 'xdir', 'reverse')
xlabel('temperature (\circC)')
title('cooling','color', foreColor,'fontname','Arial','FontAngle','italic')
if strcmp(groupName, 'Berlin LTS caviar')
    for i = 1:2
        subplot(r,c,i)
        set(gca, 'xtick', 15:5:35)
        xlim([13, 37])
    end
end

save_figure(fig, [figDir 'normalized behavior for timing temp tuning curve'],fig_type);

%% ANALYSIS AND FIGURE: correlation between cumulative speed and sleep per fly

clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% 1) (Y) CUMULATIVE SLEEP vs (X) CUMULATIVE MOVEMENT PER FLY
fps = fly(1).parameters.FPS;
[dist, sleep] = deal(nan(num.trials,2)); % empty struct for fly speed & sleep values
% 
sexList = {'m','f'};
for i = 1:num.trials
    for sex = 1:2
        
        % total distance traveled
        speed = fly(i).(sexList{sex}).speed; % in mm/s
        dist(i,sex) = (sum(speed,'omitnan')/fps)/1000; % meters traveled over full experiment
        % total sleep
        raw_sleep = fly(i).(sexList{sex}).sleep; % logical of sleep over time
        sleep(i,sex) = (sum(raw_sleep)/fps)/60; % total sleep in minutes
        
    end
end


% Comparison of total sleep vs movement between males and females
sz = 75;
fig = getfig('',1,[400 647]); hold on
    hold on
    for sex = 1:2
        scatter(dist(:,sex),sleep(:,sex), sz, data.color(sex,:), 'filled')
    end
    xlabel('distance traveled (m)')
    ylabel('total sleep (min)')
formatFig(fig,blkbgd);

fig_title = 'total distance vs total sleep';
save_figure(fig, [figDir fig_title],fig_type);

% Quick stats on male vs female total sleep quantity
buff = 0.2;
fig = getfig('',1,[400 560]); hold on
for sex = 1:2
    x = linspace(sex-buff, sex+buff, num.trials);
    x = shuffle_data(x);
    y = sleep(:,sex);
    scatter(x,y,sz,data.color(sex,:), 'filled', 'MarkerFaceAlpha', 0.8)
    y_avg = mean(y, 'omitnan');
    plot([sex-(buff*1.5), sex+(buff*1.5)], [y_avg, y_avg], 'color',...
        data.color(sex,:), 'linewidth', 1.5)
end
formatFig(fig, blkbgd);
xlim([0, 3])
set(gca, 'xtick', 1:2, 'xticklabel', {'M','F'})
set(gca, 'xcolor', 'none')
ylabel('total sleep per fly (min)')
[h,p,~,stats] = ttest(sleep(:,1), sleep(:,2));
fprintf('Sleep: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% display a significance bar and star if the two are different...
fig_title = 'total sleep by sex';
save_figure(fig, [figDir fig_title],fig_type);


% Quick stats on male vs female total movememt quantity
buff = 0.2;
fig = getfig('',1,[400 560]); hold on
for sex = 1:2
    x = linspace(sex-buff, sex+buff, num.trials);
    x = shuffle_data(x);
    y = dist(:,sex);
    scatter(x,y,sz,data.color(sex,:), 'filled', 'MarkerFaceAlpha', 0.8)
    y_avg = mean(y, 'omitnan');
    plot([sex-(buff*1.5), sex+(buff*1.5)], [y_avg, y_avg], 'color',...
        data.color(sex,:), 'linewidth', 1.5)
end
formatFig(fig, blkbgd);
xlim([0, 3])
set(gca, 'xtick', 1:2, 'xticklabel', {'M','F'})
set(gca, 'xcolor', 'none')
ylabel('total distance traveled per fly (m)')
[h,p,~,stats] = ttest(sleep(:,1), sleep(:,2));
fprintf('Distance: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% display a significance bar and star if the two are different...
fig_title = 'total distance by sex';
save_figure(fig, [figDir fig_title],fig_type);

%% ANALYSIS: SLEEP FIGURES
% What is the relationship between avg speed and total distance traveled?

% what is the avg speed of flies in the period of time before they
% sleep? And does that correlate to how long they sleep during that bout?
clearvars('-except',initial_var{:})
fps = fly(1).fps;
foreColor = formattingColors(blkbgd); % get background colors
sleepData = struct; % new structure for the sleep data...

% 1) find the start time for each sleep bout
% 2) find the avg speed in the time bin leading up to the sleep
% 3) look at a range of timebins preceeding the sleep onset
% 4) check the relationship between avg speed vs sleep duration

% 0) find the onset of each sleep bout

sleep_buffer = fps; % buffer for sleep continuity
sleep = struct;
sleep(1).duration = []; % empty to start concatenating the new data
sleep(2).duration = []; % empty to start concatenating the new data

for trial = 1:num.trials
    for sex = 1:2
        x = data.sleep(:,sex,trial);
        sleepOnset = find(diff(x)==1)+1; % sleep start
        sleepOffset = find(diff(x)==-1)+1; % sleep end
        if ~isempty(sleepOffset)
            % find sleep bouts back to back that fall into the same bout with the
            % time between bouts under the sleep buffer 1/2 second limit
            merge_end = [((sleepOnset(2:end)-sleepOffset(1:end-1))<sleep_buffer);false];
            merge_loc = [false; merge_end(1:end-1)];
            if any(merge_loc)
                sleepOnset(merge_loc) = [];
                sleepOffset(merge_end) = [];
            end
            % % Quick visual check of the sleep start and stops
            % figure;
            % hold on 
            % plot(x)
            % v_line(sleepOnset,'green')
            % v_line(sleepOffset,'red')

            % sleep duration in frames->sec->mins
            sleepDuration = ((sleepOffset-sleepOnset)/fps)/60;
            fly(trial).data(sex).sleep.ON_OFFidx = [sleepOnset, sleepOffset];
            fly(trial).data(sex).sleep.duration = sleepDuration;
            sleep(sex).duration = [sleep(sex).duration; sleepDuration];
            sleep(sex).nBouts(trial) = length(sleepDuration);
        else
            fly(trial).data(sex).sleep.ON_OFFidx = [];
            fly(trial).data(sex).sleep.duration = [];
        end
    end
end
        
% FIGURE: quick look at the duration of sleep bouts between the sexes
% Quick stats on male vs female sleep duration
sz = 75;
buff = 0.2;
fig = getfig('',1,[400 560]); hold on
for sex = 1:2
    y = sleep(sex).duration;
    x = linspace(sex-buff, sex+buff, length(y));
    x = shuffle_data(x);
    scatter(x,y,sz,data.color(sex,:), 'filled', 'MarkerFaceAlpha', 0.8)
    y_avg = mean(y, 'omitnan');
    plot([sex-(buff*1.5), sex+(buff*1.5)], [y_avg, y_avg], 'color',...
        data.color(sex,:), 'linewidth', 1.5)
end
formatFig(fig, blkbgd);
xlim([0, 3])
set(gca, 'xtick', 1:2, 'xticklabel', {'M','F'})
set(gca, 'xcolor', 'none')
ylabel('sleep duration per bout (min)')
[h,p,~,stats] = ttest2(sleep(M).duration, sleep(F).duration);
fprintf('Sleep Duration: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% display a significance bar and star if the two are different...
fig_title = 'sleep bout duration by sex';
save_figure(fig, [figDir fig_title],fig_type);

% FIGURE: look at the number of sleep bouts in male and female flies
sz = 75;
buff = 0.2;
fig = getfig('',1,[400 560]); hold on
for sex = 1:2
    y = sleep(sex).nBouts;
    x = sex*ones([num.trials, 1]);
    x = linspace(sex-buff, sex+buff, length(y));
    x = shuffle_data(x);
    scatter(x,y,sz,data.color(sex,:), 'filled', 'MarkerFaceAlpha', 0.8,'xjitter', 'density')
    y_avg = mean(y, 'omitnan');
    plot([sex-(buff*1.5), sex+(buff*1.5)], [y_avg, y_avg], 'color',...
        data.color(sex,:), 'linewidth', 1.5)
end
formatFig(fig, blkbgd);
xlim([0, 3])
set(gca, 'xtick', 1:2, 'xticklabel', {'M','F'})
set(gca, 'xcolor', 'none')
ylabel('number of sleep bouts per fly')
[h,p,~,stats] = ttest2(sleep(M).nBouts, sleep(F).nBouts);
fprintf('Num Sleep Bouts: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% display a significance bar and star if the two are different...
fig_title = 'number of sleep bouts by sex';
save_figure(fig, [figDir fig_title],fig_type);

% FIGURE: sleep instances per temperature regime type
% pull the sleep data
types = {'WT', 'WS', 'CT', 'CS'}; % different temperature regimes
tb = data.tempbin;
for sex = 1:2
    % initialize empty variables for the data
    [y_avg, y_sem] = deal(nan([4,1]));
    y = nan([num.trials, 4]);
    for i = 1:4 % for each of the temp regime types
        roi = tb.(types{i});
        y_raw = squeeze(data.sleep(roi,sex,:));
        y(:,i) = (mean(y_raw,'omitnan')./fps)*100; % avg sleep (% of time)
        y_avg(i) = mean(y(:,i),'omitnan');
        y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
    end
    switch sex
        case 1
            m.y = y;
            m.y_raw = y_raw;
            m.y_avg = y_avg;
            m.y_sem = y_sem;
        case 2
            f.y = y;
            f.y_raw = y_raw;
            f.y_avg = y_avg;
            f.y_sem = y_sem;
    end
end

% plot the sleep data
% Create a bar plot for the average sleep data across temperature regimes
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
sz = 40;
sex_buff = 0.2;
[xM, xF] = deal(1:4);
xM = xM-sex_buff;
xF = xF+sex_buff;

figure; hold on;
% Male Sleep:
    bar(xM,m.y_avg, 'FaceColor', data.color(M,:),'barwidth', 0.4);
    errorbar(xM, m.y_avg, m.y_sem, 'color', foreColor, 'linestyle', 'none', 'LineWidth', 1.5);
    x = repmat(xM,[num.trials, 1]);
    scatter(x(:), m.y(:), sz, foreColor, 'filled', 'xjitter',...
        'density','MarkerFaceAlpha', 0.75,'xjitterwidth', 0.25,'markeredgecolor',data.color(M,:));

    bar(xF,f.y_avg, 'FaceColor', data.color(F,:),'barwidth', 0.4);
    errorbar(xF, f.y_avg, f.y_sem, 'color', foreColor, 'linestyle', 'none', 'LineWidth', 1.5);
    x = repmat(xF,[num.trials, 1]);
    scatter(x(:), f.y(:), sz, foreColor, 'filled', 'xjitter',...
        'density','MarkerFaceAlpha', 0.75,'xjitterwidth', 0.25,...
        'markeredgecolor', data.color(F,:));

    set(gca, 'XTickLabel', labels,'xtick', 1:4, 'XTickLabelRotation',40);
    ylabel('Sleep per minute (%)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
save_figure(gcf, [figDir 'sleep per temp regime'], fig_type);
    
% Stats: 

% Safe vs threat warm male
[h,p,~,stats] = ttest2(m.y(:,1), m.y(:,2)); % safe vs unsafe Warm Male
fprintf('Male unsafe vs safe Warm: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% Safe vs threat cold male
[h,p,~,stats] = ttest2(m.y(:,3), m.y(:,4)); % safe vs unsafe Warm Male
fprintf('Male unsafe vs safe Cold: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% Safe vs threat warm female
[h,p,~,stats] = ttest2(f.y(:,1), f.y(:,2)); % safe vs unsafe Warm Male
fprintf('Female unsafe vs safe Warm: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% Safe vs threat cold female
[h,p,~,stats] = ttest2(f.y(:,3), f.y(:,4)); % safe vs unsafe Warm Male
fprintf('Female unsafe vs safe Cold: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);


% TODO: collapse across hot and cold thermal threat regions
% FIGURE: sleep instances per temperature regime type
% pull the sleep data
tb = data.tempbin;
for sex = 1:2
    % initialize empty variables for the data
    [y_avg, y_sem] = deal(nan([2,1]));
    y = nan([num.trials, 2]);
    for i = 1:2 % for each of the temp regime types
        switch i
            case 1 % threat
                roi = replaceNaN(tb.CT, false) | replaceNaN(tb.WT, false);
            case 2 % safe
                roi = replaceNaN(tb.CS, false) | replaceNaN(tb.WS, false);
        end
        y_raw = squeeze(data.sleep(roi,sex,:));
        y(:,i) = (mean(y_raw,'omitnan')./fps)*100; % avg sleep (% of time)
        y_avg(i) = mean(y(:,i),'omitnan');
        y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
    end
    switch sex
        case 1
            m.y = y;
            m.y_raw = y_raw;
            m.y_avg = y_avg;
            m.y_sem = y_sem;
        case 2
            f.y = y;
            f.y_raw = y_raw;
            f.y_avg = y_avg;
            f.y_sem = y_sem;
    end
end

% plot the sleep data
% Create a bar plot for the average sleep data across temperature regimes
labels = {'threat', 'safe'}; % different temperature regimes
sz = 40;
sex_buff = 0.2;
[xM, xF] = deal([1,2]);
xM = xM-sex_buff;
xF = xF+sex_buff;

figure; hold on;
% Male Sleep:
    bar(xM,m.y_avg, 'FaceColor', data.color(M,:),'barwidth', 0.4);
    errorbar(xM, m.y_avg, m.y_sem, 'color', foreColor, 'linestyle', 'none', 'LineWidth', 1.5);
    x = repmat(xM,[num.trials, 1]);
    scatter(x(:), m.y(:), sz, foreColor, 'filled', 'xjitter',...
        'density','MarkerFaceAlpha', 0.75,'xjitterwidth', 0.25,'markeredgecolor',data.color(M,:));
% Female Sleep:
    bar(xF,f.y_avg, 'FaceColor', data.color(F,:),'barwidth', 0.4);
    errorbar(xF, f.y_avg, f.y_sem, 'color', foreColor, 'linestyle', 'none', 'LineWidth', 1.5);
    x = repmat(xF,[num.trials, 1]);
    scatter(x(:), f.y(:), sz, foreColor, 'filled', 'xjitter',...
        'density','MarkerFaceAlpha', 0.75,'xjitterwidth', 0.25,...
        'markeredgecolor', data.color(F,:));

    set(gca, 'XTickLabel', labels,'xtick', 1:2);
    ylabel('Sleep per minute (%)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
save_figure(gcf, [figDir 'sleep in safe vs threat'], fig_type);
    
% Stats: 

% Safe vs threat male
[h,p,~,stats] = ttest2(m.y(:,1), m.y(:,2)); 
fprintf('Male unsafe vs safe: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% Safe vs threat female
[h,p,~,stats] = ttest2(f.y(:,1), f.y(:,2)); 
fprintf('Female unsafe vs safe: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);


%% TODO: What are the rates of sleep given the immediately preceeding period of time for behavior?



%% TODO: What are the time delays from the start of a temperature change to each behavior?


%% TODO: Timecourse of behavior states for a given fly in a given trial
% What is the distribution of behavioral states within a temperature regime
% -- try a chord plot for each 
















