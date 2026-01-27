


%% TEMP TUNING CURVE
clearvars('-except',initial_var{:})

% Select the type of information to plot: 
[title_str,pName,y_dir,y_lab,nullD,scaler,dType,dir_end,sexSep,ylimits] = PlotParamSelectionHR('Heating and Cooling');
plot_err = true;

if isempty(title_str)
    return
end
fig_dir = [figDir, dir_end];
% set figure folder
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
 
% ------------- Plotting parameters -------------

% select temp protocol specific plotting features
autoYLim = false; % Y LIMITS
if any(isnan(ylimits)) % in case new data is added without specs for axis limits
    autoYLim = false;
end

if strcmp(groupName,'Berlin LTS caviar') % X LIMITS
    xlimits = [13, 37];
    autoXLim = false;
    xPos = [14.5, 25];  % for the shaded threat region on the plot
else
    autoXLim = true;
end

kolor = Color('vaporwavepurple'); % set color selection
r = 1; % figure rows
c = 2; % heating and cooling separated plots
LW = 2; % plotting line width
FA = 0.35; % SEM shading face alpha level

% ------------- DATA AND PLOTTING -------------

% pull larger data type group
yy = data.(pName);
if sexSep % data separated per fly or per group of flies
    y_all = [squeeze(yy(:,M,:)), squeeze(yy(:,F,:))];
else 
    y_all = yy;
end
x  = data.tempbin.temps; % temp bins
nTemps = length(x); % number of temperature bins
types  = {'cooling', 'warming'};

% Extract and Plot data:
fig = getfig('',1);
for ii = 1:2
    subplot(r,c,ii); hold on
        Idx = data.tempbin.(types{ii});
        rawY = nan([nTemps, 2]); % first col = avg, second = sem

        % extract the tuning information across the temp bins
        for tt = 1:nTemps 
            % all the cooling data across the flies that fits this temp bin
            y = yy(Idx(:,tt),:); 

            % fill the temp bin data into the appropriate structure
            raw = mean(y,'omitnan').*scaler; % find cooling fly data 
            rawY(tt,1) = mean(raw,'omitnan');
            rawY(tt,2) = sem(raw);
        end

        % plot the data: 
         plot_error_fills(plot_err, x, rawY(:,1), rawY(:,2), kolor, fig_type, FA);
         plot(x,rawY(:,1),'color', kolor, 'LineWidth', LW)
end
         
% formatting: 
formatFig(fig, blkbgd,[r,c]);
matchAxis(fig, true);
for ii = 1:2
    subplot(r,c,ii) 
    title(types{ii},'color', foreColor)
    xlabel('temperature (\circC)')
    if ~autoXLim
        xlim(xlimits)
    end
    if ~autoYLim
        ylim(ylimits)
    end

    % subplot specific adjustments
    if ii==1 % cooling
        set(gca, 'XDir','reverse')
        ylabel(y_lab)
    end
    if ii==2 % warming
        set(gca, 'YColor', 'none')
    end
end

% add shaded area for 'threat' temp region
for tt = 1:2
    subplot(r, c, tt)
    ylimits = ylim;
    pos = [xPos(1,tt), ylimits(1), range(xPos), range(ylimits)]; % [lower-left X, lower-left Y, X-width, Y-height]
    h = rectangle('Position', pos, ...
              'FaceColor', foreColor, ...   % RGB color
              'FaceAlpha', 0.1, ...
              'EdgeColor', 'none');
end

% Save the Figure
save_figure(fig, [fig_dir title_str ' tuning curve'])

%%


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