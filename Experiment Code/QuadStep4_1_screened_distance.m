
%% ANALYSIS: Screen out the area over the food well for comparison of tracking with glare errors

% Make folder for saving comparison data:
comp_saveDir = [ saveDir 'Food screen Comparision\'];
if ~exist(comp_saveDir,'dir')
    mkdir(comp_saveDir)
end
initial_vars{end+1} = 'comp_saveDir';

% Screen out flies that are within certain radius of the fly food (screen glare tracking errors from food)
for exp = 1:num.exp
    [screened_dist, FC_all, FC_screen] = deal(nan(size(grouped(exp).dist.all))); %set empty matrix for the screened distance to populate
    foodLocs = data(exp).T.foodLoc; % well with food or fictive food for all trials
    for trial = 1:num.trial(exp)
        % Pull fly coordinate position data for this trial from raw sleap tracking positions 
        % (those that have been screened into this arena, anyhow)
        X = data(exp).data(trial).data.x_loc;
        Y = data(exp).data(trial).data.y_loc;
        
        % Get center position for food well
        food_center = data(exp).data(trial).data.wellcenters(:,foodLocs(trial));
        % Find distance to well center for ALL flies within the arena
        dist2well = (sqrt((X-food_center(1)).^2 + (Y-food_center(2)).^2))./pix2mm; %convert to mm

        % full number of flies BEFORE screening out food well
        temp = sum(dist2well>0,2);
        FC_all(1:size(dist2well,1),trial) = temp; temp = [];

        % screen points within/on the food well
        loc = dist2well<=2.5; % 2.5mm since the wells are 5mm in diameter 
        dist2well(loc) = nan; % set tracked points within the screen to nan
        % calculate the mean distance to the food well if we screen out the food well
        screened_dist(1:size(dist2well,1),trial) = mean(dist2well,2,'omitnan');

        % number of flies AFTER screening out food well
        temp = sum(dist2well>0,2);
        FC_screen(1:size(dist2well,1),trial) = temp; temp = [];

    end
    
    % calculate the averages and errors for the new screened data:
    grouped(exp).screen_d.all = screened_dist;
    grouped(exp).screen_d.avg = mean(screened_dist,2,'omitnan');
    grouped(exp).screen_d.err = std(screened_dist,0,2,'omitnan');
    % all fly count avg & err
    grouped(exp).screen_d.FC_all.all = FC_all;
    grouped(exp).screen_d.FC_all.avg = mean(FC_all,2,'omitnan');
    grouped(exp).screen_d.FC_all.err = std(FC_all,0,2,'omitnan');
    % screen fly count avg & err
    grouped(exp).screen_d.FC_screen.all = FC_screen;
    grouped(exp).screen_d.FC_screen.avg = mean(FC_screen,2,'omitnan');
    grouped(exp).screen_d.FC_screen.err = std(FC_screen,0,2,'omitnan');

end
    
% Cluster the distances by temperature %TODO: update this for distance from
% food of the flies....
for i = 1:num.exp  
    temps = unique(data(i).G(1).TR.temps);
    rateIdx = data(i).G(1).TR.rateIdx;
    tempIdx = data(i).G(1).TR.tempIdx;
    % find rate index
    heatRate = find(data(i).G(1).TR.rates>0);
    coolRate = find(data(i).G(1).TR.rates<0);
    try 
        holdRate = find(data(i).G(1).TR.rates==0);
        ntypes = 3;
    catch
        ntypes = 2;
    end
    
    for temp = 1:length(temps)
        for type = 1:ntypes
            switch type
                case 1 %heating
                    g_name = 'increasing';
                    idxSex = heatRate;
                case 2 %cooling
                    g_name = 'decreasing';
                    idxSex = coolRate;
                case 3 %holding
                    g_name = 'holding';
                    idxSex = holdRate;
            end
            % take avg distance for each binned temp & rate
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            % distances
            grouped(i).screen_d.(g_name).avg(temp,1) = mean(mean(grouped(i).screen_d.all(loc,:),1,'omitnan'),'omitnan'); %avg 
            grouped(i).screen_d.(g_name).avg(temp,2) = std(mean(grouped(i).screen_d.all(loc,:),1,'omitnan'),'omitnan'); %err 
            grouped(i).screen_d.(g_name).all(temp,:) = mean(grouped(i).screen_d.all(loc,:),1,'omitnan');
 
            % full fly count
            grouped(i).screen_d.FC_all.(g_name).avg(temp,1) = mean(mean(grouped(i).screen_d.FC_all.all(loc,:),1,'omitnan'),'omitnan'); %avg 
            grouped(i).screen_d.FC_all.(g_name).avg(temp,2) = std(mean(grouped(i).screen_d.FC_all.all(loc,:),1,'omitnan'),'omitnan'); %err 
            grouped(i).screen_d.FC_all.(g_name).all(temp,:) = mean(grouped(i).screen_d.FC_all.all(loc,:),1,'omitnan');

            % screened fly count
            grouped(i).screen_d.FC_screen.(g_name).avg(temp,1) = mean(mean(grouped(i).screen_d.FC_screen.all(loc,:),1,'omitnan'),'omitnan'); %avg 
            grouped(i).screen_d.FC_screen.(g_name).avg(temp,2) = std(mean(grouped(i).screen_d.FC_screen.all(loc,:),1,'omitnan'),'omitnan'); %err 
            grouped(i).screen_d.FC_screen.(g_name).all(temp,:) = mean(grouped(i).screen_d.FC_screen.all(loc,:),1,'omitnan');

        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        grouped(i).screen_d.temp_all(temp,1) = mean(mean(grouped(i).screen_d.all(loc,:),1,'omitnan'),'omitnan'); %avg 
        grouped(i).screen_d.temp_all(temp,2) = std(mean(grouped(i).screen_d.all(loc,:),1,'omitnan'),'omitnan'); %err
        % full fly count
        grouped(i).screen_d.FC_all.temp_all(temp,1) = mean(mean(grouped(i).screen_d.FC_all.all(loc,:),1,'omitnan'),'omitnan'); %avg 
        grouped(i).screen_d.FC_all.temp_all(temp,2) = std(mean(grouped(i).screen_d.FC_all.all(loc,:),1,'omitnan'),'omitnan'); %err 
        % screened fly count
        grouped(i).screen_d.FC_screen.temp_all(temp,1) = mean(mean(grouped(i).screen_d.FC_screen.all(loc,:),1,'omitnan'),'omitnan'); %avg 
        grouped(i).screen_d.FC_screen.temp_all(temp,2) = std(mean(grouped(i).screen_d.FC_screen.all(loc,:),1,'omitnan'),'omitnan'); %err 
        
    end
    grouped(i).screen_d.temps = temps;
end
clearvars('-except',initial_vars{:})

%% FIGURES: distance to food screened vs full arena (group by group)
% food over time between all and screened data per experiment
clearvars('-except',initial_vars{:})
plot_err = true;
autoLim = true;
% Y limit ranges
dist_lim = [5,35];       %distance
dt_lim = [10,32];        %distance-temp
auto_time = true;      % automated time axis limits
time_lim = [0,400];     %time limit (x-axis)
nMax =  num.exp;%
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = 3:c:r*c; %binned distance alignment

LW = 0.75;
sSpan = 180;
dataString = cell([1,2]);

% FIGURES:
for i = 1:nMax
    fig = getfig('',true);
    x = grouped(i).time;
    kolor_1 = foreColor; % full data
    kolor_2 = Color('dodgerblue'); % screened data

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor_1)

    %distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).dist.avg,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor_1)
        y2 =smooth(grouped(i).screen_d.avg,'moving',sSpan);
        plot(x,y2,'LineWidth',LW,'Color',kolor_2)
        if ~autoLim
            ylim(dist_lim)
        end

    %temp dependent distance
    subplot(r,c,sb(3).idx); hold on
        x = grouped(i).dist.distavgbytemp(:,1);
        y = grouped(i).dist.distavgbytemp(:,2);
        y_err = grouped(i).dist.distavgbytemp_err(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        plot(x,y,'color',kolor_1,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor_1,  fig_type, 0.35);

        x = grouped(i).screen_d.temps;
        y = grouped(i).screen_d.temp_all(:,1);
        y_err = grouped(i).screen_d.temp_all(:,2)./num.trial(i); %flip to SEM like the unscreened data
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        plot(x,y,'color',kolor_2,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor_2,  fig_type, 0.35);

    dataString{1} = 'All data';
    dataString{2} = 'Screened data';

    
    % FORMATING AND LABELS
    formatFig(fig,blkbgd,[r,c],sb);
    % temp
    subplot(r,c,sb(1).idx)
    ylabel('\circC')
    set(gca,"XColor",backColor)
    if ~auto_time
        xlim(time_lim)
    end
    % distance
    subplot(r,c,sb(2).idx)
    ylabel('proximity to food (mm)')
    xlabel('time (min)')
    set(gca,'ydir','reverse')
    if ~auto_time
        xlim(time_lim)
    end
    % temp-distance relationship
    subplot(r,c,sb(3).idx)
    ylabel('proximity to food (mm)')
    xlabel('temp (\circC)')
    if ~autoLim
        ylim(dt_lim)
    end
    h_line(18.1,'grey',':',1) %36.2
    set(gca,'ydir','reverse')
    %
    legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 7)
    
    % save figure
    save_figure(fig,[comp_saveDir grouped(i).name  ' timecourse summary'],fig_type,true,true,'-r80');

end

%% FIGURE: Comparison of tracked flies over time between screened and unscreened data

clearvars('-except',initial_vars{:})
plot_err = true;
autoLim = true;
% Y limit ranges
dist_lim = [5,35];       %distance
dt_lim = [10,32];        %distance-temp
auto_time = true;      % automated time axis limits
time_lim = [0,400];     %time limit (x-axis)
nMax =  num.exp;%
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = 3:c:r*c; %binned distance alignment

LW = 0.75;
sSpan = 180;
dataString = cell([1,2]);

% FIGURES:
for i = 1:nMax
    fig = getfig('',true);
    x = grouped(i).time;
    kolor_1 = foreColor; % full data
    kolor_2 = Color('dodgerblue'); % screened data

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor_1)

    %number of flies
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(i).screen_d.FC_all.avg,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor_1)
        y =smooth(grouped(i).screen_d.FC_screen.avg,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor_2)
        if ~autoLim
            ylim(dist_lim)
        end
        h_line(15,'r','-',2)

    %temp dependent fly count
    subplot(r,c,sb(3).idx); hold on
        x = grouped(i).screen_d.temps;
        y = grouped(i).screen_d.FC_all.temp_all(:,1);
        y_err = grouped(i).screen_d.FC_all.temp_all(:,2)./num.trial(i);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        plot(x,y,'color',kolor_1,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor_1,  fig_type, 0.35);

        x = grouped(i).screen_d.temps;
        y = grouped(i).screen_d.FC_screen.temp_all(:,1);
        y_err = grouped(i).screen_d.FC_screen.temp_all(:,2)./num.trial(i);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        plot(x,y,'color',kolor_2,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor_2,  fig_type, 0.35);

    dataString{1} = 'All data';
    dataString{2} = 'Screened data';

    
    % FORMATING AND LABELS
    formatFig(fig,blkbgd,[r,c],sb);
    % temp
    subplot(r,c,sb(1).idx)
    ylabel('\circC')
    set(gca,"XColor",backColor)
    if ~auto_time
        xlim(time_lim)
    end
    % fly count over time
    subplot(r,c,sb(2).idx)
    ylabel('avg fly count')
    xlabel('time (min)')
    if ~auto_time
        xlim(time_lim)
    end
    % temp-fly count relationship
    subplot(r,c,sb(3).idx)
    ylabel('avg fly count')
    xlabel('temp (\circC)')
    if ~autoLim
        ylim(dt_lim)
    end
    h_line(15,'grey',':',1) 
    %
    legend(dataString,'textcolor',foreColor, 'location', 'northeast', 'box', 'off','fontsize', 7)
    
    % save figure
    save_figure(fig,[comp_saveDir grouped(i).name  ' flycount over time summary'],fig_type,true,true,'-r80');

end

%% FIGURE: Histogram of fly count overlay -- before and after screening
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd); %get background colors
kolor_1 = foreColor;
kolor_2 = Color('dodgerblue');


for i = 1:num.exp
    fig = getfig('', 1);

    histo_data = grouped(i).screen_d.FC_all.all(:);
    histo_screened = grouped(i).screen_d.FC_screen.all(:);

    h = histogram(histo_data,'FaceColor',kolor_1,'FaceAlpha',0.75) ;
    hold on
    histogram(histo_screened,h.BinEdges,'FaceColor',kolor_2,'FaceAlpha',0.75) ;
    v_line(15,'r','-',2)
    xlabel('fly count')
    ylabel('frame count')

    formatFig(fig, blkbgd);
    legend({'Raw','Screened'},'textcolor', foreColor, 'location', 'northeast', 'box', 'off','fontsize',12)
    title_str = strrep(grouped(i).name,'_',' ');
    title(title_str,'color',foreColor)
    
    save_figure(fig,[comp_saveDir title_str  ' fly count histogram'],fig_type,true,true,'-r80');
end


%% FIGURE & STATS: cumulative hysteresis for each genotype / trial
clearvars('-except',initial_vars{:})
LW = 0.75;
buff = 0.2;
SZ = 50;
r = 1; %rows
c = 3; %columns
plot_err = false;
plotSig = true; %plot significance stars
[foreColor,backColor] = formattingColors(blkbgd);
kolor_1 = Color('white'); % full data
kolor_2 = Color('dodgerblue'); % screened data

for i = 1:num.exp

    % FIGURE:
    fig = getfig('',true);
    
    % Hystersis
    subplot(r,c,1)
        hold on
        % ----- unscreened -------
        %increasing
        x = grouped(i).increasing.temps;
        y = grouped(i).increasing.avg;
        y_err = grouped(i).increasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];
        plot_error_fills(plot_err, x, y, y_err, kolor_1,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW+0.5,'Color',kolor_1,'linestyle','-')
        %decreasing
        x = grouped(i).decreasing.temps;
        y = grouped(i).decreasing.avg;
        y_err = grouped(i).decreasing.err;
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];
        plot_error_fills(plot_err, x, y, y_err, kolor_1,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW+.5,'Color',kolor_1,'linestyle','--','HandleVisibility','off');
        % ----- screened -------
        %increasing
        x = grouped(i).screen_d.temps;
        y = grouped(i).screen_d.increasing.avg(:,1);
        y_err = grouped(i).screen_d.increasing.avg(:,2);
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];
        plot_error_fills(plot_err, x, y, y_err, kolor_2,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW+0.5,'Color',kolor_2,'linestyle','-')
        %decreasing
        x = grouped(i).screen_d.temps;
        y = grouped(i).screen_d.decreasing.avg(:,1);
        y_err = grouped(i).screen_d.decreasing.avg(:,2);
        loc = isnan(y) | isnan(y_err);% remove nans
        y(loc) = []; x(loc) = []; y_err(loc) = [];
        plot_error_fills(plot_err, x, y, y_err, kolor_2,  fig_type, 0.2);
        plot(x,y,'LineWidth',LW+.5,'Color',kolor_2,'linestyle','--','HandleVisibility','off');

        ylabel('proximity to food (mm)')
        xlabel('temp (\circC)')
        set(gca, 'ydir', 'reverse')

    % Pull difference in distance heating-cooling
    subplot(r,c,2)
        hold on
        
        x = repmat(grouped(i).decreasing.temps,[1,num.trial(i)]);
        y = grouped(i).decreasing.all-grouped(i).increasing.all;
        plot(mean(x,2),mean(y,2),'color',kolor_1,'LineWidth',2)
        
        x = repmat(grouped(i).screen_d.temps',[1,num.trial(i)]);
        y = grouped(i).screen_d.decreasing.all - grouped(i).screen_d.increasing.all;
        plot(mean(x,2),mean(y,2),'color',kolor_2,'LineWidth',2)

        h_line(0,foreColor,':',1)
        xlabel('temp (\circC)')
        ylabel('distance difference (mm)')

% Cumulative difference in proximity
    subplot(r,c,3)
        hold on
        ii = 1;
        y = grouped(i).decreasing.all-grouped(i).increasing.all;
        plotY = sum(y,1,'omitnan');
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        scatter(x,plotY,SZ,kolor_1,"filled","o")
        plot([ii-buff,ii+buff],[mean(plotY),mean(plotY)],'color',foreColor,'LineWidth',2)
        
        ii = 2;
        y = grouped(i).screen_d.decreasing.all - grouped(i).screen_d.increasing.all;
        plotY = sum(y,1,'omitnan');
        x = shuffle_data(linspace(ii-buff,ii+buff,num.trial(i)));
        scatter(x,plotY,SZ,kolor_2,"filled","o")
        plot([ii-buff,ii+buff],[mean(plotY),mean(plotY)],'color',foreColor,'LineWidth',2)

        xlim([0.5,num.exp+0.5])
        h_line(0,foreColor,':',1)
        ylabel('cumulative difference (mm)')
        
        formatFig(fig,blkbgd,[r,c]);
        set(gca,'XTick',[],'xcolor',backColor)
        xlabel('Group','color',foreColor)

% % STATS: are the means of any groups different from zero?
% [p, mlt, id] = deal([]);
% for ii = 1:num.exp
%     i = expOrder(ii);
%     y = grouped(i).decreasing.all-grouped(i).increasing.all;
%     plotY = sum(y,1,'omitnan');
%     [~,p(ii)] = ttest(plotY);
%     group_name{ii} = expNames{i};
%     %multicompare
%     mlt = autoCat(mlt, plotY',false);
%     id = autoCat(id,i*ones(length(plotY),1),false);
% end
% %Bonferonni correction:
% alpha = 0.05;
% m = num.exp;
% p_limit = alpha/m;
% h = p<=p_limit;
% stats_tbl = table(group_name',h',p','VariableNames',{'group','significant','p value'});
% disp(stats_tbl)
% 
% % add significance stars to the figure:
% if plotSig
%     y_pos = rangeLine(fig,1);
%     subplot(r,c,3); hold on
%     for ii = 1:num.exp
%         if h(ii)
%             scatter(ii,y_pos,100,foreColor,'*')
%         end
%     end
% end
% 
% % Multicompare across the groups for significance
% % STATS:
% % TODO -- update all other stats to reflect this vv
% % determine which groups differ from each other
% [~,~,stats] = anova1(mlt(:),id(:),'off');
% alpha = 0.05; %significance level
% [c,~,~,~] = multcompare(stats,alpha,'off');
% % bonferonni multiple comparisons correction
% m = size(c,1); %number of hypotheses
% sigThreshold = alpha/m;
% %find p-values that fall under the threshold
% significantHypotheses = c(:,6)<=sigThreshold;
% fprintf('\n\nPosition hysteresis cross group comparison statistics\n\n')
% [Group1,Group2,P_Value] = deal([]);
% idx = 0;
% for i = 1:length(significantHypotheses)
%     if significantHypotheses(i)
%         idx = idx+1;
%         Group1{idx,1} = expNames{c(i,1)};
%         Group2{idx,1} = expNames{c(i,2)};
%         P_Value(idx,1) = c(i,6);
%     end
% end
% sig_comp = table(Group1,Group2,P_Value);
% disp(sig_comp)

% save figure
    save_figure(fig,[comp_saveDir grouped(i).name  ' hysteresis summary'],fig_type,true,true,'-r80');

end

%%  
%% FIGURE: Distance from food: basic over-lap of time-trials and temperature protocols 
clearvars('-except',initial_vars{:})
plot_err = true;
autoLim = false;
% Y limit ranges
dist_lim = [5,35];       %distance
dt_lim = [14,32];        %distance-temp
auto_time = true;      % automated time axis limits
time_lim = [0,400];     %time limit (x-axis)
nMax =  num.exp;%
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = [1,2]; %temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = 3:c:r*c; %binned distance alignment

LW = 0.75;
sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true);
for i = 1:nMax
%     i = expOrder(ii);
    x = grouped(i).time;
    kolor = grouped(i).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(i).temp;
        plot(x,y,'LineWidth',2,'Color',kolor)

    %distance
    subplot(r,c,sb(2).idx); hold on
        y =smooth(grouped(i).screen_d.avg,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        if ~autoLim
            ylim(dist_lim)
        end

    %temp dependent distance
    subplot(r,c,sb(3).idx); hold on
        x = grouped(i).screen_d.temps;
        y = grouped(i).screen_d.temp_all(:,1);
        y_err = grouped(i).screen_d.temp_all(:,2)./num.trial(i); %flip to SEM like the unscreened data
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];

        plot(x,y,'color',kolor,'linewidth',LW+1)
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
        dataString{i} = grouped(i).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",backColor)
if ~auto_time
    xlim(time_lim)
end
% distance
subplot(r,c,sb(2).idx)
ylabel('proximity to food (mm)')
xlabel('time (min)')
set(gca,'ydir','reverse')
if ~auto_time
    xlim(time_lim)
end
% temp-distance relationship
subplot(r,c,sb(3).idx)
ylabel('proximity to food (mm)')
xlabel('temp (\circC)')
if ~autoLim
    ylim(dt_lim)
end
h_line(18.1,'grey',':',1) %36.2
set(gca,'ydir','reverse')
%
legend(dataString,'textcolor', foreColor, 'location', 'southeast', 'box', 'off','fontsize', 5)


% save figure
save_figure(fig,[comp_saveDir expGroup ' distance to food over time'],fig_type);


% fig_type = '-pdf'; autoLim = false;


%% Fig for easy overlay with other data later TODO
clearvars('-except',initial_vars{:})
blkbgd = false;
fig_type = '-pdf'; 
buff = 0.1;
sz = 50;
autoLim = false;
dist_lim = [10,32];
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

dataString = cell([1,num.exp]);

% Find the max temp and min temp of all the experiments
temp_min = 16; 
temp_max = 26;
temp_bin = 0.5;
cMapRange = temp_min:temp_bin:temp_max;
ntemps = length(cMapRange);
% color map information
g1 = floor(ntemps/2);
g2 = ntemps-g1;
g1_cMap = Color('deepskyblue','grey',g1); %deepskyblue
g2_cMap = Color('grey','red',g2);
cMap = [g1_cMap;g2_cMap];

% FIGURE:
fig = getfig('',true,[565 649]);
for ii = 1:num.exp
    hold on
    i = expOrder(ii);
    % get x and y values
    y = grouped(i).screen_d.temp_all(:,1);
    x = shuffle_data(linspace(i-buff,i+buff, length(y)));
    % x = i*ones(size(y));
    loc = isnan(y);
    x(loc) = [];
    y(loc) = [];
    

    % get color map for the temperatures in this experiment
    CC = grouped(i).screen_d.temps;
    color_idx = discretize(CC,cMapRange);
    kolor = cMap(color_idx,:);
    kolor(loc,:) = [];

    scatter(x,y,sz,kolor,'filled')

    % dataString{ii} = strrep(data(i).foodNames{1},'_',' ');
    dataString{ii} = strrep(grouped(i).name,'_',' ');

end

% FORMATING AND LABELS
formatFig(fig,blkbgd);
xlim([0,num.exp+1])
ax = gca;
set(ax,'YDir','reverse')
if ~autoLim
    ylim(dist_lim)
end
h_line(18.1,foreColor)
ylabel('proximity to food (mm)')

set(ax,'XTick',1:num.exp)%,'XTickLabel',dataString,'XTickLabelRotation',0

% Setup color code map:
colormap(cMap)
h = colorbar;
h.Color = foreColor;
N = findLargestDivisor(temp_min, temp_max, 5);
temp_label = (temp_min:N:temp_max);
h.Ticks  = linspace(0,1,length(temp_label));
h.TickLabels = temp_label;
h.Label.String = 'Temp (\circC)';
h.Label.Color = foreColor;

% Update group labels to match color scheme 
xTicks = ax.XTick;% Get current tick positions and labels
xTickLabels = string(ax.XTickLabel); % Convert to string array for easier manipulation
ax.XTickLabel = [];% Hide original tick labels
% Create new tick labels with specified colors
for ii = 1:length(xTicks)
    i = expOrder(ii);
    text(xTicks(i), ax.YLim(2), xTickLabels(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', grouped(i).color,'FontSize',16);
end
ax.XColor = backColor;
save_figure(fig,[comp_saveDir expGroup ' distance to food vertical'],fig_type);

%% FIGURE: heating and cooling separated vertical temp colored distance to food

clearvars('-except',initial_vars{:})
% blkbgd = true;  fig_type = '-png'; 
fig_type = '-pdf'; blkbgd = false;
buff = 0.1;
sz = 50;
autoLim = false;
dist_lim = [10,32];
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

dataString = cell([1,num.exp]);

% Find the max temp and min temp of all the experiments
temp_min = 16; 
temp_max = 26;
temp_bin = 0.5;
cMapRange = temp_min:temp_bin:temp_max;
ntemps = length(cMapRange);
% color map information
g1 = floor(ntemps/2);
g2 = ntemps-g1;
g1_cMap = Color('deepskyblue','grey',g1); %deepskyblue
g2_cMap = Color('grey','red',g2);
cMap = [g1_cMap;g2_cMap];

% FIGURE:
fig = getfig('',true,[565 649]);
for ii = 1:num.exp
    hold on
    i = expOrder(ii);
    % get color map for the temperatures in this experiment
    CC = grouped(i).screen_d.temps;
    color_idx = discretize(CC,cMapRange);

    % cooling
    y = grouped(i).screen_d.decreasing.avg(:,1);
    x = shuffle_data(linspace(ii-buff*1.5,ii-buff/2, length(y)))';
    loc = isnan(y);
    x(loc) = [];
    y(loc) = [];
    kolor = cMap(color_idx,:);
    kolor(loc,:) = [];
    scatter(x,y,sz,kolor,'LineWidth',1)

    % warming
    y = grouped(i).screen_d.increasing.avg(:,1);
    x = shuffle_data(linspace(ii+buff/2,ii+buff*1.5, length(y)))';
    loc = isnan(y);
    x(loc) = [];
    y(loc) = [];
    kolor = cMap(color_idx,:);
    kolor(loc,:) = [];
    scatter(x,y,sz,kolor,'filled')

    % dataString{ii} = strrep(data(i).foodNames{1},'_',' ');
    dataString{ii} = strrep(grouped(i).name,'_',' ');

end

% FORMATING AND LABELS
formatFig(fig,blkbgd);
xlim([0,9])
ax = gca;
set(ax,'YDir','reverse')
if ~autoLim
    ylim(dist_lim)
end
h_line(18.1,foreColor)
ylabel('proximity to food (mm)')

set(ax,'XTick',1:num.exp)%,'XTickLabel',dataString,'XTickLabelRotation',0

% Setup color code map:
colormap(cMap)
h = colorbar;
h.Color = foreColor;
N = findLargestDivisor(temp_min, temp_max, 5);
temp_label = (temp_min:N:temp_max);
h.Ticks  = linspace(0,1,length(temp_label));
h.TickLabels = temp_label;
h.Label.String = 'Temp (\circC)';
h.Label.Color = foreColor;

% Update group labels to match color scheme 
xTicks = ax.XTick;% Get current tick positions and labels
xTickLabels = string(ax.XTickLabel); % Convert to string array for easier manipulation
ax.XTickLabel = [];% Hide original tick labels
% Create new tick labels with specified colors
for ii = 1:length(xTicks)
    i = expOrder(ii);
    text(xTicks(ii), ax.YLim(2), xTickLabels(ii), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', grouped(i).color,'FontSize',16);
end
ax.XColor = backColor;

save_figure(fig,[comp_saveDir expGroup ' distance to food vertical warm cool split'],fig_type);

%% FIGURE:  distance to food as scatter plot (for simple view)
clearvars('-except',initial_vars{:})
% blkbgd = true;  fig_type = '-png'; 
fig_type = '-pdf'; blkbgd = false;
buff = 0.2;
sz = 50;
autoLim = false;
dist_lim = [10,32];
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

dataString = cell([1,num.exp]);

% temperatures to compare
temp = 17:1:25;

% FIGURE:
fig = getfig('',true,[565 649]);
for ii = 1:num.exp
    hold on
    i = expOrder(ii);
    
    % find index for the selected temperatures
    for tt = 1:length(temp)
        [~, idx(tt)] = min(abs(grouped(i).screen_d.temps-temp(tt)));
        % cooling
        y = grouped(i).screen_d.decreasing.all(idx(tt),:);
        % x = shuffle_data(linspace(temp(tt)-buff, temp(tt)+buff, length(y)));
        x_roi = [tt-buff-0.1, tt-0.1];
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y)));
        scatter(x,y,sz, grouped(i).color)
        plot(x_roi,[mean(y),mean(y)], 'linewidth', 1, 'color', grouped(i).color,'LineStyle',':')
        % warming
        y = grouped(i).screen_d.increasing.all(idx(tt),:);
        % x = shuffle_data(linspace(temp(tt)-buff, temp(tt)+buff, length(y)));
        x_roi = [tt+0.1, tt+buff+0.1];
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y)));
        scatter(x,y,sz, grouped(i).color,'filled')
        plot(x_roi,[mean(y),mean(y)], 'linewidth', 1, 'color', grouped(i).color)
    end
end
    
    


% FIGURE: plot a comparison line between the heating and cooling position for each temp
fig = getfig('',true,[565 649]);
for ii = 1:num.exp
    hold on
    i = expOrder(ii);
    
    % find index for the selected temperatures
    for tt = 1:length(temp)
        [~, idx(tt)] = min(abs(grouped(i).screen_d.temps-temp(tt)));
        % cooling

        y1 = grouped(i).screen_d.decreasing.all(idx(tt),:);
        % x = shuffle_data(linspace(temp(tt)-buff, temp(tt)+buff, length(y)));
        x1_roi = [tt-buff-0.1, tt-0.1];
        x1 = shuffle_data(linspace(x1_roi(1), x1_roi(2), length(y)));
        scatter(x1,y1,sz, grouped(i).color)
        % plot(x1_roi,[mean(y2),mean(y2)], 'linewidth', 1, 'color', grouped(i).color,'LineStyle',':')
        
        % warming
        y2 = grouped(i).screen_d.increasing.all(idx(tt),:);
        % x = shuffle_data(linspace(temp(tt)-buff, temp(tt)+buff, length(y)));
        x2_roi = [tt+0.1, tt+buff+0.1];
        x2 = shuffle_data(linspace(x2_roi(1), x2_roi(2), length(y)));
        scatter(x2,y2,sz, grouped(i).color,'filled')
        % plot(x2_roi,[mean(y2),mean(y2)], 'linewidth', 1, 'color', grouped(i).color)

        plot([x1;x2], [y1;y2],'color',grouped(i).color)
    
    end
end

%% FIXED FIGURE [static vs dynamic temp] WORKING

% make a comparison for fixed temperature/temp-rate points to compare the static
% temperature protocols and the dynamic temperature holds
clearvars('-except',initial_vars{:})
% blkbgd = true;  fig_type = '-png'; 
fig_type = '-pdf'; blkbgd = false;
buff = 0.2;
sz = 50; 
autoLim = false;
dist_lim = [10,32];
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

temp_comparisions = [15, 17, 20, 23, 25]; % these are the temps with fixed, held points -- might do just 17:25 for easy comp.

% find the average temp

% TODO HERE: just plot the average for each data group at the time-locked
% version (grouped(i).increasing)...



% FIGURE: plot a comparison line between the heating and cooling position for each temp
fig = getfig('',true,[565 649]);
for ii = 1:num.exp
    hold on
    i = expOrder(ii);
    
    % find index for the selected temperatures
    for tt = 1:length(temp)
        [~, idx(tt)] = min(abs(grouped(i).screen_d.temps-temp(tt)));
        % cooling

        y1 = grouped(i).screen_d.decreasing.all(idx(tt),:);
        % x = shuffle_data(linspace(temp(tt)-buff, temp(tt)+buff, length(y)));
        x1_roi = [tt-buff-0.1, tt-0.1];
        x1 = shuffle_data(linspace(x1_roi(1), x1_roi(2), length(y)));
        scatter(x1,y1,sz, grouped(i).color)
        % plot(x1_roi,[mean(y2),mean(y2)], 'linewidth', 1, 'color', grouped(i).color,'LineStyle',':')
        
        % warming
        y2 = grouped(i).screen_d.increasing.all(idx(tt),:);
        % x = shuffle_data(linspace(temp(tt)-buff, temp(tt)+buff, length(y)));
        x2_roi = [tt+0.1, tt+buff+0.1];
        x2 = shuffle_data(linspace(x2_roi(1), x2_roi(2), length(y)));
        scatter(x2,y2,sz, grouped(i).color,'filled')
        % plot(x2_roi,[mean(y2),mean(y2)], 'linewidth', 1, 'color', grouped(i).color)

        plot([x1;x2], [y1;y2],'color',grouped(i).color)
    
    end
end









































