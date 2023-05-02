

% NULL DISTRIBUTION ANALYSIS /  PROCESSING

%% Near inifinite place null distribution
clearvars('-except',initial_vars{:})

% NULL DISTRIBUTION (this should be identical for all trials/experiments in the same arena)
% Generate the null distribution of distances to the food well:
pixel_buffer = 50; % edge of arena pixel buffer
i = 1;
trial = 1;

% Set axis limits for the selected arena
c1 = data(i).data(trial).data.centre(1);
c2 = data(i).data(trial).data.centre(2);
r = data(i).data(trial).data.r;
xlimit = [c1-(r+pixel_buffer),c1+(r+pixel_buffer)];
ylimit = [c2-(r+pixel_buffer),c2+pixel_buffer+r];
foodWellLoc = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));

nbins = 10000;

% find the 'auto bin' lines
xedge = linspace(xlimit(1),xlimit(2),nbins+1);
yedge = linspace(ylimit(1),ylimit(2),nbins+1);
X = []; Y = [];
idx = 1;
for y_loc = 1:nbins
    for x_loc = 1:nbins
        X(idx) = xedge(x_loc);
        Y(idx) = yedge(y_loc);
        idx = idx + 1;
    end
end
% screen out the units outside the arena circle
temp_dist = sqrt((X-c1).^2 + (Y-c2).^2);
loc = temp_dist>r;
X(loc) = [];
Y(loc) = [];

% find food well  distance
null_distance = sqrt((X-foodWellLoc(1)).^2 + (Y-foodWellLoc(2)).^2)./pix2mm;

binSpacing = 0.10; %mm spacing for bins
binedges = 0:binSpacing:55; 

% find probability of occupancy for each bin from the null distribution
N = discretize(null_distance,binedges);
N_tot = length(null_distance);
prob = zeros(length(binedges),1);
for i = 1:length(binedges)
    prob(i) = (sum(N==i))/N_tot;
end

null_dist = struct;
null_dist.distance = null_distance;
null_dist.binwidth = binSpacing;
null_dist.binedges = binedges';
null_dist.prob = prob;

% sampling density
arena_area = pi*((r/pix2mm)^2); %mm^2
sampling_density = length(null_distance)/arena_area;

initial_vars{end+1} = 'null_dist';

% Save Null_Distribution
save([baseFolder 'Fundamentals\Distance_to_well_distribution.mat'],'null_dist','-v7.3');



%% Find probability of occupancy for each distance
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);
LW = 3;

fig = getfig('',1); hold on
    yyaxis left
    histogram(null_dist.distance,null_dist.binedges,'FaceColor',Color('grey'),'FaceAlpha',1)
    yyaxis right
    plot(null_dist.binedges,null_dist.prob,'color', Color('teal'),'linewidth',LW,'linestyle','-')
    formatFig(fig,blkbgd);
    axis tight
    yyaxis left
    ylabel('count','FontSize',24)
    yyaxis right
    xlabel('distance to well (mm)','FontSize',24)
    ylabel('occupancy probability','FontSize',24)
    yyaxis left
    set(gca,'Ycolor',foreColor)

% Save null distribution probability
save_figure(fig,[baseFolder 'Fundamentals\Final distance probabilty null distribution'],fig_type);

%% FIGURE: Cumulative distribution of sleeping distances to food
clearvars('-except',initial_vars{:})
LW = 2;
[~,backColor] = formattingColors(blkbgd);

% FIGURE:
fig = getfig('',1,[687 680]);  
hold on
% NULL DISTRIBUTION:
null_CDF = cdfplot(null_dist.distance);
null_CDF.Color = Color('teal');
null_CDF.LineWidth = 3; 

formatFig(fig,blkbgd);
xlabel('distance to well (mm)','FontSize',20)
ylabel('Empirical Cumulative Distribution','FontSize',20)
set(gca,'TickDir','out')
set(gca,'GridColor',backColor)
title('')

% save figure
save_figure(fig,[baseFolder 'Fundamentals\Final CDF of null distribution'],fig_type);

%% FIGURE: Weighted distance to food timecourse
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);

fig_dir = [saveDir 'Null Distribution\'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end 

% scale distance by probability wieght *divide by probability?
binedges = null_dist.binedges;
prob = null_dist.prob;
scaled_dist = struct;
for i = 1:num.exp
    scaledY = struct;
    all_dist = [];
    for trial = 1:num.trial(i)      
        y = data(i).data(trial).occupancy.dist2wells(:,data(i).T.foodLoc(trial));
        new_y = nan(size(y));
        for tt = 1:length(binedges)-1 %loop through the 551 bins
            loc = (y >= binedges(tt)) & (y < binedges(tt+1));
            new_y(loc) = y(loc)/prob(tt);
        end
         scaledY(trial).dist = new_y;
         all_dist = autoCat(all_dist,new_y,false);
    end
    scaled_dist(i).data = scaledY;
    scaled_dist(i).all_dist = all_dist;
end

LW = 0.25;
sSpan = 1;

fig = getfig('',1);
for i = 1:num.exp
        subplot(2,1,1); hold on 
            y = scaled_dist(i).all_dist;
            y_avg = mean(y,2,'omitnan');
            x = grouped(i).time;
            plot(x,smooth(y_avg,sSpan,'moving'),'color',grouped(i).color,'linewidth', LW)
        subplot(2,1,2); hold on 
            y = grouped(i).dist.avg;
            plot(x,smooth(y,sSpan,'moving'),'color',grouped(i).color,'linewidth', LW)
end
formatFig(fig,blkbgd,[2,1]);
% save figure
save_figure(fig,[fig_dir  'norm and wieghed distances'],fig_type);

LW = 0.25;
sSpan = 1;
for i = 1:num.exp
    fig = getfig('',1);
    hold on 
        x = grouped(i).time;
        % OG
        y = grouped(i).dist.avg;
        plot(x,smooth(y,sSpan,'moving'),'color',foreColor,'linewidth', LW)
        % scaled
        yyaxis right
        y = scaled_dist(i).all_dist;
        y_avg = mean(y,2,'omitnan');
        plot(x,smooth(y_avg,sSpan,'moving'),'color',Color('red'),'linewidth', LW)        
    % Formatting
    formatFig(fig,blkbgd);
    yyaxis right
    set(gca,'ydir','reverse','YColor', 'r')
    ylabel('Weighted distance')
    yyaxis left
    set(gca,'ycolor',foreColor,'ydir', 'reverse')
    ylabel('Proximity to food (mm)')
    xlabel('time (min)')
    % save figure
    save_figure(fig,[fig_dir  expNames{i} ' norm and wieghed distances'],fig_type);

end

%%  FIGURES: temp binned weighted distances to food 

% Cluster the sleeping flies by temperature
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
            %fraction of flies sleeping
            loc = rateIdx==idxSex & tempIdx==temp; %rate and temp align
            scaled_dist(i).(g_name)(temp,1) = mean(mean(scaled_dist(i).all_dist(loc,:),2,'omitnan'),'omitnan'); %avg 
            scaled_dist(i).(g_name)(temp,2) = std(mean(scaled_dist(i).all_dist(loc,:),1,'omitnan'),'omitnan');%./num.trial(i); %err  
        end
        % Clustered by temp (regardless of heating/cooling)
        loc = tempIdx==temp; %temp align only
        scaled_dist(i).temp_all(temp,1) = mean(mean(scaled_dist(i).all_dist(loc,:),2,'omitnan'),'omitnan'); %avg 
        scaled_dist(i).temp_all(temp,2) = std(mean(scaled_dist(i).all_dist(loc,:),1,'omitnan'),'omitnan')./num.trial(i);% %err
    end
    scaled_dist(i).temps = temps;
end

plot_err = true;
equalLim = true;
LW = 1.5;
r = 1;
c = 2;
nMax = num.exp;

% FIGURE:
fig = getfig('',true); 
% AVG
subplot(r,c,1)
hold on
for i = 1:nMax
    kolor = grouped(i).color;
    x =scaled_dist(i).temps;
    y = scaled_dist(i).temp_all(:,1);
    y_err = scaled_dist(i).temp_all(:,2);
    loc = isnan(y)|isnan(y_err);
    x(loc) = [];
    y(loc) = [];
    y_err(loc) = [];

    plot(x,y,'color',kolor,'linewidth',LW+1)
    if plot_err
        fill_data = error_fill(x, y, y_err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
        set(h, 'facealpha', 0.35)
    end
end
set(gca, 'ydir', 'reverse')
xlabel('Temperature (\circC)')
ylabel('weighted proximity to well')
% SEP HEAT / COOL
subplot(r,c,2)
hold on
for i = 1:nMax
    kolor = grouped(i).color;   
    for type = 1:2 %increasing | decreasing 
        switch type
            case 1
                section_type = 'increasing';
                l_style = '-';
            case 2
                section_type = 'decreasing';
                l_style = '--';
        end
        x = scaled_dist(i).temps;
        y = scaled_dist(i).(section_type)(:,1);
        y_err = scaled_dist(i).(section_type)(:,2);
        loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
%         if plot_err
%             fill_data = error_fill(x, y, y_err);
%             h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
%             set(h, 'facealpha', 0.2)
%         end
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',l_style)
    end
end
set(gca, 'ydir', 'reverse')
xlabel('Temperature (\circC)')
ylabel('weighted proximity to well')
    
% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
if equalLim
    fig = matchAxis(fig,true);
end
% ylim(num_temp_lim)

% save figure
save_figure(fig,[fig_dir  'Wieghted proximity by temp'],fig_type);



% FIGURE: OVERLAY OF WIEGHTED VS ABSOLUTE DISTANCE
plot_err = true;
LW = 1.5;

% FIGURES:
for i = 1:num.exp
    fig = getfig('',true,[551 680]); 
    hold on
    
    for type = 1:2 %increasing | decreasing 
        switch type
            case 1
                section_type = 'increasing';
                l_style = '-';
            case 2
                section_type = 'decreasing';
                l_style = '--';
        end
        % OG
        yyaxis left
        kolor = 'w';
        x = grouped(i).(section_type).temps;
        y = grouped(i).(section_type).avg;
        y_err = grouped(i).(section_type).err;
         loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',l_style)
        
        % WIEGHTED
        yyaxis right
        kolor = 'r';
       x = scaled_dist(i).temps;
        y = scaled_dist(i).(section_type)(:,1);
        y_err = scaled_dist(i).(section_type)(:,2);
         loc = isnan(y)|isnan(y_err);
        x(loc) = [];
        y(loc) = [];
        y_err(loc) = [];
        if plot_err
            fill_data = error_fill(x, y, y_err);
            h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
            set(h, 'facealpha', 0.2)
        end
        plot(x,y,'color',kolor,'linewidth',LW+1,'linestyle',l_style)
    end
    formatFig(fig,blkbgd);
    yyaxis left 
    set(gca, 'ycolor','w','ydir','reverse')
    ylabel('proximity to well (mm)')
    yyaxis right 
    set(gca, 'ycolor','r','ydir','reverse')
    ylabel('weighted proximity to well')
    
    title(expNames{i},'color','w')
    xlabel('temp (\circC)')
    % save figure
    save_figure(fig,[fig_dir  expNames{i} ' wieghted proximity by temp'],fig_type);

end

%% Video comparision of null distribution population
clearvars('-except',initial_vars{:})

if strcmpi(questdlg('Make video???'),'Yes')
    i = 1;
    trial = 1;
    vid = 1;
    frame = 1;
    
    SZ  = 25;
    row = 5; col = 10;
    sb(1).idx = [1:5,11:15, 21:25, 31:35, 41:45];
    sb(2).idx = [17:20, 27:30, 37:40];
    % sb(3).idx = 3;
    
    % Plotting data:
    binedges = null_dist.binedges*pix2mm;
    binedges_mm = null_dist.binedges;
    prob = null_dist.prob;
    arena_r = data(i).data(trial).data.r;
    arena_c = data(i).data(trial).data.centre;
    well_c = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));
    n_bins = length(binedges);
    
    % Video formatting
    vid_path = 'G:\My Drive\Jeanne Lab\DATA\Fundamentals\distribution size.avi';
    video_fps = 40; 
    
    % pull info for the first trial:
    dataDate = data(i).T.Date{trial};
    vid_name = data(i).T.ExperimentID{trial};
    vidDir = [baseFolder dataDate '/' vid_name '_'];
    videoPath = [vidDir num2str(vid) '.avi'];
    movieInfo = VideoReader(videoPath); %read in video
    
    % Set axis limits for the selected arena
    x = data(i).data(trial).data.centre(1);
    y = data(i).data(trial).data.centre(2);
    r = data(i).data(trial).data.r;
    xlimit = [x-(r+50),x+(r+50)];
    ylimit = [y-(r+50),y+50+r];
    
    % Open video object to which you will write new frames
    v = VideoWriter(vid_path, 'Uncompressed AVI');
    v.FrameRate = video_fps;
    open(v);
    
    % Plot image of video frame
    fig = getfig('',1,[1064 687]); set(fig,'color','k');
    for tt = 1:n_bins
        subplot(row,col,sb(1).idx); hold on
        set(gca,'color','k');
            currentImg = rgb2gray(read(movieInfo,frame));
            imshow(currentImg);
            xlim(xlimit); ylim(ylimit);  
            %plot full arena circle & well center
            viscircles(arena_c',arena_r,'Color','w'); %plot full arena circle
            scatter(well_c(1),well_c(2),SZ,"yellow",'filled');
            % plot current ring size
            viscircles(well_c',binedges(tt),'Color','r'); %plot full arena circle
        subplot(row,col,sb(2).idx); hold on
            plot(binedges_mm(1:tt),prob(1:tt),'color','r', 'linewidth', 2)
            xlim([binedges_mm(1),60]) %binedges_mm(end)])
            ylim([0 0.003])
            xlabel('distance (mm)')
            ylabel('probability')
            set(gca,'color','k','ycolor','w','xcolor','w','FontSize',18);
             % subplot(row,col,sb(2).idx); hold on
    
         % save the current figure as new frame in the video
        f = getframe(fig);
        writeVideo(v, f)
    end
    
    close(v)

end

%% Near inifinite place null distribution distance to CENTER OF ARENA
clearvars('-except',initial_vars{:})

% NULL DISTRIBUTION (this should be identical for all trials/experiments in the same arena)
% Generate the null distribution of distances to the food well:
pixel_buffer = 50; % edge of arena pixel buffer
i = 1;
trial = 1;
nbins = 10000;

% Set axis limits for the selected arena
c1 = data(i).data(trial).data.centre(1);
c2 = data(i).data(trial).data.centre(2);
r = data(i).data(trial).data.r;
xlimit = [c1-(r+pixel_buffer),c1+(r+pixel_buffer)];
ylimit = [c2-(r+pixel_buffer),c2+pixel_buffer+r];

% find the 'auto bin' lines
xedge = linspace(xlimit(1),xlimit(2),nbins+1);
yedge = linspace(ylimit(1),ylimit(2),nbins+1);
X = []; Y = [];
idx = 1;
for y_loc = 1:nbins
    for x_loc = 1:nbins
        X(idx) = xedge(x_loc);
        Y(idx) = yedge(y_loc);
        idx = idx + 1;
    end
end
% screen out the units outside the arena circle
temp_dist = sqrt((X-c1).^2 + (Y-c2).^2);
loc = temp_dist>r;
X(loc) = [];
Y(loc) = [];

% find food well  distance
null_distance_2center = sqrt((X-c1).^2 + (Y-c2).^2)./pix2mm;

binSpacing = 0.10; %mm spacing for bins
binedges = 0:binSpacing:55; 

% find probability of occupancy for each bin from the null distribution
N = discretize(null_distance_2center,binedges);
N_tot = length(null_distance_2center);
prob = zeros(length(binedges),1);
for i = 1:length(binedges)
    prob(i) = (sum(N==i))/N_tot;
end

null_dist_2center = struct;
null_dist_2center.distance = null_distance_2center;
null_dist_2center.binwidth = binSpacing;
null_dist_2center.binedges = binedges';
null_dist_2center.prob = prob;

% Save Null_Distribution
save([baseFolder 'Fundamentals\Distance_to_center_distribution.mat'],'null_dist_2center','-v7.3');


% Average 'eccentricity' given the distribution:
mean(null_distance_2center)


%% Find probability of occupancy for each distance
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);
LW = 3;

fig = getfig('',1); hold on
    yyaxis left
    histogram(null_dist_2center.distance,null_dist_2center.binedges,'FaceColor',Color('grey'),'FaceAlpha',1)
    yyaxis right
    plot(null_dist_2center.binedges,null_dist_2center.prob,'color', Color('teal'),'linewidth',LW,'linestyle','-')
    formatFig(fig,blkbgd);
    axis tight
    yyaxis left
    ylabel('count','FontSize',24)
    yyaxis right
    xlabel('eccentricity (mm)','FontSize',24)
    ylabel('occupancy probability','FontSize',24)
    yyaxis left
    set(gca,'Ycolor',foreColor)
    v_line(mean(null_distance_2center),'r','-')

    xlim([0,35])
% Save null distribution probability
save_figure(fig,[baseFolder 'Fundamentals\Eccentricity distribution'],fig_type);












