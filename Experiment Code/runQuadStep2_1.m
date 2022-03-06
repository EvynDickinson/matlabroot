


% Load in the selected data file: 
load(inputPath);
disp(['Starting figures for ' folder ' ' expName ' Arena ' arenaSel])

%% --------Find well count + distance from each well + temperature-------------------
% Find number of flies that are near each well for each frame
for vid = 1:nvids
    radii = 165; %110; %distance must be less than this number to count for a well ROI 
    dims = size(data(vid).occupancy_matrix);
    % loop for all wells 
    for well = 1:4
        b = wellcenters(:,well); 
        % reshape data points from whole video for optimized maths
        x_loc = reshape(data(vid).x_loc,numel(data(vid).y_loc),1); 
        y_loc = reshape(data(vid).y_loc,numel(data(vid).y_loc),1); 
        % within well range?
        euDist = sqrt((b(1)-x_loc).^2 + (b(2)-y_loc).^2); %euclidian distance from well center
        well_loc = reshape((euDist<=radii),[dims(2),dims(1)]);
        welldata(well).loc = well_loc;
        % how many within the circle?
        welldata(well).count = sum(well_loc,2);
        data(vid).well_counts(:,well) = welldata(well).count;
    end
    % temperature alignment
    logROI(1) = find(tempLog(:,1)==expData.tempLogStart(vid,3));
    logROI(2) = find(tempLog(:,1)==expData.tempLogEnd(vid,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    x = round(linspace(1, dims(2), length(tempCourse)));
    % upsample the temperature log:
    fullTempList = interp1(x,tempCourse,1:dims(2),'spline');   
    data(vid).tempLog = fullTempList;
end

clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

%% ------------------------------ Summary figure -------------------------------------
nrow = 5; ncol = 4;
subplotInd(1).idx = 5:7; % temperature
subplotInd(2).idx = [9:11,13:15,17:19]; % occupation
subplotInd(3).idx = 1:3; % fly count
subplotInd(4).idx = 4:4:20; % histogram

% group data across videos:
[plotX, plotY] = deal([]);
for vid = 1:nvids
    plotX = [plotX, data(vid).tempLog]; % temperature
    plotY = [plotY; data(vid).well_counts]; % flies per well
end
plotZ = occupancy.flycount;
plotY = plotY./nflies;
sSpan = 180;
LW = 2;
time = occupancy.time;

fig = getfig(''); 
 % tracking accuracy proxy (# flies)
 subplot(nrow,ncol,subplotInd(3).idx)
    y = smooth(plotZ,sSpan);
    roi = 2:length(y)-1;
    plot(time(roi), y(roi), 'linewidth', LW, 'color', Color('grey'))
    hline(nflies, 'w--')
    ylabel('fly count')
 
 % temperature over time
 subplot(nrow,ncol,subplotInd(1).idx)
    y = smooth(plotX,sSpan);
    roi = 2:length(y)-1;
    plot(time(roi), y(roi), 'linewidth', LW, 'color', 'w')
    ylabel('temp (\circC)')
    ylim([5,26])
    
 % occupation probability
 subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        y = smooth(plotY(:,well),sSpan);
        roi = 2:length(y)-1;
        plot(time(roi), y(roi), 'linewidth', LW, 'color', kolor);
    end
    xlabel('time (min)'); ylabel('occupation probability')
    
 % fly count histogram
 subplot(nrow,ncol,subplotInd(4).idx)
    yyaxis left
    set(gca,'YColor', 'k')
    yyaxis right
    h = histogram(flyCount); vline(nflies, 'w--'); 
    set(h, 'facecolor', Color('grey'))
    xlabel('Number tracked flies'); ylabel('Frame count')
%     labelHandles = findall(gca, 'type', 'text', 'handlevisibility', 'off');
%     set(labelHandles,'FontSize', 18);
%     set(gca,'fontsize',15,'FontWeight','normal');
    
formatFig(fig, true, [nrow, ncol], subplotInd);
subplot(nrow,ncol,subplotInd(2).idx)
l = legend(strrep(wellLabels,'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','position', [0.6463 0.5633 0.0963 0.1126]);% [0.8780 0.8119 0.0963 0.1126])%
subplot(nrow,ncol,subplotInd(1).idx)
set(gca, 'XColor', 'k')
subplot(nrow,ncol,subplotInd(3).idx)
set(gca, 'XColor', 'k')
titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
title(titleName,'color', 'w')

% Save image
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
    end
end  
save_figure(fig, [analysisDir expName arenaSel ' summary figure'], '-png', autoSave);


clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

%% -------Simple visualization: relationship between temp & well occupation-----------
if essentialfigs == false
    nrow = 4;
    ncol = 1;
    subplotInd(1).idx = 1;
    subplotInd(2).idx = 2:4;
    % group data across videos:
    [plotX, plotY] = deal([]);
    for vid = 1:nvids
        plotX = [plotX, data(vid).tempLog];
        plotY = [plotY; data(vid).well_counts];
    end
    plotY = plotY./nflies;
    sSpan = 180;
    LW = 1;
    time = occupancy.time;

    fig = getfig(''); 
        subplot(nrow,ncol,subplotInd(1).idx)
        y = smooth(plotX,sSpan);
        plot(time(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'w')
        ylabel('temp (\circC)')
        ylim([5,27])

        subplot(nrow,ncol,subplotInd(2).idx)
        hold on
        % error fills
        for well = 1:4
            kolor = pullFoodColor(wellLabels{well}); % plotting color for food
            y_avg(:,well) = smooth(plotY(:,well),sSpan);
            y_err = movstd(plotY(:,well),sSpan);
            fill_data = error_fill(time, y_avg(:,well), y_err);
            h(well) = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h(well), 'facealpha', 0.2)
        end
        % average line
        for well = 1:4
            kolor = pullFoodColor(wellLabels{well});
            plot(time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
        end
        xlabel('time (min)'); ylabel('occupation probability')

    formatFig(fig, true, [nrow, ncol], subplotInd);
    l = legend([{'';'';'';''};strrep(wellLabels,'_','-')]);
    set(l, 'color', 'k', 'textcolor', 'w','position', [0.7947 0.6462 0.0963 0.1126])
    for well = 1:4
        set(get(get(h(well),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    subplot(nrow,ncol,subplotInd(1).idx)
    set(gca, 'XColor', 'k')
    titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
    title(titleName,'color', 'w')


    % Save image
    save_figure(fig, [analysisDir expName arenaSel ' well occupation timcourse'], '-png', autoSave);

    initial_vars = [initial_vars; 'time'];
    clearvars('-except',initial_vars{:})
    fprintf('\nNext\n')
end

%% ----------------------Organize info from experiment--------------------------------
% group data across videos:
[plotX, plotY] = deal([]);
for vid = 1:nvids
    plotX = [plotX, data(vid).tempLog];
    plotY = [plotY; data(vid).well_counts];
end
occupancy.temp = plotX;
occupancy.occ = plotY./nflies;
occupancy.count = plotY;
% add total well occupancy to occupancy structure
occupancy.allwellOcc = sum(occupancy.occ,2);


% Measure of eccentricity
%how far are the flies from the center of the arena, on average?
% auto generate the arena center from the middlepoint of the 4 wells?
x1 = wellcenters(1,1:2:4);
y1 = wellcenters(2,1:2:4);
x2 = wellcenters(1,2:2:4);
y2 = wellcenters(2,2:2:4);
[xi,yi] = polyxpoly(x1,y1,x2,y2);
wellcenters(:,5) = [xi;yi];
% 
% %visualize the arena center
% fig = getfig; set(fig, 'color', 'k')
% imshow(demoImg); axis tight square
% hold on
% scatter(xi,yi, 45, 'y', 'filled')

% find distance from center for each fly:
EDist = [];
for vid = 1:nvids
    N = [];
    for frame = 1:length(data(vid).tempLog)
        Y = data(vid).y_loc(frame,:);
        X = data(vid).x_loc(frame,:);
        centre = wellcenters(:,5);
        D = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5); %distance from center of arena
        D(isnan(D)) = [];
        N(frame,:) = [mean(D), std(D)];
    end
    data(vid).eccentricity = N;
    EDist = [EDist;N];
end
occupancy.eccentricity = EDist;

% Plot the results of the eccentricity ** could be cool to do this overlaid
% onto the arena image like a polar plot but with temp being color
% coordinated ** 
% 
% % % Try polar plot vectorization:
% function [r,theta]=cart2polar(x,y)
%     r=sqrt(x.^2+y.^2);
%     theta=atan(y./x);
% end

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% -----------------------Calculate degree of clustering------------------------------
% pull info for each video
LDist = []; LD_err = [];
for vid = 1:nvids
    flyDistance = []; D_err = [];
    for frame = 1:length(data(vid).tempLog)
        test = [data(vid).x_loc(frame,:)',data(vid).y_loc(frame,:)'];
        D = pdist(test);
        D(isnan(D)) = [];
        flyDistance(frame) = mean(D);
        D_err(frame) = std(D);
    end
    data(vid).flyDistance = flyDistance;
    LDist = [LDist,flyDistance]; 
    LD_err = [LD_err,D_err];
end

%add inter-fly-distance to the occupancy time course structure
occupancy.IFD = LDist;
occupancy.IFD_err = LD_err;

% demo the clustering proxy
if nvids>8
    divisor = round(nvids/6);
    vidList = 1:divisor:nvids;
else
    vidList = 1:nvids;
end
nrow = 2; ncol = length(vidList); ii = 0; 
[~,minidx] = deal([]);

% VISUALIZE a demo of the clustering accuracy
if essentialfigs == false
    fig = getfig(''); set(fig, 'pos',[120 331 1244 368], 'color', 'k');
    for vid = vidList
        ii = ii+1;
        movieInfo = VideoReader([baseFolder vidFolder '\' expName '_' num2str(vid) '.avi']); %read in video
        headData = squeeze(data(vid).tracks(:,1,:,:));

        % Most clustered:
        [M,minidx(ii)] = min(data(vid).flyDistance);
        img = read(movieInfo,minidx(ii));
        subplot(nrow, ncol, ii)
        imshow(img); hold on
        axis tight square
        set(gca, 'visible', 'off')
        x = data(vid).x_loc(minidx(ii),:); x(isnan(x)) = [];
        y = data(vid).y_loc(minidx(ii),:); y(isnan(y)) = [];
    %     x = squeeze(headData(minidx(ii), 1, :)); x(isnan(x)) = [];
    %     y = squeeze(headData(minidx(ii), 2, :)); y(isnan(y)) = [];
        scatter(x,y, 10, 'y', 'filled')
        title(num2str(M))
        % overlay 'size bar' for min dist:
        plot([100,100+M], [20,20], 'linewidth', 0.5, 'color', 'r')

        % Least clustered:
        [M,maxidx(ii)] = max(data(vid).flyDistance);
        img = read(movieInfo,maxidx(ii));
        subplot(nrow, ncol, ii+length(vidList))
        imshow(img); hold on
        axis tight square
        set(gca, 'visible', 'off')
        x = data(vid).x_loc(maxidx(ii),:); x(isnan(x)) = [];
        y = data(vid).y_loc(maxidx(ii),:); y(isnan(y)) = [];
    %     x = squeeze(headData(maxidx(ii), 1, :)); x(isnan(x)) = [];
    %     y = squeeze(headData(maxidx(ii), 2, :)); y(isnan(y)) = [];
        scatter(x,y, 10, 'y', 'filled')
        title(num2str(M))
        % overlay 'size bar' for min dist:
        plot([100,100+M], [20,20], 'linewidth', 0.5, 'color', 'r')
    end
    labelHandles = findall(gcf, 'type', 'text', 'handlevisibility', 'off');
    set(labelHandles,'FontSize', 15, 'color', 'w')

    save_figure(fig, [analysisDir expName arenaSel ' linear clustering demo'], '-pdf', autoSave);
end

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% ------------------------------Movement analysis------------------------------------
%bin frame and check bin occupation changes across frames as proxy for
%movement --> won't give an accurate 'speed' etc. but it might give
%movement.
tic
nbins = 100;
spaceChange = [];

for vid = 1:nvids   
    N = [];
    for frame = 1:length(data(vid).tempLog)
        X = data(vid).x_loc(frame,:); X(isnan(X)) = [];
        Y = data(vid).y_loc(frame,:); Y(isnan(Y)) = [];
        N(:,:,frame) = histcounts2(X,Y,nbins);
    end
    if vid>1 % add last frame of previous vid to prevent frame loss 
        N = cat(3, lastFrame, N);
    end
    lastFrame =  N(:,:,end);
    binDiff = diff(N,1,3);
    binDiff(binDiff<0) = 0;
    BinChange = squeeze(sum(sum(binDiff,1),2));
    % save pertinent information:
    data(vid).BinChange = BinChange;
    spaceChange = [spaceChange; BinChange];
end
occupancy.movement = spaceChange;
toc

% preview the movement average
if essentialfigs==false
    fig = getfig;  hold on
    plot(occupancy.time(2:end), occupancy.movement,'color', Color('teal'))
    plot(occupancy.time(2:end), smooth(occupancy.movement,180),...
        'linewidth', 2, 'color', 'w')
    xlabel('time (min)'), ylabel('movement (a.u.)')
    formatFig(fig, true);
    save_figure(fig, [analysisDir expName arenaSel ' movement over time'], '-pdf', autoSave);
end

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% -----------------------Time course feature comparison------------------------------
% plot parameters:
nrow = 8; ncol = 1;
sSpan = 180; %1 minute filter length
sbpts(1).idx = 1;
sbpts(2).idx = 2:4;
sbpts(3).idx = 5:6;
sbpts(4).idx = 7:8;
LW = 1.5;

% FIGURE:
fig = getfig; set(fig, 'pos', [157 86 1232 878])
% TEMPERATURE
subplot(nrow,ncol,sbpts(1).idx)
plot(occupancy.time,smooth(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
ylabel('(\circ)')
title('temperature')
ax = gca;
x_lim = ax.XLim;

% OCCUPANCY
subplot(nrow,ncol,sbpts(2).idx)
hold on
y = [];
for well = 1:4
   y = [y, smooth(occupancy.occ(:,well),sSpan)];
end
h = area(occupancy.time,y);
for well = 1:4
    h(well).FaceColor = pullFoodColor(wellLabels{well});
end
xlim(x_lim)
set(gca, 'tickdir', 'out')
l = legend(strrep(wellLabels,'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k',...
'position', [0.7972 0.7194 0.1039 0.0792]);% [0.8780 0.8119 0.0963 0.1126])%
ylabel('occupancy probability')
title('individual well occupancy')

% 
% hold on
% for well = 1:4
%     kolor = pullFoodColor(wellLabels{well});
%     plot(time,smooth(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
% end
% plot(time,smooth(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
% ylabel('occupancy probability')
% title('individual well occupancy')
% legend(wellLabels)
% l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
% set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
%     'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%

% CLUSTERING
subplot(nrow,ncol,sbpts(3).idx); hold on
y_avg = smooth(occupancy.IFD,sSpan);
y_err = smooth(occupancy.IFD_err,sSpan);
x = occupancy.time;
fill_data = error_fill(x, y_avg, y_err);
h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(occupancy.time,smooth(occupancy.IFD,sSpan),'linewidth', LW, 'color', 'w')
ylabel('pixels')
title('inter-fly-distance')

% MOVEMENT
subplot(nrow,ncol,sbpts(4).idx); hold on
y_avg = smooth(occupancy.movement,sSpan);
y_err = movstd(occupancy.movement,sSpan);
x = occupancy.time(2:end);
fill_data = error_fill(x, y_avg, y_err);
h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(x,y_avg, 'linewidth', LW, 'color', 'w')  
ylabel('(a.u.)')
title('movement')
xlabel('time (min)')

formatFig(fig,true,[nrow, ncol], sbpts);
for ii = 1:3
    subplot(nrow,ncol,sbpts(ii).idx)
    set(gca, 'XColor', 'k')
end

% save and export figure
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
    end
end  
save_figure(fig, [analysisDir expName arenaSel ' timecourse features'], '-png',autoSave);

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% ---------------------Position & Movement Summary Figure----------------------------
% **TODO remove total well occupancy
% figure parameters:
nrow = 10; ncol = 1;
sSpan = 180; %1 minute filter length
sbpts(1).idx = 1;    % temp
sbpts(2).idx = 2:5;  % occupancy
sbpts(3).idx = 6:8;  % clustering & eccentricity
sbpts(4).idx = 9:10; % movement
LW = 1.5;

% FIGURE:
fig = getfig; set(fig, 'pos', [157 86 1232 878])
subplot(nrow,ncol,sbpts(1).idx)
plot(occupancy.time,smooth(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
ylabel('Temp(\circ)')
% title('temperature')

subplot(nrow,ncol,sbpts(2).idx)
hold on
for well = 1:4
    kolor = pullFoodColor(wellLabels{well});
    plot(occupancy.time,smooth(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
end
% plot(occupancy.time,smooth(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
ylabel('occupancy probability')
% title('individual well occupancy')
legend(wellLabels)
% l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
l = legend(strrep(wellLabels,'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
    'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%

% CLUSTERING
subplot(nrow,ncol,sbpts(3).idx); hold on
y_avg = smooth(occupancy.IFD,sSpan);
y_err = smooth(occupancy.IFD_err,sSpan);
x = occupancy.time;
fill_data = error_fill(x, y_avg, y_err);
h(1) = fill(fill_data.X, fill_data.Y, get_color('teal'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(occupancy.time,smooth(occupancy.IFD,sSpan),'linewidth', LW, 'color', Color('teal'))
% ECCENTRICITY
y_avg = smooth(occupancy.eccentricity(:,1),sSpan);
y_err = smooth(occupancy.eccentricity(:,2),sSpan);
x = occupancy.time;
fill_data = error_fill(x, y_avg, y_err);
h(2) = fill(fill_data.X, fill_data.Y, get_color('orange'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(x,y_avg, 'linewidth', LW, 'color', Color('orange'))  
ylabel('pixels')
l2 = legend({'','inter fly distance', '', 'eccentricity'});
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(l2, 'color', 'k', 'textcolor', 'w','FontSize', 8,'edgecolor', 'k',...
    'position', [0.1350 0.4638 0.1039 0.0370]);% 

% MOVEMENT
subplot(nrow,ncol,sbpts(4).idx); hold on
y_avg = smooth(occupancy.movement,sSpan);
y_err = movstd(occupancy.movement,sSpan);
x = occupancy.time(2:end);
fill_data = error_fill(x, y_avg, y_err);
h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
set(h, 'facealpha', 0.2)
plot(x,y_avg, 'linewidth', LW, 'color', 'w')  
ylabel('movement (a.u.)')
% title('movement')
xlabel('time (min)')

formatFig(fig,true,[nrow, ncol], sbpts);
for ii = 1:3
    subplot(nrow,ncol,sbpts(ii).idx)
    set(gca, 'XColor', 'k')
end

% save and export figure
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
    end
end
save_figure(fig, [analysisDir expName arenaSel ' timecourse features'], '-png',autoSave);

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% ----------------Average distance between each fly and each food source-------------
%(takes a full 2 mins to run)

% find distance from center for each fly:
tic
FDist = [];
idx = 0;
for vid = 1:nvids
    for frame = 1:length(data(vid).tempLog)
        idx = idx+1;
        for well = 1:4
            test = [wellcenters(:,well)'; data(vid).x_loc(frame,:)',data(vid).y_loc(frame,:)'];
            D = squareform(pdist(test));
            D = D(2:end,1);
            D(isnan(D)) = [];
            FDist(well).N(idx,:) = [median(D), std(D)];
        end
    end
end
toc

occupancy.dist2wells = FDist;

% occupancy.eccentricity = EDist;
nrow = 4;
ncol = 1;
subplotInd(1).idx = 1;
subplotInd(2).idx = 2:4;
% group data across videos:
plotX = occupancy.time;
sSpan = 180;
LW = 1;

fig = getfig(''); 
    subplot(nrow,ncol,subplotInd(1).idx)
    y = smooth(occupancy.temp,sSpan);
    plot(plotX(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'w')
    ylabel('temp (\circC)')
    ylim([5,27])
    
    subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    % error fills
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well}); % plotting color for food
        y_avg(:,well) = smooth(FDist(well).N(:,1),sSpan);
        y_err(:,well) = smooth(FDist(well).N(:,2),sSpan);
        fill_data = error_fill(plotX, y_avg(:,well), y_err(:,well));
        h(well) = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
        set(h(well), 'facealpha', 0.2)
    end
    % average line
    for well = 1:4
        kolor = pullFoodColor(wellLabels{well});
        plot(occupancy.time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
    end
    xlabel('time (min)'); ylabel('avg distance between fly and food (a.u.)')
    
formatFig(fig, true, [nrow, ncol], subplotInd);
l = legend([{'';'';'';''};strrep(wellLabels,'_','-')]);
set(l, 'color', 'k', 'textcolor', 'w','position', [0.7947 0.6462 0.0963 0.1126])
for well = 1:4
    set(get(get(h(well),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
subplot(nrow,ncol,subplotInd(1).idx)
set(gca, 'XColor', 'k')
titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
title(titleName,'color', 'w')
 
% save and export figure
if autoSave==true
    export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
else
    if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
    end
end
save_figure(fig, [analysisDir expName arenaSel ' avg distance to well'], '-png', autoSave);

clearvars('-except',initial_vars{:})


%% ----------------------- Save experiment data thus far -----------------------------
if autoSave == true
    clearvars('-except',initial_vars{:})
    initial_vars = unique(initial_vars);
    save([analysisDir expName arenaSel ' timecourse data'])
    copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')
    fprintf('Experiment data saved\n')
    results = 'Saved data';
else
    switch questdlg('Save experiment analysis?')
        case 'Yes'
            clearvars('-except',initial_vars{:})
            initial_vars = unique(initial_vars);
            save([analysisDir expName arenaSel ' timecourse data'])
            copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')
            fprintf('Experiment data saved\n')
            results = 'Saved data';
        case 'No'
            results = 'No data saved';
            return
        case 'Cancel'
            results = 'No data saved';
            return
    end
end
