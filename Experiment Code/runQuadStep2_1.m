
function results = runQuadStep2_1(inputPath,autoSave,essentialfigs)
% results = runQuadStep2(inputPath,autoSave,essentialfigs)
%
%
% ES Dickinson

% Load in the selected data file: 
warning off
load(inputPath);
disp(['Starting figures for ' folder ' ' expName])

nArenas = 4;
initial_vars{end+1} = 'nArenas';
initial_vars{end+1} = 'figDir';
figDir = [baseFolder folder '/'];
for arena = 1:nArenas
    arenaData(arena).figDir = [figDir arenaData(arena).name '/'];
end
tPoints = getTempTurnPoints(expData.parameters.protocol);
threshHigh = tPoints.threshHigh;
threshLow = tPoints.threshLow;
binSpace = 1; %temp degree bin widths
initial_vars = [initial_vars(:)', 'arena', 'threshHigh', 'threshLow', 'binSpace'];

%% ------------------------------ summary figure -------------------------------------
nrow = 5; ncol = 4;
subplotInd(1).idx = 5:7; % temperature
subplotInd(2).idx = [9:11,13:15,17:19]; % occupation
subplotInd(3).idx = 1:3; % fly count
subplotInd(4).idx = 4:4:20; % histogram

for arena = 1:nArenas
    
    % group data across videos:
    plotZ = T.flyCount(:,arena);
    plotY = arenaData(arena).occ_P;
    sSpan = 180;
    LW = 2;
    time = T.time;
    wellLabels = arenaData(arena).wellLabels;
    arenaSel = arenaData(arena).name;

    fig = getfig(''); 
     % tracking accuracy proxy (# flies)
     subplot(nrow,ncol,subplotInd(3).idx)
        y = smooth(plotZ,sSpan);
        roi = 2:length(y)-1;
        plot(time(roi), y(roi), 'linewidth', LW, 'color', Color('grey'))
        hline(arenaData(arena).nflies, 'w--')
        ylabel('fly count')

     % temperature over time
     subplot(nrow,ncol,subplotInd(1).idx)
        y = smooth(T.temperature,sSpan);
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
        h = histogram(arenaData(arena).flyCount); vline(arenaData(arena).nflies, 'w--'); 
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
    titleName = strrep([folder ' ' expName arenaSel], '_',' ');
    title(titleName,'color', 'w')

    % Save image
    expPDF = [arenaData(arena).figDir folder ' ' expName ' ' arenaSel ' summary.pdf'];
    if autoSave==true
        export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
    else
        if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
            export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
        end
    end  
    save_figure(fig, [arenaData(arena).figDir expName ' ' arenaSel ' summary figure'], '-png', autoSave);
end

clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

%% -------Simple visualization: relationship between temp & well occupation-----------
if essentialfigs == false
  for arena = 1:nArenas
    nrow = 4;
    ncol = 1;
    subplotInd(1).idx = 1;
    subplotInd(2).idx = 2:4;
    sSpan = 180;
    LW = 1;
    wellLabels = arenaData(arena).wellLabels;
    time = T.time;
    occ = arenaData(arena).occ_P;
    arenaSel = arenaData(arena).name;
    
    % MAKE THE FIGURE:
    fig = getfig(''); 
        % TEMPERATURE data
        subplot(nrow,ncol,subplotInd(1).idx)
        y = smooth(T.temperature,sSpan);
        plot(T.time(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'w')
        ylabel('temp (\circC)')
        ylim([5,27])
        % OCCUPANCY
        subplot(nrow,ncol,subplotInd(2).idx)
        hold on
        % error fills
        for well = 1:4
            kolor = pullFoodColor(wellLabels{well}); % plotting color for food
            y_avg(:,well) = smooth(occ(:,well),sSpan);
            y_err = movstd(occ,sSpan);
            fill_data = error_fill(time, y_avg(:,well), y_err(:,well));
            h(well) = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
            set(h(well), 'facealpha', 0.3)
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
    titleName = strrep([folder ' ' expName ' ' arenaSel], '_',' ');
    title(titleName,'color', 'w')

    % Save image
    save_figure(fig, [arenaData(arena).figDir expName ' ' arenaSel ' well occupation timcourse'], '-png', autoSave);
  end
  clearvars('-except',initial_vars{:})
  fprintf('\nNext\n')
end

%% Calculate several parameters: eccentricity, movement, clustering

% organize data for forward compatability
for arena = 1:nArenas
    tic
    disp(['Starting ' arenaIdx{arena}'])
    occupancy.temp = T.temperature;
    occupancy.occ = arenaData(arena).occ_P;
    occupancy.count = arenaData(arena).occ_N;
    occupancy.allwellOcc = sum(occupancy.occ,2);

    % ECCENTRICITY
    N = [];
    x1 = arenaData(arena).wellcenters(1,1:2:4);
    y1 = arenaData(arena).wellcenters(2,1:2:4);
    x2 = arenaData(arena).wellcenters(1,2:2:4);
    y2 = arenaData(arena).wellcenters(2,2:2:4);
    [xi,yi] = polyxpoly(x1,y1,x2,y2);
    arenaData(arena).wellcenters(:,5) = [xi;yi];
    % %visualize the arena center
    % fig = getfig; set(fig, 'color', 'k')
    % imshow(demoImg); axis tight square
    % hold on
    % scatter(xi,yi, 45, 'y', 'filled')
    % find distance from center for each fly: 
    D = sqrt(((arenaData(arena).x_loc-xi).^2 + (arenaData(arena).y_loc-yi).^2))./pix2mm; %distance from center of arena
    N(:,1) = mean(D,2,'omitnan');
    N(:,2) = std(D,0,2,'omitnan');
    occupancy.eccentricity = N;
    
    % CLUSTERING
    [LDist,LD_err] = deal([]); 
    for frame = 1:T.frame(end)
        inputData = [arenaData(arena).x_loc(frame,:)',arenaData(arena).y_loc(frame,:)'];
        inputData(isnan(inputData(:,1)),:) = [];
        D = pdist(inputData);
        LDist(frame) = mean(D);
        LD_err(frame) = std(D);
    end
    occupancy.IFD = LDist;
    occupancy.IFD_err = LD_err;
     
    % MOVEMENT: 
    nbins = 100;
    N = [];
    for frame = 1:T.frame(end)
        X = arenaData(arena).x_loc(frame,:); X(isnan(X)) = [];
        Y = arenaData(arena).y_loc(frame,:); Y(isnan(Y)) = [];
        N(:,:,frame) = histcounts2(X,Y,nbins);
    end
    binDiff = diff(N,1,3);
    binDiff(binDiff<0) = 0;
    occupancy.movement = squeeze(sum(sum(binDiff,1),2));
    
    % save data into occupancy structure:
    occupancy.dist2wells = arenaData(arena).dist2well;
    arenaData(arena).occupancy = occupancy;
    toc
end

for arena = 1:nArenas
    arenaData(arena).genotype = expData.parameters.(['Arena' arenaIdx{arena}]).genotype;
end

save([analysisDir 'half processed data.mat'])

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Temp rate heat map: 

% Get the temp rate, temp, and distance from food
for arena = 1:nArenas
    [TR,T_rates] = deal([]); % temp rate data will go here then get saved into Arena Data

    % temperature
    temp = T.temperature;
    time = T.time(1:end-1);
    plotData(:,1) = temp(1:end-1); 

    % rate of temperature change
    dT = diff(smooth(temp,180)); %180 = 1 minute smoothing kernal
    dt = diff(T.time(:)); 
    plotData(:,2) = dT./dt; 
    % figure; plot(T.time(1:end-1),plotData(:,2)) %to look at the calculated temp rate

    % distance from food
    for well = 1:4
        [kolor,num] = pullFoodColor(arenaData(arena).wellLabels{well});
        if ~(num==3) % aka not an empty well
            y = arenaData(arena).dist2well(:,well);
            foodwell = well;
            arenaData(arena).foodwell = foodwell;
            plotData(:,3) = y(1:end-1);
            break
        end
    end

    % find the mean temp rate during each ramp period:
    tPoints = getTempTurnPoints(expData.parameters.protocol); %accomodates multiple temp protocols within the data group
    keepLoc = false(1,size(plotData,1));
    for ii = 1:length(tPoints.transitions)-1
        roi = tPoints.transitions(ii):tPoints.transitions(ii+1);
        distD = median(plotData(roi,2));
        plotData(roi,4) = round(distD*ones(range(roi)+1,1),2);
        keepLoc(roi) = true;
    end
    plotData(~keepLoc,:) = nan; % exclude data outside prescribed ROI regions
%     figure; 
%     plot(T.time(1:end-1),plotData(:,2))
%     hold on 
%     time = T.time(1:end-1);
%     scatter(time(keepLoc),plotData(keepLoc,2))

    % Temp-rate identification and sorting: 
    buffSize = 0.05;
    if any(ismember(strsplit(expData.parameters.protocol,'_'),'sweeps')) %account for high freq jitter
        buffSize = 0.1;
    end
    for ii = 1:tPoints.nRates
        edges(ii,:) = [tPoints.rates(ii)-buffSize, tPoints.rates(ii)+buffSize];
    end  
    rateData = plotData(:,4);
    nRates = tPoints.nRates;
    rateIdx = discretize(rateData,nRates);
    for ii = 1:nRates
        TD = rateData(rateIdx==ii);
        T_rates(ii) = round(mean(TD,'omitnan'),2);
        % do these match the assumed rates?
        idx = find(T_rates(ii) > edges(:,1) & T_rates(ii) < edges(:,2));
        if isempty(idx)
            warndlg('Temp rate not within expected range')
            return
        end
        % set the rate value to the uniform cross-group value:
        rateData(rateIdx==ii) = tPoints.rates(ii);
        T_rates(ii) = tPoints.rates(ii);
%        figure; hold on
%        plot(rateData)
%        plot(plotData(:,4))
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
    [HM,FC,FD,dist_mat,FC_mat] = deal([]);
    for col = 1:nTemps
        for row = 1:nRates
        % Pull the data that matches this category
        loc = (rateIdx==row) & (tempIdx==col);

        % fly count data: (for tracking assessment)
        fc = flyCount(loc,arena);
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
        end
    end
    TR.heatmap = HM;
    TR.dist_mat = dist_mat;
    TR.FC_mat = FC_mat;
    TR.nRates = length(TR.rates);
    arenaData(arena).TR = TR;

    % FIGURE: show the temp rate break down 
    if arena == 1
    fig = figure; yyaxis left; plot(time,tempIdx,'Color', Color('DodgerBlue')); ylabel('Temp bins');
        yyaxis right; plot(time,rateIdx,'color', Color('orange')); ylabel('Heating vs. cooling'); 
        formatFig(fig,true); yyaxis right; set(gca, 'YColor', Color('orange')); 
        yyaxis left; set(gca, 'YColor', Color('DodgerBlue'));xlabel('Time (min)')
        save_figure(fig, [analysisDir 'Temp rate assignment'],'-png',true);
    end

    % FIGURE: temperature rate aligment 
    title_str = [folder ' ' arenaData(arena).name ' ' arenaData(arena).genotype];
    fig_file = [arenaData(arena).figDir expName ' temp temp_rate dist'];
    if ~isfile([fig_file '.png']) %skip already generated figures
        fig = figure; set(fig, 'pos', [283 226 1019 622])
        % temp
        subplot(3,1,1)
            plot(time, plotData(:,1),'color','w')
            v_line(time(tPoints.transitions),'y',':')
            ylabel('Temp (\circC)')
            title(title_str)
        % temp rate
        subplot(3,1,2); hold on
            plot(time, plotData(:,2),'color','w')
            plot(time, plotData(:,4),'color','r','linewidth', 2)
            v_line(time(tPoints.transitions),'y',':')
            ylabel('dT/dt (\circC/min)')
        % distance from food
        subplot(3,1,3)
            plot(time, plotData(:,3),'color',kolor)
            v_line(time(tPoints.transitions),'y',':')
            xlabel('time (min)')
            ylabel('distance (mm)')
            yyaxis right
            yy = arenaData(arena).occ_P(1:end-1,foodwell);
            yy(~keepLoc) = nan;
            plot(time, smoothdata(yy,'movmean',180),'color', 'w')
            ylabel('occupancy (prob)')
        fig = formatFig(fig, true, [3,1]);
            subplot(3,1,3)
            yyaxis left 
            set(gca,'YColor', kolor)
            save_figure(fig, fig_file, '-png', true);
    end
    
    % PLOT DISTANCE AS FUNCTION OF TEMP AND RATE OF TEMP
    fig_file = [arenaData(arena).figDir expName ' temp temp_rate dist heatmap'];
    if ~isfile([fig_file '.png']) %skip already generated figures
    %   colormap default
        fig = figure; set(fig, 'pos', [560 127 983 417]);
            hold on
            imAlpha=ones(size(HM));
            imAlpha(isnan(HM))=0;
            imagesc(HM,'AlphaData',imAlpha);
            set(gca,'color',0*[1 1 1]);
            axis tight
            title(title_str)
            % Axes formatting
            ax = gca;
            fig = formatFig(fig, true);
            XtickNum = ax.XTick;
            try ax.XTickLabel = t_roi(XtickNum);
            catch
                ax.XTick = 1:length(t_roi);
                ax.XTickLabel = t_roi(ax.XTick);
            end
%             YtickNum = ax.YTick;
            set(gca, 'ytick', 1:nRates,'YTickLabel',T_rates)
            ylabel('\DeltaT/dt (\circC/min)')
            xlabel('Temp (\circC)')
            % Colorbar formatting
            cbh = colorbar(); 
            cbh.Label.String = 'Distance from food (mm)';
            cbh.Color = Color('white');
            % flip colormap around to make yellow closer to food
            cmap = colormap;
            set(gca, 'colormap', flip(cmap))
        save_figure(fig, fig_file, '-png', true);
    end
    
    % FIGURE: Line plots of each rate comparison 
    fig_file = [arenaData(arena).figDir 'food distance hysteresis'];
    if ~isfile([fig_file '.png'])
        LS = {'--','-.','-'}; %cooling|stationary|heating
        nrows = 1;
        ncols = 3;
        sb(1).idx = 1; %distance vs temp
        sb(2).idx = 2; %fly count vs temp
        sb(3).idx = 3; %fly count across temp rates
        fig = figure; set(fig, 'pos', [17 49 1354 603])
        % DISTANCE
        subplot(nrows,ncols,sb(1).idx)
        hold on
        for rr = 1:TR.nRates
            curr_rate = TR.rates(rr);
            kolor = pullFoodColor(curr_rate);
            if curr_rate>0
                lstyle = LS{3};
            elseif curr_rate<0
                lstyle = LS{1};
            elseif curr_rate==0
                continue
                lstyle = LS{2};
            end
            x = TR.temps;
            y = dist_mat.avg(rr,:);
            % err = dist_mat.err(rr,:);
            plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
        end
    
        ylabel('Distance to food (mm)')
        xlabel('Temperature (\circC)')
        
        % FLY COUNT by rate
        subplot(nrows,ncols,sb(2).idx)
        hold on
        for rr = 1:nRates
            curr_rate = TR.rates(rr);
            kolor = pullFoodColor(curr_rate);
            if curr_rate>0
                lstyle = LS{3};
            elseif curr_rate<0
                lstyle = LS{1};
            elseif curr_rate==0
                continue
                lstyle = LS{2};
            end
            x = TR.temps;
            y = FC_mat.avg(rr,:);
            % err = dist_mat.err(rr,:);
            plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
        end
    %     h_line(T.NumFlies(trial),'r',':',1)
        ylabel('Fly Count (#)')
        xlabel('Temperature (\circC)')
        title(title_str)
        fig = formatFig(fig, true,[nrows,ncols],sb);
        
        % FLY COUNT across all temps
        subplot(nrows,ncols,sb(3).idx)
        idx = discretize(T.temperature(:),TR.temps);
        boxplot(flyCount(:,arena), idx,'Plotstyle', 'traditional','colors', 'w')
        h_line(arenaData(arena).nflies, 'Teal', '-',3)
        xlabel('Temperature (\circC)')
        ylabel('Tracking Fly Count')
        
        set(gca, 'color', 'k','box', 'off','XColor', 'w', 'YColor',...
                 'w','XTickLabels', TR.temps,...
                 'FontSize', 14)
        save_figure(fig, fig_file, '-png',true);
    end

clearvars('-except',initial_vars{:}) 
end

%% Build a data structure for each arena from occupancy, T, and TR data

for arena = 1:nArenas
    AD = T(:,1:6); % shared variables across all trials
    AD = addvars(AD,T.flyCount(:,arena),'NewVariableNames','flyCount');
    AD = addvars(AD,arenaData(arena).occ_P,'NewVariableNames','occ_P');
    AD = addvars(AD,arenaData(arena).occ_N,'NewVariableNames','occ_N');
    AD = addvars(AD,arenaData(arena).occupancy.eccentricity,'NewVariableNames','eccentricity');
    AD = addvars(AD,[arenaData(arena).occupancy.movement;nan],'NewVariableNames','movement');
    AD = addvars(AD,[arenaData(arena).occupancy.IFD',arenaData(arena).occupancy.IFD_err'],...
        'NewVariableNames','IFD');
    AD = addvars(AD,arenaData(arena).occupancy.dist2wells,'NewVariableNames','dist2wells');    
    AD = addvars(AD,arenaData(arena).occupancy.dist2wells(:,arenaData(arena).foodwell),...
        'NewVariableNames','dist2food');
    AD = addvars(AD,[arenaData(arena).TR.rateIdx; nan],'NewVariableNames','temprateIdx');
    temp = arenaData(arena).TR.data;
    AD = addvars(AD,[temp(:,4); nan],'NewVariableNames','temprate');
    % add degree-binned temp column:
    temp = AD.temperature;
    AD = addvars(AD,round(temp),'NewVariableNames','binTemp');
    % save data
    arenaData(arena).T = AD;
end  
clearvars('-except',initial_vars{:})
disp('next')

%% FIGURE: Temperature dependence of distance to food
LS = {'--','-.','-'}; %cooling|stationary|heating
for arena = 1:nArenas
    title_str = [folder ' ' expName ' ' arenaData(arena).genotype];
    fig_file = [arenaData(arena).figDir 'fly count hysteresis'];
    LW = 1.5;
    % pull data
    AD = arenaData(arena).T;
    nRates = arenaData(arena).TR.nRates;
    nTemps = arenaData(arena).TR.nTemps;
    temp_list = arenaData(arena).TR.temps;
    rate_list = arenaData(arena).TR.rates;
    temp = AD.binTemp;
    rate = AD.temprate;

    fig = figure; set(fig, 'pos', [80 218 1318 740])
    % Temp vs distance (all directions and rates)
    subplot(1,2,1); hold on
    bin_edges = threshLow:0.25:threshHigh;
    for well = 1:4
        plotData = [];
        x = AD.temperature;
        y = AD.dist2wells(:,well);
        idx = discretize(x,bin_edges);
        for ii = 1:length(bin_edges)
            loc = idx==ii;
            plotData(ii,:) = [mean(x(loc)),mean(y(loc))];
        end
        if well==arenaData(arena).foodwell
            kolor = 'green';
        else
            kolor = 'grey';
        end
        plot(plotData(:,1),plotData(:,2),'color', Color(kolor),'linewidth', LW)
    end
    ylabel('Distance (mm)')
    xlabel('Temperature (\circC)')
    title('All wells')
    
    % temp rate isolated
    subplot(1,2,2)
    hold on
    for rr = 1:nRates
        curr_rate = arenaData(arena).TR.rates(rr);
        kolor = pullFoodColor(curr_rate);
        if curr_rate>0
            lstyle = LS{3};
        elseif curr_rate<0
            lstyle = LS{1};
        elseif curr_rate==0
            continue
            lstyle = LS{2};
        end
        x = arenaData(arena).TR.temps;
        y = arenaData(arena).TR.dist_mat.avg(rr,:);
        % err = dist_mat.err(rr,:);
        plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
    end
    ylabel('Distance to food (mm)')
    xlabel('Temperature (\circC)')
    title('Only food well')

    % format figure
    fig = formatFig(fig, true, [1,2]);
    subplot(1,2,1)
    legend(strrep(arenaData(arena).wellLabels,'_',' '),'textcolor','w','location', 'best','box','off')
    subplot(1,2,2)
    for rr = 1:nRates
        str{rr} = [num2str(arenaData(arena).TR.rates(rr)) '\circC/min'];
    end
    legend(str,'textcolor','w','box','off','location', 'best')

save_figure(fig,fig_file,'-png',true);
end
clearvars('-except',initial_vars{:}) 

%% Save Data structures:
if autoSave == true
    for arena = 1:nArenas
        data = arenaData(arena);
        data.date = folder;
        data.r = r;
        data.radii = radii;
        save([arenaData(arena).figDir expName arenaIdx{arena} ' timecourse data.mat'],...
            'data','expData','expName','tempLog');
    end
    fprintf('Experiment data saved\n')
    results = 'Saved data';
%     copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')    
else
    switch questdlg('Save experiment analysis?')
        case 'Yes'
            for arena = 1:nArenas
                data = arenaData(arena);
                data.date = folder;
                data.r = r;
                data.radii = radii;
                save([arenaData(arena).figDir expName arenaIdx{arena} ' timecourse data.mat'],...
                    'data','expData','expName','tempLog');
            end
    fprintf('Experiment data saved\n')
    results = 'Saved data';
%     copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')    
        case 'No'
            results = 'No data saved';
            return
        case 'Cancel'
            results = 'No data saved';
            return
    end
end



    

    
    
    
 %% -----------------------Calculate degree of clustering------------------------------
% 
% % demo the clustering proxy
% if nvids>8
%     divisor = round(nvids/6);
%     vidList = 1:divisor:nvids;
% else
%     vidList = 1:nvids;
% end
% nrow = 2; ncol = length(vidList); ii = 0; 
% [~,minidx] = deal([]);
% 
% % VISUALIZE a demo of the clustering accuracy
% if essentialfigs == false
%     fig = getfig(''); set(fig, 'pos',[120 331 1244 368], 'color', 'k');
%     for vid = vidList
%         ii = ii+1;
%         movieInfo = VideoReader([baseFolder vidFolder '\' expName '_' num2str(vid) '.avi']); %read in video
%         headData = squeeze(data(vid).tracks(:,1,:,:));
% 
%         % Most clustered:
%         [M,minidx(ii)] = min(data(vid).flyDistance);
%         img = read(movieInfo,minidx(ii));
%         subplot(nrow, ncol, ii)
%         imshow(img); hold on
%         axis tight square
%         set(gca, 'visible', 'off')
%         x = data(vid).x_loc(minidx(ii),:); x(isnan(x)) = [];
%         y = data(vid).y_loc(minidx(ii),:); y(isnan(y)) = [];
%     %     x = squeeze(headData(minidx(ii), 1, :)); x(isnan(x)) = [];
%     %     y = squeeze(headData(minidx(ii), 2, :)); y(isnan(y)) = [];
%         scatter(x,y, 10, 'y', 'filled')
%         title(num2str(M))
%         % overlay 'size bar' for min dist:
%         plot([100,100+M], [20,20], 'linewidth', 0.5, 'color', 'r')
% 
%         % Least clustered:
%         [M,maxidx(ii)] = max(data(vid).flyDistance);
%         img = read(movieInfo,maxidx(ii));
%         subplot(nrow, ncol, ii+length(vidList))
%         imshow(img); hold on
%         axis tight square
%         set(gca, 'visible', 'off')
%         x = data(vid).x_loc(maxidx(ii),:); x(isnan(x)) = [];
%         y = data(vid).y_loc(maxidx(ii),:); y(isnan(y)) = [];
%     %     x = squeeze(headData(maxidx(ii), 1, :)); x(isnan(x)) = [];
%     %     y = squeeze(headData(maxidx(ii), 2, :)); y(isnan(y)) = [];
%         scatter(x,y, 10, 'y', 'filled')
%         title(num2str(M))
%         % overlay 'size bar' for min dist:
%         plot([100,100+M], [20,20], 'linewidth', 0.5, 'color', 'r')
%     end
%     labelHandles = findall(gcf, 'type', 'text', 'handlevisibility', 'off');
%     set(labelHandles,'FontSize', 15, 'color', 'w')
% 
%     save_figure(fig, [analysisDir expName arenaSel ' linear clustering demo'], '-pdf', autoSave);
% end
% 
% clearvars('-except',initial_vars{:})
% fprintf('Next\n')
% 
% %% ------------------------------Movement analysis------------------------------------
% %bin frame and check bin occupation changes across frames as proxy for
% %movement --> won't give an accurate 'speed' etc. but it might give
% %movement.
% tic
% nbins = 100;
% spaceChange = [];
% 
% for vid = 1:nvids   
%     N = [];
%     for frame = 1:length(data(vid).tempLog)
%         X = data(vid).x_loc(frame,:); X(isnan(X)) = [];
%         Y = data(vid).y_loc(frame,:); Y(isnan(Y)) = [];
%         N(:,:,frame) = histcounts2(X,Y,nbins);
%     end
%     if vid>1 % add last frame of previous vid to prevent frame loss 
%         N = cat(3, lastFrame, N);
%     end
%     lastFrame =  N(:,:,end);
%     binDiff = diff(N,1,3);
%     binDiff(binDiff<0) = 0;
%     BinChange = squeeze(sum(sum(binDiff,1),2));
%     % save pertinent information:
%     data(vid).BinChange = BinChange;
%     spaceChange = [spaceChange; BinChange];
% end
% occupancy.movement = spaceChange;
% toc
% 
% % preview the movement average
% if essentialfigs==false
%     fig = getfig;  hold on
%     plot(occupancy.time(2:end), occupancy.movement,'color', Color('teal'))
%     plot(occupancy.time(2:end), smooth(occupancy.movement,180),...
%         'linewidth', 2, 'color', 'w')
%     xlabel('time (min)'), ylabel('movement (a.u.)')
%     formatFig(fig, true);
%     save_figure(fig, [analysisDir expName arenaSel ' movement over time'], '-pdf', autoSave);
% end
% 
% clearvars('-except',initial_vars{:})
% fprintf('Next\n')
% 
% %% -----------------------Time course feature comparison------------------------------
% % plot parameters:
% nrow = 8; ncol = 1;
% sSpan = 180; %1 minute filter length
% sbpts(1).idx = 1;
% sbpts(2).idx = 2:4;
% sbpts(3).idx = 5:6;
% sbpts(4).idx = 7:8;
% LW = 1.5;
% 
% % FIGURE:
% fig = getfig; set(fig, 'pos', [157 86 1232 878])
% % TEMPERATURE
% subplot(nrow,ncol,sbpts(1).idx)
% plot(occupancy.time,smooth(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
% ylabel('(\circ)')
% title('temperature')
% ax = gca;
% x_lim = ax.XLim;
% 
% % OCCUPANCY
% subplot(nrow,ncol,sbpts(2).idx)
% hold on
% y = [];
% for well = 1:4
%    y = [y, smooth(occupancy.occ(:,well),sSpan)];
% end
% h = area(occupancy.time,y);
% for well = 1:4
%     h(well).FaceColor = pullFoodColor(wellLabels{well});
% end
% xlim(x_lim)
% set(gca, 'tickdir', 'out')
% l = legend(strrep(wellLabels,'_','-'));
% set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k',...
% 'position', [0.7972 0.7194 0.1039 0.0792]);% [0.8780 0.8119 0.0963 0.1126])%
% ylabel('occupancy probability')
% title('individual well occupancy')
% 
% % 
% % hold on
% % for well = 1:4
% %     kolor = pullFoodColor(wellLabels{well});
% %     plot(time,smooth(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
% % end
% % plot(time,smooth(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
% % ylabel('occupancy probability')
% % title('individual well occupancy')
% % legend(wellLabels)
% % l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
% % set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
% %     'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%
% 
% % CLUSTERING
% subplot(nrow,ncol,sbpts(3).idx); hold on
% y_avg = smooth(occupancy.IFD,sSpan);
% y_err = smooth(occupancy.IFD_err,sSpan);
% x = occupancy.time;
% fill_data = error_fill(x, y_avg, y_err);
% h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% plot(occupancy.time,smooth(occupancy.IFD,sSpan),'linewidth', LW, 'color', 'w')
% ylabel('pixels')
% title('inter-fly-distance')
% 
% % MOVEMENT
% subplot(nrow,ncol,sbpts(4).idx); hold on
% y_avg = smooth(occupancy.movement,sSpan);
% y_err = movstd(occupancy.movement,sSpan);
% x = occupancy.time(2:end);
% fill_data = error_fill(x, y_avg, y_err);
% h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% plot(x,y_avg, 'linewidth', LW, 'color', 'w')  
% ylabel('(a.u.)')
% title('movement')
% xlabel('time (min)')
% 
% formatFig(fig,true,[nrow, ncol], sbpts);
% for ii = 1:3
%     subplot(nrow,ncol,sbpts(ii).idx)
%     set(gca, 'XColor', 'k')
% end
% 
% % save and export figure
% if autoSave==true
%     export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
% else
%     if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
%         export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
%     end
% end  
% save_figure(fig, [analysisDir expName arenaSel ' timecourse features'], '-png',autoSave);
% 
% clearvars('-except',initial_vars{:})
% fprintf('Next\n')
% 
% %% ---------------------Position & Movement Summary Figure----------------------------
% % **TODO remove total well occupancy
% % figure parameters:
% nrow = 10; ncol = 1;
% sSpan = 180; %1 minute filter length
% sbpts(1).idx = 1;    % temp
% sbpts(2).idx = 2:5;  % occupancy
% sbpts(3).idx = 6:8;  % clustering & eccentricity
% sbpts(4).idx = 9:10; % movement
% LW = 1.5;
% 
% % FIGURE:
% fig = getfig; set(fig, 'pos', [157 86 1232 878])
% subplot(nrow,ncol,sbpts(1).idx)
% plot(occupancy.time,smooth(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
% ylabel('Temp(\circ)')
% % title('temperature')
% 
% subplot(nrow,ncol,sbpts(2).idx)
% hold on
% for well = 1:4
%     kolor = pullFoodColor(wellLabels{well});
%     plot(occupancy.time,smooth(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
% end
% % plot(occupancy.time,smooth(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
% ylabel('occupancy probability')
% % title('individual well occupancy')
% legend(wellLabels)
% % l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
% l = legend(strrep(wellLabels,'_','-'));
% set(l, 'color', 'k', 'textcolor', 'w','FontSize', 10,'edgecolor', 'k',...
%     'position', [0.1552 0.6918 0.0963 0.1126] );% [0.8780 0.8119 0.0963 0.1126])%
% 
% % CLUSTERING
% subplot(nrow,ncol,sbpts(3).idx); hold on
% y_avg = smooth(occupancy.IFD,sSpan);
% y_err = smooth(occupancy.IFD_err,sSpan);
% x = occupancy.time;
% fill_data = error_fill(x, y_avg, y_err);
% h(1) = fill(fill_data.X, fill_data.Y, get_color('teal'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% plot(occupancy.time,smooth(occupancy.IFD,sSpan),'linewidth', LW, 'color', Color('teal'))
% % ECCENTRICITY
% y_avg = smooth(occupancy.eccentricity(:,1),sSpan);
% y_err = smooth(occupancy.eccentricity(:,2),sSpan);
% x = occupancy.time;
% fill_data = error_fill(x, y_avg, y_err);
% h(2) = fill(fill_data.X, fill_data.Y, get_color('orange'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% plot(x,y_avg, 'linewidth', LW, 'color', Color('orange'))  
% ylabel('pixels')
% l2 = legend({'','inter fly distance', '', 'eccentricity'});
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(l2, 'color', 'k', 'textcolor', 'w','FontSize', 8,'edgecolor', 'k',...
%     'position', [0.1350 0.4638 0.1039 0.0370]);% 
% 
% % MOVEMENT
% subplot(nrow,ncol,sbpts(4).idx); hold on
% y_avg = smooth(occupancy.movement,sSpan);
% y_err = movstd(occupancy.movement,sSpan);
% x = occupancy.time(2:end);
% fill_data = error_fill(x, y_avg, y_err);
% h = fill(fill_data.X, fill_data.Y, get_color('white'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% plot(x,y_avg, 'linewidth', LW, 'color', 'w')  
% ylabel('movement (a.u.)')
% % title('movement')
% xlabel('time (min)')
% 
% formatFig(fig,true,[nrow, ncol], sbpts);
% for ii = 1:3
%     subplot(nrow,ncol,sbpts(ii).idx)
%     set(gca, 'XColor', 'k')
% end
% 
% % save and export figure
% if autoSave==true
%     export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
% else
%     if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
%         export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
%     end
% end
% save_figure(fig, [analysisDir expName arenaSel ' timecourse features'], '-png',autoSave);
% 
% clearvars('-except',initial_vars{:})
% fprintf('Next\n')
% 
% %% ----------------Average distance between each fly and each food source-------------
% %(takes a full 2 mins to run)
% 
% % find distance from center for each fly:
% tic
% FDist = [];
% idx = 0;
% for vid = 1:nvids
%     for frame = 1:length(data(vid).tempLog)
%         idx = idx+1;
%         for well = 1:4
%             test = [wellcenters(:,well)'; data(vid).x_loc(frame,:)',data(vid).y_loc(frame,:)'];
%             D = squareform(pdist(test));
%             D = D(2:end,1);
%             D(isnan(D)) = [];
%             FDist(well).N(idx,:) = [median(D), std(D)];
%         end
%     end
% end
% toc
% 
% occupancy.dist2wells = FDist;
% 
% % occupancy.eccentricity = EDist;
% nrow = 4;
% ncol = 1;
% subplotInd(1).idx = 1;
% subplotInd(2).idx = 2:4;
% % group data across videos:
% plotX = occupancy.time;
% sSpan = 180;
% LW = 1;
% 
% fig = getfig(''); 
%     subplot(nrow,ncol,subplotInd(1).idx)
%     y = smooth(occupancy.temp,sSpan);
%     plot(plotX(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'w')
%     ylabel('temp (\circC)')
%     ylim([5,27])
%     
%     subplot(nrow,ncol,subplotInd(2).idx)
%     hold on
%     % error fills
%     for well = 1:4
%         kolor = pullFoodColor(wellLabels{well}); % plotting color for food
%         y_avg(:,well) = smooth(FDist(well).N(:,1),sSpan);
%         y_err(:,well) = smooth(FDist(well).N(:,2),sSpan);
%         fill_data = error_fill(plotX, y_avg(:,well), y_err(:,well));
%         h(well) = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%         set(h(well), 'facealpha', 0.2)
%     end
%     % average line
%     for well = 1:4
%         kolor = pullFoodColor(wellLabels{well});
%         plot(occupancy.time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
%     end
%     xlabel('time (min)'); ylabel('avg distance between fly and food (a.u.)')
%     
% formatFig(fig, true, [nrow, ncol], subplotInd);
% l = legend([{'';'';'';''};strrep(wellLabels,'_','-')]);
% set(l, 'color', 'k', 'textcolor', 'w','position', [0.7947 0.6462 0.0963 0.1126])
% for well = 1:4
%     set(get(get(h(well),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% end
% subplot(nrow,ncol,subplotInd(1).idx)
% set(gca, 'XColor', 'k')
% titleName = strrep([folder ' ' expName ' Arena ' arenaSel], '_',' ');
% title(titleName,'color', 'w')
%  
% % save and export figure
% if autoSave==true
%     export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
% else
%     if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
%         export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-rgb','-append');
%     end
% end
% save_figure(fig, [analysisDir expName arenaSel ' avg distance to well'], '-png', autoSave);
% 
% clearvars('-except',initial_vars{:})
% 
% 
% %% ----------------------- Save experiment data thus far -----------------------------
% if autoSave == true
%     clearvars('-except',initial_vars{:})
%     initial_vars = unique(initial_vars);
%     save([analysisDir expName arenaSel ' timecourse data'])
%     copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')
%     fprintf('Experiment data saved\n')
%     results = 'Saved data';
% else
%     switch questdlg('Save experiment analysis?')
%         case 'Yes'
%             clearvars('-except',initial_vars{:})
%             initial_vars = unique(initial_vars);
%             save([analysisDir expName arenaSel ' timecourse data'])
%             copyfile(expPDF,'G:\My Drive\Jeanne Lab\DATA\Analysis')
%             fprintf('Experiment data saved\n')
%             results = 'Saved data';
%         case 'No'
%             results = 'No data saved';
%             return
%         case 'Cancel'
%             results = 'No data saved';
%             return
%     end
% end
