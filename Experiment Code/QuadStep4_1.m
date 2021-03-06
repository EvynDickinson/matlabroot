

%% EXPANDING RAMP PROTOCOL ANALYSIS:
% Do individual ramps show hysteresis?

trial = 2;

tempPoints = getTempTurnPoints(T.TempProtocol{1});
% 
% figure; 
% plot(G(trial).TR.tempIdx)
% vline(tempPoints.transitions)

temp = G(trial).TR.tempIdx;
tempList = G(trial).TR.temps;
dist = G(trial).TR.data(:,3);
G1 = [2 3; 5,7; 5,7; 8,9; 8,9];    %cooling
G2 = [3,5; 3,5; 7,8; 7,8; 9,10];   %heating
row = 1;
col = 5;


fig = figure; set(fig, 'pos',[29 458 1746 504])
for ramp = 1:5
    subplot(row,col,ramp)
    
    % DATA SELECTION
    roiDOWN = tempPoints.transitions(G1(ramp,1)):tempPoints.transitions(G1(ramp,2));
    roiUP = tempPoints.transitions(G2(ramp,1)):tempPoints.transitions(G2(ramp,2));

    hold on
    % COOLING
    y_avg = [];
    cList = unique(temp(roiDOWN));
    cList(isnan(cList)) = [];
    y = dist(roiDOWN);
    x = temp(roiDOWN);
    % scatter(x,y,30,Color('dodgerblue'),'filled')
    for ii = 1:length(cList)
        y_avg(ii) = mean(y(x==cList(ii)));
        y_err = std(y(x==cList(ii)));
        t = tempList(cList(ii));
        scatter(t,y_avg(ii),50,Color('dodgerblue'),'filled')
        plot([t,t], [y_avg(ii)-y_err,y_avg(ii)+y_err], 'linewidth', 1.5, 'color', Color('dodgerblue'))
    end
    plot(tempList(cList),y_avg,'linewidth', 1.5, 'color', Color('dodgerblue'))

    % HEATING
    y_avg = [];
    cList = unique(temp(roiUP));
    cList(isnan(cList)) = [];
    y = dist(roiUP);
    x = temp(roiUP);
    % scatter(x,y,30,Color('red'),'filled')
    for ii = 1:length(cList)
        y_avg(ii) = mean(y(x==cList(ii)));
        y_err = std(y(x==cList(ii)));
        t = tempList(cList(ii));
        scatter(t,y_avg(ii),50,Color('red'),'filled')
        plot([t,t], [y_avg(ii)-y_err,y_avg(ii)+y_err], 'linewidth', 1.5, 'color', Color('red'))
    end
    plot(tempList(cList),y_avg,'linewidth', 1.5, 'color', Color('red'))
    
    % FORMAT
    xlabel('Temp (\circC)')
    ylabel('Distance from well (mm)')
end

formatFig(fig,true,[row,col]);
for ramp = 1:5
    subplot(row,col,ramp)
    set(gca, 'fontsize', 18)
end

save_figure(fig, [figDir 'hysteresis by ramp ' T.Genotype{trial} ' ' num2str(trial)], '-png');


%% MOVEMENT ANALYSIS: 

clearvars('-except',vars{:})

ID_mat = nan(ntrials,3);

% what are the temp protocols: 
tempList = unique(T.TempProtocol);
nTempLists = size(tempList,1);
for ii = 1:nTempLists
    ID_mat(strcmp(T.TempProtocol,tempList(ii)),1) = ii;
end

% what are the genotypes:
genoList = unique(T.Genotype);
ngenoList = size(genoList,1);
for ii = 1:ngenoList
    ID_mat(strcmp(T.Genotype,genoList(ii)),2) = ii;
end

% what are the foods/well options:
foodList = unique(T.foodCat);
nfoodList = size(foodList,1);
for ii = 1:nfoodList
    ID_mat(strcmp(T.foodCat,foodList(ii)),3) = ii;
end

% Find unique instances across all trials:
ID_List = string(ID_mat(:,1));
for ii = 2:3
    ID_List = [ID_List+string(ID_mat(:,ii))];
end
IDs = unique(ID_List);
nIDs = size(IDs,1);

% Make an empty structure in which to group the data
plotData = struct;
for ii = 1:nIDs
    plotData(ii).h = []; %heating
    plotData(ii).c = []; %cooling
    plotData(ii).r = []; %temp rate list
end

% Compare across all conditions (overlays all temp rates...)
% CList = {'BlueViolet', 'DeepPink','Orange','Lime','DodgerBlue','Teal','Red',...
%          'Turquoise', 'DarkRed', 'Indigo', 'Plum'};
colors = {'BlueViolet',...
         'Indigo',...
         'Plum',...
         'Thistle',...
         'Teal',...
         'Turquoise',...
         'Aquamarine',...
         'Red',...
         'DarkRed',...
         'Orange',...
         'Gold'};
CList = colors([1,5,10,8,2,3,6,11,9,4,7]);
     
% group together trials with the same identity
for trial = 1:ntrials
    id = find(strcmp(ID_List(trial),IDs));
    kolor = Color(CList{id});
    plotData(id).color = kolor;
    
    % collapse all temp rates into one:
    for ii = 1:G(trial).TR.nRates
        plotData(id).r = [plotData(id).r, G(trial).TR.rates(ii)];
%         temp = G(trial).TR.heatmap(ii,:); % DISTANCE
        temp = G(trial).TR.movement.avg(ii,:); % MOVEMENT
        if G(trial).TR.rates(ii)<0 % Cooling data
            plotData(id).c = [plotData(id).c; temp];
        elseif G(trial).TR.rates(ii)>0 % Heating data
            plotData(id).h = [plotData(id).h; temp];
        end
    end
end


%% FIGURE: overlay movement tracks for all trials in the structure color coded by heating and cooling
% plot all trials individually
sSpan = 1; %3-degree smoothing
id = 1;

fig = figure;
    hold on
    % Cooling
    for id = 1:nIDs
        x = G(1).TR.temps; %all trials should have same temp ids
        y = plotData(id).c;
        for ii = 1:size(y,1)
            plot(x,smooth(y(ii,:),sSpan),'color', Color('dodgerblue')) %plotData(id).color) %
        end
    end
    % Heating
    for id = 1:nIDs
        x = G(1).TR.temps; %all trials should have same temp ids
        y = plotData(id).h;
        for ii = 1:size(y,1)
            plot(x,smooth(y(ii,:),sSpan),'color', 'r') %plotData(id).color) %
        end
    end
% xlim([7,22])
xlabel('Temperature (\circC)')
% ylabel('Distance')
ylabel('Movement')
formatFig(fig, true);
h_line(5,'w',1)
v_line([14,22],'gold',1)
set(gca, 'fontsize', 18)

save_figure(fig, [figDir 'movement hysteresis loops all trials'], '-png');


%%

% Characterize movement hysteresis:
SZ = 75;
LW = 2;
sSpan = 8;

for trial = 1:ntrials
    
    cool = G(trial).TR.movement.avg(1,:); % heating
    heat = G(trial).TR.movement.avg(2,:); % heating
    temp = G(trial).TR.temps;
    % cutoff temp and heat cool etc by location
    loc = isnan(cool);
    cool(loc) = [];
    heat(loc) = [];
    temp(loc) = [];
    
    % find peak movement in the cooling phase
    [cool_max,Cmax_loc] = max(smooth(cool,sSpan));
    [cool_min,Cmin_loc] = min(smooth(cool,sSpan));
    % find peak movement in the heating phase
    [heat_max,Hmax_loc] = max(smooth(heat,sSpan));
    [heat_min,Hmin_loc] = min(smooth(heat,sSpan));
    % quarter-movement line
    quarter_line = cool_min + (cool_max-cool_min)*0.25;
    crossPoint = find(smooth(heat,sSpan)>=quarter_line);
    crossPoint = crossPoint(1)-1;
    
    fig = figure; hold on
        %heating
        plot(temp,smooth(heat,sSpan),'color', Color('orangered'),'linewidth', LW) 
        plot(temp,heat,'color', Color('orangered'),'linewidth', LW,'linestyle', ':') 
        %cooling
        plot(temp,smooth(cool,sSpan),'color', Color('dodgerblue'),'linewidth', LW) 
        plot(temp,cool,'color', Color('dodgerblue'),'linewidth', LW,'linestyle', ':') 
        % plot the peak speed in the cooling phase
        scatter(temp(Cmax_loc),cool_max,SZ,Color('yellow'),'filled','d')
        scatter(temp(Cmin_loc),cool_min,SZ,Color('yellow'),'filled','d')
%         % plot the peak speed in the heating phase
%         scatter(temp(Hmax_loc),heat_max,SZ,Color('orangered'),'filled','d')
%         scatter(temp(Hmin_loc),heat_min,SZ,Color('orangered'),'filled','d')
        % quarter-movement line
        h_line(quarter_line, 'yellow',1)
        v_line(temp(crossPoint), 'yellow',1)
        
        % labels and formatting
        xlabel('Temperature (\circC)')
        ylabel('Movement')
        formatFig(fig,true);
        set(gca, 'fontsize', 18)
        
        save_figure(fig, [figDir 'movement hysteresis loop recovery ' num2str(trial)], '-png');

        
% plot the hysteresis...
mov_diff = cool-heat;

fig = figure; hold on
    plot(temp,mov_diff, 'Color', Color('teal'),'linewidth', LW,'linestyle', ':') 
    plot(temp,smooth(mov_diff,sSpan), 'Color', Color('teal'),'linewidth', LW) 
    h_line(0,'grey',1)
    xlabel('Temperature (\circC)')
    ylabel('cool-heat difference')
    formatFig(fig,true);
    set(gca, 'fontsize', 18)
        

end



%%




    
% Average across trials and 'uncode' the parameter names
for id = 1:nIDs
    plotData(id).c_avg = mean(plotData(id).c,1,'omitnan');
    plotData(id).h_avg = mean(plotData(id).h,1,'omitnan');
    % decode parameter names
    a = char(IDs(id));
    plotData(id).name = [tempList{str2double(a(1))} ' ' genoList{str2double(a(2))} ' ' foodList{str2double(a(3))}];
end
    
% FIGURE: Plot the temp-rate food proximity tuning curves
temp = G(1).TR.temps; %all trials should have same temp ids
LW = 1.5;
lStr = [];
fig = figure; set(fig, 'pos', [210 121 977 660])
    hold on
    ii = 1;
    for id = 1:nIDs
        plot(temp, plotData(id).c_avg,'color', plotData(id).color,'linewidth',LW,'linestyle', '--')
        plot(temp, plotData(id).h_avg,'color', plotData(id).color,'linewidth',LW,'linestyle', '-')
        % legend string
        if nIDs<5
            lStr{ii} = [plotData(id).name ' cooling'];
            lStr{ii+1} = [plotData(id).name ' heating'];
            ii = ii+2;
        else
            lStr{ii} = plotData(id).name;
            lStr{ii+1} = '';
            ii = ii+2;
        end
    end
    %labels and formatting
    xlabel('Temperature (\circC)')
    ylabel('Movement (mm)')
    formatFig(fig, true);
    set(gca,'fontsize', 18)
    xlim([floor(temp(1)-1),temp(end)+ceil(range(temp)*0.33)])
    %legend
    lStr = strrep(lStr,'_',' ');
    legend(lStr,'textcolor', 'w', 'box', 'off','fontsize', 8)
    
save_figure(fig, [figDir 'movement tuning overlay'], '-png');

% TODO: get measure of movement hysteresis and recovery from cooling in
% flies


























    
    
    
    
    
    
    
    
    
    










