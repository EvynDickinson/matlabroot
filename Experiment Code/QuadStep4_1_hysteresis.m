
%% HYSTERESIS ANALYSIS

%% PULL DATA 

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
        temp = G(trial).TR.heatmap(ii,:); % DISTANCE
%         temp = G(trial).TR.movement.avg(ii,:); % MOVEMENT
        if G(trial).TR.rates(ii)<0 % Cooling data
            plotData(id).c = [plotData(id).c; temp];
        elseif G(trial).TR.rates(ii)>0 % Heating data
            plotData(id).h = [plotData(id).h; temp];
        end
    end
end

newvars = who; newvars{end+1} = 'newvars';

%% FIGURE: overlay distance lines for all trials in the structure color coded by heating and cooling
% plot all trials individually
sSpan = 1; %3-degree smoothing
LW = 2;
fig = figure;
    hold on
    % Cooling
    for id = 1:nIDs
        x = G(1).TR.temps; %all trials should have same temp ids
        y = plotData(id).c;
        for ii = 1:size(y,1)
            plot(x,smooth(y(ii,:),sSpan),'color', Color('dodgerblue'),'linewidth',LW) %plotData(id).color) %
        end
    end
    % Heating
    for id = 1:nIDs
        x = G(1).TR.temps; %all trials should have same temp ids
        y = plotData(id).h;
        for ii = 1:size(y,1)
            plot(x,smooth(y(ii,:),sSpan),'color', 'r','linewidth', LW) %plotData(id).color) %
        end
    end
% xlim([7,22])
xlabel('Temperature (\circC)')
ylabel('Distance (mm)')
formatFig(fig, true);
% h_line(5,'w',1)
% v_line([14,22],'gold',1)
set(gca, 'fontsize', 18)

save_figure(fig, [figDir 'distance hysteresis loops all trials'], '-png');

% % POSTER CODE FIGURE: overlay distance lines for all trials in the structure color coded by heating and cooling
% % plot all trials individually
% [groupH, groupC] = deal([]);
% sSpan = 1; % 3-degree smoothing
% LW = 0.5;
% fig = figure; set(fig, 'pos',[-768 538 341 645])
%     hold on
%     % Cooling
%     for id = 1:nIDs
%         x = G(1).TR.temps; %all trials should have same temp ids
%         y = plotData(id).c;
%         groupC = [groupC; y];
%         for ii = 1:size(y,1)
%             plot(x,smooth(y(ii,:),sSpan),'color', Color('dodgerblue'),'linewidth',LW) %plotData(id).color) %
%         end
%     end
%     % Heating
%     for id = 1:nIDs
%         x = G(1).TR.temps; %all trials should have same temp ids
%         y = plotData(id).h;
%         groupH = [groupH; y];
%         for ii = 1:size(y,1)
%             plot(x,smooth(y(ii,:),sSpan),'color', 'r','linewidth', LW) %plotData(id).color) %
%         end
%     end
%     % Plot average lines
%     plot(x,mean(groupC),'color', Color('navy'),'linewidth', 2.5) %plotData(id).color) %
%     plot(x,mean(groupH),'color', Color('darkred'),'linewidth', 2.5) %plotData(id).color) %
%     
%     v_line(12,'black','--',1)
%     
% % xlim([7,22])
% xlabel('temperature (\circC)')
% ylabel('distance from food (mm)')
% formatFig(fig);
% % h_line(5,'w',1)
% % v_line([14,22],'gold',1)
% set(gca, 'fontsize', 15)
% 
% save_figure(fig, [figDir 'distance hysteresis loops all trials'], '-pdf');



%% Plot the hysteresis values for all trials


% % POSTER CODE FIGURE: overlay distance lines for all trials in the structure color coded by heating and cooling
% % plot all trials individually
% [groupH, groupC] = deal([]);
% sSpan = 1; % 3-degree smoothing
% LW = 0.5;
% fig = figure; set(fig, 'pos',[-768 538 341 645])
%     hold on
%     % Cooling
%     for id = 1:nIDs
%         x = G(1).TR.temps; %all trials should have same temp ids
%         y = plotData(id).c;
%         groupC = [groupC; y];
%         for ii = 1:size(y,1)
%             plot(x,smooth(y(ii,:),sSpan),'color', Color('dodgerblue'),'linewidth',LW) %plotData(id).color) %
%         end
%     end
%     % Heating
%     for id = 1:nIDs
%         x = G(1).TR.temps; %all trials should have same temp ids
%         y = plotData(id).h;
%         groupH = [groupH; y];
%         for ii = 1:size(y,1)
%             plot(x,smooth(y(ii,:),sSpan),'color', 'r','linewidth', LW) %plotData(id).color) %
%         end
%     end
%     % Plot average lines
%     plot(x,mean(groupC),'color', Color('navy'),'linewidth', 2.5) %plotData(id).color) %
%     plot(x,mean(groupH),'color', Color('darkred'),'linewidth', 2.5) %plotData(id).color) %
%     
%     v_line(12,'black','--',1)
%     
% % xlim([7,22])
% xlabel('temperature (\circC)')
% ylabel('distance from food (mm)')
% formatFig(fig);
% % h_line(5,'w',1)
% % v_line([14,22],'gold',1)
% set(gca, 'fontsize', 15)


%% Characterize food distance hysteresis:
clearvars('-except',newvars{:})

for trial = 1:ntrials
    
    cool = G(trial).TR.heatmap(1,:); % heating
    heat = G(trial).TR.heatmap(2,:); % heating
    temp = G(trial).TR.temps;
    % cutoff temp and heat cool etc by location
    loc = isnan(cool);
    cool(loc) = [];
    heat(loc) = [];
    temp(loc) = [];

    % Hysteresis
    mov_diff = cool-heat;
    
    % save data
    G(trial).mov_diff = [temp',mov_diff'];
    G(trial).color = pullGenotypeColor(T.Genotype{trial});
end

% FIGURE: hysteresis overlay across all trials colored by genotype
buff = 0.2;
SZ = 75;
LW = 2;
sSpan = 8;
row = 1;
col = 2;

% pull the 'hysteresis index'
for ii = 1:ngenoList
    hdata(ii).s = [];
end 
for trial = 1:ntrials
    idx = find(strcmp(T.Genotype{trial},genoList));
    hdata(idx).s = [hdata(idx).s; sum(G(trial).mov_diff(:,2))];
end


fig = figure; set(fig,'pos',[174 207 872 537])
% hystereis by time
subplot(row,col,1); hold on
    for trial = 1:ntrials
        temp = G(trial).mov_diff;
%         plot(temp(:,1),temp(:,2),'color', G(trial).color,'linewidth', LW) 
        plot(temp(:,1),temp(:,2),'color', 'w','linewidth', LW) 
    end
    xlabel('temperature (\circC)')  
    ylabel('cool-warm distance hystersis (mm)')
    h_line(0,'white',':',1)
% hysteresis scatter plot
subplot(row,col,2); hold on 
    for ii = 1:ngenoList
%         kolor = pullGenotypeColor(genoList{ii});
        kolor = 'w';
        x = shuffle_data(linspace(ii-buff,ii+buff,length(hdata(ii).s)));
        scatter(x,hdata(ii).s,SZ,kolor,'filled')
        plot([ii-(buff*1.25),ii+(buff*1.25)],[mean(hdata(ii).s),mean(hdata(ii).s)],'color',kolor,'linewidth',LW)
    end
    ylabel('total distance difference (mm)')
    xlabel('Genotype')
    set(gca,'XTick',1:ngenoList,'XTickLabel',strrep(genoList,'_',' '),'XTickLabelRotation',15)
    h_line(0,'white',':',1)
    
    ylim([-200,300])
formatFig(fig,true,[row,col]);
save_figure(fig, ['G:\My Drive\Presentations\SRT May 2022\distance hysteresis scatter'], '-pdf');
% save_figure(fig, [figDir 'distance hysteresis scatter'], '-png');

% % % POSTER CODE
% fig = figure; set(fig,'pos',[-980 374 314 861])
% % hystereis by time
% subplot(2,1,1); hold on
%     for trial = 1:ntrials
%         temp = G(trial).mov_diff;
%         plot(temp(:,1),temp(:,2),'color', G(trial).color,'linewidth', LW) 
%     end
%     xlabel('temperature (\circC)')
%     ylabel('cool-warm distance hystersis (mm)')
%     h_line(0,'black',':',1)
%     axis tight
% % hysteresis scatter plot
% subplot(2,1,2); hold on 
%     for ii = 1:ngenoList
%         kolor = pullGenotypeColor(genoList{ii});
%         x = shuffle_data(linspace(ii-buff,ii+buff,length(hdata(ii).s)));
%         scatter(x,hdata(ii).s,SZ,kolor,'filled')
%         plot([ii-(buff*1.25),ii+(buff*1.25)],[mean(hdata(ii).s),mean(hdata(ii).s)],'color',kolor,'linewidth',LW)
%     end
%     ylabel('total distance difference (mm)')
%     xlabel('Genotype')
%     set(gca,'XTick',1:ngenoList,'XTickLabel',strrep(genoList,'_',' '),'XTickLabelRotation',15)
%     h_line(0,'black',':',1)
% formatFig(fig,false,[2,1]);
% save_figure(fig, [figDir 'distance hysteresis scatter'], '-pdf');




%% FIGURE -- single trial hysteresis outlines
clearvars('-except',vars{:})

SZ = 75;
LW = 2;
sSpan = 8;
row = 1;
col = 2;

trial = 24;

cool = G(trial).TR.heatmap(1,:); % heating
heat = G(trial).TR.heatmap(2,:); % heating
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
% Hysteresis
mov_diff = cool-heat;

fig = figure;  set(fig,'pos',[292 299 911 606])
subplot(row,col,1)
hold on
    %heating
    plot(temp,smooth(heat,sSpan),'color', Color('orangered'),'linewidth', LW) 
    plot(temp,heat,'color', Color('orangered'),'linewidth', LW,'linestyle', ':') 
    %cooling
    plot(temp,smooth(cool,sSpan),'color', Color('dodgerblue'),'linewidth', LW) 
    plot(temp,cool,'color', Color('dodgerblue'),'linewidth', LW,'linestyle', ':') 
    h_line(quarter_line, 'yellow',1)

    % labels and formatting
    xlabel('Temperature (\circC)')
    ylabel('Distance to food (mm)')
    title(strrep(T.Genotype{trial},'_',' '))

% plot the hysteresis...
subplot(row,col,2); hold on
    % fig = figure; hold on
    plot(temp,mov_diff, 'Color', Color('teal'),'linewidth', LW,'linestyle', ':') 
    plot(temp,smooth(mov_diff,sSpan), 'Color', Color('teal'),'linewidth', LW) 
    h_line(0,'grey',1)
    xlabel('Temperature (\circC)')
    ylabel('cool-heat difference')
    title(['Sum: ' num2str(sum(mov_diff))])
    formatFig(fig,true,[row,col]);

save_figure(fig, [figDir 'distance hysteresis loop recovery ' num2str(trial)], '-png');



