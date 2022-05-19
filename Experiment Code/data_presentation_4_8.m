



%% Overlay various temperature protocols

baseFolder = getCloudPath;
tempList = {'03.24.2022','02.01.2022','02.02.2022','02.03.2022','02.04.2022',...
            '02.05.2022','02.06.2022','01.12.2022','01.20.2022','01.27.2022',...
            '01.25.2022','01.06.2022','01.11.2022','03.06.2022','03.05.2022',...
            '03.09.2022','02.22.2022','02.20.2022','03.23.2022','03.31.2022',...
            '04.02.2022'};
        



for ii = 1:length(tempList)

    filePath = [baseFolder tempList{ii}];

    list_dirs = dir([filePath '/*dataMat.mat']);
    
   % Load data:
    tempLog = readmatrix([filePath '/' list_dirs.name(1:end-12) '_RampLog']);
    expData = load([filePath '/' list_dirs.name]);
    
    % Pull tempdata from important places only:
    logROI(1) = find(tempLog(:,1)==expData.tempLogStart(1,3));
    logROI(2) = find(tempLog(:,1)==expData.tempLogEnd(expData.parameters.numVids,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    
    data(ii).temp = tempCourse;
    
end



% Plot time courses
fig = figure; hold on
for ii = 1:length(tempList)
    y = data(ii).temp;
    plot(x(1:length(y)), y,'linewidth', 1.5)
end
axis tight
formatFig(fig,true);
xlabel('time (hr)')
ylabel('temperature (\circC)')
set(gca,'fontsize', 20)
ylim([5 35])
save_figure(fig,[baseFolder(1:end-5) '/Presentations/Data Presentation 04.08.2022/all temp ramps'],'-png');

% xlimsssss
x_lim = xlim;
x = 1:x_lim(2);
x = ((x*5)./60)./60;

%%


fig = figure; hold on
y = data(12).temp;
plot(x(1:length(y)), y,'linewidth', 1.5,'color',Color('magenta'))

y = data(13).temp;
plot(x(1:length(y)), y,'linewidth', 1.5,'color',Color('teal'))

y = data(5).temp;
plot(x(1:length(y)), y,'linewidth', 1.5,'color',Color('yellow'))

xlim([0.0014 14.2597])
ylim([5 35])
formatFig(fig,true);
xlabel('time (hr)')
ylabel('temperature (\circC)')
set(gca,'fontsize', 20)

save_figure(fig,[baseFolder(1:end-5) '/Presentations/Data Presentation 04.08.2022/temp ramp 3'],'-png');






%% FIGURE: hysteresis index for Velocity Ramps A-E
[r1,r2,r3] = deal([]);
for trial = 1:ntrials
    
    r1(trial) = sum(G(trial).heatmap(1,:)-G(trial).heatmap(7,:),'omitnan');
    r2(trial) = sum(G(trial).heatmap(2,:)-G(trial).heatmap(6,:),'omitnan');
    r3(trial) = sum(G(trial).heatmap(3,:)-G(trial).heatmap(5,:),'omitnan');
end

buff = 0.3;
SZ = 75;
LW = 1.5;
fig = figure; hold on
    %rate 0.5
    x = shuffle_data(linspace(1-buff,1+buff,ntrials));
    scatter(x,r1,SZ,Color('orange'),'filled')
    plot([1-buff,1+buff],[mean(r1),mean(r1)],'color', Color('orange'),'linewidth', LW)
    %rate 0.25
    x = shuffle_data(linspace(2-buff,2+buff,ntrials));
    scatter(x,r2,SZ,Color('darkviolet'),'filled')
    plot([2-buff,2+buff],[mean(r2),mean(r2)],'color', Color('darkviolet'),'linewidth', LW)
    %rate 0.1
    x = shuffle_data(linspace(3-buff,3+buff,ntrials));
    scatter(x,r3,SZ,Color('turquoise'),'filled')
    plot([3-buff,3+buff],[mean(r3),mean(r3)],'color', Color('turquoise'),'linewidth', LW)
    %labels and formatting
    h_line(0,'grey', ':', 2)
    xlabel('temp rate (\circC/min)')
    ax = gca;
    ax.XTick = 1:3;
    ax.XTickLabels = {'0.5','0.25','0.1'};

    formatFig(fig,true);

save_figure(fig, [figDir 'hysteresis summary'], '-png');
% 
% % POSTER CODE:
% buff = 0.3;
% SZ = 75;
% LW = 2;
% fig = figure; set(fig, 'pos',[-768 538 341 645]); hold on
%     %rate 0.5
%     x = shuffle_data(linspace(1-buff,1+buff,ntrials));
%     scatter(x,r1,SZ,Color('orange'),'filled')
%     plot([1-buff,1+buff],[mean(r1),mean(r1)],'color', Color('orange'),'linewidth', LW)
%     %rate 0.25
%     x = shuffle_data(linspace(2-buff,2+buff,ntrials));
%     scatter(x,r2,SZ,Color('darkviolet'),'filled')
%     plot([2-buff,2+buff],[mean(r2),mean(r2)],'color', Color('darkviolet'),'linewidth', LW)
%     %rate 0.1
%     x = shuffle_data(linspace(3-buff,3+buff,ntrials));
%     scatter(x,r3,SZ,Color('teal'),'filled')
%     plot([3-buff,3+buff],[mean(r3),mean(r3)],'color', Color('teal'),'linewidth', LW)
%     %labels and formatting
%     h_line(0,'black', ':', 2)
%     xlabel('temp rate (\circC/min)')
%     ax = gca;
%     ax.XTick = 1:3;
%     ax.XTickLabels = {'0.5','0.25','0.1'};
%     ylabel('cumulative hysteresis (mm)')
%     formatFig(fig,false);
%     xlim([0.4,3.6])
%     set(gca, 'fontsize', 15)
% 
% save_figure(fig, [figDir 'hysteresis summary'], '-pdf');


%% Velocity ramps A-E temperature hysteresis plots by TEMP RATE 

HM = [];
% Group across the trials by temp and temp rate
for trial = 1:ntrials
    HM(:,:,trial) = G(trial).heatmap;
end
% Average across trials
avgHM = mean(HM,3,'omitnan');
avgHM(4,:) = [];   % remove the holding temp data
avgHM(:,end) = []; % remove the nan at the highest temp row
rates = G(1).rates;
rates(4) = [];
temp = G(1).temps(1:end-1);

LW = 2;
CList = {'orange','darkviolet','teal','teal','darkviolet','orange'};
lStyle = {'--','--','--','-','-','-'};

fig = figure; set(fig, 'pos',[-768 538 341 645]); hold on
for ii = 1:length(rates)
   plot(temp,avgHM(ii,:),'color',Color(CList{ii}),'linestyle',lStyle{ii},'linewidth',LW) 
end
xlabel('temperature (\circC)')
ylabel('distance from food (mm)')
formatFig(fig);
set(gca, 'fontsize', 15)
save_figure(fig, [figDir 'distance vs temp by rate summary'], '-pdf');


%% Velocity ramps A-E temperature hysteresis plots by HEAT/COOL
[HM_cool,HM_warm] = deal([]);
% Group across the trials by temp and temp rate
for trial = 1:ntrials
    HM_cool = [HM_cool; G(trial).heatmap(1:3,:)];
    HM_warm = [HM_warm; G(trial).heatmap(5:7,:)];
end

% Average across trials
avgCool = mean(HM_cool);
errCool = std(HM_cool)./sqrt(ntrials);
avgWarm = mean(HM_warm);
errWarm = std(HM_warm)./sqrt(ntrials);

temp = G(1).temps;

LW = 2;

fig = figure; set(fig, 'pos',[-931 749 689 429]); hold on
    % warming
    plot(temp,avgWarm,'color',Color('darkred'),'linewidth',2) 
    plot(temp,avgWarm+errWarm,'color',Color('darkred'),'linewidth',0.5)
    plot(temp,avgWarm-errWarm,'color',Color('darkred'),'linewidth',0.5)
    % cooling
    plot(temp,avgCool,'color',Color('navy'),'linewidth',2) 
    plot(temp,avgCool+errCool,'color',Color('navy'),'linewidth',0.5)
    plot(temp,avgCool-errCool,'color',Color('navy'),'linewidth',0.5)
    
xlabel('temperature (\circC)')
ylabel('distance from food (mm)')
formatFig(fig);
set(gca, 'fontsize', 15)
save_figure(fig, [figDir 'distance vs temp by temp direction'], '-pdf');



%% Single trial plot data
trial = 1;

nrow = 4;
ncol = 1;
subplotInd(1).idx = 1;
subplotInd(2).idx = 2:4;
sSpan = 180;
LW = 1.5;
wellLabels = data(trial).wellLabels;
temp = data(trial).occupancy.temp;
time = data(trial).occupancy.time;
dist =  data(trial).dist2wells;

% MAKE THE FIGURE:
lbl = [];
fig = getfig(''); set(fig,'pos', [-917 254 628 846]);
    % TEMPERATURE data
    subplot(nrow,ncol,subplotInd(1).idx)
    y = smooth(temp,sSpan);
    plot(time(2:end-1), y(2:end-1), 'linewidth', LW, 'color', 'k')
    ylabel('temp (\circC)')
%         ylim([5,27])
    % OCCUPANCY
    subplot(nrow,ncol,subplotInd(2).idx)
    hold on
    % error fills
    for well = 1:4
        if strcmp(wellLabels{well},'Empty')
            kolor = Color('grey');
            lbl{well} = 'empty';
        else 
            kolor = Color('purple');
            lbl{well} = T.foodCat{trial};
        end
        
        y_avg(:,well) = smooth(dist(:,well),sSpan);
%         y_err(:,well) = movstd(dist(:,well),sSpan);
%         fill_data = error_fill(time, y_avg(:,well), y_err(:,well));
%         h(well) = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%         set(h(well), 'facealpha', 0.3)
        plot(time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
    end
%     % average line
%     for well = 1:4
%         kolor = pullFoodColor(wellLabels{well});
%         plot(time,y_avg(:,well), 'linewidth', LW, 'color', kolor);
%     end
    xlabel('time (min)'); ylabel('distance to well (mm)')
% FORMATTING
formatFig(fig, false, [nrow, ncol], subplotInd);
l = legend(lbl);
set(l, 'box','off','textcolor', 'k')
subplot(nrow,ncol,subplotInd(1).idx)
set(gca, 'XColor', 'w')

% Save image
save_figure(fig, [figDir 'trial ' num2str(trial) ' distance temp plot'], '-pdf');



save_figure(fig, [arenaData(arena).figDir expName ' ' arenaSel ' well occupation timcourse'], '-png', autoSave);

clearvars('-except',initial_vars{:})
fprintf('\nNext\n')

















