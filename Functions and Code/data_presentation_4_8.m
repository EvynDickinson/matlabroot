



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






