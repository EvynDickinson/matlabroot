

rootDir = 'G:\My Drive\Jeanne Lab\DATA\Literature Temperature Rate\';

literature_temprates = [5.25, 1 16.6 6000 85 7.5 168 10.2 10.5 27.5 96, 300, 120, 180, 120];
x_range = [0.8,1.2];

figure; 
x = shuffle_data(linspace(0.8,1.2,length(temprates)));
scatter(x,temprates,50,Color('gold'),'filled')
set(gca,'yscale','log')
xlim([0.7,1.3])
hold on
plot(x_range, [0.75,0.75],'r:')


%%

% fastest drop: 32C in 27 minutes = 1.185 C/min (1943 spearfish, SD)
%

highs = [87 88 61 102 75 100 100 107 87 61 95 91 88 56 76 96 97 63 70 83 67 87 89 74 73 97 101 105 91 80 92];
lows = [30 31 28 66 23 37 71 42 12 28 39 60 60 25 11 70 68 25 24 17 25 36 39 21 22 30 52 52 33 13 70];

highs_C = (highs-32).*(5/9);
lows_C = (lows-32).*(5/9);

tempRange = (highs_C-lows_C);

tempRates_24 = tempRange/(24*60);
tempRates_12 = tempRange/(12*60);
tempRates_6 = tempRange/(6*60);
tempRates_3 = tempRange/(3*60);
tempRates_1 = tempRange/(60);
tempRates_mins = tempRange/(40);

% kolors = Color('magenta','lightpink',4);
CList = [];
kolors = {'Lavender','Thistle','Plum','MediumPurple','BlueViolet','Indigo'};
for i = 1:length(kolors)
    CList(i,:) = Color(kolors{i});
end

fig = figure; set(fig,'pos',[-834 480 586 674])
hold on
xlim([0.7,1.3])
formatFig(fig,true);
set(gca,'xcolor','k')
ylabel('Temperature Rate (\circC/min)')

x = shuffle_data(linspace(0.8,1.2,length(tempRange)));
scatter(x,tempRates_24,50,CList(1,:),'filled')
% save_figure(fig,[rootDir 'TempRates_24'],'-png',true);

x = shuffle_data(linspace(0.8,1.2,length(tempRange)));
scatter(x,tempRates_12,50,CList(2,:),'filled')
% save_figure(fig,[rootDir 'TempRates_12'],'-png',true);

x = shuffle_data(linspace(0.8,1.2,length(tempRange)));
scatter(x,tempRates_6,50,CList(3,:),'filled')
% save_figure(fig,[rootDir 'TempRates_6'],'-png',true);

x = shuffle_data(linspace(0.8,1.2,length(tempRange)));
scatter(x,tempRates_3,50,CList(4,:),'filled')
% save_figure(fig,[rootDir 'TempRates_3'],'-png',true);

x = shuffle_data(linspace(0.8,1.2,length(tempRange)));
scatter(x,tempRates_1,50,CList(5,:),'filled')
% save_figure(fig,[rootDir 'TempRates_1'],'-png',true);

x = shuffle_data(linspace(0.8,1.2,length(tempRange)));
scatter(x,tempRates_mins,50,CList(6,:),'filled')
% save_figure(fig,[rootDir 'TempRates_minutes'],'-png',true);

plot(x_range, [0.75,0.75],'m:','linewidth', 1.5)
% save_figure(fig,[rootDir 'TempRates_1 with cheyenne storm line'],'-png',true);

plot(x_range, [1.18,1.18],'color', 'r','linestyle','-','linewidth', 1.5)
% save_figure(fig,[rootDir 'TempRates_1 with world record'],'-png',true);

x = shuffle_data(linspace(0.8,1.2,length(literature_temprates)));
scatter(x,literature_temprates,50,Color('gold'),'filled')
set(gca,'yscale','log')
% save_figure(fig,[rootDir 'TempRates with literature rates'],'-png',true);


plot(x_range, [0.1,0.1],'color', Color('lime'),'linestyle','-','linewidth', 1.5)
% plot(x_range, [0.16,0.16],'color', Color('lime'),'linestyle','-','linewidth', 1.5)
plot(x_range, [0.5,0.5],'color', Color('lime'),'linestyle','-','linewidth', 1.5)
% plot(x_range, [0.25,0.25],'color', Color('lime'),'linestyle','-','linewidth', 1.5)
% save_figure(fig,[rootDir 'TempRates with experiment rates'],'-png',true);

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        