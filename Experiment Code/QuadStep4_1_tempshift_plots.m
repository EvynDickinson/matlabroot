

%% Correlation plots

clear; close all; clc
baseFolder = getCloudPath;  

list_dirs = dir([baseFolder 'Grouped Data Structures\']);
list_dirs = {list_dirs(:).name};
list_dirs(1:2) = [];
dirIdx = listdlg('ListString', list_dirs, 'SelectionMode', 'single','ListSize',[300,450]);
expGroup = list_dirs{dirIdx}; %name of experiment groups selected
saveDir = [baseFolder 'Grouped Data Structures\' expGroup '\'];
load([saveDir expGroup ' temp distance correlation ramps only.mat']);
disp([expGroup ' loaded'])


num.exp = length(plotData);
blkbgd = true;
[foreColor,backColor] = formattingColors(blkbgd);
buff = 0.2;
SZ = 50;
LW = 1.5;

pixWidth = 60; % additional pixel size for image for each extra experiment group
figSize = [pixWidth + (pixWidth*num.exp),590];

% get correlation data
for i = 1:num.exp
    disp(plotData(i).groupName)
end

%% Manually check the temperature shift-order

expOrder = [3,1,2];
colors = {'dodgerblue','gold','red'};

ylimits = [-0.8,0.1];

%display proper order:
for exp = expOrder
    disp(plotData(exp).groupName)
end

% correlation coefficients
fig = getfig('',true,figSize); 
hold on
 for ii = 1:num.exp
   
   i = expOrder(ii);
   kolor = Color(colors{ii});
   xlow = ii-buff-0.1;
   xhigh = ii+buff+0.1;
   y = plotData(i).rho;
   y_avg = mean(plotData(i).rho);
   x = shuffle_data(linspace(ii-buff,ii+buff,length(y)));
   
   scatter(x,y,SZ,kolor,'filled')
   plot([xlow,xhigh],[y_avg,y_avg],'color',kolor,'linewidth',LW)
 end
 xlim([0.5,num.exp+.5])       
 ylabel('temp-distance corr. coef.')
 h_line(0,foreColor,':',1)    
 formatFig(fig,blkbgd);    
 set(gca,'xcolor',backColor)
 ylim(ylimits)
 
% save figure
save_figure(fig,[saveDir expGroup ' temp distance correlation ramps only'],'-png',false,false);  
save_figure(fig,[saveDir expGroup ' temp distance correlation ramps only'],'-pdf',false,true);  







