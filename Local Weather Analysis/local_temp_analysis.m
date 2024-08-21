

%% Compile local data logs for a general analysis

clear; clc; close all

% Load fast logging temperature data
folder = getCloudPath;
folder = [folder(1:end-5) 'Temp logging/'];
fileList = dir([folder '*.mat']);
fileList = {fileList(:).name};
idx = listdlg("ListString",fileList,'PromptString','Select the temp log to load', 'ListSize',[300,300],'SelectionMode','multiple');
if isempty(idx)
    return
end
data = [];
num.trial = length(idx);
for i = 1:num.trial
    fileName = fileList{idx(i)};
    name = fileName(14:end-4);
    temp = load([folder fileName]);
    data(i).name = name;
    data(i).filename = fileName;
    data(i).T = temp.T;
end
clear temp

tr = [];
for i = 1:num.trial
    tr = [tr; data(i).T.tempRate];
end
tr(isnan(tr)) = [];

%% Large temp rate histogram
kolor = Color('teal');

% calculate the empircal distribution of the temp rate data
[fx,x] = ecdf(tr);

% Calculate mean and standard deviation
mu = mean(tr,'omitnan');
sigma = std(tr);

% % Calculate the 90% confidence interval boundaries
% lower_bound = mu - 1.645 * sigma;
% upper_bound = mu + 1.645 * sigma;
% 
% % Calculate the 95% confidence interval boundaries
% lower_bound = mu - 1.96 * sigma;
% upper_bound = mu + 1.96 * sigma;

% Calculate the 50% confidence interval boundaries
lower_bound = mu - 0.647 * sigma;
upper_bound = mu + 0.647 * sigma;

LW = 2;

fig = figure; hold on
    ylim([0,1])
    v_line([upper_bound, lower_bound],'r',':',LW)
    % v_line(mu,'r',':',LW)
    plot(x, fx,'color', 'k', 'linewidth', LW)
    xlabel('rate of temperature change (\circC/min)')
    ylabel('cumulative probability')
    formatFig(fig, false);
    xlim([-25,25])

save_figure(fig, [folder '/Figures/Cumulative Probability Density'],'-pdf');

%% Fit distribution of temperature rate

pd = fitdist(tr,'Normal');
x = -25:25;
y = pdf(pd,x);

% Calculate the 50% confidence interval boundaries
lower_bound = pd.mu - 0.647 * pd.sigma;
upper_bound = pd.mu + 0.647 * pd.sigma;

fig = figure; hold on
    yyaxis right
    h = histogram(tr, 'FaceAlpha',1,'FaceColor',Color('grey'));
    yyaxis left
    plot(x,y, 'Color', 'k', 'LineWidth',3)
    axis tight
    formatFig(fig, false);
    yyaxis right
    set(gca, 'ycolor', 'w')
    xlim([-25,25])
    v_line([upper_bound, lower_bound],'r',':',2)
    xlabel('temperature rate (\circC/min)')
    yyaxis left
    ylabel('Probability Density')

save_figure(fig, [folder '/Figures/Probability Density'],'-pdf');

