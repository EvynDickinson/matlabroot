

%% Estimate the sample sizes needed for proper power?
clear; clc

alpha_level = 0.05; % significance level (~p-value)
power = 0.80;       % Desired power
testType = 't';    % Single-sided T-Test
PN_SR = 100; % ~PN spike rate 100Hz ORN spike rate in intact flies
PN_std = 15; % estimated STD for the PN spike rate
PN_noLN = 140; % ~PN spike rate for 100Hz ORN when lateral inhibition is removed from antennal segmentation
% allow for multiple comparisons corrections later
numComparisons = 9;
alphaPerTest = alpha_level / numComparisons;


% Function call to compute the required sample size per group
n = sampsizepwr(testType, [PN_SR PN_std], PN_noLN, power,[],'Alpha', alphaPerTest);

% Display the result
fprintf('Required sample size per group: %d\n', ceil(n));


% plot the power over different sample sizes
nn = 1:20;
pwrout = sampsizepwr(testType, [PN_SR PN_std], PN_noLN, [],nn,'Alpha', alphaPerTest);

fig = figure; hold on
plot(nn,pwrout,'b-','linewidth', 2)
h_line(0.8,'r',':')
v_line(n,'r', ':')
title(['Power versus Sample Size with ' num2str(numComparisons) ' comparisons'])
xlabel('Sample Size')
ylabel('Power')
formatFig(fig, false);


%% Estimate the sample sizes needed for proper power?
clear; clc

alpha_level = 0.05; % significance level (~p-value)
power = 0.80;       % Desired power
testType = 't';    % Single-sided T-Test
PN_SR = 60; % ~PN spike rate 100Hz ORN spike rate in intact flies
PN_std = 7; % estimated STD for the PN spike rate
PN_noLN = 80; % ~PN spike rate for 100Hz ORN when lateral inhibition is removed from antennal segmentation
% allow for multiple comparisons corrections later
numComparisons = 9;
alphaPerTest = alpha_level / numComparisons;


% Function call to compute the required sample size per group
n = sampsizepwr(testType, [PN_SR PN_std], PN_noLN, power,[],'Alpha', alphaPerTest);

% Display the result
fprintf('Required sample size per group: %d\n', ceil(n));


% plot the power over different sample sizes
nn = 1:20;
pwrout = sampsizepwr(testType, [PN_SR PN_std], PN_noLN, [],nn,'Alpha', alphaPerTest);

fig = figure; hold on
plot(nn,pwrout,'Color', Color('green'),'linewidth', 2)
h_line(0.8,'r',':')
v_line(n,'r', ':')
% title(['Power versus Sample Size with ' num2str(numComparisons) ' comparisons'])
xlabel('Sample Size')
ylabel('Power')
formatFig(fig, false);


%% Quick hypothesis testing:
clear

folder = getCloudPath;
figDir = [folder, 'Electrophysiology Modeling\'];

% raw start data
T = 15:5:35;
PN_base = [0.1, 0.5, 0.7, 0.8, 0.85];
TRN_base = [1, 0.25, 0, 0.25, 1];

% upsample the data
N = 100;
temp = linspace(15,35,N);
PN = interp1(T,PN_base, temp,'spline');
TRN = interp1(T,TRN_base, temp,'spline');

% determine the effect of subtractive inhibition with rectification
PN_TRN = PN - TRN;
PN_TRN(PN_TRN<0) = 0; % rectify

% Plot base figure:
LW = 2;
fig = figure; hold on
    plot(temp, PN, 'color', 'k','linewidth', LW) % upsampled PN activity
    plot(temp, TRN, 'color', 'r','linewidth', LW) % upsampled  TRN activity
    % PN with TRN activity
    plot(temp, PN_TRN, 'color',Color('dodgerblue'),'linewidth', LW) % upsampled
% Labels
xlabel('Temperature (\circC)')
ylabel('Relative activity (a.u.)')
formatFig(fig, false);

save_figure(fig,[figDir 'Strong TRN inhibition'],'-png',false,false);



save_figure(fig,[figDir 'Range of TRN inhibition'],'-png',false,false);




% Plot all lines figure:
LW = 0.5;
n = 30;
gains = linspace(0.1,1,n); % test a range of TRN activity levels
CList = Color('Cyan', 'Darkblue', n);
ZList = Color('yellow', 'darkred',n);
Z = [];
fig = figure; hold on
    % TRN activity
    for i = 1:n
        Z(:,i) = TRN*gains(i);
        plot(temp, Z(:,i), 'color', ZList(i,:),'linewidth', LW) % upsampled PN activity
    end
    % new PN activity
    for i = 1:n
        Y = PN' - Z(:,i);
        Y(Y<0) = 0;% rectify
        plot(temp, Y, 'color',CList(i,:),'linewidth', LW) % upsampled
    end
    plot(temp, PN, 'color', 'k','linewidth', 2) % upsampled PN activity

% Labels
xlabel('Temperature (\circC)')
ylabel('Relative activity (a.u.)')
formatFig(fig, false);
save_figure(fig,[figDir 'Range of TRN inhibition'],'-png',false,false);
