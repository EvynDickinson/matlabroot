

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









