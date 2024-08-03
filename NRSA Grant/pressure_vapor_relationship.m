

%% Antoine Equation
% https://webbook.nist.gov/chemistry/name-ser/
clear

ko = 272.15; % kelvin offset
t = 15:35; % temperature range to plot
T = t+ko; % temperature in Kelvin

odors = {'2-butanone', 'ethyl acetate', 'acetic acid','3-octanol','ethyl butyrate','methyl acetate', 'benzaldehyde'};
coeffc = [3.9894	1150.207	-63.904;...
            4.22809	1245.702	-55.189;...
            4.68206	1642.54	-39.764;...
            4.8465	1663.322	-97.47;...
            4.33187	1509.443	-45.284;...
            4.20364	1164.426	-52.69;...
            5.21496	2337.539	-5.103] ;
nOdors = length(odors);
CList = colormap(nOdors);
cmap(nOdors)

fold_change = [];
c = 3;
r = 1;
sb(1).idx = 1:2;
sb(2).idx = 3;

fig = getfig('',1);
subplot(r,c,sb(1).idx)
hold on
for i = 1:nOdors
    A = coeffc(i, 1);
    B = coeffc(i,2);
    C = coeffc(i,3);

    P = 10.^(A-(B./(T + C))); %pressure in BAR
    fold_change(i) = (P(end)-P(1))/P(1);
    plot(t, P, 'DisplayName', odors{i},'linewidth', 1.5,'color', CList(i,:));
end
xlabel('Temperature (\circC)');
ylabel('Vapor Pressure (bar)');
hold off;
set(gca, 'tickdir', 'out')
legend show;
legend('box', 'off','location', 'northwest','FontSize',10)
 
subplot(r,c,sb(2).idx)
scatter(1:nOdors, fold_change,100,CList,'filled')
xlim([0,nOdors+1])
ax = gca;
set(ax, 'xtick', 1:nOdors,'xticklabel', odors)
ylabel('fold change from 15\circC to 35\circC')
formatFig(fig, false,[r,c],sb);
