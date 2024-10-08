

%% Antoine Equation: gives Vapor Pressure value for each odor at specified temps
% https://webbook.nist.gov/chemistry/name-ser/
clear

ko = 272.15; % kelvin offset
t = 15:35; % temperature range to plot
T = t+ko; % temperature in Kelvin

odors = {'2-butanone', 'ethyl acetate', 'acetic acid','methyl acetate',...
                'isopentyl acetate','benzaldehyde','3-octanol','ethyl butyrate','proprionic acid'};
coeffc = [3.9894	1150.207	-63.904;...
            4.22809	1245.702	-55.189;...
            4.68206	1642.54	-39.764;...
            4.20364	1164.426	-52.69;...
            5.08047	1932.043	-28.698;...
            5.21496	2337.539	-5.103;...
            4.8465	1663.322	-97.47;...
            4.33187	1509.443	-45.284;...
            4.74558	1679.869	-59.832] ;

nOdors = length(odors);
CList = jet(nOdors);

fold_change = [];
VP = [];
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
    % save the vapor pressure values
    VP(:,i) = P;
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

save_figure(fig,[bFolder 'Electrophysiology Modeling/Temperature pressure curves'],'-png',true);

%% Use Henry's Law for Ideal Dilute Solutions to determine the change in odor concentration across temps
% assuming ideal conditions, so solubitily = 1 across this range of temps
C = 1;
dilutions = 10e-1; % dilution of solution

pVap = VP*0.986923; % convert pvap from bar to atm
% 1) Determine Henry's Constant for each temp
H = (pVap./C);
% 2) Predict headspace concentration
Conc = H .* dilutions;

% determine the range of concentration change
[delta_C, dil_c] = deal([]);
idx = 0;
for ii = 2 : 6
    idx = idx + 1;
    d = 10^-ii; % dilution of solution
    Cm = H .* d;
    delta_C(idx,:) = 100 .* ((Cm(end,:) -  Cm(1,:)) ./ Cm(end,:));
    dil_c(idx,:) = d * ones(1,nOdors);
end

% plot the solution concentration vs prediction odor headspace concentration for one
% odor for all the temperatures
nT = floor(length(t)/2);
CList = [Color('dodgerblue', 'grey', nT); Color('grey', 'red', length(t)-nT)];
bFolder = getCloudPath;
for odor =  1:nOdors

    fig = figure;
    hold on
        % odor = 1; % odor ID
        % determine the range of concentration change
        [X, Y] = deal([]);
        for ii = 6:-1:2 % start dilutions
            d = 10^-ii; % dilution of solution
            Cm = H(:,odor) .* d;
            Y(:,ii) = Cm;
            X(:,ii) = d*ones(length(t),1);
        end
        % plot each concentration tuning curve for a given temp
        for i = 1:length(t)
            plot(X(i,:), Y(i,:), 'color', CList(i,:), 'linewidth', 1, 'marker', 'o')
        end
        set(gca, 'xscale', 'log')
        set(gca, 'yscale', 'log')
        title(odors{odor})
        ylabel('[Odor headspace]')
        xlabel('[Odor solution]')
        formatFig(fig, false);
        % legend('location', 'northwest')
        colormap(CList);
        h = colorbar;
        % h.Label.String = 'temperature (\circC)';
        h.Ticks = [0,1];
        h.TickLabels = {[num2str(t(1)) '(\circC)'], [num2str(t(end)) '(\circC)']};

        save_figure(fig,[bFolder 'Electrophysiology Modeling/' odors{odor} ' concentration curves'],'-png',true);

end















