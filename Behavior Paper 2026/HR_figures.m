
saveDir = createFolder([figDir, 'Paper Figures/']);
initial_var = add_var(initial_var, 'saveDir');

%% FIGURE: normalized 0 to max temp tuning curves for the four essential behaviors
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

paramList = {'CI', 'sleep', 'foodQuad', 'OutterRing', 'jump'};
colors = getPaperColors(paramList); % colors determined by universal color scheme

ntemps = length(data.tempbin.temps);
nParams = length(paramList);
TC = struct;
for i = 1:nParams
    [TC.(paramList{i}).avg, TC.(paramList{i}).std] = deal(nan(ntemps,2)); % cooling then warming for columns
end

% roi = 1:640000; % eliminate the last bit of the experiment since they fall off the food...for the LTS trials

% preformat the z-score for the data that will be extracted and used: 
Zdata = [];
Zdata(:,1) = (sum(data.CI,2)./num.trials); % courtship index
for i = 2:nParams
    Zdata(:,i) = sum(sum(data.(paramList{i}),2),3)./(num.trials*2);
end
Zdata(isnan(Zdata)) = 0;
% for i = 1:nParams
%     Zdata(:,i) = smooth(Zdata(:,i),200,'moving');
%     Zdata(:,i) = rescale(Zdata(:,i));
% end
% mVal = max(Zdata);
% Zdata = Zdata./mVal; % normalize on a [0-1] scale

% % dirty time course plot of the zscores
% figure;
% hold on
% for i = 1:4
%     % y = smooth(Zdata(:,i),5000,'moving');
%     y = Zdata(:,i);
%     plot(y(1:100:end))
% end

for t = 1:ntemps
    for type = 1:2 % cooling then warming
        switch type 
            case 1 
                ROI = data.tempbin.cooling(:,t);
            case 2
                ROI = data.tempbin.warming(:,t);
        end
        
        for i = 1:nParams % courtship, sleep, food, escape
            y = Zdata(ROI,i);
            y_avg = mean(y,'omitnan');
            TC.(paramList{i}).avg(t,type) = y_avg;
            if ~isempty(y)
                y_err = std(y, 0,1,'omitnan');
                TC.(paramList{i}).std(t,type) = y_err;
            end
        end
    end
end

% Plot
plotErr = false;
sSpan = 8;
r = 1;
c = 2; 

x = data.tempbin.temps';

% fig = getfig('', 1,[774 680]);
% for param = 1: nParams
%     raw_y = TC.(paramList{param}).avg;
%     for i = 1:size(raw_y,2)
%         raw_y(:,i) = smooth(raw_y(:,i), sSpan, 'moving');
%     end
%     scaleY = rescale(raw_y);
%     for type = 1:2
%         subplot(r,c,type); hold on
%         kolor = Color(colors{param});
%         y = scaleY(:,type);
%         plot(x, y, 'color', kolor, 'linewidth',2)
%     end
% end

fig = getfig('', 0,[ 868 806]]);
for type = 1:2
    subplot(r,c,type); hold on
    for param = 1: nParams
        raw_y = smooth(TC.(paramList{param}).avg(:,type), sSpan, 'moving');
        scaleY = rescale(raw_y);
        kolor = colors(param, :);
        plot(x, scaleY, 'color', kolor, 'linewidth',2)
    end
end

% --------- Formatting -------------
formatFig(fig, blkbgd, [r,c]);
matchAxis(fig, true);
subplot(r,c,2)
set(gca, 'ycolor', 'none')
xlabel('temperature (\circC)')
title('warming','color', foreColor,'fontname', 'Arial','FontAngle','italic')
subplot(r,c,1)
ylabel('flies (norm %)')
set(gca, 'xdir', 'reverse')
xlabel('temperature (\circC)')
title('cooling','color', foreColor,'fontname','Arial','FontAngle','italic')
if strcmp(groupName, 'Berlin LTS caviar')
    for i = 1:2
        subplot(r,c,i)
        set(gca, 'xtick', 15:5:35)
        xlim([13, 37])
        temp_lims = [14.5, 35.5];
    end
end

% update y axis to fit temp region indicator
ylims = [-0.05, 1.05];
for ii = 1:2
    subplot(r,c,ii)
    ylim(ylims)
end

% add time arrows 
for ii = 1:2
    subplot(r,c,ii)
    addTimeArrow(gca, foreColor)
end

% plot a line for the 'safe' vs 'threat' zones: 
y = rangeLine(fig, 0, false);
LW = 5;
x_less = [temp_lims(1), 25];
x_more = [25, temp_lims(2)];
subplot(r,c,2)
    plot(x_more, [y,y], 'Color', Color('darkred'), 'LineWidth',LW, 'HandleVisibility','off') % threat
    plot(x_less, [y,y], 'Color', Color('grey'), 'LineWidth',LW, 'HandleVisibility','off') % safe
subplot(r,c,1)
    plot(x_less, [y,y], 'Color', Color('darkred'), 'LineWidth',LW, 'HandleVisibility','off') % threat
    plot(x_more, [y,y], 'Color', Color('grey'), 'LineWidth',LW, 'HandleVisibility','off') % safe


% add arrows for the time at half-max? (TODO 5.1)

save_figure(fig, [saveDir 'normalized behavior for timing temp tuning curve'],fig_type);

%% FIGURE: scatter plot of normalized half-max timepoint for each behavior

