

% NULL DISTRIBUTION ANALYSIS /  PROCESSING

%% Near inifinite place null distribution
clearvars('-except',initial_vars{:})

% NULL DISTRIBUTION (this should be identical for all trials/experiments in the same arena)
% Generate the null distribution of distances to the food well:
pixel_buffer = 50; % edge of arena pixel buffer
i = 1;
trial = 1;

% Set axis limits for the selected arena
c1 = data(i).data(trial).data.centre(1);
c2 = data(i).data(trial).data.centre(2);
r = data(i).data(trial).data.r;
xlimit = [c1-(r+pixel_buffer),c1+(r+pixel_buffer)];
ylimit = [c2-(r+pixel_buffer),c2+pixel_buffer+r];
foodWellLoc = data(i).data(trial).data.wellcenters(:,data(i).T.foodLoc(trial));

binSpacing = logspace(1,4,4);

nulls = struct;

for tt = 1: length(binSpacing)
    nbins = binSpacing(tt);
    % find the 'auto bin' lines
    xedge = linspace(xlimit(1),xlimit(2),nbins+1);
    yedge = linspace(ylimit(1),ylimit(2),nbins+1);
    X = []; Y = [];
    idx = 1;
    for y_loc = 1:nbins
        for x_loc = 1:nbins
            X(idx) = xedge(x_loc);
            Y(idx) = yedge(y_loc);
            idx = idx + 1;
        end
    end
    % screen out the units outside the arena circle
    temp_dist = sqrt((X-c1).^2 + (Y-c2).^2);
    loc = temp_dist>r;
    X(loc) = [];
    Y(loc) = [];
    
    % find food well  distance
    nulls(tt).dist = sqrt((X-foodWellLoc(1)).^2 + (Y-foodWellLoc(2)).^2)./pix2mm;
    disp(tt)
end



%% FIGURE: Cumulative distribution of sleeping distances to food
% clearvars('-except',initial_vars{:})
% [foreColor,backColor] = formattingColors(blkbgd);
LW = 2;
[~,backColor] = formattingColors(blkbgd);

% PLOT NULL DISTRIBUTION:
fig_dir = [baseFolder 'Fundamentals\'];

colorList = Color('White','teal', length(binSpacing));

% FIGURE:
fig = getfig('',1,[687 680]);  
hold on
for tt = 1: length(binSpacing)
    % NULL DISTRIBUTION:
    null_CDF = cdfplot(nulls(tt).dist);
    null_CDF.Color = colorList(tt,:);
    null_CDF.LineWidth = 2; 
end
formatFig(fig,blkbgd);
xlabel('distance to well (mm)','FontSize',20)
ylabel('Empirical Cumulative Distribution','FontSize',20)
set(gca,'TickDir','out')
set(gca,'GridColor',backColor)

% save figure
save_figure(fig,[fig_dir  'CDF of null distribution'],fig_type);

binWidth = 1; %mm spacing for bins
binedges = 0:binWidth:55; 

% find probability of occupancy for each bin from the null distribution
fig = getfig('',1,[687 680]); hold on
for tt = 1: length(binSpacing)

    N = discretize(nulls(tt).dist,binedges);
    N_tot = length(nulls(tt).dist);
    prob = zeros(length(binedges),1);
    for i = 1:length(binedges)
        prob(i) = (sum(N==i))/N_tot;
    end

    plot(binedges,prob,'color', colorList(tt,:),'linewidth',LW,'linestyle','-')

end
% plot(binedges,smooth(prob,sSpan,'moving'),'color', Color('cyan'),'linewidth',LW+2,'linestyle','-')
formatFig(fig,blkbgd);
xlabel('distance to well (mm)','FontSize',20)
ylabel('occupancy probability','FontSize',20)
set(gca,'TickDir','out')


% Save null distribution probability
save_figure(fig,[fig_dir  'Distance probabilty null distribution'],fig_type);


% find probability of occupancy for each bin from the null distribution
fig = getfig('',1);
for tt = 1: length(binSpacing)
subplot(2,2,tt)
    N = discretize(nulls(tt).dist,binedges);
    N_tot = length(nulls(tt).dist);
    prob = zeros(length(binedges),1);
    for i = 1:length(binedges)
        prob(i) = (sum(N==i))/N_tot;
    end
    plot(binedges,prob,'color', colorList(tt,:),'linewidth',LW,'linestyle','-')
    xlabel('mm')
    ylabel('occ prob')
    title(num2str(binSpacing(tt)))
end
formatFig(fig,blkbgd,[2,2]);

% Save null distribution probability
save_figure(fig,[fig_dir  'Separated distance probabilty null distribution'],fig_type);



% find probability of occupancy for each bin from the null distribution
fig = getfig('',1);
for tt = 1: length(binSpacing)
subplot(2,2,tt)
    yyaxis left
    histogram(nulls(tt).dist,binedges,'FaceColor',foreColor,'FaceAlpha',0.8)
    hold on
    yyaxis right
    N = discretize(nulls(tt).dist,binedges);
    N_tot = length(nulls(tt).dist);
    prob = zeros(length(binedges),1);
    for i = 1:length(binedges)
        prob(i) = (sum(N==i))/N_tot;
    end
    plot(binedges,prob,'color', colorList(tt,:),'linewidth',LW,'linestyle','-')
end
formatFig(fig,blkbgd,[2,2]);
for tt = 1: length(binSpacing)
    subplot(2,2,tt)
    yyaxis left
    set(gca,'YColor',backColor)
end

% Save null distribution probability
save_figure(fig,[fig_dir  'Separated null distribution'],fig_type);





