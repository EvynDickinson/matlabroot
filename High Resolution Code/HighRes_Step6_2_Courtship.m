

%% FIGURE: Quantity of courtship during each thermal threat period by fly
% manually adjust for other temperature protocols outside the LTS
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% build indexes for warm threat, cool threat, warm safe, cool safe
data.tempbin.WT = data.tempbin.h_idx & data.temp>25; % warm threat
data.tempbin.WS = data.tempbin.c_idx & data.temp>25; % warm safe
data.tempbin.CT = data.tempbin.c_idx & data.temp<25; % cool threat
data.tempbin.CS = data.tempbin.h_idx & data.temp<25; % cool safe

% pull the CI data
[y_avg, y_sem] = deal(nan([4,1]));
y = nan([num.trials, 4]);

types = {'WT', 'WS', 'CT', 'CS'};
tb = data.tempbin;
for i = 1:4 % for each of the temp regime types
    roi = tb.(types{i});
    y_raw = data.CI(roi,:); % raw data
    y(:,i) = mean(y_raw,'omitnan');
    y_avg(i) = mean(y(:,i),'omitnan');
    y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
end

% plot the CI data
% Create a bar plot for the average CI data across temperature regimes
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
sz = 75;
figure;
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:4, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:4,[num.trials, 1]);
    scatter(x(:), y(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels);
    ylabel('Courtship Index (CI)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
save_figure(gcf, [figDir 'CI per temp regime'], fig_type);
    
% run statistical comparisions between safe and threatening for hot/cold
% respectively 

% Warm threat vs warm safe
[h,p,~,stats] = ttest(y(:,1), y(:,2)); 
fprintf('CI warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

% cool threat vs cool safe
[h,p,~,stats] = ttest(y(:,3), y(:,4)); 
fprintf('CI cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);


% Comparison binned across all thermal threats vs safe temps (regardless of
% warm vs cold absolute temp)

% pull the CI data
[y_avg, y_sem] = deal(nan([2,1]));
y = nan([num.trials, 2]);

tb = data.tempbin;
for i = 1:2 % for each of the temp regime types
    switch i
        case 1
            roi = tb.(types{1}) | tb.(types{3}); % threat
        case 2
            roi = tb.(types{2}) | tb.(types{4}); % safe
    end
    y_raw = data.CI(roi,:); % raw data
    y(:,i) = mean(y_raw,'omitnan');
    y_avg(i) = mean(y(:,i),'omitnan');
    y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
end

% plot the CI data
% Create a bar plot for the average CI data across temperature regimes
labels = {'Threat', 'Safe'};
sz = 75;
figure;
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:2, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:2,[num.trials, 1]);
    scatter(x(:), y(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels);
    ylabel('Courtship Index (CI)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
save_figure(gcf, [figDir 'CI safe vs threat'], fig_type);
    
% run statistical comparisions between safe and threatening for hot/cold
% respectively 

% Warm threat vs warm safe
[h,p,~,stats] = ttest(y(:,1), y(:,2)); 
fprintf('CI safe vs threat: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

%% FIGURE: Quantity of all Courtship Measures during each thermal threat period by fly
% manually adjust for other temperature protocols outside the LTS
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% build indexes for warm threat, cool threat, warm safe, cool safe
data.tempbin.WT = data.tempbin.h_idx & data.temp>25; % warm threat
data.tempbin.WS = data.tempbin.c_idx & data.temp>25; % warm safe
data.tempbin.CT = data.tempbin.c_idx & data.temp<25; % cool threat
data.tempbin.CS = data.tempbin.h_idx & data.temp<25; % cool safe

% pull the CI data
[y_avg, y_sem] = deal(nan([4,1]));
y = nan([num.trials, 4]);

types = {'WT', 'WS', 'CT', 'CS'};
tb = data.tempbin;
for i = 1:4 % for each of the temp regime types
    roi = tb.(types{i});
    data.wing_ext_all(isnan(data.wing_ext_all)) = false;
    data.chase_all(isnan(data.chase_all)) = false;
    data.circling_all(isnan(data.circling_all)) = false;
    y_raw = data.wing_ext_all | data.chase_all | data.circling_all;
    y_raw = y_raw(roi,:); % raw data for the selected regime
    y(:,i) = mean(y_raw,'omitnan');
    y_avg(i) = mean(y(:,i),'omitnan');
    y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
end

% plot the CI data
% Create a bar plot for the average CI data across temperature regimes
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
sz = 75;
figure;
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:4, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:4,[num.trials, 1]);
    scatter(x(:), y(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels);
    ylabel('Courtship Index All');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
save_figure(gcf, [figDir 'CI_all per temp regime'], fig_type);
    
% run statistical comparisions between safe and threatening for hot/cold
% respectively 

% Warm threat vs warm safe
[h,p,~,stats] = ttest(y(:,1), y(:,2)); 
fprintf('CI warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

% cool threat vs cool safe
[h,p,~,stats] = ttest(y(:,3), y(:,4)); 
fprintf('CI cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);


% Comparison binned across all thermal threats vs safe temps (regardless of
% warm vs cold absolute temp)

% pull the CI data
[y_avg, y_sem] = deal(nan([2,1]));
y = nan([num.trials, 2]);

tb = data.tempbin;
for i = 1:2 % for each of the temp regime types
    switch i
        case 1
            roi = tb.(types{1}) | tb.(types{3}); % threat
        case 2
            roi = tb.(types{2}) | tb.(types{4}); % safe
    end
    data.wing_ext_all(isnan(data.wing_ext_all)) = false;
    data.chase_all(isnan(data.chase_all)) = false;
    data.circling_all(isnan(data.circling_all)) = false;
    y_raw = data.wing_ext_all | data.chase_all | data.circling_all;
    y_raw = y_raw(roi,:); % raw data for the selected regime
    y(:,i) = mean(y_raw,'omitnan');
    y_avg(i) = mean(y(:,i),'omitnan');
    y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
end

% plot the CI data
% Create a bar plot for the average CI data across temperature regimes
labels = {'Threat', 'Safe'};
sz = 75;
figure;
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:2, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:2,[num.trials, 1]);
    scatter(x(:), y(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels);
    ylabel('Courtship Index All');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
save_figure(gcf, [figDir 'CI safe vs threat'], fig_type);
    
% run statistical comparisions between safe and threatening for hot/cold
% respectively 

% Warm threat vs warm safe
[h,p,~,stats] = ttest(y(:,1), y(:,2)); 
fprintf('CI_all safe vs threat: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);


%% FIGURE: Distribution of courtship behaviors by temperature regime
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% pull the CI data
[y_avg, y_sem] = deal(nan([4,1]));
y = nan([num.trials, 4]);

types = {'WT', 'WS', 'CT', 'CS'};
tb = data.tempbin;
for i = 1:4 % for each of the temp regime types
    roi = tb.(types{i});
    y_raw = replaceNaN(data.wing_ext_all,false) | ...
        replaceNaN(data.chase_all,false) | replaceNaN(data.circling_all,false);
    y_raw = y_raw(roi,:); % raw data for the selected regime
    y(:,i) = mean(y_raw,'omitnan');
    y_avg(i) = mean(y(:,i),'omitnan');
    y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
end

% plot the CI data
% Create a bar plot for the average CI data across temperature regimes
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
sz = 75;
figure;
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:4, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:4,[num.trials, 1]);
    scatter(x(:), y(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels);
    ylabel('Courtship Index All');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
save_figure(gcf, [figDir 'CI_all per temp regime'], fig_type);

%% FIGURE: Distribution of each courtship behavior across temperature regimes
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% pull the CI data
[y_avg, y_sem] = deal(nan([4,1]));
[yWi, yCh, yCi] = deal(nan([num.trials, 4]));

types = {'WT', 'WS', 'CT', 'CS'};
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
colors = {'teal', 'purple', 'gold'};
behaviors = {'wing_ext_all','chase_all', 'circling_all'};
tb = data.tempbin;
x = 1:length(types);

fig = getfig('',1,[504 620]);
    hold on
    for i = 1:length(behaviors)
        kolor = Color(colors{i});
        y = nan([num.trials, 4]);
        % extract the behavior data for each temp condition to plot
        for t = 1:4 % for each of the temp regime types
            roi = tb.(types{t});
            y(:,t) = mean(data.(behaviors{i})(roi,:),'omitnan').*100; % convert to percent of time
            y_avg(t) = mean(y(:,t),'omitnan');
            y_sem(t) = std(y(:,t),'omitnan')/sqrt(num.trials);
        end
        % plot the data: 
        fill_data = error_fill(x, y_avg, y_sem);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
        set(h, 'facealpha', 0.2)
        plot(y_avg, 'o-', 'Color', kolor, 'MarkerFaceColor', kolor,...
            'LineWidth', 1.5,'linestyle', '--');
    end
% formatting and labels
formatFig(fig, blkbgd);
set(gca, 'xtick', x, 'XTickLabel', labels)
ylabel('avg time in behavior (%)')
legend(strrep(behaviors,'_',' '),'box', 'off')
xlim([0.5,length(types)+0.5])

save_figure(gcf, [figDir 'Courtship behaviors across temp regions'], fig_type);


%% Courtship Attempts per Encounters
% define when a male came within courthship proximity to the female
% on how many of those encounters did he try to court?
% encounters: male fly within 3 body lengths of the female fly
clearvars('-except',initial_var{:})

IFD = 7; % mm distance that counts as an encounter

% frame locations of any courtship behavior:
courtship_locs = replaceNaN(data.wing_ext_all,false) | ...
        replaceNaN(data.chase_all,false) | replaceNaN(data.circling_all,false);

% initialize parameters for saving the data
encounters = struct;

% extract the data for each trial
for trial = 1:num.trials
    % find list of locations the male was within courtship distance
    possible_locs = fly(trial).T.IFD <= IFD; %(mm (flies within 1 BL = 2.5mm so this is ~ <3BLs)

    % create a list of total encounters 
    en_ON_loc = find(diff(possible_locs)==1); % start of encounter
    en_OFF_loc = find(diff(possible_locs)==-1); % end of encounter
    % check if we start with a possible encounter: 
    if possible_locs(1)
        en_ON_loc = [true; en_ON_loc]; % add start at frame 1 if they are within distance
    end
    if possible_locs(end) % add end location if they are still within range in last frame
        en_OFF_loc = [en_OFF_loc; true];
    end
    en_locs = [en_ON_loc,en_OFF_loc(1:length(en_ON_loc))]; 

    total_encounters = size(en_locs,1); % how many encounter periods did the flies have
    
    % find the number of the encounters that included a courtship attempt:
    pos_en = false([total_encounters,1]);
    for i = 1:total_encounters
        roi = (en_locs(i,1):en_locs(i,2)); % frames for this encounter
        pos_en(i) = any(courtship_locs(roi));
    end
    en_locs = [en_locs,pos_en]; % add y/n encounter success to location index

    % save the data into the encounters structure for further analysis
    encounters(trial).locs = en_locs;
    encounters(trial).tot_en = total_encounters;
    encounters(trial).succ_en = sum(pos_en);
    encounters(trial).succ_per = (sum(pos_en)/total_encounters)*100;
    encounters(trial).possible_locs = possible_locs;

end

% quick figure of the total courtship per possible encounter across the
% full experiment, regardless of temperature

foreColor = formattingColors(blkbgd);
fig = getfig('',1,[230 620]); hold on
    y = [encounters(:).succ_per];
    x = ones([1, num.trials]); %1:num.trials;
    scatter(x,y,75,foreColor, 'filled','xjitter', 'density','MarkerFaceAlpha',0.75)
    y_avg = median(y, 'omitnan');
    plot([0.5, 1.5],[y_avg, y_avg], 'color', Color('grey'), 'LineWidth',1.5)
    % TODO: add a grey mean line
    xlabel('trial')
    ylabel('courtship attempts of encounters (%)')
    formatFig(fig, blkbgd);
    xlim([0.25, 1.75])
    set(gca, 'xcolor', 'none')
    y_avg = num2str(round(y_avg,1));
    y_max = min([length(y_avg),3]);
    title([num2str(IFD) ' mm | ' y_avg(1:y_max) '%'],'color', foreColor)

save_figure(gcf, [figDir 'Courtship attempts of encounters total dist ' num2str(IFD)], fig_type);

% How does the courtship-per-encounter rate change across different
% temperature regimes?


types = {'WT', 'WS', 'CT', 'CS'};
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};

tb = data.tempbin;
[tot_encounters,tot_courtship,court_rate] = deal(nan([num.trials,4]));

for i = 1:4 % for each of the temp regime types
    roi = find(tb.(types{i}));
    % which encounter locations are within the temp regime
    for trial = 1:num.trials
        % encounters within this temperature regime
        encounter_locs = ismember(encounters(trial).locs(:,1),roi);
        tot_encounters(trial,i) = sum(encounter_locs);
        tot_courtship(trial,i) = sum(encounters(trial).locs(encounter_locs,3));
        court_rate(trial,i) = (tot_courtship(trial,i)/tot_encounters(trial,i))*100;
    end
end

y_avg = median(court_rate,1,'omitnan');
y_sem = std(court_rate,0,1,'omitnan')./sqrt(num.trials);
% plot the encounter rate data (bar graph)
% Create a bar plot for the average CI data across temperature regimes
sz = 75;
fig = getfig('', 1, [436 620]);
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:4, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:4,[num.trials, 1]);
    scatter(x(:), court_rate(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels,'XTickLabelRotation',30);
    ylabel('Courtship per encounters (%)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
    ylim([-3,70])
save_figure(gcf, [figDir 'courtship per encounter across temp regimes'], fig_type);

% run simple stats:
% Warm threat vs warm safe
[h,p,~,stats] = ttest(court_rate(:,1), court_rate(:,2)); 
fprintf('court rate warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

% cool threat vs cool safe
[h,p,~,stats] = ttest(court_rate(:,3), court_rate(:,4)); 
fprintf('court rate cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);







%% TODO: what is the correlation between male chase and male vs female speed? 





%% TODO: What are the locations of the different courtship behaviors? 




%% TODO: 