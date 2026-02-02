

%% FIGURE: Quantity any type of courtship during each temperature region by fly
% manually adjust for other temperature protocols outside the LTS
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

all_types =  {'WT', 'WS', 'CT', 'CS', 'SS'};
all_labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe','Static Safe'};

% Select the type of courtship data to plot:
court_options = {'CI', 'CI_all', 'wing_ext', 'wing_ext_all', 'court_chase', 'chase_all', 'circling_1sec', 'circling_all'};
court_string = {'Courtship Index', 'Unrestricted CI', 'Wing Extension', 'All Wing Extension', 'Chase', 'All Chase', 'Circling', 'All Circling'};
idx = listdlg("PromptString","Select the type of data to plot","SelectionMode","single",...
    "ListSize",[200, 300], 'ListString',court_options);
param = court_options{idx};
param_str = court_string{idx};

% ADJUST HERE TO WORK FOR DIFFERENT DYNAMIC PROTOCOLS (1/16/2026) | make a
% universal temp protocol loading function
if contains(groupName, 'F LRR 25-17')
    type_sel = 3:5;
    temp_type = 1;
elseif contains(groupName, 'LTS')
    type_sel = 1:4;
    temp_type = 2;
end

types = all_types(type_sel);
labels = all_labels(type_sel);

% pull the CI data
nTypes = length(types);
[y_avg, y_sem] = deal(nan([nTypes,1]));
y = nan([num.trials, nTypes]);

tb = data.tempbin;
for i = 1:nTypes % for each of the temp regime types
    roi = tb.(types{i});
    y_raw = data.(param)(roi,:); % raw data
    y(:,i) = mean(y_raw,'omitnan');
    y_avg(i) = mean(y(:,i),'omitnan');
    y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
end

% plot the CI data
% Create a bar plot for the average CI data across temperature regimes
sz = 75;
fig = getfig('All threat vs safe', 1, [904 579]);
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:nTypes, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:nTypes,[num.trials, 1]);
    scatter(x(:), y(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels);
    ylabel(param_str);
    title_str = [param_str ' per temp regions'];
    % title(title_str);
    formatFig(gcf, blkbgd);
save_figure(fig, [figDir title_str], fig_type);

% run STATISTICAL comparisions between safe and threatening for hot/cold
switch temp_type
    case 1 % F LRR 25-17
        % {'CT', 'CS', 'SS'} comparisons between these
        n = 3; % number of comparisons
        alpha = 0.05/n;
        % Cold threat vs Cold Safe 
        [~,p,~,stats] = ttest(y(:,1), y(:,2)); 
        h = p<=alpha; % Bonferonni MCC for two comparisons
        fprintf('cold threat vs safety: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
        % Cold threat vs Static Safe
        [~,p,~,stats] = ttest(y(:,1), y(:,3)); 
        h = p<=alpha; % Bonferonni MCC
        fprintf('Cold threat vs Static Safe: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
        % Cold threat vs Static Safe
        [~,p,~,stats] = ttest(y(:,2), y(:,3)); 
        h = p<=alpha; % Bonferonni MCC
        fprintf('Cold Safe vs Static Safe: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

    case 2 % LTS 15-35
        % {'WT' VS 'WS' | 'CT' vs 'CS'} comparisons between these
        % Warm threat vs safe & cold threat vs safe
        n = 2;
        alpha = 0.05/n; % Bonferonni lowering of significance level for a comparison
        [~,p,~,stats] = ttest(y(:,1), y(:,2));
        h = p<=alpha; % Bonferonni MCC
        fprintf('warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
        [~,p,~,stats] = ttest(y(:,3), y(:,4)); 
        h = p<=alpha; % Bonferonni MCC
        fprintf('cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
end

% Comparison binned across all thermal threats vs safe temps (regardless of
% warm vs cold absolute temp)
if temp_type==2 % LTS
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
        y_raw = data.(param)(roi,:); % raw data
        y(:,i) = mean(y_raw,'omitnan');
        y_avg(i) = mean(y(:,i),'omitnan');
        y_sem(i) = std(y(:,i),'omitnan')/sqrt(num.trials);
    end
    
    % plot the CI data
    % Create a bar plot for the average CI data across temperature regimes
    labels = {'Threat', 'Safe'};
    sz = 75;
    fig = getfig('All threat vs safe', 1, [469 579]);
        bar(y_avg, 'FaceColor', foreColor);
        hold on;
        errorbar(1:2, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
        % scatter plot points
        x = repmat(1:2,[num.trials, 1]);
        scatter(x(:), y(:), sz, Color('grey'), 'filled', 'xjitter','density','MarkerFaceAlpha', 0.75);
        set(gca, 'XTickLabel', labels);
        ylabel(param_str);
        % title('Courtship Behavior During Thermal Threats');
        formatFig(gcf, blkbgd);
    save_figure(gcf, [figDir param_str ' safe vs threat'], fig_type);
        
    % run statistical comparisions between safe and threatening for hot/cold
    
    % Warm threat vs warm safe
    [h,p,~,stats] = ttest(y(:,1), y(:,2)); 
    fprintf('CI safe vs threat: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

end

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


%% ANALYSIS: Courtship Attempts per Encounters
% define when a male came within courthship proximity to the female
% on how many of those encounters did he try to court?
% encounters: male fly within 3 body lengths of the female fly
clearvars('-except',initial_var{:})
initial_var{end+1} = 'encounters';

IFD = 12; % mm distance that counts as an encounter

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
        pos_en(i) = any(data.CI_all(roi));
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

%% FIGURE: Rate of encounters across temperature regions
clearvars('-except',initial_var{:})

% % ------------------------------------------------------------------------
% How does the courtship-per-encounter rate change across different
% temperature regimes?

types = {'WT', 'WS', 'CT', 'CS'};
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
nTypes = length(types);

tb = data.tempbin;
[tot_encounters,tot_courtship,court_rate, encounter_rate] = deal(nan([num.trials,4]));

for i = 1:nTypes % for each of the temp regime types
    roi = find(tb.(types{i})); % find the frames within the region
    time_period = length(roi); % how many frames possible could there be for encounters (to normalize in case there are difference in the time periods of the diff temp regimess
    % disp(time_period)
    % which encounter locations are within the temp regime
    for trial = 1:num.trials
        % encounters within this temperature regime
        encounter_locs = ismember(encounters(trial).locs(:,1),roi);
        tot_encounters(trial,i) = sum(encounter_locs);
        tot_courtship(trial,i) = sum(encounters(trial).locs(encounter_locs,3));
        court_rate(trial,i) = (tot_courtship(trial,i)/tot_encounters(trial,i))*100;
    end
    encounter_rate(:,i) = (tot_encounters(:,i)/time_period).*100;
end

% PLOT THE PERCENTAGE OF FRAMES WITH AN ENCOUNTER PER TEMP REGIME
y_avg = median(encounter_rate,1,'omitnan');
y_sem = std(encounter_rate,0,1,'omitnan')./sqrt(num.trials);
% plot the encounter rate data (bar graph)
% Create a bar plot for the average CI data across temperature regimes
sz = 75;
fig = getfig('', 1, [436 620]);
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:nTypes, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:nTypes,[num.trials, 1]);
    scatter(x(:), encounter_rate(:), sz, Color('grey'), 'filled', ...
        'xjitter','density','xjitterwidth', 0.3, 'MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels,'XTickLabelRotation',30);
    ylabel('Rate of encounters (% of total time)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
    % ylim([-3,70])
save_figure(gcf, [figDir 'rate of encounters across temp regimes'], fig_type);

% run simple stats:
% Warm threat vs warm safe
[h,p,~,stats] = ttest(encounter_rate(:,1), encounter_rate(:,2)); 
fprintf('encounter rate warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% cool threat vs cool safe
[h,p,~,stats] = ttest(encounter_rate(:,3), encounter_rate(:,4)); 
fprintf('encounter rate cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

%%  FIGURE: (todo) PERCENT OF FRAMES WITH AN ENCOUNTER COLLAPSED TO SAFE AND THREAT
clearvars('-except',initial_var{:})

%  ------------------------------------------------------------------------
y_avg = median(encounter_rate,1,'omitnan');
y_sem = std(encounter_rate,0,1,'omitnan')./sqrt(num.trials);
% plot the encounter rate data (bar graph)
% Create a bar plot for the average CI data across temperature regimes
sz = 75;
fig = getfig('', 1, [436 620]);
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:nTypes, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:nTypes,[num.trials, 1]);
    scatter(x(:), encounter_rate(:), sz, Color('grey'), 'filled', ...
        'xjitter','density','xjitterwidth', 0.3, 'MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels,'XTickLabelRotation',30);
    ylabel('Rate of encounters (% of total time)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
    % ylim([-3,70])
save_figure(gcf, [figDir 'rate of encounters across temp regimes'], fig_type);

% ------------------------------------------------------------------------
% PLOT THE COURTSHIP RATE FOR EACH TEMP REGIME
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


%%  FIGURE: (TODO) duration of courtship encounters across each temperature regime
clearvars('-except',initial_var{:})
% ------------------------------------------------------------------------

fps = fly(1).fps;
% encounters(trial).locs = en_locs;
types = {'WT', 'WS', 'CT', 'CS'};
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
nTypes = length(types);

tb = data.tempbin;
[enc_duration] = deal(nan([num.trials,4]));

% TODO: figure out the timing here -- recheck with encounters running only
% for selective courtship metrics to see if the duration increases (it
% should...)

for i = 1:nTypes % for each of the temp regime types
    roi = find(tb.(types{i})); % find the frames within the region
    % find which encounters start within the region and their avg duration
    for trial = 1:num.trials
        % encounters within this temperature regime
        loc = ismember(encounters(trial).locs(:,1),roi) & encounters(trial).locs(:,3)==1;
        % encounter duration in seconds
        enc_duration(trial,i) = mean(encounters(trial).locs(loc,2)-...
                                encounters(trial).locs(loc,1),'omitnan')/fps;
    end
end

% PLOT THE PERCENTAGE OF FRAMES WITH AN ENCOUNTER PER TEMP REGIME
y_avg = median(encounter_rate,1,'omitnan');
y_sem = std(encounter_rate,0,1,'omitnan')./sqrt(num.trials);
% plot the encounter rate data (bar graph)
% Create a bar plot for the average CI data across temperature regimes
sz = 75;
fig = getfig('', 1, [436 620]);
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:nTypes, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:nTypes,[num.trials, 1]);
    scatter(x(:), encounter_rate(:), sz, Color('grey'), 'filled', ...
        'xjitter','density','xjitterwidth', 0.3, 'MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels,'XTickLabelRotation',30);
    ylabel('Rate of encounters (% of total time)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
    % ylim([-3,70])
save_figure(gcf, [figDir 'rate of encounters across temp regimes'], fig_type);

% run simple stats:
% Warm threat vs warm safe
[h,p,~,stats] = ttest(encounter_rate(:,1), encounter_rate(:,2)); 
fprintf('encounter rate warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% cool threat vs cool safe
[h,p,~,stats] = ttest(encounter_rate(:,3), encounter_rate(:,4)); 
fprintf('encounter rate cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);

% ------------------------------------------------------------------------



%% Fly speed for each temperature region
clearvars('-except',initial_var{:})

[foreColor, ~] = formattingColors(blkbgd); % get background colors

types = {'WT', 'WS', 'CT', 'CS'};
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
nTypes = length(types);

% extract the speeds for each fly within the temp region
tb = data.tempbin;
[fly_speed] = deal(nan([num.trials*2,4]));
for i = 1:nTypes % for each of the temp regime types
    roi = find(tb.(types{i})); % find the frames within the region
    all_speed = squeeze(mean(data.speed(roi,:,:),1,'omitnan'));
    fly_speed(:,i) = all_speed(:); % combine across male and female flies here
end
   

% PLOT THE avg speed for each region
y_avg = median(fly_speed,1,'omitnan');
y_sem = std(fly_speed,0,1,'omitnan')./sqrt(num.trials*2);
% plot the encounter rate data (bar graph)
% Create a bar plot for the average CI data across temperature regimes
sz = 75;
fig = getfig('', 1, [436 620]);
    bar(y_avg, 'FaceColor', foreColor);
    hold on;
    errorbar(1:nTypes, y_avg, y_sem, 'color', Color('grey'), 'linestyle', 'none', 'LineWidth', 1.5);
    % scatter plot points
    x = repmat(1:nTypes,[num.trials*2, 1]);
    scatter(x(:), fly_speed(:), sz, Color('grey'), 'filled', ...
        'xjitter','density','xjitterwidth', 0.3, 'MarkerFaceAlpha', 0.75);
    set(gca, 'XTickLabel', labels,'XTickLabelRotation',30);
    ylabel('Avg fly speed (mm/s)');
    % title('Courtship Behavior During Thermal Threats');
    formatFig(gcf, blkbgd);
    % ylim([-3,70])
save_figure(gcf, [figDir 'fly speed per temp regime'], fig_type);

% run simple stats:
% Warm threat vs warm safe
[h,p,~,stats] = ttest(fly_speed(:,1), fly_speed(:,2)); 
fprintf('speed warm temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);
% cool threat vs cool safe
[h,p,~,stats] = ttest(fly_speed(:,3), fly_speed(:,4)); 
fprintf('speed cold temps: h=%d, p=%.4f, t=%.3f, df=%d\n', h, p, stats.tstat, stats.df);


%% COURTSHIP TUNING CURVES

clearvars('-except',initial_var{:})

currFigDir = createFolder([figDir, 'courtship tuning curves/']);

[foreColor, ~] = formattingColors(blkbgd); % get background colors
r = 1; 
c = 2; 
LW = 1;
kolor = Color('dodgerBlue'); % foreColor
plot_err = true;
if strcmp(groupName,'Berlin LTS caviar')
    xlimits = [13, 37];
    xPos = [14.5, 25]; 
end
ylabel_str = '% flies';
typeNames = {'cooling', 'warming'};

fields = {'CI', 'CI_all', 'circling_1sec', 'circling_all', 'court_chase', 'chase_all', 'wing_ext', 'wing_ext_all'};
scaler = 100; % convert fractions to percentages of flies

idx = listdlg("PromptString",'Select data types', 'SelectionMode','multiple', 'ListSize', [200, 300], 'ListString',fields);
nParams = length(idx);
legend_str = cell(1,nParams);
colorList = {'DarkOrchid','Gold','dodgerblue','magenta', 'turquoise','lime','red','Orange'};


fig = getfig('Tuning Curve', 1);
for p = 1:nParams
    kolor = Color(colorList{p});
    % EXTRACT THE DATA
    sel_field = fields{idx(p)};
    legend_str{p} =  sel_field;
    double_field = size(data.(sel_field),2)==2; % if there are two fields for both sexes
    % maxRoi = 730000; % cutoff for food quality loss in the smaller plate
    % find averages for each fly for each temperature bin
    temps = data.tempbin.temps;
    nTemps = length(temps);
    [pData.cooling.avg,pData.cooling.sem, pData.warming.avg,pData.warming.sem]  = deal(nan([nTemps,1]));    
    for type = 1:2 % heating and cooling
        type_name = typeNames{type};
        for t = 1:nTemps
            roi = data.tempbin.(type_name)(:,t);
            % roi(maxRoi:end) = false;
            if double_field
                raw = [squeeze(data.(sel_field)(roi,1,:)), squeeze(data.(sel_field)(roi,2,:))];
            else 
                raw = data.(sel_field)(roi,:);
            end
            processed =  mean(raw,1,'omitnan');
            pData.(type_name).avg(t) = mean(processed,'omitnan');
            pData.(type_name).sem(t) = std(processed,'omitnan')/sqrt(length(processed));
        end
    end

    % PLOT THE DATA
    for type = 1:2 % cooling and warming
        subplot(r, c, type)
        hold on
        x = temps;
        y = pData.(typeNames{type}).avg .* scaler;
        y_err = pData.(typeNames{type}).sem .* scaler;
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.4);
        plot(x,y,'color',kolor,'linewidth',LW+1)
    end
end

% FORMATTING
for type = 1:2
    subplot(r, c, type)
    % formatting
    if type == 1
        set(gca, 'xdir', 'reverse')
        ylabel(ylabel_str)
    else
        set(gca, 'ycolor', 'none')
    end
    xlabel('temp (\circC)')
    xlim(xlimits)
    % ylim(ylimits)
    title(typeNames{type},'color', foreColor)
end
matchAxis(fig, true); % match the y axes between heating and cooling
formatFig(fig, blkbgd,[r,c]);

for type = 1:2
    subplot(r, c, type)
    ylimits = ylim;
    pos = [xPos(1,type), ylimits(1), 10, range(ylimits)]; % [lower-left X, lower-left Y, X-width, Y-height]
    h = rectangle('Position', pos, ...
              'FaceColor', foreColor, ...   % RGB color
              'FaceAlpha', 0.2, ...
              'EdgeColor', 'none');
end
legend(strrep(legend_str,'_',' '), 'textcolor', foreColor, 'box', 'off','fontsize', 12,'location', 'northwest')
            
save_figure(fig, [currFigDir strjoin(legend_str,' ') ' courtship tuning curve'])




%% TODO: what is the correlation between courtship-all & courtship by temperature?
% z-score them?
% or time-series correlation
% or temp by temp bin comparison?


%% TODO: what is the correlation between male chase and male vs female speed? 





%% TODO: What are the locations of the different courtship behaviors? 




%% TODO: 