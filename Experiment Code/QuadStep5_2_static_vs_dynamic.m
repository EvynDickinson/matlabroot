


%% [MANUAL LOADING SECTION]
% load data that already has been saved to plot across static and dynamic
% trials without having to load the massive set of all of them
clear; clc

paths = getPathNames; % get the appropriate file path names
baseFolder = getDataPath(5,0,'Select where you want to find the grouped data structures');
% baseFolder = getCloudPath;
structFolder = [baseFolder paths.group_comparision];

% Pick type of data to load from 4_2 where we export processed data subtypes
param_options = {'ring', 'inner75', 'fullquad', 'innerquad', 'quadring', 'circle10', 'circle7', 'circle5', 'fliesonfood', 'sleep'}; 
idx = listdlg("PromptString",'select data type to load', 'ListString', param_options, 'SelectionMode', 'single');
try param = param_options{idx};
catch
    return
end

% Pick static & dynamic data sets
folderName = selectFolder(structFolder,true,'Select dynamic and static data to load');
nGroups = length(folderName);

% Check if the selected data sets have the appropriate data saved
file_check = false([1,nGroups]);
for i = 1:nGroups
    file_check(i) = isfile([structFolder folderName{i} '/' param ' data.mat']);
end
if any(~file_check)
    {folderName{~file_check}}' %#ok<CCAT1,NOPTS>
    warndlg('Presaved data not found')
    return
end

% Load data for the selected types
data = struct;
for i = 1:nGroups
    dummy = load([structFolder folderName{i} '/' param ' data.mat']);
    data(i).y = dummy.y;
end

% Plot data
% [MANUAL CHECK]
disp(folderName)
d = 1; % dynamic
s = 2; % static
data(s).foodPairs = [1,2; 3,4; 5,6; 7,8; 9,10; 11,12; 13,14; 15,16]; % food vs empty for the static trials
data(d).foodPairs = [1,2]; % food vs empty for the dynamic trials

blkbgd = false; % black background or not

initial_vars = who; 
initial_vars{end+1} = 'initial_vars';
initial_vars{end+1} = 'static_names';
initial_vars{end+1} = 'static_temps';
initial_vars{end+1} = 'LTS_temp';

%% Load fictive temperature protocol: 
% check if it's a hold and if it's dummy LTS then load that data here: 
if ~exist('LTS_temp', 'var')
    drivePath = getCloudPath;   
    drivePath = drivePath(1:end-5);
    load([drivePath, 'LTS 15-35 temp data.mat']); 
end

%% Plot the data 
clearvars('-except',initial_vars{:})

% scatter plot of static vs dynamic sep by warming and cooling for the loaded parameter

tempList = [15 20 25 30 33]; % what temps to plot
nTemps = length(tempList);
% full_temp_list =  [15 17 20 23 25 27 30 33 35];

% Plotting Parameters:
foreColor = formattingColors(blkbgd); %get background colors

food = true; % plot the food data trials
if food
    expIdxD = data(d).foodPairs(1);
    expIdxS = 1;
else 
    expIdxS = 2;
    expIdxD = data(d).foodPairs(2);
end

% extract the temperature information from the static group trial names:
static_names = {data(s).y.name};
static_temps = nan([1,length(static_names)]);
for i = 1:length(static_names)
    temp = str2double(regexp(static_names{i}, '\d+', 'match')); % pull the number from the name
    data(s).y(i).temp = temp;
    static_temps(i) = temp;
end

% Plotting Parameters
buff = 0.3; % x scatter buffer size
SZ = 35; % scatter point size
sideBuff = 0.5; % side shift from temp to fit both comparisons
Lbuff = 2.5; % scaling on distance for avg line
LW = 2;
r = 1;
c = 2;

% PLOT THE FIGURE:
fig = getfig('',1);
%subplots for heating vs cooling
for type = 1:2 % heating then cooling
    subplot(r,c,type); hold on
    switch type
        case 1 % heating
            type_str = 'increasing';
        case 2 % cooling
            type_str = 'decreasing';
    end
    % ____ plot the dynamic trial data first ____ 
    for i = 1:nTemps
        x = tempList(i);
        y = data(d).y(expIdxD).(param);
        % find the temp index for the desired temperatures
        [~, temp_loc] = min(abs(x - y.temps));
        % plot data 
        kolor = Color('teal');
        X = x-sideBuff;
        y_points = y.(type_str).raw(temp_loc,:);
        y_avg = mean(y_points);
        xx = shuffle_data(linspace(X-buff, X+buff,length(y_points)));
        scatter(xx,y_points,SZ, kolor, 'filled','o')
        plot([X-Lbuff*buff, X+Lbuff*buff], [y_avg, y_avg],"Color",kolor,'LineWidth', LW)
    end

    % ____ plot the static trial data  ____ 
    for i = 1:nTemps
        % find the experiment that matches the temperature
        x = tempList(i);
        sIdx = find((abs(x - static_temps))==0);
        sIdx = sIdx(any(data(s).foodPairs(:,expIdxS)==sIdx,1)); % find from the food or empty trials, as it were
    
        % find the fictive temp bin
        [~, temp_loc] = min(abs(x - LTS_temp.temp_list));
        fictive_temp_match = LTS_temp.temp_list(temp_loc);
        if abs(fictive_temp_match - x)>1
            h = warndlg('Large temp discrepency between fictive temp and real temp');
            uiwait(h)
        end
        switch type
            case 1 % warming
                rateIDX = find(LTS_temp.temp_rates>0);
            case 2 % cooling 
                rateIDX = find(LTS_temp.temp_rates<0);
        end
        ROI = LTS_temp.loc(rateIDX,temp_loc).frames;
        % extract the averages for these frames: 
        y = data(s).y(sIdx).(param).all(ROI,:);
        y_points = mean(y,1,'omitnan');
        y_avg = mean(y_points);
        X = x+sideBuff;
        xx = shuffle_data(linspace(X-buff, X+buff,length(y_points)));
    
        % plot data
        kolor = foreColor;
        scatter(xx,y_points,SZ, kolor, 'filled','o')
        plot([X-Lbuff*buff, X+Lbuff*buff], [y_avg, y_avg],"Color",kolor,'LineWidth', LW)
    end
    title(type_str)
end

% format the figure:

% GOTTA LOOK AT THIS WITH THE NO FOOD TRIALS -- DOES THE PRESENCE OF FOOD
% BUFFER THE THREAT OF THE FOOD?
% COULD ALSO BE THE CHANGE IN FOOD VALUE AT HIGHER TEMPS AS IT DRIES OUT
% AND THEREFORE PROVIDES A DIFFERENT EXPERIENCE AND IS LESS INTERESTING TO
% THE FLIES AND THUS DOESN'T KEEP THEM IN THE FOOD REGION... NEED TO FIND A
% WAY TO CONTROL FOR THE FOOD QUALITY CHANGING OVER TIME....













% 
% % autoLim = false;
% % xlim_auto = false; % change the time range for the x axis
% % temp_lims = [12, 35]; % x limit default values
% % occ_lims = [0 100]; % y limit default values
% % food_color = Color('gold');
% % empty_color = foreColor;
% % empty_type = 2; % high occupancy null data
% % LW = 2; % avg line width
% % eLW = 2; % error bar width
% % SZ = 25; % scatter point size 
% % buff = 0.3; % scatter plot buffer size
% % r = 1; % rows
% % c = 2; % columns
% % foreColor = formattingColors(blkbgd);
% 
% % % Stats setup: 
% % SD = []; % stats data structure
% % [SD.w, SD.c] = deal(nan([max(num.trial),num.exp])); % empty mat for the trial average data to fill
% % [SD.w_temps,SD.c_temps] = deal(nan([num.exp,1]));
% 
% % FIGURE:
% fig = getfig('',true);
% ylimits = [];
% for i = 1:num.exp
%     % get temperature for the experiment group: 
%     x = strrep(data(i).temp_protocol,'Hold','');
%     x = str2double(strrep(x,'C',''));
% 
%     % determine which data to grab depending on food or no food: 
%     switch data(i).emptytrial
%         case true  % empty trial 
%             kolor = empty_color;
%             subfield = 'high'; % low or high occupancy quadrant to plot
%         case false % food present
%             kolor = food_color;
%             subfield = 'food';
%     end
% 
%     % pull the subfield data structure from the grouped data set
%     if quad_regions % sub regions (requires '.food' or '.low' extension etc)
%         yy = grouped(i).(pName).(subfield);
%     else
%         yy = grouped(i).(pName); % no subregions in the metric (e.g., ring)
%     end
% 
%     for type = 1:2 % heating and cooling
%         subplot(r,c,type)
%         hold on
%         if FT % fictive temp
%             % for now, use the built in LTS 15-35 but need to adjust this for
%             % future use to be more dynamic...TODO IMPORTANT 7.23
%             [~, temp_loc] = min(abs(x - LTS_temp.temp_list));
%             fictive_temp_match = LTS_temp.temp_list(temp_loc);
%             % disp(['Target temp vs fictive temp: ' num2str(x) ' vs ' num2str(fictive_temp_match)])
%             % throw warning if the temps are really far from each other
%             if abs(fictive_temp_match - x)>1
%                 h = warndlg('Large temp discrepency between fictive temp and real temp');
%                 uiwait(h)
%             end
%             switch type
%                 case 1 % warming
%                     rateIDX = find(LTS_temp.temp_rates>0);
%                 case 2 % cooling 
%                     rateIDX = find(LTS_temp.temp_rates<0);
%             end
%             ROI = LTS_temp.loc(rateIDX,temp_loc).frames;
%             sbpt_str = {'''warming''', '''cooling'''};
%         else % duration! 
%            ROI = timeROI;
%            sbpt_str = {['avg first ' num2str(duration_time) ' hours'], ['avg first ' num2str(duration_time) ' hours']};
%         end
% 
%         % plot data: 
%         y = mean(yy.all(ROI,:),1,'omitnan');
%         y_avg = mean(y);
%         y_err = std(y, 0,2)/sqrt(length(y));
%         X = shuffle_data(linspace(x-buff, x+buff, length(y)));
%         scatter(X, y, SZ, kolor, 'filled')
%         plot([x-buff, x+buff],[y_avg, y_avg], 'color', kolor, 'linewidth', LW)
%         errorbar(x,y_avg,y_err,'color', kolor, 'linewidth', eLW)
% 
%         % Save stats data: 
%         switch type
%             case 1 % warming
%                 SD.w(1:length(y),i) = y;
%                 SD.w_temps(i) = x;
%             case 2 % cooling
%                 SD.c(1:length(y),i) = y;
%                 SD.c_temps(i) = x;
%         end
%     end
%     ylimits = [ylimits, ylim]; %#ok<AGROW>
% end
% ylimits = [min(ylimits), max(ylimits)];
% 
% % Formatting
% formatFig(fig, blkbgd,[r,c]);
% for i = 1:2 % 
%     subplot(r,c,i)
%     title(sbpt_str{i})
%     xlabel('temp \circC')
%     ylabel(y_lab)
%     if ~autoLim 
%         ylim(ylimits)
%     end
%     if ~xlim_auto
%         xlim(temp_lims)
%     end
%     h_line(nullD, 'grey','--', 1)
% end
% 
% % Add statistics here...compare across temps within food and no food trials
% % ANOVA with post-hoc multicompare
% 
% % Are there differences between the food and no food conditions for each temperature? 
% [pC, pW, C_TP, W_TP] = deal([]);
% for i = 1:size(foodPairs,1)
%     food = foodPairs(i,1);
%     empty = foodPairs(i,2);
%     [~,pC(i)] = ttest(SD.c(:,food), SD.c(:,empty));
%     [~,pW(i)] = ttest(SD.w(:,food), SD.w(:,empty));
%     C_TP(i) = SD.c_temps(food); % working here
%     W_TP(i) = SD.w_temps(food); % working here
% end
% %Bonferonni correction:
% alpha = 0.05;
% m = num.exp;
% p_limit = alpha/m;
% hC = pC<=p_limit; % 'cooling' comparison btwn food and empty
% hW = pW<=p_limit; % 'warming' comparison btwn food and empty
% 
% % plot the significance stars on the data: 
% subplot(r, c, 1) % warming
%     ySig = rangeLine(fig,1,true);
%     XX = W_TP;
%     YY = ySig*ones(size(W_TP));
%     XX(~hW) = []; % remove non-significant comparisons
%     YY(~hW) = []; % remove non-significant comparisons
%     scatter(XX,YY, 200, foreColor,'filled', 'pentagram')
% subplot(r, c, 2) % cooling
%     ySig = rangeLine(fig,1,true);
%     XX = C_TP;
%     YY = ySig*ones(size(C_TP));
%     XX(~hC) = []; % remove non-significant comparisons
%     YY(~hC) = []; % remove non-significant comparisons
%     scatter(XX,YY, 200, foreColor,'filled', 'pentagram')
% 
% % Save the figure
% fig_folder = createFolder([saveDir 'Figures/']);
% save_figure(fig,[fig_folder title_str ' occ scatter tuning food vs no food ' temp_type],fig_type);
% 
% % statistical comparisons across the different temperatures -- set up the
% % groups and identifiers
% 
% temps = string(num2str(SD.w_temps));
% temps = repmat(temps',size(SD.w,1),1);
% temps = reshape(temps,[numel(temps),1]);
% foodList = repmat("Food",num.exp,1);
% foodList([data(:).emptytrial]) = "Empty";
% foodList = repmat(foodList',size(SD.w,1),1);
% foodList = reshape(foodList,[numel(foodList),1]);
% warm = reshape(SD.w,[numel(SD.w),1]);
% cool = reshape(SD.c,[numel(SD.c),1]);
% warm_idx = ~isnan(warm);
% cool_idx = ~isnan(cool);
% wTL = temps(warm_idx);
% wFL = foodList(warm_idx);
% cTL = temps(cool_idx);
% cFL = foodList(cool_idx);
% warm = warm(warm_idx);
% cool = cool(cool_idx);
% 
% % TODO: need to pull these out for each individual item....
% w_aov = anova({wTL, wFL}, warm,FactorNames=["Temperature" "Food_Status"]);
% c_aov = anova({cTL, cFL}, cool,FactorNames=["Temperature" "Food_Status"]);
% 
% m = multcompare(w_aov,["Temperature" "Food_Status"], CriticalValueType="bonferroni");
% 
% % TODO: make this into a comparison matrix to show which are different between the conditions
% 
% p = anovan(warm,{wTL, wFL},'model','interaction','varnames',{'Temperature','Food_Status'});
% 
% [~,~,stats] = anovan(warm,{wTL, wFL},'model','interaction','varnames',{'Temperature','Food_Status'});
% 
% 
% 


























