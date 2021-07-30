


%% Load excel file
clear
% fileRoot = 'F:\JeanneRotation\Evyn\';
fileRoot = 'E:\Evyn\';

try
    xlFile = 'E:\Evyn\Experiment Summary.xlsx';
    [~,~,excelfile] = xlsread(xlFile,'Sheet1');
catch 
    xlFile = 'F:\JeanneRotation\Evyn\Experiment Summary.xlsx';
    [~,~,excelfile] = xlsread(xlFile,'Sheet1');
end
% Excel column for the headers:
Excel.headers = excelfile(1,:); %xltitles
Excel.date = find(strcmpi('Date',Excel.headers) == 1);
Excel.flynum = find(strcmpi('Fly Num',Excel.headers) == 1);
Excel.Struct = find(strcmpi('Struct',Excel.headers) == 1);

% Select structure to analyze:
structList = unique(excelfile(2:end,Excel.Struct));
idx = listdlg('PromptString', 'Select group', 'ListString', structList);
flygroup = structList{idx};
disp(['Selected ' flygroup])
loc = find(strcmpi(excelfile(:,Excel.Struct),flygroup));
num.flies = length(loc);

% Load the fly data:
for ifly = 1:num.flies
    folderdate = excelfile{loc(ifly),Excel.date};
    flynum = excelfile{loc(ifly),Excel.flynum};
    flyID = generate_flyID(folderdate, num2str(flynum));
    FLY = load([fileRoot folderdate '\' flyID '_analysisData.mat']);
    % remove select fields
    fields = fieldnames(FLY);
    fullfields = {'ans','strtNum','I','A','vidName','y','flyList','v_temp','vidNum'};
    fieldList = fields(ismember(fields,fullfields));
    FLY = rmfield(FLY,fieldList);
    fly(ifly) = FLY;
end
clear FLY loc idx ans

figDir = [fileRoot 'matlabroot\Jeanne Rotatation\Figures\'];

initial_vars = who; initial_vars{end+1} = 'initial_vars';
clearvars('-except',initial_vars{:})


%% Fig: Scatter plot peak df_f
fileTag = '-png';               %image save type
stats_text = true;              %line with p-vals and sig stars
color_opt = true;               %black vs white background
LW = 1;                         %line width
SZ = 50;                        %scatter point size
l_style = '-';                  %line style for connecting line
paired = true;                  %connecting line
G1.color = Color('purple');     %Solvent
G2.color = Color('teal');       %CS-
G3.color = Color('orange');     %CS+
markerList = {'o','s','d','<','p','h','^','v','>','+','*','x'}; % marker list: 
% subplotInd(1).idx = [1:3];
% subplotInd(2).idx = [4,5];
% nrows = 1;
% ncols = 5;
subplotInd(1).idx = [1:3,6:8,11:13];
subplotInd(2).idx = [4,5,9,10,14,15];
% subplotInd(3).idx = [4,5];
nrows = 3;
ncols = 5;

if color_opt==true
    baseColor = 'white';
else
    baseColor = 'black';
end
param = fly(1).param;

% range for maxF
fps = 100/3;
ROI = round(param.odorON*fps):round((param.odorON+5)*fps);
sSpan = round(fps/2);

% find the max df_f
for ifly = 1:num.flies 
    temp = fly(ifly).solvent.naive;
    for ii = 1:size(temp,2)
        G1.N(ii,ifly) = mean(max(smooth(fly(ifly).solvent.naive(ROI,ii),sSpan)));
        G1.T(ii,ifly) = mean(max(smooth(fly(ifly).solvent.test(ROI,ii),sSpan)));
    end

    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G2.N(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.naive(ROI,ii),sSpan)));
        G2.T(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.test(ROI,ii),sSpan)));
    end
        
    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G3.N(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.naive(ROI,ii),sSpan)));
        G3.T(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.test(ROI,ii),sSpan)));
    end
end
G1.naive = mean(G1.N,1);
G2.naive = mean(G2.N,1);
G3.naive = mean(G3.N,1);
G1.test = mean(G1.T,1);
G2.test = mean(G2.T,1);
G3.test = mean(G3.T,1);

% Find the difference in the rate of change...
G1.diff = G1.test-G1.naive;
G2.diff = G2.test-G2.naive;
G3.diff = G3.test-G3.naive;
G1.avg_diff = mean(G1.diff);
G2.avg_diff = mean(G2.diff);
G3.avg_diff = mean(G3.diff);

% Find the absolute percent change in df/f:
G1.perchange = nanmean(abs(((G1.T-G1.N)./G1.N).*100),1);
G2.perchange = nanmean(abs(((G2.T-G2.N)./G2.N).*100),1);
G3.perchange = nanmean(abs(((G3.T-G3.N)./G3.N).*100),1);
G1.avgPC = median(G1.perchange);
G2.avgPC = median(G2.perchange);
G3.avgPC = median(G3.perchange);

% plot the peak
fig = getfig('Avg df_f',1);
subplot(nrows,ncols,subplotInd(1).idx);
hold on
% Solvent
len = length(G1.naive);
for ifly = 1:len
    scatter([1,2],[G1.naive(ifly),G1.test(ifly)],SZ,G1.color,'filled',markerList{ifly});
end

% CS-
len = length(G2.naive);
for ifly = 1:len
    scatter([3,4],[G2.naive(ifly),G2.test(ifly)],SZ,G2.color,'filled',markerList{ifly});
end

% CS+
len = length(G3.naive);
for ifly = 1:len
    scatter([5,6],[G3.naive(ifly),G3.test(ifly)],SZ,G3.color,'filled',markerList{ifly});
end

% connecting lines
if paired==true
    plot([1,2], [G1.naive;G1.test],'color', G1.color,'linewidth',LW,'linestyle',l_style)
    plot([3,4], [G2.naive;G2.test],'color', G2.color,'linewidth',LW,'linestyle',l_style)
    plot([5,6], [G3.naive;G3.test],'color', G3.color,'linewidth',LW,'linestyle',l_style)
end

% quick stats:
[~,G1.p] = ttest(G1.naive, G1.test);
[~,G2.p] = ttest(G2.naive, G2.test);
[~,G3.p] = ttest(G3.naive, G3.test);
p_val = [G1.p,G2.p,G3.p];
odors = {param.solvent, param.csNEG_odor, param.csPOS_odor};
% multiple comparisons test:
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = 1:3
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\nOdor ' odors{p_loc(idx)} ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
end
fprintf('\nStats Done\n')

% labels
ylabel('\DeltaF/F')
title({param.cross;'Naive vs Test'})
xlim([0,7])
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.5,3.5,5.5]);
set(ax, 'XTickLabels', {param.solvent,['CS- ' param.csNEG_odor],['CS+ ' param.csPOS_odor]})


subplot(nrows,ncols,subplotInd(2).idx);
hold on
offset = 0.5/num.flies;
for ii = 1:num.flies
    x = 1+(ii-1)*offset;
    scatter(x,G1.diff(ii),SZ,G1.color,'filled',markerList{ii});
end
for ii = 1:num.flies
    x = 2+(ii-1)*offset;
    scatter(x,G2.diff(ii),SZ,G2.color,'filled',markerList{ii});
end
for ii = 1:num.flies
    x = 3+(ii-1)*offset;
    scatter(x,G3.diff(ii),SZ,G3.color,'filled',markerList{ii});
end
%plot avg line
plot([0.9,1.6],[G1.avg_diff,G1.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
plot([1.9,2.6],[G2.avg_diff,G2.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
plot([2.9,3.6],[G3.avg_diff,G3.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
xlim([0,4])

% labels
% title('Change in \DeltaF/F from conditioning')
ylabel('Change in \DeltaF/F')
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.25,2.25,3.25])
set(ax, 'XTickLabels', {'Solvent', 'CS-','CS+'})

% % ANOVA of the absolute change in df/f:
% [p,tbl,stats] = anova1([G1.diff',G2.diff',G3.diff'],{'Solvent', 'CS-','CS+'})

[~,G1.p] = ttest(G1.naive, G1.test);
[~,G2.p] = ttest(G2.naive, G2.test);
[~,G3.p] = ttest(G3.naive, G3.test);
p_val = [G1.p,G2.p,G3.p];
odors = {param.solvent, param.csNEG_odor, param.csPOS_odor};
% multiple comparisons test:
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = 1:3
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\nOdor ' odors{p_loc(idx)} ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
end
fprintf('\nStats Done\n')


% ====== change in df/f =========
fig = formatFig(fig,color_opt,[nrows,ncols],subplotInd);
subplot(nrows,ncols,subplotInd(1).idx);

% Plot option
if stats_text == true 
   ax = gca;
   ub = ax.YTick;
   rng = range([ub(1),ub(end)])*0.03; %significance star offset
   pk = ub(end)+mean(diff(ub))+rng;
   % Solvent: 
   x = [1,2];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(round(G1.p,3))];
   t = text(1.25, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if sig(2)==true
       T = text(1.5, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
   % Unconditioned: 
   x = [3,4];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(round(G2.p,3))];
   t = text(3.25, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if sig(1)==true
       T = text(3.5, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
   % Conditioned: 
   x = [5,6];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(round(G3.p,3))];
   t = text(5.25, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if sig(2)==true
       T = text(5.5, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
end

save_figure(fig, [figDir flygroup ' df_f comparison'], fileTag);
clearvars('-except',initial_vars{:})


%% Fig: naive->test timecourse overlay
G1.color = 'teal';      % CSneg color
G2.color = 'orange';    % CSpos color
LW = 1.5;
color_opt = true;
fileTag = '-png';
lines = false;          % single line for each fly trial avg
shading = true;
norm_axes = true;
if color_opt==true
    baseColor = 'white';
else 
    baseColor = 'black';
end    
param = fly(1).param;
time = fly(1).time;
nrows = 2;
ncols = 1;    
    
for ifly = 1:num.flies
    G1.pre(:,ifly) = mean(fly(ifly).CSneg.naive(:,1:end),2);
    G1.post(:,ifly) = mean(fly(ifly).CSneg.test(:,1:end),2);
    
    G2.pre(:,ifly) = mean(fly(ifly).CSpos.naive(:,1:end),2);
    G2.post(:,ifly) = mean(fly(ifly).CSpos.test(:,1:end),2);
end

fig = getfig('Odors Responses',1);
% CS-
subplot(nrows,ncols,1); hold on
    % Naive
    y = G1.pre;
    err = std(y,0,2);
    avg = mean(y,2);
    if lines==true
        plot(time,y, 'color', Color('grey'),'linewidth', 0.5,'linestyle', ':')
    end
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(baseColor), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(baseColor),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(baseColor),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(baseColor),'linewidth', LW)
    % Conditioned
    y = G1.post;
    err = std(y,0,2);
    avg = mean(y,2);
    if lines==true
        plot(time,y, 'color', Color(G1.color),'linewidth', 0.5,'linestyle',':')
    end
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G1.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(G1.color),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(G1.color),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(G1.color),'linewidth', LW)
    % labels
    if norm_axes==false
        y1 = rangeLine(fig);
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
    end
    ymax(1,:) = ylim;
    xlabel('time (s)')
    ylabel('\DeltaF/F')
    title({['CS-  ' param.csNEG_odor];''})  

% CS-
subplot(nrows,ncols,2); hold on
    % Naive
    y = G2.pre;
    err = std(y,0,2);
    avg = mean(y,2);
    if lines==true
        plot(time,y, 'color', Color('grey'),'linewidth', 0.5,'linestyle', ':')
    end
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(baseColor), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(baseColor),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(baseColor),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(baseColor),'linewidth', LW)
    % Conditioned
    y = G2.post;
    err = std(y,0,2);
    avg = mean(y,2);
    if lines==true
        plot(time,y, 'color', Color(G2.color),'linewidth', 0.5,'linestyle',':')
    end
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G2.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(G2.color),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(G2.color),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(G2.color),'linewidth', LW)
    % labels
    if norm_axes==false
        y1 = rangeLine(fig);
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
    end
    ymax(2,:) = ylim;
    xlabel('time (s)')
    ylabel('\DeltaF/F')
    title({['CS+  ' param.csPOS_odor];''})    
    

if norm_axes==true
    ub = max(ymax(:,2));
    lb = min(ymax(:,1));
    for n = 1:(nrows*ncols)
        subplot(nrows,ncols,n)
        ylim([lb,ub])
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [ub,ub], 'color',Color(baseColor),'linewidth', 3)
    end
end

fig = formatFig(fig, color_opt,[nrows,ncols]);
save_figure(fig, [figDir flygroup ' df_f overlay'], fileTag);

clearvars('-except',initial_vars{:})
       

%% Summary : time course + change in df_f
G1.color = 'teal';              % CSneg color
G2.color = 'orange';            % CSpos color
color_opt = true;
fileTag = '-png';               %image save time
shading = false;                %error shading
norm_axes = true;               %normalize axes on time course figs
stats_text = true;              %line with p-vals and sig stars
LW = 1;                         %linewidth for plots
SZ = 50;                        %scatter point size
paired = true;                  %connecting line
markerList = {'o','s','d','<','p','h','^','v','>','+','*','x'}; % marker list: 

if color_opt==true
    baseColor = 'white';
else 
    baseColor = 'black';
end    
param = fly(1).param;
time = fly(1).time;

nrows = 2;
ncols = 5;    
subplotInd(1).idx = [1:3];
subplotInd(2).idx = [6:8];
subplotInd(3).idx = [4,5,9,10];

fig = getfig('Odors Responses',1);
% ========== Time course ==========
for ifly = 1:num.flies
    G1.pre(:,ifly) = mean(fly(ifly).CSneg.naive(:,1:end),2);
    G1.post(:,ifly) = mean(fly(ifly).CSneg.test(:,1:end),2);
    
    G2.pre(:,ifly) = mean(fly(ifly).CSpos.naive(:,1:end),2);
    G2.post(:,ifly) = mean(fly(ifly).CSpos.test(:,1:end),2);
end
% CS-
subplot(nrows,ncols,subplotInd(1).idx); hold on
    % Naive
    y = G1.pre;
    err = std(y,0,2);
    avg = mean(y,2);
    plot(time,y, 'color', Color('grey'),'linewidth', 0.5,'linestyle', ':')
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(baseColor), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(baseColor),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(baseColor),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(baseColor),'linewidth', LW)
    % Conditioned
    y = G1.post;
    err = std(y,0,2);
    avg = mean(y,2);
    plot(time,y, 'color', Color(G1.color),'linewidth', 0.5,'linestyle',':')
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G1.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(G1.color),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(G1.color),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(G1.color),'linewidth', LW)
    % labels
    if norm_axes==false
        y1 = rangeLine(fig);
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
    end
    ymax(1,:) = ylim;
    xlabel('time (s)')
    ylabel('\DeltaF/F')
    title({['CS-  ' param.csNEG_odor];''})  

% CS-
subplot(nrows,ncols,subplotInd(2).idx); hold on
    % Naive
    y = G2.pre;
    err = std(y,0,2);
    avg = mean(y,2);
    plot(time,y, 'color', Color('grey'),'linewidth', 0.5,'linestyle', ':')
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(baseColor), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(baseColor),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(baseColor),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(baseColor),'linewidth', LW)
    % Conditioned
    y = G2.post;
    err = std(y,0,2);
    avg = mean(y,2);
    plot(time,y, 'color', Color(G2.color),'linewidth', 0.5,'linestyle',':')
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G2.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(time, avg+err, 'color', Color(G2.color),'linewidth', 0.5)
        plot(time, avg-err, 'color', Color(G2.color),'linewidth', 0.5)
    end
    plot(time, avg, 'color', Color(G2.color),'linewidth', LW)
    % labels
    if norm_axes==false
        y1 = rangeLine(fig);
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
    end
    ymax(2,:) = ylim;
    xlabel('time (s)')
    ylabel('\DeltaF/F')
    title({['CS+  ' param.csPOS_odor];''})    
    

if norm_axes==true
    ub = max(ymax(:,2));
    lb = min(ymax(:,1));
    for n = 1:2
        subplot(nrows,ncols,subplotInd(n).idx)
        ylim([lb,ub])
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [ub,ub], 'color',Color(baseColor),'linewidth', 3)
    end
end

% ========== Scatter plot ==========
% range for maxF
fps = 100/3;
ROI = round(param.odorON*fps):round((param.odorON+5)*fps);
sSpan = round(fps/2);

% find the max df_f
for ifly = 1:num.flies 
    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G1.N(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.naive(ROI,ii),sSpan)));
        G1.T(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.test(ROI,ii),sSpan)));
    end
        
    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G2.N(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.naive(ROI,ii),sSpan)));
        G2.T(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.test(ROI,ii),sSpan)));
    end
end
G1.naive = mean(G1.N,1);
G2.naive = mean(G2.N,1);
G1.test = mean(G1.T,1);
G2.test = mean(G2.T,1);

% plot the peak
subplot(nrows,ncols,subplotInd(3).idx)
hold on
% CS-
len = length(G1.naive);
for ifly = 1:len
    scatter([1,2],[G1.naive(ifly),G1.test(ifly)],SZ,Color(G1.color),'filled',markerList{ifly});
end
% CS+
len = length(G2.naive);
for ifly = 1:len
    scatter([3,4],[G2.naive(ifly),G2.test(ifly)],SZ,Color(G2.color),'filled',markerList{ifly});
end
% connecting lines
if paired==true
    plot([1,2], [G1.naive;G1.test],'color', Color(G1.color),'linewidth',LW)
    plot([3,4], [G2.naive;G2.test],'color', Color(G2.color),'linewidth',LW)
end

% quick stats:
[~,G1.p] = ttest(G1.naive, G1.test);
[~,G2.p] = ttest(G2.naive, G2.test);
p_val = [G1.p,G2.p];
odors = {param.csNEG_odor, param.csPOS_odor};
% multiple comparisons test:
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = 1:length(p_err)
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\nOdor ' odors{p_loc(idx)} ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
end
fprintf('\nStats Done\n')

% labels
ylabel('\DeltaF/F')
title({param.cross;'Naive vs Test'})
xlim([0.5,4.5])
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.5,3.5])
set(ax, 'XTickLabels', {['CS- ' param.csNEG_odor],['CS+ ' param.csPOS_odor]})

fig = formatFig(fig,color_opt,[nrows,ncols],subplotInd);

% Plot option
if stats_text == true 
   ax = gca;
   ub = ax.YTick;
   rng = range([ub(1),ub(end)])*0.03; %significance star offset
   pk = ub(end)+mean(diff(ub))+rng;
   % Unconditioned: 
   x = [1,2];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(round(G1.p,3))];
   t = text(1.1, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if sig(1)==true
       T = text(1.25, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
   % Conditioned: 
   x = [3,4];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(round(G2.p,3))];
   t = text(3.1, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if sig(2)==true
       T = text(3.25, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
end

save_figure(fig, [figDir flygroup ' summary'], fileTag);

clearvars('-except',initial_vars{:})


%% Plot the df/f for each airpuff over time
fileTag = '-png';               %image save type
stats_text = true;              %line with p-vals and sig stars
color_opt = true;               %black vs white background
LW = 1;                         %line width
SZ = 50;                        %scatter point size
l_style = '-';                  %line style for connecting line
paired = true;                  %connecting line
G1.color = Color('purple');     %Solvent
G2.color = Color('teal');       %CS-
G3.color = Color('orange');     %CS+

if color_opt==true
    baseColor = 'white';
    highColor = 'black';
else
    baseColor = 'black';
    highColor = 'white';
end
param = fly(1).param;

% range for maxF
fps = 100/3;
ROI = round(param.odorON*fps):round((param.odorON+5)*fps);
sSpan = round(fps/2);
% find the max df_f
for ifly = 1:num.flies 
    temp = fly(ifly).solvent.naive;
    for ii = 1:size(temp,2)
        G1.N(ii,ifly) = mean(max(smooth(fly(ifly).solvent.naive(ROI,ii),sSpan)));
        G1.T(ii,ifly) = mean(max(smooth(fly(ifly).solvent.test(ROI,ii),sSpan)));
    end

    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G2.N(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.naive(ROI,ii),sSpan)));
        G2.T(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.test(ROI,ii),sSpan)));
    end
        
    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G3.N(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.naive(ROI,ii),sSpan)));
        G3.T(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.test(ROI,ii),sSpan)));
    end
end
G1.naive = mean(G1.N,1);
G2.naive = mean(G2.N,1);
G3.naive = mean(G3.N,1);
G1.test = mean(G1.T,1);
G2.test = mean(G2.T,1);
G3.test = mean(G3.T,1);

markerList = {'o','s','d','<','p','h','^','v','>','+','*','x'};


% plot the peak
fig = getfig('Trial based df_f',1);
hold on
idx = 1;

% Solvent
len = size(G1.N,1); %number of trials
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G1.N(:,ii),'linewidth', 1, 'color',Color(baseColor),...
         'marker', markerList{ii},'markerfacecolor',Color(baseColor))
end
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G1.T(:,ii),'linewidth', 1, 'color',G1.color,'marker', markerList{ii},'markerfacecolor',G1.color)
end
idx = idx+len;

% CS-
len = size(G2.N,1); %number of trials
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G2.N(:,ii),'linewidth', 1, 'color',Color(baseColor),...
         'marker', markerList{ii},'markerfacecolor',Color(baseColor))
end
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G2.T(:,ii),'linewidth', 1, 'color',G2.color,'marker', markerList{ii},'markerfacecolor',G2.color)
end
idx = idx+len;

% CS-
len = size(G3.N,1); %number of trials
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G3.N(:,ii),'linewidth', 1, 'color',Color(baseColor),...
         'marker', markerList{ii},'markerfacecolor',Color(baseColor))
end
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G3.T(:,ii),'linewidth', 1, 'color',G3.color,'marker', markerList{ii},'markerfacecolor',G3.color)
end
idx = idx+len;

% labels
ylabel('\DeltaF/F')
title({param.cross;'Naive vs Test'})
xlim([0,idx])
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.5,5,10])
set(ax, 'XTickLabels', {param.solvent,['CS- ' param.csNEG_odor],['CS+ ' param.csPOS_odor]})

fig = formatFig(fig,color_opt);

save_figure(fig, [figDir flygroup ' Trial-based df-f'], fileTag);
clearvars('-except',initial_vars{:})


%% Plot the peak df/f for TRAINING sessions && SLOPE over time
fileTag = '-png';               %image save type
stats_text = true;              %line with p-vals and sig stars
color_opt = true;               %black vs white background
LW = 1;                         %line width
SZ = 50;                        %scatter point size
l_style = '-';                  %line style for connecting line
paired = true;                  %connecting line
G2.color = Color('teal');       %CS-
G3.color = Color('orange');     %CS+
markerList = {'o','s','d','<','p','h','^','v','>','+','*','x'};
sbpltInd(1).idx = [1:3];
sbpltInd(2).idx = [4:5];
nrows = 1;
ncols = 5;

if color_opt==true
    baseColor = 'white';
    highColor = 'black';
else
    baseColor = 'black';
    highColor = 'white';
end
param = fly(1).param;

% range for maxF
fps = 100/3;
ROI = round(param.odorON*fps):round((param.odorON+5)*fps);
sSpan = round(fps/2);
% find the max df_f
for ifly = 1:num.flies 
    temp = fly(ifly).CSneg.training;
    for ii = 1:size(temp,2)
        G2.T(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.training(ROI,ii),sSpan)));
    end
        
    temp = fly(ifly).CSneg.training;
    for ii = 1:size(temp,2)
        G3.T(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.training(ROI,ii),sSpan)));
    end
end
G2.test = mean(G2.T,1);
G3.test = mean(G3.T,1);

% ========= plot the df/f data ============
fig = getfig('Training session df_f',1);
subplot(nrows,ncols,sbpltInd(1).idx)
hold on
idx = 1;

% CS-
len = size(G2.T,1); %number of trials
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G2.T(:,ii),'linewidth', 1, 'color',G2.color,'marker', markerList{ii},'markerfacecolor',G2.color)
end
idx = idx+len;

% CS-
len = size(G3.T,1); %number of trials
for ii = 1:num.flies
    x = idx:idx+len-1;
    plot(x,G3.T(:,ii),'linewidth', 1, 'color',G3.color,'marker', markerList{ii},'markerfacecolor',G3.color)
end
idx = idx+len;

% labels
ylabel('\DeltaF/F')
title({param.cross;'Training trials'})
xlim([0,idx])
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.5,5,10])
set(ax, 'XTickLabels', {['CS- ' param.csNEG_odor],['CS+ ' param.csPOS_odor]})

% =========== Slope analysis ============
subplot(nrows,ncols,sbpltInd(2).idx)
    hold on
    % find the difference in df/f for each of the odor conditions 
    G2.diff = sum(diff(G2.T,1),1);
    G3.diff = sum(diff(G3.T,1),1);
    G2.avg_diff = mean(G2.diff);
    G3.avg_diff = mean(G3.diff);

    offset = 0.5/num.flies;
    for ii = 1:num.flies
        x = 1+(ii-1)*offset;
        scatter(x,G2.diff(ii),SZ,G2.color,'filled',markerList{ii});
    end
    for ii = 1:num.flies
        x = 2+(ii-1)*offset;
        scatter(x,G3.diff(ii),SZ,G3.color,'filled',markerList{ii});
    end
    %plot avg line
    plot([0.9,1.6],[G2.avg_diff,G2.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
    plot([1.9,2.6],[G3.avg_diff,G3.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
    xlim([0,3])

    % labels
    title('Change in \DeltaF/F during training')
    ylabel('slope of (\DeltaF/F)')
    xlabel('Odor')
    ax = gca;
    set(ax, 'XTick', [1.25,2.25])
    set(ax, 'XTickLabels', {'CS-','CS+'})

fig = formatFig(fig,color_opt,[nrows,ncols],sbpltInd);

% Stats comparison
[~,pval] = ttest(G2.diff, G3.diff);
disp(['CS+ vs CS- p-value: ' num2str(pval)])

% Plot option
if stats_text == true 
   subplot(nrows,ncols,sbpltInd(2).idx)
   ax = gca;
   ub = ax.YTick;
   rng = range([ub(1),ub(end)])*0.03; %significance star offset
   pk = ub(end)+mean(diff(ub))+rng;
   % Unconditioned: 
   x = [0.9,2.1];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(pval)];
   t = text(1.25, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if pval<0.05
       T = text(1.5, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
end

save_figure(fig, [figDir flygroup ' Training sessions'], fileTag);
clearvars('-except',initial_vars{:})


%% Signal decay
fileTag = '-png';               %image save type
stats_text = true;              %line with p-vals and sig stars
color_opt = true;               %black vs white background
LW = 1;                         %line width
SZ = 50;                        %scatter point size
l_style = '-';                  %line style for connecting line
paired = true;                  %connecting line
G1.color = Color('purple');     %Solvent
G2.color = Color('teal');       %CS-
G3.color = Color('orange');     %CS+
markerList = {'o','s','d','<','p','h','^','v','>','+','*','x'}; % marker list: 
sbpltInd(1).idx = [1:3];
sbpltInd(2).idx = [4:5];
nrows = 1;
ncols = 5;

if color_opt==true
    baseColor = 'white';
else
    baseColor = 'black';
end
param = fly(1).param;

% range for maxF
fps = 100/3;
ROI = round(param.odorON*fps):round((param.odorON+5)*fps);
sSpan = round(fps/2);

% find the max df_f
for ifly = 1:num.flies 
    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        data = fly(ifly).CSneg.naive(:,ii);
        G2.N(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.naive(ROI,ii),sSpan)));
        G2.N_dur(ifly) = Odor_response_duration(data,fps,ROI);
        data = fly(ifly).CSneg.test(:,ii);
        G2.T(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.test(ROI,ii),sSpan)));
        G2.T_dur(ifly) = Odor_response_duration(data,fps,ROI);
    end
        
    temp = fly(ifly).CSpos.naive;
    for ii = 1:size(temp,2)
        data = fly(ifly).CSpos.naive(:,ii);
        G3.N(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.naive(ROI,ii),sSpan)));
        G3.N_dur(ifly) = Odor_response_duration(data,fps,ROI);
        data = fly(ifly).CSpos.test(:,ii);
        G3.T(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.test(ROI,ii),sSpan)));
        G3.T_dur(ifly) = Odor_response_duration(data,fps,ROI);
    end
end

G2.naive = mean(G2.N_dur,1);
G3.naive = mean(G3.N_dur,1);
G2.test = mean(G2.T_dur,1);
G3.test = mean(G3.T_dur,1);
 
% Find the difference in the rate of change...
G2.diff = G2.test-G2.naive;
G3.diff = G3.test-G3.naive;
G2.avg_diff = mean(G2.diff);
G3.avg_diff = mean(G3.diff);

% ======== Plot the Duration data ========
% plot the peak
fig = getfig('Response duration',1);
subplot(nrows,ncols,sbpltInd(1).idx)
hold on
% CS-
len = length(G2.naive);
for ifly = 1:len
    scatter([1,2],[G2.naive(ifly),G2.test(ifly)],SZ,G2.color,'filled',markerList{ifly});
end

% CS+
len = length(G3.naive);
for ifly = 1:len
    scatter([3,4],[G3.naive(ifly),G3.test(ifly)],SZ,G3.color,'filled',markerList{ifly});
end

% connecting lines
if paired==true
    plot([1,2], [G2.naive;G2.test],'color', G2.color,'linewidth',LW,'linestyle',l_style)
    plot([3,4], [G3.naive;G3.test],'color', G3.color,'linewidth',LW,'linestyle',l_style)
end

% quick stats:
[~,G2.p] = ttest(G2.naive, G2.test);
[~,G3.p] = ttest(G3.naive, G3.test);
p_val = [G2.p,G3.p];
odors = {param.csNEG_odor, param.csPOS_odor};
% multiple comparisons test:
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = 1:2
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\nOdor ' odors{p_loc(idx)} ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
end
fprintf('\nStats Done\n')

% labels
ylabel('\DeltaF/F response duration (s)')
title({param.cross;'Naive vs Test'})
xlim([0,5])
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.5,3.5])
set(ax, 'XTickLabels', {['CS- ' param.csNEG_odor],['CS+ ' param.csPOS_odor]})

% =========== Slope analysis ============
subplot(nrows,ncols,sbpltInd(2).idx)
    hold on
    % find the difference in df/f for each of the odor conditions 
    G2.diff = sum(diff(G2.T,1),1);
    G3.diff = sum(diff(G3.T,1),1);
    G2.avg_diff = mean(G2.diff);
    G3.avg_diff = mean(G3.diff);

    offset = 0.5/num.flies;
    for ii = 1:num.flies
        x = 1+(ii-1)*offset;
        scatter(x,G2.diff(ii),SZ,G2.color,'filled',markerList{ii});
    end
    for ii = 1:num.flies
        x = 2+(ii-1)*offset;
        scatter(x,G3.diff(ii),SZ,G3.color,'filled',markerList{ii});
    end
    %plot avg line
    plot([0.9,1.6],[G2.avg_diff,G2.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
    plot([1.9,2.6],[G3.avg_diff,G3.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
    xlim([0.5,3])

    % labels
    ylabel('\Delta response duration (s)')
    xlabel('Odor')
    ax = gca;
    set(ax, 'XTick', [1.25,2.25])
    set(ax, 'XTickLabels', {'CS-','CS+'})

fig = formatFig(fig,color_opt,[nrows,ncols],sbpltInd);


% Stats for pt 1:
subplot(nrows,ncols,sbpltInd(1).idx)
if stats_text == true 
   ax = gca;
   ub = ax.YTick;
   rng = range([ub(1),ub(end)])*0.03; %significance star offset
   pk = ub(end)+mean(diff(ub))+rng;
   % Unconditioned: 
   x = [1,2];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(round(G2.p,3))];
   t = text(1.25, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if sig(1)==true
       T = text(1.5, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
   % Conditioned: 
   x = [3,4];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(round(G3.p,3))];
   t = text(3.25, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if sig(2)==true
       T = text(3.5, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
end


% Stats comparison
[~,pval] = ttest(G2.diff, G3.diff);
disp(['CS+ vs CS- p-value: ' num2str(pval)])

% Plot option
if stats_text == true 
   subplot(nrows,ncols,sbpltInd(2).idx)
   ax = gca;
   ub = ax.YTick;
   rng = range([ub(1),ub(end)])*0.03; %significance star offset
   pk = ub(end)+mean(diff(ub))+rng;
   x = [0.9,2.6];
   plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
   texttag = ['p = ' num2str(pval)];
   t = text(1.5, pk-rng, texttag);
   set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
   if pval<0.05
       T = text(1.75, pk+rng, '*');
       set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
   end
end

save_figure(fig, [figDir flygroup ' response duration'], fileTag);
clearvars('-except',initial_vars{:})


%% Relative change in df/f
fileTag = '-png';               %image save type
stats_text = true;              %line with p-vals and sig stars
color_opt = true;               %black vs white background
LW = 1;                         %line width
SZ = 50;                        %scatter point size
l_style = '-';                  %line style for connecting line
paired = true;                  %connecting line
G1.color = Color('purple');       %CS-
G2.color = Color('teal');       %CS-
G3.color = Color('orange');     %CS+
markerList = {'o','s','d','<','p','h','^','v','>','+','*','x'};
sbpltInd(1).idx = [1:3];
sbpltInd(2).idx = [4:5];
nrows = 1;
ncols = 5;

if color_opt==true
    baseColor = 'white';
    highColor = 'black';
else
    baseColor = 'black';
    highColor = 'white';
end
param = fly(1).param;

% range for maxF
fps = 100/3;
ROI = round(param.odorON*fps):round((param.odorON+5)*fps);
sSpan = round(fps/2);
% find the max df_f
for ifly = 1:num.flies 
    temp = fly(ifly).solvent.naive;
    for ii = 1:size(temp,2)
        G1.N(ii,ifly) = mean(max(smooth(fly(ifly).solvent.naive(ROI,ii),sSpan)));
        G1.T(ii,ifly) = mean(max(smooth(fly(ifly).solvent.test(ROI,ii),sSpan)));
    end

    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G2.N(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.naive(ROI,ii),sSpan)));
        G2.T(ii,ifly) = mean(max(smooth(fly(ifly).CSneg.test(ROI,ii),sSpan)));
    end
        
    temp = fly(ifly).CSneg.naive;
    for ii = 1:size(temp,2)
        G3.N(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.naive(ROI,ii),sSpan)));
        G3.T(ii,ifly) = mean(max(smooth(fly(ifly).CSpos.test(ROI,ii),sSpan)));
    end
end
G1.naive = mean(G1.N,1);
G2.naive = mean(G2.N,1);
G3.naive = mean(G3.N,1);
G1.test = mean(G1.T,1);
G2.test = mean(G2.T,1);
G3.test = mean(G3.T,1);

% Find the difference in the rate of change...
G1.diff = G1.test-G1.naive;
G2.diff = G2.test-G2.naive;
G3.diff = G3.test-G3.naive;
G1.avg_diff = mean(G1.diff);
G2.avg_diff = mean(G2.diff);
G3.avg_diff = mean(G3.diff);

% ========= plot the df/f data ============
fig = getfig('Relative rate of change in df_f',1);
hold on
offset = 0.5/num.flies;
for ii = 1:num.flies
    x = 1+(ii-1)*offset;
    scatter(x,G1.diff(ii),SZ,G1.color,'filled',markerList{ii});
end
for ii = 1:num.flies
    x = 2+(ii-1)*offset;
    scatter(x,G2.diff(ii),SZ,G2.color,'filled',markerList{ii});
end
for ii = 1:num.flies
    x = 3+(ii-1)*offset;
    scatter(x,G3.diff(ii),SZ,G3.color,'filled',markerList{ii});
end
%plot avg line
plot([0.9,1.6],[G1.avg_diff,G1.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
plot([1.9,2.6],[G2.avg_diff,G2.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
plot([2.9,3.6],[G3.avg_diff,G3.avg_diff],'color', Color(baseColor), 'linewidth', 1, 'linestyle', ':')
xlim([0,4])

% labels
title('Change in \DeltaF/F from conditioning')
ylabel('slope of (\DeltaF/F)')
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.25,2.25,3.25])
set(ax, 'XTickLabels', {'Solvent', 'CS-','CS+'})

fig = formatFig(fig,color_opt);

% % Stats comparison
% [~,pval] = ttest(G2.diff, G3.diff);
% disp(['CS+ vs CS- p-value: ' num2str(pval)])
% 
% % Plot option
% if stats_text == true 
%    subplot(nrows,ncols,sbpltInd(2).idx)
%    ax = gca;
%    ub = ax.YTick;
%    rng = range([ub(1),ub(end)])*0.03; %significance star offset
%    pk = ub(end)+mean(diff(ub))+rng;
%    % Unconditioned: 
%    x = [0.9,2.1];
%    plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
%    texttag = ['p = ' num2str(pval)];
%    t = text(1.25, pk-rng, texttag);
%    set(t,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',10);
%    if pval<0.05
%        T = text(1.5, pk+rng, '*');
%        set(T,'Color', Color(baseColor), 'FontName', 'Arial', 'FontSize',25);
%    end
% end

save_figure(fig, [figDir flygroup ' Absolute change in df_f'], fileTag);
clearvars('-except',initial_vars{:})








































