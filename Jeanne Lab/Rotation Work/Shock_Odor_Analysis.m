

%% Load TIFF image stack to analyze
clear
warning('off','all') % Suppress all the tiff warnings
folderdate = '10.8.20';
% addpath('G:\JeanneRotation\Evyn\matlabroot\Jeanne Rotatation')
% fileRoot = ['F:\JeanneRotation\Evyn\' folderdate '\'];
fileRoot = ['E:\Evyn\' folderdate '\'];

A = questdlg('Open previously generated file?','','Yes', 'No', 'Cancel', 'No');
switch A
    case 'No'
        [excelfile, Excel, xlFile] = load_ExperimentSummary;
        % Select fly to analyze
        a = excelfile(strcmpi(excelfile(:,Excel.date),folderdate),:);
        for ii = 1:length(a)
            b{ii} = num2str(a{ii,Excel.flynum});
        end
        flyList = unique(b);
        idx = listdlg('PromptString', 'Select fly', 'ListString', flyList);
        flynum = flyList{idx}; 
        flyID = generate_flyID(folderdate, flynum);
        loc = strcmpi(b,flynum);
        flyDir = a(loc,:);
        num.vids = size(flyDir,1);
        clear a b loc
        
        % need to add in a selection here for the appriate fly number!! 
        vidList = dir([fileRoot '*.tif*']);
        for ii = 1:length(vidList)
            TiffFiles{ii} = vidList(ii).name;
        end
        
        for ii = 1:num.vids
            vidNum(ii) = flyDir{ii,Excel.session};
            v_temp = regexp(TiffFiles, regexptranslate('wildcard', ['*' num2str(vidNum(ii)) ').tif']));
            for n = 1:length(v_temp)
                if v_temp{n}==1
                   vidName{ii} =  TiffFiles{n};
                   continue;
                end
            end
        end
        
        % import video for all selected files
        h = waitbar(0,'Loading videos');
        for ii = 1:num.vids
            disp(['Loading: ' vidName{ii}])
            % load stack:
            tstack = Tiff(fullfile(fileRoot,vidName{ii}),'r');
            [I,J] = size(tstack.read()); % get the pixel sizes for frame
            nframes = length(imfinfo(fullfile(fileRoot,vidName{ii})));
            IM = zeros(I,J,nframes);
            IM(:,:,1)  = tstack.read();
            for n = 2:nframes
                tstack.nextDirectory()
                IM(:,:,n) = tstack.read();
            end
            % Extract the intensity information for the given file
            vidmean = mean(mean(IM,1),2); % find the mean intensity per frame
            y = squeeze(vidmean);
            vid(ii).y = y;
%             vidpath = [fileRoot flyID '_' num2str(ii)];
%             save(vidpath, 'IM')
%             vidList(ii).vidpath = vidpath;
            % add full video file to vid structure
            clear IM tstack
            % Update waitbar
            waitbar(ii/num.vids,h)
        end
        %close waitbar
        close(h)
        disp('done!')
        clear ii n J idx h

        % Save the loaded image data:
        initial_vars = who; initial_vars{end+1} = 'initial_vars';
        save([fileRoot flyID '_videoData'])
        disp('Data saved!')
    case 'Yes'
        temp = dir([fileRoot '*videoData.mat']);
        for ii = 1:length(temp)
            tempFiles{ii} = temp(ii).name;
        end
        idx = listdlg('PromptString', 'Select data', 'ListString', tempFiles, 'ListSize', [200 100]);
        load(fullfile(fileRoot, temp(idx).name))
    case 'Cancel'
        return
end
        
clearvars('-except',initial_vars{:})


% Organize data into useable structures
% Load matlab data:
try 
    fly = load([fileRoot flyID '_data.mat']);
    param = fly.param;
    param.vid_del = 2;
catch
    disp('No matlab data found to load')
    return
end
initial_vars{end+1} = 'fly';
initial_vars{end+1} = 'CSneg';
initial_vars{end+1} = 'CSpos';
initial_vars{end+1} = 'solvent';
initial_vars{end+1} = 'param';
initial_vars{end+1} = 'time';
initial_vars{end+1} = 'trial';
initial_vars{end+1} = 'empty';
num = fly.num;

% Sort data into groups based on experiment phase:
emp = strcmpi(flyDir(:,Excel.experiment), 'Empty');
pre = strcmpi(flyDir(:,Excel.experiment), 'Control');
stim_shock = strcmpi(flyDir(:,Excel.experiment), 'CS+ training');
stim = strcmpi(flyDir(:,Excel.experiment), 'CS- training');
post = strcmpi(flyDir(:,Excel.experiment), 'Test');
% group all vidoes (assumes vids are the same length)
for n = 1:length(vid)
    intensity(:,n) = vid(n).y;
end

% find change in florescence for all trials:
odorON = param.ODOR_start - param.vid_del;
odorOFF = odorON + param.ODOR_dur;
param.odorON = odorON; param.odorOFF = odorOFF;
buff = 2; % extra start buffer for calculating base florescence
viddur = param.trial_duration-buff;
poststim = viddur-odorON;
num.fps = size(intensity,1)/viddur;
time = linspace(-param.odorON,poststim,size(intensity,1));
controlROI = round(buff*num.fps):round(odorON*num.fps);
F0 = mean(intensity(controlROI,:),1);
df_f = ((intensity-F0)./F0)*100;
trial.intensity = intensity;
trial.F0 = F0;
trial.df_f = df_f;
trial.controlROI = controlROI;
trial.time = time;
% figure;
% plot(df_f)

% ------ Empty ------
try 
    param.empty = true;
    empty.odor = 'empty';
    empty.loc = strcmpi(flyDir(:,Excel.odor),empty.odor) & emp;
    empty.df_f = df_f(:,empty.loc);
catch
    empty = [];
    param.empty = false;
    disp('No empty vial trials found')
end

% ------ Solvent ------
solvent.odor = param.solvent;
solvent.pre_loc = strcmpi(flyDir(:,Excel.odor),solvent.odor) & pre;
solvent.naive = df_f(:,solvent.pre_loc);
solvent.post_loc = strcmpi(flyDir(:,Excel.odor),solvent.odor) & post;
solvent.test = df_f(:,solvent.post_loc);

% ------ CSpos ------ 
CSpos.odor = param.csPOS_odor;
CSpos.pre_loc = strcmpi(flyDir(:,Excel.odor),CSpos.odor) & pre;
CSpos.naive = df_f(:,CSpos.pre_loc);
CSpos.train_loc = strcmpi(flyDir(:,Excel.odor),CSpos.odor) & stim_shock;
CSpos.training = df_f(:,CSpos.train_loc);
CSpos.post_loc = strcmpi(flyDir(:,Excel.odor),CSpos.odor) & post;
CSpos.test = df_f(:,CSpos.post_loc);

% ------ CSneg ------ 
CSneg.odor = param.csNEG_odor;
CSneg.pre_loc = strcmpi(flyDir(:,Excel.odor),CSneg.odor) & pre;
CSneg.naive = df_f(:,CSneg.pre_loc);
CSneg.train_loc = strcmpi(flyDir(:,Excel.odor),CSneg.odor) & stim;
CSneg.training = df_f(:,CSneg.train_loc);
CSneg.post_loc = strcmpi(flyDir(:,Excel.odor),CSneg.odor) & post;
CSneg.test = df_f(:,CSneg.post_loc);


clearvars('-except',initial_vars{:})

% save all data:
save([fileRoot flyID '_analysisData'])

%% Fig: Empty vial pre-post comparision
if param.empty == true
    G1.color = 'black';
    G2.color = 'purple';
    LW = 1;
    color_opt = true;
    fileTag = '-png';
    fig_title = [param.cross ' Empty vial control'];
    setY = true;
    y_range = [-5,25];

    if color_opt==true
        baseColor = 'white';
        highColor = 'black';
        if strcmpi(G1.color,'black')
            G1.color = baseColor;
        end
        if strcmpi(G2.color,'black')
            G1.color = baseColor;
        end
    else 
        baseColor = 'black';
        highColor = 'white';
    end

    fig = getfig('Empty Controls',1);
    hold on
    % pre
    plot(time, empty.df_f(:,1), 'color', Color(G1.color), 'linewidth', LW)
    % post
    plot(time, empty.df_f(:,2), 'color', Color(G2.color), 'linewidth', LW)
    % odor & labels
    if setY==true
        ylim(y_range)
    end
    y1 = rangeLine(fig);
    odorOFF = param.ODOR_dur;
    plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
    xlabel('time (s)')
    ylabel('\DeltaF/F')
    title({fig_title; strrep(flyID,'_','-')})
    l = legend({'first','last'});
    set(l,'Color', Color(highColor), 'TextColor', Color(baseColor),'EdgeColor',Color(baseColor));

    fig = formatFig(fig, color_opt);

    % save figure option here...
    save_figure(fig, [fileRoot flyID ' empty vial time course'], fileTag);
    clearvars('-except',initial_vars{:})
end

%% Fig: Solvent pre-post comparision
G1.color = 'black';
G2.color = 'purple';
LW = 1;
color_opt = true;
fileTag = '-png';
fig_title = [param.cross ' ' param.solvent ' control'];
setY = true;
y_range = [-5,25];

if color_opt==true
    baseColor = 'white';
    highColor = 'black';
    if strcmpi(G1.color,'black')
        G1.color = baseColor;
    end
    if strcmpi(G2.color,'black')
        G1.color = baseColor;
    end
else 
    baseColor = 'black';
    highColor = 'white';
end

fig = getfig('Solvent Controls',1);
hold on
% pre
plot(time, solvent.naive, 'color', Color(G1.color), 'linewidth', LW)
% post
plot(time, solvent.test, 'color', Color(G2.color), 'linewidth', LW)
% odor & labels
if setY==true
    ylim(y_range)
end
y1 = rangeLine(fig);
odorOFF = param.ODOR_dur;
plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
xlabel('time (s)')
ylabel('\DeltaF/F')
title({fig_title; strrep(flyID,'_','-')})
% l = legend({'first','last'});
% set(l,'Color', Color(highColor), 'TextColor', Color(baseColor),'EdgeColor',Color(baseColor));

fig = formatFig(fig, color_opt);

% save figure option here...
save_figure(fig, [fileRoot flyID ' solvent time course'], fileTag);
clearvars('-except',initial_vars{:})


%% Fig: Conditioned odors comparions
G1.color = 'teal';  % CSneg color
G2.color = 'orange'; %CSpos color
LW = 1;
color_opt = true;
fileTag = '-png';
alltrials = true;
shading = false;
norm_axes = true;
if color_opt==true
    baseColor = 'white';
else 
    baseColor = 'black';
end

nrows = 2;
ncols = 2;
% color shades for alltrials
greyscale = Color(baseColor,'grey',num.baseTrials);

fig = getfig('Odors Responses',1);
% Naive CS-
subplot(nrows,ncols,1)
    hold on
    y = CSneg.naive;
    if alltrials==true
        for n = 1:num.baseTrials
        plot(time, y(:,n), 'color', greyscale(n,:), 'linewidth', 0.5)
        end
    end
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G1.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
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
    title({['CS- naive : ' param.csNEG_odor];''})   
    
% Test CS-
subplot(nrows,ncols,3)
    hold on
    y = CSneg.test;
    if alltrials==true
        for n = 1:num.testTrials
        plot(time, y(:,n), 'color', greyscale(n,:), 'linewidth', 0.5)
        end
    end
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G1.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    end
    plot(time, avg, 'color', Color(G1.color),'linewidth', LW)
    % labels
    if norm_axes==false
        y1 = rangeLine(fig);
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
    end
    ymax(3,:) = ylim;
    xlabel('time (s)')
    ylabel('\DeltaF/F')
    title({['CS- test : ' param.csNEG_odor];''})  
  

% Naive CS+
subplot(nrows,ncols,2)
    hold on
    y = CSpos.naive;
    if alltrials==true
        for n = 1:num.baseTrials
        plot(time, y(:,n), 'color', greyscale(n,:), 'linewidth', 0.5)
        end
    end
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G2.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
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
    title({['CS+ naive : ' param.csPOS_odor];''})   
    
% Test CS-
subplot(nrows,ncols,4)
    hold on
    y = CSpos.test;
    if alltrials==true
        for n = 1:num.testTrials
        plot(time, y(:,n), 'color', greyscale(n,:), 'linewidth', 0.5)
        end
    end
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G2.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    end
    plot(time, avg, 'color', Color(G2.color),'linewidth', LW)
    % labels
    if norm_axes==false
        y1 = rangeLine(fig);
        odorOFF = param.ODOR_dur;
        plot([0, odorOFF], [y1,y1], 'color',Color(baseColor),'linewidth', 3)
    end
    ymax(4,:) = ylim;
    xlabel('time (s)')
    ylabel('\DeltaF/F')
    title({['CS+ test : ' param.csPOS_odor];''})  
     
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
save_figure(fig, [fileRoot flyID ' all CS trials'], fileTag);

clearvars('-except',initial_vars{:})


%% Fig: naive->test overlaid w/ scatter avg points
G1.color = 'teal';  % CSneg color
G2.color = 'orange'; %CSpos color
LW = 1;
color_opt = true;
fileTag = '-png';
shading = true;
norm_axes = true;
if color_opt==true
    baseColor = 'white';
else 
    baseColor = 'black';
end    

nrows = 1;
ncols = 2;    
    
fig = getfig('Odors Responses',1);
% CS-
subplot(nrows,ncols,1)
    hold on
    % Naive
    y = CSneg.naive;
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(baseColor), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    end
    plot(time, avg+err, 'color', Color(baseColor),'linewidth', 0.5)
    plot(time, avg-err, 'color', Color(baseColor),'linewidth', 0.5)
    plot(time, avg, 'color', Color(baseColor),'linewidth', LW)
    % Test
    y = CSneg.test;
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G1.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    end
    plot(time, avg+err, 'color', Color(G1.color),'linewidth', 0.5)
    plot(time, avg-err, 'color', Color(G1.color),'linewidth', 0.5)
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
    title({['CS- : ' param.csNEG_odor];''})  
    
subplot(nrows,ncols,2)
hold on
    % Naive
    y = CSpos.naive;
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(baseColor), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    end
    plot(time, avg+err, 'color', Color(baseColor),'linewidth', 0.5)
    plot(time, avg-err, 'color', Color(baseColor),'linewidth', 0.5)
    plot(time, avg, 'color', Color(baseColor),'linewidth', LW)
    % Test
    y = CSpos.test;
    err = std(y,0,2);
    avg = mean(y,2);
    if shading==true
        fill_data = error_fill(time, avg, err);
        h = fill(fill_data.X, fill_data.Y, get_color(G2.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    end
    plot(time, avg+err, 'color', Color(G2.color),'linewidth', 0.5)
    plot(time, avg-err, 'color', Color(G2.color),'linewidth', 0.5)
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
    title({['CS+ : ' param.csPOS_odor];''})
     
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
save_figure(fig, [fileRoot flyID ' overlaid CS pre post'], fileTag);

clearvars('-except',initial_vars{:})

    
%% Fig: Scatter plot peak df_f
fileTag = '-png';               %image save type
stats_text = true;              %line with p-vals and sig stars
color_opt = true;               %black vs white background
LW = 1;                         %line width
SZ = 50;                        %scatter point size
l_style = '--';                  %line style for connecting line
paired = true;                  %connecting line
G1.color = Color('purple');     %Solvent
G2.color = Color('teal');       %CS-
G3.color = Color('orange');     %CS+

if color_opt==true
    baseColor = 'white';
else
    baseColor = 'black';
end

% find the max df_f
G1.naive = max(solvent.naive);
G1.test = max(solvent.test);

G3.naive = max(CSpos.naive);
G3.test = max(CSpos.test);    

G2.naive = max(CSneg.naive);
G2.test = max(CSneg.test);    

% plot the peak
fig = getfig('Avg df_f',1);
hold on
% Solvent
len = length(G1.naive);
x1 = 1*ones(1,len);
x2 = 2*ones(1,len);
scatter([x1,x2],[G1.naive,G1.test],SZ,G1.color,'filled')
% CS-
len = length(G2.naive);
x1 = 3*ones(1,len);
x2 = 4*ones(1,len);
scatter([x1,x2],[G2.naive,G2.test],SZ,G2.color,'filled')
% CS-
len = length(G3.naive);
x1 = 5*ones(1,len);
x2 = 6*ones(1,len);
scatter([x1,x2],[G3.naive,G3.test],SZ,G3.color,'filled')
% connecting lines
if paired==true
    plot([1,2], [G1.naive;G1.test],'color', G1.color,'linewidth',LW,'linestyle',l_style)
    plot([3,4], [G2.naive;G2.test],'color', G2.color,'linewidth',LW,'linestyle',l_style)
    plot([5,6], [G3.naive;G3.test],'color', G3.color,'linewidth',LW,'linestyle',l_style)
end

% quick stats:
[~,G2.p] = ttest(G2.naive, G2.test);
[~,G3.p] = ttest(G3.naive, G3.test);
p_val = [G2.p,G3.p];
odors = {param.csNEG_odor,param.csPOS_odor};
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
ylabel('\DeltaF/F')
title({param.cross;'Naive vs Test'})
xlim([0,7])
xlabel('Odor')
ax = gca;
set(ax, 'XTick', [1.5,3.5,5.5])
set(ax, 'XTickLabels', {param.solvent,param.csNEG_odor,param.csPOS_odor})

fig = formatFig(fig,color_opt);

% Plot option
if stats_text == true 
    ax = gca;
    ub = ax.YTick;
    rng = range([ub(1),ub(end)])*0.03; %significance star offset
    pk = ub(end)+mean(diff(ub))+rng;
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

save_figure(fig, [fileRoot flyID ' peak df_f'], fileTag);
clearvars('-except',initial_vars{:})












































