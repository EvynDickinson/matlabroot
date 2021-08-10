clear
close all
clc

%% Headless fly leg movement analysis
% file_root = '/Volumes/Evyn SSD/Evyn UW work/matlabroot';
file_root = 'C:\matlabroot';
fig_dir = [file_root '/Interneuron Lines/Figures/'];

% name of data set
answer = questdlg('Load which data?', '', 'Onball', 'Offball', 'Cancel');
switch answer
    case 'Onball'
        file_name = '13B-20847-csChrimson-7V-onball';
        num.fly = 7;
    case 'Offball'
        file_name = '13B-20847-csChrimson-7V-offball';
        num.fly = 4;
end; clear answer

% load tracking data
load([file_root, '/', file_name ' tracking data.mat'])
% load response data
load([file_root, '/', file_name ' T1 response data.mat'])
% behavior labels
load([file_root, '/', file_name ' behavior class.mat'])

 
behavior = 'Stationary';
Joints = {'CoFe', 'FeTi', 'TiTa'};
condlist = [1:15, 22];
controlConds = [1, 8, 15, 22];


% Convert starting behaviors into numbers and then a filter:
for kk = 1:num.fly
    for cond = 1:28
        for rep = 1:3
            state = group(kk).behavior{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            switch state
                case {'stationary', 'Stationary'}
                    STATE = 1;
                case {'walking', 'Walking'}
                    STATE = 2;
                case {'grooming', 'Grooming'}
                    STATE = 3;
                case {'other', 'Other'}
                    STATE = 4;
                case {'-'}
                    STATE = 5;
            end
            group(kk).STATE(cond, rep) = STATE;
        end
    end
    group(kk).walking = (group(kk).STATE==2);
    group(kk).stationary = (group(kk).STATE==1);
    group(kk).grooming = (group(kk).STATE==3);
    group(kk).other = (group(kk).STATE==4);
end

type = {'stim', 'control'};
rep = 1;
behavior = 'Stationary';
Joints = {'CoFe', 'FeTi', 'TiTa'};
condlist = [1:15, 22];
controlConds = [1, 8, 15, 22];
%Build empty data structure & fill with joint angle data;
for ifly = 1:num.fly
    % build empty structures to fill later with joint positions
    for iJ = 1:3
        data(ifly).(Joints{iJ}).all = NaN(337,28); 
        data(ifly).(Joints{iJ}).raw = NaN(337,28);
    end
    for cond = condlist
        % adjust the first joint position
         idx = 0;
         for iFrame = 90:150 % tracking range
             temp = tracking(ifly).Left_Front(cond,rep).frame(iFrame).pos;
             if ~isempty(temp)  % control period
                idx = idx+1;
                fixPoint(idx,:) = temp(1,:);
             end
             temp = [];
         end
         temp = mean(fixPoint,1);
         Angles = NaN(426,length(Joints));
         % find the Joint positions for all the frames
        for iFrame = 90:426 % tracking range
            input(iFrame).pos = tracking(ifly).Left_Front(cond,rep).frame(iFrame).pos;
            if ~isempty(input(iFrame).pos) % add the fixed position point
                input(iFrame).pos(1,:) = temp;
                Angles(iFrame,:) = calculateAngles(input(iFrame).pos);
            end
        end
        % fill in the gaps with interpolated data:
        xstart = 1:337;
        for iJ = 1:length(Joints)
            vstart = Angles(90:426,iJ);
            nanloc = isnan(vstart);
            xq = find(nanloc==true);
            v = vstart(~nanloc);
            x = xstart(~nanloc);
            vq = interp1(x,v,xq); 
            allV = vstart; 
            allV(nanloc) = vq;
            data(ifly).(Joints{iJ}).all(:,cond) = allV;
            data(ifly).(Joints{iJ}).raw(:,cond) = vstart;
        end       
    end
end


% find the averages for the stationary flies:
for ifly = 1:num.fly
  for iJ = 1:length(Joints)
    % Stim
    filter = group(ifly).stationary(:,rep)';
    filter(controlConds) = false;
    input = data(ifly).(Joints{iJ}).all;
    temp = input(:,filter);
    filter = isnan(temp(1,:));
    temp(:,filter) = [];
    data(ifly).(Joints{iJ}).stim.all = temp;
    % Control
    filter = (logical(1:28)); filter(:) = false;
    filter(controlConds) = true;
    input = data(ifly).(Joints{iJ}).all;
    temp = input(:,filter);
    filter = group(ifly).stationary(controlConds,rep)';
    temp(:,~filter) = [];
    data(ifly).(Joints{iJ}).control.all = temp;
    
    % % Find the average across all trials for a given fly% %
    % find the average
    data(ifly).(Joints{iJ}).stim.avg = nanmean(data(ifly).(Joints{iJ}).stim.all,2);
    data(ifly).(Joints{iJ}).control.avg = nanmean(data(ifly).(Joints{iJ}).control.all,2);
    % find the error
    data(ifly).(Joints{iJ}).stim.err = sem(data(ifly).(Joints{iJ}).stim.all,2,2);
    data(ifly).(Joints{iJ}).control.err = sem(data(ifly).(Joints{iJ}).control.all,2,2);
  end
end

%find change in angle rather than absolute degree
Joints = {'CoFe', 'FeTi', 'TiTa'};
type = {'stim', 'control'};
for iJ = 1:3
  for tt = 1:2
    for ifly = 1:num.fly 
        offset = data(ifly).(Joints{iJ}).(type{tt}).all(1,:);
        data(ifly).(Joints{iJ}).(type{tt}).startangle = offset;
        data(ifly).(Joints{iJ}).(type{tt}).change = data(ifly).(Joints{iJ}).(type{tt}).all-offset;
        data(ifly).(Joints{iJ}).(type{tt}).avgchange = nanmean(data(ifly).(Joints{iJ}).(type{tt}).change,2);
    end
  end
end

% Extract all the joint angles over time -- line/fly
for itype = 1:2
   for iJ = 1:3
      for ifly = 1:num.fly
        inputdata.(Joints{iJ}).(type{itype}).all(:,ifly) = data(ifly).(Joints{iJ}).(type{itype}).avgchange;
      end
      inputdata.(Joints{iJ}).(type{itype}).avg = mean(inputdata.(Joints{iJ}).(type{itype}).all,2);
      inputdata.(Joints{iJ}).(type{itype}).err = sem(inputdata.(Joints{iJ}).(type{itype}).all,2);
   end
end  

clear behavior allV cond filter fixPoint idx ifly iFrame iJ itype 
clear nanloc offset STATE temp vq x xq xstart vstart v tt state kk input

initialVars = who;
initialVars{end+1} = 'initialVars';

%% Find the probability of flexion for a given initial joint angle
% color code with behavior:
stimON = 60; %frame 60 stim turns on
breakpoint = 210; %stimON + (0.5*fps); 
ROI = stimON:breakpoint;
SZ = 200;
fig = getfig('',1);
hold on 

% plot the data
for ifly = 1:num.fly
    % load all trials for a given fly:
    raw = data(ifly).FeTi.all; % totally raw data
    % exclude moving trials & control conditions
    loc = group(ifly).stationary(:,1);
    raw(:,~loc) = nan;
    raw(:,controlConds) = nan;
    % FILTER & PLOT by behavior response:
    
    % Find the max flexion within the first 500ms after stimulus onset:
    initialJoint = raw(stimON,:);
    deltaJoint = raw(stimON,:) - min(raw(ROI,:));
    
    % response locations:
    G1 = strcmpi(ResponseData(ifly).behavior(:,1), 'nothing'); 
    G2 = strcmpi(ResponseData(ifly).behavior(:,1), 'grooming');
    G3 = strcmpi(ResponseData(ifly).behavior(:,1), 'flex');
    G4 = ~(G1|G2|G3);% other:

    % Plot each group:
    scatter(initialJoint(G1), deltaJoint(G1), SZ,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color('grey'));
    scatter(initialJoint(G2), deltaJoint(G2), SZ,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color('YellowGreen'));
    scatter(initialJoint(G3), deltaJoint(G3), SZ,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color('LightPink'));
    scatter(initialJoint(G4), deltaJoint(G4), SZ,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color('Navy'));
    
    %Print the N
    n(ifly) = sum(~isnan(initialJoint));
    fprintf(['\n Fly ' num2str(ifly) ' #trials: ' num2str(n(ifly)) '\n'])
end
fprintf(['\n Total #trials: ' num2str(sum(n)) '\n'])
xlabel('Initial joint angle (\circ)')
ylabel('Max FeTi joint flexion in first 500 ms (\Delta\circ)')
figname = ['Max flexion in ' file_name];
title(figname)
set(gca, 'TickDir', 'out')
xlim([0,180])

save_figure(fig, [fig_dir, figname]);

clearvars('-except',initialVars{:})

%% Plot the probability of flexion for the original binned groups: 
stimON = 60; %frame 60 stim turns on
stimOFF = 276; %stimON + (0.72*fps); 

[y,initialAngle, flexed] = deal([]);
% pull the raw data for each fly from the stim group:
for ifly = 1:num.fly
   raw = data(ifly).FeTi.all; % totally raw data
   % exclude moving trials & control conditions
   loc = group(ifly).stationary(:,1);
   % get trials that flexed:
   temp = double(strcmpi(ResponseData(ifly).behavior(:,1),'flex'));
   temp(~loc) = nan; % use nan for trials that weren't stationary
   temp(controlConds) = nan;
   % angle of the joint at the onset of stimulus
   initialAngle(:,ifly) = raw(stimON,:);
   flexed(:,ifly) = temp;
   
end   
clear raw temp loc

% plot the probability of flexion for binned groups:
Edges = [30, 60, 90, 120, 180];
bins = discretize(initialAngle,Edges);
% find total number in each group:
temp = bins;
temp(isnan(flexed)) = nan;

% print number of trials
n = sum(~isnan(temp),1);
disp(['#trials by fly: ' num2str(n)])
disp(['#trials total: ' num2str(sum(n))])

for nbin = 1:length(Edges)-1
   % find total number of trials in a bin
   tot(nbin) = sum(sum(temp==nbin));
   % find number of trials that flexed
   flexcount(nbin) = sum(sum(temp==nbin & flexed==1));
   % find number of trials that DIDN'T flex
   Noflexcount(nbin) = sum(sum(temp==nbin & flexed==0));
   % pull prob for graph
   y(nbin) = (flexcount(nbin)/tot(nbin))*100;
end

fig = getfig('',1); hold on
x = 1:length(Edges)-1;
bar(x,y)
% scatter(nbin, y, 100, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
xlim([0,nbin+1])
ylim([0,100])
xlabel(['Bin edges: ' num2str(Edges)])
ylabel('Probability of flexion')
figname = [file_name ' prob of flexion'];
title({figname; ['Prob: ' num2str(y)]})
save_figure(fig, [fig_dir, figname]);

clearvars('-except',initialVars{:})

%% Cumulative probability of flexion: 

stimON = 60; %frame 60 stim turns on
stimOFF = 276; %stimON + (0.72*fps); 


[x,y,initialAngle, flexed] = deal([]);
% pull the raw data for each fly from the stim group:
for ifly = 1:num.fly
   raw = data(ifly).FeTi.all; % totally raw data
   % exclude moving trials & control conditions
   loc = group(ifly).stationary(:,1);
   % get trials that flexed:
   temp = double(strcmpi(ResponseData(ifly).behavior(:,1),'flex'));
   temp(~loc) = nan; % use nan for trials that weren't stationary
   temp(controlConds) = nan;
   % angle of the joint at the onset of stimulus
   initialAngle(:,ifly) = raw(stimON,:);
   flexed(:,ifly) = temp;
end   
clear raw temp loc

% remove instances of non movers / control conditions from initial angle:
temp = initialAngle;
temp(isnan(flexed)) = nan;
angle = temp;
temp = flexed;
temp(isnan(angle)) = nan;
flex = temp;
loc = (isnan(flex) & isnan(angle));
flex(loc) = [];
angle(loc) = [];

% sort and organize the data:
[~,I] = sort(angle);
x = angle(I)';
y = flex(I)';

% find the cumulative probability for each trial (ascending order)
for ii = 1:length(angle)
   z(ii,1) = (sum(y(1:ii))/length(angle))*100;    
end
% add a data point for 0 and 180:
cprob = [0,0; x,z; 180,z(end)];

fig = getfig('',1); hold on
plot(cprob(:,1), cprob(:,2), 'color', 'k', 'linewidth',1)
xlim([0,180])
ylim([0,100])
xlabel('Joint Angle')
ylabel('C-prob of flexion')

figname = [file_name ' C-prob of flexion'];
title(figname)
save_figure(fig, [fig_dir, figname]);

% fit1 = cprob(:,1);
% fit2 = cprob(:,2);

clearvars('-except',initialVars{:})

%% Bar graph of response data 
responselist = [];
for ifly = 1:num.fly
   a = ResponseData(ifly).behavior(:,1);
   responselist = [responselist, a];
end
options = unique(responselist);

labelname = [];
stimConds = [2:7,9:14];
% find proportions of responses
for ii = 1:length(options)
    checklist(ii).logic = strcmpi(options{ii},responselist);
    checklist(ii).name = options{ii};
    checklist(ii).stim = checklist(ii).logic(stimConds,:);
    checklist(ii).control = checklist(ii).logic(controlConds,:);
    % find the proportions:
    raw = checklist(ii).stim;
    a = sum(sum(raw));
    checklist(ii).stimstats = a/(numel(raw))*100;
    plotdata.stim(ii) = checklist(ii).stimstats;
    raw = checklist(ii).control;
    a = sum(sum(raw));
    checklist(ii).controlstats = a/(numel(raw))*100;
    plotdata.control(ii) = checklist(ii).controlstats;
    labelname = [labelname ', ' options{ii}];
end

% Blues:
StimInColor = {'006DDB', '198BFF', '50A7FF', '8AC4FF', 'B6DAFF'};
% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f', 'b2b2b2', 'cccccc'};
% LineStyles
Lstyles = {'-', '--', ':', '-.'};

for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end


%plot the data
fig = getfig; hold all
%control
x = 0.85:1:length(options)-.15;
y = plotdata.control;
bar(x,y, 0.3, 'FaceColor', cntlColors(1,:),'EdgeColor','none')
% stim
x = 1.15:1:length(options)+.15;
y = plotdata.stim;
bar(x,y, 0.3, 'FaceColor', stimColors(2,:),'EdgeColor','none')

% labels
set(gca,'TickDir','out');
ylim([0 100])
xlabel(labelname)
ylabel('% of trials')
figname = [file_name ' T1 responses graph'];
title({figname; ''})

save_figure(fig, [fig_dir, figname]);

clearvars('-except',initialVars{:})

%% Starting joint angles of all included trials
stimON = 60; % loc of start of stim
sz = 50; % size of scatter points
% colors for the diff joints:
% Blue/greens for 13B:
StimInColor = {'00dba0', '006DDB', '00c8db'}; 
% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f'};
for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end

fig = getfig;  
% PULL THE DATA
for iJ = 1:3
    stimdata = []; controldata = [];
    for ifly = 1:num.fly
        a = data(ifly).(Joints{iJ}).stim.all(stimON,:)';
        stimdata = [stimdata; a];
        b = data(ifly).(Joints{iJ}).control.all(stimON,:)';
        controldata = [controldata; b];
    end
    subplot(3,1,iJ)
    hold all
    h = histogram(stimdata,0:10:180);
    xlim([0,180])
    ylim([0, 20])
    title(Joints{iJ})
%     ylim([0,50])
%     
%     figure(Cfig)
%     subplot(3,1,iJ)
%     hold all
%     h = histogram(controldata,18);
%     xlim([0,180])
%     ylim([0,50])
end
subplot(3,1,1)
figname = [file_name, ' start joint angle histo'];
title({figname; 'CoFe'})

save_figure(fig, [fig_dir, figname]);

%     % PLOT stim data
%     n = length(stimdata);
%     x = 2.5:1/(n-1):3.5;  
%     y = stimdata(randi([1,n],1,n));
%     scatter(x, y, sz, stimColors(iJ,:), 'filled')
%     % PLOT control data
%     n = length(controldata);
%     x = 2.5:1/(n-1):3.5;  
%     y = controldata(randi([1,n],1,n));
%     x = 1:1/(length(controldata)-1):2;    
%     scatter(x, y, sz, cntlColors(iJ,:), 'filled')
% end
% LABEL
% xlim([0 6])
% ylim([0 180])
% ylabel('Joint angles')
% xlabel('Control vs Stim')
% figname = [file_name ' Start Joint Angles'];
% set(gca,'TickDir','out');
% title({figname; ''})
% 
% save_figure(fig, [fig_dir, figname])

clearvars('-except',initialVars{:})

%% figure of T2 movement vs no movement trials in joint angles

% Find winshield responses
for kk = 1:num.fly
    for cond = 1:28
        for rep = 1:3
            STATE = [];
            state = ResponseData(kk).behavior{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            switch state
                case {'windshield', 'Windshield'}
                    STATE = 1;
                case {'leg_lift', 'leglift'}
                    STATE = 2;
                case {'stepping', 'Stepping'}
                    STATE = 3;
                case {'other', 'Other'}
                    STATE = 4;
                case {'-', 'Nothing'}
                    STATE = 5;
            end
            if isempty(STATE)
                STATE = 6;
            end
            ResponseData(kk).STATE(cond, rep) = STATE;
        end
    end
    ResponseData(kk).windshield = (ResponseData(kk).STATE==1);
    ResponseData(kk).leglift = (ResponseData(kk).STATE==2);
    ResponseData(kk).stepping = (ResponseData(kk).STATE==3);
    ResponseData(kk).nothing = (ResponseData(kk).STATE==5);
end

clear STATE state kk 
rep = 1;
stimlist = [2:7,9:14];
% Find proportion of windshield 
for ifly = 1:num.fly
    temp = ResponseData(ifly).windshield(:,rep); %windshield only
    filter = group(ifly).stationary(:,rep); %stationary only
    loc = all([temp, filter],2); %both stationary and windshield
    
    % STIM probability / locations
    total = sum(filter(stimlist)); % total stationary
    WS = sum(loc(stimlist));
    inputdata.response.stim.total(ifly) = total;
    inputdata.response.stim.WS(ifly) = WS;
    inputdata.response.stim.percentWS(ifly) = (WS/total)*100;
    
    MT = logical(1:28); MT(controlConds) = false; % stim locations
    a = all([MT',loc],2);% stim locations that also did the windshield
    inputdata.response.stim.wsLoc(:,ifly) = a;
    MT = logical(1:28); MT(:) = false; MT(stimlist) = true;
    a = all([MT',filter,~temp],2); % stim locations that did nothing
    inputdata.response.stim.nonLoc(:,ifly) = a; 
    
    % CONTROL probability / locations
    total = sum(filter(controlConds)); % total stationary
    WS = sum(loc(controlConds));
    inputdata.response.control.total(ifly) = total;
    inputdata.response.control.WS(ifly) = WS;
    inputdata.response.control.percentWS(ifly) = (WS/total)*100;

    MT = logical(1:28); MT(:) = false; MT(controlConds)=true; % control locations
    a = all([MT',loc],2);% control locations that also did the windshield
    inputdata.response.control.wsLoc(:,ifly) = a;
    a = all([MT',filter,~temp],2); % control locations that did nothing
    inputdata.response.control.nonLoc(:,ifly) = a; 
end

% Blue/greens for 13B:
StimInColor = {'00dba0', '006DDB', '00c8db'}; %color per joint
% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f'};
for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end


% PLOT the absolute value for all three joints
fig = getfig; 
x = 0:1/300:1.12; sSpan = 3;
for iJ = 1:3 %1:length(Joints)
    subplot(1,3,iJ); hold all
    plotdata.non = []; plotdata.WS = [];
    for ifly = 1:num.fly
       % non WS cases
       filter = inputdata.response.stim.nonLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
           plotdata.non = [plotdata.non, b];
       end
       
       % WS cases
       filter = inputdata.response.stim.wsLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
            plotdata.WS = [plotdata.WS, b];
       end
    end
    
    % AVERAGES & ERR
    % no WS
    avg = smooth(mean(plotdata.non,2),5);
    err = smooth(sem(plotdata.non,2),5);
    plot(x, avg, 'color', stimColors(iJ,:), 'linewidth', 2)
    plot(x, avg+err, 'color', stimColors(iJ,:), 'linewidth', 1)
    plot(x, avg-err, 'color', stimColors(iJ,:), 'linewidth', 1)

    % WS
    avg = smooth(nanmean(plotdata.WS,2),5);
    err = smooth(sem(plotdata.WS,2),5);
    plot(x, avg, 'color', Color('lightpink'), 'linewidth', 2)
    plot(x, avg+err, 'color', Color('lightpink'), 'linewidth', 1)
    plot(x, avg-err, 'color', Color('lightpink'), 'linewidth', 1)
    
    % plot light lines / add labels
    plot([0.2, 0.92], [100, 100], 'linewidth', 1, 'color', 'g')
    ylim([0 180])
    xlabel('time (s)')
    ylabel([Joints{iJ} ' angle (deg)'])
    set(gca,'TickDir','out');
end
subplot(1,3,2)
figname = [file_name ' headless WS vs no WS by joint'];
title({figname; 'pinks = WS, blue = noWS'}) %(1:end-18)

save_figure(fig, [fig_dir, figname]);


% PLOT the CHANGE in angle for all three joints
fig = getfig; hold all
Pinks = Color('purple', 'lightpink', 3);
x = 0:1/300:1.12; sSpan = 3;
for iJ = 1:3 %1:length(Joints)
    plotdata.non = []; plotdata.WS = [];
    for ifly = 1:num.fly
       % non WS cases
       filter = inputdata.response.stim.nonLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
           plotdata.non = [plotdata.non, b];
       end
       
       % WS cases
       filter = inputdata.response.stim.wsLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
            plotdata.WS = [plotdata.WS, b];
       end
    end
    
    % AVERAGES & ERR
    % no WS
    offset = plotdata.non(60,:);
    input = plotdata.non-offset;
    avg = smooth(nanmean(input,2),5);
    err = smooth(sem(input,2),5);
%     plot(x, avg, 'color', stimColors(iJ,:), 'linewidth', 2)
%     plot(x, avg+err, 'color', stimColors(iJ,:), 'linewidth', 1)
%     plot(x, avg-err, 'color', stimColors(iJ,:), 'linewidth', 1)
    
    fill_data = error_fill(x, avg,err);
    h = fill(fill_data.X, fill_data.Y, stimColors(iJ,:), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(x, avg, 'color', stimColors(iJ,:), 'linewidth', 3)

    % WS
    offset = plotdata.WS(60,:);   
    input = plotdata.WS-offset;
    avg = smooth(nanmean(input,2),5);
    err = smooth(sem(input,2),5);
%     plot(x, avg, 'color', Pinks(iJ,:), 'linewidth', 2)
%     plot(x, avg+err, 'color', Pinks(iJ,:), 'linewidth', 1)
%     plot(x, avg-err, 'color', Pinks(iJ,:), 'linewidth', 1)
    
    fill_data = error_fill(x, avg,err);
    h = fill(fill_data.X, fill_data.Y, Pinks(iJ,:), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(x, avg, 'color', Pinks(iJ,:), 'linewidth', 3)
    
    % plot light lines / add labels
    plot([0.2, 0.92], [30, 30], 'linewidth', 1, 'color', 'g')
    hline(0,'k:')
end

set(gca,'TickDir','out');
xlabel('Time (sec)')
ylabel('Change in joint angle (deg)')

title({[file_name ' headless']; 'pinks = WS, blue = noWS'}) %(1:end-18)

save_figure(fig, [fig_dir, file_name ' change in joint angles WS vs no windshield']);


clearvars('-except',initialVars{:})

%% Joint angles time course for flexed vs no flex in 13B neuron activation
rep = 1;
% Find winshield responses
for kk = 1:num.fly
    for cond = 1:28
            STATE = [];
            state = ResponseData(kk).behavior{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            switch state
                case {'flex', 'Flex'}
                    STATE = 1;
                case {'nothing', 'Nothing'}
                    STATE = 2;
                case {'step', 'Step'}
                    STATE = 3;
                case {'other', 'Other','-'}
                    STATE = 4;
            end
            if isempty(STATE)
                STATE = 6;
            end
            ResponseData(kk).STATE(cond, rep) = STATE;

    end
    ResponseData(kk).flex = (ResponseData(kk).STATE==1);
    ResponseData(kk).nothing = (ResponseData(kk).STATE==2);
    ResponseData(kk).step = (ResponseData(kk).STATE==3);
    ResponseData(kk).other = (ResponseData(kk).STATE==5);
end

clear STATE state kk 
stimlist = [2:7,9:14];
% Find proportion of windshield 
for ifly = 1:num.fly
    temp = ResponseData(ifly).flex(:,rep); %windshield only
    filter = group(ifly).stationary(:,rep); %stationary only
    loc = all([temp, filter],2); %both stationary and windshield
    
    % STIM probability / locations
    total = sum(filter(stimlist)); % total stationary
    WS = sum(loc(stimlist));
    inputdata.response.stim.total(ifly) = total;
    inputdata.response.stim.WS(ifly) = WS;
    inputdata.response.stim.percentWS(ifly) = (WS/total)*100;
    
    MT = logical(1:28); MT(controlConds) = false; % stim locations
    a = all([MT',loc],2);% stim locations that also did the windshield
    inputdata.response.stim.wsLoc(:,ifly) = a;
    MT = logical(1:28); MT(:) = false; MT(stimlist) = true;
    a = all([MT',filter,~temp],2); % stim locations that did nothing
    inputdata.response.stim.nonLoc(:,ifly) = a; 
    
    % CONTROL probability / locations
    total = sum(filter(controlConds)); % total stationary
    WS = sum(loc(controlConds));
    inputdata.response.control.total(ifly) = total;
    inputdata.response.control.WS(ifly) = WS;
    inputdata.response.control.percentWS(ifly) = (WS/total)*100;

    MT = logical(1:28); MT(:) = false; MT(controlConds)=true; % control locations
    a = all([MT',loc],2);% control locations that also did the windshield
    inputdata.response.control.wsLoc(:,ifly) = a;
    a = all([MT',filter,~temp],2); % control locations that did nothing
    inputdata.response.control.nonLoc(:,ifly) = a; 
end

% Blue/greens for 13B:
StimInColor = {'00dba0', '006DDB', '00c8db'}; %color per joint
% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f'};
for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end


% PLOT the absolute value for all three joints
fig = getfig; 
x = 0:1/300:1.12; sSpan = 3;
for iJ = 1:3 %1:length(Joints)
    subplot(1,3,iJ); hold all
    plotdata.non = []; plotdata.WS = [];
    for ifly = 1:num.fly
       % non WS cases
       filter = inputdata.response.stim.nonLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
           plotdata.non = [plotdata.non, b];
       end
       
       % WS cases
       filter = inputdata.response.stim.wsLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
            plotdata.WS = [plotdata.WS, b];
       end
    end
    
    % AVERAGES & ERR
    % no WS
    avg = smooth(mean(plotdata.non,2),5);
    err = smooth(sem(plotdata.non,2),5);
    plot(x, avg, 'color', stimColors(iJ,:), 'linewidth', 2)
    plot(x, avg+err, 'color', stimColors(iJ,:), 'linewidth', 1)
    plot(x, avg-err, 'color', stimColors(iJ,:), 'linewidth', 1)

    % WS
    avg = smooth(nanmean(plotdata.WS,2),5);
    err = smooth(sem(plotdata.WS,2),5);
    plot(x, avg, 'color', Color('lightpink'), 'linewidth', 2)
    plot(x, avg+err, 'color', Color('lightpink'), 'linewidth', 1)
    plot(x, avg-err, 'color', Color('lightpink'), 'linewidth', 1)
    
    % plot light lines / add labels
    plot([0.2, 0.92], [100, 100], 'linewidth', 1, 'color', 'g')
    ylim([0 180])
    xlabel('time (s)')
    ylabel([Joints{iJ} ' angle (deg)'])
    set(gca,'TickDir','out');
end
subplot(1,3,2)
figname = [file_name ' flexed vs no repsonse'];
title({figname; 'pinks = flex, blue = no response'}) %(1:end-18)

save_figure(fig, [fig_dir, figname]);


% PLOT the CHANGE in angle for all three joints
fig = getfig; hold all
Pinks = Color('purple', 'lightpink', 3);
x = 0:1/300:1.12; sSpan = 3;
for iJ = 1:3 %1:length(Joints)
    plotdata.non = []; plotdata.WS = [];
    for ifly = 1:num.fly
       % non WS cases
       filter = inputdata.response.stim.nonLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
           plotdata.non = [plotdata.non, b];
       end
       
       % WS cases
       filter = inputdata.response.stim.wsLoc(:,ifly);
       b = data(ifly).(Joints{iJ}).all(:,filter);
       if ~isempty(b)
            plotdata.WS = [plotdata.WS, b];
       end
    end
    
    % AVERAGES & ERR
    % no WS
    offset = plotdata.non(60,:);
    input = plotdata.non-offset;
    avg = smooth(nanmean(input,2),5);
    err = smooth(sem(input,2),5);
    plot(x, avg, 'color', stimColors(iJ,:), 'linewidth', 2)
    plot(x, avg+err, 'color', stimColors(iJ,:), 'linewidth', 1)
    plot(x, avg-err, 'color', stimColors(iJ,:), 'linewidth', 1)
    
%     fill_data = error_fill(x, avg,err);
%     h = fill(fill_data.X, fill_data.Y, stimColors(iJ,:), 'EdgeColor','none');
%     set(h, 'facealpha', 0.2)
%     plot(x, avg, 'color', stimColors(iJ,:), 'linewidth', 3)

    % WS
    offset = plotdata.WS(60,:);   
    input = plotdata.WS-offset;
    avg = smooth(nanmean(input,2),5);
    err = smooth(sem(input,2),5);
    plot(x, avg, 'color', Pinks(iJ,:), 'linewidth', 2)
    plot(x, avg+err, 'color', Pinks(iJ,:), 'linewidth', 1)
    plot(x, avg-err, 'color', Pinks(iJ,:), 'linewidth', 1)
    
%     fill_data = error_fill(x, avg,err);
%     h = fill(fill_data.X, fill_data.Y, Pinks(iJ,:), 'EdgeColor','none');
%     set(h, 'facealpha', 0.2)
%     plot(x, avg, 'color', Pinks(iJ,:), 'linewidth', 3)
%     
    % plot light lines / add labels
    plot([0.2, 0.92], [30, 30], 'linewidth', 1, 'color', 'g')
    hline(0,'k:')
end
ylim([-60, 30])
set(gca,'TickDir','out');
xlabel('Time (sec)')
ylabel('Change in joint angle (deg)')

title({[file_name ' headless']; 'pinks = flex, blue = no response'}) %(1:end-18)

save_figure(fig, [fig_dir, file_name ' change in joint angles flex vs nothing']);

clear avg cond condlist err filter h input loc MT offset Pinks 
clear sSpan StimInColor stimlist temp total WS x a b fig figname


% FeTi joint angle by inital joint angle
% Edges = [40 80 100 120 180]; nbins = 4;
Edges = [20 60 90 120 180]; nbins = 4;
% make a struct with all the trials together
temp = struct('control', [], 'stim', [], 'raw_index', [], 'flex_index', []);
for ifly = 1:num.fly
    input = data(ifly).FeTi.all;
    filter = group(ifly).stationary(:,rep); % only include stationary flies in the analysis
    input(:,~filter) = NaN;
    % find the portion of flex within each range group...
    b = discretize(input(1,:), Edges);
    b(controlConds) = NaN;
    
    filter = ResponseData(ifly).flex; % filter out everything but a flex
    a = input(:,controlConds); % isolate the control conditions *no response filter*
    input(:,~filter) = NaN;
    input(:,controlConds) = NaN; % filter out the control conditions

    temp.stim = [temp.stim, input];
    temp.control = [temp.control, a];
    temp.raw_index = [temp.raw_index, b];
    b(~filter) = NaN;
    temp.flex_index = [temp.flex_index, b];
end; clear a b filter

inputdata.ALL.FeTi = temp;

% find the relative portions of flex vs no flex per joint angle:
for ii = 1:nbins
    flex.raw(ii) = length(find(temp.raw_index == ii));
    flex.flexed(ii) = length(find(temp.flex_index == ii));
end
flex.percent = (flex.flexed./flex.raw)*100;





% bin the joint angle data
[Y,~] = discretize(inputdata.ALL.FeTi.stim(1,:),Edges); %stim
Ycontrol = discretize(inputdata.ALL.FeTi.control(1,:),Edges); %control
for ii = 1:nbins
    % -- stim 
      loc = []; loc = (Y==ii); 
      if sum(loc)>1
        inputdata.ALL.bin(ii).data = inputdata.ALL.FeTi.stim(:,loc);
        inputdata.ALL.bin(ii).avg = nanmean(inputdata.ALL.bin(ii).data,2);
        inputdata.ALL.bin(ii).err = sem(inputdata.ALL.bin(ii).data,2);
        % change in angle over time
        inputdata.ALL.bin(ii).offset = inputdata.ALL.bin(ii).avg(1);
        inputdata.ALL.bin(ii).change = inputdata.ALL.bin(ii).avg-inputdata.ALL.bin(ii).offset;
      else
        inputdata.ALL.bin(ii).change = NaN(1,337);
      end
      
        % -- control
      loc = []; loc = (Ycontrol==ii);
      if sum(loc)>2  
        inputdata.ALL.cntl(ii).data = inputdata.ALL.FeTi.control(:,loc);
        inputdata.ALL.cntl(ii).avg = nanmean(inputdata.ALL.cntl(ii).data,2);
        inputdata.ALL.cntl(ii).err = sem(inputdata.ALL.cntl(ii).data,2);
        % change in angle over time
        inputdata.ALL.cntl(ii).offset = inputdata.ALL.cntl(ii).avg(1);
        inputdata.ALL.cntl(ii).change = inputdata.ALL.cntl(ii).avg-inputdata.ALL.cntl(ii).offset;
      else
        inputdata.ALL.cntl(ii).change = NaN(1,337);
      end
end



% colorList = {'Navy', 'Blue', 'DodgerBlue', 'LightSkyBlue',...
%              'Black', 'DimGrey', 'DarkGrey', 'LightGrey'};
% Blues:
StimInColor = {'006DDB', '198BFF', '50A7FF', '8AC4FF', 'B6DAFF'};
% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f', 'b2b2b2', 'cccccc'};
% LineStyles
Lstyles = {'-', '--', ':', '-.'};

for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end

% average change in joint angle over time -- color = initial angle

x = 0:1/300:1.12;
fig = getfig;
hold all
for ii = 1:nbins
    avg = inputdata.ALL.bin(ii).change;
    if sum(avg)>=0 || sum(avg) <0
        err = inputdata.ALL.bin(ii).err;
        plot(x, avg, 'color', stimColors(ii,:),...
                'linewidth', 1, 'linestyle', Lstyles{ii})
        plot(x, avg+err, 'color', stimColors(ii,:),...
                'linewidth', 1, 'linestyle', Lstyles{ii}) 
        plot(x, avg-err, 'color', stimColors(ii,:),...
            'linewidth', 1, 'linestyle', Lstyles{ii}) 
    else 
        disp(ii)
    end
end
set(gca,'TickDir','out');
hline(0, 'k'); xlim([0 1.2]); ylim([-80 40])
% vline([0.2, 0.92],'g-')
plot([0.2, 0.92], [8, 8], 'linewidth', 3, 'color', 'g')
plot([0.2, 0.92], [10 10], 'linewidth', 2, 'color', stimColors(1,:))
plot([0.2, 0.92], [12 12], 'linewidth', 2, 'color', stimColors(2,:))
plot([0.2, 0.92], [14 14], 'linewidth', 2, 'color', stimColors(3,:))
plot([0.2, 0.92], [16 16], 'linewidth', 2, 'color', stimColors(4,:))
xlabel('time (sec)')
ylabel('Change in FeTi joint angle (deg)')
title(['Femur-Tibia joint angles during ' file_name ' activation'])

textnote = ['percent that flexed: ' num2str(flex.percent)];
annotation('textbox', [0.2, 0.15, 0.1, 0.1], 'String', textnote)
textnote = ['N flexed: ' num2str(flex.flexed)];
annotation('textbox', [0.2, 0.1, 0.1, 0.1], 'String', textnote)
textnote = ['N total: ' num2str(flex.raw)];
annotation('textbox', [0.2, 0.05, 0.1, 0.1], 'String', textnote)

figname = [fig_dir file_name ' FeTi change in joint angle FLEX only'];
saveas(fig, figname);
save_figure(fig, figname);

clearvars('-except',initialVars{:})

%% Plot the average joint angle over time for the different joints
% Blue/greens for 13B:
StimInColor = {'00dba0', '006DDB', '00c8db'}; 
% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f'};
for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end

% PLOT ALL THREE JOINT ANGLES OVER TIME
x = 0:1/300:1.12;
fig = getfig; hold all
for iJ = 1:3
    %CONTROL
%     fill_data = error_fill(x, inputdata.(Joints{iJ}).control.avg,...
%                            inputdata.(Joints{iJ}).control.err);
%     h = fill(fill_data.X, fill_data.Y, cntlColors(iJ,:), 'EdgeColor','none');
%     set(h, 'facealpha', 0.2)
    avg = inputdata.(Joints{iJ}).control.avg;
    err = inputdata.(Joints{iJ}).control.err;
    plot(x, avg, 'color', cntlColors(iJ,:), 'linewidth', 3)
    plot(x, avg+err, 'color', cntlColors(iJ,:), 'linewidth', 1)
    plot(x, avg-err, 'color', cntlColors(iJ,:), 'linewidth', 1)
end
for iJ = 1:3
    %STIM
%     fill_data = error_fill(x, inputdata.(Joints{iJ}).stim.avg,...
%                            inputdata.(Joints{iJ}).stim.err);
%     h = fill(fill_data.X, fill_data.Y, stimColors(iJ,:), 'EdgeColor','none');
%     set(h, 'facealpha', 0.2)
    avg = inputdata.(Joints{iJ}).stim.avg;
    err = inputdata.(Joints{iJ}).stim.err;
    plot(x, avg, 'color', stimColors(iJ,:), 'linewidth', 3)
    plot(x, avg+err, 'color', stimColors(iJ,:), 'linewidth', 1)
    plot(x, avg-err, 'color', stimColors(iJ,:), 'linewidth', 1)
end
plot([0.2, 0.92], [30, 30], 'linewidth', 1, 'color', 'g')
% hline(0, 'k'); xlim([0 1.2]); 
ylim([-50, 30])
set(gca,'TickDir','out');
xlabel('Time (sec)')
ylabel('Joint angle (deg)')

title({[file_name ' headless']; ''}) %(1:end-18)

save_figure(fig, [fig_dir, file_name ' avg joint angles']);
clearvars('-except',initialVars{:})

%% Deep dive into FeTi joint angle over time %% 
% make a struct with all the trials together
shading_opt = true;
fileTag = '-png';

temp = struct('control', [], 'stim', []);
for ifly = 1:length(data)
    input = data(ifly).FeTi.all;
    a = input(:,1);
    b = input(:,8);
    input(:,[1 8]) = []; % delete the two control traces

    temp.control = [temp.control, a, b];
    temp.stim = [temp.stim, input];
end
inputdata.ALL.FeTi = temp;
%    % filter out flies that are moving before the stimulus starts
% if strcmpi(file_name, '13B-20847-csChrimson-7V-offball tracking data.mat')
%     errlist = range(inputdata.ALL.FeTi.stim(1:8,:),1);
%     figure; hold all; scatter(1:length(errlist), errlist)
%     loc = (errlist>=6);
%     inputdata.ALL.FeTi.stim(:,loc) = [];
%     % remove the fly that was grooming:
%     inputdata.ALL.FeTi.control(:,1) = [];
% end

% bin the joint angle data

E = [30 60 90 120 180]; nbins = 4;
Y = discretize(inputdata.ALL.FeTi.stim(1,:),E); %stim
Ycontrol = discretize(inputdata.ALL.FeTi.control(1,:),E); %control
for ii = 1:nbins
    % -- stim
      loc = []; loc = (Y==ii); 
      if sum(loc)>2
        inputdata.ALL.bin(ii).data = inputdata.ALL.FeTi.stim(:,loc);
        inputdata.ALL.bin(ii).avg = nanmean(inputdata.ALL.bin(ii).data,2);
        inputdata.ALL.bin(ii).err = sem(inputdata.ALL.bin(ii).data,2);
        % change in angle over time
        inputdata.ALL.bin(ii).offset = inputdata.ALL.bin(ii).avg(1);
        inputdata.ALL.bin(ii).change = inputdata.ALL.bin(ii).avg-inputdata.ALL.bin(ii).offset;
      else
        inputdata.ALL.bin(ii).change = NaN(1,337);
      end
      
        % -- control
      loc = []; loc = (Ycontrol==ii);
      if sum(loc)>2  
        inputdata.ALL.cntl(ii).data = inputdata.ALL.FeTi.control(:,loc);
        inputdata.ALL.cntl(ii).avg = nanmean(inputdata.ALL.cntl(ii).data,2);
        inputdata.ALL.cntl(ii).err = sem(inputdata.ALL.cntl(ii).data,2);
        % change in angle over time
        inputdata.ALL.cntl(ii).offset = inputdata.ALL.cntl(ii).avg(1);
        inputdata.ALL.cntl(ii).change = inputdata.ALL.cntl(ii).avg-inputdata.ALL.cntl(ii).offset;
      else
        inputdata.ALL.cntl(ii).change = NaN(1,337);
      end
end
% 
% colorList = {'Navy', 'Blue', 'DodgerBlue', 'LightSkyBlue',...
%              'Black', 'DimGrey', 'DarkGrey', 'LightGrey'};
% Blues:
% StimInColor = {'006DDB', '198BFF', '50A7FF', '8AC4FF', 'B6DAFF'};
% orange-brown
StimInColor = {'000000', '7C5A3F', 'F7931E', 'FFB677'};

% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f', 'b2b2b2', 'cccccc'};
% LineStyles
Lstyles = {':', '--', '-.', '--'};

for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end

% average change in joint angle over time -- color = initial angle

x = 0:1/300:1.12;
LW = 1;
fig = getfig('',1);
hold all
for ii = 1:nbins
    y = inputdata.ALL.bin(ii).change;
    if sum(~isnan(y))==0; continue; end
    err = inputdata.ALL.bin(ii).err;
    % shading
    if shading_opt==true
        fill_data = error_fill(x, y, err);  
        h = fill(fill_data.X, fill_data.Y, stimColors(ii,:), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    end
    % outlines and main trace
    plot(x, y, 'color', stimColors(ii,:),...
            'linewidth', LW, 'linestyle', Lstyles{ii})
    plot(x, y+err, 'color', stimColors(ii,:),'linewidth', 0.5)
    plot(x, y-err, 'color', stimColors(ii,:),'linewidth', 0.5)  
    
end
set(gca,'TickDir','out');
% hline(0, 'k');
xlim([0 1.2]);
ax = gca;
ax.XTick = [0,0.6,1.2];
% Y axis adjusts:
%offball:
% ylim([-60 20])
% ax.YTick = [-60:20:20];
%onball:
ylim([-30 10])
ax.YTick = -30:10:10;

y1 = rangeLine(fig);
plot([0.2, 0.92], [y1, y1], 'linewidth', 3, 'color', 'g')
xlabel('time (sec)')
ylabel('Change in FeTi joint angle (deg)')
title(['Femur-Tibia joint angles during ' file_name ' activation'])

figname = [fig_dir file_name ' FeTi change in joint angle'];
saveas(fig, figname);
save_figure(fig, figname, fileTag);


% % plotted individually
% figure; hold all
% for ii = 1:nbins
%     plot(xstart, inputdata.ALL.bin(ii).avg, 'color',  get_color(colorList{ii}), 'linewidth', 3)
%     plot(xstart, inputdata.ALL.cntl(ii).avg, 'color', get_color(colorList{ii+4}), 'linewidth', 3)
% end  
% vline([60,317],'g-')
% set(gca,'TickDir','out');
clearvars('-except',initialVars{:})

%% Bootstrapping statistics for the joint flexion data
% ---------------------------------------------------------------------
% 1) calculate the difference in speed from stimulus onset to 200ms post
% stimulus offset (diff in area under the curv) between SH and IN
% ---------------------------------------------------------------------
% compare between each light length (0ms, 90ms, 720ms) in the split half
% and the interneuron (eg 10B, 13B or 9A)

% ---------------------------------------------------------------------
% 2) Mix and randomly pull (10,000 times) samples of the trials (with
% replacement) and assign to control | stim and then find the diff in speed
% ---------------------------------------------------------------------
% CHANGE THIS ALL TO BE ABSOLUTE CHANGE IN JOINT ANGLE -- NOT THE ABSOLUTE
% JOINT ANGLE 

% Blue/greens for 13B:
StimInColor = {'00dba0', 'f37b4f', '00c8db'}; 
% Greyscale
ControlInColor = {'000000', '4c4c4c', '7f7f7f'};
for ii = 1:length(StimInColor)
    stimColors(ii,:) =  hex2rgb(StimInColor{ii});     
    cntlColors(ii,:) =  hex2rgb(ControlInColor{ii});   
end

stats = [];
cConds = [1,8,15,22];
sConds = [2:7, 9:14];
sROI = 217:276; % 200ms period at the end of the stimulation period
cROI = 1:60; % 200ms pre stimstart

rep = 1;
for iJ = 1:3
    temp = struct('control', [], 'stim', []);
    for ifly = 1:length(data)
        raw = data(ifly).(Joints{iJ}).all;
        %eliminate nonstationary trials
        filter = group(ifly).stationary(:,rep);
        raw(:,~filter) = NaN;
        a = raw(:,cConds);
        b = raw(:,sConds);

        temp.control = [temp.control, a];
        temp.stim = [temp.stim, b];
    end
    % remove the NaN trials
    loc = sum(isnan(temp.control),1) > 0; % nan anywhere in the list
    temp.control(:,loc) = [];
    loc = sum(isnan(temp.stim),1) > 0; % nan anywhere in the list
    temp.stim(:,loc) = [];
    % adjust everything to be normalized change in joint angle:
    offset = temp.control(1,:);
    temp.control = temp.control-offset;
    offset = temp.stim(1,:);
    temp.stim = temp.stim-offset;
    
    % save parameters
    stats(iJ).param.stimROI = sROI;
    stats(iJ).param.controlROI = cROI;
    stats(iJ).param.jointName = Joints{iJ};
    stats(iJ).param.control_trialsN = size(temp.control,2);
    stats(iJ).param.stim_trialsN = size(temp.stim,2);
    
    % calculate the CONTROL change in joint angle
    numflies = stats(iJ).param.control_trialsN;
    a = temp.control;
    [diff, err] = statsChange(a,cROI,sROI);
    % 'save' the data into the stats struct
    stats(iJ).OG.C_diff = diff;
    stats(iJ).OG.C_err = err;
    stats(iJ).OG.C_raw = a;
    
    % calculate the STIM change in joint angle
    numflies = stats(iJ).param.stim_trialsN;
    a = temp.stim;
    [diff, err] = statsChange(a,cROI,sROI);
    % 'save' the data into the stats struct
    stats(iJ).OG.S_diff = diff;
    stats(iJ).OG.S_err = err;
    stats(iJ).OG.S_raw = a;
    
    % diff between stim and control
    stats(iJ).OG.diff = stats(iJ).OG.S_diff-stats(iJ).OG.C_diff;
end

% figure of change in joint angle for each joint: 
% plot the change in value for all conditions:
sz = 50;
fig = getfig; idx = 0; hold all
for iJ = 1:3
    idx = idx+1;
    % pull up the diff + err values for both SH and IN:
    C_err = stats(iJ).OG.C_err;
    S_err = stats(iJ).OG.S_err;
    diff = stats(iJ).OG.diff;
    % find the err:
    err = sqrt(C_err^2+S_err^2)/sqrt(num.fly);
    
    scatter(idx, diff, sz, stimColors(iJ,:))
    errorbar(idx, diff, err, 'color', stimColors(iJ,:))
end
xlim([0, idx+1])
set(gca,'TickDir','out');
xlabel('Joint')
ylabel('Difference in joint angle before-after stim')
figname = [file_name ' diff in joint angle before-after'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);




% bootstrap the trials and make a distribution of possibilities
tic
fig = getfig;
N = 10E4;
for iJ = 1:3
    % combine the control (no laser) and stim (laser) data
    mixed_data = [stats(iJ).OG.S_raw, stats(iJ).OG.C_raw];
    test = [];
    cNum = stats(iJ).param.control_trialsN;
    sNum = stats(iJ).param.stim_trialsN;
    
  for n = 1:N

    % draw 'new' data:
    randLoc = randperm(cNum+sNum);
    C_loc = randLoc(1:cNum);
%     randLoc = randperm(cNum+sNum); % with replacement
    S_loc = randLoc(cNum+1:end);
    
    % calc the change in speed for the control:
    a = [];
    a = mixed_data(:,C_loc);
    diff = statsChange(a,cROI,sROI);
    test.C(n) = diff;
    
    % calc the change in speed for IN:
    a = [];
    a = mixed_data(:,S_loc);
    diff = statsChange(a,cROI,sROI);
    test.S(n) = diff;
    % calc the diff between SH and IN:
    test.diff(n) = test.S(n)-test.C(n);
  end
  % 'save' the data into the test struct
    stats(iJ).distrb = test;
    rdistrb = test.diff;
    subplot(1,3,iJ); hold all
    histogram(rdistrb)
    vline(stats(iJ).OG.diff, 'r-')
    p = sum(abs(rdistrb)>=abs(stats(iJ).OG.diff))/length(rdistrb);
    disp(p);
    title({['Joint: ' Joints{iJ}]; ['p = ' num2str(p)]})
    stats(iJ).p = p;
end
toc


save_figure(fig, [fig_dir, file_name, ' 10E4 bootstrap distrb change in joint angle']);
save([fig_dir, file_name, ' 10E4 bootstrap distrb change in joint angle'], 'stats')

% multiple comparisons test:
for iJ = 1:3
    p_val(iJ) = stats(iJ).p;
end
 p_err = p_val;
 [P_err,loc] = sort(p_err);
 fprintf('\n All:')
for idx = 1:length(p_val)
    q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
    r = (P_err((idx)) > (idx/length(P_err))*.05);
    fprintf(['\nInsignificant change: ' num2str(q) ' vs ' num2str(r) ' Joint: ' Joints{loc(idx)}])
end
fprintf('\n Done\n')


% find the number of trials:
fprintf('\n N''s by trial')
fprintf(['\ncontrol trials: ' num2str(cNum)])
fprintf(['\nstimulus trials: ' num2str(sNum)])
% number of flies:
fprintf('\n N''s by fly')
fprintf(['\nFlies total: ' num2str(num.fly) '\n'])


clear a alltrials ans b cond conds cROI diff duration err fig filter
clear idx ii ifly IN_flies IN_loc INavg INcntl INnum INraw INstim loc
clear mixed_data loc n N numflies p p_err P_err p_val post post_err pre
clear pre_err pval q r randLoc raw rdistrb ROI SH_flies SH_loc SHavg
clear SHcntl SHnum SHstim sROI stp strt sz tEnd tRange x S_loc sConds sNum
clear stats test temp rep offset iJ cNum cConds C_loc

%% Bootstrapping statistics for binned FeTi joint angles
% ---------------------------------------------------------------------

stats = [];
cConds = [1,8,15,22];
sConds = [2:7, 9:14];
sROI = 217:276; % 200ms period at the end of the stimulation period
cROI = 1:60; % 200ms pre stimstart

% gather all of the raw data -- filtered by only stationary flies and those
% that 'flexed' during the experiment
rep = 1;
% Edges = [35 80 100 120 180]; nbins = 4;
Edges = [20 60 90 120 180]; nbins = 4;
temp.all = [];
for ifly = 1:num.fly
   a = data(ifly).FeTi.all; % load the raw data
   filter = group(ifly).stationary(:,rep); % only include stationary flies in the analysis 
   a(:,~filter) = NaN; %removed moving trials
   filter = strcmp(ResponseData(ifly).behavior(:,1),'flex');
   a(:,~filter) = NaN; %removed nonflexing trials
   a(:,cConds) = NaN; % remove the control trials
   
   temp.all = [temp.all, a];
end
% remove the NaN instances
filter = isnan(temp.all(1,:));
temp.raw = temp.all(:,~filter);
% find the control / stim regions and remove trials with missing points
a = [temp.raw(cROI,:); temp.raw(sROI,:)];
filter = sum(isnan(a),1)>0;
temp.raw(:,filter) = [];
% sort the data into the joint angle bins:
bIdx = discretize(temp.raw(1,:), Edges); %index for bins by start angle
% [aa, idx] = discretize(temp.raw, nbins);

% adjust all to change in angle:
offset = temp.raw(1,:);
temp.raw = temp.raw-offset;



ANVdata = mean(temp.raw(sROI,:));
bIdx;
[~,~,anvSTATS] = anova1(ANVdata,bIdx);
[c,~,~,~] = multcompare(anvSTATS);

stats.data = temp;
stats.anovadata = ANVdata;
stats.anovabins = bIdx;
stats.anova = c;
stats.anvSTATS = anvSTATS;



% multiple comparisons test:
for ii = 1:size(c,1)
    namecomp{ii} = [num2str(c(ii,1)) ' vs ' num2str(c(ii,2))];
    p_val(ii) = c(ii,6);
end
 p_err = p_val;
 [P_err,loc] = sort(p_err);
 fprintf('\n All:')
for idx = 1:length(p_val)
    q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
    r = (P_err((idx)) > (idx/length(P_err))*.05);
    fprintf(['\nInsignificant change: ' num2str(q) ' was ' num2str(r) ' bins: ' namecomp{loc(idx)}])
end
fprintf('\n Done\n')


% Create a table with the data and variable names
T = table(namecomp', c(:,6), 'VariableNames', { 'comparison', 'pval'} ); disp(T)
% Write data to text file
writetable(T, [fig_dir, file_name, ' ANOVA of FeTi joint angle chnages.txt'])
save([fig_dir, file_name, ' ANOVA of FeTi joint angle chnages'], 'stats')


Edges
% find the number of trials:
fprintf('\n N''s by trial')
for bb = 1:nbins
   a = sum(bIdx == bb); 
   fprintf(['\n Bin ' num2str(bb) ' trials: ' num2str(a)])
end
% number of flies:
fprintf('\n N''s by fly')
fprintf(['\nFlies total: ' num2str(num.fly) '\n'])


clear a alltrials ans b cond conds cROI diff duration err fig filter
clear idx ii ifly IN_flies IN_loc INavg INcntl INnum INraw INstim loc
clear mixed_data loc n N numflies p p_err P_err p_val post post_err pre
clear pre_err pval q r randLoc raw rdistrb ROI SH_flies SH_loc SHavg
clear SHcntl SHnum SHstim sROI stp strt sz tEnd tRange x S_loc sConds sNum
clear stats test temp rep offset iJ cNum cConds C_loc ANVdata anvSTATS
clear avg bb bIdx c Edges nbins

%% intact fly stuff below: 

% %% Intact fly time course figures
% % load the data set with first secttion of MIDSTIM_Step_7
% % file_root = '/Volumes/Evyn SSD/Evyn UW work/matlabroot';
% file_root = 'C:\matlabroot';
% fig_dir = [file_root '/Interneuron Lines/Figures/'];
% 
% 
% type = 'walking'; % stationary
% for cond = 1:7 
%   for ifly = 1:num.fly  
%   
%     % select the data and generate a filter for the initial behavior state
%     filter = [group(ifly).(type)(cond,:), group(ifly).(type)(cond+7,:)];  
%     cw = [fly(ifly).Control.speed(cond).data(1:end-1,:); fly(ifly).Stim.speed(cond).data];
%     ccw = [fly(ifly).Control.speed(cond+7).data(1:end-1,:); fly(ifly).Stim.speed(cond+7).data];
%     input = [cw, ccw];
%     input(:,~filter) = nan; %nix the trials that didn't fit the movement type
%     % if there aren't 2 trials min, all go NaN
%     if sum(filter)<2
%         input(:,filter) = nan;
%     end
%     %stats on the fly's data
%     data(cond).(type).raw(ifly).data = input;
%     data(cond).(type).raw(ifly).avg = nanmean(input,2);
%     data(cond).(type).raw(ifly).err = sem(input,2);
%     % add averages to the group data
%     data(cond).(type).avg(:,ifly) = data(cond).(type).raw(ifly).avg;
%     data(cond).(type).err(:,ifly) = data(cond).(type).raw(ifly).err;
%   end
%   % group avg
%   data(cond).(type).AVG = nanmean(data(cond).(type).avg,2);
%   data(cond).(type).ERR = sem(data(cond).(type).avg,2);
% 
% end
% 
% % COLOR SELECTION %
% structures = {'10B-04751-gal4xUAS-csChrimson', '10B-04751-gal4xUAS-gtACR1',...
%               'SH-gal4xUAS-csChrimson', 'SH-gal4xUAS-gtACR1'};
% for ii = 1:length(structures)
%     if strcmpi(structure_name, structures{ii})
%         temp = ii;
%     end
% end
% switch temp
%     case 1
%         % Oranges (10Bxchrimson)
%         % InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'};
%         InColor = {'000000', '000000', '9c3009', '000000', 'f15a24', '000000', 'f9bca7'}; % spread out the colors
%     case 2
%         % reds for 10B gtacr1
%         InColor = {'000000', '000000', '6e1d23', '000000', '8e1f23', '000000', 'ce2323'}; % spread out the colors
%     case 3
%         % redwater marqeadon (Split half chrimson color)
%         InColor = {'000000', '000000', '634b4b', '000000', '866868', '000000', 'a98585'}; % spread out the colors
%     case 4
%         % dull dreams SHxgtacr1
%         InColor = {'000000', '000000', '77738d', '000000', 'b296a0', '000000', 'f6c2c2'}; % spread out the colors
% end
% % % Teals:
% % InColor = {'000000', '005364', '006f85', '008ba7', '46abbf', '8ac9d6', 'bde9e9'};
% 
% for ii = 1:length(InColor)
%     colorList(ii,:) =  hex2rgb(InColor{ii});     
% end
% 
% 
% % PLOT THE DATA %
% y = [];
% condLength = [0, 0.03, 0.06, 0.09, 0.180, 0.320, 0.720];
% x =  -2:1/30:2;
% sSpan = 3;
% fig = getfig;
% hold all
% %CONTROL
% y.avg = smooth(data(1).(type).AVG,sSpan);
% y.err = smooth(data(1).(type).ERR,sSpan);
% % fill_data = error_fill(x,  y.avg, y.err);
% % h = fill(fill_data.X, fill_data.Y, Color('Black'), 'EdgeColor','none');
% % set(h, 'facealpha', 0.2)
% plot(x, y.avg, 'color', Color('Black'), 'linewidth', 1)
% plot(x, y.avg+y.err, 'color', Color('Black'), 'linewidth', 0.5)
% plot(x, y.avg-y.err, 'color', Color('Black'), 'linewidth', 0.5)
% %STIM % 60, 180, 720 = 3,5,7;
% for cond = [3,5,7]
%     y.avg = smooth(data(cond).(type).AVG,sSpan);
%     y.err = smooth(data(cond).(type).ERR,sSpan);
% %     fill_data = error_fill(x, y.avg, y.err);
% %     h = fill(fill_data.X, fill_data.Y, colorList(cond,:), 'EdgeColor','none');
% %     set(h, 'facealpha', 0.2)
%     plot(x, y.avg, 'color', colorList(cond,:), 'linewidth', 1)
%     plot(x, y.avg+y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
%     plot(x, y.avg-y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
% end
% %axes & labels
% ymax = 1.6;
% offset = .03;
% ylim([0 ymax])
% LX = [0,condLength(3); 0, condLength(5); 0, condLength(7)];
% for ii = 1:3
%     LY = [ymax-(offset*ii), ymax-(offset*ii)];
%     plot(LX(ii,:), LY, 'Color', 'g', 'linewidth', 3)
% end
% % vline([0, 0, 0.6, 0.18, 0.72], 'g')
% xlabel('time (sec)')
% ylabel('Speed (cm/s)')
% title([filename(1:end-4) ' headless ' type])
% set(gca,'TickDir','out');
% % save the data
% figname = [fig_dir, filename(1:end-4) ' ' type ' speed timecourse'];
% % savefig(fig, figname);
% save_figure(fig, figname); 
% 
% 
% 
% 

% %% Quantifying the startle effects of the intact walking / stationary flies 
% % type = 'stationary';
% 
% type = 'walking';
% switch type
%     case 'walking'
%         min_speed = 0.3;
%     case 'stationary'
%         min_speed = 0.0;
% end
% 
% 
% condFrames = round(condLength*num.fps);
% windowEnd = condFrames + 60 + round(0.2*num.fps);
% windowEnd(1) = windowEnd(end);
% condFrames = round(condLength*num.fps);
% sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% % smooth of ...6 seconds
% for ifly = 1:num.fly
% %     fig1 = getfig;
%     for cond = 1:7
% %     subplot(2,4,cond)
%         % plot each individual track per lfy
% %         hold all
%         for rep = 1:6
%             vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
%                                   'gaps',[],'singlespacing',[],'duration',[]);
%             raw = data(cond).(type).raw(ifly).data(:,rep); 
%             % if the data fits the behavior type, analyze
%             if ~isnan(raw(1))
%                 temp = smooth(raw, 'rlowess');
%                 vStart = nanmean(temp(54:60)); %value prestim
%                 % %only plot trials with running before the stim ONLY FOR WALKING
%                 if vStart >= min_speed
%                     %find locations that the fly isn't moving post stim start
%                     vLow = find(temp(61:121) <= min_speed); 
%                     if ~isempty(vLow)
%                         vFreeze(rep).pos = vLow+60;
% %                         scatter(x(vFreeze(rep).pos), temp(vFreeze(rep).pos), 20, colorList{rep}, 'filled') %plot as points
% %                         plot(x, temp, 'color',  colorList{rep}, 'linewidth', 1)
% 
%                         % calculate the specific bits that qualify for 'freezing'
%                         vFreeze(rep).spacing = diff(vFreeze(rep).pos);
%                         vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
%                         %account for no gap situations
%                         try % gaps present
%                             vFreeze(rep).end = vFreeze(rep).gaps(1);
%                             vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1);           
%                             vFreeze(rep).start = vFreeze(rep).singlespacing(1);
%                             if vFreeze(rep).start>vFreeze(rep).end
%                                 vFreeze(rep).start = vFreeze(rep).end;
%                             end
%                         catch % no gaps in the 'frozen' data points
%                             vFreeze(rep).start = 1;
%                             vFreeze(rep).end = length(vFreeze(rep).pos);
%                         end
%                         vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
%                         % plot only the initial freeze data:
%                         vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
%                         if vFreeze(rep).idx(1)<windowEnd(cond)
%                              %last possible start point for freezing
% %                             scatter(x(vFreeze(rep).idx), temp(vFreeze(rep).idx), 40, 'r') %plot as points
%                         else 
%                             vFreeze(rep).idx = [];
%                             vFreeze(rep).start = [];
%                             vFreeze(rep).end = [];
%                         end 
%                     end
%                 end
%             end
%         end
%         FreezeData(cond).all(ifly).data = vFreeze;
% %         y = data(cond).(type).raw(ifly).avg;
% %         plot(x, y, 'linewidth', 2)
% %         vline([0, condLength(cond)], 'k')
% %         hline(0.3, 'k')
% %         title(['Cond ' num2str(cond)])
% %         set(gca,'TickDir','out');
%     end
% %     title({filename(1:end-4); ['Fly ' num2str(ifly)]; 'Freezing points <0.3'; ['Cond ' num2str(cond)]})
% %     savefig(fig1, [fig_dir, filename(1:end-4) ' ' type ' freezing points fly ' num2str(ifly)]);
% %     close(fig1)
% end
% 
% 
% %% Quantify the freezing data PER INSTANCE
% % FreezeData is the starting data structure to work from
% % colorList = [];
% % InColor = {'000000', '005364', '006f85', '008ba7', '46abbf', '8ac9d6', 'bde9e9'}; %teals
% % InColor = {'35063e','650021', '840000', 'b00149', 'cb0162', 'f7879a', 'ffb7ce'}; %maroons
% % Oranges (10B)
% % 
% % InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'};
% % for ii = 1:length(InColor)
% %     colorList(ii,:) =  hex2rgb(InColor{ii});     
% % end
% 
% % Find the data across all flies for each condition % 
% MT = NaN(6,1);
% for cond = 1:7    
%     for ifly = 1:num.fly
%         inputdata = FreezeData(cond).all(ifly).data;   
%         output = struct('freezeresponse', MT, 'delay', MT, 'duration', MT);
%         for rep = 1:6
%            idx = inputdata(rep).idx;
%            % did the fly freeze on that trial?
%            if isempty(idx)
%              output.freezeresponse(rep) = false;
%            else
%              output.freezeresponse(rep) = true;  
%              output.delay(rep) = (idx(1)-60)/num.fps;
%              output.duration(rep) = length(idx)/num.fps;
%            end
%         end
%         output.totalcount = sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]);
%         output.freezecount = sum(output.freezeresponse);
%         output.avgdelay = nanmean(output.delay);
%         output.avgduration = nanmean(output.duration);
%         
%         FreezeData(cond).data(ifly) = output;
%         FreezeData(cond).totalcount(ifly) = output.totalcount;
%         FreezeData(cond).freezecount(ifly) = output.freezecount;
%         FreezeData(cond).delay(:,ifly) = output.delay;
%         FreezeData(cond).duration(:,ifly) = output.duration;
%     end
%     FreezeData(cond).total = sum(FreezeData(cond).totalcount);
%     FreezeData(cond).freeze = sum(FreezeData(cond).freezecount);
%     FreezeData(cond).freq = FreezeData(cond).freeze/FreezeData(cond).total;
%     % Delay stats:
%     loc = ~isnan(FreezeData(cond).delay);
%     a = reshape(FreezeData(cond).delay(loc),FreezeData(cond).freeze,1);
%     FreezeData(cond).Delay.all = a;
%     FreezeData(cond).Delay.avg = mean(a);
%     FreezeData(cond).Delay.err = sem(a);
%     % Duration stats:
%     loc = ~isnan(FreezeData(cond).duration);
%     a = reshape(FreezeData(cond).duration(loc),FreezeData(cond).freeze,1);
%     FreezeData(cond).Duration.all = a;
%     FreezeData(cond).Duration.avg = mean(a);
%     FreezeData(cond).Duration.err = sem(a);
% end
% 
% 
% 
% % -- Plot the frequency of freezing -- %
% kolor = colorList(5,:);
% y = []; 
% fig = getfig;x = condLength;
% for cond = 1:7
%     y(cond) = FreezeData(cond).freq*100;
% end
% scatter(x, y, 70, kolor, 'filled',...
%     'markeredgecolor', 'none', 'markerfacecolor', 'flat')
% xlim([-.1, .8]); ylim([0, 100]);
% 
% xlabel('Stimulation length (sec)')
% ylabel('Frequency of freezing (%)')
% title([filename(1:end-4) ' ' type ' freeze frequency'])
% set(gca,'TickDir','out');
% % save the data
% figname = [fig_dir, filename(1:end-4) ' ' type ' freezing tuning curve'];
% % savefig(fig, figname);
% save_figure(fig, figname);
% 
% % -- Plot the freezing duration tuning curve -- %
% y = []; 
% fig = getfig;
% hold all
% for cond = 1:7
%     y(cond) = FreezeData(cond).Duration.avg;
%     err(cond) = FreezeData(cond).Duration.err;
% end
% scatter(x, y, 70, kolor, 'filled',...
%     'markeredgecolor', 'none', 'markerfacecolor', 'flat')
% errorbar(x,y,err, 'color', kolor, 'linestyle', 'none')
% xlim([-.1, .8]); ylim([0, 0.8])
% 
% xlabel('Stimulation length (sec)')
% ylabel('Freezing duration (sec)')
% title([filename(1:end-4) ' ' type ' freeze duration'])
% set(gca,'TickDir','out');
% % save the data
% figname = [fig_dir, filename(1:end-4) ' ' type ' freeze duration tuning'];
% savefig(fig, figname);
% save_figure(fig, figname);
% 
% 
% % -- Plot the freezing delay tuning curve -- %
% y = []; 
% fig = figure; set(fig, 'color', 'w');
% hold all
% for cond = 1:7
%     y(cond) = FreezeData(cond).Delay.avg;
%     err(cond) = FreezeData(cond).Delay.err;
% end
% scatter(x, y, 70, kolor, 'filled',...
%     'markeredgecolor', 'none', 'markerfacecolor', 'flat')
% errorbar(x,y,err, 'color', kolor, 'linestyle', 'none')
% xlim([-.1, .8]); % ylim([0, 0.5])
% xlabel('Stimulation length (sec)')
% ylabel('Freezing delay (sec)')
% title([filename(1:end-4) ' ' type ' freeze delay'])
% set(gca,'TickDir','out');
% % save the data
% figname = [fig_dir, filename(1:end-4) ' ' type ' freeze delay tuning'];
% savefig(fig, figname);
% save_figure(fig, figname);
% 
% 
% 
% 
% 
% %% Quantify the freezing data (data.Freez) PER FLY
% 
% % need to run the above code for each of the flies that are in the
% % structure--maybe assign freeze to a new data structure for ease with
% % handling.
% 
% % Start analysis with the 'FreezeData' structure
% % for each fly: find the duration of freezing for each trial and the avg
% % for each fly -- of those that paused and then those that didn't
% % Teals:
% colorList = [];
% InColor = {'000000', '005364', '006f85', '008ba7', '46abbf', '8ac9d6', 'bde9e9'}; %teals
% InColor = {'35063e','650021', '840000', 'b00149', 'cb0162', 'f7879a', 'ffb7ce'}; %maroons
% 
% for ii = 1:length(InColor)
%     colorList(ii,:) =  hex2rgb(InColor{ii});     
% end
% 
% MT = NaN(6,1);
% for cond = 1:7    
%     for ifly = 1:num.fly
%         inputdata = FreezeData(cond).all(ifly).data;   
%         output = struct('freezeresponse', MT, 'delay', MT, 'duration', MT);
%         for rep = 1:6
%            idx = inputdata(rep).idx;
%            % did the fly freeze on that trial?
%            if isempty(idx)
%              output.freezeresponse(rep) = false;
%            else
%              output.freezeresponse(rep) = true;  
%              output.delay(rep) = (idx(1)-60)/num.fps;
%              output.duration(rep) = length(idx)/num.fps;
%            end
%         end
%         output.totalcount = sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]);
%         output.freezecount = sum(output.freezeresponse);
%         output.avgdelay = nanmean(output.delay);
%         output.avgduration = nanmean(output.duration);
%         FreezeData(cond).data(ifly) = output;
%         FreezeData(cond).totalcount(ifly) = output.totalcount;
%         FreezeData(cond).freezecount(ifly) = output.freezecount;
%         FreezeData(cond).avgdelay(ifly) = output.avgdelay;
%         FreezeData(cond).avgduration(ifly) = output.avgduration;
%     end
% end
% % --- avg duration or avg delay ------ %
% sz = 70;
% clear MT output input
% field = 'avgduration' ;% field = 'avgdelay';  avgdelay  
% % plot the data? 
% fig = getfig; hold all
% xlength = 0.2; xlim([0,8])
% for cond = 1:7
%     loc = ~(FreezeData(cond).freezecount==0);
%     idx = sum(loc);
%     x = cond:xlength/idx:cond+(xlength-xlength/idx);
%     y = FreezeData(cond).(field)(loc);
%     scatter(x, y, sz, colorList(cond,:), 'filled')
%     avg = nanmean(FreezeData(cond).(field));
%     err = sem(FreezeData(cond).(field),2);
% %     errorbar(cond+(xlength/2), avg, err, 'Color', 'r')
%     plot([cond, cond+xlength], [avg,avg], 'Color', colorList(cond,:), 'linewidth', 3)
% end
% title([structure_name ' Freezing Behavior: ' field])
% xlabel('Condition')
% ylabel([field ' time (sec)'])
% ylim([0,2])
% set(gca,'TickDir','out');
% 
% figname = [fig_dir, structure_name ' ' type ' freezing_' field];
% savefig(fig, figname);
% macsave(fig, figname);
% 
% % ---- freezing percent ---- %
% % Exploring the proportion of flies that exhibit freezing
% fig = getfig; hold all
% xlength = 0.2; xlim([0,8])
% for cond = 1:7
%     y = FreezeData(cond).freezecount./FreezeData(cond).totalcount;
%     y(isnan(y)) = [];
%     idx = length(y);
%     x = cond:xlength/idx:cond+(xlength-xlength/idx);
%     scatter(x, y, sz, colorList(cond,:), 'filled')
%     avg = nanmean(y);
%     err = sem(y,2);
% %     errorbar(cond+(xlength/2), avg, err, 'Color', 'r')
%     plot([cond, cond+xlength], [avg,avg], 'Color', colorList(cond,:), 'linewidth', 3)
% end
% title({[structure_name ' Frequency of freezing behavior'];' ';' '})
% xlabel('Condition')
% ylabel('Percent of trials freezing (%)')
% figname  = [fig_dir, structure_name ' ' type ' percent freezing'];
% 
% savefig(fig, figname)
% macsave(fig, figname)
% 
% %% STATIONARY FLY ANALYSIS
% 
% 
% 
% type = 'stationary';
% for cond = 1:7 
%   for ifly = 1:num.fly  
%   
%     % select the data and generate a filter for the initial behavior state
%     filter = [group(ifly).(type)(cond,:), group(ifly).(type)(cond+7,:)];  
%     cw = [fly(ifly).Control.speed(cond).data(1:end-1,:); fly(ifly).Stim.speed(cond).data];
%     ccw = [fly(ifly).Control.speed(cond+7).data(1:end-1,:); fly(ifly).Stim.speed(cond+7).data];
%     input = [cw, ccw];
%     input(:,~filter) = nan; %nix the trials that didn't fit the movement type
%     %stats on the fly's data
%     data(cond).(type).raw(ifly).data = input;
%     data(cond).(type).raw(ifly).avg = nanmean(input,2);
%     data(cond).(type).raw(ifly).err = sem(input,2);
%     % add averages to the group data
%     data(cond).(type).avg(:,ifly) = data(cond).(type).raw(ifly).avg;
%     data(cond).(type).err(:,ifly) = data(cond).(type).raw(ifly).err;
%   end
%   % group avg
%   data(cond).(type).AVG = nanmean(data(cond).(type).avg,2);
%   data(cond).(type).ERR = sem(data(cond).(type).avg,2);
% 
% end
% 
% % COLOR SELECTION %
% structures = {'10B-04751-gal4xUAS-csChrimson', '10B-04751-gal4xUAS-gtACR1',...
%               'SH-gal4xUAS-csChrimson', 'SH-gal4xUAS-gtACR1'};
% for ii = 1:length(structures)
%     if strcmpi(structure_name, structures{ii})
%         temp = ii;
%     end
% end
% switch temp
%     case 1
%         % Oranges (10Bxchrimson)
%         % InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'};
%         InColor = {'000000', '000000', '9c3009', '000000', 'f15a24', '000000', 'f9bca7'}; % spread out the colors
%     case 2
%         % reds for 10B gtacr1
%         InColor = {'000000', '000000', '6e1d23', '000000', '8e1f23', '000000', 'ce2323'}; % spread out the colors
%     case 3
%         % redwater marqeadon (Split half chrimson color)
%         InColor = {'000000', '000000', '634b4b', '000000', '866868', '000000', 'a98585'}; % spread out the colors
%     case 4
%         % dull dreams SHxgtacr1
%         InColor = {'000000', '000000', '77738d', '000000', 'b296a0', '000000', 'f6c2c2'}; % spread out the colors
% end
% % % Teals:
% % InColor = {'000000', '005364', '006f85', '008ba7', '46abbf', '8ac9d6', 'bde9e9'};
% 
% for ii = 1:length(InColor)
%     colorList(ii,:) =  hex2rgb(InColor{ii});     
% end
% 
% 
% % PLOT THE DATA %
% condLength = [0, 0.03, 0.06, 0.09, 0.180, 0.320, 0.720];
% x =  -2:1/30:2; y = [];
% sSpan = 3;
% fig = figure; set(fig, 'Color', 'w');
% hold all
% %CONTROL
% y.avg = smooth(data(1).(type).AVG,sSpan);
% y.err = smooth( data(1).(type).ERR,sSpan);
% fill_data = error_fill(x,  y.avg, y.err);
% h = fill(fill_data.X, fill_data.Y, Color('Black'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% plot(x, y.avg, 'color', Color('Black'), 'linewidth', 1)
% %STIM % 60, 180, 720 = 3,5,7;
% for cond = [3,5,7]
%     y.avg = smooth(data(cond).(type).AVG,sSpan);
%     y.err = smooth( data(cond).(type).ERR,sSpan);
%     fill_data = error_fill(x, y.avg, y.err);
%     h = fill(fill_data.X, fill_data.Y, colorList(cond,:), 'EdgeColor','none');
%     set(h, 'facealpha', 0.2)
%     plot(x, y.avg, 'color', colorList(cond,:), 'linewidth', 1)
% end
% %axes & labels
% ymax = 1.6;
% offset = .03;
% ylim([0 ymax])
% LX = [0,condLength(3); 0, condLength(5); 0, condLength(7)];
% for ii = 1:3
%     LY = [ymax-(offset*ii), ymax-(offset*ii)];
%     plot(LX(ii,:), LY, 'Color', 'g', 'linewidth', 3)
% end
% % hline(0.3, 'k-')
% % vline(-0.2, 'r')
% % vline([0, 0, 0.6, 0.18, 0.72], 'g')
% xlabel('time (sec)')
% ylabel('Speed (cm/s)')
% title([filename(1:end-4) ' headless ' type])
% set(gca,'TickDir','out');
% % save the data
% figname = [fig_dir, filename(1:end-4) ' ' type ' speed timecourse'];
% savefig(fig, figname);
% macsave(fig, figname);
% 
% 
% %% Stationary starts on a cond-by-cond level 
% 
% for start_cond = [3,5,7]
% stim_conds = start_cond:7:28;
% newdata = 1;
% % pointrange: 500ms period after light offset
% range_list = [0,0; 1,16; 2,17; 3,18; 6,21; 12,27; 22,37];
% pointrange = range_list(start_cond,1):range_list(start_cond,2); %long light
% % Stationary flies that start walking after activation|silencing
% % save an individual image and then the data, then, once all four crosses
% % are run, make the figure with all the appropriate data
% type = {'control', 'light'};
% min_speed = 0.3;
% 
% num.divide = 12;
% ind = 0;
% Pdata = struct;
% if newdata == 1
%   for kk = 1:num.fly
%     %control information
%     Pdata.control.speed(kk,1:num.divide) = nan; %fill the trial spots with NaN since moving flies won't count 
% 
%     idx = 0;
%     for cond = [1,8,15,22] %control stimuli
%         for rep = 1:num.reps
%             idx = idx + 1;
%             if group(kk).stationary(cond,rep) == true
%                 Pdata.control.speed(kk,idx) = nanmean(fly(kk).Stim.speed(cond).data(pointrange, rep));
%             end
%         end
%     end
%     Pdata.control.stationary(kk) = (sum(~isnan(Pdata.control.speed(kk,:))));
%     Pdata.control.idx(kk,:) = Pdata.control.speed(kk,:)>=min_speed;
%     
%     Pdata.control.fraction(kk) = sum(Pdata.control.idx(kk,:))/Pdata.control.stationary(kk)*100;
%     Pdata.control.increase(kk) = nanmean(Pdata.control.speed(kk,Pdata.control.idx(kk,:)));
% 
%     
%     % stimulus information
%     Pdata.light.speed(kk,1:num.divide) = nan;
%     idx = 0;
%     for cond = stim_conds %
%         for rep = 1:num.reps
%             idx = idx + 1;
%             if group(kk).stationary(cond,rep) == true
%                 Pdata.light.speed(kk,idx) = mean(fly(kk).Stim.speed(cond).data(pointrange, rep));
%             end
%         end
%     end
%     Pdata.light.stationary(kk) = (sum(~isnan(Pdata.light.speed(kk,:))));
%     Pdata.light.idx(kk,:) = Pdata.light.speed(kk,:)>=min_speed;
%     Pdata.light.fraction(kk) = sum(Pdata.light.idx(kk,:))/Pdata.light.stationary(kk)*100;
%     Pdata.light.increase(kk) = nanmean(Pdata.light.speed(kk,Pdata.light.idx(kk,:)));
% 
%     Pdata.plot(kk,1) = Pdata.control.fraction(kk);
%     Pdata.plot(kk,2) = Pdata.light.fraction(kk);
% 
%   end
%     
%   for tt = 1:2
%     Pdata.(type{tt}).totalstationary = sum(Pdata.(type{tt}).stationary);
%     Pdata.(type{tt}).totalmoved = sum(sum(Pdata.(type{tt}).idx));
%     Pdata.(type{tt}).percent = Pdata.(type{tt}).totalmoved/Pdata.(type{tt}).totalstationary*100;
%     Pdata.bargraph(tt) = Pdata.(type{tt}).percent;
%   end
% 
% end
% 
% 
% fig = figure; set(fig, 'color', 'w')
% bar([1,2], Pdata.bargraph, 'BarWidth', 1)
% ylim([0,100])
% xlabel(['total stationary: ' num2str(Pdata.control.totalstationary),...
%         ' vs ' num2str(Pdata.light.totalstationary)])
%     
% title({structure_name; 'fly starts control vs stim'; ['cond ' num2str(start_cond)]})
% set(gca,'TickDir','out');
% figure_name = [figures_dir, structure_name, ' Fly starts from stationary cond '  num2str(start_cond)];
% % save(figure_name, 'Pdata')
% savefig(fig, figure_name)
% macsave(fig, figure_name);
% end 
% 
% 
% %% stationary start probability for the flies
% %All of the conditions run together on the same graph
%  
% temp = [];
% for start_cond = 2:7
%     stim_conds = start_cond:7:28;
%     newdata = 1;
%     % pointrange: 500ms period after light offset
%     range_list = [0,0; 1,16; 2,17; 3,18; 6,21; 12,27; 22,37];
%     pointrange = range_list(start_cond,1):range_list(start_cond,2); %long light
%     % Stationary flies that start walking after activation|silencing
%     % save an individual image and then the data, then, once all four crosses
%     % are run, make the figure with all the appropriate data
%     type = {'control', 'light'};
%     min_speed = 0.3;
% 
%     num.divide = 12;
%     ind = 0;
%     Pdata = struct;
%     if newdata == 1
%       for kk = 1:num.fly
%         %control information
%         Pdata.control.speed(kk,1:num.divide) = nan; %fill the trial spots with NaN since moving flies won't count 
% 
%         idx = 0;
%         for cond = [1,8,15,22] %control stimuli
%             for rep = 1:num.reps
%                 idx = idx + 1;
%                 if group(kk).stationary(cond,rep) == true
%                     Pdata.control.speed(kk,idx) = nanmean(fly(kk).Stim.speed(cond).data(pointrange, rep));
%                 end
%             end
%         end
%         Pdata.control.stationary(kk) = (sum(~isnan(Pdata.control.speed(kk,:))));
%         Pdata.control.idx(kk,:) = Pdata.control.speed(kk,:)>=min_speed;
% 
%         Pdata.control.fraction(kk) = sum(Pdata.control.idx(kk,:))/Pdata.control.stationary(kk)*100;
%         Pdata.control.increase(kk) = nanmean(Pdata.control.speed(kk,Pdata.control.idx(kk,:)));
% 
% 
%         % stimulus information
%         Pdata.light.speed(kk,1:num.divide) = nan;
%         idx = 0;
%         for cond = stim_conds %
%             for rep = 1:num.reps
%                 idx = idx + 1;
%                 if group(kk).stationary(cond,rep) == true
%                     Pdata.light.speed(kk,idx) = mean(fly(kk).Stim.speed(cond).data(pointrange, rep));
%                 end
%             end
%         end
%         Pdata.light.stationary(kk) = (sum(~isnan(Pdata.light.speed(kk,:))));
%         Pdata.light.idx(kk,:) = Pdata.light.speed(kk,:)>=min_speed;
%         Pdata.light.fraction(kk) = sum(Pdata.light.idx(kk,:))/Pdata.light.stationary(kk)*100;
%         Pdata.light.increase(kk) = nanmean(Pdata.light.speed(kk,Pdata.light.idx(kk,:)));
% 
%         Pdata.plot(kk,1) = Pdata.control.fraction(kk);
%         Pdata.plot(kk,2) = Pdata.light.fraction(kk);
% 
%       end
% 
%       for tt = 1:2
%         Pdata.(type{tt}).totalstationary = sum(Pdata.(type{tt}).stationary);
%         Pdata.(type{tt}).totalmoved = sum(sum(Pdata.(type{tt}).idx));
%         Pdata.(type{tt}).percent = Pdata.(type{tt}).totalmoved/Pdata.(type{tt}).totalstationary*100;
%         Pdata.bargraph(tt) = Pdata.(type{tt}).percent;
%       end
% 
%     end
%     % pull out the data to use for plotting multiple light lengths
%     temp.data(start_cond,:) = Pdata.bargraph;
%     temp.stimN(start_cond) = Pdata.light.totalstationary;
%     temp.cntlN(start_cond) = Pdata.control.totalstationary;
% end
% 
% % plot the data
% x = [1 2; 4 5; 7 8; 10 11; 13 14; 16 17; 19 20];
% fig = figure; set(fig, 'color', 'w'); hold all
% for ii = 1:6
%     bar(x(ii,:), temp.data(ii+1,:), 'BarWidth', 1)    
% end
% ylim([0,100])
% xlabel('Condition')
% ylabel('Movement initiation (%)')
%         
% title({structure_name; 'fly starts control vs stim'})
% set(gca,'TickDir','out');
% figure_name = [figures_dir, structure_name, ' ALL fly starts from stationary'];
% % save(figure_name, 'Pdata')
% savefig(fig, figure_name)
% macsave(fig, figure_name);
% 
% 
% 
% %% color work
% %Mat Colored Rainbow
% colorList = [];
% InColor = {'000000', 'd8a810', '4f66a5', 'b85b5b', '347455', 'a1468f', '950c15'};
% %Bright Rainbow
% InColor = {'0504AA';  'A2CFFE'; '90E4C1'; '0B8B87'; '5CAC2D';...
%              'EEDC5B'; 'FF9408'; 'FF724C';  'FB2943'; '7B002C'};   
% %Maroon colors             
% InColor = {'35063e','650021', '840000', 'b00149', 'cb0162', 'f7879a', 'ffb7ce'};
% % Teals:
% InColor = {'000000', '005364', '006f85', '008ba7', '46abbf', '8ac9d6', 'bde9e9'};
% 
% for ii = 1:length(InColor)
%     colorList(ii,:) =  hex2rgb(InColor{ii});     
% end
% 
% figure; hold all
% for ii = 1:length(colorList)
%     plot(1:4,1+ii:4+ii, 'color', colorList(ii,:), 'linewidth', 20)  
% end

