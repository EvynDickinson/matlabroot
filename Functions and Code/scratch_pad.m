

% find the diff between the change in running speed in SH and IN
stats = [];
sROI = 66:72; % 200ms period 200ms post stim start
cROI = 55:61; % 200ms pre stimstart
for cond = 1:7
    stats(cond).param.stimROI = sROI;
    stats(cond).param.controlROI = cROI;
   % SPLIT HALF difference in speed:
    numflies = sum(~isnan(Cdata(cond).walking.avg(1,:)));
    a = Cdata(cond).stats.raw; % all trials for this condition
    ROI = a(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    a(:,filter) = []; %remove non walkers
    [diff, err]= statsChange(a,cROI,sROI);
    % 'save' the data into the stats struct
    stats(cond).OG.SH_diff = diff;
    stats(cond).OG.SH_err = err;
    
    % INTERNEURON difference in speed:
    numflies = sum(~isnan(Sdata(cond).walking.avg(1,:)));
    a = Sdata(cond).stats.raw; % all trials for this condition
    ROI = a(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    a(:,filter) = []; %remove non walkers
    % calc the change in speed for that condition:
    [diff, err]= statsChange(a,cROI,sROI);
    % 'save' the data into the stats struct
    stats(cond).OG.IN_diff = diff;
    stats(cond).OG.IN_err = err;
    % find the difference between the two changes in speed:
    stats(cond).OG.diff = stats(cond).OG.IN_diff-stats(cond).OG.SH_diff;
end

% plot the change in value for all conditions:
sz = 50;
fig = getfig; idx = 0; hold all
for cond = [1, 4, 7]
    idx = idx+1;
    % pull up the diff + err values for both SH and IN:
    SH_err = stats(cond).OG.SH_err;
    IN_err = stats(cond).OG.IN_err;
    diff = stats(cond).OG.diff;
    % find the err:
    err = sqrt(SH_err^2+IN_err^2)/sqrt(min([CONTROL.num.fly, STIM.num.fly]));
    
    scatter(idx, diff, sz, kolor(2).K(4,:))
    errorbar(idx, diff, err, 'color', kolor(2).K(4,:))
end
xlim([0, idx+1])



% bootstrap the trials and make a distribution of possibilities
tic
fig = getfig;
N = 10E4;
for cond = 1:7
    % SPLIT HALF difference in speed:
    SH_flies = sum(~isnan(Cdata(cond).walking.avg(1,:)));
    a = Cdata(cond).stats.raw; % all trials for this condition
    ROI = a(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    a(:,filter) = []; %remove non walkers
    
    IN_flies = sum(~isnan(Sdata(cond).walking.avg(1,:)));
    b = Sdata(cond).stats.raw; % all trials for this condition
    ROI = b(30:60,:); %screen for non walkers
    filter = median(ROI)<=0.3;
    b(:,filter) = []; %remove non walkers
    
    mixed_data = [a,b];
     test = [];
  for n = 1:N

    % draw 'new' data:
    randLoc = randperm(SH_flies+IN_flies);
    SH_loc = randLoc(1:SH_flies);
    IN_loc = randLoc(SH_flies+1:end);
    
    % calc the change in speed for SH:
    a = [];
    a = mixed_data(:,SH_loc);
    diff = statsChange(a,cROI,sROI);
    test.SH(n) = diff;
    
    % calc the change in speed for IN:
    a = [];
    a = mixed_data(:,IN_loc);
    diff = statsChange(a,cROI,sROI);
    test.IN(n) = diff;
    % calc the diff between SH and IN:
    test.diff(n) = test.IN(n)-test.SH(n);
  end
  % 'save' the data into the test struct
    stats(cond).distrb = test;
    rdistrb = test.diff;
    subplot(2,4,cond); hold all
    histogram(rdistrb)
    vline(stats(cond).OG.diff, 'r-')
    p = sum(abs(rdistrb)>=abs(stats(cond).OG.diff))/length(rdistrb);
    disp(p);
    title({['cond: ' num2str(cond)]; ['p = ' num2str(p)]})
    stats(cond).p = p;
end
toc


%% Delay to freezing in walking flies (by fly)  -- COMBINE FOR MULTIPLE DATA SETS

% Delay data
stppoint = 1.5;
windowEnd = (60+stppoint*num.fps)*ones(1,7);
% stppoint = 0.5;
% windowEnd = condFrames + 60 + round(stppoint*num.fps);
% windowEnd(1) = windowEnd(end);


sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for cond = 1:7
    for ifly = 1:num.fly
        for rep = 1:6
            vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                  'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
            raw = data(cond).(type).raw(ifly).data(:,rep); 
            % if the data fits the behavior type, analyze
            if ~isnan(raw(1))
                temp = smooth(raw, 'rlowess');
                vStart = nanmean(raw(54:60)); %value prestim
%                 vChange = NaN;
                
                % only include trials with greater than min speed
                if vStart >= min_speed
                    % Change in speed
                    vStim = nanmean(raw(76:82)); %value prestim
                    vChange = vStim-vStart;
                    vFreeze(rep).speedchange = vChange;
                    %find locations that the fly isn't moving post stim start
                    vLow = find(temp(61:121) <= min_speed); 
                    if ~isempty(vLow)
                        vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                        % calculate the specific bits that qualify for 'freezing'
                        vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                        vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                        % find the start and end points of the 'freeze'
                        try % gaps present
                            vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                            vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                            vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                            if vFreeze(rep).start>vFreeze(rep).end
                                vFreeze(rep).start = vFreeze(rep).end;
                            end
                        catch % no gaps in the 'frozen' data points
                            vFreeze(rep).start = 1;
                            vFreeze(rep).end = length(vFreeze(rep).pos);
                        end
                        vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                        vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                        if vFreeze(rep).idx(1)>windowEnd(cond)
                            vFreeze(rep).idx = [];
                            vFreeze(rep).start = [];
                            vFreeze(rep).end = [];
                        end 
                    end
                end
            end
        end
        FreezeData(cond).all(ifly).data = vFreeze;
    end
end


MT = NaN(6,1);
for cond = 1:7   
    %empty data here -- so we don't factor in zeros
    FreezeData(cond).avgchange = NaN(1,num.fly);
    FreezeData(cond).avgdelay = NaN(1,num.fly);
    for ifly = 1:num.fly
        inputdata = FreezeData(cond).all(ifly).data;   
        output = struct('freezeresponse', MT, 'delay', MT, 'duration', MT, 'change', MT);
        for rep = 1:6
           idx = inputdata(rep).idx;
           % did the fly freeze on that trial?
           if isempty(idx)
             output.freezeresponse(rep) = false;
           else
             output.freezeresponse(rep) = true;  
             output.delay(rep) = (idx(1)-60)/num.fps;
             output.duration(rep) = length(idx)/num.fps;
             output.change(rep) = inputdata(rep).speedchange;
           end
        end
        tot = (sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]));
        if tot >= 2 %must have a possibility of 2 trials for each fly to include
            output.freezecount = sum(output.freezeresponse);
            output.totalcount = sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]);
            output.avgdelay = nanmean(output.delay);
            output.avgduration = nanmean(output.duration);
            output.avgchange = nanmean(output.change);
        else
            output.freezecount = NaN;
            output.totalcount = NaN;
            output.avgdelay = NaN;
            output.avgduration = NaN;
            output.avgchange = NaN;
        end    
            FreezeData(cond).data(ifly) = output;
            FreezeData(cond).totalcount(ifly) = output.totalcount;
            FreezeData(cond).freezecount(ifly) = output.freezecount;
            FreezeData(cond).percent(ifly) = FreezeData(cond).freezecount(ifly)/FreezeData(cond).totalcount(ifly);
            FreezeData(cond).avgdelay(ifly) = output.avgdelay;
            FreezeData(cond).avgduration(ifly) = output.avgduration;
            FreezeData(cond).avgchange(ifly) = output.avgchange;
    end
end


% -- line plot delay to freeze -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(FreezeData(cond).avgdelay);
    y.err(cond) = sem(FreezeData(cond).avgdelay,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([0, 1]);
title({[structure_name ' Freezing Delay. End point: stimend+' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Time to freeze initiation (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, structure_name ' ' type ' freezing delay linegraph ' num2str(stppoint)];

% save the data
figname  = [fig_dir, structure_name ' ' type ' freezing lineplot endpoint ' num2str(stppoint)];
save_figure(fig, figname);





% Time delay for freezing
fig = getfig; hold all
xlength = 0.2; xlim([0,8])
for cond = 1:7
    y = FreezeData(cond).avgdelay;
    y(isnan(y)) = [];
    idx = length(y);
    x = cond:xlength/idx:cond+(xlength-xlength/idx);
    scatter(x, y, sz, Kolor, 'filled')
    avg = nanmean(y);
    err = sem(y,2);
%     errorbar(cond+(xlength/2), avg, err, 'Color', 'r')
    plot([cond-(xlength/2), cond+(xlength*1.5)], [avg,avg], 'Color', 'k', 'linewidth', 3)
%     errorbar(cond, avg, err, 'color', 'k')
end
title({[structure_name ' Freezing Delay. End point: stimend+' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Time to freeze initiation (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, structure_name ' ' type ' freezing delay endpoint ' num2str(stppoint)];
ylim([0, 1.5])

save_figure(fig, figname);







%% Heatmap of speed for all flies within a trial
Scombo = []; 
MT = NaN(3,121);
%Stim data
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:STIM.num.fly
       a = Sdata(cond).walking.raw(ifly).data;
       loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Scombo = [Scombo; MT;All_sorted]; %
end
% control data
Ccombo = [];
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:CONTROL.num.fly
       a = Cdata(cond).walking.raw(ifly).data;
       loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Ccombo = [Ccombo;MT; All_sorted]; %
end

fig = getfig; 
subplot(2,1,1)
    imagesc(Scombo)%(:,30:90)
    clims =  [0, 3];
    caxis(clims);
    hold all
    colormap hot
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
%     xlabel('time')
    ylabel('Trials by Cond')  
    title([comp ' interneuron data'])
subplot(2,1,2)
    imagesc(Ccombo)% (:,30:90)
    clims =  [0, 3];
    caxis(clims);
    hold all
    colormap hot
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
    xlabel('time')
    ylabel('Trials by Cond')  
    title([comp ' SH data'])

    
figname  = [fig_dir, comp ' WALKING trials heatmap of speed'];

save_figure(fig, figname);

% fig = getfig;
% for cond = 1:7
%     ALL = []; All_sorted = [];
%     for ifly = 1:STIM.num.fly
%        a = Sdata(cond).walking.raw(ifly).data;
%        loc = isnan(a(1,:)); a(:,loc) = [];
%        ALL = [ALL; a'];
%     end
%     % sort by initial speed
%     avg_speed = nanmean(ALL(:,50:60),2);
%     [~, idx] = sort(avg_speed);
%     All_sorted = ALL(idx,:);
%     %PLOT the data
%     subplot(2,4,cond)
%     imagesc(All_sorted(:,30:90))
%     clims =  [0, 3];
%     caxis(clims);
%     hold all
%     colormap hot
%     % Adjust colorbar
%     c = colorbar;
%     c.Label.String = 'Speed (cm/s)';
%     xlabel('time')
%     ylabel(['Trials in Cond ' num2str(cond)])
% end


%% heatmap of all flies -- not just the walking flies

Scombo = []; 
MT = NaN(3,121);
%Stim data
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:STIM.num.fly
       a = Sdata(cond).ALL(ifly).data;
%        loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Scombo = [Scombo; MT;All_sorted];% 
end
% control data
Ccombo = [];
for cond = 1:7
    ALL = []; All_sorted = [];
    for ifly = 1:CONTROL.num.fly
       a = Cdata(cond).ALL(ifly).data;
%        loc = isnan(a(1,:)); a(:,loc) = [];
       ALL = [ALL; a'];
    end
    % sort by initial speed
    avg_speed = nanmean(ALL(:,50:60),2);
    [~, idx] = sort(avg_speed);
    All_sorted = ALL(idx,:);
    %PLOT the data
    Ccombo = [Ccombo; MT;All_sorted];% 
end

fig = getfig; 
subplot(2,1,1)
    imagesc(Scombo) % (:,30:90)
    clims =  [0, 3];
    caxis(clims);
    hold all
    colormap hot
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
%     xlabel('time')
    ylabel('Trials by Cond')  
    title([comp ' interneuron data'])
subplot(2,1,2)
    imagesc(Ccombo) %(:,30:90)
    clims =  [0, 3];
    caxis(clims);
    hold all
    colormap hot
    % Adjust colorbar
    c = colorbar;
    c.Label.String = 'Speed (cm/s)';
    xlabel('time')
    ylabel('Trials by Cond')  
    title([comp ' SH data'])
    
figname  = [fig_dir, comp ' ALL trials heatmap of speed'];

save_figure(fig, figname);


%% Reset the joint angles based on the adjusted select points
type = {'stim', 'control'};
rep = 1;
behavior = 'Stationary';
num.fly = 3;
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
    


    




  for iJ = 1:3 %for each joint measured
    data(ifly).(Joints{iJ}).all = NaN(337,28); 
    data(ifly).(Joints{iJ}).raw = NaN(337,28);
    for cond = condlist
        idx = 0;
        % screen for desired start behaviors
        if strcmp(group(ifly).behavior{cond},behavior) % behavior matched
            for iFrame = 90:426 %specific time period analyzed (200ms pre : 200ms post)
  
                input = tracking(ifly).Left_Front(cond,1).frame(iFrame).(Joints{iJ});
                idx = idx+1;
                if sum(input) == 0 %carry NaNs through
                    input = NaN;
                end
                data(ifly).(Joints{iJ}).raw(idx,cond) = input;
            end
            % fill in the gaps with interpolated data:
            xstart = 1:337;
            vstart = data(ifly).(Joints{iJ}).raw(:,cond);
            nanloc = isnan(vstart);
            xq = find(nanloc==true);
            v = vstart(~nanloc);
            x = xstart(~nanloc);
            vq = interp1(x,v,xq); 
            allV = vstart; 
            allV(nanloc) = vq;
            data(ifly).(Joints{iJ}).all(:,cond) = allV;
        end
    end
    % filter out flies that are moving before the stimulus starts
    if strcmpi(file_name, '13B-20847-csChrimson-7V-offball tracking data.mat')
        if ifly == 1 %&& iJ == 2 % fly 1, cond 1|11 -> fly moving
            data(ifly).(Joints{iJ}).all(:,1) = NaN(337,1);
            data(ifly).(Joints{iJ}).all(:,11) = NaN(337,1);
        end
    end
    % Separate into Stim vs Control
    temp = logical(1:28); 
    filter = temp; filter(controlConds) = false; % eliminate the controls
    data(ifly).(Joints{iJ}).stim.all = data(ifly).(Joints{iJ}).all(:,filter);
    filter = temp; filter(:) = false; filter(controlConds) = true; %eliminate the stim conditions
    data(ifly).(Joints{iJ}).control.all = data(ifly).(Joints{iJ}).all(:,filter);    

    % % Find the average across all trials for a given fly% %
    % find the average
    data(ifly).(Joints{iJ}).stim.avg = nanmean(data(ifly).(Joints{iJ}).stim.all,2);
    data(ifly).(Joints{iJ}).control.avg = nanmean(data(ifly).(Joints{iJ}).control.all,2);
    % find the error
    data(ifly).(Joints{iJ}).stim.err = sem(data(ifly).(Joints{iJ}).stim.all,2,2);
    data(ifly).(Joints{iJ}).control.err = sem(data(ifly).(Joints{iJ}).control.all,2,2);
    
  end




%% compare delay and walking speed...

hlimit = 0.5;
itype = {'CONTROL', 'STIM'};
sz = 60;
MT = NaN(6,1);    
RNGE = 54:61;
fig = getfig; hold all; idx = 0;
num.fly = size(FreezeData(1).CONTROL,2);
for dd = 1:2 %SH vs IN
    if dd == 1
       data = Cdata;
       klr = 'k';
    else
        data = Sdata;
        klr = 'r';
    end
    for cond = 1:7 %[1, 4:6]  
        idx = idx+1;
        subplot(2,7,idx); hold all; ylim([0 2]); xlim([0 2.5])
        % load the delay and duration data into 'output'
        output = struct('delay', MT, 'duration', MT, 'speed', MT);
        for ifly = 1:num.fly
            inputdata = FreezeData(cond).(itype{dd})(ifly).data; 
            output.duration = MT; output.delay = MT; output.speed = MT;
            for rep = 1:6
                if ~isempty(inputdata(rep).duration)
                    output.duration(rep) = inputdata(rep).duration/fps;  
                    output.delay(rep) = (inputdata(rep).pos(1)-60)/fps;
                    % avg speed pre-stim onset
                    a = nanmean(data(cond).walking.raw(ifly).data(RNGE,rep));
                    output.speed(rep) = a;
                end
            end
            if sum(~isnan(output.speed))>1
                output.allduration(:,ifly) = output.duration;
                output.alldelay(:,ifly) = output.delay;
                output.allspeed(:,ifly) = output.speed;
            else
                output.allduration(:,ifly) = MT;
                output.alldelay(:,ifly) = MT;
                output.allspeed(:,ifly) = MT;
            end
        end
        % group averages for the delay and duration

        loc = ~isnan(output.allduration);
        %duration
        output.dur_list = output.allduration(loc);
        output.dur_avg = mean(output.dur_list);
        output.dur_err = nanstd(output.dur_list);
        %delay
        output.delay_list = output.alldelay(loc);
        output.delay_avg = mean(output.delay_list);
        output.delay_err = nanstd(output.delay_list);
        %speed
        output.speed_list = output.allspeed(loc);
        output.speed_avg = mean(output.speed_list);
        output.speed_err = nanstd(output.speed_list);
        % how many fall into the right delay group?
        TOTAL = length(output.delay_list);
        withinLim = sum(output.delay_list<=hlimit);
        percentLim = (withinLim/TOTAL)*100;
        title({[num2str(percentLim) '% under ' num2str(hlimit) 's'];''})

        % all points
        x = output.speed_list;
        y = output.delay_list;
        
        scatter(x,y,sz,kolor(dd).K(cond,:), 'filled')
%       
        h = lsline;
        h.Color = klr;
      
        set(gca,'TickDir','out');
        hline(condLength(cond), 'g')
        hline(hlimit, 'm')
        ylim([0 2]); xlim([0 2.5])
    end 
end
subplot(2,7,1)
xlabel('Speed at stim onset (cm/s)')
ylabel('Delay to freeze onset (s)')
figname  = [fig_dir, comp ' speed vs delay percentages'];
% title({[comp ' speed vs delay'];'';''})

save_figure(fig, figname);


%% load multiple FreezeData structures to compare across genotypes

pathroot = 'C:\matlabroot\Interneuron Lines\Freeze Data';
complist = {'10B vs SH Chrimson', '13B vs SH Chrimson', '9A vs SH Chrimson'};
for igen = 1:3
    iname = complist{igen};
    FD(igen) = load([pathroot, '/', iname, ' FreezeData']);
end; clear iname;

for igen = 1:3
    InColor = {'9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8';...% 10B oranges
               '00005C', '003AA8', '006DDb', '198bff', '50a7ff', '8ac4ff', 'b6daff';...% 13B blues
               '16005F', '490092', '6400c8', '7f00ff', '9a36ff', 'B66dfe', 'd1a3ff'}; % 9A purples
    for cond = 1:7
        kolor(2).K(cond,:) = hex2rgb(InColor{igen,cond});
    end
    % SH greyscale controls
    kolor(1).K = Color('black', 'lightgrey', 7);
    % add to structure
    FD(igen).color = kolor;
end

%% plot the freeze duration for each condition:

itype = {'CONTROL', 'STIM'};
sz = 60;
MT = NaN(6,1);    
fig = getfig; hold all

for igen = 1:3
    FreezeData = FD(igen).FreezeData;
    num.fly = size(FreezeData(cond).CONTROL,2);
    for dd = 1:2 %SH vs IN
        for cond = [1, 4:6]  
            % load the delay and duration data into 'output'
            output = struct('delay', MT, 'duration', MT);
            for ifly = 1:num.fly
                inputdata = FreezeData(cond).(itype{dd})(ifly).data; 
                output.duration = MT; output.delay = MT;
                for rep = 1:6
                    if ~isempty(inputdata(rep).duration)
                        output.duration(rep) = inputdata(rep).duration/fps;  
                        output.delay(rep) = (inputdata(rep).pos(1)-60)/fps;
                    end
                end
                output.allduration(:,ifly) = output.duration;
                output.alldelay(:,ifly) = output.delay;
            end
            % group averages for the delay and duration

            loc = ~isnan(output.allduration);
            output.dur_list = output.allduration(loc);
            output.dur_avg = mean(output.dur_list);
            output.dur_err = nanstd(output.dur_list);
            output.delay_list = output.alldelay(loc);
            output.delay_avg = mean(output.delay_list);
            output.delay_err = nanstd(output.delay_list);

            % plot the data for STIM/CONTROL:
            if cond == 1
                scatter(output.delay_avg, output.dur_avg, sz*3, 'g', 'filled')
            end
%             % all points
%             scatter(output.delay_list, output.dur_list, sz, FD(igen).color(dd).K(cond,:), 'filled')
            % avg points
            errorbar(output.delay_avg,output.dur_avg,output.dur_err,output.dur_err,output.delay_err,output.delay_err,...
                'color', FD(igen).color(dd).K(cond,:))
            scatter(output.delay_avg, output.dur_avg, sz*2, FD(igen).color(dd).K(cond,:), 'filled')
        end 
    end 
end
xlabel('Freeze delay (s)')
ylabel('Freeze duration (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, 'CsChrimson freeze delay vs duration'];
title({'CsChrimson freeze delay vs duration';'';''})
save_figure(fig, figname);


% delay/duration:


fig = getfig; hold all
for igen = 1:3
    FreezeData = FD(igen).FreezeData;
    num.fly = size(FreezeData(cond).CONTROL,2);
    for dd = 1:2 %SH vs IN
        for cond = [1, 4:6]  
            % load the delay and duration data into 'output'
            output = struct('delay', MT, 'duration', MT);
            for ifly = 1:num.fly
                inputdata = FreezeData(cond).(itype{dd})(ifly).data; 
                output.duration = MT; output.delay = MT;
                for rep = 1:6
                    if ~isempty(inputdata(rep).duration)
                        output.duration(rep) = inputdata(rep).duration/fps;  
                        output.delay(rep) = (inputdata(rep).pos(1)-60)/fps;
                    end
                end
                output.allduration(:,ifly) = output.duration;
                output.alldelay(:,ifly) = output.delay;
            end
            % group averages for the delay and duration

            loc = ~isnan(output.allduration);
            output.dur_list = output.allduration(loc);
            output.dur_avg = mean(output.dur_list);
            output.dur_err = nanstd(output.dur_list);
            output.delay_list = output.alldelay(loc);
            output.delay_avg = mean(output.delay_list);
            output.delay_err = nanstd(output.delay_list);
            output.avg_metric = output.dur_avg/output.delay_avg;
            output.metric = output.dur_list./output.delay_list;
            
            y = output.metric;
            x = condLength(cond)*ones(length(y),1);
            scatter(x, y, sz, FD(igen).color(dd).K(cond,:), 'filled')
            
            % avg
            scatter(condLength(cond), output.avg_metric, sz*2, 'g', 'filled')
%            
%             % plot the data for STIM/CONTROL:
%             if cond == 1
%                 scatter(output.delay_avg, output.dur_avg, sz*3, 'g', 'filled')
%             end
%             % all points
%             scatter(output.delay_list, output.dur_list, sz, FD(igen).color(dd).K(cond,:), 'filled')
%             avg points
%             errorbar(output.delay_avg,output.dur_avg,output.dur_err,output.dur_err,output.delay_err,output.delay_err,...
%                 'color', FD(igen).color(dd).K(cond,:))
%             scatter(output.delay_avg, output.dur_avg, sz*2, FD(igen).color(dd).K(cond,:), 'filled')
        end 
    end 
end
















% -- line plot delay to freeze -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(FreezeData(cond).avgdelay);
    y.err(cond) = sem(FreezeData(cond).avgdelay,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([0, 1]);
title({[structure_name ' Freezing Delay. End point: stimend+' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Time to freeze initiation (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, structure_name ' ' type ' freezing delay linegraph ' num2str(stppoint)];

% save the data
figname  = [fig_dir, structure_name ' ' type ' freezing lineplot endpoint ' num2str(stppoint)];
save_figure(fig, figname);




















































%% Probability density function
type = 'walking';
min_speed = 0.3;

% percent freezing
stppoint = 0.500;
windowEnd = (60+stppoint*fps)*ones(1,7);

sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for dd = 1:2
    switch dd
        case 1 %SH control
            data = Cdata;
            group = CONTROL.group;
            num = CONTROL.num;
        case 2 %INT stim
            data = Sdata;
            group = STIM.group;
            num = STIM.num;
    end
    for cond = 1:7
        for ifly = 1:num.fly
            for rep = 1:6
                vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                      'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
                raw = data(cond).(type).raw(ifly).data(:,rep); 
                % if the data fits the behavior type, analyze
                if ~isnan(raw(1))
                    temp = smooth(raw, sSpan, 'rlowess');
                    vStart = nanmean(raw(54:60)); %value prestim
    %                 vChange = NaN;
                    % only include trials with greater than min speed
                    if vStart >= min_speed
                        % Change in speed
                        vStim = nanmean(raw(76:82)); %value prestim
                        vChange = vStim-vStart;
                        vFreeze(rep).speedchange = vChange;
                        %find locations that the fly isn't moving post stim start
                        vLow = find(temp(61:121) <= min_speed); 
                        if ~isempty(vLow)
                            vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                            % calculate the specific bits that qualify for 'freezing'
                            vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                            vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                            % find the start and end points of the 'freeze'
                            try % gaps present
                                vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                                vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                                vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                                if vFreeze(rep).start>vFreeze(rep).end
                                    vFreeze(rep).start = vFreeze(rep).end;
                                end
                            catch % no gaps in the 'frozen' data points
                                vFreeze(rep).start = 1;
                                vFreeze(rep).end = length(vFreeze(rep).pos);
                            end
                            vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                            vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                            if vFreeze(rep).idx(1)>windowEnd(cond)
                                vFreeze(rep).idx = [];
                                vFreeze(rep).start = [];
                                vFreeze(rep).end = [];
                            end 
                        end
                    end
                end
            end
            tot = (sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]));
            switch dd
                case 1
                    FreezeData(cond).CONTROL(ifly).data = vFreeze;
                    FreezeData(cond).CONTROL(ifly).totwalking = tot;
                case 2
                    FreezeData(cond).STIM(ifly).data = vFreeze;
                    FreezeData(cond).STIM(ifly).totwalking = tot;
            end
        end
    end
end
       

% Plot the data together:


% Calculate the probabilty of a fly stopping...
MT = max([CONTROL.num.fly, STIM.num.fly]);
MT = NaN(MT,1);

for cond = 1:7
    for dd = 1:2
        temp = []; output = [];
        switch dd
            case 1
                temp = FreezeData(cond).CONTROL;
                num = CONTROL.num;
                output.raw = zeros(6,num.fly);
            case 2
                temp = FreezeData(cond).STIM;
                num = STIM.num;
                output.raw = zeros(6,num.fly);
        end
        for ifly = 1:num.fly
            for rep = 1:6
                if ~isempty(temp(ifly).data(rep).duration)
                    % find the start place
                    output.raw(rep,ifly) = temp(ifly).data(rep).pos(1);
                end
            end
            output.totwalking(ifly) = temp(ifly).totwalking;
        end
        % pull together data across flies
        output.totalpeak = sum(output.totwalking);
        a = (output.raw>0); % find all instances of freezing
        output.totalfreeze =  sum(sum(a));
        output.maxfreeze = output.totalfreeze/output.totalpeak;
        output.inc = output.maxfreeze/output.totalfreeze;
        output.loc = output.raw(a);
        
        % find the probability
        a = sort(output.loc);
        a = a-60;
        x = []; y = [];
        for ii = 1:length(a)
           x(ii) = a(ii); 
           y(ii) = ii*output.inc;
        end
        output.x = [0, x, 60];
        output.y = [0, y, output.maxfreeze];
        switch dd
            case 1
                cntl = output;           
            case 2
                stim = output;
        end
    end  
    %     fig = getfig;
    %     hold all
    %     plot(cntl.x, cntl.y, 'color', 'k')
    %     plot(stim.x, stim.y, 'color', 'r')
    FreezeData(cond).C_output = cntl;
    FreezeData(cond).S_output = stim;
end
   
controlColor = Color('black', 'lightgrey', 6);
controlColor = [0 0 0; controlColor];
% InColor = {'000000', '9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'};
InColor = {'003400', '004D00', '006700', '178017', '319A31', '4AB34A', '64CD64'}; % gtacr1 greens
for ii = 1:length(InColor)
    ColorList(ii,:) = hex2rgb(InColor{ii});
end


% Plot the data together:
LW = 2;
fig = getfig; hold all

for cond = 2:7
    plot(FreezeData(cond).C_output.x, FreezeData(cond).C_output.y, 'color', controlColor(cond,:), 'linewidth', LW);
    plot(FreezeData(cond).S_output.x, FreezeData(cond).S_output.y, 'color', ColorList(cond,:), 'linewidth', LW);
end
xlim([0, 60])
ylim([0, 0.9])
set(gca,'TickDir','out');
xlabel('Time (fps)')
ylabel('Percent of freezes (%)')
figname = [comp ' freezing probabilty curve'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);







% interpolated version

controlColor = Color('black', 'lightgrey', 7);
controlColor = [0 0 0; controlColor];
InColor = {'f37b4f', '9c3009', 'c33c0c', 'f15a24', 'f37b4f', 'f69b7b', 'f9bca7', 'fbd5c8'};
for ii = 1:length(InColor)
    ColorList(ii,:) = hex2rgb(InColor{ii});
end

 
fig = getfig; hold all
for dd = 1:2
    switch dd
        case 1 
            kolor = controlColor;
        case 2
            kolor = ColorList;
    end
    % derivative plot
    for cond = 1:7
        switch dd
            case 1
                x = [FreezeData(cond).C_output.x(2:end), 61]; % 
                y = [FreezeData(cond).C_output.y(2:end), FreezeData(cond).C_output.maxfreeze];% 
            case 2
                x = [FreezeData(cond).S_output.x(2:end), 61]; % 
                y = [FreezeData(cond).S_output.y(2:end), FreezeData(cond).S_output.maxfreeze];% 
        end
       if cond == 1
            LW = 3;
       else
            LW = 1;
       end
        xq = 1:61; xq(x) = [];
        ix = find(diff(x)==0); % where is there a repeat number
        % run each section...
        sloc = 1; output = []; raw =[];
        raw = [x',y'];
    %     figure; hold all;
        for ii = 1:(length(ix)+1)
            try
                loc = (ix(ii));
            catch
                loc = length(x); 
            end
            x1 = x(sloc:loc);
            y1 = y(sloc:loc);
            tmin = xq>=x1(1); tmax = (xq<=x1(end));
            xq1 = xq(all([tmin; tmax]));
            yq = interp1(x1,y1,xq1); %,'spline'
    %         plot(x1,y1); plot(xq1, yq)      
            output(ii).yq = yq;
            output(ii).xq1 = xq1;
            sloc = loc+1;
            raw = [raw; output(ii).xq1', output(ii).yq'];
        end
        % sort the order 
        [~, isort] = sort(raw(:,1));
        temp = [];
        temp = raw(isort,:); temp = [0,0; temp];
        raw = smooth(temp(:,2),1);
%         % plot the data
%         subplot(1,2,1); hold all
%         a = smooth(diff([0;raw]), 3);
%         plot(temp(:,1), a , 'color', kolor(cond,:)) %
% %         plot(temp(:,1), temp(:,2), 'color', kolor(cond,:))
%         subplot(1,2,2); hold all
        plot(temp(:,1), raw, 'color', kolor(cond,:), 'linewidth', LW)
    end
end
xlim([0, 61])
    
set(gca,'TickDir','out');
xlabel('Time (fps)')
ylabel('Cumulative freezing (%)')
figname = [comp ' cumulative distribution of freezing curve'];
title({figname;' ';' '})

save_figure(fig, [fig_dir, figname]);





%% Percent of walking flies that freeze (by fly) 
type = 'walking';
min_speed = 0.5;

% percent freezing
stppoint = 0.400;
windowEnd = (60+stppoint*num.fps)*ones(1,7);

condFrames = round(condLength*num.fps);
sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for cond = 1:7
    for ifly = 1:num.fly
        for rep = 1:6
            vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                  'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
            raw = data(cond).(type).raw(ifly).data(:,rep); 
            % if the data fits the behavior type, analyze
            if ~isnan(raw(1))
                temp = smooth(raw, 'rlowess');
                vStart = nanmean(raw(54:60)); %value prestim
%                 vChange = NaN;
                
                % only include trials with greater than min speed
                if vStart >= min_speed
                    % Change in speed
                    vStim = nanmean(raw(76:82)); %value prestim
                    vChange = vStim-vStart;
                    vFreeze(rep).speedchange = vChange;
                    %find locations that the fly isn't moving post stim start
                    vLow = find(temp(61:121) <= min_speed); 
                    if ~isempty(vLow)
                        vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                        % calculate the specific bits that qualify for 'freezing'
                        vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                        vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                        % find the start and end points of the 'freeze'
                        try % gaps present
                            vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                            vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                            vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                            if vFreeze(rep).start>vFreeze(rep).end
                                vFreeze(rep).start = vFreeze(rep).end;
                            end
                        catch % no gaps in the 'frozen' data points
                            vFreeze(rep).start = 1;
                            vFreeze(rep).end = length(vFreeze(rep).pos);
                        end
                        vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                        vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                        if vFreeze(rep).idx(1)>windowEnd(cond)
                            vFreeze(rep).idx = [];
                            vFreeze(rep).start = [];
                            vFreeze(rep).end = [];
                        end 
                    end
                end
            end
        end
        FreezeData(cond).all(ifly).data = vFreeze;
    end
end


MT = NaN(6,1);
for cond = 1:7   
    %empty data here -- so we don't factor in zeros
    FreezeData(cond).avgchange = NaN(1,num.fly);
    FreezeData(cond).avgdelay = NaN(1,num.fly);
    for ifly = 1:num.fly
        inputdata = FreezeData(cond).all(ifly).data;   
        output = struct('freezeresponse', MT, 'delay', MT, 'duration', MT, 'change', MT);
        for rep = 1:6
           idx = inputdata(rep).idx;
           % did the fly freeze on that trial?
           if isempty(idx)
             output.freezeresponse(rep) = false;
           else
             output.freezeresponse(rep) = true;  
             output.delay(rep) = (idx(1)-60)/num.fps;
             output.duration(rep) = length(idx)/num.fps;
             output.change(rep) = inputdata(rep).speedchange;
           end
        end
        tot = (sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]));
        if tot >= 2 %must have a possibility of 2 trials for each fly to include
            output.freezecount = sum(output.freezeresponse);
            output.totalcount = sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]);
            output.avgdelay = nanmean(output.delay);
            output.avgduration = nanmean(output.duration);
            output.avgchange = nanmean(output.change);
        else
            output.freezecount = NaN;
            output.totalcount = NaN;
            output.avgdelay = NaN;
            output.avgduration = NaN;
            output.avgchange = NaN;
        end    
            FreezeData(cond).data(ifly) = output;
            FreezeData(cond).totalcount(ifly) = output.totalcount;
            FreezeData(cond).freezecount(ifly) = output.freezecount;
            FreezeData(cond).percent(ifly) = FreezeData(cond).freezecount(ifly)/FreezeData(cond).totalcount(ifly);
            FreezeData(cond).avgdelay(ifly) = output.avgdelay;
            FreezeData(cond).avgduration(ifly) = output.avgduration;
            FreezeData(cond).avgchange(ifly) = output.avgchange;
    end
end

% -- line plot the frequency of freezing -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(FreezeData(cond).percent);
    y.err(cond) = sem(FreezeData(cond).percent,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([0, 1]);
set(gca,'TickDir','out');
title([filename(1:end-4) ' ' type ' freeze frequency'])
xlabel('Stimulation length (sec)')
ylabel('Frequency of freezing (%)')

% save the data
figname  = [fig_dir, structure_name ' ' type ' freezing lineplot endpoint ' num2str(stppoint)];
save_figure(fig, figname);



%---scatter plot Percent of flies that freeze
sz = 70;
Kolor = colorList(5,:);
fig = getfig; hold all
xlength = 0.2; xlim([0,8])
for cond = 1:7
    y = FreezeData(cond).percent;
    y(isnan(y)) = [];
    idx = length(y);
    x = cond:xlength/idx:cond+(xlength-xlength/idx);
    scatter(x, y, sz, Kolor, 'filled')
    avg = nanmean(y);
    err = sem(y,2);
%     errorbar(cond+(xlength/2), avg, err, 'Color', 'r')
    plot([cond-(xlength/2), cond+(xlength*1.5)], [avg,avg], 'Color', 'k', 'linewidth', 3)
%     errorbar(cond, avg, err, 'color', 'k')
end
title({[structure_name ' Frequency of freezing behavior. Window end: ' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Percent of trials freezing (%)')
figname  = [fig_dir, structure_name ' ' type ' percent freezing endpoint ' num2str(stppoint)];
set(gca,'TickDir','out');

save_figure(fig, figname);

%% Delay to freezing in walking flies (by fly) 

condFrames = round(condLength*num.fps);
% Delay data
stppoint = 1.5;
windowEnd = (60+stppoint*num.fps)*ones(1,7);
% stppoint = 0.5;
% windowEnd = condFrames + 60 + round(stppoint*num.fps);
% windowEnd(1) = windowEnd(end);

condFrames = round(condLength*num.fps);
sSpan = 6;% 1/30 is the frame rate, so if we want a rolling 200ms avg, we want a
% smooth of ...6 seconds
for cond = 1:7
    for ifly = 1:num.fly
        for rep = 1:6
            vFreeze(rep) = struct('idx',[],'start',[],'end',[],'pos',[],'spacing',[],...
                                  'gaps',[],'singlespacing',[],'duration',[],'speedchange', []);
            raw = data(cond).(type).raw(ifly).data(:,rep); 
            % if the data fits the behavior type, analyze
            if ~isnan(raw(1))
                temp = smooth(raw, 'rlowess');
                vStart = nanmean(raw(54:60)); %value prestim
%                 vChange = NaN;
                
                % only include trials with greater than min speed
                if vStart >= min_speed
                    % Change in speed
                    vStim = nanmean(raw(76:82)); %value prestim
                    vChange = vStim-vStart;
                    vFreeze(rep).speedchange = vChange;
                    %find locations that the fly isn't moving post stim start
                    vLow = find(temp(61:121) <= min_speed); 
                    if ~isempty(vLow)
                        vFreeze(rep).pos = vLow+60; %position of 'freeze data' in matrix
                        % calculate the specific bits that qualify for 'freezing'
                        vFreeze(rep).spacing = diff(vFreeze(rep).pos);
                        vFreeze(rep).gaps = find(vFreeze(rep).spacing>1);
                        % find the start and end points of the 'freeze'
                        try % gaps present
                            vFreeze(rep).end = vFreeze(rep).gaps(1); % stop when the data pops into moving speeds
                            vFreeze(rep).singlespacing = find(vFreeze(rep).spacing == 1); %consecutive freezing points         
                            vFreeze(rep).start = vFreeze(rep).singlespacing(1); %start of 2+ points below moving
                            if vFreeze(rep).start>vFreeze(rep).end
                                vFreeze(rep).start = vFreeze(rep).end;
                            end
                        catch % no gaps in the 'frozen' data points
                            vFreeze(rep).start = 1;
                            vFreeze(rep).end = length(vFreeze(rep).pos);
                        end
                        vFreeze(rep).duration = vFreeze(rep).end - vFreeze(rep).start +1; 
                        vFreeze(rep).idx = vFreeze(rep).pos(vFreeze(rep).start:vFreeze(rep).end);
                        if vFreeze(rep).idx(1)>windowEnd(cond)
                            vFreeze(rep).idx = [];
                            vFreeze(rep).start = [];
                            vFreeze(rep).end = [];
                        end 
                    end
                end
            end
        end
        FreezeData(cond).all(ifly).data = vFreeze;
    end
end


MT = NaN(6,1);
for cond = 1:7   
    %empty data here -- so we don't factor in zeros
    FreezeData(cond).avgchange = NaN(1,num.fly);
    FreezeData(cond).avgdelay = NaN(1,num.fly);
    for ifly = 1:num.fly
        inputdata = FreezeData(cond).all(ifly).data;   
        output = struct('freezeresponse', MT, 'delay', MT, 'duration', MT, 'change', MT);
        for rep = 1:6
           idx = inputdata(rep).idx;
           % did the fly freeze on that trial?
           if isempty(idx)
             output.freezeresponse(rep) = false;
           else
             output.freezeresponse(rep) = true;  
             output.delay(rep) = (idx(1)-60)/num.fps;
             output.duration(rep) = length(idx)/num.fps;
             output.change(rep) = inputdata(rep).speedchange;
           end
        end
        tot = (sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]));
        if tot >= 2 %must have a possibility of 2 trials for each fly to include
            output.freezecount = sum(output.freezeresponse);
            output.totalcount = sum([group(ifly).(type)(cond,:),group(ifly).(type)(cond+7,:)]);
            output.avgdelay = nanmean(output.delay);
            output.avgduration = nanmean(output.duration);
            output.avgchange = nanmean(output.change);
        else
            output.freezecount = NaN;
            output.totalcount = NaN;
            output.avgdelay = NaN;
            output.avgduration = NaN;
            output.avgchange = NaN;
        end    
            FreezeData(cond).data(ifly) = output;
            FreezeData(cond).totalcount(ifly) = output.totalcount;
            FreezeData(cond).freezecount(ifly) = output.freezecount;
            FreezeData(cond).percent(ifly) = FreezeData(cond).freezecount(ifly)/FreezeData(cond).totalcount(ifly);
            FreezeData(cond).avgdelay(ifly) = output.avgdelay;
            FreezeData(cond).avgduration(ifly) = output.avgduration;
            FreezeData(cond).avgchange(ifly) = output.avgchange;
    end
end


% -- line plot delay to freeze -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(FreezeData(cond).avgdelay);
    y.err(cond) = sem(FreezeData(cond).avgdelay,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([0, 1]);
title({[structure_name ' Freezing Delay. End point: stimend+' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Time to freeze initiation (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, structure_name ' ' type ' freezing delay linegraph ' num2str(stppoint)];

% save the data
figname  = [fig_dir, structure_name ' ' type ' freezing lineplot endpoint ' num2str(stppoint)];
save_figure(fig, figname);





% Time delay for freezing
fig = getfig; hold all
xlength = 0.2; xlim([0,8])
for cond = 1:7
    y = FreezeData(cond).avgdelay;
    y(isnan(y)) = [];
    idx = length(y);
    x = cond:xlength/idx:cond+(xlength-xlength/idx);
    scatter(x, y, sz, Kolor, 'filled')
    avg = nanmean(y);
    err = sem(y,2);
%     errorbar(cond+(xlength/2), avg, err, 'Color', 'r')
    plot([cond-(xlength/2), cond+(xlength*1.5)], [avg,avg], 'Color', 'k', 'linewidth', 3)
%     errorbar(cond, avg, err, 'color', 'k')
end
title({[structure_name ' Freezing Delay. End point: stimend+' num2str(stppoint)];' ';' '})
xlabel('Condition')
ylabel('Time to freeze initiation (s)')
set(gca,'TickDir','out');
figname  = [fig_dir, structure_name ' ' type ' freezing delay endpoint ' num2str(stppoint)];
ylim([0, 1.5])

save_figure(fig, figname);


%% Change in speed after stimulus onset

condFrames = round(condLength*num.fps);
% Delay data
stppoint = 1.5;
windowEnd = (60+stppoint*num.fps)*ones(1,7);
% stppoint = 0.5;
% windowEnd = condFrames + 60 + round(stppoint*num.fps);
% windowEnd(1) = windowEnd(end);

condFrames = round(condLength*num.fps);


controlRange = 54:60;
stimRange = 76:82;
% %0.5
% controlRange = 45:60;
% stimRange = 75:90;

for cond = 1:7
    for ifly = 1:num.fly
            cntl = data(cond).walking.raw(ifly).data(controlRange,:);
            stim = data(cond).walking.raw(ifly).data(stimRange,:);
            speed(cond).raw(:,ifly) = nanmean(stim,1) - nanmean(cntl,1);
            speed(cond).flyavg(ifly) = nanmean(speed(cond).raw(:,ifly));
%             % percent change
%             speed(cond).per_raw(:,ifly) = (nanmean(stim,1) - nanmean(cntl,1))./nanmean(cntl,1);
%             speed(cond).per_avg(ifly) = nanmean(speed(cond).per_raw(:,ifly));
    end
    speed(cond).avg = nanmean(speed(cond).flyavg);
    speed(cond).err = sem(speed(cond).flyavg,2);
end
    
%absolute change   
fig = getfig; hold all;
ylim([-1.2,0.4])
xlength = 0.2; xlim([0,8])
for cond = 1:7
    inpt = speed(cond).flyavg;
    loc = ~isnan(inpt);
    idx = sum(loc);
    x = cond:xlength/idx:cond+(xlength-xlength/idx);
    y = inpt(loc);
    scatter(x, y, sz, Kolor, 'filled')
    avg = nanmean(inpt);
    plot([cond, cond+xlength], [avg,avg], 'Color', 'k', 'linewidth', 3)
end
title([structure_name ' change in speed from stim '])
xlabel('Condition')
ylabel('Change in speed (cm/sec)')
hline(0,'k:')

set(gca,'TickDir','out');

figname  = [fig_dir, structure_name ' change in speed after stimulus'];
save_figure(fig, figname);



% -- line plot delay to freeze -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(speed(cond).flyavg);
    y.err(cond) = sem(speed(cond).flyavg,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([-.80, 0]);
title([structure_name ' change in speed from stim '])
xlabel('Condition')
ylabel('Change in speed (cm/sec)')
hline(0,'k:')

figname  = [fig_dir, structure_name ' change in speed line graph'];

% save the data;
save_figure(fig, figname);

%% Time to half-speed drop in walking flies
num.fps = 30;
stppoint = 0.2;
condFrames = round(condLength*num.fps);
windowEnd = condFrames + 60 + round(stppoint*num.fps);
windowEnd(1) = windowEnd(end);

controlRange = 55:61;
% stimRange = 61:88;
% x = cond:xlength/idx:cond+(xlength-xlength/idx);
for ii = 1:2
    switch ii
        case 1 %control
            num = CONTROL.num;
        case 2 %stim (10B)
            
    end
for cond = 1:7
    stimRange = [61:windowEnd(cond)];
    for ifly = 1:num.fly
            cntl = nanmean(Cdata(cond).walking.raw(ifly).data(controlRange,:),2);
            stim = nanmean(Cdata(cond).walking.raw(ifly).data(stimRange,:),2);
            if ~isnan(cntl)
                % min value:
                [yMin, loc] = (min(stim));
                MinLoc = loc+60;
                HLoc = round(loc/2)+60;
                Cdata(cond).walking.MinLoc(ifly) = MinLoc;
                Cdata(cond).walking.HLoc(ifly) = HLoc;

%                 %PLOT preview
%                 fig = getfig; hold all
%                 plot(x, nanmean(data(cond).walking.raw(ifly).data,2))
%                 vline(x(HLoc))
%                 vline([x(stimRange(1)), x(stimRange(end))], 'b')
%                 vline([x(controlRange(1)), x(controlRange(end))], 'm')
%                 scatter(x(stimRange), nanmean(data(cond).walking.raw(ifly).data(stimRange,:),2))
%                 scatter(x(MinLoc), yMin, 60, 'r', 'filled')
            else
                Cdata(cond).walking.MinLoc(ifly) = NaN;
                Cdata(cond).walking.HLoc(ifly) = NaN;
            end
    end
end
end
   

%PLOT condition average
fig = getfig;
    hold all
    cond = 1;
    plot(x, data(cond).walking.AVG, 'k')
    loc = round(nanmean(data(cond).walking.HLoc));
    scatter(x(loc), data(cond).walking.AVG(loc), 60, 'r', 'filled')
    
    for cond = 3:2:7
    plot(x, data(cond).walking.AVG, 'b')
    loc = round(nanmean(data(cond).walking.HLoc));
    scatter(x(loc), data(cond).walking.AVG(loc), 60, 'b', 'filled')
    end



%absolute change   
fig = getfig; hold all;
ylim([-1.2,0.4])
xlength = 0.2; xlim([0,8])
for cond = 1:7
    inpt = speed(cond).flyavg;
    loc = ~isnan(inpt);
    idx = sum(loc);
    x = cond:xlength/idx:cond+(xlength-xlength/idx);
    y = inpt(loc);
    scatter(x, y, sz, Kolor, 'filled')
    avg = nanmean(inpt);
    plot([cond, cond+xlength], [avg,avg], 'Color', 'k', 'linewidth', 3)
end
title([structure_name ' change in speed from stim '])
xlabel('Condition')
ylabel('Change in speed (cm/sec)')
hline(0,'k:')

set(gca,'TickDir','out');

figname  = [fig_dir, structure_name ' change in speed after stimulus'];
save_figure(fig, figname);



% -- line plot delay to freeze -- %
Kolor = colorList(5,:);
x = []; y = [];
x = condLength;
for cond = 1:7
    y.avg(cond) = nanmean(speed(cond).flyavg);
    y.err(cond) = sem(speed(cond).flyavg,2);
end
fig = getfig; hold all
scatter(x, y.avg, 70, Kolor, 'filled',...
    'markeredgecolor', 'none', 'markerfacecolor', 'flat')
errorbar(x, y.avg, y.err, 'linestyle', 'none', 'color', Kolor)
xlim([-.1, .8]); ylim([-.80, 0]);
title([structure_name ' change in speed from stim '])
xlabel('Condition')
ylabel('Change in speed (cm/sec)')
hline(0,'k:')

figname  = [fig_dir, structure_name ' change in speed line graph'];

% save the data;
save_figure(fig, figname);




%% zoom in of drop



condLength = [0, 0.03, 0.06, 0.09, 0.180, 0.320, 0.720];
x =  -2:1/30:2;
sSpan = 3;
for cond = 1:2:7
    fig = getfig; hold all
    for ii = 1:2 
        y = [];
        switch ii
            case 1 % control
                InColor = {'000000', '000000', '000000', '000000', '000000', '000000', '000000'};
%                  InColor = {'301818', '000000', '634b4b', '000000', '866868', '000000', 'a98585'};
                 y.avg = smooth(Cdata(cond).(type).AVG,sSpan);
                 y.err = smooth(Cdata(cond).(type).ERR,sSpan);
            case 2 % stim (10b, etc)
                InColor = {'f15a24', 'f15a24', 'f15a24', 'f15a24', 'f15a24', 'f15a24', 'f15a24'}; 
%                 InColor = {'000000', '000000', '9c3009', '000000', 'f15a24', '000000', 'f9bca7'}; 
                y.avg = smooth(Sdata(cond).(type).AVG,sSpan);
                y.err = smooth(Sdata(cond).(type).ERR,sSpan);
        end
        % COLOR SELECTION %
        for ii = 1:length(InColor)
            colorList(ii,:) =  hex2rgb(InColor{ii});     
        end
%         fill_data = error_fill(x,  y.avg, y.err);
%         h = fill(fill_data.X, fill_data.Y, colorList(cond,:), 'EdgeColor','none');
%         set(h, 'facealpha', 0.2)
        plot(x, y.avg, 'color', colorList(cond,:), 'linewidth', 1)
        plot(x, y.avg+y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
        plot(x, y.avg-y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
    end
    %axes & labels
    ymax = 1.6;
    offset = .03;
    ylim([0 ymax])
    vline([0, condLength(cond)], 'k')
    xlabel('time (sec)')
    ylabel('Speed (cm/s)')
    title([comp ' headless ' type ' cond: ' num2str(cond)])
    set(gca,'TickDir','out');
    figname = [fig_dir, comp ' speed timecourse cond ' num2str(cond)];
    % savefig(fig, figname);
    save_figure(fig, figname); 
end





















%% Intact fly time course figures
% load the data set with first secttion of MIDSTIM_Step_7
% file_root = '/Volumes/Evyn SSD/Evyn UW work/matlabroot';
file_root = 'C:\matlabroot';
fig_dir = [file_root '/Interneuron Lines/Figures/'];
comp = '10B vs SH Chrimson';
type = 'walking'; % stationary

% first time get control data:
for ii = 1:2
    switch ii
        case 1 % CONROL
            group = CONTROL.group;
            fly = CONTROL.fly; 
            num = CONTROL.num;
        case 2 % STIM
            group = STIM.group;
            fly = STIM.fly; 
            num = STIM.num;
    end
    % Pull out the desired data
    for cond = 1:7 
        for ifly = 1:num.fly  
            % select the data and generate a filter for the initial behavior state
            filter = [group(ifly).(type)(cond,:), group(ifly).(type)(cond+7,:)];  
            cw = [fly(ifly).Control.speed(cond).data(1:end-1,:); fly(ifly).Stim.speed(cond).data];
            ccw = [fly(ifly).Control.speed(cond+7).data(1:end-1,:); fly(ifly).Stim.speed(cond+7).data];
            input = [cw, ccw];
            input(:,~filter) = nan; %nix the trials that didn't fit the movement type
            % if there aren't 2 trials min, all go NaN
            if sum(filter)<2
                input(:,filter) = nan;
            end
            %stats on the fly's data
            data(cond).(type).raw(ifly).data = input;
            data(cond).(type).raw(ifly).avg = nanmean(input,2);
            data(cond).(type).raw(ifly).err = sem(input,2);
            % add averages to the group data
            data(cond).(type).avg(:,ifly) = data(cond).(type).raw(ifly).avg;
            data(cond).(type).err(:,ifly) = data(cond).(type).raw(ifly).err;
        end
        % group avg
        data(cond).(type).AVG = nanmean(data(cond).(type).avg,2);
        data(cond).(type).ERR = sem(data(cond).(type).avg,2);
    end
    switch ii
        case 1 % CONROL
            Cdata = data;
        case 2 % STIM
            Sdata = data;
    end
end
clear filter group fly num cw ccw input data

condLength = [0, 0.03, 0.06, 0.09, 0.180, 0.320, 0.720];
x =  -2:1/30:2;
sSpan = 3;
for cond = 1:2:7
    fig = getfig; hold all
    for ii = 1:2 
        y = [];
        switch ii
            case 1 % control
                InColor = {'000000', '000000', '000000', '000000', '000000', '000000', '000000'};
%                  InColor = {'301818', '000000', '634b4b', '000000', '866868', '000000', 'a98585'};
                 y.avg = smooth(Cdata(cond).(type).AVG,sSpan);
                 y.err = smooth(Cdata(cond).(type).ERR,sSpan);
            case 2 % stim (10b, etc)
                InColor = {'f15a24', 'f15a24', 'f15a24', 'f15a24', 'f15a24', 'f15a24', 'f15a24'}; 
%                 InColor = {'000000', '000000', '9c3009', '000000', 'f15a24', '000000', 'f9bca7'}; 
                y.avg = smooth(Sdata(cond).(type).AVG,sSpan);
                y.err = smooth(Sdata(cond).(type).ERR,sSpan);
        end
        % COLOR SELECTION %
        for ii = 1:length(InColor)
            colorList(ii,:) =  hex2rgb(InColor{ii});     
        end
%         fill_data = error_fill(x,  y.avg, y.err);
%         h = fill(fill_data.X, fill_data.Y, colorList(cond,:), 'EdgeColor','none');
%         set(h, 'facealpha', 0.2)
        plot(x, y.avg, 'color', colorList(cond,:), 'linewidth', 1)
        plot(x, y.avg+y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
        plot(x, y.avg-y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
    end
    %axes & labels
    ymax = 1.6;
    offset = .03;
    ylim([0 ymax])
    vline([0, condLength(cond)], 'k')
    xlabel('time (sec)')
    ylabel('Speed (cm/s)')
    title([comp ' headless ' type ' cond: ' num2str(cond)])
    set(gca,'TickDir','out');
    figname = [fig_dir, comp ' speed timecourse cond ' num2str(cond)];
    % savefig(fig, figname);
    save_figure(fig, figname); 
end





% PLOT THE DATA %


fig = getfig;
hold all
%CONTROL
y.avg = smooth(data(1).(type).AVG,sSpan);
y.err = smooth(data(1).(type).ERR,sSpan);
% fill_data = error_fill(x,  y.avg, y.err);
% h = fill(fill_data.X, fill_data.Y, Color('Black'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
plot(x, y.avg, 'color', Color('Black'), 'linewidth', 1)
plot(x, y.avg+y.err, 'color', Color('Black'), 'linewidth', 0.5)
plot(x, y.avg-y.err, 'color', Color('Black'), 'linewidth', 0.5)
%STIM % 60, 180, 720 = 3,5,7;
for cond = [3,5,7]
    y.avg = smooth(data(cond).(type).AVG,sSpan);
    y.err = smooth(data(cond).(type).ERR,sSpan);
%     fill_data = error_fill(x, y.avg, y.err);
%     h = fill(fill_data.X, fill_data.Y, colorList(cond,:), 'EdgeColor','none');
%     set(h, 'facealpha', 0.2)
    plot(x, y.avg, 'color', colorList(cond,:), 'linewidth', 1)
    plot(x, y.avg+y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
    plot(x, y.avg-y.err, 'color', colorList(cond,:), 'linewidth', 0.5)
end
%axes & labels
ymax = 1.6;
offset = .03;
ylim([0 ymax])
LX = [0,condLength(3); 0, condLength(5); 0, condLength(7)];
for ii = 1:3
    LY = [ymax-(offset*ii), ymax-(offset*ii)];
    plot(LX(ii,:), LY, 'Color', 'g', 'linewidth', 3)
end
% vline([0, 0, 0.6, 0.18, 0.72], 'g')
xlabel('time (sec)')
ylabel('Speed (cm/s)')
title([filename(1:end-4) ' headless ' type])
set(gca,'TickDir','out');
% save the data
figname = [fig_dir, filename(1:end-4) ' ' type ' speed timecourse'];
% savefig(fig, figname);
save_figure(fig, figname); 







%% Headless fly leg movement analysis
file_root = '/Volumes/Evyn SSD/Evyn UW work/matlabroot';
file_name = '13B-20847-csChrimson-7V-offball tracking data.mat';
load([file_root, '/', file_name])

Joints = {'CoFe', 'FeTi', 'TiTa'};

%Build empty data structure & fill with joint angle data;
for ifly = 1:length(tracking)
  for iJ = 1:3 %for each joint measured
    for cond = 1:14
        idx = 0;
        data(ifly).(Joints{iJ}).raw = NaN(337,28);
        for iFrame = 90:426 %specific time period analyzed
            input = tracking(ifly).Left_Front(cond,1).frame(iFrame).(Joints{iJ});
            idx = idx+1;
            if sum(input) == 0 %carry NaNs through
                input = NaN;
            end
            data(ifly).(Joints{iJ}).raw(idx,cond) = input;
        end
        % fill in the gaps with interpolated data:
        xstart = 1:337;
        vstart = data(ifly).(Joints{iJ}).raw(:,cond);
        nanloc = isnan(vstart);
        xq = find(nanloc==true);
        v = vstart(~nanloc);
        x = xstart(~nanloc);
        vq = interp1(x,v,xq); 
        allV = vstart;
        allV(nanloc) = vq;
        data(ifly).(Joints{iJ}).all(:,cond) = allV;
    end
    % % Find the average across all trials for a given fly% %
    % isolate the stim trials
    data(ifly).(Joints{iJ}).stim.all = data(ifly).(Joints{iJ}).all;
    data(ifly).(Joints{iJ}).stim.all(:,1) = [];
    data(ifly).(Joints{iJ}).stim.all(:,8) = [];
    % isolate the control trials
    data(ifly).(Joints{iJ}).control.all = ...
        [data(ifly).(Joints{iJ}).all(:,1), data(ifly).(Joints{iJ}).all(:,8)];
    % find the average
    data(ifly).(Joints{iJ}).stim.avg = nanmean(data(ifly).(Joints{iJ}).stim.all,2);
    data(ifly).(Joints{iJ}).control.avg = nanmean(data(ifly).(Joints{iJ}).control.all,2);
    % find the error
    data(ifly).(Joints{iJ}).stim.err = sem(data(ifly).(Joints{iJ}).stim.all,2,2);
    data(ifly).(Joints{iJ}).control.err = sem(data(ifly).(Joints{iJ}).control.all,2,2);
    
  end
end

% DEMO of the data from 'data' within a fly for all trials
% figure;
% hold all
% % stim
% plot(xstart,data(ifly).(Joints{iJ}).stim.all, 'r', 'linewidth', 0.5)
% plot(xstart,data(ifly).(Joints{iJ}).stim.avg, 'r', 'linewidth', 2)
% % control
% plot(xstart,data(ifly).(Joints{iJ}).control.all, 'k', 'linewidth', 0.5)
% plot(xstart,data(ifly).(Joints{iJ}).control.avg, 'k', 'linewidth', 2)
% vline([60,317],'g-')
% 
% fill_data = error_fill(xstart, data(ifly).(Joints{iJ}).stim.avg, data(ifly).(Joints{iJ}).stim.err);
% h = fill(fill_data.X, fill_data.Y, get_color('red'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)
% 
% fill_data = error_fill(xstart, data(ifly).(Joints{iJ}).control.avg, data(ifly).(Joints{iJ}).control.err);
% h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
% set(h, 'facealpha', 0.2)



type = {'stim', 'control'};
% Plot a figure of all the joint angles over time -- line/fly
for itype = 1:2
   for iJ = 1:3
      for ifly = 1:length(tracking) 
        inputdata.(Joints{iJ}).(type{itype}).all(:,ifly) = data(ifly).(Joints{iJ}).(type{itype}).avg;
      end
      inputdata.(Joints{iJ}).(type{itype}).avg = mean(inputdata.(Joints{iJ}).(type{itype}).all,2);
      inputdata.(Joints{iJ}).(type{itype}).err = sem(inputdata.(Joints{iJ}).(type{itype}).all,2);
   end
end

colorList = {'Navy', 'DodgerBlue', 'SkyBlue', 'Black', 'DimGrey', 'DarkGrey'};

fig = getfig('',1); 
for iJ = 1:3
%     subplot(3,1,iJ); 
    hold all
    %STIM
    fill_data = error_fill(xstart, inputdata.(Joints{iJ}).stim.avg,...
                           inputdata.(Joints{iJ}).stim.err);
    h = fill(fill_data.X, fill_data.Y, get_color(colorList{iJ}), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(xstart, inputdata.(Joints{iJ}).stim.avg, 'color', get_color(colorList{iJ}), 'linewidth', 3)
    %CONTROL
    fill_data = error_fill(xstart, inputdata.(Joints{iJ}).control.avg,...
                           inputdata.(Joints{iJ}).control.err);
    h = fill(fill_data.X, fill_data.Y, get_color(colorList{iJ+3}), 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(xstart, inputdata.(Joints{iJ}).control.avg, 'color', get_color(colorList{iJ+3}), 'linewidth', 3)
    vline([60,317],'g-')
end
xlabel('Time')
ylabel('Joint angle (deg)')
title([file_name(1:end-18) ' headless'])







figure;
hold all
for cond = [2:7,9:14]
    plot(xstart, data(ifly).(Joints{iJ}).all(:,cond),'r','linewidth', 1)
end
for cond = [1,8]
    plot(xstart, data(ifly).(Joints{iJ}).all(:,cond),'k','linewidth', 1)
end

vline([60,317],'g-')



clear data
xstart = 1:337;
vstart = data(ifly).(Joints{iJ}).all(:,cond);
nanloc = isnan(vstart);
xq = find(nanloc==true);
v = vstart(~nanloc);
x = xstart(~nanloc);
vq = interp1(x,v,xq,'spline'); 
allV = vstart;
allV(nanloc) = vq;

figure;
hold all
scatter(x,v,'k') % og points
scatter(xq,vq,'r') %new points
plot(xstart, allV, 'k', 'linewidth', 2)

% plot([x;xq], [v;vq])
vline([60,317],'g-')



x = [1 2 5];
v = [5 6 9];
xq = [3,4];


vq = interp1(x,v,xq);
disp(vq)















%% Find delay to movement time in stationary flies!
% how long after the stim onset does the fly reach 0.3 cm/sec of movement?
% 
















%% COMPRESS_AVI Compresses avi files 

% input folder
dir_input =  'E:\Basler Trig\Intersegmental CS activation\';

% output folder
dir_root = 'G:\My Drive\Evyn\Intersegmental CS activation\';
mkdir(dir_root)

% Get list of folders
list_folders = dir(dir_input);

nVideos = 0;
% for each folder within the input folder find the num vids and create
% replacement folders
for iSession = 3:length(list_folders)
    list_videos = dir([dir_input,list_folders(iSession).name,'\*.avi']);
    nVideos = nVideos+length(list_videos);
    dir_output = [dir_root, list_folders(iSession).name '\']; 
    mkdir(dir_output)
end

% Create waitbar
h = waitbar(0,'Compressing videos');


% Iterate through sessions
nVideo = 1;
for iSession = 3:length(list_folders) 
    % Get list of videos in current session
    list_videos = dir([dir_input,list_folders(iSession).name,'\*.avi']);

    % Designate session folder in output directory
    dir_output = [dir_root, list_folders(iSession).name '\']; 
    
    % Iterate through videos
    for iVideo = 1:length(list_videos)
        % Create objects to read and write the video, and open the AVI file
        % for writing
        disp(['Compressing ',list_videos(iVideo).name,' ...'])
        reader = VideoReader([dir_input,list_folders(iSession).name,'\',list_videos(iVideo).name]);        
        writer = VideoWriter([dir_output,list_videos(iVideo).name],'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
        writer.FrameRate = reader.FrameRate;
        writer.Quality = 75;
        open(writer);
        
        % Read and write each frame
        while hasFrame(reader)
            img = readFrame(reader);
            writeVideo(writer,img);
        end
        writeVideo(writer,img);
        close(writer);           
        disp('Done!')
        
        % Update waitbar
        waitbar(nVideo/nVideos,h)
        nVideo = nVideo+1;
    end
end  
    
close(h)
clc





















%%




% create the raw data and add to fly structure:
for ifly = 1:length(fly)
figures_dir = ['C:\matlabroot\Tony Paper Work\' fly(ifly).structure_name '\'];
% structure_name = fly(ifly).structure_name;

    for pp = 1:1
        clear trial LA
        switch pp
            case 1
                param = 'speed';
                InputData = fly(ifly).Speed;
            case 2
                param = 'rotvelocity';
                InputData = fly(ifly).Rotvel;
        end
        shorttime = 60:(60+round(.09*num.fps)+(duration*num.fps));
        longtime = 60:(60+round(.72*num.fps)+(duration*num.fps));
        timerange(1).data = shorttime;
        timerange(2).data = longtime;
        controltime = [60-(.2*num.fps):59];
        mean_speed.control = nanmean(nanmean(InputData(1).y(controltime, :)));
        % find the change in speed from baseline
        trial.speed.short = InputData(2).y(shorttime, :)-mean_speed.control;
        trial.speed.long = InputData(3).y(longtime, :)-mean_speed.control;
        %remove nans:
        trial.speed.short(:,isnan(trial.speed.short(1,:))) = [];
        trial.speed.long(:,isnan(trial.speed.long(1,:))) = [];
        trial.num.short = size(trial.speed.short,2);
        trial.num.long = size(trial.speed.long,2);
        fly(ifly).bootstrap = trial;
    end
end

%% Run the bootstrap on the distribution of speeds:
% Group flydata for all data
tic
icount = 0;
N = 10E3;
hfig = getfig;
type = {'short', 'long'};
for ifly = 1:6
    for itype = 1:2 %short|long
      icount = icount+1;
      err = [];
       if ifly <=3 
           iflycontrol = 7;
       else 
           iflycontrol = 8;
       end
    
    data = [];

%     f = getfig('', 1); hold all
    idx = 0;
    data = [fly(iflycontrol).bootstrap.speed.(type{itype}), fly(ifly).bootstrap.speed.(type{itype})];

    Asize = size(data,2);
    Csize = fly(iflycontrol).bootstrap.num.(type{itype});
    Ssize = fly(ifly).bootstrap.num.(type{itype});
%     
%     poolobj = gcp;
%     addAttachedFiles(poolobj,{'sort'})
%     
    

for n = 1:N
        test = []; sorted = []; a = [];
        loc = randi(Asize,[1,Csize]);
        a = data(:,loc);
        test.control.raw = a;
%         test.control.x.raw = sort(reshape(a,[numel(a),1]));
        a = reshape(a,[numel(a),1]);
         test.control.x.raw = sort(a);
        test.control.y.raw = (1:numel(a))/numel(a);

        loc = randi(size(data,2),[1,Ssize]);
        a = data(:,loc);
        test.stim.raw = a;
        test.stim.x.raw = sort(reshape(a,[numel(a),1]));
        test.stim.y.raw = (1:numel(a))/numel(a);

        %find the min value and max value in the two distrutions 
        minval = min([test.stim.x.raw(1), test.control.x.raw(1)]);
        maxval = max([test.stim.x.raw(end), test.control.x.raw(end)]);

        % add ones and zeros to the end of the vectors to bring to same x value
        test.stim.x.raw = [minval; test.stim.x.raw; maxval];
        test.stim.y.raw = [0; test.stim.y.raw'; 1];
        test.control.x.raw = [minval; test.control.x.raw; maxval];
        test.control.y.raw = [0; test.control.y.raw'; 1];

        % merge to find common x values:
        a = union(test.control.x.raw, test.stim.x.raw);
        sorted.control.y = zeros(length(a),1);
        for ii = 1:length(a)
            x = find(test.control.x.raw>= a(ii)); %all x values larger than the one one in question
            sorted.control.y(ii,1) = test.control.y.raw(x(1));
            xx = find(test.stim.x.raw>= a(ii)); %all x values larger than the one one in question
            sorted.stim.y(ii,1) = test.stim.y.raw(xx(1));
        end
        
        test.diff = abs(sorted.control.y-sorted.stim.y);
        [bval, mdiff] = max(test.diff);
        rawdiff = (sorted.control.y-sorted.stim.y);
        test.DIFF = rawdiff(mdiff);
        
        
%         %update figure
%         if  mod(n,100) == 0
%             idx = idx+1;
%             subplot(10,10,idx)
%             plot(test.control.x.raw, test.control.y.raw, 'b')
%             plot(test.stim.x.raw, test.stim.y.raw, 'y')
%             xlabel('Speed')
%             ylabel('Cumulative distribution')
%             plot(a,sorted.control.y, 'color', Color('grey'), 'LineWidth', 2)
%             plot(a,sorted.stim.y, 'color', Color('darkorange'), 'LineWidth', 2)
%             vline(a(mdiff))
%             yyaxis right
%             plot(a, test.diff, 'color', Color('brown'), 'LineWidth', 1)
%             ylabel('Difference')
%         end

        %save the data
         rdistr(n) = test.DIFF;
    end 
    
    % FIND THE ORIGINAL DIFFERENCE:
    test = []; sorted = []; a = [];
    a = fly(iflycontrol).bootstrap.speed.(type{itype});
    test.control.raw = a;
    test.control.x.raw = sort(reshape(a,[numel(a),1]));
    test.control.y.raw = (1:numel(a))/numel(a);

    a = fly(ifly).bootstrap.speed.(type{itype});
    test.stim.raw = a;
    test.stim.x.raw = sort(reshape(a,[numel(a),1]));
    test.stim.y.raw = (1:numel(a))/numel(a);

    %find the min value and max value in the two distrutions 
    minval = min([test.stim.x.raw(1), test.control.x.raw(1)]);
    maxval = max([test.stim.x.raw(end), test.control.x.raw(end)]);

    % add ones and zeros to the end of the vectors to bring to same x value
    test.stim.x.raw = [minval; test.stim.x.raw; maxval];
    test.stim.y.raw = [0; test.stim.y.raw'; 1];
    test.control.x.raw = [minval; test.control.x.raw; maxval];
    test.control.y.raw = [0; test.control.y.raw'; 1];

    % merge to find common x values:
    a = union(test.control.x.raw, test.stim.x.raw);
    sorted.control.y = zeros(length(a),1);
    for ii = 1:length(a)
        x = find(test.control.x.raw>= a(ii)); %all x values larger than the one one in question
        sorted.control.y(ii,1) = test.control.y.raw(x(1));
        xx = find(test.stim.x.raw>= a(ii)); %all x values larger than the one one in question
        sorted.stim.y(ii,1) = test.stim.y.raw(xx(1));
    end

    test.diff = abs(sorted.control.y-sorted.stim.y);
    [test.DIFF, mdiff] = max(test.diff);

%     %update figure
%     fig = getfig('',1); hold all
%         plot(test.control.x.raw, test.control.y.raw, 'b')
%         plot(test.stim.x.raw, test.stim.y.raw, 'y')
%         xlabel('Speed')
%         ylabel('Cumulative distribution')
%         plot(a,sorted.control.y, 'color', Color('grey'), 'LineWidth', 2)
%         plot(a,sorted.stim.y, 'color', Color('darkorange'), 'LineWidth', 2)
%         vline(a(mdiff))
%         yyaxis right
%         plot(a, test.diff, 'color', Color('brown'), 'LineWidth', 1)
%         ylabel('Difference')

    OG.diff = test.DIFF;
    p = sum(abs(rdistr)>=abs(OG.diff))/length(rdistr);
    pval(icount) = p;
    fly(ifly).bootstrap.rdistr = rdistr;
    fly(ifly).bootstrap.OGdiff = test;

    figure(hfig)
    subplot(3,4, icount)
    histogram(rdistr)
    vline(OG.diff, 'r-')
    title({[fly_cross_list{ifly} ' ' type{itype}]; ['p = ' num2str(p)]})
    end
end

toc






 %activation
 p_err = p_val([1,2,5,6,9,10]);
 P_err = sort(p_err);
 fprintf('\n Activation:')
for idx = 1:6
q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
r = (P_err((idx)) > (idx/length(P_err))*.05);
fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
end

%silencing
p_err = p_val([3 4 7 8 11 12]);
P_err = sort(p_err);
fprintf('\n\n Silencing:')
for idx = 1:6
q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
r = (P_err((idx)) > (idx/length(P_err))*.05);
fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
end


save_figure(hfig, ['C:\matlabroot\Tony Paper Work\Speed comparison KS distributition one tailed']);

'C:\matlabroot\Tony Paper Work\'







%%




space_adj = -0.3;

color_list = {'gold', 'orangered', 'skyblue', 'midnightblue', ...
              'pink', 'mediumvioletred','lightgreen', 'darkgreen'};
         
type = {'short', 'long'};          
          
fig = getfig;
subplot(2,2,1) %speed activation
hold all
xlim([0, 9])
ylim([-100, 180])
idx = 0;
for ifly = [7, 1:2:6]
    for tt = 1:2 %short and long light
       idx = idx+1;
        switch tt
            case 1 %short
                x = idx-space_adj;
            case 2
                x = idx+space_adj;
        end
        numtrials = length(fly(ifly).LA_speed.(type{tt}));
        X = (x-0.1):(0.2/numtrials):(x+0.1)-(0.2/numtrials);
            scatter(X,fly(ifly).LA_speed.(type{tt}), 25, Color(color_list{idx}))
        scatter(x, fly(ifly).LA_speed.avg.(type{tt}), 30,'k', 'filled')
        errorbar(x, fly(ifly).LA_speed.avg.(type{tt}), fly(ifly).LA_speed.sem.(type{tt}), 'Color', 'k')
    end
end  
hline(0, 'k:');set(gca,'TickDir','out');
title('Activation Speed changes in running flies')


subplot(2,2,3) %speed silencing
hold all
xlim([0, 9])
ylim([-100, 180])
idx = 0;
for ifly = [8, 2:2:6]
    for tt = 1:2 %short and long light
       idx = idx+1;
        switch tt
            case 1 %short
                x = idx-space_adj;
            case 2
                x = idx+space_adj;
        end
        numtrials = length(fly(ifly).LA_speed.(type{tt}));
        X = (x-0.1):(0.2/numtrials):(x+0.1)-(0.2/numtrials);
            scatter(X,fly(ifly).LA_speed.(type{tt}), 25, Color(color_list{idx}))
        scatter(x, fly(ifly).LA_speed.avg.(type{tt}), 30, 'k', 'filled')
        errorbar(x, fly(ifly).LA_speed.avg.(type{tt}), fly(ifly).LA_speed.sem.(type{tt}), 'Color', 'k')
    end
end  
hline(0, 'k:');set(gca,'TickDir','out');
title('Silencing Speed changes in running flies')


subplot(2,2,2) %speed activation
hold all
xlim([0, 9])
% ylim([-100, 180])
idx = 0;
for ifly = [7, 1:2:6]
    for tt = 1:2 %short and long light
       idx = idx+1;
        switch tt
            case 1 %short
                x = idx-space_adj;
            case 2
                x = idx+space_adj;
        end
        numtrials = length(fly(ifly).LA_rotvel.(type{tt}));
        X = (x-0.1):(0.2/numtrials):(x+0.1)-(0.2/numtrials);
            scatter(X,fly(ifly).LA_rotvel.(type{tt}), 25, Color(color_list{idx}))
        scatter(x, fly(ifly).LA_rotvel.avg.(type{tt}), 30,'k', 'filled')
        errorbar(x, fly(ifly).LA_rotvel.avg.(type{tt}), fly(ifly).LA_rotvel.sem.(type{tt}), 'Color', 'k')
    end
end  
hline(0, 'k:');set(gca,'TickDir','out');
title('Activation Rotvel changes in running flies')

subplot(2,2,4) %rotvel silencing
hold all
xlim([0, 9])
% ylim([-100, 180])
idx = 0;
for ifly = [8, 2:2:6]
    for tt = 1:2 %short and long light
       idx = idx+1;
        switch tt
            case 1 %short
                x = idx-space_adj;
            case 2
                x = idx+space_adj;
        end
        numtrials = length(fly(ifly).LA_rotvel.(type{tt}));
        X = (x-0.1):(0.2/numtrials):(x+0.1)-(0.2/numtrials);
            scatter(X,fly(ifly).LA_rotvel.(type{tt}), 25, Color(color_list{idx}))
        scatter(x, fly(ifly).LA_rotvel.avg.(type{tt}), 30, 'k', 'filled')
        errorbar(x, fly(ifly).LA_rotvel.avg.(type{tt}), fly(ifly).LA_rotvel.sem.(type{tt}), 'Color', 'k')
    end
end  
hline(0, 'k:');set(gca,'TickDir','out');
title('Silencing RotVel changes in running flies')

% 
% save_figure(fig, ['C:\matlabroot\Tony Paper Work\light effect fig ' num2str(duration)])


%% Change in movement statistics between fly line comparisons
%isolate numbers for bootstrapping
N = 10E3;
% datatype = 'rotvel';
datatype = 'speed'; 
icount = 0;
hfig = getfig;
type = {'short', 'long'};
for ifly = 1:6
  for itype = 1:2 %short|long
      icount = icount+1;
      err = [];
       if ifly <=3 
           iflycontrol = 7;
       else 
           iflycontrol = 8;
       end
      c_data = fly(iflycontrol).(['LA_' datatype]).(type{itype});
      s_data = fly(ifly).(['LA_' datatype]).(type{itype});
      % remove nans:
      c_data(isnan(c_data)) = [];
      s_data(isnan(s_data)) = [];
      % find data lengths
      ncontrol = length(c_data);
      nstim = length(s_data);
      err.all = [c_data, s_data];
      idx = 0;
      fig = getfig('', 1);
      for n = 1:N
%           test = [];
          loc = randperm(ncontrol+nstim);
          rdata = err.all(loc);
          test(n).control = mean(rdata(1:ncontrol));
          test(n).stim = mean(rdata(ncontrol+1:end));
          test(n).diff = (test(n).stim-test(n).control);
          rdistrb(n) = test(n).diff;
         % graphical test
          if mod(n,100) == 0
                idx = idx+1;
                subplot(10,10,idx)
                set(gca,'YTickLabel',[],'XTickLabel',[]); 
                hold all 
                scatter(1:1/ncontrol:2-1/ncontrol, rdata(1:ncontrol))
                scatter(3:1/nstim:4-1/nstim, rdata(ncontrol+1:end))
                xlim([0.5, 4.5])
          end
      end
      OG.diff = fly(ifly).(['LA_' datatype]).avg.(type{itype})-fly(iflycontrol).(['LA_' datatype]).avg.(type{itype});
      p = sum(abs(rdistrb)>=abs(OG.diff))/length(rdistrb);
      
    figure(hfig)
    subplot(3,4,icount)
    histogram(rdistrb)
    vline(OG.diff, 'r-')
    title({[fly_cross_list{ifly} ' ' type{itype}]; ['p = ' num2str(p)]})
    if strcmpi(datatype, 'rotvel') == true
        fly(ifly).(['rot_diff_' type{itype}]).rdistrb = rdistrb;
        fly(ifly).(['rot_diff_' type{itype}]).rawdata = test;
        fly(ifly).(['rot_diff_' type{itype}]).p_val = p;
        fprintf(['\n ' fly_cross_list{ifly} ' ' type{itype} ': p=' num2str(fly(ifly).(['rot_diff_' type{itype}]).p_val)])
    else
        fly(ifly).(['diff_' type{itype}]).rdistrb = rdistrb;
        fly(ifly).(['diff_' type{itype}]).rawdata = test;
        fly(ifly).(['diff_' type{itype}]).p_val = p;
        fprintf(['\n ' fly_cross_list{ifly} ' ' type{itype} ': p=' num2str(fly(ifly).(['diff_' type{itype}]).p_val)])
    end
     p_val(icount) = p;
  end
end



 
 %activation
 p_err = p_val([1,2,5,6,9,10]);
 P_err = sort(p_err);
 fprintf('\n Activation:')
for idx = 1:6
q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
r = (P_err((idx)) > (idx/length(P_err))*.05);
fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
end

%silencing
p_err = p_val([3 4 7 8 11 12]);
P_err = sort(p_err);
fprintf('\n\n Silencing:')
for idx = 1:6
q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
r = (P_err((idx)) > (idx/length(P_err))*.05);
fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
end


save_figure(hfig, [figures_dir, datatype ' comparison bootstrap distributition ' num2str(duration)]);










































%% Load Fly Structure (already built with "MIDSTIM_Step_4_groupflies.m")
clc; close all; clear all
 
min_speed = 0.3;

[filename, directory] = uigetfile('*.mat', 'Select your Fly Structure');
fullname = fullfile(directory, filename); 
load(fullname);
try fly = FLY; catch
end
if isfield(fly(1), 'parameters')
    for kk = 1:length(fly)
        fly(kk).param = fly(kk).parameters;
    end
    fly = rmfield(fly,'parameters');
end
kk = 1;
parameters = fly(1).param;
fprintf(['\n Loaded ' filename '\n'])
setGlobalx(parameters)

% load labels & condition parameters
[num, Fictrac, Matlab, labels] = NUM(fly(1).param);  
num.fly = length(fly);
Type = {'Stim', 'Control'};

% Create a Figures Folder for this structure:
figures_dir = [directory 'MN line figures\' filename(1:end-4) '\'];
if ~isfolder(figures_dir)
    mkdir(figures_dir)
end 

cond_figures_dir = [figures_dir 'Conditions\'];
if ~isfolder(cond_figures_dir)
    mkdir(cond_figures_dir)
end 
clc
% ELIMINATE OUTLIERS
fly = eliminate_speed_outliers(fly, 5);

% LOAD BEHAVIOR CLASSIFICATION DATA:
switch questdlg('Load behavior classification data?')
    case 'Yes'
        load([directory, structure_name, ' behavior class'])
end

% TRAJECTORY PATHS AND HEATMAPS
ROI = 1:60;
switch questdlg('Load former trajectory info?')
    case 'Yes' 
        load([figures_dir, 'Trajectory Data'])
    case {'No', 'Cancel'}
        return
end


%% Time-Course SPEED|ROT VEL Analysis grouped by initial state:(finished)
param = {'speed', 'rotvelocity'};

for condition = 1:2 %speed|rotvelocity
   
    fig = getfig('',1);
    coloropts = {'black', 'mediumblue', 'darkcyan', 'coral'}; %stationary, walking, grooming, other
    idx = 0;
    for icond = 1:num.conds %control for now
        idx = idx+1;
        for ii = 1:4
            InputData(ii).x = [];
            InputData(ii).y = [];
            InputData(ii).Ind = 1;
            InputData(ii).Color = Color(coloropts{ii});
        end

        for kk = 1:num.fly
          for rep = 1:num.reps
            % filter:
            filter = group(kk).STATE(icond,rep);
            % data:
            x = (-2:1/30:2)';
            y = [fly(kk).Control.(param{condition})(icond).data(2:end, rep);...
                 fly(kk).Stim.(param{condition})(icond).data(:, rep)];
            InputData(filter).x(:,InputData(filter).Ind) = x;
            InputData(filter).y(:,InputData(filter).Ind) = y;
            % set the next position
            InputData(filter).Ind = InputData(filter).Ind+1; 
          end
        end

        subplot(4,7,idx)
        fig = TimeCourseFigure(InputData, fig);
        getaxes(fig, 10);
        vline(0, 'k')
        title(getcond(icond))
    end
    subplot(4,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})

    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Behavior sorted time course Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Behavior sorted time course Rot Velocity']);
    end

end


%% Combined CW and CCW state-dependent time course graphs:(finished)

param = {'speed', 'rotvelocity'};
for condition = 1:2 %speed|rotvelocity

    fig = getfig('',1);
    coloropts = {'black', 'mediumblue', 'darkcyan', 'coral'}; %stationary, walking, grooming, other
    idx = 0;
    for icond = [1:7, 15:21] %control for now
        idx = idx+1;
        for ii = 1:4
            InputData(ii).x = [];
            InputData(ii).y = [];
            InputData(ii).Ind = 1;
            InputData(ii).Color = Color(coloropts{ii});
        end

        for tt = 1:2
            adj = 1;
            switch tt % CW and CCW
                case 1
                    cond = icond;
                    if condition == 2
                        adj = -1;
                    end
                case 2
                    cond = icond+7;
            end
            for kk = 1:num.fly
              for rep = 1:num.reps
                % filter:
                filter = group(kk).STATE(cond,rep);
                % data:
                x = (-2:1/30:2)';
                y = [fly(kk).Control.(param{condition})(cond).data(2:end, rep);...
                     fly(kk).Stim.(param{condition})(cond).data(:, rep)];
                InputData(filter).x(:,InputData(filter).Ind) = x;
                InputData(filter).y(:,InputData(filter).Ind) = y.*adj;
                % set the next position
                InputData(filter).Ind = InputData(filter).Ind+1; 

              end
            end
        end

        subplot(2,7,idx)
        fig = TimeCourseFigure(InputData, fig);
        getaxes(fig, 10);
        vline(0, 'k')
        title(getcond(icond))
    end
    subplot(2,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})

    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Combined behavior sorted time course Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Combined behavior sorted time course Rot Velocity']);
    end


end

%% Speed-dependent time-course of running flies (finished)

clear InputData
maxframes = num.fly*num.reps;
nframes = 0.2*num.fps;
BehaveDist.AllSpeeds = [];
% get speeds during the 200ms pre-stim:
for icond = 1:num.conds
    BehaveDist.Speeds(icond,1:maxframes) = NaN;
    temp = [];
    for kk = 1:num.fly
        filter = group(kk).walking(icond,:);
        y = mean(fly(kk).Control.speed(icond).data(end-nframes:end,:));
        fly(kk).behavior.screenspeed(icond, :) = y;
        fly(kk).behavior.group = group(kk);
        temp = [temp, y(filter)];
    end
    BehaveDist.Speeds(icond,1:length(temp)) = temp;
    BehaveDist.AllSpeeds = [BehaveDist.AllSpeeds, temp];
end

% histogram to divide by initial speeds:
nbins = 4;
% figure;
% histogram(BehaveDist.AllSpeeds,nbins); title(['Speed Distribution within ' structure_name])
% ylabel('Count'); xlabel('Mean speed during 200ms pre-stim (cm/s)')
% get the bin edges:
[~,E] = discretize(BehaveDist.AllSpeeds,nbins);
maxspeed = max(BehaveDist.AllSpeeds);

% Make the figures for both flies:

param = {'speed', 'rotvelocity'};
for condition = 1:2
    % Throw data into the order needed for generating a figure
    fig = getfig('Speed Based Analysis', 1);
    idx = 0;
    for icond = [1:7, 15:21]
        idx = idx+1;
        for tt = 1:2 % CW|CCW
            adj = 1;
            switch tt % CW and CCW
                case 1
                    cond = icond;
                    if condition == 2
                        adj = -1;
                    end
                case 2
                    cond = icond+7;
            end
            coloropts = {'Thistle', 'Orchid', 'BlueViolet', 'Indigo'};
%             coloropts = Color('Plum','Indigo', nbins);
            for ii = 1:nbins
                InputData(ii).Ind = 1;
                InputData(ii).x = (-2:1/num.fps:2)';
                InputData(ii).y = [];
                InputData(ii).Color = Color(coloropts{ii});
%                 InputData(ii).Color = coloropts(ii,:);
            end
            for kk = 1:num.fly
                for irep = 1:num.reps
                   meanspeed = fly(kk).behavior.screenspeed(cond,irep);
                    % create a sort filter for the speeds:
                    ii = discretize(meanspeed,E);
                    %filter data for walking only:
                    if group(kk).STATE(cond,irep) == 2
                        InputData(ii).y(:,InputData(ii).Ind) = ...
                            ([fly(kk).Control.(param{condition})(cond).data(2:61,irep);...
                             fly(kk).Stim.(param{condition})(cond).data(:,irep)]).*adj;
                        InputData(ii).Ind = InputData(ii).Ind+1;
                    end
                end
            end
        end
        % Generate a figure
        subplot(2,7,idx)
        fig = TimeCourseFigure(InputData, fig,1);
        getaxes(fig, 10);
        switch condition
            case 1 %speed
                ylim([0, maxspeed+0.25])
            case 2 %rotvel
        end
        vline(0, 'g')
        title(getcond(icond))
        vline((fly(1).param.conds_matrix(icond).opto), 'g')
    end
    subplot(2,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})
    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Speed dependent change in walking Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Speed dependent change in walking Rot Velocity']);
    end
end
    


%% Loaded Vs Unloaded Time Course -- speed and rotational velocity (finished)
% Variables: 
clear InputData
param = {'speed', 'rotvelocity'};
coloropts = [Color('black'); Color('grey'); Color('saddlebrown'); Color('goldenrod')];

for condition = 1:2 % speed|rotvelocity
    % Throw data into the order needed for generating a figure
    fig = getfig('', 1);
    
    idx = 0;
    for icond = [1:7, 15:21] %combine CW&CCW
        idx = idx+1; % condition count (i.e. 1-7)      
        % Set the Input Data Structure:
        for ii = 1:4 % all types of loading
            InputData(ii).Ind = 1;
            InputData(ii).x = (-2:1/num.fps:2)';
            InputData(ii).y = [];
            InputData(ii).Color = coloropts(ii,:);
        end
        for tt = 1:2 % CW|CCW
            adj = 1;
            switch tt % CW and CCW
                case 1
                    cond = icond;
                    if condition == 2
                        adj = -1;
                    end
                case 2
                    cond = icond+7;
            end
            for kk = 1:num.fly %for all flies within the condition:
                for irep = 1:num.reps
                   istate = group(kk).STATE(icond,irep);
                   iphase = group(kk).PHASE(icond,irep);
                    % send data to appropriate grouping based on state and phase
                   switch iphase %loaded or unloaded
                       case 1 %loaded
                           if istate == 1 || istate == 3 %stationary or grooming
                               ii = 1;
                           elseif istate == 2 %walking
                               ii = 3;
                           end
                       case 2 %unloaded
                           if istate == 1 || istate == 3 %stationary or grooming
                               ii = 2;
                           elseif istate == 2 %walking
                               ii = 4;
                           end
                   end
                    % load the data into the InputData structure's
                    % appropriate field:
                    if group(kk).STATE(cond,irep) < 4
                        InputData(ii).y(:,InputData(ii).Ind) = ...
                            ([fly(kk).Control.(param{condition})(cond).data(2:61,irep);...
                             fly(kk).Stim.(param{condition})(cond).data(:,irep)]).*adj;
                        InputData(ii).Ind = InputData(ii).Ind+1;
                    end
                end
            end
        end
        % Generate a figure
        subplot(2,7,idx)
        fig = TimeCourseFigure(InputData, fig, 1);
        getaxes(fig, 10);
%         switch condition
%             case 1 %speed
%                 ylim([0, maxspeed+0.25])
%             case 2 %rotvel
%         end
        vline(0, 'g')
        title(getcond(icond))
        vline((fly(1).param.conds_matrix(icond).opto), 'g')
    end
    subplot(2,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})
    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Loaded vs Unloaded time course Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Loaded vs Unloaded time course Rot Velocity']);
    end
end
    


for cond = 1:27
    
    
end






