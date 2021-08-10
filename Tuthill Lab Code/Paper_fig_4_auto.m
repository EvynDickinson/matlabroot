clear all
% load data
fly_cross_list = {'81A07-gal4xUAS-csChrimson', '81A07-gal4xUAS-gtACR1', ...
                  '22A08-gal4xUAS-csChrimson', '22A08-gal4xUAS-gtACR1', ...
                  '35C09-gal4xUAS-csChrimson', '35C09-gal4xUAS-gtACR1', ...
                  'BDP-gal4xUAS-csChrimson', 'BDP-gal4xUAS-gtACR1'};


for ifly = 1:length(fly_cross_list)
    fly_cross = fly_cross_list{ifly};
    
    figures_dir = ['C:\matlabroot\Tony Paper Work\' fly_cross '\'];
%     load([figures_dir, fly_cross ' Change in speed']) % need to find this
%     fly(ifly).LA_speed = LA;
%     
%     load([figures_dir, fly_cross ' Change in rotvelocity']) % need to find this
%     fly(ifly).LA_rotvel = LA;
    
    load([figures_dir, fly_cross ' Percent Change in speed']) % need to find this
    fly(ifly).LA_speed = LA;
    
    load([figures_dir, fly_cross ' Percent Change in rotvelocity']) % need to find this
    fly(ifly).LA_rotvel = LA;
    
    
    load([figures_dir, fly_cross ' speed time course figure'])
    fly(ifly).Speed = InputData;
    
    load([figures_dir, fly_cross ' rotvel time course figure'])
    fly(ifly).Rotvel = InputData;
    
    load([figures_dir, fly_cross ' Fly starts from stationary']) % need to find this
    fly(ifly).data = data;
    load([figures_dir, fly_cross ' Fly starts from stationary 720ms'])
    fly(ifly).ldata = data;
    load([figures_dir, fly_cross ' Average Trajectory walking selective'])
    fly(ifly).stim = stim;
    fly(ifly).ctrl = ctrl;
    
    clear LA InputData data stim ctrl
    
    load(['C:\matlabroot\MN line figures\', fly_cross, '\Trajectory Data'])
    fly(ifly).traj = traj;
    load(['C:\matlabroot\behavior class\', fly_cross,  ' behavior class'])
    % Convert behvariors into numbers and then a filter:
    for kk = 1:num.fly
    for cond = 1:num.conds
        for rep = 1:num.reps
            state = group(kk).behavior{cond, rep};
            phase = group(kk).phase{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            if ~ischar(phase)
                phase = cell2mat(phase);
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
            end
            switch phase
                case {'stance', 'Stance'}
                    PHASE = 1;
                case {'swing', 'Swing'}
                    PHASE = 2;
            end 
            group(kk).STATE(cond, rep) = STATE;
            group(kk).PHASE(cond, rep) = PHASE;
        end
    end
    group(kk).walking = (group(kk).STATE==2);
    group(kk).stationary = (group(kk).STATE==1);
    group(kk).grooming = (group(kk).STATE==3);
    group(kk).other = (group(kk).STATE==4);
    group(kk).stance = (group(kk).PHASE==1);
    group(kk).swing = (group(kk).PHASE==2);
end

    fly(ifly).group = group;
    fly(ifly).structure_name = fly_cross;
    fly(ifly).parameters = parameters;
    fly(ifly).num = num;
end
% 
% for ii =1:3
% reset = fly(7).Speed(ii).y';
% fly(7).Speed(ii).y = [];
% fly(7).Speed(ii).y = reset;
% end
% 
% for ii =1:3
% reset = fly(7).Rotvel(ii).y';
% fly(7).Rotvel(ii).y = [];
% fly(7).Rotvel(ii).y = reset;
% end

%%
duration = 0.2;


for ifly = 1:length(fly)
figures_dir = ['C:\matlabroot\Tony Paper Work\' fly(ifly).structure_name '\'];
% structure_name = fly(ifly).structure_name;

for pp = 1:2
    clear trial LA
    switch pp
        case 1
            param = 'speed';
            InputData = fly(ifly).Speed;
        case 2
            param = 'rotvelocity';
            InputData = fly(ifly).Rotvel;
    end
    shorttime = [60:(60+round(.09*num.fps)+(duration*num.fps))];
    longtime = [60:(60+round(.72*num.fps)+(duration*num.fps))];
    timerange(1).data = shorttime;
    timerange(2).data = longtime;
    numtrials = size(InputData(1).y,2);
        % find the percent change in speed for each fly in the group -- if
        % no control data is present--disgard that fly's data.
        for kk = 1:numtrials %per fly
            %raw data
            trial(kk).short.control = InputData(1).y(timerange(1).data,kk);
            trial(kk).short.stim = InputData(2).y(timerange(1).data,kk);
            
            trial(kk).long.control = InputData(1).y(timerange(2).data,kk);
            trial(kk).long.stim = InputData(3).y(timerange(2).data,kk);
            %percent change SHORT stim
            a = mean(trial(kk).short.stim);
            b = mean(trial(kk).short.control);
            trial(kk).short.change = (a-b)/b*100;
            LA.short(kk) = trial(kk).short.change;
            %percent change LONG stim
            a = mean(trial(kk).long.stim);
            b = mean(trial(kk).long.control);
            trial(kk).long.change = (a-b)/b*100;
            LA.long(kk) = trial(kk).long.change;
        end
         %check for outlier in 6
        if pp == 1
            LA.short(abs(LA.short)>200) = nan;
            LA.long(abs(LA.long)>200) = nan;
        else
            LA.short(abs(LA.short)>1000) = nan;
            LA.long(abs(LA.long)>1000) = nan;
        end
        % FLY LINE AVERAGES:
        LA.avg.short = nanmean(LA.short);
        LA.avg.long = nanmean(LA.long);
        LA.sem.short = nanstd(LA.short);%/sqrt(numtrials);
        LA.sem.long = nanstd(LA.long);%/sqrt(numtrials);
        
        fig = figure;
            set(fig, 'color', 'w', 'pos', [-841,315,560,420]);
            hold all
            x = 1:1/numtrials:2-1/numtrials;
            scatter(x,LA.short)
            x = 3:1/numtrials:4-1/numtrials;
            scatter(x,LA.long)

            scatter(1.5, LA.avg.short, 'k', 'filled')
            scatter(3.5, LA.avg.long, 'k', 'filled')

            errorbar(1.5, LA.avg.short, LA.sem.short, 'Color', 'k')
            errorbar(3.5, LA.avg.long, LA.sem.long, 'Color', 'k')
            xlim([0, 5])
            xlabel('Short light, long light')
        title({fly_cross_list{ifly}; ['Percent Change in ' param]; ''})
        set(gca,'TickDir','out');
        ylabel(['percent change in ' param])
        figure_name = [figures_dir, fly_cross_list{ifly}, ' Percent Change in ' param];
        switch pp
            case 1
                fly(ifly).LA_speed = LA;
            case 2
                fly(ifly).LA_rotvel = LA;
        end
        
        save(figure_name, 'LA')
        save_figure(fig, figure_name);
     
end
fprintf(['\n Finished: ' fly_cross_list{ifly}])
end



%%  Group the percent change in speed data into one figure

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
save_figure(fig, ['C:\matlabroot\Tony Paper Work\light effect fig ' num2str(duration)])


    
%%  Initiation of movement Bar Graphs:
% data = short stim ; ldata = long stim;

color_list = {'black','orangered','black', 'midnightblue', 'black', 'mediumvioletred','black','darkgreen'};

fig = getfig; 
subplot(1,2,1)
ylim([0,100])
data = [];
idx = 0;
for ifly = [7,1:2:6]
   for tt = 1:2 
    idx = idx+1;
    data(idx) = fly(ifly).data.bargraph(tt);
%     data(idx) = fly(ifly).ldata.bargraph(tt);
   end
end

hold all
b = bar(1:8, data, 'BarWidth', 0.8,'FaceColor','flat');
% for ii = 1:8
%    b.CData(ii,:) = Color(color_list{ii}) ;
% end
set(gca,'TickDir','out');
title('activation')
subplot(1,2,2)
data = [];
idx = 0;
for ifly = [8,2:2:7]
   for tt = 1:2 
    idx = idx+1;
    data(idx) = fly(ifly).data.bargraph(tt);
%     data(idx) = fly(ifly).ldata.bargraph(tt);
   end
end
ylim([0,100])
hold all
b = bar(1:8, data, 'BarWidth', 0.8,'FaceColor','flat');
% for ii = 1:8
%    b.CData(ii,:) = Color(color_list{ii}) ;
%     
% end
set(gca,'TickDir','out');
title('silencing')

figure_name = 'C:\matlabroot\Tony Paper Work\Fly starts from stationary 90ms';
% figure_name = 'C:\matlabroot\Tony Paper Work\Fly starts from stationary 720ms';
save_figure(fig, figure_name);




%% Initiation of Movement Statistics:
C_ifly = 8; %control fly genotype
% isolate numbers for bootstrapping for gtACR1 data
fig = getfig;
type = {'short', 'long'};
idx = 0;
for tt = 1:2 %short and long laser
    for ifly = 2:2:6
        idx = idx+1;
        switch tt 
            case 1
                % Short data
                Ca = fly(C_ifly).data.light.totalstationary;
                Cm = fly(C_ifly).data.light.totalmoved;
                Sa = fly(ifly).data.light.totalstationary;
                Sm = fly(ifly).data.light.totalmoved;
            case 2
                % Long data:
                Ca = fly(C_ifly).ldata.light.totalstationary;
                Cm = fly(C_ifly).ldata.light.totalmoved;
                Sa = fly(ifly).ldata.light.totalstationary;
                Sm = fly(ifly).ldata.light.totalmoved; 
        end
        
        % Stim and unstim are vectors of ones and zeros, in this case of 40% vs 25%
        stim = zeros(Sa,1); stim(randperm(Sa,Sm)) = 1;
        unstim = zeros(Ca,1); unstim(randperm(Ca,Cm)) = 1;
        Dpctg = sum(stim)/numel(stim) - sum(unstim)/numel(unstim);

        noiserate = 0.1;
        N = 10E4;
        alltrials = [stim;unstim];
        noisenum = round(numel(alltrials)*noiserate);
        dpctg_reps = 1:N;
        for n = 1:N
           inds = randperm(numel(alltrials));
           %            rand(length(Rinds),1)>=0.5;
           ralltrials = alltrials;
           ralltrials(1:noisenum) = rand(noisenum,1)>=0.5;

           stim_draw = ralltrials(inds(1:numel(stim)));
           unstim_draw = ralltrials(inds(numel(stim)+1:end));
           
           dpctg_reps(n) = sum(stim_draw)/numel(stim_draw) - sum(unstim_draw)/numel(unstim_draw);
        end
        
        % Plot distributions
        subplot(2,3,idx)
        f = gcf;
        a = histogram(dpctg_reps);
        a.Parent.NextPlot = 'add';
        plot([1 1]*Dpctg,a.Parent.YLim,'r');

        p = sum(abs(dpctg_reps)>=abs(Dpctg))/length(dpctg_reps);
        title({[fly_cross_list{ifly} ' ' type{tt}];[' p=' num2str(p)]})
        fprintf(['\n ' fly_cross_list{ifly} ' p=' num2str(p)])

        fly(ifly).stationarystarts_err.(type{tt}) = p;
        p_err(idx) = p;

    end
end

% sort(p_err)

% for idx = 1:6
% 
% q = (p_err(sort(idx)) > 0.05/(length(p_err) +1 - idx));
% fprintf(['\ninsig: ' num2str(q)])
% end
sort(p_err)
P_err = sort(p_err);
for idx = 1:6

q = (P_err((idx)) > 0.05/(length(p_err) +1 - idx));
r = (P_err((idx)) > (idx/length(p_err))*.05);


fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
end


%% Initiation of Movement Statistics: ACTIVATION
C_ifly = 7; %control fly genotype
% isolate numbers for bootstrapping for gtACR1 data
fig = getfig;
type = {'short', 'long'};
idx = 0;
for tt = 1:2 %short and long laser
    for ifly = 1:2:6
        idx = idx+1;
        switch tt 
            case 1
                % Short data
                Ca = fly(C_ifly).data.light.totalstationary;
                Cm = fly(C_ifly).data.light.totalmoved;
                Sa = fly(ifly).data.light.totalstationary;
                Sm = fly(ifly).data.light.totalmoved;
            case 2
                % Long data:
                Ca = fly(C_ifly).ldata.light.totalstationary;
                Cm = fly(C_ifly).ldata.light.totalmoved;
                Sa = fly(ifly).ldata.light.totalstationary;
                Sm = fly(ifly).ldata.light.totalmoved; 
        end
        
        % Stim and unstim are vectors of ones and zeros, in this case of 40% vs 25%
        stim = zeros(Sa,1); stim(randperm(Sa,Sm)) = 1;
        unstim = zeros(Ca,1); unstim(randperm(Ca,Cm)) = 1;
        Dpctg = sum(stim)/numel(stim) - sum(unstim)/numel(unstim);
noiserate = 0.1;
        N = 10E4;
        alltrials = [stim;unstim];
        noisenum = round(numel(alltrials)*noiserate);
        dpctg_reps = 1:N;
        for n = 1:N
           inds = randperm(numel(alltrials));
           %            rand(length(Rinds),1)>=0.5;
           ralltrials = alltrials;
           ralltrials(1:noisenum) = rand(noisenum,1)>=0.5;

           stim_draw = ralltrials(inds(1:numel(stim)));
           unstim_draw = ralltrials(inds(numel(stim)+1:end));
           
           dpctg_reps(n) = sum(stim_draw)/numel(stim_draw) - sum(unstim_draw)/numel(unstim_draw);
        end

        % Plot distributions
        subplot(2,3,idx)
        f = gcf;
        a = histogram(dpctg_reps);
        a.Parent.NextPlot = 'add';
        plot([1 1]*Dpctg,a.Parent.YLim,'r');

        p = sum(abs(dpctg_reps)>=abs(Dpctg))/length(dpctg_reps);
        title({[fly_cross_list{ifly} ' ' type{tt}];[' p=' num2str(p)]})
        fprintf(['\n ' fly_cross_list{ifly} ' p=' num2str(p)])

        fly(ifly).stationarystarts_err.(type{tt}) = p;
        p_err(idx) = p;

    end
end

sort(p_err)
P_err = sort(p_err);

for idx = 1:6

q = (P_err((idx)) > 0.05/(length(p_err) +1 - idx));
r = (P_err((idx)) > (idx/length(p_err))*.05);


fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
end













    
  %% Generate figures:

% duration = 0.2;
% 
% for ifly = 7:8%1:length(fly)
% figures_dir = ['C:\matlabroot\Tony Paper Work\' fly(ifly).structure_name '\'];
% % structure_name = fly(ifly).structure_name;
% 
% for pp = 1:2
%     switch pp
%         case 1
%             param = 'speed';
%             InputData = fly(ifly).Speed;
%         case 2
%             param = 'rotvelocity';
%             InputData = fly(ifly).Rotvel;
%     end
%     %find the offsets:
%     offset.range = 57:59;
% %     if ifly == 7
% %         offset.control = nanmean(nanmean(InputData(1).y(:,offset.range)));
% %         LA(1).offset = offset.control-nanmean(nanmean(InputData(2).y(:,offset.range))); %short offset
% %         LA(2).offset = offset.control-nanmean(nanmean(InputData(3).y(:,offset.range))); %long offset
% %     else
%         offset.control = nanmean(nanmean(InputData(1).y(offset.range,:)));
%         LA(1).offset = offset.control-nanmean(nanmean(InputData(2).y(offset.range,:))); %short offset
%         LA(2).offset = offset.control-nanmean(nanmean(InputData(3).y(offset.range,:))); %long offset
% %     end
%     % time ranges for data:
%     shorttime = [60:(60+round(.09*num.fps)+(duration*num.fps))];
%     longtime = [60:(60+round(.72*num.fps)+(duration*num.fps))];
%     timerange(1).data = shorttime;
%     timerange(1).length = length(shorttime);
%     timerange(2).data = longtime;
%     timerange(2).length = length(longtime);
%     
%     %add data ranges to new structure
% %     if ifly == 7
% %         LA(1).control = InputData(1).y(:,timerange(1).data)'; %short control
% %         LA(2).control = InputData(1).y(:,timerange(2).data)'; %long control
% %         LA(1).stim = InputData(2).y(:,timerange(1).data)'; %short light stim
% %         LA(2).stim = InputData(3).y(:,timerange(2).data)'; %long light stim
% %     else
%         LA(1).control = InputData(1).y(timerange(1).data,:); %short control
%         LA(2).control = InputData(1).y(timerange(2).data,:); %long control
%         LA(1).stim = InputData(2).y(timerange(1).data,:); %short light stim
%         LA(2).stim = InputData(3).y(timerange(2).data,:); %long light stim
% %     end
% 
%     for tt = 1:2 %short light|long light
%         LA(tt).controlavg = nanmean(sum(LA(tt).control)/timerange(tt).length);
%         LA(tt).controlerr = nanstd(sum(LA(tt).control))/fly(ifly).num.fly;
%         LA(tt).stimavg = nanmean(sum(LA(tt).stim)/timerange(tt).length)+LA(tt).offset;
%         LA(tt).stimerr = nanstd(sum(LA(tt).stim))/fly(ifly).num.fly;
% %         LA(tt).diff = LA(tt).stimavg/LA(tt).controlavg;
%         LA(tt).diff = (LA(tt).stimavg)-LA(tt).controlavg/LA(tt).controlavg;
%         LA(tt).differr = sqrt((LA(tt).controlerr/LA(tt).controlavg)^2+(LA(tt).stimerr/LA(tt).stimavg)^2);
%     end
% 
%     fig = getfig;
%     hold all
%     for tt = 1:2
%        scatter(tt, LA(tt).diff, 100, 'filled', 'k')
%        errorbar(tt, LA(tt).diff , LA(tt).differr/num.fly, 'Color', 'k')
%     end
%     xlim([0, 3])
% %     ylim([-2, 2])
%     xlabel('Short light, long light')
%     title({fly_cross_list{ifly}; ['Change in ' param]; ''})
%     ylabel(['fraction change in ' param])
%     set(gca,'TickDir','out');
%     figure_name = [figures_dir, fly_cross_list{ifly}, ' Percent Change in ' param];
%      switch pp
%         case 1
%             fly(ifly).LA_speed = LA;
%         case 2
%             fly(ifly).LA_rotvel = LA;
%      end
%     save(figure_name, 'LA')
%     save_figure(fig, figure_name);
% %     close (fig)
%    
% end
% end
% 
% 
% space_adj = -0.3;
% % color_list = {'midnightblue', 'skyblue', 'midnightblue', 'skyblue',...
% %               'mediumvioletred', 'pink', 'mediumvioletred', 'pink',...
% %               'darkgreen', 'lightgreen', 'darkgreen', 'lightgreen'...
% %               'orangered', 'gold'};
% % color_list = {'orangered', 'gold', 'midnightblue', 'skyblue', ...
% %               'mediumvioletred', 'pink','darkgreen', 'lightgreen'};
% color_list = {'gold', 'orangered', 'skyblue', 'midnightblue', ...
%               'pink', 'mediumvioletred','lightgreen', 'darkgreen'};
%           
% 
% fig = getfig;
% subplot(2,2,1) %speed activation
% hold all
% xlim([0, 9])
% % ylim([-0.4,0.6])
% idx = 0;
% vline(0.5, 'r-')
% for ifly = [7, 1:2:6]
%     for tt = 1:2 %short and long light
%       idx = idx+1;  
%         switch tt
%             case 1
%                 x = idx-space_adj;
%             case 2
%                 x = idx+space_adj;
%         end
%         
%        scatter(x, fly(ifly).LA_speed(tt).diff, 100, Color(color_list{idx}), 'filled', 'MarkerEdgeColor', 'none')
%        errorbar(x,  fly(ifly).LA_speed(tt).diff ,  fly(ifly).LA_speed(tt).differr/ fly(ifly).num.fly,...
%                 'Color', Color(color_list{idx}))
%     end
%     
%     vline(idx+0.5, 'r-')
% end  
% hline(0, 'k:');set(gca,'TickDir','out');
% title('Activation Speed changes in running flies')
% 
% 
% subplot(2,2,3) %speed silencing
% hold all
% xlim([0, 9])
% % ylim([-.4, 0.6])
% idx = 0;
% vline(0.5, 'r-')
% for ifly = [8, 2:2:6]
%     for tt = 1:2 %short and long light
%       idx = idx+1;  
%         switch tt
%             case 1
%                 x = idx-space_adj;
%             case 2
%                 x = idx+space_adj;
%         end
%         
%        scatter(x, fly(ifly).LA_speed(tt).diff, 100, Color(color_list{idx}), 'filled', 'MarkerEdgeColor', 'none')
%        errorbar(x,  fly(ifly).LA_speed(tt).diff ,  fly(ifly).LA_speed(tt).differr/ fly(ifly).num.fly,...
%                 'Color', Color(color_list{idx}))
%     end
%     
%      vline(idx+0.5, 'r-')
% end  
% hline(0, 'k:');set(gca,'TickDir','out');
% title('Silencing Speed changes in running flies')
% 
% 
% subplot(2,2,2) %rot velocity
% hold all
% xlim([0, 9])
% % ylim([0, 1.4])
% vline(0.5, 'r-')
% idx = 0;
% for ifly = [7, 1:2:6]
%     for tt = 1:2 %short and long light
%          idx = idx+1;  
%        switch tt
%             case 1
%                 x = idx-space_adj;
%             case 2
%                 x = idx+space_adj;
%         end
%        scatter(x, fly(ifly).LA_rotvel(tt).diff, 100, Color(color_list{idx}), 'filled', 'MarkerEdgeColor', 'none')
%        errorbar(x,  fly(ifly).LA_rotvel(tt).diff, fly(ifly).LA_rotvel(tt).differr/fly(ifly).num.fly,...
%                'Color', Color(color_list{idx}))
%     end
%     vline(idx+0.5, 'r-')
% end
% hline(0, 'k:');set(gca,'TickDir','out');
% title('activation rotvelocity changes in turning flies')
% 
% subplot(2,2,4) %rot velocity
% hold all
% xlim([0, 9])
% % ylim([0, 1.4])
% vline(0.5, 'r-')
% idx = 0;
% for ifly = [8, 2:2:6]
%     for tt = 1:2 %short and long light
%          idx = idx+1;  
%        switch tt
%             case 1
%                 x = idx-space_adj;
%             case 2
%                 x = idx+space_adj;
%         end
%        scatter(x, fly(ifly).LA_rotvel(tt).diff, 100, Color(color_list{idx}), 'filled', 'MarkerEdgeColor', 'none')
%        errorbar(x,  fly(ifly).LA_rotvel(tt).diff, fly(ifly).LA_rotvel(tt).differr/fly(ifly).num.fly,...
%                'Color', Color(color_list{idx}))
%     end
%     vline(idx+0.5, 'r-')
% end
% hline(0, 'k:');set(gca,'TickDir','out');
% title('silencing rotvelocity changes in turning flies')
% 
% % 
% save_figure(fig, ['C:\matlabroot\Tony Paper Work\light effect fig ' num2str(duration)])
  
    