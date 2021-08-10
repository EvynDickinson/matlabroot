
figures_dir = ['C:\matlabroot\Tony Paper Work\'];
    
%% Initiation of Movement Statistics: -- silencing gtACR1 data
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
save_figure(fig, [figures_dir, ' stationary starts STATISTICS silencing']);

%% Initiation of Movement Statistics: ACTIVATION - csChrimson data
C_ifly = 7; %control fly genotype
% isolate numbers for bootstrapping for gtACR1 data
fig = getfig('',1);
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

save_figure(fig, [figures_dir, ' stationary starts STATISTICS activation']);

%% Change in movement statistics between fly line comparisons
%isolate numbers for bootstrapping
datatype = 'Rotvel';
% datatype = 'Speed'; 
tic
duration = 0.2; % -- post stim length duration
% time ranges for data:
shorttime = [60:(60+round(.09*num.fps)+(duration*num.fps))];
longtime = [60:(60+round(.72*num.fps)+(duration*num.fps))];
% timerange(1).data = shorttime;
% timerange(1).length = length(shorttime);
% timerange(2).data = longtime;
% timerange(2).length = length(longtime);
% Isolate the time course data vectors for long activation
hfig = getfig; icount = 0;
for ifly = 1:6
    for stimtype = 2:3
        icount = icount+1;
    if ifly <=3 
        iflycontrol = 7;
    else 
        iflycontrol = 8;
    end

    % Control Data for BDP fly
%     stimtype = 2; %3 = long, 2= short
    switch stimtype
        case 2
            lightlength = '90ms';
            timerange = shorttime;
        case 3
            lightlength = '720ms';
            timerange = longtime;
    end

    % Create struct for raw data
    if strcmpi(datatype, 'Rotvel') == true
        err.BDP.control = fly(iflycontrol).Rotvel(1).y(timerange,:);
        err.BDP.stim = fly(iflycontrol).Rotvel(stimtype).y(timerange,:);
        err.MN.control = fly(ifly).Rotvel(1).y(timerange,:);
        err.MN.stim = fly(ifly).Rotvel(stimtype).y(timerange,:);
    else
        err.BDP.control = fly(iflycontrol).Speed(1).y(timerange,:);
        err.BDP.stim = fly(iflycontrol).Speed(stimtype).y(timerange,:);
        err.MN.control = fly(ifly).Speed(1).y(timerange,:);
        err.MN.stim = fly(ifly).Speed(stimtype).y(timerange,:);
    end
    %remove nans:
    LINE = {'BDP', 'MN'};
    TYPE = {'control', 'stim'};
    for iline = 1:2 %BDP|MN
        for itype = 1:2 %Control|stim
            a = err.(LINE{iline}).(TYPE{itype});
            loc = ~isnan(a(1,:));
            err.(LINE{iline}).(TYPE{itype}) = a(:,loc);
        end
        %number of trials per condition
        err.(LINE{iline}).Cnum = size(err.(LINE{iline}).control,2);
        err.(LINE{iline}).Snum = size(err.(LINE{iline}).stim,2);
    end
    % Get the original data here:
%     offset1 = nanmean(err.MN.control(1,:));
%     offset2 = nanmean(err.MN.stim(1,:));
    OG.MN_diff = (nanmean(nanmean(err.MN.control)))-...
                 (nanmean(nanmean(err.MN.stim)));
%     offset1 = nanmean(err.BDP.control(1,:));
%     offset2 = nanmean(err.BDP.stim(1,:));
    OG.BDP_diff = (nanmean(nanmean(err.BDP.control)))-...
                  (nanmean(nanmean(err.BDP.stim)));
    OG.diff = OG.BDP_diff-OG.MN_diff;

    % combined data--all possible bits to choose from:
    % control error = BDP err + MN err
    err.control_ALL = [err.BDP.control, err.MN.control];
    err.stim_ALL = [ err.BDP.stim, err.MN.stim,];
    % select a random sampling for the control & the stim

    color_idx = {Color('black'); Color('orange'); Color('grey'); Color('green')};


    N = 10E3;

%     fig = getfig(fly_cross_list{ifly}, 1);
    idx = 0;
    for n = 1:N       
        
        % controls
        loc = randperm(err.BDP.Cnum+err.MN.Cnum);
        
        a = err.control_ALL(:,loc(1:err.BDP.Cnum));
        test(n).BDP.avg_control = mean(mean(a,1)) ;%- nanmean(a(1,:));
        
        b = err.control_ALL(:,loc(err.BDP.Cnum+1:end));
        test(n).MN.avg_control = mean(mean(b,1));% - nanmean(b(1,:));

        % stimulus
        loc = randperm(err.BDP.Snum+err.MN.Snum);

        c = err.stim_ALL(:,loc(1:err.BDP.Snum));
        test(n).BDP.avg_stim = mean(mean(c,1));% - nanmean(c(1,:));%-offset;
        
        d = err.stim_ALL(:,loc(err.BDP.Snum+1:end));
        test(n).MN.avg_stim = mean(mean(d,1)) ;% - nanmean(d(1,:));%-offset;

        % graphical test
%         if mod(n,100) == 0
%             idx = idx+1;
%             subplot(10,10,idx)
%             set(gca,'YTickLabel',[],'XTickLabel',[]); 
%             hold all 
%             plot(mean(a,2), 'color', color_idx{1,:})
%             plot(mean(b,2), 'color', color_idx{3,:})
%             plot(mean(c,2), 'color', color_idx{2,:})
%             plot(mean(d,2), 'color', color_idx{4,:})
%         end
    end



    for n = 1:N
        % calculate the difference between 'in-genotype' control vs stim
        test(n).MN_diff = test(n).MN.avg_control-test(n).MN.avg_stim ;
        test(n).BDP_diff = test(n).BDP.avg_control - test(n).BDP.avg_stim;
        % difference between the BDP and MN:
        test(n).diff = test(n).BDP_diff-test(n).MN_diff;
        rdistrb(n) = test(n).diff;
    end

    p = sum(abs(rdistrb)>=abs(OG.diff))/length(rdistrb);
%     figure(fig); subplot(10,10,5)
%     title({[fly_cross_list{ifly} ' ' lightlength]; ['p = ' num2str(p)]})

    figure(hfig)
    subplot(3,4,icount)
    histogram(rdistrb)
    vline(OG.diff, 'r-')
    title({[fly_cross_list{ifly} ' ' lightlength]; ['p = ' num2str(p)]})
    if strcmpi(datatype, 'Rotvel') == true
        fly(ifly).(['rot_diff_' lightlength]).rdistrb = rdistrb;
        fly(ifly).(['rot_diff_' lightlength]).rawdata = test;
        fly(ifly).(['rot_diff_' lightlength]).p_val = p;
        fprintf(['\n ' fly_cross_list{ifly} ' ' lightlength ': p=' num2str(fly(ifly).(['rot_diff_' lightlength]).p_val)])
    else
        fly(ifly).(['diff_' lightlength]).rdistrb = rdistrb;
        fly(ifly).(['diff_' lightlength]).rawdata = test;
        fly(ifly).(['diff_' lightlength]).p_val = p;
        fprintf(['\n ' fly_cross_list{ifly} ' ' lightlength ': p=' num2str(fly(ifly).(['diff_' lightlength]).p_val)])
    end
    
    p_val(icount) = p;
    end

end
fprintf('\n')
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


save_figure(hfig, [figures_dir, datatype ' comparison bootstrap distributition ' num2str(duration)]);



%% In-house bootstrapping of the effects:
err = [];
color_idx = {Color('black'); Color('orange'); Color('grey'); Color('green')};
%isolate numbers for bootstrapping
datatype = 'Rotvel';
% datatype = 'Speed';
tic
duration = 0.5; % -- post stim length duration
% time ranges for data:
shorttime = 60:(60+round(.09*num.fps)+(duration*num.fps));
longtime = 60:(60+round(.72*num.fps)+(duration*num.fps));
% timerange(1).data = shorttime;
% timerange(1).length = length(shorttime);
% timerange(2).data = longtime;
% timerange(2).length = length(longtime);
% Isolate the time course data vectors for long activation
hfig = getfig; icount = 0;

for ifly = 1:8

% Create struct for raw data
if strcmpi(datatype, 'Rotvel') == true
    err.short.control = fly(ifly).Rotvel(1).y(shorttime,:);
    err.short.stim = fly(ifly).Rotvel(2).y(shorttime,:);
    err.long.control = fly(ifly).Rotvel(1).y(longtime,:);
    err.long.stim = fly(ifly).Rotvel(3).y(longtime,:);
else
   err.short.control = fly(ifly).Speed(1).y(shorttime,:);
    err.short.stim = fly(ifly).Speed(2).y(shorttime,:);
    err.long.control = fly(ifly).Speed(1).y(longtime,:);
    err.long.stim = fly(ifly).Speed(3).y(longtime,:);
end

TYPE = {'control', 'stim'};
lightlength = {'short', 'long'};
for ilength = 1:2 %short and long
    icount = icount+1;
    for itype = 1:2 %Control|stim
        a = err.(lightlength{ilength}).(TYPE{itype});
        loc = ~isnan(a(1,:));
        err.(lightlength{ilength}).(TYPE{itype}) = a(:,loc);
    end
    %number of trials per condition
    Cnum = size(err.(lightlength{ilength}).control,2);
    Snum = size(err.(lightlength{ilength}).stim,2);

    % original data
    OG.(lightlength{ilength}).diff = (nanmean(nanmean(err.(lightlength{ilength}).control)))-...
                 (nanmean(nanmean(err.(lightlength{ilength}).stim)));

    err.(lightlength{ilength}).ALL =...
        [err.(lightlength{ilength}).control, err.(lightlength{ilength}).stim];

    % select a random sampling for the control & the stim
    N = 10E3;
    fig = getfig([fly_cross_list{ifly} ' ' lightlength{ilength}], 1);
    idx = 0;
    for n = 1:N 
        loc = randperm(Cnum+Snum);
        a = err.(lightlength{ilength}).ALL(:,loc(1:Cnum));
        test(n).(lightlength{ilength}).avg_control = mean(mean(a,1)) ; 
        b = err.(lightlength{ilength}).ALL(:,loc(Cnum+1:end));
        test(n).(lightlength{ilength}).avg_stim = mean(mean(b,1));

        % graphical test
        if mod(n,100) == 0
            idx = idx+1;
            subplot(10,10,idx)
            set(gca,'YTickLabel',[],'XTickLabel',[]); 
            hold all 
            plot(mean(a,2), 'color', color_idx{1,:})
            plot(mean(b,2), 'color', color_idx{3,:})
        end
    end
    for n = 1:N
        % calculate the difference between 'in-genotype' control vs stim
        test(n).(lightlength{ilength}).diff = ...
                    test(n).(lightlength{ilength}).avg_control-...
                    test(n).(lightlength{ilength}).avg_stim ;
        rdistrb(n) = test(n).(lightlength{ilength}).diff;
    end

    p = sum(abs(rdistrb)>=abs(OG.(lightlength{ilength}).diff))/N;
    figure(fig); subplot(10,10,5)
    title({[fly_cross_list{ifly} ' ' lightlength{ilength}]; ['p = ' num2str(p)]})

    figure(hfig)  
    subplot(4,4,icount)
    histogram(rdistrb)
    vline(OG.(lightlength{ilength}).diff, 'r-')
    title({[fly_cross_list{ifly} ' ' lightlength{ilength}]; ['p = ' num2str(p)]})


    if strcmpi(datatype, 'Rotvel') == true
        fly(ifly).timecourse.(['rot_diff_' lightlength{ilength}]).rdistrb = rdistrb;
        fly(ifly).timecourse.(['rot_diff_' lightlength{ilength}]).rawdata = test;
        fly(ifly).timecourse.(['rot_diff_' lightlength{ilength}]).p_val = p;
    else
        fly(ifly).timecourse.(['diff_' lightlength{ilength}]).rdistrb = rdistrb;
        fly(ifly).timecourse.(['diff_' lightlength{ilength}]).rawdata = test;
        fly(ifly).timecourse.(['diff_' lightlength{ilength}]).p_val = p;
    end
    fprintf(['\n ' fly_cross_list{ifly} ' ' lightlength{ilength} ': p=' num2str(fly(ifly).timecourse.(['diff_' lightlength{ilength}]).p_val)])
    p_val(icount) = p;

end
end


%activation
 p_err = p_val;
 P_err = sort(p_err);
 fprintf('\n All:')
for idx = 1:icount
q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
r = (P_err((idx)) > (idx/length(P_err))*.05);
fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
end




% %activation
%  p_err = p_val([1,2,5,6,9,10]);
%  P_err = sort(p_err);
%  fprintf('\n Activation:')
% for idx = 1:6
% q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
% r = (P_err((idx)) > (idx/length(P_err))*.05);
% fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
% end
% 
% %silencing
% p_err = p_val([3 4 7 8 11 12]);
% P_err = sort(p_err);
% fprintf('\n\n Silencing:')
% for idx = 1:6
% q = (P_err((idx)) > 0.05/(length(P_err) +1 - idx));
% r = (P_err((idx)) > (idx/length(P_err))*.05);
% fprintf(['\ninsig: ' num2str(q) ' vs ' num2str(r)])
% end


save_figure(hfig, [figures_dir, datatype ' inhouse comparison bootstrap distributition ' num2str(duration)]);


%%

% O.G. calc of difference between stim and control
for ifly = 1:length(fly)
figures_dir = ['C:\matlabroot\Tony Paper Work\' fly(ifly).structure_name '\'];
% structure_name = fly(ifly).structure_name;

for pp = 1:2 %speed|rotvel
    switch pp
        case 1
            param = 'speed';
            InputData = fly(ifly).Speed;
        case 2
            param = 'rotvelocity';
            InputData = fly(ifly).Rotvel;
    end
    %find the offsets:
    offset.range = 57:59; %200ms preceeding the stimulus
    if ifly == 7
        offset.control = nanmean(nanmean(InputData(1).y(:,offset.range)));
        LA(1).offset = offset.control-nanmean(nanmean(InputData(2).y(:,offset.range))); %short offset
        LA(2).offset = offset.control-nanmean(nanmean(InputData(3).y(:,offset.range))); %long offset
    else
        
        offset.control = nanmean(nanmean(InputData(1).y(offset.range,:)));
        LA(1).offset = offset.control-nanmean(nanmean(InputData(2).y(offset.range,:))); %short offset
        LA(2).offset = offset.control-nanmean(nanmean(InputData(3).y(offset.range,:))); %long offset
    end
    % time ranges for data:
    shorttime = [60:(60+round(.09*num.fps)+(duration*num.fps))];
    longtime = [60:(60+round(.72*num.fps)+(duration*num.fps))];
    timerange(1).data = shorttime;
    timerange(1).length = length(shorttime);
    timerange(2).data = longtime;
    timerange(2).length = length(longtime);
    
    %add data ranges to new structure
    if ifly == 7
        LA(1).control = InputData(1).y(:,timerange(1).data)'; %short control
        LA(2).control = InputData(1).y(:,timerange(2).data)'; %long control
        LA(1).stim = InputData(2).y(:,timerange(1).data)'; %short light stim
        LA(2).stim = InputData(3).y(:,timerange(2).data)'; %long light stim
    else
        LA(1).control = InputData(1).y(timerange(1).data,:); %short control
        LA(2).control = InputData(1).y(timerange(2).data,:); %long control
        LA(1).stim = InputData(2).y(timerange(1).data,:); %short light stim
        LA(2).stim = InputData(3).y(timerange(2).data,:); %long light stim
    end

    for tt = 1:2 %short light|long light
        LA(tt).controlavg = nanmean(sum(LA(tt).control)/timerange(tt).length);
        LA(tt).controlerr = nanstd(sum(LA(tt).control))/fly(ifly).num.fly;
        LA(tt).stimavg = nanmean(sum(LA(tt).stim)/timerange(tt).length)+LA(tt).offset;
        LA(tt).stimerr = nanstd(sum(LA(tt).stim))/fly(ifly).num.fly;
        LA(tt).diff = LA(tt).stimavg/LA(tt).controlavg;
        LA(tt).differr = sqrt((LA(tt).controlerr/LA(tt).controlavg)^2+(LA(tt).stimerr/LA(tt).stimavg)^2);
    end

    fig = getfig;
    hold all
    for tt = 1:2
       scatter(tt, LA(tt).diff, 100, 'filled', 'k')
       errorbar(tt, LA(tt).diff , LA(tt).differr/num.fly, 'Color', 'k')
    end
    xlim([0, 3])
    ylim([0, 2])
    xlabel('Short light, long light')
    title({fly_cross_list{ifly}; ['Change in ' param]; ''})
    ylabel(['fraction change in ' param])
    set(gca,'TickDir','out');
    figure_name = [figures_dir, fly_cross_list{ifly}, ' Change in ' param];
     switch pp
        case 1
            fly(ifly).LA_speed = LA;
        case 2
            fly(ifly).LA_rotvel = LA;
     end
    save(figure_name, 'LA')
    save_figure(fig, figure_name);
%    close (fig)
   
end
end



% activation - long


for ifly = 1:3

Ca = fly(7).ldata.light.totalstationary;
Cm = fly(7).ldata.light.totalmoved;
Sa = fly(ifly).ldata.light.totalstationary;
Sm = fly(ifly).ldata.light.totalmoved;

% Stim and unstim are vectors of ones and zeros, in this case of 40% vs 25%
stim = zeros(Sa,1); stim(randperm(Sa,Sm)) = 1;
unstim = zeros(Ca,1); unstim(randperm(Ca,Cm)) = 1;
Dpctg = sum(stim)/numel(stim) - sum(unstim)/numel(unstim);

N = 10E4;
alltrials = [stim;unstim];
dpctg_reps = 1:N;
for n = 1:N
   inds = randperm(numel(alltrials));
   stim_draw = alltrials(inds(1:numel(stim)));
   unstim_draw = alltrials(inds(numel(stim)+1:end));
   dpctg_reps(n) = sum(stim_draw)/numel(stim_draw) - sum(unstim_draw)/numel(unstim_draw);
end


Plot distributions
figure, f = gcf;
a = histogram(dpctg_reps);
a.Parent.NextPlot = 'add';
plot([1 1]*Dpctg,a.Parent.YLim,'r');
p = sum(dpctg_reps>=Dpctg)/length(dpctg_reps);
fprintf(['\n ' fly_cross_list{ifly} ' p='])
disp(p)


end


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
        
%         save(figure_name, 'LA')
%         save_figure(fig, figure_name);
        close(fig)
end
fprintf(['\n Finished: ' fly_cross_list{ifly}])
end



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
    











%% Change in movement statistics 
% shows the results of bootstrapping with baseline adjustments...doesn't
% look right though...
% %isolate numbers for bootstrapping
% tic
% duration = 0.5; % -- post stim length duration
% % time ranges for data:
% shorttime = [60:(60+round(.09*num.fps)+(duration*num.fps))];
% longtime = [60:(60+round(.72*num.fps)+(duration*num.fps))];
% % timerange(1).data = shorttime;
% % timerange(1).length = length(shorttime);
% % timerange(2).data = longtime;
% % timerange(2).length = length(longtime);
% % Isolate the time course data vectors for long activation
% 
% for ifly = 1:6
%     for stimtype = 2:3
%     if ifly <=3 
%         iflycontrol = 7;
%     else 
%         iflycontrol = 8;
%     end
% 
%     % Control Data for BDP fly
% %     stimtype = 2; %3 = long, 2= short
%     switch stimtype
%         case 2
%             lightlength = '90ms';
%             timerange = shorttime;
%         case 3
%             lightlength = '720ms';
%             timerange = longtime;
%     end
% 
%     % Create struct for raw data
%     err.BDP.control = fly(iflycontrol).Speed(1).y(timerange,:);
%     err.BDP.stim = fly(iflycontrol).Speed(stimtype).y(timerange,:);
%     err.MN.control = fly(ifly).Speed(1).y(timerange,:);
%     err.MN.stim = fly(ifly).Speed(stimtype).y(timerange,:);
%     %remove nans:
%     LINE = {'BDP', 'MN'};
%     TYPE = {'control', 'stim'};
%     for iline = 1:2 %BDP|MN
%         for itype = 1:2 %Control|stim
%             a = err.(LINE{iline}).(TYPE{itype});
%             loc = ~isnan(a(1,:));
%             err.(LINE{iline}).(TYPE{itype}) = a(:,loc);
%         end
%         %number of trials per condition
%         err.(LINE{iline}).Cnum = size(err.(LINE{iline}).control,2);
%         err.(LINE{iline}).Snum = size(err.(LINE{iline}).stim,2);
%     end
%     % Get the original data here:
%     offset1 = nanmean(err.MN.control(1,:));
%     offset2 = nanmean(err.MN.stim(1,:));
%     OG.MN_diff = (nanmean(nanmean(err.MN.control))-offset1)-...
%                  (nanmean(nanmean(err.MN.stim)-offset2));
%     offset1 = nanmean(err.BDP.control(1,:));
%     offset2 = nanmean(err.BDP.stim(1,:));
%     OG.BDP_diff = (nanmean(nanmean(err.BDP.control))-offset1)-...
%                   (nanmean(nanmean(err.BDP.stim))-offset2);
%     OG.diff = OG.BDP_diff-OG.MN_diff;
% 
%     % combined data--all possible bits to choose from:
%     % control error = BDP err + MN err
%     err.control_ALL = [err.BDP.control, err.MN.control];
%     err.stim_ALL = [ err.BDP.stim, err.MN.stim,];
%     % select a random sampling for the control & the stim
% 
%     color_idx = {Color('black'); Color('orange'); Color('grey'); Color('green')};
% 
% 
%     N = 10E3;
% 
%     fig = getfig(fly_cross_list{ifly}, 1);
%     idx = 0;
%     for n = 1:N
%          offset1 = nanmean(err.MN.control(1,:));
%          offset2 = nanmean(err.MN.stim(1,:));
%         
%         
%         % controls
%         loc = randperm(err.BDP.Cnum+err.MN.Cnum);
%         
%         a = err.control_ALL(:,loc(1:err.BDP.Cnum));
%         test(n).BDP.avg_control = mean(mean(a,1)) - nanmean(a(1,:));
%         
%         b = err.control_ALL(:,loc(err.BDP.Cnum+1:end));
%         test(n).MN.avg_control = mean(mean(b,1)) - nanmean(b(1,:));
% 
%         % stimulus
%         loc = randperm(err.BDP.Snum+err.MN.Snum);
% 
%         c = err.stim_ALL(:,loc(1:err.BDP.Snum));
%         test(n).BDP.avg_stim = mean(mean(c,1)) - nanmean(c(1,:));%-offset;
%         
%         d = err.stim_ALL(:,loc(err.BDP.Snum+1:end));
%         test(n).MN.avg_stim = mean(mean(d,1)) - nanmean(d(1,:));%-offset;
%         
%         
%         
%         % graphical test
%         if mod(n,100) == 0
%             idx = idx+1;
%             subplot(10,10,idx)
%             set(gca,'YTickLabel',[],'XTickLabel',[]); 
%             hold all 
%             plot(mean(a,2) - nanmean(a(1,:)), 'color', color_idx{1,:})
%             plot(mean(b,2) - nanmean(b(1,:)), 'color', color_idx{3,:})
%             plot(mean(c,2) - nanmean(c(1,:)), 'color', color_idx{2,:})
%             plot(mean(d,2) - nanmean(d(1,:)), 'color', color_idx{4,:})
%         end
%     end
% 
% 
% 
%     for n = 1:N
%         % calculate the difference between 'in-genotype' control vs stim
%         test(n).MN_diff = test(n).MN.avg_control-test(n).MN.avg_stim ;
%         test(n).BDP_diff = test(n).BDP.avg_control - test(n).BDP.avg_stim;
%         % difference between the BDP and MN:
%         test(n).diff = test(n).BDP_diff-test(n).MN_diff;
%         rdistrb(n) = test(n).diff;
%     end
% 
%     p = sum(abs(rdistrb)>=abs(OG.diff))/length(rdistrb);
%     figure(fig); subplot(10,10,5)
%     title({[fly_cross_list{ifly} ' ' lightlength]; ['p = ' num2str(p)]})
% 
%     figure;
%     histogram(rdistrb)
%     vline(OG.diff, 'r-')
%     title({[fly_cross_list{ifly} ' ' lightlength]; ['p = ' num2str(p)]})
%     fly(ifly).(['diff_' lightlength]).rdistrb = rdistrb;
%     fly(ifly).(['diff_' lightlength]).rawdata = test;
% 
% end
% end
% 
% toc