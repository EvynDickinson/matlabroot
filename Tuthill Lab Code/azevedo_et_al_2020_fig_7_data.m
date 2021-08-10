

%% Figure 7 %% 

%% Panel B - headless fly joint angles
% BDP
load('C:\matlabroot\Tony Paper Work\data\BDP-gal4xUAS-csChrimson-headless-offball angle data.mat')
for itrial = 1:20
    Fig7PanelB.BDP_Chrimson(:,itrial) = data(itrial).intrpdata;
end
% FAST MN 81A07
load('C:\matlabroot\Tony Paper Work\data\81A07xgtACR1-headless-offball angle data.mat')
for itrial = 1:15
    Fig7PanelB.MN_81A07_Chrimson(:,itrial) = data(itrial).intrpdata;
end
% INTERMEDIATE MN 22A08
load('C:\matlabroot\Tony Paper Work\data\22A08xgtACR1-headless-offball angle data.mat')
for itrial = 1:14
    Fig7PanelB.MN_22A08_Chrimson(:,itrial) = data(itrial).intrpdata;
end
% SLOW MN 35C09
load('C:\matlabroot\Tony Paper Work\data\35C09xgtACR1-headless-offball angle data.mat')
for itrial = 1:17
    Fig7PanelB.MN_35C09_Chrimson(:,itrial) = data(itrial).intrpdata;
end

% save the data into a clean new structure:

save('G:\My Drive\MN packaged data\Fig_7_Panel_B', 'Fig7PanelB');



%% Figure 7 Panel E -- intact walking flies speed trajectories
% save('G:\My Drive\MN packaged data\Fig_7_Panel_E', 'Fig7PanelE');


% load the data individually
load('BDP-gal4xUAS-gtACR1.mat')
cross = 'BDP_gtACR1';

load('BDP-gal4xUAS-csChrimson.mat')
cross = 'BDP_csChrimson';

load('D:\MN line figures\DATA\81A07 data\81A07-gal4xUAS-csChrimson.mat')
cross = 'MN_81A07_csChrimson';

load('D:\MN line figures\DATA\81A07 data\81A07-gal4xUAS-gtACR1.mat')
cross = 'MN_81A07_gtACR1';

load('D:\MN line figures\DATA\22A08 data\22A08-gal4xUAS-csChrimson.mat')
cross = 'MN_22A08_csChrimson';

load('D:\MN line figures\DATA\22A08 data\22A08-gal4xUAS-gtACR1.mat')
cross = 'MN_22A08_gtACR1';

load('D:\MN line figures\DATA\35C09 data\35C09-gal4xUAS-csChrimson.mat')
cross = 'MN_35C09_csChrimson';

load('D:\MN line figures\DATA\35C09 data\35C09-gal4xUAS-gtACR1.mat')
cross = 'MN_35C09_gtACR1';


% manually drop the raw data into this:
cross = 'BDP_gtACR1';
numfly = length(FLY);
% pull the speed data out for trials 1, 4, 7 in the CW and CCW direction
load(['C:\matlabroot\behavior class\' structure_name ' behavior class']);
 

% Convert behvariors into numbers and then a filter:
for kk = 1:numfly
    for cond = 1:28
        for rep = 1:3
            state = group(kk).behavior{cond, rep};
            phase = group(kk).phase{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            if ~ischar(phase)
                phase = cell2mat(phase);
            else
                PHASE = 0;
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
                case {'-'}
                    PHASE = 0;
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

type = 'walking';
fly = FLY;

for cond = 1:7 
    for ifly = 1:numfly
    % select the data and generate a filter for the initial behavior state
    filter = [group(ifly).(type)(cond,:), group(ifly).(type)(cond+7,:)];  
    cw = [fly(ifly).Control.speed(cond).data(1:end-1,:); fly(ifly).Stim.speed(cond).data];
    ccw = [fly(ifly).Control.speed(cond+7).data(1:end-1,:); fly(ifly).Stim.speed(cond+7).data];
    input = [cw, ccw];
    data(cond).ALL(ifly).data = input;
    input(:,~filter) = nan; %nix the trials that didn't fit the movement type
    % if there aren't 2 trials min, all go NaN
    if sum(filter)<1 %2 ADJUSTING THIS NUMBER FOR A TRIAL
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



% pull out the walking data per fly and then that for each trial   
% control:
    a = data(1).walking.avg;
    Fig7PanelE.(cross).control_walking_data.fly_avg_speed = a;
    output = [];
    for ifly = 1:numfly
        b = data(1).walking.raw(ifly).data;
        output = [output, b];
    end
    loc = isnan(output(1,:));
    output(:,loc)=[];
    Fig7PanelE.(cross).control_walking_data.all_trials = output;
% 90ms stim
    a = data(4).walking.avg;
    Fig7PanelE.(cross).stim90ms_walking_data.fly_avg_speed = a;
    output = [];
    for ifly = 1:numfly
        b = data(4).walking.raw(ifly).data;
        output = [output, b];
    end
    loc = isnan(output(1,:));
    output(:,loc)=[];
    Fig7PanelE.(cross).stim90ms_walking_data.all_trials = output;
% 720ms stim
    a = data(7).walking.avg;
    Fig7PanelE.(cross).stim720ms_walking_data.fly_avg_speed = a;
    output = [];
    for ifly = 1:numfly
        b = data(7).walking.raw(ifly).data;
        output = [output, b];
    end
    loc = isnan(output(1,:));
    output(:,loc)=[];
    Fig7PanelE.(cross).stim720ms_walking_data.all_trials = output;


 clear a b ccw cw cond cross data filter fly FLY group ifly input kk loc 
 clear numfly output phase PHASE rep state STATE structure_name type   
    
    
%% Panel H - stationary intact flies that start moving

% save('G:\My Drive\MN packaged data\Fig_7_Panel_H', 'Fig7PanelH');
   
 
cross = 'BDP_CsChrimson';
struct_name = 'BDP-gal4xUAS-csChrimson';

cross = 'BDP_gtACR1';
struct_name = 'BDP-gal4xUAS-gtACR1';

cross = 'MN_81A07_gtACR1';
struct_name = '81A07-gal4xUAS-gtACR1';

cross = 'MN_81A07_CsChrimson';
struct_name = '81A07-gal4xUAS-csChrimson';

cross = 'MN_22A08_gtACR1';
struct_name = '22A08-gal4xUAS-gtACR1';

cross = 'MN_81A07_CsChrimson';
struct_name = '81A07-gal4xUAS-csChrimson';

cross = 'MN_35C09_gtACR1';
struct_name = '35C09-gal4xUAS-gtACR1';

cross = 'MN_35C09_CsChrimson';
struct_name = '35C09-gal4xUAS-csChrimson';
 
load(['C:\matlabroot\Tony Paper Work\' struct_name '\' struct_name ' Fly starts from stationary.mat']);

Fig7PanelH.(cross).stim90ms.no_laser.percent_moved = data.control.percent;
Fig7PanelH.(cross).stim90ms.no_laser.total_moved = data.control.totalmoved;
Fig7PanelH.(cross).stim90ms.no_laser.total_possible = data.control.totalstationary;

Fig7PanelH.(cross).stim90ms.laser.percent_moved = data.light.percent;
Fig7PanelH.(cross).stim90ms.laser.total_moved = data.light.totalmoved;
Fig7PanelH.(cross).stim90ms.laser.total_possible = data.light.totalstationary;
 
load(['C:\matlabroot\Tony Paper Work\' struct_name '\' struct_name ' Fly starts from stationary 720ms.mat']);

Fig7PanelH.(cross).stim720ms.no_laser.percent_moved = data.control.percent;
Fig7PanelH.(cross).stim720ms.no_laser.total_moved = data.control.totalmoved;
Fig7PanelH.(cross).stim720ms.no_laser.total_possible = data.control.totalstationary;
 
Fig7PanelH.(cross).stim720ms.laser.percent_moved = data.light.percent;
Fig7PanelH.(cross).stim720ms.laser.total_moved = data.light.totalmoved;
Fig7PanelH.(cross).stim720ms.laser.total_possible = data.light.totalstationary;
 

 
 
 
 
 
 
 
 
 
 
    
    