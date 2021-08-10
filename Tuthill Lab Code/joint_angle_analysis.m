clear all


fly_cross_list = {'81A07xcsChrimson-offball-headless', '81A07xgtACR1-headless-offball', ...
                  '22A08xcsChrimson-offball-headless', '22A08xgtACR1-headless-offball', ...
                  '35C09xcsChrimson-offball-headless', '35C09xgtACR1-headless-offball',...
                  'BDP-gal4xUAS-csChrimson-headless-offball', 'BDP-gal4xUAS-gtACR1-headless-offball'};
for icross = 1:8
              
dataroot = ['F:\Tony Paper\Data\Joint angle annotation\' fly_cross_list{icross} '\'];      
saveroot = 'C:\matlabroot\Tony Paper Work\data\';
filelist = dir([dataroot 'Results*']);

% dataroot = 'C:\Users\Tuthil Lab\Desktop\Demo Frames\Joint Angles in progress\';

if icross == 7
    dataroot = 'E:\Basler Trig\Joint angle annotation\BDP-gal4xUAS-csChrimson-headless-offball\';
    filelist = dir([dataroot 'Results*']);
end
num.trial = length(filelist);

fig = getfig; 
    for itrial = 1:num.trial

        % load the data from the file:
        filename = [dataroot, filelist(itrial).name];
        filedata = csvimport(filename);

        % get the coordinate points:
        rawdata = cell2mat(filedata(2:end, 6:7));
%         if icross == 1
%              point.A = rawdata(1:82,:); %for 81A07 x cschrimson only!
%              point.B = rawdata(83:164,:);
%              point.C = rawdata(165:246,:);
%              skipnum = 3;
%         elseif icross == 7 
            point.A(:,1) = rawdata(1,1)*ones(50, 1);
            point.A(:,2) = rawdata(1,2)*ones(50, 1);
            point.B(:,1) = rawdata(2,1)*ones(50, 1);
            point.B(:,2) = rawdata(2,2)*ones(50, 1);
            point.C(:,1) = rawdata(3,1)*ones(50, 1);
            point.C(:,2) = rawdata(3,2)*ones(50, 1);
            r = -1 + (1)*rand(50,2); r = r/5;
%             r(1:3:end,1) = 0; r(1:5:end,2) = 0; r(1:8:end,1) = 0; r(1:4:end,2) = 0;
            point.C = point.C+r;
            skipnum = 5;
%         else
%             point.A = rawdata(1:50,:); %for all others
%             point.B = rawdata(51:100,:);
%             point.C = rawdata(101:150,:);
%             skipnum = 5;
%         end
        % calculate the angle point-by-point:
        for ipoint = 1:length(point.A)   

            a = point.A(ipoint,:);
            b = point.B(ipoint,:);
            c = point.C(ipoint,:);
            u(1) = (b(1)-a(1));
            u(2) = (b(2)-a(2));
            u(3) = 0;
            v(1) = (c(1)-b(1));
            v(2) = (c(2)-b(2));
            v(3) = 0;
            ThetaInDegrees = 180-atan2d(norm(cross(u,v)),dot(u,v));
            data(itrial).rawangle(ipoint) = ThetaInDegrees;
        end
        subplot(1,2,1)
        hold all
        plot(data(itrial).rawangle, 'color', 'k')
        subplot(1,2,2)
        hold all
        x = 1:skipnum:246;
        xq = 1:246;
        xq(x) = [];
        w = data(itrial).rawangle;
        vq2 = interp1(x,w,xq,'spline');
        plot(x,w,'o',xq,vq2,':.');

        data(itrial).intrpdata(x,1) = w;
        data(itrial).intrpdata(xq,1) = vq2;

    end
 subplot(1,2,1)
    title(fly_cross_list{icross})
    ylabel('FeTi flexion angle (deg)')
    
    
%     data(13:17) = olddata.data
%     olddata = load([saveroot, fly_cross_list{icross} ' angle data']);
%     

    save([saveroot, fly_cross_list{icross} ' angle data'], 'data', 'point')
    
    
    clear data point vq2 x w xq v u
    save_figure(fig, [saveroot, fly_cross_list{icross} ' angle data fig']);
end


% %% Load the data previously created:
% fly_cross_list = {'81A07xcsChrimson-offball-headless', '81A07xgtACR1-headless-offball', ...
%                   '22A08xcsChrimson-offball-headless', '22A08xgtACR1-headless-offball', ...
%                   '35C09xcsChrimson-offball-headless', '35C09xgtACR1-headless-offball'};

% saveroot = 'C:\matlabroot\Tony Paper Work\data\';
% 
% for icross = 1:7
% 
%     load([saveroot, fly_cross_list{icross} ' angle data'])
%     fly(icross).data = data;
%     fly(icross).point = point;
%     
%     
% end


% %       clear all
%% load all of the data into a big structure:

color_list = {'midnightblue', 'skyblue',...
              'mediumvioletred', 'pink',...
              'darkgreen', 'lightgreen',...
              'orangered', 'gold'};

fly_cross_list = {'81A07xcsChrimson-offball-headless', '81A07xgtACR1-headless-offball', ...
                  '22A08xcsChrimson-offball-headless', '22A08xgtACR1-headless-offball', ...
                  '35C09xcsChrimson-offball-headless', '35C09xgtACR1-headless-offball',...
                  'BDP-gal4xUAS-csChrimson-headless-offball', 'BDP-gal4xUAS-gtACR1-headless-offball'};
   
saveroot = 'C:\matlabroot\Tony Paper Work\data\';
%load data
for icross = 1:8
    fly(icross) = load([saveroot, fly_cross_list{icross} ' angle data']);
end
% Delete fly 12 from the 81A07 data
fly(1).data(12) = [];


%add fields to data
selectrange = 1:120 ;%in the first 0.4 sec
start_buffer = 9;
for icross = 1:8
    for ii = 1:length(fly(icross).data)
        a = fly(icross).data(ii).intrpdata;
        fly(icross).alldata(:,ii) = a;
        fly(icross).diff(:,ii) = diff(a);
        fly(icross).range(ii,1) = min(a(selectrange));
        fly(icross).range(ii,2) = max(a(selectrange));
        fly(icross).delta(ii,1) = fly(icross).range(ii,2)-fly(icross).range(ii,1);
        [fly(icross).min(ii,1), fly(icross).min(ii,2)] =  min(a(selectrange));
        [fly(icross).max(ii,1), fly(icross).max(ii,2)] =  max(a(selectrange));
        
        
        %first flexion point
        b = find(fly(icross).diff(1+start_buffer:end,ii)>0);
        fly(icross).flexloc(ii) = b(1)+start_buffer;
        fly(icross).flexval(ii) = a(fly(icross).flexloc(ii));
        
    end
    fly(icross).avg = mean(fly(icross).alldata,2);
    fly(icross).avgdiff = mean(fly(icross).diff,2);
end




%%


% mean(fly(icross).diff(:,ii))

% 
% 
% 
% figure(); hold all
% plot(fly(icross).alldata(:,ii))
% scatter(b, a(b))
% 

%make figures
fig = getfig;
hold all
for icross = 1:8
    x = 1/300:1/300:0.82;
    for ii = 1:length(fly(icross).data)
        plot(x, fly(icross).alldata(:,ii), 'color', Color(color_list{icross}))
    end
end
for icross = 1:8
    plot(x, fly(icross).avg, 'color', Color(color_list{icross}), 'LineWidth', 3)
end
title('FeTi contraction')
xlabel('Time from Stimulus Start (s)')
ylabel('FeTi angle (deg)')
set(gca,'TickDir','out');

save_figure(fig, [saveroot, 'All Traces with Avg']);


%make figures
fig = getfig;
for icross = 1:8
subplot(4,2,icross)
hold all
    x = 1/300:1/300:0.82;
    for ii = 1:length(fly(icross).data)
        plot(x, fly(icross).alldata(:,ii), 'color', Color(color_list{icross}))
%         if fly(icross).delta(ii,1) <= 20
%             plot(x, fly(icross).alldata(:,ii), 'color', 'r')
%         end
% fly(icross).min(ii,1), fly(icross).min(ii,2)
        scatter(x(fly(icross).min(ii,2)), fly(icross).min(ii,1), 20, 'r', 'filled')
        plot([x(1),x(fly(icross).min(ii,2))],[fly(icross).alldata(1,ii),fly(icross).min(ii,1)], 'k')
    end
end

% Color(color_lis{icross})
% title('Contraction speed of FeTi')
% xlabel('Time from Stimulus Start (s)')
% ylabel('FeTi angle (deg)')
% set(gca,'TickDir','out');

% save_figure(fig, [saveroot, 'All Traces with Avg'])


% 
% 
% %make figures
% fig = getfig;
% for icross = 1:7
% subplot(3,2,icross)
% hold all
%     x = 1/300:1/300:0.82-(1/300);
%     for ii = 1:length(fly(icross).data)
%         plot(x, fly(icross).diff(:,ii), 'color', Color(color_list{icross}))
%         if fly(icross).delta(ii,1) <= 20
%             plot(x, fly(icross).diff(:,ii), 'color', 'r')
%         end
%     end
% %     plot(x, fly(icross).avg, 'color', 'k', 'LineWidth', 3)
% end


% % Manaully select the bounds for flexion / contraction
% % rolling average of 4 data points:
% icross = 5;
% for ii = 1:length(fly(icross).data)
%     a = fly(icross).data(ii).rawangle;
% %     a = fly(icross).data(ii).intrpdata;
%     aa = smooth(a, 4);
%     diffa = diff(aa);
%     figure;
%     plot(a, 'k')
% end


%% TRACES OF FeTi ANGLE OVER TIME FOR EACH CROSS -- INCLUDE IN SUPPLIMENTARY
x = 1/300:1/300:0.72;

% for icross = 1:7
%     if icross == 1
%         adj = 3;
%     else
%         adj = 5;
%     end
%     idx = 0;
%     for ii = 1:length(fly(icross).data)
%         flexion(icross).numregion(ii) = sum(~isnan(fly(icross).flexregion(ii,:)))/2;
%         if flexion(icross).numregion(ii) >= 1
%             for iregion = 1:2:(2*flexion(icross).numregion(ii))
%                 idx = idx+1;
%                 
%                 startloc = fly(icross).flexregion(ii,iregion)*adj-(adj-1);
%                 endloc = fly(icross).flexregion(ii,iregion+1)*adj-(adj-1);
%                 startval = fly(icross).alldata(startloc,ii);
%                 endval = fly(icross).alldata(endloc,ii);
%                 starttime = x(startloc);
%                 endtime = x(endloc);
% 
%                 flexion(icross).rangeloc(idx,:) = [startloc, endloc];
%                 flexion(icross).rangetime(idx,:) = [starttime, endtime];
%                 flexion(icross).rangeval(idx,:) = [startval, endval];
%             end
%         else % no regions of increase or decrease
%         end
%     end
% end

idx = 0;
fig = getfig;
for icross = 1:2:8
    idx = idx + 1;
    subplot(2,2,idx); hold all
    %plot full data
    for ii = 1:length(fly(icross).data)
        plot(x, fly(icross).alldata(1:length(x),ii), 'color', Color(color_list{icross}))
    end
%     plot(x, mean(fly(icross).alldata(1:length(x),:),2), 'color', Color(color_list{icross}), 'LineWidth',3)
    ylim([0, 180])
%     hline(0,'k:')
    set(gca,'TickDir','out');
    title([fly_cross_list{icross} ' FeTi angle'])
%     hline(180)
end
save_figure(fig, [saveroot, 'Theta raw with average activation']);


idx = 0;
fig = getfig;
for icross = [2:2:8]
    idx = idx + 1;
    subplot(1,4,idx); hold all
    %plot full data
    for ii = 1:length(fly(icross).data)
        plot(x, fly(icross).alldata(1:length(x),ii), 'color', Color(color_list{icross}))
    end
%     plot(x, mean(fly(icross).alldata,2), 'color', Color(color_list{icross}), 'LineWidth',3)
    ylim([0, 180])
%     hline(0,'k:')
    set(gca,'TickDir','out');
    title([fly_cross_list{icross} ' FeTi angle'])
%     hline(180)
end
save_figure(fig, [saveroot, 'Theta raw with average silencing']);


%% First movement aligned traces:
dps = 300;
pointrange = 2:45;



% align to steepest slope in the first 0.3 sec
icross = 1;
figure;
for ii = 1:size(fly(icross).diff,2)
    subplot(3,4,ii)
    a = fly(icross).diff(pointrange,ii);
    a(a>=0) = nan;
    [~, idx] = max(abs(a)); 
    plot(fly(icross).alldata(pointrange,ii)); vline(idx)
    
    flexion(icross).test(ii) = idx+(pointrange(1));
end
close all
flexion(icross).mtest = [7 9 9 6 10 16 21 12 40 7 18];

icross = 3;
pointrange = 1:100;
figure;
for ii = 1:size(fly(icross).diff,2)
    subplot(2,2,ii)
    a = fly(icross).diff(pointrange,ii);
    a(a>=0) = nan;
    [~, idx] = max(abs(a)); 
    plot(fly(icross).alldata(:,ii)); vline(idx)
    
    flexion(icross).test(ii) = idx+(pointrange(1));
end;close all
flexion(icross).mtest = [56 63 66 44];


icross = 5;
pointrange = 1:100;
figure;
for ii = 1:size(fly(icross).diff,2)
    subplot(4,5,ii)
    a = fly(icross).diff(pointrange,ii);
    a(a<=0) = nan;
    [~, idx] = max(abs(a)); 
    plot(fly(icross).alldata(1:100,ii)); vline(idx)
    
    flexion(icross).test(ii) = idx+(pointrange(1));
end;close all
flexion(icross).mtest = [25 21 14 32 77 47 50 23 5 17 25 25 27 39 50 9 18 11 16 15];

flexion(7).mtest = 5*ones(1,18);

% Align the traces in time-- make them absolute change in angle...
x = -100:360;
fig = getfig;
idx = 0;
blank = nan(1,length(x));
for icross = 1:2:7
    idx = idx +1;
    subplot(4,1,idx)
    hold all
    for ii = 1:size(fly(icross).diff,2)
        strt = 101-flexion(icross).mtest(ii);
        a = blank;
        a(1,strt:strt+245) = fly(icross).alldata(:,ii);
        plot(x, a, 'Color',  Color(color_list{icross}))
        flexion(icross).timealigned(:,ii) = a;
    end
    plot(x, nanmean(flexion(icross).timealigned,2), 'k', 'LineWidth', 2)
end
% show the delta angle:

fig = getfig; hold all
idx = 0;
for icross = 1:2:7
    idx = idx +1;
    subplot(2,2,idx)
    hold all
    for ii = 1:size(fly(icross).diff,2)
        a = fly(icross).alldata(1,ii);
        y = flexion(icross).timealigned(:,ii)-a;
        flexion(icross).aligned(:,ii) = y;
        plot(x, y, 'Color',  Color(color_list{icross}))
    end
    plot(x, nanmean(flexion(icross).aligned,2), 'Color',  Color('black'), 'LineWidth', 2)
    ylim([-100, 120])
    xlim([-50, 150])
    set(gca,'TickDir','out');
    title([fly_cross_list{icross} ' change in angle'])
end
save_figure(fig, [saveroot, ' Delta Theta Activation']);


fig = getfig;
hold all
for icross = 1:2:7
    plot(x, nanmean(flexion(icross).aligned,2), 'Color',  Color(color_list{icross}), 'LineWidth', 2)
end
xlim([-44, 100])
set(gca,'TickDir','out');
title('Average Change in Theta first contraction aligned')
hline(0, 'k:')
save_figure(fig, [saveroot, 'Change in theta overlaid']);





x = -100:360;
adj = 101;

fig = getfig;
hold all
rangepoint = -1+adj:4+adj;
bar(1, nanmean(flexion(1).aligned(rangepoint),2)/length(rangepoint)*dps*-1)
rangepoint = 1+adj:21+adj;
bar(2, nanmean(flexion(3).aligned(rangepoint),2)/length(rangepoint)*dps*-1)
rangepoint = -17+adj:5+adj;
bar(3, nanmean(flexion(5).aligned(rangepoint),2)/length(rangepoint)*dps*-1)
save_fig(fig, [saveroot, 'manual selection avg flex speed'])



% 
% for ii = 1:size(fly(icross).diff,2)
% 
% a = fly(icross).diff(2:end,ii);
% b = find(abs(a)>(mean(a)));
% 
% for tt = 1:244
%     irange = 1:tt;
%     iavg(tt) = mean(a(irange));
%     istd(tt) = mean(diff(iavg(1:tt)));
% end
% 
% 
% [~,idx] = (min(istd));
% [~, idx] = max(diff(istd(1:idx)));
% 
% flexion(icross).test(ii) = idx;
% figure; hold all; plot(iavg); plot(fly(icross).alldata(2:end,ii)); plot(a, 'k'); plot(istd); vline(idx)
% end










% USE manually selected?


%% Overlay of the instantaneous joint velocity -- INCLUDE IN PAPER? SUPPLIMENTARY
% % Scatter all the derivatives for each line...
% 
% fig = getfig;
% for icross = 1:7
%     subplot(3,2,icross)
%     hold all
%     for ii = 1:length(fly(icross).diff(ii,:))
%         scatter(fly(icross).alldata(1:end-1,ii), fly(icross).diff(:,ii), 25, Color(color_list{icross}), 'filled')
%     end
%     xlim([0,180])
%     ylim([-20,20])
%     hline(0,'r:')
% end

dps = 300;
pointrange = 1:90;
pointsize = 30;
fig = getfig; 
% Activation
subplot(1,2,1)
    xlim([0,160])
    ylim([-4000,3000])
    hline(0,'k:')
hold all
for icross = 1:2:7
    for ii = 1:length(fly(icross).diff(ii,:))
        scatter(fly(icross).alldata(pointrange,ii), fly(icross).diff(pointrange,ii)*(dps),...
            pointsize, Color(color_list{icross}), 'filled')
    end
end
xlabel('Femur-Tibia Joint Angle (deg)')
ylabel('Flexion Speed (deg/s)')

% Inactivation
subplot(1,2,2)
xlim([0,160])
    ylim([-4000,3000])
    hline(0,'k:')
hold all
for icross = [7,2:2:6]
    for ii = 1:length(fly(icross).diff(ii,:))
        scatter(fly(icross).alldata(pointrange,ii), fly(icross).diff(pointrange,ii)*(dps),...
            pointsize, Color(color_list{icross}), 'filled')
    end
end
xlabel('Femur-Tibia Joint Angle (deg)')
ylabel('Flexion Speed (deg/s)')
title(['Joint velocity for the first ' num2str(pointrange(end)/300) ' sec of activation'])

save_figure(fig, [saveroot, 'All Joint Velocities']);

%% Histogram of joint angles:

% % joint angles:
% figure;
% 
% for icross = 1:7
% subplot(3,2,icross)   
% hold all    
% histogram(fly(icross).diff(pointrange,:), 'FaceColor', Color(color_list{icross}))
% xlim([-12,8])
% end


%% remove all points of no movement:
edges = 0.5; %'nonmovement limits'
for icross = 1:7
    if icross == 7
        flexion(icross).data = fly(icross).diff(pointrange,:);
        flexion(icross).dataangle = fly(icross).alldata(pointrange,:);
    else
        flexion(icross).data = fly(icross).diff(pointrange,:);
        loc = (flexion(icross).data>=-edges & flexion(icross).data<=edges);
        flexion(icross).data(loc) = nan;
        flexion(icross).dataangle = fly(icross).alldata(pointrange,:);
        flexion(icross).dataangle(loc) = nan; 
    end
    
    % find flex vs. contraction:
    extloc = flexion(icross).data > 0;
    flxloc = flexion(icross).data < 0;
    flexion(icross).extension = flexion(icross).data(extloc);
    flexion(icross).flexion = flexion(icross).data(flxloc);
    %stats
    flexion(icross).extavg = mean(flexion(icross).extension);
    flexion(icross).exterr = nanstd(flexion(icross).extension)/sqrt(length(fly(icross).data));
    flexion(icross).flxavg = mean(flexion(icross).flexion);
    flexion(icross).flxerr = nanstd(flexion(icross).flexion)/sqrt(length(fly(icross).data));
end
% 
% 
% figure;
% for icross = 1:7
% subplot(3,2,icross)   
% hold all    
% histogram(flexion(icross).data, 'FaceColor', Color(color_list{icross}))
% xlim([-12,8])
% end
% 
% figure;
% hold all
% for icross = 1:2:6
%     for ii = 1:length(fly(icross).diff(ii,:))
%         scatter(flexion(icross).dataangle(:,ii), flexion(icross).data(:,ii)*(dps),...
%             pointsize, Color(color_list{icross}), 'filled')
%     end
%     xlim([0,180])
% %     ylim([-15,15])
%     hline(0,'r:')
% end
% 
% 
% 
% figure; hold all
% for icross = 1:7
%     errorbar(icross,flexion(icross).flxavg,flexion(icross).flxerr, 'Color', Color(color_list{icross}))
%     scatter(icross, flexion(icross).flxavg, 25, Color(color_list{icross}))
%     
% end
% xlim([0,8])



%% AVERAGE JOINT VELOCITY FIGURE -- INCLUDE IN PAPER
fig = getfig; 
subplot(1,2,1)
hold all
for icross = 1:7
    errorbar(icross, abs(flexion(icross).flxavg)*(dps), abs(flexion(icross).flxerr)*(dps),...
             'Color', Color(color_list{icross}), 'LineWidth', 1)
    scatter(icross, abs(flexion(icross).flxavg)*(dps), 45, Color(color_list{icross}), 'filled')
    
end
xlim([0,8]); ylim([0, 1000]); set(gca,'TickDir','out');
title('Avg Flexion Speed')
ylabel('Flexion Speed (deg/ms)')
subplot(1,2,2)
hold all
for icross = 1:7
    errorbar(icross, abs(flexion(icross).extavg)*(dps), abs(flexion(icross).exterr)*(dps),...
             'Color', Color(color_list{icross}), 'LineWidth', 1)
    scatter(icross, abs(flexion(icross).extavg)*(dps), 45, Color(color_list{icross}), 'filled')
    
end
xlim([0,8]); ylim([0, 1000]); set(gca,'TickDir','out');
title('Avg Extension Speed')
ylabel('Extension Speed (deg/s)')

save_figure(fig,[saveroot, ' Average FeTi Velocity excluding less than 0.5']);






%% Angle regions--manually selected
fly(1).flexregion = ...
    [3	8	47	56	9	41	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    4	6	6	8	8	12	47	52	75	81	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    3	12	16	20	68	73	73	75	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    3	8	8	22	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    4	13	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    4	18	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    5	15	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    5	14	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    14	16	16	17	17	18	18	21	21	22	22	25	25	27	27	28	28	30	30	33	37	38	38	43;...
    nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    6	8	10	15	15	17	17	21	21	23	23	26	26	27	27	33	33	35	35	40	40	42	48	52];
% 37	38	38	42	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan

fly(2).flexregion = ...
    [7	11	11	13	22	24	24	28	28	26	36	41	41	46;...
    4	19	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    2	15	15	18	18	22	22	27	27	31	34	36	36	46;...
    2	7	10	13	19	22	nan	nan	nan	nan	nan	nan	nan	nan;...
    7	8	8	9	18	19	34	42	nan	nan	nan	nan	nan	nan;...
    4	5	26	30	30	31	31	32	37	40	46	49	nan	nan];

fly(3).flexregion = ...
    [9	15	nan	nan	nan	nan;...
    14	21	39	41	41	42;...
    14	21	nan	nan	nan	nan;...
    10	14	25	26	26	27];

fly(4).flexregion = [17	20	30	34	38	40];
    

fly(5).flexregion = ...
    [6	7	7	13	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    6	12	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    2	3	3	10	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    4	6	8	12	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    10	43	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    11	13	13	14	14	21	23	25	33	43	nan	nan	nan	nan	nan	nan;...
    12	19	24	30	30	43	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    2	42	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    2	7	7	14	17	20	24	26	26	29	29	33	33	36	36	39;...
    3	8	8	14	30	34	34	39	nan	nan	nan	nan	nan	nan	nan	nan;...
    4	5	7	16	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    3	8	8	19	19	24	24	28	33	38	38	41	nan	nan	nan	nan;...
    4	9	12	31	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    7	11	12	26	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    5	13	22	24	24	28	28	43	nan	nan	nan	nan	nan	nan	nan	nan;...
    6	25	35	43	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    3	7	7	21	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    5	16	21	34	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    2	11	11	28	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan	nan;...
    5	7	7	8	8	9	13	17	21	25	nan	nan	nan	nan	nan	nan];


fly(6).flexregion = ...
    [8	14	14	19;...
    3	39	nan	nan;...
    10	30	nan	nan;...
    13	28	nan	nan;...
    15	42	nan	nan];



































