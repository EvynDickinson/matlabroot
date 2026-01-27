

%% Speed across regions timecourse

% pull out when flies in outer ring = 1
% pull out speed

region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'};
color_list = {'MetroPurple','MetroRed','MetroOrange'};

fig = getfig;
hold on
for ii = 1:length(region_name)
    current_var = data.(region_name{ii});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    speed(~loc) = nan;
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    m = mean(speed_all,2,'omitnan');    
    plot(data.time,smooth(m,300,'moving'),"Color",Color(color_list{ii}))
end

formatFig(fig,blkbgd);
legend(region_name)
ylabel('average speed (mm/s)')
xlabel('time (min)')

%% Speed across regions within temp regimes scatter plot

% pull out when flies are in each region
% pull out speed for each temp regime time period

region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'};
regime_name = {'WT','WS','CT','CS','SS'};
color_list = {'MetroPurple','MetroRed','MetroOrange'};

for ii = 1:length(region_name)
    current_var = data.(region_name{ii});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    speed(~loc) = nan;
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    for regime = 1:length(regime_name)
    current_reg = data.tempbin.(regime_name{regime});
    loc2 = replaceNaN(current_reg,0);
    speedR = speed_all;
    speedR(~loc2) = nan;

end

