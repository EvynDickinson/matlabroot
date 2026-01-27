

%% Speed in outer ring

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

formatFig(fig,blkbgd)
legend(region_name)
ylabel('average speed (mm/s)')
xlabel('time (min)')

