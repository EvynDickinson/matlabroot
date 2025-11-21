
function fig = getfig(name,secondscreen,new_size)
% fig = getfig(name,secondscreen,figsize)


% Get monitor positions
monitorPositions = get(0, 'MonitorPositions');

% Check if there is more than one monitor
if size(monitorPositions, 1) < 2
    single_monitor = true;
else
    single_monitor = false;
    if strcmp(getenv('COMPUTERNAME'), '')||strcmp(getenv('COMPUTERNAME'), 'DENALI') %mac laptop -- display on computer as second display
        idx = find(monitorPositions(:,1)==1);
        secondscreen_pos = monitorPositions(idx(1),1:2)+[1,50]; %50 point offset for home bar on the screen
    else % other computers with 2 monitors
        idx = find(monitorPositions(:,1)<0);
        secondscreen_pos = monitorPositions(idx(1),1:2)+[1,50]; %50 point offset for home bar on the screen
        if strcmp(version('-release'),'2025b') % updates for new version of figure positioning
            secondscreen_pos(2) = 50;
        end
    end
end


 % [-1069 15 1064 680]

primary_pos = [50,50]; %position from left, position from bottom
primary_size = [1450, 900];
second_size = [1064 680];

% defaults
fig_pos = primary_pos; %hold over from previous verions
figsize = primary_size;

% set name option
if nargin==0
    name = 'Fancy Figure';
end

% set second screen options
if nargin>1
    if ~secondscreen || single_monitor% primary screen selected
        fig_pos = primary_pos;
        figsize = primary_size;
    elseif secondscreen 
        fig_pos = secondscreen_pos;
        % switch getenv('COMPUTERNAME')
        %     case 'DENALI'
        %         fig_pos = [2022 161];
        %     case 'TOGIAK'
        %         fig_pos = [-1044 261];
        %     case 'EVYNPC'
        %         fig_pos = [-1070 220];
        %      case '' %mac laptop
        %         fig_pos = [-1044 261];
        %     case 'ACADIA'
        %         fig_pos = [-1074 562];
        % end
        figsize = second_size;
    end
end

% if there is a specific size, update that here:
if nargin==3
    figsize = new_size;
end

fig = figure; 
    set(fig, 'color', 'w', 'pos', [fig_pos figsize],'name', name) 
    box off


