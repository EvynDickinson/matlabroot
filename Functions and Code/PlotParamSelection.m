

function [title_str, pName,y_dir,y_lab,nullD,scaler,dType,fig_dir,ext] = PlotParamSelection(plotType,location_only)
% TODO: update this for the new region occupancy structures...6.16.25 ESD
% [pName,y_dir,y_lab,nullD,scaler] = PlotParamSelection(plotType);
% plotType = true --> select the type of data to plot (e.g., avg, single
% trial, separated heating and cooling
% {'Distance to Food', 'Food Occupancy', 'Food Circle Occupancy', 'Quadrant Occupancy', 'Ring Occupancy','Speed'};
% title_str = Name of the selected parameter
% pName = parameter name in the grouped structure
% y_dir = y axis plotting direction
% y_lab  = y axis label
% nullD = null distribution value
% scaler = multiplication scaler for the parameter 
% dType = display type (e.g., single trial, avg, sep H &C)
% fig_dir = fig directory ending for subfolder of this plot type
% 
% ES Dickinson, 2024

% Select the type of information to plot: 
if nargin>1 && location_only
    % paramList = { 'Food Occupancy', 'Food Circle Occupancy', 'Quadrant Occupancy', 'Ring Occupancy'};
    paramList = { 'ring', 'inner75', 'fullquad','quadring', 'innerquad','circle10', 'circle7', 'circle5'};
else
    % paramList = { 'Food Occupancy', 'Food Circle Occupancy', 'Quadrant Occupancy', 'Ring Occupancy','Proximity to Food','Speed','Sleep'};
    paramList = { 'ring', 'inner75', 'fullquad', 'innerquad','circle10', 'circle7', 'circle5','fliesonfood','sleep','speed'};
end

idx = listdlg('ListString', paramList,'PromptString', 'Select the type of data you want to plot:','ListSize',[200,200]);
if isempty(idx)
    disp('No choice selected')
    return
end
title_str = paramList{idx};
switch title_str
    case 'ring'
        pName = 'ring';
        ext = false; % no sub regions
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 25;
        scaler = 1;
    case 'inner75'
        pName = 'inner75';
        ext = false; % no sub regions
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 75;
        scaler = 1;
    case 'fullquad'
        pName = 'fullquad';
        ext = true; % extension for sub region required
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 75;
        scaler = 1;
     case 'quadring'
        pName = 'quadring';
        ext = true; % extension for sub region required
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 6.25;
        scaler = 1;
    case 'innerquad'
        pName = 'innerquad';
        ext = true; % extension for sub region required
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 18.75;
        scaler = 1;
    case 'circle10'
        pName = 'circle10';
        ext = true; % extension for sub region required
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 10;
        scaler = 1;
    case 'circle7'
        pName = 'circle7';
        ext = true; % extension for sub region required
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 7;
        scaler = 1;
    case 'circle5'
        pName = 'circle5';
        ext = true; % extension for sub region required
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 5;
        scaler = 1;
    % case 'Proximity to Food'
    %     pName = 'dist';
    %     y_dir = 'reverse';
    %     y_lab = [title_str ' (mm)'];
    %     nullD = 18.1;
    %     scaler = 1;
    %     % % Y limit ranges
    %     % dist_lim = [10,35];       %distance
    %     % dt_lim = [14, 32];        %distance-temp
    % case 'Food Occupancy'
    %     pName = 'occ';
    %     y_dir = 'normal';
    %     y_lab = [title_str ' (%)'];
    %     nullD = 14.4;
    %     scaler = 100;
    % case 'Food Circle Occupancy'
    %     pName = 'foodcircle';
    %     y_dir = 'normal';
    %     y_lab = [title_str ' (%)'];
    %     nullD = 10;
    %     scaler = 1;
    % case 'Quadrant Occupancy'
    %     pName = 'quadrant';
    %     y_dir = 'normal';
    %     y_lab = [title_str ' (%)'];
    %     nullD = 25;
    %     scaler = 1;
    % case 'Ring Occupancy'
    %     pName = 'ring';
    %     y_dir = 'normal';
    %     y_lab = [title_str ' (%)'];
    %     nullD = 25;
    %     scaler = 1;
    case 'speed'
        pName = 'speed';
        ext = false;
        y_dir = 'normal';
        y_lab = [title_str ' (mm/s)'];
        nullD = nan;
        scaler = 1;
    case 'sleep'
        pName = 'sleep';
        ext = false;
        y_dir = 'normal';
        y_lab = 'Sleeping flies (%)';
        nullD = nan;
        scaler = 1;
    case 'fliesonfood'
        pName = 'fliesonfood';
        ext = false;
        y_dir = 'normal';
        y_lab = 'flies on food (#)';
        nullD = nan;
        scaler = 1;
end

if plotType
    qList = {'Single trial lines', 'Average', 'Heating and Cooling'};
    idx = listdlg('ListString', qList,'PromptString', ['How do you want to plot the ' title_str ' data:'],'ListSize',[300,200]);
    if isempty(idx)
        disp('No choice selected')
        return
    end
    typeString = qList{idx};
    switch typeString
        case 'Single trial lines'
            dType = 1;
            fig_dir = '/all trial lines/';
        case 'Average'
            dType = 2;
            fig_dir = '/time course/';
        case 'Heating and Cooling'
            dType = 3;
            fig_dir = '/time course H and C/';
         case ''
            disp('')
            return
    end
else
    dType = 1;
    fig_dir = '';
end


% % Select the type of information to plot: 
% if nargin>1 && location_only
%     paramList = { 'Food Occupancy', 'Food Circle Occupancy', 'Quadrant Occupancy', 'Ring Occupancy'};
% else
%     paramList = { 'Food Occupancy', 'Food Circle Occupancy', 'Quadrant Occupancy', 'Ring Occupancy','Proximity to Food','Speed','Sleep'};
% end
% 
% idx = listdlg('ListString', paramList,'PromptString', 'Select the type of data you want to plot:','ListSize',[200,200]);
% if isempty(idx)
%     disp('No choice selected')
%     return
% end
% title_str = paramList{idx};
% switch title_str
%     case 'Proximity to Food'
%         pName = 'dist';
%         y_dir = 'reverse';
%         y_lab = [title_str ' (mm)'];
%         nullD = 18.1;
%         scaler = 1;
%         % % Y limit ranges
%         % dist_lim = [10,35];       %distance
%         % dt_lim = [14, 32];        %distance-temp
%     case 'Food Occupancy'
%         pName = 'occ';
%         y_dir = 'normal';
%         y_lab = [title_str ' (%)'];
%         nullD = 14.4;
%         scaler = 100;
%     case 'Food Circle Occupancy'
%         pName = 'foodcircle';
%         y_dir = 'normal';
%         y_lab = [title_str ' (%)'];
%         nullD = 10;
%         scaler = 1;
%     case 'Quadrant Occupancy'
%         pName = 'quadrant';
%         y_dir = 'normal';
%         y_lab = [title_str ' (%)'];
%         nullD = 25;
%         scaler = 1;
%     case 'Ring Occupancy'
%         pName = 'ring';
%         y_dir = 'normal';
%         y_lab = [title_str ' (%)'];
%         nullD = 25;
%         scaler = 1;
%     case 'Speed'
%         pName = 'speed';
%         y_dir = 'normal';
%         y_lab = [title_str ' (mm/s)'];
%         nullD = nan;
%         scaler = 1;
%     case 'Sleep'
%         pName = 'sleep';
%         y_dir = 'normal';
%         y_lab = 'Sleeping flies (%)';
%         nullD = nan;
%         scaler = 1;
% end
% 
% if plotType
%     qList = {'Single trial lines', 'Average', 'Heating and Cooling'};
%     idx = listdlg('ListString', qList,'PromptString', ['How do you want to plot the ' title_str ' data:'],'ListSize',[300,200]);
%     if isempty(idx)
%         disp('No choice selected')
%         return
%     end
%     typeString = qList{idx};
%     switch typeString
%         case 'Single trial lines'
%             dType = 1;
%             fig_dir = '/all trial lines/';
%         case 'Average'
%             dType = 2;
%             fig_dir = '/time course/';
%         case 'Heating and Cooling'
%             dType = 3;
%             fig_dir = '/time course H and C/';
%          case ''
%             disp('')
%             return
%     end
% else
%     dType = 1;
%     fig_dir = '';
% end

