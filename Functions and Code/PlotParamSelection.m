

function [title_str,pName,y_dir,y_lab,nullD,scaler,dType,fig_dir,ext] = PlotParamSelection(plotType,location_only,multiselect)
% [title_str, pName,y_dir,y_lab,nullD,scaler,dType,fig_dir,ext] = PlotParamSelection(plotType,location_only,multiselect)
% 
% PURPOSE
% User selects a low resolution data type that returns parameters specific 
% to that type of data which are related to plotting the data
%
% INPUTS (all optional)
%   'plotType' : select the type of data to plot (avg, single trial, separated heating and cooling)
%       true = user selects the type of data to plot
%       false = the type of data is not returned
%       (default : average)
%   'location_only' : restrict the plotting data types to only those that
%       are based on spatial location of the flies
%       true = only spatially based parameters are desired (eg. not speed or sleep)
%       false = all data types are available
%    'mutliselect' : user ability to select multiple data types
%       true = more than one parameter can be selected for plotting
%       false = returns only a single parameter 
%
% OUTPUTS
%    'title_str' : Name of the selected parameter
%       (e.g. 'Distance to Food', 'Food Occupancy', 'Food Circle', 'Occupancy')
%    'pName' : parameter name in the grouped structure
%    'y_dir' : y axis plotting direction ('normal' vs 'reverse')
%    'y_lab' : y axis label 
%       (e.g. 'Distance to Food (mm)', 'Food Occupancy (% flies)')
%    'nullD' : null distribution value (e.g. 25 for 25% of the arena)
%    'scaler' : multiplication scaler for the parameter (usually 1 or 100)
%    'dType' : display type (e.g., single trial, avg, sep H &C)
%    'fig_dir' : fig directory ending for subfolder of this plot type [string]
%    'ext' : are there subgroups of the data that can be plotted [logical]
%       (e.g., 'innerquad' has 4 extension subquadrants for each location)
%
% ES DICKINSON, 2024

%%
% allow multiple parameters to be selected at a time or not
if nargin > 2 && multiselect
    selectionMode = 'multiple';
    prompt_str = 'Select MULTPLE types for comparison:';
else
    selectionMode = 'single';
    prompt_str = 'Select a data type:';
end

% Select the type of information to plot: 
if nargin>1 && location_only
    paramList = {'fullquad','quadring', 'innerquad','circle10', ...
                 'circle7', 'circle5', 'ring', 'inner75'};
else
    paramList = {'fullquad', 'innerquad', 'quadring','circle10', 'circle7',...
                 'circle5', 'fliesonfood','ring','speed','sleep','inner75'};
end

idx = listdlg('ListString', paramList,'PromptString', prompt_str,...
              'ListSize',[200,200],'SelectionMode',selectionMode);
if isempty(idx)
    disp('No choice selected')
    return
else 
    numParams = length(idx);
end

% initialize empty parameters (need to do this for the cases where there
% are multiple selection options available for the data)
[title_str, pName, ext, y_dir, y_lab, nullD, scaler] = deal([]);

% load data-specific parameters into each group
for i = 1:numParams
    title_str{i} = paramList{idx(i)};
    switch title_str{i}
        case 'ring'
            pName{i} = 'ring';
            ext(i) = false; % no sub regions
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 25;
            scaler(i) = 1;
        case 'inner75'
            pName{i} = 'inner75';
            ext(i) = false; % no sub regions
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 75;
            scaler(i) = 1;
        case 'fullquad'
            pName{i} = 'fullquad';
            ext(i) = true; % extension for sub region required
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 25;
            scaler(i) = 1;
         case 'quadring'
            pName{i} = 'quadring';
            ext(i) = true; % extension for sub region required
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 6.25;
            scaler(i) = 1;
        case 'innerquad'
            pName{i} = 'innerquad';
            ext(i) = true; % extension for sub region required
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 18.75;
            scaler(i) = 1;
        case 'circle10'
            pName{i} = 'circle10';
            ext(i) = true; % extension for sub region required
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 10;
            scaler(i) = 1;
        case 'circle7'
            pName{i} = 'circle7';
            ext(i) = true; % extension for sub region required
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 7;
            scaler(i) = 1;
        case 'circle5'
            pName{i} = 'circle5';
            ext(i) = true; % extension for sub region required
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 5;
            scaler(i) = 1;
        case 'speed'
            pName{i} = 'speed';
            ext(i) = false;
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (mm/s)'];
            nullD(i) = nan;
            scaler(i) = 1;
        case 'sleep'
            pName{i} = 'sleep';
            ext(i) = false;
            y_dir{i} = 'normal';
            y_lab{i} = 'Sleeping flies (%)';
            nullD(i) = nan;
            scaler(i) = 1;
        case 'fliesonfood'
            pName{i} = 'fliesonfood';
            ext(i) = false;
            y_dir{i} = 'normal';
            y_lab{i} = 'flies on food (#)';
            nullD(i) = nan;
            scaler(i) = 1;
    end
end

% reset to not be in a cell if only one parameter is selected
if numParams==1
    i = 1;
    pName = pName{i};
    y_dir = y_dir{i};
    y_lab = y_lab{i};
    title_str = title_str{i};
end

% Select the type of data that will be plotted (e.g. avg, single, H&C)
if ~exist('plotType','var')
    plotType = false;
end
if plotType
    qList = { 'Heating and Cooling','Average', 'Single trial lines'};
    idx = listdlg('ListString', qList,'PromptString', 'How do you want to plot the data:','ListSize',[300,200]);
    if isempty(idx)
        disp('No choice selected')
        return
    end
    typeString = qList{idx};
    switch typeString
        case 'Single trial lines'
            dType = 1;
            fig_dir = 'all trial lines/';
        case 'Average'
            dType = 2;
            fig_dir = 'time course/';
        case 'Heating and Cooling'
            dType = 3;
            fig_dir = 'time course H and C/';
         case ''
            disp('')
            return
    end
else
    dType = 2;
    fig_dir = '';
end

%% 


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

