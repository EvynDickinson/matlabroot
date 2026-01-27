

function [title_str,pName,y_dir,y_lab,nullD,scaler,dType,fig_dir,sexSep,ylimits] = PlotParamSelectionHR(plotType,dataType,multiselect)
% [title_str, pName,y_dir,y_lab,nullD,scaler,dType,fig_dir,ext,ylimits] = PlotParamSelectionHR(plotType,dataType,multiselect)
% 
% PURPOSE
% User selects a HIGH resolution data type that returns parameters specific 
% to that type of data which are related to plotting the data
%
% INPUTS (all optional)
%     'plotType' : select the type of data to plot (avg, single trial, separated heating and cooling)
%           true = user selects the type of data to plot
%           false = the type of data is not returned
%           (default : average)
%     'dataType' : What subsection (or all) data types only [all, location,courtship]
%           'location' = only spatially based parameters are desired (eg. not speed or sleep)
%           'all' = all data types [default value]
%           'courtship' = only courtship behaviors
%     'mutliselect' : user ability to select multiple data types
%           true = more than one parameter can be selected for plotting
%           false = returns only a single parameter (default value)
%
% OUTPUTS
%     'title_str' : Name of the selected parameter
%           (e.g. 'Distance to Food', 'Food Occupancy', 'Food Circle', 'Occupancy')
%     'pName' : parameter name in the grouped structure
%     'y_dir' : y axis plotting direction ('normal' vs 'reverse')
%     'y_lab' : y axis label 
%           (e.g. 'Distance to Food (mm)', 'Food Occupancy (% flies)')
%     'nullD' : null distribution value (e.g. 25 for 25% of the arena)
%     'scaler' : multiplication scaler for the parameter (usually 1 or 100)
%     'dType' : display type (e.g., single trial, avg, sep H &C)
%     'fig_dir' : fig directory ending for subfolder of this plot type [string]
%     'sexSep' : is the data separated by sex (i.e. is the data structure
%           [time, sex, trial] format or [time, trial] format
%     'ylimits' : suggested y limits for this parameter
%
% ES DICKINSON, 2026

%%
% allow multiple parameters to be selected at a time or not
if exist('multiselect', 'var') && multiselect 
    selectionMode = 'multiple';
    prompt_str = 'Select MULTPLE types for comparison:';
else
    selectionMode = 'single';
    prompt_str = 'Select a data type:';
end

% Select the type of information to plot: 
if ~exist('dataType', 'var') 
    dataType = 'all';  % default to all options
end
switch dataType
    case 'all'
        paramList = {'Outer Ring', 'Sleep',  'Food Quadrant', 'Speed',...
          'Inner Food Quadrant', 'Distance to Food', 'Fly On Food', 'Turning', 'Courtship Index',...
          'Courtship Index All', 'Circling', 'Circling All', 'Chase', 'Chase All',...
          'Food Circle', 'Eccentricity','Wing Extension', 'Wing Extension All', 'Inter Fly Distance'};
    case 'location' % behavior is determined by location in the arena (percent based ones) 
        paramList = {'Outer Ring','Food Quadrant','Inner Food Quadrant', ...
          'Fly On Food','Food Circle',};
    case 'courtship' % courtship specific parameters
        paramList = {'Courtship Index','Courtship Index All','Inter Fly Distance'...
           'Circling', 'Circling All', 'Chase', 'Chase All','Wing Extension', 'Wing Extension All'};
end
idx = listdlg('ListString', paramList,'PromptString', prompt_str,...
              'ListSize',[200,400],'SelectionMode',selectionMode);
if isempty(idx)
    disp('No choice selected')
    return
else 
    numParams = length(idx);
end

% initialize empty parameters (need to do this for the cases where there
% are multiple selection options available for the data)
[title_str, pName, y_lab, nullD, scaler] = deal([]); %nan(numParams,1)
ylimits = nan(numParams,2);
sexSep = true([numParams,1]); % default is true (that each fly has their own data)
y_dir = repmat({'normal'},[numParams,1]); % set default to normal direction

% load data-specific parameters into each group
for i = 1:numParams
    title_str{i} = paramList{idx(i)};
    switch title_str{i}
        case 'Outer Ring'
            pName{i} = 'OutterRing'; % yes, I have bad spelling sometimes lol
            scaler(i) = 100;
            y_lab{i} = 'edge occupancy (% flies)';
            ylimits(i,:) = [0, 50];
            nullD(i) = 25;
        case 'Speed'
            pName{i} = 'speed';
            scaler(i) = 1;
            y_lab{i} = 'speed (cm/s)';
            ylimits(i,:) = [0 18];
        case 'Food Quadrant'
            pName{i} = 'foodQuad';
            scaler(i) = 100;
            y_lab{i} = 'food quadrant (% flies)';
            ylimits = [0 90];
        case 'Inner Food Quadrant'
            pName{i} = 'innerFoodQuad';
            scaler(i) = 100;
            y_lab{i} = 'inner food quadrant (% flies)';
            ylimits = [0 90];
            nullD(i) = 18.75;
        case 'Food Circle'
            pName{i} = 'foodcircle';
            scaler(i) = 100;
            y_lab{i} = 'food circle (% flies)';
            ylimits = [0 65];
            nullD(i) = 7; %TODO check if this is the correct null dist percent
        case 'Sleep'
            pName{i} = 'sleep';
            scaler(i) = 100;
            y_lab{i} = 'sleeping (% flies)';
            ylimits = [0 60];
        case 'Courtship Index'
            pName{i} = 'CI';
            scaler(i) = 100;
            y_lab{i} = 'restrictive courtship index (% male flies)';
            ylimits = [0 60];
            sexSep(i) = false; 
        case 'Courtship Index All'
            pName{i} = 'CI_all';
            scaler(i) = 100;
            y_lab{i} = 'courtship index all (% male flies)';
            ylimits(i,:) = [0, 50];
            sexSep(i) = false; 
        case 'Eccentricity'
            pName{i} = 'Eccentricity';
            scaler(i) = 1;
            y_lab{i} = 'distance from center (mm)'; % TODO check if this is distance from edge
            % ylimits(i,:) = [0, 50];
        case 'Distance to Food'
            pName{i} = 'dist2food';
            scaler(i) = 100;
            y_lab{i} = 'distance to food (mm)';
            ylimits(i,:) = [0, 35];
            y_dir{i} = 'reverse';
        case 'Flies on Food'
            pName{i} = 'FlyOnFood';
            scaler(i) = 100;
            y_lab{i} = 'Flies on Food (% flies)';
            % ylimits(i,:) = [0, 35];
            % nullD(i) = 25; % TODO : could actually calculate spatial distribution of the food relative to the arena size
        case 'Circling'
            pName{i} = 'circling_1sec';
            scaler(i) = 100;
            y_lab{i} = 'circling behavior (% male flies)';
            % ylimits(i,:) = [0, 35];
            sexSep(i) = false; 
        case 'Circling All'
            pName{i} = 'circling_all';
            scaler(i) = 100;
            y_lab{i} = 'all circling behavior (% male flies)';
            % ylimits(i,:) = [0, 35];
            sexSep(i) = false; 
        case 'Chase'
            pName{i} = 'court_chase';
            scaler(i) = 100;
            y_lab{i} = 'chase behavior (% male flies)';
            % ylimits(i,:) = [0, 35];
            sexSep(i) = false; 
        case 'Chase All'
            pName{i} = 'chase_all';
            scaler(i) = 100;
            y_lab{i} = 'all chase behavior (% male flies)';
            % ylimits(i,:) = [0, 35];
            sexSep(i) = false; 
        case 'Wing Extension'
            pName{i} = 'wing_ext';
            scaler(i) = 100;
            y_lab{i} = 'wing extension behavior (% male flies)';
            % ylimits(i,:) = [0, 35];
            sexSep(i) = false;         
        case 'Wing Extension All'
            pName{i} = 'wing_ext_all';
            scaler(i) = 100;
            y_lab{i} = 'wing extension behavior (% male flies)';
            % ylimits(i,:) = [0, 35];
            sexSep(i) = false; 
        case 'Inter Fly Distance'
            pName{i} = 'IFD';
            scaler(i) = 1;
            y_lab{i} = 'distance between flies (mm)';
            % ylimits(i,:) = [0, 35];
            sexSep(i) = false; 
        case 'Turning'
            pName{i} = 'turning speed';
            scaler(i) = 1;
            y_lab{i} = 'turning speed (mm/s)'; % TODO check that this is correct
            % ylimits(i,:) = [0, 35];
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
if plotType % if choice of options is desired
    qList = { 'Heating and Cooling','Average', 'Single trial lines'};
    idx = listdlg('ListString', qList, PromptString='How do you want to plot the data:',...
        ListSize=[300,200], SelectionMode=selectionMode);
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
else % default option is average lines
    dType = 2;
    fig_dir = '';
end

