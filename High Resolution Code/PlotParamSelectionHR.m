

function [title_str,pName,y_dir,y_lab,nullD,scaler,dType,fig_dir,ext,ylimits] = PlotParamSelectionHR(plotType,dataType,multiselect)
% [title_str, pName,y_dir,y_lab,nullD,scaler,dType,fig_dir,ext,ylimits] = PlotParamSelectionHR(plotType,dataType,multiselect)
% 
% PURPOSE
% User selects a HIGH resolution data type that returns parameters specific 
% to that type of data which are related to plotting the data
%
% INPUTS (all optional)
%   'plotType' : select the type of data to plot (avg, single trial, separated heating and cooling)
%       true = user selects the type of data to plot
%       false = the type of data is not returned
%       (default : average)
%   'dataType' : What subsection (or all) data types only
%       'location' = only spatially based parameters are desired (eg. not speed or sleep)
%       'all' = all data types [default value]
%       'courtship' = only courtship behaviors
%    'mutliselect' : user ability to select multiple data types
%       true = more than one parameter can be selected for plotting
%       false = returns only a single parameter (default value)
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
%    'ylimits' : suggested y limits for this parameter
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
if exist('dataType', 'var') 
    switch dataType
        case 'all'
            paramList = {'Outer Ring', 'Sleep',  'Food Quadrant', 'Speed',...
              'Inner Food Quadrant', 'Distance to Food', 'Fly On Food', ...
              'Circling', 'Circling All', 'Chase', 'Chase All',...
              'Food Circle', 'Eeccentricity', 'CI','Wing Extension', 'Wing Extension All'};
        case 'location'
            paramList = {'Outer Ring', 'foodQuad', 'innerFoodQuad',...
              'FlyOnFood', 'foodcircle'};
        case 'courtship'
            paramList = {'CI','CI_all', 'circling_1sec', 'circling_all',...
              'court_chase', 'chase_all', 'wing_ext', 'wing_ext_all'};
    end
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
[title_str, pName, ext, y_dir, y_lab, nullD, scaler] = deal([nan(numParams,1)]);
ylimits = nan(numParams,2);

% load data-specific parameters into each group
for i = 1:numParams
    title_str{i} = paramList{idx(i)};
    switch title_str{i}
        case 'Outer Ring'
            pName{i} = 'OutterRing';
            scaler(i) = 100;
            y_lab{i} = 'edge occupancy (% flies)';
            ylimits(i,:) = [0, 50];
            nullD(i) = 25;
            ext(i) = false; % no extended subregions
        case 'Speed'
            pName{i} = 'speed';
            scaler(i) = 1;
            ylabel_str = 'speed (cm/s)';
            ylimits(i,:) = [0 18];
        case 'foodQuad'
            scaler(i) = 100;
            ylabel_str = 'food quadrant (% flies)';
            ylimits = [0 90];
        case 'innerFoodQuad'
            scaler(i) = 100;
            ylabel_str = 'inner food quadrant (% flies)';
            ylimits = [0 90];
        case 'foodcircle'
            scaler(i) = 100;
            ylabel_str = 'food circle (% flies)';
            ylimits = [0 65];
        case 'sleep'
            scaler(i) = 100;
            ylabel_str = 'sleeping (% flies)';
            ylimits = [0 60];
        case 'CI'
            scaler(i) = 100;
            ylabel_str = 'sleeping (% flies)';
            ylimits = [0 60];
        case 'Eccentricity'
            pName{i} = 'OutterRing';
            scaler(i) = 100;
            y_lab{i} = 'edge occupancy (% flies)';
            ylimits(i,:) = [0, 50];
            nullD(i) = 25;
            ext(i) = false; 
        case 
            pName{i} = 'OutterRing';
            scaler(i) = 100;
            y_lab{i} = 'edge occupancy (% flies)';
            ylimits(i,:) = [0, 50];
            nullD(i) = 25;
            ext(i) = false; 
        case 
            pName{i} = 'OutterRing';
            scaler(i) = 100;
            y_lab{i} = 'edge occupancy (% flies)';
            ylimits(i,:) = [0, 50];
            nullD(i) = 25;
            ext(i) = false; 
        case 
            pName{i} = 'OutterRing';
            scaler(i) = 100;
            y_lab{i} = 'edge occupancy (% flies)';
            ylimits(i,:) = [0, 50];
            nullD(i) = 25;
            ext(i) = false; 


            
            pName{i} = 'ring';
            ext(i) = false; % no sub regions
            y_dir{i} = 'normal';
            y_lab{i} = [title_str{i} ' (%)'];
            nullD(i) = 25;
            scaler(i) = 1;




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

