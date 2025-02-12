

function [title_str, pName,y_dir,y_lab,nullD,scaler,dType,fig_dir] = PullParamProperties(title_str,typeString)
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

switch title_str
    case 'Proximity to Food'
        pName = 'dist';
        y_dir = 'reverse';
        y_lab = [title_str ' (mm)'];
        nullD = 18.1;
        scaler = 1;
        % % Y limit ranges
        % dist_lim = [10,35];       %distance
        % dt_lim = [14, 32];        %distance-temp
    case 'Food Occupancy'
        pName = 'occ';
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 14.4;
        scaler = 100;
    case 'Food Circle Occupancy'
        pName = 'foodcircle';
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 10;
        scaler = 1;
    case 'Quadrant Occupancy'
        pName = 'quadrant';
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 25;
        scaler = 1;
    case 'Ring Occupancy'
        pName = 'ring';
        y_dir = 'normal';
        y_lab = [title_str ' (%)'];
        nullD = 25;
        scaler = 1;
    case 'Speed'
        pName = 'speed';
        y_dir = 'normal';
        y_lab = [title_str ' (mm/s)'];
        nullD = nan;
        scaler = 1;
    case 'Sleep'
        pName = 'sleep';
        y_dir = 'normal';
        y_lab = 'Sleeping flies (%)';
        nullD = nan;
        scaler = 1;
end


    switch typeString
        case 'Single trial lines'
            dType = 1;
            fig_dir = '/tuning curves all trial lines/';
        case 'Average'
            dType = 2;
            fig_dir = '/tuning curves/';
        case 'Heating and Cooling'
            dType = 3;
            fig_dir = '/tuning curves H and C/';
         case ''
            disp('')
            return
    end
end




