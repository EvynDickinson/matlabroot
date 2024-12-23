
clear; clc
[excelfile, Excel, XL] = load_QuadBowlExperiments;

% find all the experiment start times to polar plot:
timeList = [];
timeList(:,1) = [excelfile{2:end, Excel.zeitgebertime}];

histogram(timeList)

% Convert the 24 hours into a circular radian value
theta = (24-(timeList/24)) * 2 * pi; 

% show the experiment start times as a polar plot
fig = getfig('Polar Histogram', 1);
    h = polarhistogram(theta);
    set(h,'FaceColor',Color('teal'), 'EdgeColor', 'k','FaceAlpha',0.8)
    % set labels and formats
    ax = gca;
    set(ax, 'ThetaZeroLocation', 'top')
    set(ax, 'ThetaTick',0:45:359,'ThetaTickLabel', {'0','21','18','15','12','9','6','3'})

% OVERLAY NIGHT SHADING
ax_cart = axes('Position', ax.Position, 'Color', 'none');  % Create an overlay Cartesian axes

% Step 5: Shade the left half (90° to 270°)
hold(ax_cart, 'on');

% Define the angles and radius for shading the left half
shading_theta = linspace(pi/2, 3*pi/2, 100);  % Angles from 90° to 270°
shading_r = max(ax.RLim) * ones(size(shading_theta));  % Use maximum radius of polar plot

% Convert polar coordinates to Cartesian for the shading
[x_shade, y_shade] = pol2cart(shading_theta, shading_r);

% Add the center point to complete the patch
x_shade = [0 x_shade 0];  % Include the origin (0,0)
y_shade = [0 y_shade 0];

% Create the filled patch in grey
fill(ax_cart, x_shade, y_shade, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Grey with transparency

% Adjust Cartesian axis limits to match polar plot
axis(ax_cart, 'equal');  % Ensure equal scaling on both axes
r_max = max(ax.RLim);  % Get the maximum radius from the polar plot
ax_cart.XLim = [-r_max r_max];
ax_cart.YLim = [-r_max r_max];

% Remove the Cartesian axis ticks
ax_cart.XTick = [];
ax_cart.YTick = [];
ax_cart.XColor = 'none';
ax_cart.YColor = 'none';

% Bring the polar plot to the front  %TODO: see if there is a way to do
% this without blocking the new shading
uistack(ax, 'top');


% TODO: update the label sizes and add some titles etc, 
% TODO: try making the binsizes for each hour
