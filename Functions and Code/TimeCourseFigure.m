
function [fig, InputData] = TimeCourseFigure(InputData, fig, smoothing)
% [fig, InputData] = TimeCourseFigure(InputData, fig, smoothing)
% InputData needs to include:
% x - x data == time data (e.g. -2:2)
% y - y data == speed or rotational velocity over the course of a trial
% Color == color values (rgb) for each type of data in InputData
% smoothing = number for moving average...1 = no smoothing
% EXAMPLE:
% cond = 1; %control for now
% coloropts = {'black', 'mediumblue', 'darkcyan', 'coral'}; %stationary, walking, grooming, other
% for ii = 1:4 %ii = a type of data to plot together, ex: control & stim
%     InputData(ii).x = [];
%     InputData(ii).y = [];
%     InputData(ii).Ind = 1; % 1 means that data won't graph becuase there
%     isn't any data for that instance
%     InputData(ii).Color = Color(coloropts{ii});
%     InputData(ii).Style = ':'; 
% end
% 
% for kk = 1:num.fly
%   for rep = 1:num.reps
%     % filter:
%     filter = group(kk).STATE(cond,rep);
%     % data:
%     x = (-2:1/30:2)';
%     y = [fly(kk).Control.speed(cond).data(2:end, rep); fly(kk).Stim.speed(cond).data(:, rep)];
%     InputData(filter).x(:,InputData(filter).Ind) = x;
%     InputData(filter).y(:,InputData(filter).Ind) = y;
%     % set the next position
%     InputData(filter).Ind = InputData(filter).Ind+1; 
%   end
% end

if nargin == 2
    smoothing = 5;
end

[~, ind] = max(size(InputData(1).y));
if ind == 2
    dir = 1;
elseif ind == 1
    dir = 2;
end


hold all
for ii = 1:length(InputData)
    InputData(ii).xavg = nanmean(InputData(ii).x,2);
    InputData(ii).yavg = nanmean(InputData(ii).y,dir);
    InputData(ii).yerr = sem(InputData(ii).y,dir);
    InputData(ii).ysmooth = smooth(InputData(ii).yavg, smoothing);
    if isfield(InputData, 'Style')
        LStyle = InputData(ii).Style;
    else LStyle = '-';
    end
    
    if InputData(ii).Ind > 1
%         meanerr = nanmean(InputData(ii).yerr);
%         fprintf(['  Err: ' num2str(meanerr) '\n'])
        fill_data = error_fill(InputData(ii).xavg, InputData(ii).ysmooth, InputData(ii).yerr);
                h = fill(fill_data.X, fill_data.Y, InputData(ii).Color, 'EdgeColor','none');
                set(h, 'facealpha', 0.2)
        plot(InputData(ii).xavg, InputData(ii).ysmooth, 'LineStyle', LStyle, 'color', InputData(ii).Color,...
                 'LineWidth', 1);
    end
    
    
        
end     













