
function fill_data = error_fill(xdata, ydata, error)
% fill_data = error_fill(xdata, ydata, error)
% 
% This takes x,y and y's error data and sorts it into the proper format for 
% using the fill function to plot error filled areas around a line plot
% Returns data for making a filled shaped of an error region
% ** Is executed within the function 'plot_error_fill.m'  **
% 
% INPUTS
%    'xdata' : [vector with x-axis data]
%    'ydata' : [vector with y-axis data]
%    'error' : [SEM or STD associated with ydata, in vector format]
% 
% OUTPUT
%    'fill_data' [struct with fields 'X' and 'Y' for the fill function]
%
% EXAMPLE
%    fill_data = error_fill(xdata, ydata, error) [this function]
%    h = fill(fill_data.X, fill_data.Y, Color(color), 'EdgeColor','none');
%    set(h, 'facealpha', 0.2)
%     
% ES Dickinson, University of Washington, Jan 2019

%%  all data must be set out along the long dimension--in columns, not
% within a row
% X DATA

if size(xdata,2)==1
    X = xdata';
else
    X = xdata;
end

% Y DATA
if size(ydata,2)==1
    Y.data = ydata';
else
    Y.data = ydata;
end

% Y ERROR
if size(error,2)==1
    Y.error = error';
else
    Y.error = error;
end

% Error Fill Areas
y1_1 = (Y.data+Y.error); % error high
y1_2 = (Y.data-Y.error); % error low
fill_data.X = [fliplr(X), (X)];
fill_data.Y = [fliplr(y1_1), (y1_2)];

end



