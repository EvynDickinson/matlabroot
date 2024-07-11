function v_line(x, line_color, line_style, line_width)
% v_line(x, color, style, linewidth)
% 
% Allows for multiple x inputs for vertical lines
% Each can have same or diff line colors and styles
%
% ES Dickinson
% Yale University, 2020
%
% Adaptaton of 'vline' by Brandon Kuczenski
% brandon_kuczenski@kensingtonlabs.com

% %type of color input allowable:
% a = 'w';
% A = 'white';
% aa = {'w', 'b'};
% AA = {'w', 'blue'};
% aaa = [0 1 0];
% AAA = [0 1 0; 0 0 0];

% type of linestyle inputs:
% ':'
% {':','--','-'}


% turn hold to on:
g = ishold(gca);
hold on
%get ylimits
y = get(gca,'ylim');


% ADD LINES

% defaults:
LW = 0.5;
linetype = '--';
RGB = 'r';

% user inputs:
for ii = 1:length(x)
    % x value for vline
    xval = x(ii); 
    
    % find line color 
    if nargin>1 
        if ischar(line_color) %single color only
            % letter code for color:
            if length(line_color)==1 %single value
                RGB = line_color;
            else
                RGB = Color(line_color); 
            end
        elseif iscell(line_color) %multiple colors in letter form
            temp = line_color{ii}; %pull the one associated with this ii
            % convert color:
            if length(temp)==1 %single value
                RGB = temp;
            else
                RGB = Color(temp); 
            end
        elseif isnumeric(line_color)
            if size(line_color,1)>1
                RBG = line_color(ii,:);
            else
                RBG = line_color;
            end
        end
    end
    
    %linestyle input
    if nargin>2 
        if ischar(line_style)
            linetype = line_style;
        elseif iscell(line_style)
            linetype = line_style{ii};
        end
    end
    
    %linewidth input
    if nargin==4
        if length(line_width)>1
            LW = line_width(ii);
        else
            LW = line_width;
        end
    end
    
    plot([xval, xval], y,'color', RGB, 'linestyle', linetype, 'linewidth', LW)
end  
    
     
% return hold status
if g == 0
    hold off
end









end
