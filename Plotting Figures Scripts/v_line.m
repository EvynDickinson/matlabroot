function v_line(x, in1, in2, in3)
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
        if ischar(in1) %single color only
            % letter code for color:
            if length(in1)==1 %single value
                RGB = in1;
            else
                RGB = Color(in1); 
            end
        elseif iscell(in1) %multiple colors in letter form
            temp = in1{ii}; %pull the one associated with this ii
            % convert color:
            if length(temp)==1 %single value
                RGB = temp;
            else
                RGB = Color(temp); 
            end
        elseif isnumeric(in1)
            if size(in1,1)>1
                RBG = in1(ii,:);
            else
                RBG = in1;
            end
        end
    end
    
    
    %linestyle input
    if nargin>2 
        if ischar(in2)
            linetype = in2;
        elseif iscell(in2)
            linetype = in2{ii};
        end
    end
    
    %linewidth input
    if nargin==4
        if length(in3)>1
            LW = in3(ii);
        else
            LW = in3;
        end
    end
    
    plot([xval, xval], y,'color', RGB, 'linestyle', linetype, 'linewidth', LW)
end  
    
     
% return hold status
if g == 0
    hold off
end









end
