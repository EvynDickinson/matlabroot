
function color_output = Color(color_1, color_2, N)
% 
% color_matrix = Color(color_1, color_2, N)
%
% SINGLE COLOR:
% Get RGB color vector from a color name based on the options in the
% 'Color Palette.xlsx' excel file
% OR 
% MULTICOLOR:
% Get a smooth RGB color vectors gradient (of length N) between two colors pulled from
% the options in the 'Color Palette.xlsx' excel file
%
% INPUTS
%   'color_1' : color name in string format, e.g. 'Grey'
%           * for color gradients, this is the start color *
%   'color_2' : end color for the color gradient, input as a string 
%
% OUTPUT
%    'color_matrix' : matrix of color values (R B G) for each color 
%           columns are the R, G, B value for a given color
%           rows are each a different color in the gradient
%
% EXAMPLE:
% e.g. [1, 0.5, 0] = Color('orange')
%
% ** requires the 'color_palette.m' file to be in the current matlab path**
% 
% ES DICKINSON, Dec 2018

%% Run the function

load('color_palette','colors')

% Only one color requested:
color_loc(:,1) = strcmpi(color_1, colors.names);
color_matrix(1,:) = colors.rgb(color_loc(:,1), :);
color_output(1,:) = color_matrix;
% Two colors requested:
if nargin > 1
    color_loc(:,2) = strcmpi(color_2, colors.names);
    color_matrix(2,:) = colors.rgb(color_loc(:,2), :);
% Color gradient requested:
    if nargin == 3
      %set up the gradient
        interval = (color_matrix(1,:)-color_matrix(2,:))/(N-1);
        for ii = 1:N-1 %RGB   
            color_output(ii+1,:) = abs(color_output(ii,:)-interval);
        end  
    end
end

%% HOW THE COLORS ARE UPDATED|CREATED:
% xlFile = 'C:\Users\evynd\Documents\matlabroot\Functions and Code\Color Palette.xlsx';      
% %load excel sheet data
% [~,~,excelfile] = xlsread(xlFile);
% 
% %Load colors:
% for ii = 1:length(excelfile)
%     a = excelfile{ii,3};
%     colors.names{ii} = excelfile{ii,1};
%     if ~isnan(a)
%         aa = strsplit(a, {',', '.'});
%         R = str2double(aa{1});
%         G = str2double(aa{2});
%         B = str2double(aa{3});
%         colors.raw(ii,:) = [R, G, B];
%     else
%         colors.raw(ii,1:3) = NaN; 
%     end
% end
% colors.rgb = colors.raw./255;
% 
% % check that no colors exceed a value of 1
% % or are less than 0
% for ii = 1:length(excelfile)
%     if any(colors.rgb(ii,:) > 1)
%         fprintf(['\n value >1 in row: ' num2str(ii)])
%     end 
%     if any(colors.rgb(ii,:) < 0)  
%         fprintf(['\n value <0 in row: ' num2str(ii)])
%     end
% end
% save('color_palette', 'colors')

end