


function color_output = get_color(color_1, color_2, N)
% 
% color_output = get_color(color_1, color_2, N)
%
% Give a color name and get out a double with 
% the three color values
% OR 
% Input two colors and an N and get a 
% matrix with a color gradient between the two 
% desired colors
%
% Input:
% 'x' [color name, string, e.g. 'Grey']
% Output:
% color_matrix = [R B G values]
%
% e.g. [1, 0.5, 0] = Color('orange')
% Find current colors in 'Color Palette' 
% on the google drive
% 
% ES Dickinson, University of Washington, Dec 2018


load('color_palette')

    

% First color requested:
color_loc = strcmpi(color_1, colors.names);
if sum(color_loc) == 0
   warndlg(['Cannot find ''' color_1 ''' in the color palette'])
   return 
end
color_matrix(1,:) = colors.rgb(color_loc, :);
color_output(1,:) = color_matrix;

% 2+ colors requested
switch nargin
    case 2
        % Two colors requested:
        color_loc = strcmpi(color_2, colors.names);
        if sum(color_loc) == 0
           warndlg(['Cannot find ''' color_2 ''' in the color palette'])
           return 
        end
        color_matrix(2,:) = colors.rgb(color_loc, :);
        color_output = color_matrix;
    case 3
        % Color gradient requested:
        color_loc = strcmpi(color_2, colors.names);
        if sum(color_loc) == 0
           warndlg(['Cannot find ''' color_2 ''' in the color palette'])
           return 
        end
        color_matrix(2,:) = colors.rgb(color_loc, :);
        % Create the gradient
        interval = (color_matrix(1,:)-color_matrix(2,:))/(N-1);
        for ii = 1:N-1 %RGB   
            color_output(ii+1,:) = color_output(ii,:)-interval;
        end  
        % Check for any values > 1 || < 0 
        for ii = 1:size(color_output,1)
            filterup = color_output(ii,:) > 1;
            filterdown = color_output(ii,:) < 0;
            if sum(filterup) > 0
                color_output(ii,filterup) = 1;
            end
            if sum(filterdown) > 0
                color_output(ii,filterdown) = 0;
            end
        end
            
end
    



%% HOW THE COLORS ARE UPDATED|CREATED:
% xlFile = 'G:\My Drive\Data\FicTrac Raw Data\Color Palette.xlsx';      
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
%     if colors.rgb(ii,rr) > 1  
%         fprintf(['\n value >1 in row: ' num2str(ii)])
%     end 
%     if colors.rgb(ii,rr) < 0  
%         fprintf(['\n value <0 in row: ' num2str(ii)])
%     end
% end
% save('color_palette', 'colors')


%% Visual Testing Example
% color_1 = 'Red';
% color_2 = 'Blue';
% N = 100; 
% 
% hfig = figure; set(hfig, 'pos', [50, 50, 800, 900], 'Color', 'w'); hold all
%     % case 1 one color
%     color_output = get_color(color_1);
%     subplot(3,1,1)
%         plot(1:5, 1:5, 'color', color_output, 'LineWidth', 4)
%         title('Color One')
%     % case 2 2 colors
%     color_output = get_color(color_1, color_2);
%     subplot(3,1,2)
%         hold all
%         plot(1:5, 1:5, 'color', color_output(1,:), 'LineWidth', 4)
%         plot(1:5, 1.5:5.5, 'color', color_output(2,:), 'LineWidth', 4)
%         title('Color One & Two')
%     % case 3 gradient of colors
%     color_output = get_color(color_1, color_2, N);
%     subplot(3,1,3)
%         hold all
%         for ii = 1:N
%             scatter(ii, ii, 50, color_output(ii,:), 'filled')
%         end
%         title('Color Gradient from Color One to Color Two')



end