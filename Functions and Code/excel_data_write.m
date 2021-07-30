
function excel_data_write(fly_ID, folder_date, total_distance_traveled)
% excel_data_write(fly_ID, folder_date, total_distance_traveled)
% Write 'Yes' in the data alignment column of Excel 'Fly Summary 2'
% Inputs:
% 'fly_ID' [string with the date and fly number, e.g. '12102018_fly_2_1]
% 'folder_date' [string with the date of the experiment, e.g. '10.14.18']
% 'total_distance_traveled' [optional addition of total path distance of fly]
% Outputs:
% adds 'Yes' to frame alignment column in the Excel file
% 
% ES Dickinson, University of Washington, Dec 2018


% Write in Excel that the data's been analyzed

[raw, Excel, xlFile] = load_flysummary;
sheet = 'Sheet1'; % excel sheet to write into       

fprintf('\n Writing to excel... \n')

flynum = fly_ID(end-2:end);
a = strcmpi(flynum, raw(:,Excel.flynum))==1; %fly number matches
b = strcmpi(folder_date, raw(:,Excel.date))==1; %fly date matches
Excel.row = find(a == 1 & b == 1); %excel row with this fly info
Excel.row = Excel.row(1);
xlRange = [Alphabet(Excel.frames) num2str(Excel.row)];

% write 'Yes' to excel file 
xlswrite(xlFile, {'Yes'}, sheet, xlRange);
if nargin == 3
    xlRange = [Alphabet(Excel.distance) num2str(Excel.row)];
    % write total distance traveled into excel file 
    xlswrite(xlFile, {num2str(total_distance_traveled)}, sheet, xlRange);
end
    
end


