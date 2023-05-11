% Num2Month conversion
% Use for converting 'date' output from MatLab (e.g.: '06-Mar-2018' ) to
% pure number (e.g. '3.6.18')
% allnumberdate = Num2Month(date)


function Folder_date = Num2Month(vid_folder)

% vid_folder = char(datetime('today'));

Month_matrix = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

Day = num2str(str2double(vid_folder(1:2)));
Month = num2str(find(strcmp(vid_folder(4:6),Month_matrix) == 1));
Year = num2str(str2double(vid_folder(10:11)));

Folder_date = [Month '.' Day '.' Year];

end
 



  
  
  
  