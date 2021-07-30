
function flyID = generate_flyID(fly_date, fly_num)
% flyID = generate_flyID(fly_date, fly_num)
% turn date and fly num into the fly ID:
% e.g. '3.28.18' & '2_0' --> '03282018_fly2_0'
% 
% 
% ES Dickinson, University of Washington, 2020


    %convert the number into fly ID date format, e.g 03282018 from 3.28.18
    c = strsplit(fly_date,'.');
    %month conversion:
    switch str2double(c(1))
        case num2cell(1:9)
            C{1} = ['0' c(1)];
        case num2cell(10:12)
            C{1} = c(1);
        otherwise
    end
    %day conversion:
    switch str2double(c(2))
        case num2cell(1:9)
            C{2} = ['0' c(2)];
        case num2cell(10:31)
            C{2} = c(2);
        otherwise
    end
    

    datestring = cell2mat([cell2mat(C{1}) cell2mat(C{2}) '20' c(3)]);
    flyID = [datestring '_fly' (fly_num)];
    
end