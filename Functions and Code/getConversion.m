
function [conversion, con_type] = getConversion(exp_date, plate, exp_type)
% [conversion, con_type] = getConversion(exp_date, plate)
% ex: exp_date = 'MM.DD.YYYY'
% plate is 1 or 2 for old vs new plate
% exp_type: 1 = low res, 2 = high res
% con_type will be the number(s) to use for the conversion -- but gives 2 if you
% don't give the number of the plate
% 
% pull the numbers and date standards for the experimental plates
% and their sizes etc 

date_switch = datetime('10.20.2023', 'InputFormat', 'MM.dd.yyyy');

i = 1; 
conversion(i).name = 'Cedar st plate 1';
conversion(i).cutoff = 739200; % this is the date that BEFORE uses this plate and size for plate 1
conversion(i).date_switch = date_switch;
conversion(i).pix2mm = 13.09; % conversion from pixels to mm
conversion(i).R_real = 35.2; % mm distance of the LEGIT outer limit (including inaccessible space)
conversion(i).R = 31.5; % mm distance that is accessible to the flies

i = 2; 
conversion(i).name = 'College st plate 1';
conversion(i).cutoff = 739200; % this is the date that AFTER uses this plate and size for plate 1
conversion(i).date_switch = date_switch;
conversion(i).pix2mm = 12.97; % conversion from pixels to mm
conversion(i).R_real = 35.2; % mm distance of the LEGIT outer limit (including inaccessible space)
conversion(i).R = 31.5; % mm distance that is accessible to the flies

i = 3; 
conversion(i).name = 'College st plate 2';
conversion(i).cutoff = 739200; % all trials AFTER this date use this size for plate 2
conversion(i).date_switch = date_switch;
conversion(i).pix2mm = 13.00; % conversion from pixels to mm
conversion(i).R_real = 33.7; % mm distance of the LEGIT outer limit (including inaccessible space)
conversion(i).R = 30; % mm distance that is accessible to the flies

i = 4; 
conversion(i).name = 'College st plate 2 high res';
conversion(i).cutoff = 739200; % all trials AFTER this date use this size for plate 2
conversion(i).date_switch = date_switch;
conversion(i).pix2mm = 30.00; % conversion from pixels to mm
conversion(i).R_real = 33.7; % mm distance of the LEGIT outer limit (including inaccessible space)
conversion(i).R = 30; % mm distance that is accessible to the flies

% arena sizes / distances 
for i = 1:4
    conversion(i).circle75 = conversion(i).R*sqrt(3/4); % radius of circle that occupies 75% of the arena (aka outer 25%)
    conversion(i).circle10 = conversion(i).R*sqrt(0.1); % radius of a circle occupying 10% of the arena
    conversion(i).circle7 = conversion(i).R*sqrt(0.07); % radius of circle that occupies 75% of the arena (aka outer 25%)
    conversion(i).circle5 = conversion(i).R*sqrt(0.05); % radius of a circle occupying 10% of the arena
end


% find the conversion type if needed: 
% (exp_date, plate, exp_type)
if nargout > 1 
   if exist('exp_date','var')
       date_opt = ischar(exp_date); % is there a real date?
       testDate = datetime(exp_date,'InputFormat', 'MM.dd.yyyy') ; % this is the date time for the current experiment
   else
       date_opt = false;
   end
   if exist('plate','var') 
       plate_opt = any(plate==[1,2]); % is the plate a real choice? 
   end
   if nargin == 3 && any(exp_type==[1,2])% is there an experiment type?
       exp_opt = true;
   else
       exp_opt = false;
   end

    con_opts = true(1,4); % start options all true
   
    % ----------- Rule OUT option 1 ------------
    % if we the date is after the building switch
    if date_opt && (testDate > date_switch) 
        con_opts(1) = false;
    end
    % if trial is high res
    if exp_opt && exp_type==2
        con_opts(1) = false;
    end
    % if the plate is not plate 1
    if plate_opt && plate==2
        con_opts(1) = false;
    end

    % ----------- Rule OUT option 2 ------------
    % if the plate is not plate 1
    if plate_opt && plate==2
        con_opts(2) = false;
    end
    % if the date is before the time switch
    if date_opt && (testDate < date_switch) 
        con_opts(2) = false;
    end 
    % if trial is high res
    if exp_opt && exp_type==2
        con_opts(2) = false;
    end

    % ----------- Rule OUT option 3 ------------
    % if the plate is not plate 2
    if plate_opt && plate==1
        con_opts(3) = false;
    end
    % if the date is before the time switch
    if date_opt && (testDate < date_switch) 
        con_opts(3) = false;
    end 
    % if trial is high res
    if exp_opt && exp_type==2
        con_opts(3) = false;
    end

    % ----------- Rule OUT option 4 ------------
    % if the plate is not plate 2
    if plate_opt && plate==1
        con_opts(4) = false;
    end
    % if the date is before the time switch
    if date_opt && (testDate < date_switch) 
        con_opts(4) = false;
    end 
    % if trial is low res
    if exp_opt && exp_type==1
        con_opts(4) = false;
    end
    
    con_type = find(con_opts);
end
    


%     % no plate specification:
%     if ~plate_opt 
%         if testDate < date_switch
%             con_type = 1; % this is the only plate that took place before the building switch date
%         else
%             con_type = 2:3; 
%         end
% 
%     % plate specification
%     else
%             switch plate
%                     case 1 % plate one (old/OG plate) 
%                         if testDate < date_switch
%                             con_type = 1;
%                         else con_type = 2;
%                         end
%                     case 2 % plate two (new plate, smaller) 
%                         con_type = 3; % technically this can also be 4 if it is high resolution...
%             end
%     end
% 
% end












