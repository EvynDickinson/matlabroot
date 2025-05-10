
function conversion = getConversion
% conversion = getConversion
% pull the numbers and date standards for the experimental plates
% and their sizes etc 

i = 1; 
conversion(i).name = 'Cedar st plate 1';
conversion(i).cutoff = 739200; % this is the date that BEFORE uses this plate and size for plate 1
conversion(i).pix2mm = 13.0944521751958; % conversion from pixels to mm
conversion(i).R_real = 35.2; % mm distance of the LEGIT outer limit (including inaccessible space)
conversion(i).R = 31.5; % mm distance that is accessible to the flies
conversion(i).circle75 = conversion(i).R*sqrt(3/4); % radius of circle that occupies 75% of the arena (aka outer 25%)
conversion(i).circle10 = conversion(i).R*sqrt(0.1); % radius of a circle occupying 10% of the arena


i = 2; 
conversion(i).name = 'College st plate 1';
conversion(i).cutoff = 739200; % this is the date that AFTER uses this plate and size for plate 1
conversion(i).pix2mm = 12.9709545178985; % conversion from pixels to mm
conversion(i).R_real = 35.2; % mm distance of the LEGIT outer limit (including inaccessible space)
conversion(i).R = 31.5; % mm distance that is accessible to the flies
conversion(i).circle75 = conversion(i).R*sqrt(3/4); % radius of circle that occupies 75% of the arena (aka outer 25%)
conversion(i).circle10 = conversion(i).R*sqrt(0.1); % radius of a circle occupying 10% of the arena

i = 3; 
conversion(i).name = 'College st plate 2';
conversion(i).cutoff = 739200; % this is the date that AFTER uses this plate and size for plate 1
conversion(i).pix2mm = 13.0032195776504; % conversion from pixels to mm
conversion(i).R_real = 33.7; % mm distance of the LEGIT outer limit (including inaccessible space)
conversion(i).R = 30; % mm distance that is accessible to the flies
conversion(i).circle75 = conversion(i).R*sqrt(3/4); % radius of circle that occupies 75% of the arena (aka outer 25%)
conversion(i).circle10 = conversion(i).R*sqrt(0.1); % radius of a circle occupying 10% of the arena




