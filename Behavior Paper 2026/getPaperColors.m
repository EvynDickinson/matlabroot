

function colors = getPaperColors(field_names)


% field_names = {'escape jump', 'escape ring', 'food quadrant', 'fly on food', 'courtship', 'sleep'};
% kolors = {'WongOrange', 'WongRed','WongBlue','WongLightBlue', 'WongPink', 'WongGreen'}; % colors for the diff behaviors

nFields = max(size(field_names));
colors = nan(nFields, 3);
for ii = 1:nFields
    switch field_names{ii}
        case 'jump'
            cName = 'WongOrange';
        case {'ring', 'OutterRing'}
            cName = 'WongRed';
        case {'innerfoodquad', 'foodquad', 'innerquad', 'foodQuad'}
            cName = 'WongBlue';
        case {'FlyOnFood', 'fliesonfood', 'fly on food'}
            cName = 'WongLightBlue';
        case 'CI'
            cName = 'WongPink';
        case 'sleep'
            cName = 'WongGreen';
        case 'speed'
             cName = 'WongYellow';
    end

    colors(ii,:) = Color(cName);
end
        

