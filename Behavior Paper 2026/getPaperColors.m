

function [colors, formatted_names] = getPaperColors(field_names)
% [colors, formatted_names] = getPaperColors(field_names)

% field_names = {'escape jump', 'escape ring', 'food quadrant', 'fly on food', 'courtship', 'sleep'};
% kolors = {'WongOrange', 'WongRed','WongBlue','WongLightBlue', 'WongPink', 'WongGreen'}; % colors for the diff behaviors

nFields = max(size(field_names));
formatted_names = cell(size(field_names));
colors = nan(nFields, 3);
for ii = 1:nFields
    switch field_names{ii}
        case 'jump'
            cName = 'WongOrange';
            n = 'escape jump';
        case {'ring', 'OutterRing'}
            cName = 'WongRed';
            n = 'escape ring';
        case {'innerfoodquad', 'foodquad', 'innerquad', 'foodQuad'}
            cName = 'WongBlue';
            n = 'food quadrant';
        case {'FlyOnFood', 'fliesonfood', 'fly on food'}
            cName = 'WongLightBlue';
            n = 'fly on food';
        case 'CI'
            cName = 'WongPink';
            n = 'courtship';
        case 'sleep'
            cName = 'WongGreen';
            n = 'sleep';
        case 'speed'
            cName = 'WongYellow';
            n = 'speed';
    end
    formatted_names{ii} = n;
    colors(ii,:) = Color(cName);
end
        

