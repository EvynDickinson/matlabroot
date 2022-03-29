

function [kolor,num] = pullFoodColor(foodString)
% [kolor,num] = pullFoodColor(foodString)
% kolor : [R G B]
% num : 1 - plant | 2 - yeast | 3 - empty
% 

num = nan; %default
switch foodString
    case {'Yeast', 'yeast'}
        colorname = 'gold';
        num = 2;
    case {'Plant', 'plant'}
        colorname = 'green';
        num = 1;
    case {'Empty', 'empty'}
        colorname = 'grey';
        num = 3;
    % SPECIFIC VARIETIES OF PLANT FOOD:
    case 'Plant_827'
        colorname = 'palegreen';
        num = 1;
    case {'Plant_91', 'Plant9_1', 'Plant9_1C'}
        colorname = 'Darkgreen';
        num = 1;
    case 'Plant9_20_A'
        colorname = 'Chartreuse';
        num = 1;
    case {'Plant9_20_B'}
        colorname = 'teal';
        num = 1;
    case 'Plant9_20_C'
        colorname = 'cyan';
        num = 1;
    case 'Plant_11_4_21'
        colorname = 'green';
        num = 1;
    % SPECIFIC VARIETIES OF YEAST FOOD:
    case 'Yeast_11_4_21'
        colorname = 'gold';
        num = 2;
    case {'Yeast9_1','Yeast_9_1'}
        colorname = 'sandybrown';
        num = 2;
    case {'Yeast_8_27','Yeast8_27'}
        colorname = 'sienna';
        num = 2;
    case 'Yeast_9_20'
        colorname = 'maroon';
        num = 2;
    % --------------------------------
    case 'Merlot'
        colorname = 'Purple';
        num = 2;
    case 'Water'
        colorname = 'DodgerBlue';
        num = 2;
    case 'Sugar_water'
        colorname = 'Blue';
        num = 2;
    case 'Water_and_ACV'
        colorname = 'Gold';
        num = 2;
    % --------------------------------
    % TEMPERATURE RATES....
    case {-0.65, 0.65}
        colorname = 'goldenrod';
    case {-0.5, 0.5}
        colorname = 'orange';
    case {-0.30, 0.30}
        colorname = 'Purple';
    case {-0.25, 0.25}
        colorname = 'DarkViolet';    
    case {-0.16, 0.16}
        colorname = 'DodgerBlue';
    case {-0.15, 0.15}
        colorname = 'SkyBlue';
    case {-0.1,0.1}
        colorname = 'Turquoise';
    case 0
        colorname = 'white';
end

try kolor = Color(colorname);
catch 
    if ismember('Yeast', foodString)
        kolor = Color('gold');
        num = 2;
    elseif ismember('Plant', foodString)
        kolor = Color('green');
        num = 1;
    elseif ismember('Glucose', foodString)
        kolor = Color('yellow');
        num = 2;
    elseif ismember('Molasses', foodString)
        kolor = Color('purple');
        num = 2;
    else 
        kolor = Color('orangered');
        num = 2;
    end
end



