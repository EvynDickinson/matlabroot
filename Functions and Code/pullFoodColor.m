

function kolor = pullFoodColor(foodString)
% kolor = pullFoodColor(foodString)

switch foodString
    case 'Yeast'
        colorname = 'gold';
    case 'Plant'
        colorname = 'green';
    case 'Empty'
        colorname = 'grey';
    % SPECIFIC VARIETIES OF PLANT FOOD:
    case 'Plant_827'
        colorname = 'palegreen';
    case {'Plant_91', 'Plant9_1'}
        colorname = 'Darkgreen';
    case 'Plant9_20_A'
        colorname = 'Chartreuse';
    case 'Plant9_20_B'
        colorname = 'teal';
    case 'Plant9_20_C'
        colorname = 'cyan';
    % SPECIFIC VARIETIES OF YEAST FOOD:
    case {'Yeast9_1','Yeast_9_1'}
        colorname = 'sandybrown';
    case {'Yeast_8_27','Yeast8_27'}
        colorname = 'sienna';
    case 'Yeast_9_20'
        colorname = 'maroon';
end

kolor = Color(colorname);
end