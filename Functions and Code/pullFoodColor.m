

function kolor = pullFoodColor(foodString)
% kolor = pullFoodColor(foodString)

switch foodString
    case 'Yeast'
        kolor = Color('gold');
    case 'Plant'
        kolor = Color('green');
    case 'Empty'
        kolor = Color('grey');
    case 'Plant_827'
        kolor = Color('palegreen');
    case 'Plant_91'
        kolor = Color('Darkgreen');
end

end