
function kolor = pullSexColor(str)
% kolor = pullSexColor(str)
%
% Pulls the color code based on the input string describing the sex
% ES Dickinson, 2022 Yale University

switch str
    case 'Female'
        k_str = 'DeepPink';
    case 'V. Female'
        k_str = 'Pink';
    case 'V. Male'
        k_str = 'LightSkyBlue';
    case 'Male'
        k_str = 'DodgerBlue';
    case 'Mixed'
        k_str = 'white';
end

kolor = Color(k_str);







