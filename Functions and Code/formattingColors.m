


function [foreColor,backColor,forestr, backstr] = formattingColors(blkbgd,textON)
% [foreColor,backColor] = formattingColors(blkbgd,textON)
% True = black background figure
% False = white background figure
% textON is true if you want the value in string form: 'w' | 'k'
% this just gives the colors for the foreground color and background color


if blkbgd==true
    foreColor =  [1 1 1]; %white
    backColor =  [0 0 0]; %black
    if nargin>1 && textON
        forestr = 'w';
        backstr = 'k';
    end
else
    foreColor = [0 0 0]; %black
    backColor = [1 1 1]; %white
     if nargin>1 && textON
        forestr = 'k';
        backstr = 'w';
    end
end


% 
% if blkbgd==true
%     foreColor = 'w';
%     backColor = 'k';
% else
%     foreColor = 'k';
%     backColor = 'w';
% end