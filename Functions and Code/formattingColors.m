


function [foreColor,backColor] = formattingColors(blkbgd)
% [foreColor,backColor] = formattingColors(blkbgd)
% True = black background figure
% False = white background figure
% this just gives the colors for the foreground color and background color


if blkbgd==true
    foreColor = 'w';
    backColor = 'k';
else
    foreColor = 'k';
    backColor = 'w';
end