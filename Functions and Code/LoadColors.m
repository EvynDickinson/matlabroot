
function [baseColor, foreColor] = LoadColors(BandWoption, ~)
% [baseColor, foreColor] = LoadColors(BandWoption, RGBflag)
% set background colors based on the image being a black or white
% background setting. e.g. white foreground colors for a black background
% uses a logical input for true (black background) or false (white
% background). 
% RGBflag : present  = colors come as [R B G] format
% RGBflag : absent = colors come as string format 'white'
%
%
%
% ES Dickinson,
% 2021 Yale University

if nargin == 2
    if BandWoption == true
        baseColor = [0 0 0];
        foreColor = [1 1 1];
    else 
        baseColor = [1 1 1];
        foreColor = [0 0 0];
    end
else
    if BandWoption == true
        baseColor = 'black';
        foreColor = 'white';
    else 
        baseColor = 'white';
        foreColor = 'black';
    end
end
