
function [nrows, ncols] = subplot_numbers(num_instances, maxCol)
% [nrows, ncols] = subplot_numbers(num_instances, maxCol)
%
% PURPOSE
% Get column and row dimensions for the desired number of 
% subplots that are the most square
% 
% INPUTS
%   'num_instances' : total number of subplots
%   'maxCol' : max number of columns
% 
% OUTPUTS
%   'nrows' : number of rows
%   'ncols' : number of columns
% 
% ES DICKINSON, 2020

%%
if nargin == 1
    maxCol = 4;
end

aa = divisors(num_instances);
a = sqrt(num_instances);
if length(aa) == 2 % no divisors
    aa = [2, 3, 4, 5];
    [~,idx] = min(abs(aa-a));
    dim1 = aa(idx);
    dim2 = ceil(num_instances/dim1);
else
    [~,idx] = min(abs(aa-a));
    dim1 = aa(idx);
    dim2 = num_instances/dim1;
end
if dim1 > maxCol || dim2 > maxCol
    dim1 = maxCol;
    dim2 = ceil(num_instances/dim1);
end

nrows = max([dim1, dim2]);
ncols = min([dim1, dim2]);


