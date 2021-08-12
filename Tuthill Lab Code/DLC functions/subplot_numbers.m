
function [nrows, ncols] = subplot_numbers(num_instances, maxCol)
% [nrows, ncols] = subplot_numbers(num_instances, maxCol)
% Get numbers for the rows and columns of a subplot that are 
% most square
% input: 
% num_instances -- how many plots are desired
% maxCol -- max number of columns
% 
% 
% ES Dickinson
% University of Washington, 2020

if nargin == 1
    maxCol = 4;
end

aa = divisors(num_instances);
a = sqrt(num_instances);
if length(aa) == 2 %no divisors
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
%     a = max([dim1, dim2]);
    dim1 = maxCol;
    dim2 = ceil(num_instances/dim1);
end
nrows = max([dim1, dim2]);
ncols = min([dim1, dim2]);

end