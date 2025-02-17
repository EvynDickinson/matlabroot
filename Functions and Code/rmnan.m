function y = rmnan(y)
% y = rmnan(y)
% remove nans from the vector y 
% shortens the vector to the length of y without nans
%
% ES Dickinson

loc = isnan(y);
y(loc) = [];


