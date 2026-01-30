

function data = empty2nan(data)
% data = empty2nan(data)
%
% Replace an empty variable with a nan
% if the variable is empty, replace it with a single nan (for the entire variable)
% 
% INPUT
%   'data' : variable of unknown type, to be determined if it is empty
%            if data variable is empty, it will be assigned nan
% 
% OUTPUT
%   'data' : variable as either a nan or the original (nonempty) variable value
%
% EXAMPLE
% e.g., data = [5 6]; data out will be [5 6]
% e.g, data = 0 x 0 empty; data out will be nan
%
% CODE 
% if isempty(data)
%     data = nan;
% end
%
% ES DICKINSON,  2026

%% 
if isempty(data)
    data = nan;
end
