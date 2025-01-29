

function data = empty2nan(data)
% data = empty2nan(data)
% if the variable is empty, replace it with a nan (for the entire variable)
% e.g., data = [5 6]; data out will be [5 6]
% e.g, data = 0 x 0 empty; data out will be nan
%
% code: 
% if isempty(data)
%     data = nan;
% end

if isempty(data)
    data = nan;
end
