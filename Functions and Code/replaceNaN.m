

function data = replaceNaN(data, newValue)
% data = replaceNaN(data, newValue)
%
% PURPOSE
% replace all nans in the structure with a new set value
%
% INPUTS
%    'data' : input matrix OR structure! this will iterate through all fields
%           and if they are numeric, will run this function on them
%    'newValue' : value to replaces the NaN locations (e.g. 'false)
%
% OUTPUTS
%     'data' : data matrix or structure with updated values
%
% ES DICKINSON, 2025

%%
if isstruct(data)
    fields = fieldnames(data);
    for i = 1:numel(fields)
        if isnumeric(data.(fields{i}))
            loc = isnan(data.(fields{i})); % location of nan values in the matrix
            data.(fields{i})(loc) = newValue;
        end
    end
else
    data(isnan(data)) = newValue;
end
