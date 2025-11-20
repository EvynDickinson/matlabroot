

function data = replaceNaN(data, newValue)
% data = replaceNaN(data, newValue)
% replace all nans in the structure with a new set value
% ES Dickinson 2025 Yale University


if isstruct(data)
    fields = fieldnames(data);
    for i = 1:numel(fields)
        if isnumeric(data.(fields{i}))
            data.(fields{i})(isnan(data.(fields{i}))) = newValue;
        end
    end
else
    data(isnan(data)) = newValue;
end
