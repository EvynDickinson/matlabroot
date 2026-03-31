

function field = trimToROI(field, t_max, exp_dur)
% field = trimToROI(field, t_max, exp_dur)
%
% PURPOSE
% Trims a matrix field along its time dimension, removing all data
% beyond the region of interest. Automatically detects which dimension
% corresponds to time by matching against the full experiment duration.
%
% INPUTS:
%   'field'    : matrix to be trimmed (up to 3D: n x m, or n x m x d)
%   't_max'    : last frame index to keep (i.e. the end of the time ROI)
%   'exp_dur'  : total experiment duration in frames — used to identify
%                which dimension is the time axis
%                e.g. if exp_dur = 873000, the dimension with size 873000
%                is treated as time and trimmed to 1:t_max
%
% OUTPUTS:
%   'field'    : trimmed matrix with time dimension cut to 1:t_max
%                (unchanged if no dimension matches exp_dur)
%
% ES DICKINSON, 2026
%%
    a = size(field);
    timeDim = find(a == exp_dur, 1);
    if isempty(timeDim)
        return
    end
    switch timeDim
        case 1; field = field(1:t_max, :, :);
        case 2; field = field(:, 1:t_max, :);
        case 3; field = field(:, :, 1:t_max);
    end
end