
function [peak_locations, peak_mag] = findPeaks(varargin)
% peak_locations = findPeaks(varargin)
% Use this to adjust and pass data to the peakfinder.m function in the proper
% format. Argument pairs should be lead with '-'.
% ES Dickinson, UW, 2020
%
% ex: max_loc = findPeaks('-x', raw, '-extrema', 1); % find max peaks
% 
%INPUTS:
%       -x - A real vector from the maxima will be found (required)
%       -sel - The amount above surrounding data for a peak to be,
%           identified (default = (max(x0)-min(x0))/4). Larger values mean
%           the algorithm is more selective in finding peaks.
%       -thresh - A threshold value which peaks must be larger than to be
%           maxima or smaller than to be minima.
%       -extrema - 1 if maxima are desired, -1 if minima are desired
%           (default = maxima, 1)
%       -includeEndpoints - If true the endpoints will be included as
%           possible extrema otherwise they will not be included
%           (default = true)
%       -interpolate - If true quadratic interpolation will be performed
%           around each extrema to estimate the magnitude and the
%           position of the peak in terms of fractional indicies. Note that
%           unlike the rest of this function interpolation assumes the
%           input is equally spaced. To recover the x_values of the input
%           rather than the fractional indicies you can do:
%           peakX = x0 + (peakLoc - 1) * dx
%           where x0 is the first x value and dx is the spacing of the
%           vector. Output peakMag to recover interpolated magnitudes.
%           See example 2 for more information.
%           (default = false)
%
%   OUTPUTS:
%       peakLoc - The indicies of the identified peaks in x0
%       peakMag - The magnitude of the identified peaks
%
%   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
%       are at least 1/4 the range of the data above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
%       that are at least sel above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local
%       maxima that are at least sel above surrounding data and larger
%       (smaller) than thresh if you are finding maxima (minima).
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
%       data if extrema > 0 and the minima of the data if extrema < 0
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema, includeEndpoints)
%       returns the endpoints as possible extrema if includeEndpoints is
%       considered true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,sel,thresh,extrema,interpolate)
%       returns the results of results of quadratic interpolate around each
%       extrema if interpolate is considered to be true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
%       local maxima as well as the magnitudes of those maxima
%


% default paths & settings: 
extrema = 1; % find maxima
thresh = []; % no manual threshold
sel = []; % no relative selection height
includeEndpoints = true;
Intrplt = true;


% parse the input parameters:
for ii = 1:nargin
    if ischar(varargin{ii}) && ~isempty(varargin{ii})
        if varargin{ii}(1) == '-' %find command descriptions
            switch lower(varargin{ii}(2:end))
                case 'x'
                    raw_data = varargin{ii+1};
                case 'sel'
                    sel = varargin{ii+1};
                case 'thresh'
                    thresh = varargin{ii+1};
                case 'extrema'
                    extrema = varargin{ii+1};
                case 'includeendpoints'
                    includeEndpoints = varargin{ii+1};
                case 'interpolate'
                    Intrplt = varargin{ii+1};
            end
        end
    end
end

% adjust the data to exist only in the positive or negative realm depending
% on the extrema selection
switch extrema
    case 1 %find maxima
        a = min(raw_data);
        if a<=0
            x = raw_data+abs(a);
            offset = abs(a);
        else
            x = raw_data;
            offset = 0;
        end
    case -1 %find minima
        a = max(raw_data);
        if a>=0
            x = raw_data-abs(a);
            offset = abs(a)*-1;
        else
            x = raw_data;
            offset = 0;
        end
end

% input the data in the actual function: peakfinder.m

[peak_locations, mag] = peakfinder(x,sel,thresh, extrema, includeEndpoints, Intrplt);

% adjust the peak_mag for the initial offset:

peak_mag = mag+offset;





end



