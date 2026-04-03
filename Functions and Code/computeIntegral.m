
function x_I = computeIntegral(x, t, dt)
% x_I = computeIntegral(x, t, dt)
%
% PURPOSE
% Computes a sliding window trapezoidal integral of a 1D signal x over a
% trailing window of t samples, returning a vector of the same length as x.
%
% INPUTS
%   'x'  : input signal
%       n-by-1 vector of evenly sampled data
%   't'  : window size in samples
%       integer, must be >= 2 for trapezoidal weighting to apply
%   'dt' : sample period in seconds
%       scalar, e.g. 1/30 for a 30 Hz signal
%
% OUTPUTS
%   'x_I' : sliding integral of x
%       n-by-1 vector, same length as x. First t-1 samples are NaN
%       (window not yet full). Units are x * seconds.
%
% NOTES
%   Endpoint weighting follows the trapezoidal rule (half-weight on first
%   and last sample of each window). Implemented via convolution for
%   efficiency — O(n log n) regardless of window size t.
%
% EXAMPLE
%   x_I = computeIntegral(x, 150, 1/30)    % 5-second window at 30 Hz
%
% ES DICKINSON, YALE, 2026


%%
kernel          = ones(t, 1);                        % boxcar kernel of length t for sliding sum
kernel([1 end]) = 0.5;                               % trapz correction: half-weight the endpoints
x_conv          = conv(x, kernel, 'full') * dt;      % sliding trapezoidal integral, scaled by sample period
n               = length(x);                         % number of input samples
x_I             = [nan(t-1, 1); x_conv(t:t+n-t)];   % pad head with NaN (window not yet full), trim to input length



