
function x_D = computeDerivative(x, t, dt)
% x_D = computeDerivative(x, t, dt)
%
% PURPOSE
% Estimates the instantaneous derivative of a 1D signal x at each sample
% using the closed-form slope of a least-squares linear fit over a trailing
% window of t samples.
%
% INPUTS
%   'x'  : input signal
%       n-by-1 vector of evenly sampled data
%   't'  : window size in samples
%       integer, must be >= 2. Larger values give smoother derivative
%       estimates at the cost of temporal resolution.
%   'dt' : sample period in seconds
%       scalar, e.g. 1/30 for a 30 Hz signal
%
% OUTPUTS
%   'x_D' : sliding derivative of x
%       n-by-1 vector, same length as x. First t-1 samples are NaN
%       (window not yet full). Units are x per second.
%
% NOTES
%   The derivative is computed as the slope of a linear fit over each
%   trailing window using a closed-form weighted sum (equivalent to
%   polyfit slope but without the overhead). Implemented via convolution
%   for efficiency — O(n log n) regardless of window size t. More robust
%   to noise than a simple two-point finite difference.
%
% EXAMPLE
%   x_D = computeDerivative(x, 150, 1/30)    % 5-second window at 30 Hz
%
% ES DICKINSON, YALE, 2026

%%
ramp   = (-(t-1)/2 : (t-1)/2)';              % centred linear ramp over window, used as slope weights
w      = ramp / (sum(ramp.^2) * dt);          % normalise weights so output is in units of x per second
x_conv = conv(x, flipud(w), 'full');          % sliding weighted sum — equivalent to slope of linear fit at each sample
n      = length(x);                           % number of input samples
x_D    = [nan(t-1, 1); x_conv(t:t+n-t)];     % pad head with NaN (window not yet full), trim to input length
