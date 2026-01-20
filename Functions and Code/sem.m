function err = sem(data, dim, dim_num_samples)
% err = sem(data, dim, dim_num_samples)
% calculate the standard error of the mean for a set of data
% 
% INPUTS
%   'data' : data matrix (m x n) from which to calculate SEM
%   'dim' : (optional) dimension of the matrix (if not a vector) overwhich to
%           calculate the error 
%           (default - 1)
%   'dim_num_samples' : (optional) number of samples within the data set for
%           conversion to SEM from STD
%           (default - assumes length of the dimension omitting nans)
%
% OUTPUT
%   'err' : vector of calculated SEM values
%
% ES DICKINSON, 2018

%%
% get defaults
if nargin == 1
    dim = 1;
    dim_num_samples = dim;
elseif nargin == 2
    dim_num_samples = dim;
end 

% get sample size default values
switch dim_num_samples
    case 2
        sampleSize = sum(~isnan(data(1,:)));
    case 1
        sampleSize = sum(~isnan(data(:,1)));
end

% calculate the error
st_err = std(data,0,dim,'omitnan');
err = st_err / sqrt(sampleSize);









