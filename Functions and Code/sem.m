function [sem_err, std_err ] = sem(data, dim, dim_num_samples, sampleSize)
% err = sem(data, dim, dim_num_samples)
%
% PURPOSE
% calculate the standard error of the mean for a set of data
% 
% INPUTS
%   'data' : data matrix (m x n) from which to calculate SEM
%   'dim' : (optional) dimension of the matrix (if not a vector) overwhich to
%           calculate the error 
%           (default = 1 for a matrix)
%           (default length for a vector)
%   'dim_num_samples' : (optional) dimension of the data set for sample
%           size to conversion from STD to SEM 
%           (default - assumes length of the dimension ** omitting nans **)
%   'sampleSize' : (optional) number of samples for converting to SEM
%
% NOTES  
% ** don't need to provide both dim_num_samples AND sampleSize -- sample
% ** size will override the dim_num_samples
%
% OUTPUT
%   'sem_err' : vector of calculated SEM values
%   'std_err' : standard deviation of the values
%
% ES DICKINSON, 2018

%%

% Find the dimension overwhich to process the statistics
if ~exist('dim', 'var') % if there isn't a supplied dimension 
    if isvector(data)
        dim = find(size(data)>1);
    else 
        dim = 1;
    end
end

% if dimension for data sample size is provided
if exist('dim_num_samples', 'var') && ~exist('sampleSize', 'var')
    sampleSize = sum(~isnan(data),dim_num_samples);
end

% Find default number of samples if none is supplied
if ~exist('dim_num_samples', 'var') && ~exist('sampleSize', 'var') % if there isn't a supplied sample number or dimension
    if isvector(data)
        sampleSize = length(data);
    else 
        sampleSize = size(data,1);  % default dimension for number of samples is 1
    end
end

% calculate the error
std_err = std(data,0,dim,'omitnan');
sem_err = std_err ./ sqrt(sampleSize);









