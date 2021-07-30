function err = sem(data, dim, dim_num_samples)
% err = sem(data, dim, dim_num_samples)
%conversion from standard deviation to standard error of the mean
%
%
%
% ED Dickinson
% University of Washington, 2018

if nargin == 1
    dim = 1;
    dim_num_samples = dim;
elseif nargin == 2
    dim_num_samples = dim;
end 



switch dim_num_samples
    case 2
        sampleSize = sum(~isnan(data(1,:)));
    case 1
        sampleSize = sum(~isnan(data(:,1)));
end

st_err = nanstd(data,0,dim);
err = st_err / sqrt(sampleSize);



% size(data,dim_num_samples)
end