

%% Run PCA on behavioral data:

% 1) Pre-process and organize the data

% select experiment
exp =1;

raw_data = [];
for trial = 1:num.trial(exp)

    % pull the data components that will be included in the analysis
    temperature = data(exp).data(trial).data.T.temperature;
    temp_rate = [data(exp).data(trial).data.TR.data(:,4);nan];
    x_loc = data(exp).data(trial).data.x_loc;
    y_loc = data(exp).data(trial).data.y_loc;
    
    %combined all components into a matrix
    for ii = 1:size(x_loc,2)
        raw_data = [raw_data; temperature, temp_rate, x_loc(:,ii), y_loc(:,ii)];
    end
end
% 2) Standardize the data

% remove locations without a temperature rate
loc = any(isnan(raw_data),2);
raw_data(loc,:) = [];

% normalize the data
std_data = zscore(raw_data);

% 3) Perform PCA

[coeff, score, latent, tsquared, explained, mu] = pca(std_data);


% 4) Analyze the results

cumulative_explained = cumsum(explained);

figure;
scatter3(score(:,1), score(:,2),score(:,3));
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('PCA of Time Series Data');
