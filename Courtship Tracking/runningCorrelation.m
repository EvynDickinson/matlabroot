
function running_corr = runningCorrelation(data, windowSize)
    % Input:
    % - data: a matrix with two columns representing the time series
    % - windowSize: the size of the sliding window (number of points)

    % Number of data points
    n = size(data, 1);
    
    % Pre-allocate running correlation array
    running_corr = NaN(n - windowSize + 1, 1);  % Result will be shorter than original data due to window

    % Loop through the data with the specified window size
    for i = 1:(n - windowSize + 1)
        % Get the current window of data
        windowData = data(i : i + windowSize - 1, :);
        
        % Calculate correlation between the two columns in the current window
        R = corrcoef(windowData(:, 1), windowData(:, 2));
        
        % Store the correlation coefficient (off-diagonal element in R)
        running_corr(i) = R(1, 2);
    end
end