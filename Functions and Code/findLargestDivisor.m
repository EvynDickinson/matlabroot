

function largestDivisor = findLargestDivisor(startRange, endRange,maxDivisorSize)
% largestDivisor = findLargestDivisor(startRange, endRange,maxDivisorSize)

    % Calculate the length of the range
    rangeLength = endRange - startRange + 1;
    
    % Calculate the largest divisor that would result in fewer than the target groups
    % We use floor(rangeLength / maxDivisorSize) because we want to ensure
    % we have fewer than the target number
    largestDivisor = floor(rangeLength / (maxDivisorSize));
    
    % Ensure the divisor is at least 1
    if largestDivisor < 1
        largestDivisor = 1;
    end
end
