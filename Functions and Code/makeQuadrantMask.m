

% Define square image dimensions
imageSize = 256;  % Size of the square image
radius = imageSize / 4;  % Radius of the circle (adjust as needed)
center = [imageSize / 2, imageSize / 2];  % Center of the circle

% Create a logical mask for the circle
[X, Y] = meshgrid(1:imageSize, 1:imageSize);
circleMask = (X - center(1)).^2 + (Y - center(2)).^2 <= radius^2;

% Divide the circle into 4 quadrants
quadrant1Mask = circleMask & X >= center(1) & Y < center(2);  % Top-right
quadrant2Mask = circleMask & X < center(1) & Y < center(2);   % Top-left
quadrant3Mask = circleMask & X < center(1) & Y >= center(2);  % Bottom-left
quadrant4Mask = circleMask & X >= center(1) & Y >= center(2); % Bottom-right

% Calculate the ROI X and Y positions for each quadrant
[Q1_Y, Q1_X] = find(quadrant1Mask);
[Q2_Y, Q2_X] = find(quadrant2Mask);
[Q3_Y, Q3_X] = find(quadrant3Mask);
[Q4_Y, Q4_X] = find(quadrant4Mask);

% Plot the circle and quadrants for visualization
figure;
imshow(circleMask); hold on;
plot(Q1_X, Q1_Y, 'r.', 'DisplayName', 'Quadrant 1 (Top-right)');
plot(Q2_X, Q2_Y, 'g.', 'DisplayName', 'Quadrant 2 (Top-left)');
plot(Q3_X, Q3_Y, 'b.', 'DisplayName', 'Quadrant 3 (Bottom-left)');
plot(Q4_X, Q4_Y, 'm.', 'DisplayName', 'Quadrant 4 (Bottom-right)');
legend('show');
title('Circle divided into Quadrants');