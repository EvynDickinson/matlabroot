
% Define square image dimensions
imageSize = 256;  % Size of the square image
radius = imageSize / 4;  % Radius of the circle
center = [imageSize / 2, imageSize / 2];  % Center of the circle

% Create a logical mask for the circle
[X, Y] = meshgrid(1:imageSize, 1:imageSize);
circleMask = (X - center(1)).^2 + (Y - center(2)).^2 <= radius^2;

% Define the diagonal quadrant boundaries (forming an "X")
% Quadrant 1: Top-left (Y > X)
% Quadrant 2: Top-right (Y > -X + 2 * center(1))
% Quadrant 3: Bottom-right (Y <= X)
% Quadrant 4: Bottom-left (Y <= -X + 2 * center(1))
quadrant1Mask = circleMask & (Y > X) & (Y > -X + 2 * center(1));   % Top-left
quadrant2Mask = circleMask & (Y > X) & (Y <= -X + 2 * center(1));  % Top-right
quadrant3Mask = circleMask & (Y <= X) & (Y <= -X + 2 * center(1)); % Bottom-right
quadrant4Mask = circleMask & (Y <= X) & (Y > -X + 2 * center(1));  % Bottom-left

% Calculate the ROI X and Y positions for each quadrant
[Q1_Y, Q1_X] = find(quadrant1Mask);
[Q2_Y, Q2_X] = find(quadrant2Mask);
[Q3_Y, Q3_X] = find(quadrant3Mask);
[Q4_Y, Q4_X] = find(quadrant4Mask);

% Plot the circle and quadrants for visualization
figure;
imshow(circleMask); hold on;
plot(Q1_X, Q1_Y, 'r.', 'DisplayName', 'Quadrant 1 (Top-left)');
plot(Q2_X, Q2_Y, 'g.', 'DisplayName', 'Quadrant 2 (Top-right)');
plot(Q3_X, Q3_Y, 'b.', 'DisplayName', 'Quadrant 3 (Bottom-right)');
plot(Q4_X, Q4_Y, 'm.', 'DisplayName', 'Quadrant 4 (Bottom-left)');
legend('show');
title('Circle Divided into Quadrants Forming an X');


%% Overlay and test the quadrant mask:
clearvars('-except',initial_vars{:})

for exp = 1:num.exp
    quad_occ = [];
    for trial = 1:num.trial(exp)
        center = data(exp).data(trial).data.centre;
        x_loc = data(exp).data(trial).data.x_loc;
        y_loc = data(exp).data(trial).data.y_loc;
        r = data(exp).data(trial).data.r;
        foodWell = data(exp).T.foodLoc(trial);

        % Adjust the X and Y coordinates relative to the new center
        adjustedX = x_loc - center(1);
        adjustedY = y_loc - center(2);
        
        % Initialize matrix to hold quadrant classification (same size as input matrices)
        quadrantMatrix = zeros(size(x_loc));
        
        % Define quadrant masks based on the new center
        Q = [];
        Q(1).Mask = (adjustedY > adjustedX) & (adjustedY <= -adjustedX);  % Top
        Q(2).Mask = (adjustedY <= adjustedX) & (adjustedY <= -adjustedX); % Bottom
        Q(3).Mask = (adjustedY <= adjustedX) & (adjustedY > -adjustedX);  % Left
        Q(4).Mask = (adjustedY > adjustedX) & (adjustedY > -adjustedX);   % Right
        
        % Determine the well locations (which determine quadrant assignment):
        adjusted_wx = data(exp).data(trial).data.wellcenters(1,foodWell) - center(1);
        adjusted_wy = data(exp).data(trial).data.wellcenters(2,foodWell) - center(2);
        
        idx_loc = false(1,4);
        % Find the food quadrant (find location with the food well coordinates included)
        idx_loc(1) = (adjusted_wy > adjusted_wx) & (adjusted_wy <= -adjusted_wx);  % top
        idx_loc(2) = (adjusted_wy <= adjusted_wx) & (adjusted_wy <= -adjusted_wx); % right
        idx_loc(3) = (adjusted_wy <= adjusted_wx) & (adjusted_wy > -adjusted_wx);  % bottom
        idx_loc(4) = (adjusted_wy > adjusted_wx) & (adjusted_wy > -adjusted_wx);   % left
        quad_loc = find(idx_loc);

        fly_loc = ~isnan(x_loc); %gives logical for all fly positions in the position matrix
        foodQuad = Q(quad_loc).Mask & fly_loc; % flies in food quad
        nflies = data(exp).T.NumFlies(trial);
        y = (sum(foodQuad,2)./nflies).*100;
        quad_occ = autoCat(quad_occ, y,false);

        %  A = quad_occ;
        %  B = data(exp).data(trial).data.occupancy.dist2wells(:,foodWell);
        %  C = data(exp).data(trial).data.occupancy.occ(:,foodWell)*100;
        % 
        % figure; hold on
        % plot(A,'color', Color('red'));
        % plot(C,'color', 'k')
        % yyaxis right 
        % plot(B, 'color', 'y')
        % R = corrcoef(A,B);
        % disp(['Distance: ' num2str(R(2))])
        %  R = corrcoef(A,C);
        % disp(['ROI: ' num2str(R(2))])
    
    end

    % save the data to a larger structure
    grouped(exp).quadrant.all = quad_occ;
    grouped(exp).quadrant.avg = mean(quad_occ,2,'omitnan');
    grouped(exp).quadrant.std = std(quad_occ,0,2,'omitnan');
end
          
% Pool the data for heating and cooling together: 
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
    all_quad = [];
    
    % Update the averages for the classic temperature bins 
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        end
        raw_c(t,:) = mean(grouped(exp).quadrant.all(c_frames,:),1,'omitnan');
        raw_h(t,:) = mean(grouped(exp).quadrant.all(h_frames,:),1,'omitnan');
    end

    % find the avg and err and save to group structure
    grouped(exp).quadrant.increasing.raw = raw_h;
    grouped(exp).quadrant.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).quadrant.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).quadrant.decreasing.raw = raw_c;
    grouped(exp).quadrant.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).quadrant.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).quadrant.temps = temps;
end

%% FIGURE: compare the arena quadrant occupancy with food well occupancy

% 1) Show the correlation between full quadrant occupancy and food well ROI occupancy
r = 2;
c = 3;
fig = getfig('',1);

for exp = 1:num.exp
subplot(r,c,exp)
  

    x = grouped(exp).time;
    quad = grouped(exp).quadrant.avg;
    foodROI = grouped(exp).occ.avg.*100;

    plot(x, quad, 'color', 'r')
    hold on
    plot(x, foodROI, 'color', 'k')
    xlabel('time (min)')
    ylabel('Occupancy (%)')
    ylim([0,90])

    R = corrcoef(quad,foodROI);
    title({grouped(exp).name; ['R = ' num2str(R(2))]})    
end
formatFig(fig, blkbgd,[r,c]);

% what are the correlations between the two ROI measures?


















