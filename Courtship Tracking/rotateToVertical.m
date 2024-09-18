

function rotatedPoints = rotateToVertical(points, p1, p2,plot_fig)
    % points: Nx2 matrix of points where each row is [x, y]
    % p1: index of the first point (or [x1, y1] coordinates)
    % p2: index of the second point (or [x2, y2] coordinates)
    

    if ~exist('plot_fig', 'var')
        plot_fig = false;
    end

    % Get the coordinates of the two points
    if isscalar(p1)
        pt1 = points(p1, :);  % Assume p1 is an index
        pt2 = points(p2, :);  % Assume p2 is an index
    else
        pt1 = p1;  % Assume p1 is a [x, y] coordinate
        pt2 = p2;  % Assume p2 is a [x, y] coordinate
    end

    % Calculate the angle between the points and the X-axis
    deltaX = pt2(1) - pt1(1);
    deltaY = pt2(2) - pt1(2);
    theta = atan2(deltaX, deltaY);  % Angle with respect to X-axis

    % Create the rotation matrix to align the points vertically
    R = [cos(theta), -sin(theta); 
         sin(theta),  cos(theta)];
     
    % % Shift points so pt1 is at the origin (for correct rotation)
    % shiftedPoints = points - pt1;
    
    % Apply the rotation matrix to all points
    rotatedPoints = (R * points')';  % Rotate and transpose back to Nx2
    
    % Shift points back to the original position of pt1
    % rotatedPoints = rotatedPoints + pt1;


    % plot fly rotation for new positions: 
    if plot_fig
        fig = getfig('',1,[560 420]);
            hold on
            plotFlySkeleton(fig, points(:,1),points(:,2),'k',true);
            plotFlySkeleton(fig, rotatedPoints(:,1),rotatedPoints(:,2),'r',true);
            % plotFlySkeleton(fig, f(frame, :,1),f(frame, :,2),Color('deeppink'),true);
            axis square equal
            h_line(0,'grey','--',0.5)
            v_line(0,'grey','--',0.5)
    end

end
