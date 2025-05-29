

function Q = findQuadLocation(x,y)
    % Q = findQuadLocation(x,y)
    % where x and y have been zeroed to the center of the arena
    % Q returns logical masks
    
    Q = struct;
    Q(1).Mask = (y > x) & (y <= -x);  % left quad
    Q(2).Mask = (y <= x) & (y <= -x); % top quad
    Q(3).Mask = (y <= x) & (y > -x);  % right quad
    Q(4).Mask = (y > x) & (y > -x);   % bottom quad

    Q(1).name = 'left quad';
    Q(2).name = 'top quad';
    Q(3).name = 'right quad';
    Q(4).name = 'bottom quad';

end