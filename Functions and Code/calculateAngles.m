

function Angles = calculateAngles(JointPositions)
% Angles = calculateAngles(JointPositions)
% joint positions: x,y matrix

numAngles = size(JointPositions,1)-2;
x = JointPositions(:,1);
y = JointPositions(:,2);
Angles = NaN(1,numAngles);

for ii = 1:numAngles
    % select the three positions
    a = [x(ii), y(ii)];
    b = [x(ii+1), y(ii+1)];
    c = [x(ii+2), y(ii+2)];
    % calculate the angle
    u(1) = (b(1)-a(1));
    u(2) = (b(2)-a(2));
    u(3) = 0;
    v(1) = (c(1)-b(1));
    v(2) = (c(2)-b(2));
    v(3) = 0;
    % set the joint angle
    Angles(ii) = 180-atan2d(norm(cross(u,v)),dot(u,v));
end


end