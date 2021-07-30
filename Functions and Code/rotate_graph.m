
function [X, Y] = rotate_graph(fig, x, y, linecolor)

min_change = 0.015;


% create position offset to (0,0) start point at t0
temp.x_off = x(1);
temp.y_off = y(1);
X = x-temp.x_off;
Y = y-temp.y_off;

X = smooth(X,5);
Y = smooth(Y,5);

% find the first step taken by the fly to use as an alignment:
x_change = abs(diff(X));
y_change = abs(diff(Y));
DistPerStep = sqrt(x_change.^2 + y_change.^2);
time = find(DistPerStep>=min_change, 1 );
if sum(time) <= 1
    time = 2;
end

% X = x-temp.x_off;
% Y = y-temp.y_off;

% find the degree of rotation value for the plot
% time = 2;
d = sqrt(X(time)^2 + Y(time)^2); %distance from t0-t1
u = [X(time), Y(time)];
v = [0, d];
theta = rad2deg(acos(dot(u,v)/d^2));

if X(time)<0 
    theta = theta*-1;
end

direction = [0,0,1];
origin = [0,0,0];
 
figure(fig)
hold on
h = plot(X, Y, 'color', linecolor);
rotate(h, direction, theta, origin)
X = get(h, 'Xdata');
Y = get(h, 'Ydata');


