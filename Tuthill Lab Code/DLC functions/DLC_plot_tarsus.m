
function [fig1, fig2] = DLC_plot_tarsus(data) 
% 
% [fig1, fig2] = DLC_plot_tarsus(data); 
% ES Dickinson
% University of Washington, 2020

leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
for leg = 1:length(leg_colors)
    C(leg,:) = Color(leg_colors{leg});
end

fig1 = getfig('',1); 
% --- view the shifted traces -----
for leg = 1:6
    a = data(leg).input;
    plot3(a(:,1), a(:,2), a(:,3), 'color', C(leg,:), 'linewidth', 1)
    hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
grid on
legend(leg_labels);
title('Shifted Tarsus position')
% ---------------------------------

sbpt = [1:2:6,2:2:6];

fig2 = getfig('',1); 
% --- view the individual traces -----
for leg = 1:6
    subplot(3,2,sbpt(leg))
    a = data(leg).ind;
    plot3(a(:,1), a(:,2), a(:,3), 'color', C(leg,:), 'linewidth', 1)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    title([leg_labels{leg} ' Tarsus position'])
end

% ---------------------------------



end