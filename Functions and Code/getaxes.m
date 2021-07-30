
function fig = getaxes(fig, weight)

if nargin == 1
    weight = 25;
end
ax = gca;
set(ax, 'FontWeight', 'Bold', 'FontSize', weight, 'XColor', 'k', 'YColor', 'k')

