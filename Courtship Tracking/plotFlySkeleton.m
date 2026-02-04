
function fig = plotFlySkeleton(fig, x, y, kolor, nodes)
% fig = plotFlySkeleton(fig, x,y,kolor,nodes)
%
%




%%

if ~exist('kolor', 'var')
    kolor = 'k';
end
skeleton = [1,2; 2,3; 2,4; 2,5];

figure(fig)
hold on
% edges
for i = 1:size(skeleton,1)
    plot(x(skeleton(i,:)),y(skeleton(i,:)), 'color', kolor,'HandleVisibility','off')
end
% nodes
if nodes 
    scatter(x, y,50, kolor, 'filled')
end
