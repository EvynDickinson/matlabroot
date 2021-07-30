function fig = getfig(name, ~)
% fig = getfig(name, leftscreenbinarychoice)


fig = figure; 
set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900]);
try
    set(fig, 'name', name);
catch
    set(fig, 'name', 'Fancy Figure');
end

if nargin == 2 %left screen image
%     set(fig, 'pos', [-1909,52,1901,938]);
    set(fig, 'pos', [-1044 261 997 724]);
end

box off

end