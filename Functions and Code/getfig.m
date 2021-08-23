function fig = getfig(name, ~)
% fig = getfig(name, leftscreenbinarychoice)

 switch getenv('COMPUTERNAME')
    case 'DENALI'
        basePos = [50, 50, 1450, 900];
        outPos = [2022 161 1232 755];
    case 'TOGIAK'
        basePos = [50, 50, 1450, 900];
        outPos = [-1044 261 997 724];
    case 'EVYNPC'
        basePos = [50, 50, 1450, 900];
        outPos = [-1044 261 997 724];
end
    


fig = figure; 
    set(fig, 'color', 'w', 'pos', basePos)
    try
        set(fig, 'name', name);
    catch
        set(fig, 'name', 'Fancy Figure');
    end
    if nargin == 2 %secondary screen positioning
        set(fig, 'pos', outPos);
    end
    box off

end