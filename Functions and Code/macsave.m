
function macsave(fig, figname)
% result = macsave(fig, figname)
% save figures as a pdf on a mac that doesn't have ghostscript

% set(fig, 'pos', [231,215,599,355])
% set(fig, 'pos', [231,215,599,355])

saveas(fig,figname,'pdf');

end

