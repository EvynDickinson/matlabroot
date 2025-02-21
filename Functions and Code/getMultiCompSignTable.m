
function [fig, cPlot, pPlot] = getMultiCompSignTable(c, expOrder, blkbgd, thresh, names)

if nargin <3
    blkbgd = true;
end
[foreColor,~] = formattingColors(blkbgd); %get background colors
n = length(unique(c(:,1))) + 1; % how many groups to plot for comparisons

if nargin <2
    expOrder = 1:n;
end
if nargin <4
    thresh = 0.05; % significance level
end
if nargin <5
    names = 1:n;
end

buff = 0.5; % spacing for plotting (1/2 column width)
sz = 150; % plot symbol size

% display the statistical differences: 
    
    cPlot = []; pPlot = [];
    for i = 1:size(c,1)
        ex1 = expOrder(c(i,1));
        ex2 = expOrder(c(i,2));
        cPlot(ex1, ex2) = (c(i,6));
        cPlot(ex2, ex1) = c(i,6);
        if c(i,6)<=thresh
            pPlot = [pPlot; ex1,ex2; ex2,ex1];
        end
    end

    lims = [buff,n+buff];
    fig = getfig('',1,[705 649]);
    % imagesc(cPlot);
    % axis square equal
    xlim(lims)
    ylim(lims)
    hold on
    scatter(pPlot(:,1), pPlot(:,2), sz, foreColor, '*')
    v_line(buff:1:n+buff,'grey','-',1)
    h_line(buff:1:n+buff,'grey','-',1)
    plot(lims,lims)

    formatFig(fig, blkbgd);
    set(gca,'XTick',1:1:n,'YTick',1:1:n,'XTickLabel',names,'YTickLabel',names)
    set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
    
