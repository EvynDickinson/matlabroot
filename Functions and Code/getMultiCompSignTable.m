


function [fig, cPlot, pPlot] = getMultiCompSignTable(c, expOrder, blkbgd, thresh, names, namematch)
% [fig, cPlot, pPlot] = getMultiCompSignTable(c, expOrder, blkbgd, thresh, names, namematch)
%
% c = the stats list structure
% expOrder = the desired order within the matrix, default is 1:n
% blkbng = T|F for the desired background color (true=black)
% thresh = significance level of p-value for plotting a significance star (default 0.05)
% names = label names for each category (default 1:n)
% namematch = T|F ; match the name order to the experiment order, (e.g., if both need
% to be ordered by the expOrder, set this to true. If false, it will pair expOrder(1)
% with name(1)) Default is false -- e.g. that the names do NOT need to rearranged
%
% ES Dickinson

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
if nargin == 6 && namematch
    nameList = names(expOrder);
else
    nameList = names;
end

buff = 0.5; % spacing for plotting (1/2 column width)
sz = 700; % plot symbol size
lw = 4;
if blkbgd
    lcolor = 'white';
else
    lcolor = 'black';
end

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
    if ~isempty(pPlot)
        scatter(pPlot(:,1), pPlot(:,2), sz, foreColor,'filled', 'pentagram')
    end
    v_line(buff:1:n+buff,lcolor,'-',lw)
    h_line(buff:1:n+buff,lcolor,'-',lw)
    %     v_line(buff:1:n+buff,'grey','-',lw)
    % h_line(buff:1:n+buff,'grey','-',lw)
    plot(lims,lims,'color', Color('gold'),'linewidth', 1.5)

    formatFig(fig, blkbgd);
    set(gca,'XTick',1:1:n,'YTick',1:1:n,'XTickLabel',nameList,'YTickLabel',nameList)
    set(gca, 'TickLength',[0,0],'Box','on','XDir', 'reverse', 'YDir','normal','XTickLabelRotation',90)
    set(gca, 'LineWidth', lw+1)  %change box lines to 1 larger than inner lines
