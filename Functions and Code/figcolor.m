

function color_string = figcolor(blkbgd)
% color_string = figcolor(blkbgd)
% returns string: '-blk' or '-wht' depending on the
% background color selection

%%
if blkbgd
    color_string = '-blk';
else 
    color_string = '-wht';
end

