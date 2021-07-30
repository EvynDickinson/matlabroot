function logData(src, evt, fid)
% Add the time stamp and the data values to data. To write data sequentially,
% transpose the matrix.

data = [evt.TimeStamps, evt.Data]' ;
fwrite(fid,data,'double');
end