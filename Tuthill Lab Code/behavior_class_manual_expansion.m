
filename = 'SH-gal4xUAS-gtACR1 behavior class.mat';

load(filename)


for ii = 6:12
    group(ii).behavior = cell(28,3);
    group(ii).phase = cell(28,3);
end

save(filename, 'group')