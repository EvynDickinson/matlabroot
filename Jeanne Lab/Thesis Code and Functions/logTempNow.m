

function tempLog = logTempNow(searchPath)
% tempLog = logTempNow(searchPath)
% 
% 
% ES Dickinson, Yale Univeristy, Aug 2021

    list_dirs = dir(searchPath); % only videos
    %I contains the index to the biggest number which is the latest file
    [~,I] = max([list_dirs(:).datenum]);
    if ~isempty(I)
        logName = [list_dirs(I).folder, '\' list_dirs(I).name];
    end
    A = readmatrix(logName);
    tempLog = [size(A,1), A(end,:)]; %index, timepoint, currtemp, settemp, %loading

end