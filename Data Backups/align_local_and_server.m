

function status = align_local_and_server()
% copy trial data files from a local drive to the server and vice versa
% current server is svalbard -- must be connected to VPN if off campus
% ESD Yale 2024

%% Compare local trial data to server data

% select the local drive to compare to the server
switch questdlg('Use standard local drive?')
    case 'Yes'
        local_drive = getBasePath;
    case 'No'
        local_drive = uigetdir;
        local_drive = [local_drive '/'];
    case 'Cancel'
        return b
end

server_drive = getServerPath; 
% TODO: add check that the server is functioning

% get the names of the folders in the single trial portion of the local drive:
localList = dir(local_drive);
localFolders = {localList(:).name};

% get the names of the folders in the single trial portion of the svalbard server:
serverList = dir(server_drive);
serverFolders = {serverList(:).name};

[localSkip, serverSkip] = deal(false);

% files NOT on the local drive:
missinglocalFolders = setdiff(serverFolders,localFolders); 
if isempty(missinglocalFolders)
    localSkip = true;
end

% files NOT on the server:
missingserverFolders = setdiff(localFolders, serverFolders); 
if isempty(missingserverFolders)
    serverSkip = true;
end

% Exit program if all files are found in both places
if serverSkip && localSkip
    disp('All files present in both locations')
    return
end

%% select folders to move between locations:

% display the number of folders that are missing and ask which to copy: 
if ~localSkip
    toLocalIdx = listdlg('PromptString', {'Select folders to copy:'; 'FROM the server'; 'TO the local drive'},...
            'ListString', missinglocalFolders,'ListSize', [250,400]);
else
    toLocalIdx = [];
    localSkip = true;
end

if ~serverSkip
    toServerIdx = listdlg('PromptString', {'Select folders to copy:'; 'FROM the local drive'; 'TO the server'},...
            'ListString', missingserverFolders,'ListSize', [250,400]);
else
    toServerIdx = [];
    serverSkip = true;
end

% Exit program if no files are selected to move between drives
n = length([toLocalIdx, toServerIdx]);
if n==0
    disp('No files selected to transfer between drives')
    return
end

N = num2str(n);

%% Copy folders to appropriate locations: 
% update bar for copying progress
f = waitbar(0,['Copying data from server to local drive...0/' N]);

idx = 0;
% Copy TO local
if ~localSkip
    for i = 1:length(toLocalIdx)
        idx = idx + 1;
        startDir = [server_drive missinglocalFolders{toLocalIdx(i)}];
        endDir = [local_drive missinglocalFolders{toLocalIdx(i)}];
        success = copyfile(startDir, endDir);
        if success
            disp(['Moved ' missinglocalFolders{toLocalIdx(i)}])
        else
            disp(['Failed to move ' missinglocalFolders{toLocalIdx(i)}])
        end
        waitbar(idx/n, f, ['Copying data from server to local drive...' num2str(idx) '/' N])
    end
end

% Copy TO server
if ~serverSkip
    waitbar(idx/n, f, ['Copying data from local drive to server...' num2str(idx) '/' N])
    for i = 1:length(toServerIdx)
        idx = idx + 1;
        startDir = [local_drive missingserverFolders{toServerIdx(i)}];
        endDir = [server_drive missingserverFolders{toServerIdx(i)}];
        success = copyfile(startDir, endDir);
        if success
            disp(['Moved ' missingserverFolders{toServerIdx(i)}])
        else
            disp(['Failed to move ' missingserverFolders{toServerIdx(i)}])
        end
        waitbar(idx/n, f, ['Copying data from local drive to server...' num2str(idx) '/' N])
    end
end
disp('Finished copying all data')
close(f)

status = true;
























