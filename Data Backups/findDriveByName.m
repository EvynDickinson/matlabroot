
% ++++++++++++++++++++++++++++++++++++
function matchingDrives = findDriveByName(targetName)
% matchingDrives = findDriveByName(targetName)
% try to find any drives that have the portable data drives:
%
% targetName = 'OnTheGoData';
%
% 

if ispc
    % Get drive letters and volume names using WMIC
    [~, result] = system('wmic volume get DriveLetter,Label /format:csv');
    
    % Split the result into lines
    lines = strsplit(result, '\n');
    
    % Initialize output
    matchingDrives = {};
    
    % Iterate over each line and check for the target volume name
    for i = 2:length(lines) % Skip the header line
        line = strtrim(lines{i});
        if isempty(line)
            continue;
        end
        
        % Split the line by commas
        parts = strsplit(line, ',');
        
        if length(parts) < 3
            continue; % Skip malformed lines
        end
        
        % Extract the drive letter and volume label
        driveLetter = strtrim(parts{2});
        volumeLabel = strtrim(parts{3});
        
        % Check if the volume name matches the target name
        if strcmpi(volumeLabel, targetName)
            matchingDrives{end+1} = driveLetter; 
        end
    end
end

% FOR MAC COMPUTERS:
if ismac
    % Get volume information using the diskutil command
    [~, result] = system('diskutil info -all');
    
    % Split the result into lines
    lines = strsplit(result, '\n');
    
    % Initialize variables
    matchingDrives = {};
    currentVolume = '';
    currentMountPoint = '';
    
    % Iterate over each line and check for the target volume name
    for i  = 1:length(lines)
        line = strtrim(lines{i});
        
        % Check if the line contains volume name information
        if startsWith(line, '   Volume Name:')
            % Extract the volume name
            currentVolume = strtrim(strrep(line, 'Volume Name:', ''));
        end
        
        % Check if the line contains mount point information
        if startsWith(line, '   Mount Point:')
            % Extract the mount point
            currentMountPoint = strtrim(strrep(line, 'Mount Point:', ''));
        end
        
        % If both volume name and mount point are found, check for a match
        if ~isempty(currentVolume) && ~isempty(currentMountPoint)
            if strcmpi(currentVolume, targetName)
                matchingDrives{end+1} = currentMountPoint; %#ok<AGROW>
            end
            % Reset for the next volume
            currentVolume = '';
            currentMountPoint = '';
        end
    end

end


end

