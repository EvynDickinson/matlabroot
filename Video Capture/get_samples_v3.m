
% get_samples_v3(120)

function get_samples_v3(nframes,fpf,vidName,tempLogName,hz)

    % If you run this function in a loop, the output files will be stored in
    % different directories
    persistent dirname
    % persistent tempLogName
    persistent count
    if isempty(count)
        count = 1;
    else
        count = count + 1;
    end
    % basedir = 'DATA\09.23.2024\';
    dirname = [vidName '_' num2str(count)];
    disp(dirname)
    mkdir(dirname);

    % Clean Up
    delete(imaqfind)

    % Create Video Object
    vi1 = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');
    src = getselectedsource(vi1);  % Get the camera source properties

    [~, vi1] = initialize_CourtshipCamera(src,vi1,hz);

    % Initialize Counter
    vi1.UserData = 1; 

    % Set Parameters
    set(vi1,'FramesAcquiredFcn',{@FrameSave,dirname,tempLogName},'FramesAcquiredFcnCount',fpf);
    set(vi1,'FramesPerTrigger',nframes,'LoggingMode','memory');
    % Dynamically adjust timeout based on frame rate
    expected_capture_time = nframes / hz + 5;  % Adding buffer time
    set(vi1,'Timeout',expected_capture_time);

    % Trigger Device
    start(vi1);

    try
        wait(vi1);
    catch ME
        warning('%s: %s', ME.identifier, ME.message);
    end
end


% Callback every fpf frames
function FrameSave(vi1,~,dirname,tempLogName)
    try
        data = getdata(vi1,vi1.FramesAcquiredFcnCount);
        filename = ['file',num2str(vi1.UserData),'.mat'];
        timestamp = datetime('now');

        % log current temp information
        A = readmatrix(tempLogName);
        tempLog = [size(A,1), A(end,:)];

        % save the video file and temp and time
        save([pwd '\' dirname '\' filename],'data','timestamp','tempLog','-v6');
        vi1.UserData = vi1.UserData + 1;
    catch ME
         warning('%s: %s', ME.identifier, ME.message);
    end
end

 % -v6 option speeds up save process to help prevent the buffer 
    % from filling during file write operation

% 
% function tempLog = logTempNow(searchPath)
% % tempLog = logTempNow(searchPath)
% % 
% % 
% % ES Dickinson, Yale Univeristy, Aug 2021
% 
%     list_dirs = dir(searchPath); % only videos
%     %I contains the index to the biggest number which is the latest file
%     [~,I] = max([list_dirs(:).datenum]);
%     if ~isempty(I)
%         logName = [list_dirs(I).folder, '\' list_dirs(I).name];
%     end
%     A = readmatrix(logName);
%     tempLog = [size(A,1), A(end,:)]; %index, timepoint, currtemp, settemp, %loading
% 
% end














