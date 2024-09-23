
% get_samples_v3(120)

function get_samples_v3(nframes,fpf)

% If you run this function in a loop, the output files will be stored in
% different directories
persistent dirname
persistent count
if isempty(count)
    count = 1;
else
    count = count + 1;
end
% basedir = 'DATA\09.23.2024\';
dirname = ['output_' num2str(count)];
disp(dirname)

mkdir(dirname);
% Clean Up
delete(imaqfind)

% Create Video Object
% vi1 = videoinput('winvideo');

vi1 = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');

% Initialize Counter
vi1.UserData = 1; 
% fpf = 60;

% Set Parameters
set(vi1,'FramesAcquiredFcn',{@FrameSave,dirname},'FramesAcquiredFcnCount',fpf);
set(vi1,'FramesPerTrigger',nframes,'LoggingMode','memory');
set(vi1,'Timeout',fpf);

% Trigger Device
start(vi1);
wait(vi1);


% Callback Every fpf Frames
function FrameSave(vi1,event,dirname)
disp(vi1.FramesAcquired);
data = getdata(vi1,vi1.FramesAcquiredFcnCount);
filename = ['file',num2str(vi1.UserData),'.mat'];
timestamp = datetime('now');
% -v6 option speeds up save process to help prevent the buffer from filling
% during file write operation
filename
save([pwd '\' dirname '\' filename],'data','timestamp','-v6');
vi1.UserData = vi1.UserData + 1;

% You can create several AVI files from the MAT files using the AVIFILE and 
% ADDFRAME functions in MATLAB,and a third-party utility such as 
% VirtualDub (<http://www.virtualdub.org/>) can create one large AVI file.

