
function write_to_excel(param)
% Write experiment information to excel 'Fly Summary' sheet about the 
% details of the experiment
% Inputs:
% 'param' [structure with the variable info to save into excel]
% Outputs:
% Line addended to the excel file summarizing the experiment for each fly
% 
% ES Dickinson, University Of Washington, Dec 2018

% WRITE EXPERIMENT INFO INTO EXCEL FILE
[raw, Excel, xlFile] = load_flysummary;

end_pos = size(raw,1)+1;
xlRange = ['A' num2str(end_pos) ':P' num2str(end_pos)];
%find time:
a = datetime;
a = char(a);
time = a(end-7:end); clear a
% Information to write:

%light lengths
lite = char(string(param.light_length(1)));
for ii = 2:length(param.light_length)
   LTE = char(string(param.light_length(ii)));
   lite = [lite ', ' LTE]; 
end

info = {param.folder_date, param.fly_num, param.cross, param.sex, '5.23.18-A', ...
        num2str(param.num_reps), num2str(param.basler_length), param.LED_light, ...
        num2str(param.LED_intensity), num2str(param.vid_on), lite, time, ...
        num2str(param.control_time), num2str(param.OL_time), cell2mat(param.temp),...
        cell2mat(param.humidity)};

% write to excel file 
xlswrite(xlFile, info, 'Sheet1', xlRange);
 
fprintf('\n Information written to Excel file \n')
close all

% Add comment into 'Notes' section in Excel
answer = questdlg('Add comment on fly behavior?', 'Notes', 'Yes', 'No', 'No');
a = strcmpi('Notes:', Excel.headers); %column number
NoteRange = [Alphabet(a) num2str(end_pos)]; %cell address
switch answer
case 'Yes'
    newNote = inputdlg('Note for Excel');
    Notes = {[param.body ', ' param.loading ', ' newNote]};
    xlswrite(xlFile, Notes, 'Sheet1', NoteRange);
    fprintf('\n Note added to Excel file \n')
case 'No'
    Notes = {[param.body ', ' param.loading]};
    xlswrite(xlFile, Notes, 'Sheet1', NoteRange);
    fprintf('\n No comments on behavior \n')
end

close all

end





