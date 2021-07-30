
function status = write_OdorShockTrial_excel(param, phase, conc, odor, session)
% status = write_OdorShockTrial_excel(param, phase, conc, odor, session)
% 
% ES Dickinson
% Yale 2020


% Write info to excel
[excelfile, ~, xlFile] = load_ExperimentSummary;
end_pos = size(excelfile,1)+1;
xlRange = ['A' num2str(end_pos) ':Q' num2str(end_pos)];
%find time:
a = char(datetime);
a = a(end-7:end);   
time = strrep(a,':','.'); clear a

info = {param.date, param.fly_num,...
        param.cross, param.sex,...
        param.age, phase, ...
        conc, odor,...
        num2str(param.ODOR_dur),  num2str(param.SHOCK_dur),...
        num2str(param.ODOR_start), num2str(param.SHOCK_start),...
        num2str(param.postStim_dur),...
        param.folder_name, param.tag,...
        session, time};

% write to excel file 
status = xlswrite(xlFile, info, 'Sheet2', xlRange);
fprintf('Information written to Excel file\n\n')

end