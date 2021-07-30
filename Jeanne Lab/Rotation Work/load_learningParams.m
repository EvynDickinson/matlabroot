function param = load_learningParams
odor_opt = {'paraffin oil','2-butanone', '3-octanol'};
conc_opt = {'0.1%','1%','2%','100%', 'other'};

param.date = Num2Month(date);    % today's date
param.cross = select_cross;        % genotype


% solvent
choice = listdlg('ListString', odor_opt, 'PromptString','Select SOLVENT',...
                'SelectionMode', 'Single', 'ListSize', [150 100]); %
param.solvent = odor_opt{choice};
choice = listdlg('ListString', conc_opt, 'PromptString',[param.solvent ' concentration?'],...
                'SelectionMode', 'Single', 'ListSize', [150 100]); %
param.solvent_conc = conc_opt{choice};
if strcmpi(param.solvent_conc,'other')
    param.solvent_conc = char(inputdlg([param.solvent ' concentration?'], 'SOLVENT CONCENTRATION %'));
    param.solvent_conc = [param.solvent_conc '%'];
end


% Conditioned Stimulus
choice = listdlg('ListString', odor_opt, 'PromptString','Select CONDITIONED odor',...
                'SelectionMode', 'Single', 'ListSize', [150 100]); %
param.csPOS_odor = odor_opt{choice};
choice = listdlg('ListString', conc_opt, 'PromptString',[param.csPOS_odor ' concentration?'],...
                'SelectionMode', 'Single', 'ListSize', [150 100]); %
param.csPOS_odor_conc = conc_opt{choice};
if strcmpi(param.csPOS_odor_conc,'other')
    param.csPOS_odor_conc = char(inputdlg([param.csPOS_odor ' concentration?'], 'CONDITIONED ODOR CONCENTRATION'));
    param.csPOS_odor_conc = [param.csPOS_odor_conc '%'];
end



% Unconditioned Stimulus
choice = listdlg('ListString', odor_opt, 'PromptString','Select UNCONDITIONED odor',...
                'SelectionMode', 'Single', 'ListSize', [150 100]); %
param.csNEG_odor = odor_opt{choice};
choice = listdlg('ListString', conc_opt, 'PromptString',[param.csNEG_odor ' concentration?'],...
                'SelectionMode', 'Single', 'ListSize', [150 100]); %
param.csNEG_odor_conc = conc_opt{choice};
if strcmpi(param.csNEG_odor_conc,'other')
    param.csNEG_odor_conc = char(inputdlg([param.csNEG_odor ' concentration?'], 'UNCONDITIONED ODOR CONCENTRATION'));
    param.csNEG_odor_conc = [param.csNEG_odor_conc '%'];
end


% fly age
age_opts = {'1 day','2 days','3 days','4 days','5 days','6 days','7 days'};
choice = listdlg('ListString', age_opts, 'PromptString','Select fly age',...
                'SelectionMode', 'Single', 'ListSize', [100 100]); %
param.age = age_opts{choice};



% Get information for today's experiment
% FormatOut = 'mmddyyyy';
% Today = datestr(datetime,FormatOut);
% param.matlab_data_file = [Today '_fly' param.fly_num];


% param.experiment = 'control';
% param.odor = '2-butanone';
% param.conc = '1%';

end