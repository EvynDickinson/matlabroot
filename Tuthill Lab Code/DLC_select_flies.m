
function FilePath = DLC_select_flies(varargin)
% FilePath = DLC_select_flies('-inpath', fileroot, '-outpath', output_root, '-export', true);
% Arguments: 
% '-inpath' --> the location of the folder housing the date folder
%             default: 'G:\My Drive\Evyn\Data\FicTrac Raw Data\'
% '-outpath' --> location to save the exported excel file
%             default: 'C:\matlabroot\DLC\'
% '-export' --> 1 = save a copy of the addresses to an excel file
% 
% ES Dickinson
% University of Washington, 2020

% % -- Select structures to analyze & get the file path info -- % %


% default paths & settings: 
fileroot = 'G:\My Drive\Evyn\Data\FicTrac Raw Data\'; %TODO
output_root = 'C:\matlabroot\DLC\'; %TODO
export_opt = false;

% parse the input parameters:
for ii = 1:nargin
    if ischar(varargin{ii}) && ~isempty(varargin{ii})
        if varargin{ii}(1) == '-' %find command descriptions
            switch lower(varargin{ii}(2:end))
                case 'inpath'
                    fileroot = varargin{ii+1};
                case 'outpath'
                    output_root = varargin{ii+1};
                case 'export'
                    export_opt = varargin{ii+1};
            end
        end
    end
end


% Pull data from fly summary:    
[excelfile, Excel] = load_flysummary; 
strtpoint = 470; %TODO

% Pull out the structure names from the file:
structure_names.excelfile = excelfile(strtpoint:end, Excel.new_struct_name);
structure_names.numberslist = excelfile(strtpoint:end, Excel.structurenum);
ind = 1; 
for ii = 1:length(structure_names.excelfile)
    if ischar(structure_names.excelfile{ii})
        structure_names.text{ind} = (structure_names.excelfile{ii});
        structure_names.num(ind) = (structure_names.numberslist{ii});
        ind = ind+1;
    end
end; clear ind
% (USEFUL)
% A = cellfun(func,C) %apply the function to each cell in the array


%find the unique structure names:
structure_names.unique = unique(structure_names.text);

% select the structure names desired: TODO make this a multi-select and
% choose multiple structures
choice = listdlg('ListString', structure_names.unique, 'PromptString','Structure choice?',...
                'SelectionMode', 'Multi', 'ListSize', [250 400]);
            
% Structure:
structure_name = structure_names.unique{choice(1)};
switch questdlg(['Use structure name: ' structure_name '?'], 'Save Data', 'Yes', 'Different Name', 'No', 'Yes')
    case 'Yes'
    case 'Different Name'
        structure_name = cell2mat(inputdlg('Name to save?'));
    case 'No'
        return
end
full_struct_name = [output_root, structure_name, '_address'];

% Pull the folder information and export to a text file:
% fid = fopen([full_struct_name, '.txt'], 'a+'); % create|open text file
idx = 1;
for icross = 1:length(choice)
    FilePath(icross).structure_name = structure_name;
    cross_name = structure_names.unique{choice(icross)};
    fprintf(['\n Selected cross: ' cross_name '\n'])
    loc_name = find(strcmpi(cross_name, excelfile(:,Excel.new_struct_name)) == 1);
    % check for a structure number:
    loc_number = cell2mat(excelfile(loc_name, Excel.structurenum));
    struct_check = (loc_number>0);
    struct_N = max(loc_number);
    % find folder dates for files with overlap in the struct name & numbers:
    folder_date = (excelfile(loc_name(struct_check), Excel.date));
    fly_nums = excelfile(loc_name(struct_check), Excel.flynum);
      
    % pull all the location information for each structure
    FilePath(icross).cross = cross_name;
    for ifly = 1:struct_N
        FilePath(icross).locations{ifly,1} = fileroot;
        FilePath(icross).locations{ifly,2} = folder_date{ifly};
        FilePath(icross).locations{ifly,3} = ['Fly ', fly_nums{ifly}];
        flyID = generate_flyID(folder_date{ifly}, fly_nums{ifly});
        FilePath(icross).locations{ifly,4} = flyID;
        A_export_locations{idx,1} = ...
            [fileroot, folder_date{ifly}, '\Fly ', fly_nums{ifly}, '\angles\',flyID];
        P_export_locations{idx,1} = ...
            [fileroot, folder_date{ifly}, '\Fly ', fly_nums{ifly}, '\pose-3d\',flyID];                  
%         fwrite(fid, [A_export_locations{idx,1}, ', '],'char'); % write to the text file  
        idx = idx+1;
    end
end
% fclose(fid); %close the text file
% write to Excel File: 

if export_opt == 1
    writecell(A_export_locations,[full_struct_name '.xls'],'Sheet','Angles')
    writecell(P_export_locations,[full_struct_name '.xls'],'Sheet','Pose-3d')
end


end


 





