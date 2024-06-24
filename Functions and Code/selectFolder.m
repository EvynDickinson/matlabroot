

function folderSelected = selectFolder(path,mulitselect)
% folderSelected = selectFolder(path)
% select a folder from an input directory with the autoadded 

if nargin <2
    mulitselect = 'Single'; % default (other option: multiple)
end

dirc = dir(path);
folderNames = ['Today', {dirc(:).name}];
folderNames(2:3) = [];
indx = listdlg('ListString', folderNames, 'SelectionMode', mulitselect);
nFolders = length(indx);
if isempty(indx)
    fprintf('\n No folder selected\n')
    folder = '';
    return
end

dir_sel = cell(1,nFolders);
for i = 1:nFolders
    if indx(i)==1 %TODAY
        dir_sel{i} = char(datetime('Today','Format', 'MM.dd.yyyy'));
    else
        dir_sel{i} = folderNames{indx(i)};
    end
end

folderSelected = dir_sel;

