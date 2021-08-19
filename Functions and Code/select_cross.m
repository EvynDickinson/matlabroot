

function fly_cross = select_cross
% 
% fly_cross = select_cross
% 
% Select the genetic cross of the fly
% List of crosses available in the text
% file 'Cross' in the matlabroot dir.
% 
% Output:
% 'fly_cross' [string name of cross]
% 
% ES Dickinson, University of Washington, 2019


file_id = fopen('Cross.txt');
A = fscanf(file_id, '%s');
Cross = strsplit(A,','); 

idx = listdlg('ListString', Cross, 'PromptString', 'Select fly cross',...
              'SelectionMode', 'Single', 'ListSize', [250, 400]);
if isempty(idx)
    fprintf('\n No cross selected\n')
    fly_cross = 'none';
else
    fly_cross = Cross{idx};
end

end