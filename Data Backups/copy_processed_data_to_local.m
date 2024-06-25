

% set up directories
curr_folder = 'G:\My Drive\Jeanne Lab\DATA\';
target_dir = 'D:\Processed Data\';

% get list of data folders
dir_list = dir(curr_folder); 
folderNames = {dir_list.name};
filtered_names = folderNames(cellfun(@(x) length(x)==10, folderNames));

% copy the folders 
 for i = 1:length(filtered_names)
     tic
     letter_list = {'A',' B', 'C', 'D'};
     for jj = 1:length(letter_list)
        specific_folder = [filtered_names{i} '\Arena ' letter_list{jj}];
        copyfile([curr_folder specific_folder], [target_dir specific_folder]) 
     end
    disp(['copied ' filtered_names{i}])
    toc
 end


%%
tic
copyfile([curr_folder specific_folder], [target_dir filtered_names{i}]) 
disp('done')
toc