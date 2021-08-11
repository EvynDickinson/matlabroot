% MIDSTIM_Step_2_sort_and_condense_data  
  
% Select folder to reorganize and compress: 
Input_dir.root = 'E:\FicTrac Raw Data';   
temp.listing = dir(Input_dir.root);
for ii = 3:length(temp.listing)-1                
    temp.folders{ii-2}= temp.listing(ii).name; 
end  
index = listdlg('ListString', temp.folders, 'PromptString',...      
                'Select folder to organize', 'SelectionMode', 'Single');
Input_dir.date = temp.folders{index};
Input_dir.rootpath = [Input_dir.root, '\', Input_dir.date, '\'];
clear temp

% Select fly trials:  
temp.listing = dir(Input_dir.rootpath); 
for ii = 3:length(temp.listing) 
    temp.folders{ii-2}= temp.listing(ii).name;  
end 
index = listdlg('ListString', temp.folders, 'PromptString',... 
                'Select flies to copy', 'SelectionMode', 'Multiple'); 
     
% Create new folders in the 'D' data storage folder and 'G' google drive
% paths.D_dir = ['D:\Evyn Data Files\' Input_dir.date '\'];
paths.D_dir = ['F:\Evyn Data Files\' Input_dir.date '\'];
paths.G_dir = ['G:\My Drive\Data\FicTrac Raw Data\' Input_dir.date '\'];

% for each fly trial, copy the data over and compress the videos
idx = 0; 
 for prep = 1:length(index)   
    idx = idx+1;
    temp.fly_dir = temp.folders{index(idx)};
   % Data Folders:
    temp.data_root = [temp.fly_dir, '\FicTrac Data'];
    Input_dir.data_root{idx} = [Input_dir.rootpath, temp.data_root];
    Output_dir.D_dir.data_root{idx} = [paths.D_dir, temp.data_root];
    Output_dir.G_dir.data_root{idx} = [paths.G_dir, temp.data_root];
    
    copyfile(Input_dir.data_root{idx}, Output_dir.D_dir.data_root{idx})
    copyfile(Input_dir.data_root{idx}, Output_dir.G_dir.data_root{idx})
    
   % Video Folders:
    
    Input_dir.vid_root{idx} = [Input_dir.rootpath, temp.fly_dir, '\Uncompressed Video'];
    temp.vid_root = [temp.fly_dir, '\Raw Video'];
    Output_dir.D_dir.vid_root{idx} = [paths.D_dir, temp.vid_root];
    Output_dir.G_dir.vid_root{idx} = [paths.G_dir, temp.vid_root];
    if ~exist(Output_dir.D_dir.vid_root{idx})
        mkdir(Output_dir.D_dir.vid_root{idx})
    end
    if ~exist(Output_dir.G_dir.vid_root{idx})
        mkdir(Output_dir.G_dir.vid_root{idx})
    end
    
    % *** Compress videos *** %
    % Create waitbar
    list_videos = dir([Input_dir.vid_root{idx},'\*.avi']);
    nVideos = length(list_videos);
    h = waitbar(0,['Compressing videos from ', Input_dir.date ' ' temp.fly_dir]);

    % Iterate through videos
    for iVideo = 1:nVideos
        % Create objects to read and write the video, and open the AVI file
        % for writing
        disp(['Compressing ', list_videos(iVideo).name,' ...'])
        reader = VideoReader([Input_dir.vid_root{idx}, '\', list_videos(iVideo).name]);        
        writer = VideoWriter([Output_dir.D_dir.vid_root{idx}, '\', list_videos(iVideo).name],'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
        writer.FrameRate = reader.FrameRate;
        writer.Quality = 75;
        open(writer);

        % Read and write each frame
        while hasFrame(reader)
            img = readFrame(reader);
            writeVideo(writer,img);
        end
        writeVideo(writer,img);
        close(writer);           
        disp('Done!')

        % Update waitbar
        waitbar(iVideo/nVideos,h)
    end

    close(h)

    % Copy the Compressed videos into 
    copyfile(Output_dir.D_dir.vid_root{idx},  Output_dir.G_dir.vid_root{idx})
    fprintf(['\n Finished ' temp.fly_dir '\n'])
    % Copy the config file to the video folder:
    copyfile([Input_dir.data_root{idx} '\Config.toml'], Output_dir.G_dir.vid_root{idx})
    copyfile([Input_dir.data_root{idx} '\Config.toml'], Output_dir.D_dir.vid_root{idx})
 end
    
fprintf('\n DONE!! \n')



% Next Step:
% MIDSTIM_Step_3_combine_fictrac_n_matlab







