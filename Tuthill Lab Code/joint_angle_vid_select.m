

% 
% fly_cross_list = {'81A07-gal4xUAS-csChrimson', '81A07-gal4xUAS-gtACR1', ...
%                   '22A08-gal4xUAS-gtACR1', ...
%                   '35C09-gal4xUAS-csChrimson', '35C09-gal4xUAS-gtACR1', ...
%                   'BDP-gal4xUAS-gtACR1'};


origState = warning;
warning('off')
strtpoint = 627; %line in excel file to start on
% Find files from Excel:    
[excelfile, Excel] = load_flysummary;             
% 
% 
% % Pull out the structure names from the file: 
% structure_names.excelfile = excelfile(strtpoint:end, Excel.structure);
% ind = 1;
% for ii = 1:length(structure_names.excelfile)
%     if ischar(structure_names.excelfile{ii})
%         structure_names.text{ind} = (structure_names.excelfile{ii});
%         ind = ind+1;
%     end
% end; clear ind
% 
% %find the unique structure names:
% structure_names.unique = unique(structure_names.text);
% indx = listdlg('ListString', structure_names.unique, 'SelectionMode', 'Multiple', 'ListSize', [300, 400]);
% 
% structure_name = structure_names.unique{indx};

names = {'13B-csChrimson-headless-offball', '13B-gtACR1-headless-offball',...
         '10B-csChrimson-headless-offball', '10B-gtACR1-headless-offball',...
         '9A-csChrimson-headless-offball', '9A-gtACR1-headless-offball',...
         '13B-csChrimson-headless-onball', '13B-gtACR1-headless-onball',...
         '10B-csChrimson-headless-onball', '10B-gtACR1-headless-onball',...
         '9A-csChrimson-headless-onball', '9A-gtACR1-headless-onball', ...
         'BDP-gal4xUAS-gtACR1-headless-offball','BDP-gal4xUAS-csChrimson-headless-offball'};


% names = {'22A08xcsChrimson-offball-headless', '22A08xgtACR1-headless-offball',...
%          '35C09xcsChrimson-offball-headless', '35C09xgtACR1-headless-offball',...
%          '81A07xgtACR1-headless-offball', 'BDP-gal4xUAS-gtACR1-headless-offball',...
%          'BDP-gal4xUAS-csChrimson-headless-offball'};

for iname = 1:14
    structure_name = names{iname};
    location = find(strcmpi(structure_name, excelfile(strtpoint:end,Excel.new_struct_name)))+strtpoint-1;
    bb = 1;
    for tt = 1:length(location)
        a = cell2mat(excelfile(location(tt), Excel.structurenum));
        if ischar(a) %no value for the structure
        else %has a number in the structure
            flynum(bb) = a; %find total number of flies
            Location.rownum(bb) = location(tt); %get the rownumber for each fly
            bb = bb+1; %increase the index for the number of flies
        end
    end
    num.flies = max(flynum); clear flynum bb a tt
    % Error Check
    if sum(num.flies) == 0
        warndlg('Check name of flies in Excel File, no matching files found')
        return
    end

    % Load the select Excel file information for the flies in the structure:
    flylist = excelfile(location,:);

%% only test the longest light lenght:
% cond 7, 14, 21, 28

starttime = 0.3;
stoptime = 1.5;

savefolder = ['F:\Joint angle annotation\' structure_name, '\'];
if ~isfolder(savefolder)
    mkdir(savefolder)
end 
% load the name of the fly videos:
for kk = 1:num.flies
    idx = 0;
    for cond = [7, 14, 21, 28]
        idx = idx +1;
        folder_date = flylist{kk,1};
        datestring = dateconverter(folder_date);
        flynum = flylist{kk,2};
        fileroot = ['D:\Evyn Data Files\' folder_date '\Fly ' flynum '\Raw Video\'];
        filename = [datestring '_fly' flynum ' R1C' num2str(cond) ' Cam-E ' getcond(cond)];
        selectdata.videoname{kk,idx} = [fileroot, filename, '.avi'];

        selectdata.newname{kk,idx} = [savefolder '\fly ' num2str(kk) ' cond ' num2str(cond) ' uncompressed.avi'];
        selectdata.compressedname{kk,idx} = [savefolder '\fly ' num2str(kk) ' cond ' num2str(cond) '.avi'];
    end
end

% Load the fly videos: 
duration = round((stoptime-starttime)*300+1);
FrameRate = 30;

for kk = 1:num.flies %each fly in the condition
    for ii = 1:4 %four conditions
        % % Create a VideoReader to use for reading frames from the file.
        flyVideo = VideoReader(selectdata.videoname{kk,ii});
        flyVideo.CurrentTime = starttime;
        cdata.mov(duration) = struct('cdata',zeros(flyVideo.Height,flyVideo.Width,3,'uint8'),...
        'colormap',[]);
        idx = 1;
        while flyVideo.CurrentTime <= stoptime
            tempframe = readFrame(flyVideo);
            newvideo(idx).data = tempframe;
           idx = idx+1;
        end
        fig = figure; set(fig, 'pos',[300,190,flyVideo.Width,flyVideo.Height], 'color','k');
        currAxes.Visible = 'off';
        hold on
        v = VideoWriter(selectdata.newname{kk,ii}, 'Uncompressed AVI');
        v.FrameRate = FrameRate;
        open(v);
        set(fig, 'color','k');
%         imshow(buffer_frame, 'border', 'tight')
        currAxes.Visible = 'off';
        % Collect image info for movie file
        f = getframe(fig);
        writeVideo(v, f)  
        clf('reset')
        for datapoint = 1:5:duration
            set(fig, 'color','k');
            %display the iamge
            imshow(newvideo(datapoint).data, 'border', 'tight')
            currAxes.Visible = 'off';
            % Collect image info for movie file
            f = getframe(fig);
            writeVideo(v, f)  
            clf('reset')
        end
        close(v)
        close all  

        % Compress the video
        reader = VideoReader(selectdata.newname{kk,ii});        
        writer = VideoWriter(selectdata.compressedname{kk,ii}, 'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
        writer.FrameRate = reader.FrameRate;
        writer.Quality = 60;
        open(writer);
        % Read and write each frame 
        while hasFrame(reader)  
            img = readFrame(reader);
            writeVideo(writer,img);
        end
        writeVideo(writer,img);
        close(writer);  
        disp(ii)
    end  
    fprintf(['\n Finished fly ' num2str(kk)])
end
fprintf(['\n Finished: ' structure_name '\n'])

        
end    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        





