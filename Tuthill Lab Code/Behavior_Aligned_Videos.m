
clear all; close all; clc 
% Video analysis program to play videos of a fly doing the same thing
% i.e. a fly that is running with leg in stance at 1.0cm/sec etc
% Data selection type: 
camera = 'B';
% behaviors = {'walking', 'stationary', 'grooming', 'other'};
behaviors = {'Walking', 'Stationary', 'Grooming', 'Other'};
% type.behavior = 'grooming';
type.phase = 'obtuse'; %'bent'; %'acute' ; %'stance';
FrameRate = 10;
 
%% Selection parameters:
 
origState = warning;
warning('off')
strtpoint = 470; %line in excel file to start on
behavior_list = {'stationary', 'walking', 'grooming', 'other', 'repeat video!', 'wrong phase last time!'};
phase_list = {'stance', 'swing', 'swing', 'swing', 'WHOOPS'};
% Find files from Excel:    
[excelfile, Excel] = load_flysummary;

% Pull out the structure names from the file: 
structure_names.excelfile = excelfile(strtpoint:end, Excel.structure);
ind = 1;
for ii = 1:length(structure_names.excelfile)
    if ischar(structure_names.excelfile{ii})
        structure_names.text{ind} = (structure_names.excelfile{ii});
        ind = ind+1;
    end
end; clear ind

%find the unique structure names:
structure_names.unique = unique(structure_names.text);
indx = listdlg('ListString', structure_names.unique, 'SelectionMode', 'Multiple', 'ListSize', [300, 400]);

structure_name = structure_names.unique{indx};

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
    
% Load the behavior classification data
% load(['C:\matlabroot\behavior class\' structure_name ' behavior class'])
load(['C:\matlabroot\behavior class\' structure_name ' Group'])


% Load the fly data structure:
load(structure_name)

% FLY(1) = [];
% num.flies = num.flies-1;
% Select the data type: *i.e. walking in swing*  
num.conds = size(group(1).behavior,1);
num.reps = size(group(1).behavior,2);
% calculate the avg speed during the behavior categorization frames:
for kk = 1:num.flies
   for COND = 1:num.conds
       group(kk).speed(COND,:) = nanmean(FLY(kk).Control.speed(COND).data(end-7:end,:));
   end
end
stoptimes = [0.8 0.8 0.8 0.8 0.9 1 1.4 0.8 0.8 0.8 0.8 0.9 1 1.4...
             0.8 0.8 0.8 0.8 0.9 1 1.4 0.8 0.8 0.8 0.8 0.9 1 1.4];
starttime = 0.40;

% starttime = 0.50;
% stoptimes = [0.65 0.65 0.65 0.65 0.65 0.65 0.65 0.8 0.8 0.8 0.8 0.9 1 1.4...
%              0.8 0.8 0.8 0.8 0.9 1 1.4 0.8 0.8 0.8 0.8 0.9 1 1.4];
extralabel = ''; 

for itype = 1:4
    type.behavior = behaviors{itype};
for cond =   1:28 %[21,28]
    clear aa a widthcoordinates heightcoordinates clear cdata vidsize coordinates selectdata
    
    stoptime = stoptimes(cond);
    fprintf(['\n Starting video ' num2str(cond) '\n'])
    
    selectdata.flynum = [];
    for kk = 1:num.flies
        filter.behavior = strcmpi(group(kk).behavior,type.behavior);
        filter.phase = strcmpi(group(kk).ANGLE,type.phase);
        for rep = 1:num.reps
            if filter.behavior(cond,rep)==true && filter.phase(cond,rep)==true 
                selectdata.flynum = [selectdata.flynum; kk, rep];
            end
        end
    end
    
%     
%     selectdata.flynum = [];
%     for kk = 1:num.flies
%         filter.behavior = strcmpi(group(kk).behavior,type.behavior);
%         filter.phase = strcmpi(group(kk).phase,type.phase);
%         for rep = 1:num.reps
%             if filter.behavior(cond,rep)==true && filter.phase(cond,rep)==true 
%                 selectdata.flynum = [selectdata.flynum; kk, rep];
%             end
%         end
%     end
    
    if size(selectdata.flynum,1)>=1
        
        num.vidfiles = size(selectdata.flynum,1);
        % load the name of the fly videos:
        for kk = 1:num.vidfiles
            flynumber = selectdata.flynum(kk,1);
            selectdata.speed(kk) = group(flynumber).speed(selectdata.flynum(kk,2));
            folder_date = flylist{flynumber,1};
            datestring = dateconverter(folder_date);
            flynum = flylist{flynumber,2};
            fileroot = ['D:\Evyn Data Files\' folder_date '\Fly ' flynum '\Raw Video\'];
            filename = [datestring '_fly' flynum ' R' num2str(selectdata.flynum(kk,2))...
                        'C' num2str(cond) ' Cam-' camera ' ' getcond(cond)];
            selectdata.videoname{kk} = [fileroot, filename, '.avi'];
    %         disp(selectdata.videoname{kk})
        end

        % Load the fly videos: 
        duration = round((stoptime-starttime)*300+1);

        % Determine the size and scaling of the video images:
        %total image size is going to be ~screensize:
        sizes.startwidth = 1450;
        sizes.startheight = 900;

        %autofill the dimensions for the videos:
        aa = divisors(num.vidfiles);
        a = sqrt(num.vidfiles);
        if length(aa) == 2 %no divisors
            aa = [2, 3, 4, 5];
            [~,idx] = min(abs(aa-a));
            dim1 = aa(idx);
            dim2 = ceil(num.vidfiles/dim1);
        else 
            [~,idx] = min(abs(aa-a));
            dim1 = aa(idx);
            dim2 = num.vidfiles/dim1;
        end
        if dim1 > 5 || dim2 > 5
            a = max([dim1, dim2]);
            dim1 = 5;
            dim2 = ceil(num.vidfiles/dim1);
        end
        sizes.numcol = max([dim1, dim2]);
        sizes.numrow = min([dim1, dim2]);
        fprintf(['\n Vid dimensions: ' num2str(sizes.numcol) 'x' num2str(sizes.numrow) '\n'])

        sizes.vidheight = round(sizes.startheight/sizes.numrow);

        % create a speed-based index
        [~,speedindx] = sort(selectdata.speed);
        % load all the video cdata into the dtructure cdata 
        % tic
        vidheight = sizes.vidheight;
        for kk = 1:num.vidfiles

            % % Create a VideoReader to use for reading frames from the file.
            flyVideo = VideoReader([selectdata.videoname{speedindx(kk)}]);
            flyVideo.CurrentTime = starttime;
            duration = round((stoptime-flyVideo.CurrentTime)*300+1);
            cdata(kk).mov(duration) = struct('cdata',zeros(flyVideo.Height,flyVideo.Width,3,'uint8'),...
            'colormap',[]);
            ii = 1;
            while flyVideo.CurrentTime <= stoptime
                tempframe = readFrame(flyVideo);
                tempframe1 = imresize(tempframe, [vidheight NaN]);
                cdata(kk).mov(ii).cdata = imresize(tempframe1, 1);
               ii = ii+1;
            end
    %         disp(kk)

        end ; clear vidheight
        % toc
        for kk = 1:num.vidfiles
            vidsize.width(kk) = size(cdata(kk).mov(1).cdata,2);
            vidsize.height(kk) = size(cdata(kk).mov(1).cdata,1);
        end

        vidsize.W = median(vidsize.width);
        vidsize.H = median(vidsize.height);

        sizes.vidwidth = vidsize.W;
        sizes.vidheight = vidsize.H;

        sizes.framewidth = vidsize.W*sizes.numcol;
        sizes.frameheight = vidsize.H*sizes.numrow;

        image_size = ([sizes.frameheight, sizes.framewidth]);
        % size and error check
        for kk = 1:num.vidfiles
            if vidsize.width(kk)>vidsize.W || vidsize.width(kk)<vidsize.W
                err(1) = kk;
                err(2) = vidsize.width(kk)-vidsize.W;
                disp('Size Mismatch')
             % adjust the size of the video to scale with the others:
               blankframe = uint8(zeros(vidsize.H, vidsize.W, 3));
               for frame = 1:duration
                   if err(2) > 0 % larger image
                        rescaledimage = imresize(cdata(err(1)).mov(frame).cdata, [NaN, vidsize.W]);
                   elseif err(2) < 0 % smaller image
                        rescaledimage = cdata(err(1)).mov(frame).cdata;
                   end
                    x = size(rescaledimage,1);
                    y = size(rescaledimage,2);
                    z = size(rescaledimage,3);
                    blankframe(1:x, 1:y, 1:z) = rescaledimage;
                    cdata(err(1)).mov(frame).cdata = [];
                    cdata(err(1)).mov(frame).cdata = blankframe;
               end
            end   
    %         disp(size(cdata(kk).mov(55).cdata))
        end  

        % find row positions:
        edges = [];
        edges(1) = 0;
        for ii = 1:sizes.numrow
            edges(ii+1) = sizes.numcol*ii;
        end

        heightcoordinates = discretize(0:num.vidfiles-1,edges);
        % find col positions:
        for kk = 1:num.vidfiles
            widthcoordinates(kk) = kk-((heightcoordinates(kk)-1)*sizes.numcol);
            if widthcoordinates(kk) == 0
                 widthcoordinates(kk) = sizes.numcol;
            end
        end
        % all coordinates
        subpos.coordinates = [widthcoordinates', heightcoordinates'];
        % pixel coordinates for each subposition:
        % concatenate the image data:
        for kk = 1:num.vidfiles
            % set the coordinates for the images:
            x = sizes.vidwidth*(subpos.coordinates(kk,1)-1)+1;
            y = sizes.vidheight*(subpos.coordinates(kk,2)-1)+1;
            coordinates(kk,:) = [x,y];
        end
        buffer_frame = uint8(zeros(image_size));  

        % Load the video data into a buffer frame:
        newvideo = struct([]);
        for frame = 1:duration

            newframe = buffer_frame;
            for kk = 1:num.vidfiles
                x = (coordinates(kk,1):coordinates(kk,1)+sizes.vidwidth-1);
                y = (coordinates(kk,2):coordinates(kk,2)+sizes.vidheight-1);
                newframe(y,x) = rgb2gray(cdata(kk).mov(frame).cdata);
            end
            newvideo(frame).data = newframe;
        end

        % size(cdata(kk).mov(frame).cdata)

        if sizes.vidwidth ~= size(cdata(1).mov(1).cdata,2) %if not equal...
            warndlg('Widths not aligned')
        end
        if sizes.vidheight ~= size(cdata(1).mov(1).cdata,1) %if not equal...
            warndlg('Heights not aligned')
        end
        if sizes.frameheight ~= (sizes.vidheight*sizes.numrow)
            warn('Total frame height incorrect')
        end
        if sizes.framewidth ~= (sizes.vidwidth*sizes.numcol)
            warn('Total frame width incorrect')
        end

        % size(newvideo(frame).data)
        % size(buffer_frame)
        % f = figure;
        % imshow(newframe)

        % Saving Info:
        analysis_directory = ['E:\Basler Trig\Behavior Aligned Videos\' structure_name '/' type.behavior '\'];
        if ~exist(analysis_directory, 'dir')
            mkdir(analysis_directory);
        end
        addpath(analysis_directory)

        video_name = [structure_name, ' ', type.behavior, ' L1 in ' type.phase ' during ' getcond(cond) extralabel];
        vid_name = [analysis_directory video_name ' uncompressed.avi'];


    %% Make the Video:

    fig = figure; set(fig, 'pos',[300,190,image_size(2),image_size(1)], 'color','k');

    currAxes.Visible = 'off';
    hold on
    v = VideoWriter(vid_name, 'Uncompressed AVI');
    v.FrameRate = FrameRate;
    open(v);
    set(fig, 'color','k');
    imshow(buffer_frame, 'border', 'tight')
    currAxes.Visible = 'off';
    % Collect image info for movie file
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')
    for datapoint = 1:duration
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
    fprintf('\n Video Saved! \n')  

    % Compress the video
    new_vid = [analysis_directory video_name];
    reader = VideoReader(vid_name);        
    writer = VideoWriter(new_vid, 'Motion JPEG AVI'); % 'Motion JPEG AVI', 'MPEG-4'
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
    else
        % no videos for the category, just skip the video
    end
end
end
fprintf('\n DONE!\n')

% 
% % USE THE 'GROUP' STRUCTURE FROM 'Behavior_Categorization.m'
% 
% % Convert behvariors into numbers and then a filter:
% for kk = 1:13
%     for cond = 1:28
%         for rep = 1:3
%             state = group(kk).behavior{cond, rep};
%             phase = group(kk).phase{cond, rep};
%             if ~ischar(state)
%                 state = cell2mat(state);
%             end
%             if ~ischar(phase)
%                 phase = cell2mat(phase);
%             end
%             switch state
%                 case {'stationary', 'Stationary'}
%                     STATE = 1;
%                 case {'walking', 'Walking'}
%                     STATE = 2;
%                 case {'grooming', 'Grooming'}
%                     STATE = 3;
%                 case {'other', 'Other'}
%                     STATE = 4;
%             end
%             switch phase
%                 case {'stance', 'Stance'}
%                     PHASE = 1;
%                 case {'swing', 'Swing'}
%                     PHASE = 2;
%             end 
%             group(kk).STATE(cond, rep) = STATE;
%             group(kk).PHASE(cond, rep) = PHASE;
%         end
%     end
%     group(kk).walking = (group(kk).STATE==2);
%     group(kk).stationary = (group(kk).STATE==1);
%     group(kk).grooming = (group(kk).STATE==3);
%     group(kk).other = (group(kk).STATE==4);
%     group(kk).stance = (group(kk).PHASE==1);
%     group(kk).swing = (group(kk).PHASE==2);
% 
% end
% 
% 



