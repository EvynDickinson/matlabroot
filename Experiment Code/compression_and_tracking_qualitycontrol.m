

vidPath = "G:\My Drive\Jeanne Lab\Video Compression Testing\testLossLess.avi";
movieInfo = VideoReader(vidPath); %read in video
demo = rgb2gray(read(movieInfo,1));

%% Compare compression methods with VirtualDub

baseFolder = 'G:\My Drive\Jeanne Lab\Video Compression Testing';
% Select the complete experiments to process
list_dirs = dir([baseFolder '\*.avi']); %only matlab files
list_dirs = {list_dirs(:).name};
nvids = length(list_dirs);

comp_label = {'Full','XVid-2.5', 'XVid-3','XVid-5','XVid-10', 'XVid-15', 'XVid-20', 'Loesless'};
[~,Idx] = sort([0, 10, 15, 2.5, 20, 3, 5, 7.5]);


X = 1:8;
Y = [229, 17, 12, 5, 3, 2, 1.5, 315];

fig = figure; set(fig, 'color', 'k')
plot(X,Y, 'LineWidth', 2, 'Color', Color('teal'), 'Marker','*')
set(gca, 'YScale', 'log')
set(gca, 'XTickLabel', comp_label)
xlabel('Model')
ylabel('Video Size (KB)')
formatFig(fig, true)



% pull the first frame for each of the videos:
for ii = 1:nvids
    vid = Idx(ii);
    movieInfo = VideoReader([baseFolder, '\' list_dirs{vid}]); %read in video
    if strcmp(movieInfo.VideoFormat, 'RGB24')
        fid = fopen([baseFolder, '\' list_dirs{vid}]);
        
        temp = fread(fid);


        demoImg(vid).img = rgb2gray(read(movieInfo,1));

        i=1;
        while hasFrame(movieInfo)
            video_data(i).cdata = readFrame(movie);
            i=i+1;
        end


    else
        demoImg(vid).img = rgb2gray(read(movieInfo,1));
    end

end

nrows = 2;
ncols = 4;

fig = getfig; hold on
for vid = 1:nvids
    subplot(nrows, ncols, vid)
    imshow(demoImg(vid).img)
end

expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
expName = expName(1:end-1);
clear expNames

movieInfo = VideoReader("G:\My Drive\Jeanne Lab\DATA\11.09.2021\Arena B\PlantYeastChoice_1.avi"); %read in video