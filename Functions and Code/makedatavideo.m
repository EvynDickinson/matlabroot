



function [] = makedatavideo(FileTag,LowI,HighI)

FileList=dir([FileTag, '*FilterAndRegisterImages.mat']);
load(FileList.name);

FileList=dir([FileTag, '*FilterAndRegisterImagesSelectROICalculateDFF.mat']);
load(FileList.name);

FileList=dir([FileTag, '*AngleForImagingFrame.mat']);
load(FileList.name);

NofRows=size(RegisteredImages,1);
NofColumns=size(RegisteredImages,2);
NofPixels=NofRows*NofColumns;
NofFrames=size(RegisteredImages,4);
framerate = 7.57; %hz, hard coded based on experience

ImageLegAngle(ImageLegAngle<0) = ImageLegAngle(find(ImageLegAngle>0, 1));
if length(ImageLegAngle)<NofFrames
   z = NaN(NofFrames-length(ImageLegAngle), 1);
   ImageLegAngle = [ImageLegAngle; z];
end

Images=squeeze(RegisteredImages(:,:,SignalChannel,:));
% set(0,'DefaultFigureVisible','off'); % Don't display figure

fig = figure;
v = VideoWriter(['G:\My Drive\Sweta to backup\2 photon data\uncompressed data videos\', FileTag, '.avi'], 'Uncompressed AVI');
v.FrameRate = 3*framerate;
% v.Quality = 100;
open(v);

for i = 1:NofFrames
%     i
    h = subplot(2, 1, 2);
    plot(DFF1(1:i))
    ylim([-0.1, max(DFF1)+0.2]);
    ylabel('DF/F')
    xlim([0, NofFrames]);
    xlabel('sec')
    set(h,'color', 'k')
    
    set(gcf,'color', 'k')
    box off
    
    ax = gca;
    xticks(ax, 0:framerate*10:NofFrames)
    ax.XTickLabel = ax.XTick./framerate;
    ax.LineWidth = 1;
    ax.XColorMode = 'manual';
    ax.XColor = 'w';
    ax.YColor = 'w';
    
    hsp1 = get(gca, 'Position');
    
    
    g = subplot(2, 1, 1);
    imshow(Images(5:end-5, 5:end-5, i), [LowI, HighI]);
    set(g, 'Xtick', [], 'YTick', []);
    box off

    hold on 
    %calculate leg diagram coordinates, plot
    x3 = 35 + 20*cos(deg2rad(180-ImageLegAngle(i)));
    y3 = 15 + 20*sin(deg2rad(180-ImageLegAngle(i)));
    line([15, 35, x3], [15, 15, y3], 'Color', 'red', 'LineWidth', 3)
    plot(x3, y3, '.r', 'MarkerSize', 20)
    hsp2 = get(gca, 'Position');
    set(gca, 'Position', [hsp1(1), 0.8*hsp2(2), hsp1(3), 1.5*hsp2(4)])
    
    f = getframe(fig);
    writeVideo(v, f)
end

close(v)


clear