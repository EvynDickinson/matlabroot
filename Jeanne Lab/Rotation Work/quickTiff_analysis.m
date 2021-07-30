% Load a Tiff stack and then look at the avg intensity change across time
% and across pixel

%% Load TIFF image stack to analyze
clear

warning('off','all') % Suppress all the tiff warnings
folderdate = '9.10.20';
fileRoot = ['E:\Evyn\' folderdate '\'];

% select specific trial:
fileList = dir([fileRoot '*.tif*']);
for ii = 1:length(fileList)
    TiffFiles{ii} = fileList(ii).name;
end
idx = listdlg('PromptString', 'Select Video', 'ListString', TiffFiles);
vidName = TiffFiles{idx};

% load stack:
tstack = Tiff(fullfile(fileRoot,vidName),'r');
[I,J] = size(tstack.read()); % get the pixel sizes for frame
nframes = length(imfinfo(fullfile(fileRoot,vidName))); % get the total number of frames
IM = zeros(I,J,nframes);
IM(:,:,1)  = tstack.read();
disp('Loading frames...')
for n = 2:nframes
    tstack.nextDirectory()
    IM(:,:,n) = tstack.read();
    if rem(n,50)==0 %indicate every 50th frame
        disp(n)
    end
end


%% quick visualization:

%plot avg intensity of entire image over time:
time = 1:nframes;
a = mean(mean(IM,1),2);
y = squeeze(a);
figure;
plot(time,y)


% find region with most change within a given time ROI:
t1 = 200;
t2 = 800;
selIM = IM(:,:,t1:t2);
IMrange = range(selIM,3);

figure;
imagesc(IMrange)
colorbar


%% Plot df/f for a given trial: * IN PROGRESS *
% ROIs
controlStart = 200;      % first frame in video to use
controlStop = 400;      % end of control baseline F region
finalFrame = 800;       % final frame in stack to include

% Background florescence
control_ROI = controlStart:controlStop; 
controlIM = IM(:,:,control_ROI);
F0 = mean(mean(mean(controlIM,3)));

% Find deltaF over F:
pdata = IM(:,:,controlStart:finalFrame);
a = mean(mean(pdata,1),2);
Ft = squeeze(a);
df_f = ((Ft-F0)/F0)*100;

% plot
fig = getfig('Change in florescence',1);
plot(df_f,'color', Color('teal'), 'linewidth', 2)

%labels etc:
ylabel('\DeltaF/F')
xlabel('Time (frames @30Hz)')


title('1% 2-butanone for 5 sec. Pan KC line')
fig = formatFig(fig, true);

save_figure(fig, [fileRoot vidName ' df-f']);
























warning('on','all')



%% Load tiff stack









































