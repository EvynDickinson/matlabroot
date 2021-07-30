

%% Load TIFF image stack to analyze
clear
warning('off','all') % Suppress all the tiff warnings
folderdate = '10.8.20';
fileRoot = ['E:\Evyn\' folderdate '\'];

% select specific trial:
fileList = dir([fileRoot '*.tif*']);
for ii = 1:length(fileList)
    TiffFiles{ii} = fileList(ii).name;
end
idx = listdlg('PromptString', 'Select Video', 'ListString', TiffFiles);
num_vids = length(idx);

% import video for all selected files
for ii = 1:num_vids
    vidName = TiffFiles{idx(ii)};
    fprintf(['\nLoading: ' vidName '\n'])
    % load stack:
    tstack = Tiff(fullfile(fileRoot,vidName),'r');
    [I,J] = size(tstack.read()); % get the pixel sizes for frame
    nframes = length(imfinfo(fullfile(fileRoot,vidName))); % get the total number of frames
    IM = zeros(I,J,nframes);
    IM(:,:,1)  = tstack.read();
    disp('frames...')
    for n = 2:nframes
        tstack.nextDirectory()
        IM(:,:,n) = tstack.read();
        if rem(n,50)==0 %indicate every 50th frame
            fprintf(['  ' num2str(n)])
        end
    end
    
    % add full video file to vid structure
    vid(ii).IM = IM;
    vid(ii).vidName = vidName;
end
disp('done!')

% ------------- FIGURE--------------%
% Overlay multiple trials single basic figure:
% Basic parameters:
shock_delay = 10-2;
odor_delay = 10-2;
stim_duration = 1;
fps = 100/3;
odor_on = odor_delay;
odor_off = (odor_delay+stim_duration);
% Basic parameters:
controlROI = 1:7; %1-sec buffer for start and before odor
controlStart = round(controlROI(1)*fps);
controlStop  = round(controlROI(end)*fps);

fig = getfig('',1);
subplot(2,1,1); hold on
for n = 1:num_vids
    IM = vid(n).IM;
    vid_length = size(IM,3)/fps;
    nframes = size(IM,3);

    %plot avg intensity of entire image over time:
    time = linspace(1,vid_length,nframes);
    a = mean(mean(IM,1),2);
    y(:,n) = squeeze(a);
end
plot(time,y, 'linewidth', 1, 'color', Color('purple'))
y1 = rangeLine(fig);
plot([odor_delay, odor_delay+stim_duration],[y1,y1], 'color', Color('blue'), 'linewidth', 3)

% Plot df/f for trials: * IN PROGRESS *
subplot(2,1,2); hold on
for n = 1:num_vids
    IM = vid(n).IM;
    %cutoff final frame (always drops to zero)
    IM(:,:,end-3:end) = [];
    vid_length = size(IM,3)/fps;
    nframes = size(IM,3);
    time = linspace(1,vid_length,nframes);
    
    % Background florescence
    ROI = controlStart:controlStop;
    controlIM = IM(:,:,ROI);
    F0 = mean(mean(mean(controlIM,3)));

    % Find change in florescence (full frame only)
    Ft = squeeze(mean(mean(IM,1),2));
    df_f(:,n) = ((Ft-F0)/F0)*100;
end

% Plot the df/f
plot(time,df_f,'color', Color('grey'), 'linewidth', 1, 'linestyle',':')
plot(time,mean(df_f,2),'color', Color('teal'), 'linewidth', 2)

%labels etc:
v_line([odor_on, odor_off],'Navy',1)
ylabel('\DeltaF/F')
xlabel('Time (s)')

% find region with most change within a given time ROI:
t1 = controlStop;
% t2 = 800;
IM = vid(1).IM;
selIM = IM(:,:,t1:end);
IMrange = range(selIM,3);

fig = getfig('',1);
imagesc(IMrange)
colorbar


%%
% save the data for quick comparisons with the other videos:

naive.vid = vid;
naive.df_f = data;
save([fileRoot, '2_butanone_NAIVE'], 'naive')
save([fileRoot, '3_octanol_NAIVE'], 'naive')

training.vid = vid;
training.df_f = data;
save([fileRoot, '2_butanone_CONDITIONING'], 'training')
save([fileRoot, '3_octanol_CONDITIONING'], 'training')

post.vid = vid;
post.df_f = data;
save([fileRoot, '2_butanone_CONDITIONING'], 'post')
save([fileRoot, '3_octanol_CONDITIONING'], 'post')

