
clearvars('-except',initial_vars{:})
close all
%  Data Structure: data(video).tracks(frame number, body nodes, X Y coordinates, indivual fly tracks)

% pull and select data
vid = 1;
roi = 1:100;
flynums = [6,10,16,22,33,41,50];
test = squeeze(data(vid).tracks(:,1,:,flynums));

% pull image
movieInfo = VideoReader([figDir,'\',expName,'_',num2str(vid),'.avi']); 
% PLOT IMAGE
fig = figure; set(fig, 'pos', [2040 499 814 733],'color', 'k'); % X-off, Y-off, width, height
nanLoc = [nan,nan];
for ii = roi
    %pull current data
    img = read(movieInfo,ii);
    plotdata = squeeze(data(vid).tracks(ii,1,:,:));
    
    imshow(img)
    axis tight square
    hold on
    
    scatter(plotdata(1,:),plotdata(2,:), 15,'w','filled')
    % plot the history of 5 random flies...
    for fly = 1:length(flynums)
        x = test(ii,1,fly);
        y = test(ii,2,fly);
        plot(test(roi(1):ii,1,fly),test(roi(1):ii,2,fly),'Color', Color('yellow'),'LineWidth',0.5)
        scatter(x,y,15,'r','filled')
        if isnan(x) 
            % save the nan locations so that all future slides can also include the points
            nanLoc = [nanLoc; test(ii-1,:,fly)];
        end
        scatter(nanLoc(:,1),nanLoc(:,2),20,Color('orange')) 
        
    end

%     % throw an orange dot for any nan locations
%     if any(plotdata())
%         
%         test(roi(1):ii,1,fly)
%     % TODO...
%     end
    pause(0.05)
end


%% 

clearvars('-except',initial_vars{:})
close all

% pull and select data
vid = 1;
roi = 1:1797;
flynums = [6,10,16,22,33,41,50];
test = squeeze(data(vid).tracks(:,1,:,flynums));

% pull video information
movieInfo = VideoReader([figDir,'\',expName,'_',num2str(vid),'.avi']); 
buffer_frame = read(movieInfo,1);
buffer_frame(:,:,:) = 1;

% Set up video parameters
fig = figure; set(fig, 'pos', [2040 499 814 733],'color', 'k'); % X-off, Y-off, width, height
nanLoc = [nan,nan];

    hold on
    v = VideoWriter(['G:\My Drive\Jeanne Lab\DATA\06.22.2022\analysis\Tracking Video2'], 'Uncompressed AVI');
    v.FrameRate = 30;
    open(v);
    imshow(buffer_frame); axis tight square
    currAxes.Visible = 'off';
    % Collect image info for movie file
    f = getframe(fig);
    writeVideo(v, f)  
    clf('reset')

    % Draw video frames...
    for ii = roi
        set(fig, 'color','k'); 
        %pull current data
        img = read(movieInfo,ii);
        plotdata = squeeze(data(vid).tracks(ii,1,:,:));
        imshow(img)
        axis tight square
        hold on
        scatter(plotdata(1,:),plotdata(2,:), 15,'w','filled')
        % plot the history of 5 random flies...
        for fly = 1:length(flynums)
            x = test(ii,1,fly);
            y = test(ii,2,fly);
            plot(test(roi(1):ii,1,fly),test(roi(1):ii,2,fly),'Color', Color('yellow'),'LineWidth',0.5)
            scatter(x,y,15,'r','filled')
            if isnan(x) 
                % save the nan locations so that all future slides can also include the points
                nanLoc = [nanLoc; test(ii-1,:,fly)];
            end
            scatter(nanLoc(:,1),nanLoc(:,2),20,Color('orange'))     
        end
    


        currAxes.Visible = 'off';
        % Collect image info for movie file
        f = getframe(fig);
        writeVideo(v, f)  
        clf('reset')
    end

    
close(v)
close all
fprintf('\n Video Saved! \n')  





    



  
   
        










        
















