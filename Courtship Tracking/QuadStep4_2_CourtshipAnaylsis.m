

%% Identify periods of presumed courtship
clearvars('-except',initial_vars{:})
vidBase = getDataPath(2,2);
% 1) Find periods  where the body distance is less than 5 mm & the speed is
% higher than something and the speed correlation is above 0 (aka can't
% just be flies sitting together at the food)
exp = 1;
closepoints = [];
IFD = [];
for trial = 1:num.trial(exp)
    y = (data(exp).data(trial).occupancy.IFD/pix2mm)';
    IFD = autoCat(IFD, y<=5,false);
    closepoints(trial) = (sum(y<=5)/length(y))*100;
end

% pull out the speed data
for trial = 1:num.trial(exp)
    speed =[];
    y = data(exp).data(trial).speed.raw;
    for i = 3:size(y,2)
        
        % 1) find locations with 3+ tracked points on the image...
        
        loc = sum(~isnan(y),2)>3; % locations with more than 3 flies tracked
        start3 = find(diff(loc)==1);
        stop3 = find(diff(loc)==-1);
        overcount = min([length(start3) length(stop3)]);
        
        % 2) plot image of the arena and the plotted points -- click on the fake one if possible
        vidDir = [vidBase, data.T.Date{trial}, '/' data.T.ExperimentID{trial}];
        for ii = 1:overcount
            % find video frames (1 before to 1 after) for the overcount
            frame = start3(ii);
            vid = data(exp).data(trial).data.T.vidNums(frame);
            vidframes = data(exp).data(trial).data.T.vidFrame(frame);
            
            % find the X-Y locs for flies on the frame
            X = data(exp).data(trial).data.x_loc(frame,:);
            Y = data(exp).data(trial).data.y_loc(frame,:);
             
            movieInfo = VideoReader([vidDir, '_' num2str(vid) '.avi']);
            img = read(movieInfo, vidframes);

            % image cropping: 
            centre = data(exp).data(trial).data.centre;
            r = data(exp).data(trial).data.r;
            % plot image
            fig = getfig;
            for rr = 1:3
                subplot(1,3,rr)
                    img = read(movieInfo, vidframes+(rr-2));
                    imshow(img)
                    xlim([centre(1)-r,centre(1)+r])
                    ylim([centre(2)-r,centre(2)+r])
                    hold on
                    scatter(X,Y,50,Color('gold'))
                    % draw 
            end

                 


end


speed(speed>=15) = nan;
time = grouped(exp).time;


r = 2; c = 1;
fig = getfig;
subplot(r,c,1)
plot(IFD)
subplot(r,c,2)
plot(speed)

% calculate the speed of the flies centroid (if they are under a specific distance)



%% Tracking cleaning

% may want to clean up the tracks manually before screening? 

temp = load("S:\Evyn\DATA\Trial Data\02.04.2023_linear_recovery_F_caviar_C\linear_recovery_F_caviar speed data.mat");

% speedTracks structure contains the x-y positions of the flies as well as
% the speed for each set of points 








































































































