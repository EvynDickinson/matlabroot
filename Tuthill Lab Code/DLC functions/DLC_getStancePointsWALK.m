

function [cloud_1, cloud_2, combinedData] = DLC_getStancePointsWALK(data,figDisplay)
% [cloud_1, cloud_2, combinedData]  = DLC_getStancePointsWALK(data);
% 
% Extract the stance frames from a walking fly trial to use
% for later fitting a sphere. Requires all six legs input.
% Input: 
% data(leg).input
% figDisplay - T|F choice to display figures
%              default 'false' *off*
% Output:
% nx3 data point matrix (x,y,z) coordinates
% cloud - from euclidian distance
% combinedData - cloud filtered by angle data
% 
% ES Dickinson,
% University of Washington, 2020



% ----------------------------------------------------------------------------
% Plot the Incoming data:
if nargin == 1
    % option: turn off figures
    set(0,'DefaultFigureVisible','off');
else
    if figDisplay == true
        set(0,'DefaultFigureVisible','on');
    else
        set(0,'DefaultFigureVisible','off');
    end
end
    
% PLOT THE DATA FOR EACH TARSUS IN SPACE (2 figures)
DLC_plot_tarsus(data); 


% ----------------------------------------------------------------------------
% Find the orientation of the fly in cartesian space
% determine the orientation of the points compared to 'zero'
for leg = 1:6
    input = data(leg).input;
    e_dist = pdist2([0,0,0],input)';
    temp(leg) = nanmean(e_dist);
end
front = mean([temp(1), temp(4)]); % avg distance for front legs
rear = mean([temp(3), temp(6)]); % avg distance for rear legs
% NOTE: if the front euclidian distance is greater, then e_dist MIN = start of stance 
% if the rear euclidian distance is greater, then e_dist MAX = start of stance
if front>rear %MIN -> start stance
    stance_val = 1;
else
    stance_val = -1;
end

% fprintf('\n Fly orientation: '); disp(stance_val)

% ----------------------------------------------------------------------------
% Find the locations of the position peaks for stance and swing starts
temp = [];
idx = 1;
for leg = 1:6
    stance = [];

%     input = data(leg).ind;
    input = data(leg).input;
    for ii = 1:3
        temp(:,ii) = smooth(input(:,ii),7);
    end
    % get euclidian distances from (0,0,0):
    e_dist = pdist2([0,0,0],temp)';
    try 
        max_loc = round(findPeaks('-x', e_dist, '-extrema', (stance_val),'-includeEndpoints', false));
        min_loc = round(findPeaks('-x', e_dist, '-extrema', -(stance_val),'-includeEndpoints', false));
    catch
        h = warndlg('No peaks detected with findPeaks, check data manually');
        uiwait(h)
        continue;
    end
    nn = min([length(min_loc), length(max_loc)]);
    % find the start of stance periods:
    for ii = 1:nn
        stance(ii,1) = max_loc(ii);
        try
            a = min_loc(min_loc>stance(ii,1));
            stance(ii,2) = a(1);
        catch
            stance(ii,:) = [];
        end
    end
    nsteps = size(stance,1);
    if nsteps == 0
        STEP(leg).qual = false;
        % remove this leg from the points
        STEP(leg).stance = nan;
        STEP(leg).nsteps = nan;
        STEP(leg).e_dist = nan;
        STEP(leg).input = nan;
        STEP(leg).ROI = nan;
        STEP(leg).dur = nan;
        continue;
    end
    
% % ------------PLOT-----------------        
%     fig = getfig('',1);
%     subplot(2,1,1) % raw and smoothed traces
%     plot3(input(:,1),input(:,2),input(:,3))
%     hold on
%     plot3(temp(:,1),temp(:,2),temp(:,3))
%     grid on
%     title(['Leg ' leg_labels{leg}])
%     legend({'Raw trace', 'smoothed'})
%     subplot(2,1,2) % euclidian distance with peaks
%     plot(e_dist)
%     vline(min_loc, 'g')
%     vline(max_loc, 'm')
%     legend({'euclidian distance'})
%     xlabel('time (fps)')
% % ---------------------------------    
       
    % Check how close the 'peak' points are in physical space:
    % in theory, the euclidian distance between the start points for swing
    % or stance should be very similar *for good walking*
    % QUALITY CONTROL
    % find dist between each set of points:  
    A = input(stance(:,1),:); % stance start
    B = input(stance(:,2),:); %stance end
    % find distance between all points for start|end of stance
    strt_dist = pdist2(A,A);
    strt_dist(strt_dist==0) = nan;
    end_dist = pdist2(B,B);
    end_dist(end_dist==0) = nan;
    %avg standard deviation between points
    strt_avg = mean(nanstd(strt_dist));
    end_avg = mean(nanstd(end_dist));
    % implement a series of quality control (10%)
    if (strt_avg/end_avg)>10 || (strt_avg/end_avg)<0.10 % 10% of other value
        % to find the appropriate region of points to select: 
        %find total length between first and last step, divide by num of steps
        %then take only 40% of the avg steps...
        dur = round(((stance(end,2)-stance(1,1))/nsteps)*0.4);
        offset = round(dur*0.20);
        % use values from points that are clumping correctly:
        %start of stance is the good landmark
        if strt_avg<end_avg 
            ROI = [stance(:,1)+offset,stance(:,1)+offset+dur];
        %end of stance is the good landmark
        else
            ROI = [stance(:,2)-offset-dur,stance(:,2)-offset];
        end
    % good quality point selection for both start and end of stance
    else
        dur = round((stance(:,2)-stance(:,1))*0.25); %middle 50%
        ROI = [stance(:,1)+dur,stance(:,2)-dur];
    end

    STEP(leg).stance = stance;
    STEP(leg).nsteps = nsteps;
    STEP(leg).e_dist = e_dist;
    STEP(leg).input = input;
    STEP(leg).ROI = ROI;
    STEP(leg).dur = dur;
    
    clear A B a
    
    % SECOND METHOD OF POINT ISOLATION:
    % run angle calc for all trials:
    roi_sz = 10;
    ofst = round(roi_sz/2);
    ntot = length(temp)-(2*roi_sz); %number of comparisions
    t_in = roi_sz:length(temp)-roi_sz; % point locs to check
    roi = [t_in',t_in'+roi_sz]; %list of ROIs
    TP = t_in'+ofst; %list of center test points indexes
    p1 = temp(roi(:,1),:);
    p2 = temp(TP,:);
    p3 = temp(roi(:,2),:);
    % find the vectors:
    U = (p1-p2);
    V = (p3-p2);
    for pp = 1:length(U)
        u = U(pp,:);
        v = V(pp,:);
        CosTheta(pp) = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
        ThetaInDegrees(pp,1) = real(acosd(CosTheta(pp)));
    end
    % find peaks in the angle data:
    min_angles = smooth(ThetaInDegrees,roi_sz);
    angle_peak = peakseek(min_angles,round(length(min_angles)/12)); %
    angle_peak(angle_peak<10) = []; % remove start edge effect
    angle_peak = angle_peak+ofst;

    % - determine stance start vs end from angle peaks -
    %piggyback on the euclidian distance measures   
    group1 = angle_peak(1:2:end);
    group2 = angle_peak(2:2:end);

    %avg stance strt&end position in space
    if  STEP(leg).nsteps==1
        strt = input(stance(:,1),:);
        fin = input(stance(:,2),:);
    else
        strt = mean(input(stance(:,1),:));
        fin = mean(input(stance(:,2),:));
    end
    % distance from each cluster to the strt&end values:
    g1_2_strt = mean(pdist2(strt,input(group1,:)));
    g1_2_end = mean(pdist2(fin,input(group1,:)));
    g2_2_strt = mean(pdist2(strt,input(group2,:)));
    g2_2_end = mean(pdist2(fin,input(group2,:)));
    % group 1 is START of stance
    if g1_2_strt < g1_2_end && g2_2_end < g2_2_strt
        stance_start = group1;
        stance_end = group2;
    else
        stance_start = group2;
        stance_end = group1;
    end
    % find paired steps:
    nn = min([length(stance_start), length(stance_end)]);
    for pp = 1:nn
        Astance(pp,1) = stance_start(pp);
        try
            a = stance_end(stance_end>Astance(pp,1));
            Astance(pp,2) = a(1);
        catch
            Astance(pp,:) = [];
        end
    end
    Ansteps = size(Astance,1);
    Adur = round((Astance(:,2)-Astance(:,1))*0.1); %middle 90%
    AROI = [Astance(:,1)+Adur,Astance(:,2)-Adur];
    
    
%     % --- view figures of the angle data ----
%     set(0,'DefaultFigureVisible','on');
%     fig = getfig('',1);
%     subplot(2,1,1)
%     plot(ofst:length(ThetaInDegrees)+ofst-1, smooth(ThetaInDegrees,roi_sz))
%     vline(angle_peak, 'r')
%     title({'Angle'})
%     subplot(2,1,2) % euclidian distance with peaks
%     plot(e_dist)
%     vline(min_loc, 'g')
%     vline(max_loc, 'm')
%     title({'euclidian distance'})
%     xlabel('time (fps)')
% 
%     LW =1;
%     fig = getfig('',1);
%     % subplot 1 : raw trace
%     plot3(input(:,1),input(:,2),input(:,3))
%     hold on; grid on
%     test = round(angle_peak);%(2:end-1)
%     scatter3(input(test,1), input(test,2), input(test,3), 100, 'r', 'filled')
%     scatter3(input(min_loc,1), input(min_loc,2), input(min_loc,3), 100, 'g','filled')
%     scatter3(input(max_loc,1), input(max_loc,2), input(max_loc,3), 100, 'm','filled')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     legend({'tarsus position', 'angle selection', 'distance selection', 'distance selection'})
%     title('Selection points via euclidian dist or \Deltaangle')
%     % ------------------------------------------

    Ang_STEP(leg).stance = Astance;
    Ang_STEP(leg).nsteps = Ansteps;
    Ang_STEP(leg).theta = ThetaInDegrees;
    Ang_STEP(leg).input = input;
    Ang_STEP(leg).ROI = AROI;
    Ang_STEP(leg).dur = Adur;
    
    
    STEP(leg).qual = true;
    if isempty(STEP(leg).stance)
        STEP(leg).qual = false;
        % remove this leg from the points
        STEP(leg).stance = nan;
        STEP(leg).nsteps = nan;
        STEP(leg).e_dist = nan;
        STEP(leg).input = nan;
        STEP(leg).ROI = nan;
        STEP(leg).dur = nan;
    end    
end

save('StepData', 'STEP','Ang_STEP');


% ------------PLOT-----------------    
fig = DLC_plot_ROIForFit(STEP);
% ---------------------------------

% set(0,'DefaultFigureVisible','on');
% % ------------PLOT-----------------    
% fig = DLC_plot_ROIForFit(Ang_STEP);
% % ---------------------------------
% set(0,'DefaultFigureVisible','off');


% ----------------------------------------------------------------------------
%%  Group all the 'stride' regions for each leg from both methods:
for leg = 1:6
    if STEP(leg).qual == false 
        [XZ(leg),YZ(leg),ZRange(leg)] = deal(nan);
        continue;
    end
    input = STEP(leg).input;
    stance = STEP(leg).stance;
    Ang_stance = Ang_STEP(leg).stance;
    % find overlap?
    %stance_1 regions
    ROI_1 = [];
    for n = 1:STEP(leg).nsteps
        a = stance(n,1):stance(n,2);
        ROI_1 = [ROI_1, a];
    end
     %stance_2 regions
    ROI_2 = [];
    for n = 1:Ang_STEP(leg).nsteps
        a = Ang_stance(n,1):Ang_stance(n,2);
        ROI_2 = [ROI_2, a];
    end
    % find overlap
    ROI = ROI_1(ismember(ROI_1,ROI_2));
    if isempty(ROI) % skip leg if no overlapping regions found
        [XZ(leg),YZ(leg),ZRange(leg)] = deal(nan);
        continue;
    end
    all_data(leg).input = input;
    all_data(leg).ROI = ROI;
    temp = input(ROI,:); 
    all_data(leg).FitPoints = temp;
    
    % find XZ & YZ slopes
    all_data(leg).XZ = mean(diff(temp(:,3))./diff(temp(:,1))); 
    all_data(leg).YZ = mean(diff(temp(:,3))./diff(temp(:,2)));
    % find the avg X,Y,Z
    XZ(leg) = all_data(leg).XZ;
    YZ(leg) = all_data(leg).YZ;
    ZRange(leg) = range(temp(:,3));
    
end

% % ----- Figure: 
% fig = getfig('',1);
% for leg = 1:6
%     input = all_data(leg).input;
%     plot3(input(:,1),input(:,2),input(:,3), 'color', 'k', 'linewidth', 1);
%     hold on
%     input = all_data(leg).FitPoints;
%     scatter3(input(:,1),input(:,2),input(:,3), 100, 'r', 'filled')
% end
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid on
% % -----------------
all_exclude = (ZRange > 0.1) & (abs(XZ) > 1) & (abs(YZ) > 1);
combinedData = [];
for leg = 1:6
    if all_exclude(leg) == true; continue; end
    try
    combinedData = [combinedData; all_data(leg).FitPoints];
    catch
    end
end

% % TODO
% set(0,'DefaultFigureVisible','on');
% % % Leg quality control image:
% fig = getfig('',1);
% for leg = 1:6
%     if all_exclude(leg) == true
%         C = 'red';
%     else
%         C = 'grey';
%     end
%     scatter3(all_data(leg).FitPoints(:,1),...
%              all_data(leg).FitPoints(:,2),...
%              all_data(leg).FitPoints(:,3), 50, Color(C))
%     hold on
% end

%% CLOUD 1 quality check - from euclidian distance

% Check for ROI leg outliers
% Find the avg ROI region point per leg (cartesian)
for leg = 1:6 %number of legs
    if STEP(leg).qual==false
        centre(leg).avg = nan;
        tester(leg,:) = [nan, nan, nan];
        continue
    end
    temp = []; centre(leg).data = [];
    input = STEP(leg).input;
    ROI = STEP(leg).ROI;
    for n = 1:STEP(leg).nsteps
        roi = ROI(n,1):ROI(n,2);
        temp = [input(roi,1), input(roi,2), input(roi,3)];
        centre(leg).data = [centre(leg).data; temp];
    end
    % find the avg X,Y,Z
    centre(leg).avg = nanmean(centre(leg).data);
    tester(leg,:) = centre(leg).avg;
end

% find the diff between the various coordinates and the center point:
std_dev = nanstd(tester); % std of each axes
av = nanmean(tester); % avg value of each axes
min_val = (av-std_dev);
max_val = (av+std_dev);

% how does each leg fall in comp to other legs?

for leg = 1:6
    if STEP(leg).qual == false
        rr(leg,:) = nan(1,3);
        check(leg,:) = false(1,3);
        continue
    end
    check(leg,:) = tester(leg,:)>max_val;
    rr(leg,:) = range(centre(leg).data);
end
check2 = rr>median(rr);
exclude = sum(check2&check,2)>0;

% set(0,'DefaultFigureVisible','on');
% % Leg quality control image:
% figure;
% for leg = 1:6
%     if exclude(leg) == true
%         C = 'red';
%     else
%         C = 'grey';
%     end
%     scatter3(centre(leg).data(:,1),centre(leg).data(:,2),centre(leg).data(:,3), 50, Color(C))
%     hold on
%     scatter3(centre(leg).avg(1),centre(leg).avg(2),centre(leg).avg(3), 50, 'k', 'filled')
% end
% av_point = mean(all_vals);
% scatter3(av_point(1), av_point(2),av_point(3), 300, 'p', 'filled')

% Organize all data points into a simple matrix
m_data = [nan,nan,nan];
% exc = false(6,1);
for ii = 1:6 %number of legs
    % Add an optional removal of the legs that have a 'Z' out of range
    if exclude(ii)==true || STEP(ii).qual==false 
        fprintf(['\n Skipped leg ' num2str(ii) '\n'])
        continue;
    end
    stance = STEP(ii).stance;
    input = STEP(ii).input;
    ROI = STEP(ii).ROI;
    newData = [];
    for n = 1:STEP(ii).nsteps
        roi = ROI(n,1):ROI(n,2);
        added = [input(roi,1), input(roi,2), input(roi,3)];
        newData = [newData; added];
    end
    m_data = [m_data;newData];
end
m_data(1,:) = [];
cloud_1 = m_data;


%% CLOUD 2 quality check - from angle data %  !!!

% Check for ROI leg outliers
% Find each step for all legs and then pull stats on those:
stepData = []; idx = 1;
for leg = 1:6 %number of legs
    temp = []; centre(leg).data = [];
    try
        input = Ang_STEP(leg).input;
    catch
       continue; 
    end
    for n = 1:Ang_STEP(leg).nsteps
        roi = Ang_STEP(leg).stance(n,1):Ang_STEP(leg).stance(n,2);
        temp = [input(roi,1), input(roi,2), input(roi,3)];
        stepData(idx).data = temp;
        XZ(idx,1) = mean(diff(temp(:,3))./diff(temp(:,1))); 
        YZ(idx,1) = mean(diff(temp(:,3))./diff(temp(:,2)));
        ZRange(idx,1) = range(temp(:,3));    
        idx = idx+1;
    end
end
    
exclude_2 = sum([(ZRange > 0.1), (abs(XZ) > 1), (abs(YZ) > 1)],2)>1;
% % ---- Figure of exclusions ----
% set(0,'DefaultFigureVisible','on');
% % Leg quality control image:
% figure;
% for idx = 1:length(XZ)
%     if exclude_2(idx) == true
%         C = 'red';
%     else
%         C = 'grey';
%     end
%     scatter3(stepData(idx).data(:,1),stepData(idx).data(:,2),stepData(idx).data(:,3), 50, Color(C))
%     hold on
% end
% % ----------------------------

% Organize all data points into a simple matrix
c_data = [];
% exc = false(6,1);
for idx = 1:length(XZ)
    if exclude_2(idx) == true
        continue;
    end
    c_data = [c_data; stepData(idx).data];
end
cloud_2 = c_data;


end


