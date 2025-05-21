%   
% download_raw_data_from_server

clear
clc

%% --------------------- Select Date & Experiment to Process ------------------------------
%load excel file:
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;
switch questdlg('Use Excel sheet to select experiment?')
    case 'Yes'
        loc1 = cellfun(@isnan,excelfile(2:end,Excel.step1));
        % loc1 = cellfun(@isnan,excelfile(2:end,Excel.processed));
        loc2 = cellfun(@ischar,excelfile(2:end,Excel.tracked));
        loc = loc1 & loc2;
        rownums = find(loc)+1;
        if isempty(rownums)
            warndlg('No available experiments') 
            return
        end
        eligible_files = excelfile([false;loc],[Excel.date, Excel.arena, Excel.expID]); 
        FileNames = join(eligible_files);
        fileIdx = listdlg('ListString', FileNames,'ListSize',[250,450]);   
        if isempty(fileIdx) % abort if there are no files selected (e.g. cancel button or esc pressed)
            return
        end 
        % get file info:
        baseFolder = getDataPath(2,0); % select raw data from a user defined location
        % baseFolder = getCloudPath;
        folder = eligible_files{fileIdx,1};
        expName = eligible_files{fileIdx,3};
        clear loc1 loc2 loc eligible_files FileNames rownums fileIdx
    case 'No'
        %get base folder pathway
        baseFolder = getDataPath(2,0);
        folder = selectFolder(baseFolder);
        folder = folder{:};
        % [baseFolder, folder] = getCloudPath(2);  %old version
%         % Select the complete experiments to process
        list_dirs = dir([baseFolder folder, '/*dataMat.mat']); %only matlab files
        list_dirs = {list_dirs(:).name};
        expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
        expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
        expName = expName(1:end-1); clear expNames 
    case 'Cancel'
        return
end

% Saving and Loading Directories:
analysisDir = [baseFolder folder '/analysis/'];
if ~isfolder(analysisDir)
    mkdir(analysisDir) 
end

% expPDF = [analysisDir folder expName ' summary.pdf'];
XLrow = find(strcmpi(excelfile(:,Excel.date), folder) & ...
                      strcmpi(excelfile(:,Excel.expID), expName)); 

% Load relevant data files (.mat, .csv, .h5)
warning off 
expData = load([baseFolder folder '/' expName ' dataMat.mat']);
tempLog = readmatrix([baseFolder folder '/' expName '_RampLog']);
nvids = expData.parameters.numVids;

%load tracking predictions     
data = struct;
disp('Loading data files...')

for vid = 1:nvids
   vidBase = [baseFolder folder '/' expName '_' num2str(vid)]; 
   fp_endings = {'.h5', '.avi.predictions.slp.h5', '.avi.predictions.analysis.h5'};
   dataIn = false; % switch for logging data
   for suf = 1:length(fp_endings)
       filePath = [vidBase fp_endings{suf}];
       if exist(filePath,'file')
           data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
           data(vid).tracks = h5read(filePath,'/tracks');
           dataIn = true; %log successful data load
            continue
       end
   end
   if ~dataIn
        disp(vidBase)
        h = warndlg('Warning: file not found');
        uiwait(h)
        return
        % If file is corrupt or permanently missing...fill with blank NaN
        % data
%         data(vid).occupancy_matrix = [];
%         data(vid).tracks = [];
   end
end; clear filePath fp_endings suf vidBase

initial_vars = who; initial_vars{end+1} = 'initial_vars';
fprintf('\nNext\n')

%% Load arena outlines and copy data to analysis folder
% Load the well location data 
well_loc_file = [baseFolder folder '/' expName ' well_locations.mat'];
try load(well_loc_file);
catch % if the location data is missing -- outline the wells here:
    movieInfo = VideoReader([baseFolder folder '/' expName '_1.avi']); %read in video
    demoImg = imadjust(rgb2gray(read(movieInfo,1)), [55/255,150/255]); % pull out a 
    getArenaOutlines(demoImg, [baseFolder folder '/' expName]) % actually outline the arenas and wells
    load(well_loc_file);
end
copyfile(well_loc_file, [analysisDir expName ' well_locations.mat']);

% Load the well location data 
well_coor_file = [baseFolder folder '/' expName ' arena coordinates.mat'];
load(well_coor_file);
copyfile(well_coor_file, [analysisDir expName ' arena coordinates.mat']);% copy the files to the new analysis folder

% Tidy variables / environment 
initial_vars = [initial_vars; 'arenaData'; 'arenaIdx'; 'r'; 'wellcenters'; 'demoImg';'radii'];
clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Get the number of flies in each arena 
% Number of flies: 
nframes = 3;
nflies = nan(1,4); 
for arena = 1:4
    arenaData(arena).nflies = excelfile{XLrow(arena),Excel.numflies};
    nflies(arena) = arenaData(arena).nflies;
end
movieInfo = VideoReader([baseFolder folder '/' expName '_1.avi']); %read in video
ii = randi(size(data(1).occupancy_matrix,2),[3,1]); %random selection of frames to count fly number
demoImg = rgb2gray(read(movieInfo,1));
demoImg = imadjust(demoImg, [55/255,150/255]);
if any(isnan(nflies)) 
    nflies = [];
    % manual count of flies
    fprintf('\nCount the number of flies in the picture by clicking them\n then hit ENTER\n')
    T = true;
   while T % get the number of flies in the arena 
       % % find the X & Y location for SLEAP tracked points for this video:
       %  headData = squeeze(data(1).tracks(:,1,:,:)); %this is ALL points in all arenas for each frame
       %  x_loc = squeeze(headData(:,1,:));
       %  y_loc = squeeze(headData(:,2,:));
       %  calc_flies = nan(1,4);
       %  for arena = 1:4
       %          centre = arenaData(arena).centre;
       %          % find points within each arena sphere
       %          fly_loc = (((x_loc-centre(1)).^2 + (y_loc-centre(2)).^2).^0.5)<=r;
       %          calc_flies(arena) = median(sum(fly_loc,2));
       %          %save the locations for the randomized frame selection
       %          X = x_loc(ii,:);   Y = y_loc(ii,:); LOC = fly_loc(ii,:);
       %          X(LOC)
       %          sum(LOC,2)
       %  end 
        
        fprintf('Number of flies counted \n    A     B     C     D\n    ---    ---    ---    ---\n')
        for jj = 1:nframes
            demoImg = rgb2gray(read(movieInfo,ii(jj))); %display the arena frame image
            demoImg = imadjust(demoImg, [55/255,150/255]);
            PT(jj).frame = readPoints(demoImg);
            for arena = 1:4
                centre = arenaData(arena).centre;
                X = PT(jj).frame(1,:);
                Y = PT(jj).frame(2,:);
                % find points within each arena sphere
                loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)<=r;
                nflies(arena,jj) = sum(loc);
            end
            disp(nflies(:,jj)')
        end
        disp('Number of flies:')
        A = nflies(1,:)'; B = nflies(2,:)'; C = nflies(3,:)'; D = nflies(4,:)';
        numFlies = [A,B,C,D];
        disp(table(A,B,C,D))
        
        % Skip Frame Speed Hack
        % if there is a row of zeros (skipped tracking frame) remove the empty frame
        skipFrameloc = sum(numFlies==0,2)==4; %skipped counting frame locations
        mismatchLoc = diff(numFlies)==0; % across frame mismatched numbers
        % if we skipped labeling a frame but the first two counts match up,
        % then auto add the number from the first two counts to the list
        if any(skipFrameloc) && all(mismatchLoc(1,:)) 
                numFlies(skipFrameloc,:) = numFlies(1,:);
        end
        % Resume normal miscounted fly error catching protocol
        count_match = (diff(numFlies,(nframes-1),1)==0);
        if all(count_match)
            nflies = [];
            for arena = 1:4
                nflies(arena) = median(numFlies(:,arena)); 
                arenaData(arena).nflies = nflies(arena);
            end
            T = false;
        else
%             disp('    A     B     C     D')
%             disp(nflies')
            switch questdlg('Nonmatching fly counts, manually select number?:')
                case 'Yes'
                    for arena = 1:4
                        arenaNums{arena} = num2str(median(numFlies(:,arena)));
                        if ~count_match(arena)
                            arenaNums{arena} = 'NaN';
                        end
                    end
                    prompt = {'Arena A','Arena B', 'Arena C', 'Arena D'}; dlgtitle = 'Input';
                    answer = inputdlg(prompt,dlgtitle,[1 25],arenaNums);
                    nflies = [];
                    for arena = 1:4
                        arenaData(arena).nflies = str2double(answer{arena});
                        nflies(arena) = arenaData(arena).nflies;
                    end
                    T = false;
                case 'No'
                    T = true;
                case 'Cancel'
                    return
            end
        end
    end
    % write number of flies into the excel sheet
    try
        for arena = 1:4
            writematrix(nflies(arena),xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.numflies) num2str(XLrow(arena))]);
        end
    catch
        h = warndlg('Close Experiment Summary excel file and then close this warning box');
        uiwait(h)
        for arena = 1:4
            writematrix(nflies(arena),xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.numflies) num2str(XLrow(arena))]);
        end
    end
end
disp('Number of flies:')
disp('    A     B     C     D')
disp(nflies)

%% Start data table that concatenates across all videos:

[vidNums, vidFrame, frame, temperature, tempWork] = deal([]);
currFrame = 0;
for vid = 1:nvids
    if ~isempty(data(vid).occupancy_matrix) %add blanks for corrupt videos
        nframes = size(data(vid).occupancy_matrix,2);
    else
        nframes = median(occupancy.frameROI(:,2)-occupancy.frameROI(:,1));
    end
    % find number of frames per vid
    occupancy.frameROI(vid,1) = currFrame+1;
    currFrame = nframes + currFrame;
    occupancy.frameROI(vid,2) = currFrame;
    % video numbers
    vidNums = [vidNums; vid*ones(nframes,1)];
    % frame num in video
    vidFrame = [vidFrame; (1:nframes)'];
    % total frame count
    frame = [frame; (occupancy.frameROI(vid,1):occupancy.frameROI(vid,2))'];
    % temperature log
    logROI(1) = find(tempLog(:,1)==expData.tempLogStart(vid,3));
    logROI(2) = find(tempLog(:,1)==expData.tempLogEnd(vid,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    x = round(linspace(1, nframes, length(tempCourse)));
    fullTempList = interp1(x,tempCourse,1:nframes,'spline');   
    data(vid).tempLog = fullTempList;
    temperature = [temperature; fullTempList']; 
    % temp plate work log
    workCourse = tempLog(logROI(1):logROI(2),4);
    x = round(linspace(1, nframes, length(workCourse)));
    fullWorkList = interp1(x,workCourse,1:nframes,'spline');   
    tempWork = [tempWork; fullWorkList']; 
end
occupancy.temp = temperature;
% Time count
time = (linspace(1, (frame(end)/expData.parameters.FPS)/60, frame(end)))';
occupancy.time = time;

% Data table with continuous variables:
T = table(frame, time, temperature, tempWork, vidNums, vidFrame);

initial_vars = [initial_vars; 'nflies'; 'demoImg','wellLabels'; 'occupancy'; 'T'];
clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Find the x-y coordinates for each fly

% how many tracks max within this whole experiment?
w = 0;
for vid = 1:nvids
    if isempty(data(vid).occupancy_matrix)
        continue
    end
    w = max([size(data(vid).occupancy_matrix,1),w]);
end
[X, Y] = deal(nan(T.frame(end),w));

% load fly locations into the X & Y matrices
for vid = 1:nvids
    if isempty(data(vid).occupancy_matrix)
        continue
    end
    % ---- fly tracked locations -------
    headData = squeeze(data(vid).tracks(:,1,:,:));
    % x-y coordinates of flies for each frame
    x_loc = squeeze(headData(:,1,:));
    y_loc = squeeze(headData(:,2,:));
    data(vid).x_loc_raw = x_loc;
    data(vid).y_loc_raw = y_loc;
    dims = size(x_loc); 
    roi_1 = occupancy.frameROI(vid,1):occupancy.frameROI(vid,2);
    roi_2 = 1:dims(2);
    X(roi_1,roi_2) = x_loc;
    Y(roi_1,roi_2) = y_loc;
end

% Add X and Y to the data table:
T = addvars(T,X,Y);
% head(T,5)

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% Determine the fly tracked points for each arena

flyCount = [];
% which points are within the sphere of each arena?
for arena = 1:4
    % Pull variables:
    X = T.X;
    Y = T.Y;
    centre = arenaData(arena).centre; % arena center for each of the four arenas
    c1 = centre(1);
    c2 = centre(2);
    
    % Find points within arena:
    loc = sqrt((X-c1).^2 + (Y-c2).^2)<=r; % tracked points within circle
    X(~loc) = nan;
    Y(~loc) = nan;
    flyNum = sum(loc,2); % how many points are found within that arena = number of flies

    % Save into data structure: 
    data(vid).x_loc = X;
    data(vid).y_loc = Y;

    % Remove some nans to make smaller matrix
    dims = size(X);
    [x_loc, y_loc] = deal(nan(dims(1),max(flyNum)));
    for ii = 1:dims(1)
        keep_loc = ~isnan(X(ii,:));
        x_loc(ii,1:flyNum(ii)) = X(ii,keep_loc);
        y_loc(ii,1:flyNum(ii)) = Y(ii,keep_loc);
    end
    arenaData(arena).x_loc = x_loc;
    arenaData(arena).y_loc = y_loc;
    
    % fly count for this arena
    flyCount = [flyCount,flyNum];
    arenaData(arena).flyCount = flyNum;
    % Save to structures
    switch arena
        case 1
            flyCount_A = flyNum;
            xA = x_loc;
            yA = y_loc;
        case 2
            flyCount_B = flyNum;
            xB = x_loc;
            yB = y_loc;
        case 3
            flyCount_C = flyNum;
            xC = x_loc;
            yC = y_loc;
        case 4
            flyCount_D = flyNum;
            xD = x_loc;
            yD = y_loc;
    end
    
end

% save the data to appropriate structures:
try 
    T.flyCount = flyCount;
    T.flyCount_A = flyCount_A;
    T.flyCount_B = flyCount_B;
    T.flyCount_C = flyCount_C;
    T.flyCount_D = flyCount_D;
    T.xA = xA; T.yA = yA;
    T.xB = xB; T.yB = yB;
    T.xC = xC; T.yC = yC;
    T.xD = xD; T.yD = yD;
catch 
    T = addvars(T,flyCount,flyCount_A,flyCount_B,flyCount_C,flyCount_D,xA,yA,xB,yB,xC,yC,xD,yD);
    T = movevars(T,'X','After','flyCount_D');
    T = movevars(T,'Y','After','X');
end

occupancy.flyCount = flyCount;

initial_vars = [initial_vars; 'flyCount'];
clearvars('-except',initial_vars{:})
fprintf('Next\n')

T = movevars(T,'X','After','flyCount_D');
T = movevars(T,'Y','After','X');
% head(T,1)
% T.Properties.VariableNames
% T = removevars(T,{'flyCount_1','flyCount_A_1','flyCount_C_1','flyCount_D_1','flyCount_B_1'});

%% FIGURES: check on the number of flies per arena

n = 1; % how many miscounted frames to look at

% Check flycount offset by arena:
for arena = 1:4
    save_loc = [baseFolder folder '/Arena ' arenaIdx{arena} '/' expName ' over and Under Tracked Images.pdf'];
    offset = flyCount(:,arena)-nflies(arena);
    [~,idx] = sort(offset);
    lowIDX = idx(1:n);          % lowest fly count frame index
    highIDX = idx(end-n+1:end); % highest fly count frame index

    % -------------- LOW COUNTS ---------------
    % load video information
    for ii = 1:n
        frame = lowIDX(ii);
        vidNum = T.vidNums(frame);
        if isempty(data(vidNum).occupancy_matrix)
            continue
        end
        vidframe = T.vidFrame(frame);
        movieInfo = VideoReader([baseFolder folder '/' expName '_' num2str(vidNum) '.avi']); %read in video
        img = read(movieInfo,vidframe);
        fig = figure;
        if ismac
            set(fig, 'pos', [534 10 972 967]);
        end
            imshow(img); set(fig,'color', 'k')
            hold on
            x = T.X(frame,:);
            y = T.Y(frame,:);
            scatter(x,y, 10, 'y')
            % draw arena circle
            kolor = arenaData(arena).color;
            centre = arenaData(arena).centre;
            viscircles(centre', r, 'color', kolor);
            % write the number of frames that are offset
            for a = 1:4
                AC = arenaData(a).centre;
                str = num2str(T.flyCount(frame,a)-nflies(a));
                text(AC(1),AC(2),str,'color','w','fontsize',14,'horizontalAlignment', 'center')
            end
            % Add title with temperature
            title([strrep(expName,'_',' ') ' ' folder ' | Temperature: ' num2str(T.temperature(frame))],...
                'color', 'w')
        % save image??
        export_fig(fig, save_loc, '-pdf', '-nocrop', '-r80' , '-painters', '-rgb','-append');
        close(fig)
    end

     % -------------- HIGH COUNTS ---------------
    % load video information
    for ii = 1:n
        frame = highIDX(ii);
        vidNum = T.vidNums(frame);
        if isempty(data(vidNum).occupancy_matrix)
            continue
        end
        vidframe = T.vidFrame(frame);
        movieInfo = VideoReader([baseFolder folder '/' expName '_' num2str(vidNum) '.avi']); %read in video
        img = read(movieInfo,vidframe);
        fig = figure;
        if ismac
            set(fig, 'pos', [534 10 972 967]);
        end
            imshow(img); set(fig,'color', 'k')
            hold on
            x = T.X(frame,:);
            y = T.Y(frame,:);
            scatter(x,y, 10, 'y')
            % draw arena circle
            kolor = arenaData(arena).color;
            centre = arenaData(arena).centre;
            viscircles(centre', r, 'color', kolor);
            % write the number of frames that are offset
            for a = 1:4
                AC = arenaData(a).centre;
                str = num2str(T.flyCount(frame,a)-nflies(a));
                text(AC(1),AC(2),str,'color','w','fontsize',14,'horizontalAlignment', 'center')
            end
            % Add title with temperature
            title([strrep(expName,'_',' ') ' ' folder ' | Temperature: ' num2str(T.temperature(frame))],...
                'color', 'w')
        export_fig(fig, save_loc, '-pdf', '-nocrop', '-r100' , '-painters', '-rgb','-append');
        close(fig)
    end
end

%% FIGURES : Visualize tracking across all four quadrants
nrows = 2;
ncols = 3;
sb(1).idx = 4; %A
sb(2).idx = 1; %B
sb(3).idx = 5; %C
sb(4).idx = 2; %D
sb(5).idx = [3,6]; % grouped fly count histogram

fig = figure; set(fig, 'pos', [99 200 895 635])
    for arena = 1:4
        subplot(nrows, ncols, sb(arena).idx)
        hold on
        kolor = arenaData(arena).color;
        histogram(arenaData(arena).flyCount,'FaceColor', kolor,'EdgeColor','w','facealpha', 1)
        v_line(arenaData(arena).nflies,'w','--',2)
        xlabel('# Flies')
        ylabel('Count')
        title(arenaData(arena).name)
    end
    subplot(nrows, ncols, sb(5).idx)
    histogram(sum(flyCount,2),'FaceColor', Color('grey'),'EdgeColor','w')
    v_line(sum(nflies),'r','-',2)
    xlabel('# Flies')
    ylabel('Count')
    title('Full plate')
    formatFig(fig, true, [nrows,ncols],sb);
    figure(fig)
save_figure(fig, [analysisDir expName ' Fly Count Histogram'],'-png',false,true,'-r80');

% Tracking with temperature
sSpan = 360;
LW = 2;
nrows = 3; ncols = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;
fig = figure; set(fig, 'pos', [298 260 508 483])
subplot(nrows, ncols, sb(1).idx)
plot(T.time, T.temperature,'color', 'w','linewidth', LW)
xlabel('Time (min)')
ylabel('Temp (\circC)')
subplot(nrows, ncols, sb(2).idx)
hold on
for arena = 1:4
    y = arenaData(arena).flyCount-arenaData(arena).nflies;
    plot(T.time, smooth(y,sSpan), 'color', arenaData(arena).color, 'linewidth', LW)
end
h_line(0,'grey',':',1)
xlabel('Time (min)')
ylabel('Fly count offset (#)')
formatFig(fig, true, [nrows, ncols],sb);
figure(fig)
save_figure(fig, [analysisDir expName ' Fly Count over time'],'-png',false,true,'-r80');

clearvars('-except',initial_vars{:})

%% Find number of flies within each well sphere & fly distance to wells
% THIS SECTION WOULD NEED TO BE RERUN FOR ALL TRIALS

% V2_UPDATE (5.21.25)
% find the plate for this trial: 
plate = excelfile{XLrow(1),Excel.plate};
[conversion, con_type] = getConversion(folder, plate, 1); 
for arena = 1:4 
    arenaData(arena).plate = plate;
    arenaData(arena).con_type = con_type;
    arenaData(arena).pix2mm = conversion(con_type).pix2mm;
end

pix2mm = conversion(con_type).pix2mm;
radii = conversion(con_type).circle10 * pix2mm;
% radii = 165; % PREVIOUS VERSION
for arena = 1:4
    % Pull fly coordinate position data for this arena 
    X = T.(['x' arenaIdx{arena}]);
    Y = T.(['y' arenaIdx{arena}]);
    
    % Get center position for each well
    centers = arenaData(arena).wellcenters;

    % Calculate distance to food and occupancy counts
    for well = 1:4
        % center of the well coordinates
        c1 = centers(1,well);
        c2 = centers(2,well);

        % Find distance to well center
        dist2well = sqrt((X-c1).^2 + (Y-c2).^2);
        
        % Find points within arena:
        loc = dist2well<=radii; % tracked points within circle
        N = sum(loc,2); % how many points tracked within that circle
        
        % store data
        arenaData(arena).dist2well(:,well) = mean(dist2well./pix2mm,2,'omitnan');
        arenaData(arena).dist2well_err(:,well) = std(dist2well./pix2mm,0,2,'omitnan');
        arenaData(arena).occ_N(:,well) = N;
        arenaData(arena).occ_P(:,well) = N/nflies(arena);
    end
end

initial_vars{end+1} = 'pix2mm';
clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% FIGURE: distance to food and well occupancy across the experiment
sSpan = 360;
LW = 1;
nrows = 6;
ncols = 1;
sb(1).idx = 1;   % tracking
sb(2).idx = 2;   % temperature
sb(3).idx = 3:4; % occupancy
sb(4).idx = 5:6; % distance to food

for arena = 1:4
    fig = figure; set(fig, 'pos', [297 48 1065 820])
    % TRACKING 
    subplot(nrows,ncols,sb(1).idx)
        plot(T.time, smooth(T.flyCount(:,arena),sSpan),'color', arenaData(arena).color,'linewidth', LW)
        ylabel('# flies') 
    % TEMPERATURE
    subplot(nrows,ncols,sb(2).idx)
        plot(T.time, T.temperature,'color', 'w','linewidth', LW)
        ylabel('Temp (\circC)')
    % OCCUPATION
    subplot(nrows,ncols,sb(3).idx)
        hold on
        for well = 1:4
            kolor = pullFoodColor(arenaData(arena).wellLabels{well});
            y = smooth(arenaData(arena).occ_P(:,well),'moving',sSpan);
            plot(T.time, y, 'linewidth', LW, 'color', kolor);
        end
        ylabel('occupation probability')

    % DISTANCE TO FOOD
    subplot(nrows,ncols,sb(4).idx)
        hold on
        for well = 1:4
            kolor = pullFoodColor(arenaData(arena).wellLabels{well});
            y = smooth(arenaData(arena).dist2well(:,well),'moving',sSpan);
            plot(T.time, y, 'linewidth', LW, 'color', kolor);
        end
        ylabel('distance to food (mm)')
        xlabel('time (min)');

    formatFig(fig, true, [nrows,ncols],sb);
    subplot(nrows,ncols,sb(3).idx)
    l = legend(strrep(arenaData(arena).wellLabels,'_','-'));
    set(l, 'box', 'off', 'textcolor', 'w','edgecolor', 'k','location', 'northwest');
    subplot(nrows, ncols, sb(1).idx)
    titleName = strrep({[folder ' ' expName ' Arena ' arenaIdx{arena}];...
                         expData.parameters.(['Arena' Alphabet(arena)]).genotype}, '_',' ');
    title(titleName,'color', 'w')

    save_loc = [baseFolder folder '/Arena ' arenaIdx{arena} '/' expName ' time course overview'];
    save_figure(fig, save_loc, '-png',true,true,'-r80');
end

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% ------------------- Save preformatted data for QuadStep2 ------------------------
disp('Saving data...')
facility = 'college';

%------ Write experiments start time and location to the excel summary file ------
searchPath = [baseFolder folder '/' expName '*_1.avi'];
videoList = dir(searchPath);
videoStartTime = videoList(1).date(end-7:end);

% write the experiment details into the excel sheet
isExcelFileOpen(xlFile); % check that file details can be written to spreadsheet
for arena = 1:4
    writecell({videoStartTime},xlFile,'Sheet','Exp List','Range', [Alphabet(Excel.starttime) num2str(XLrow(arena))])
    writecell({facility},xlFile,'Sheet','Exp List','Range', [Alphabet(Excel.facility) num2str(XLrow(arena))])
    writecell({'Y'}, xlFile, 'Sheet','Exp List','Range',[Alphabet(Excel.step1) num2str(XLrow(arena))]);
end

clearvars('-except',initial_vars{:})
save([analysisDir expName ' preformed data'],'-v7.3')
disp('Formatted data saved')
disp('Done')








%% 
% %% Find the fly counts and head tracking points for each frame
% % flyCount = nan(occupancy.frameROI(end,2),4);
% % 
% % for vid = 1:nvids    
% %     % all data for head tracked location
% %     headData = squeeze(data(vid).tracks(:,1,:,:));
% %     
% %     % x-y coordinates of flies for each frame
% %     x_loc = squeeze(headData(:,1,:));
% %     y_loc = squeeze(headData(:,2,:));
% %     data(vid).x_loc_raw = x_loc;
% %     data(vid).y_loc_raw = y_loc;
% %     nframes = size(headData,1); 
% %     
% %     % --------- Find data points that lie within each arena ---------
% %     Xdim = size(x_loc);
% %     Ydim = size(y_loc);
% %     for arena = 1:4
% %         %resize the data:
% %         X = reshape(x_loc,[numel(x_loc),1]);
% %         Y = reshape(y_loc,[numel(y_loc),1]);
% %         centre = arenaData(arena).centre;
% %         % find points within each arena sphere
% %         loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)>=r;
% %         X(loc) = NaN;
% %         Y(loc) = NaN;
% %         %resize the data:
% %         X = reshape(X,Xdim);
% %         Y = reshape(Y,Ydim);
% % 
% %         
% %     
% %         % save tracking points into structure:
% %         data(vid).(arenaIdx{arena}).x_loc = X;
% %         data(vid).(arenaIdx{arena}).y_loc = Y;
% %         FliesPerFrame = sum(~isnan(X),2);
% %         data(vid).(arenaIdx{arena}).FliesPerFrame = FliesPerFrame;
% %         % number of flies within this arena:
% %         ROI = occupancy.frameROI(vid,1):occupancy.frameROI(vid,2);
% %         flyCount(ROI, arena) = FliesPerFrame;
% %     end
% % end
% % for arena = 1:4
% %     arenaData(arena).flyCount = flyCount(:,arena);
% % end
% % occupancy.flyCount = flyCount;
% % 
% % initial_vars = [initial_vars; 'flyCount'];
% % clearvars('-except',initial_vars{:})
% % fprintf('Next\n')
% %% --------------------------- Mask the food wells & arena ---------------------------
% % maskPath = [analysisDir expName ' ArenaMask.mat'];
% % if ~isfile(maskPath)
% %     idx = find(~strcmpi('Empty', wellLabels));
% %     if ~isempty(idx)
% %         % draw out the masks:
% %         for ii = 1:length(idx)
% %             wellName = strrep(wellLabels{idx(ii)}, '_', '-');
% %             % label the ROI
% %             f = figure;
% %             imshow(demoImg)
% %             title(['Outline the ' wellName ' well'])
% %             roi = drawpolygon;
% %             mask(ii).n = roi.Position;
% %             uiwait(f)
% %         end
% %         % save the food masks:
% %         save(maskPath, 'mask')
% %     end
% % else
% %     load(maskPath)
% % end
% 
% % % Check the video size
% % if all(size(demoImg)==946)
% %     r = 400; %radius of the arena for the quadbowl
% % else
% %     warndlg('Not expected video pixel size...manually set radius')
% %     return
% %     
% %     % manually select an arena radius?
% %     
% %     % % pic preview of mask circle:
% %     % figure; hold on
% %     % imshow(demoImg)
% %     % viscircles(centre', r);
% % end
% 
% % % Screen any points outside the arena
% % x1 = wellcenters(1,1:2:4);
% % y1 = wellcenters(2,1:2:4);
% % x2 = wellcenters(1,2:2:4);
% % y2 = wellcenters(2,2:2:4);
% % [xi,yi] = polyxpoly(x1,y1,x2,y2);
% % centre = [xi;yi];
% % % 
% % % pic preview of mask circle:
% % figure; hold on
% % imshow(demoImg)
% % viscircles(centre', r);
% 
% % % Tracking matrix locations: [frame, node, xy, fly]
% % scaleFactor = 1.31; % offset from the video encoding compression
% 
% 
% 
% 
% 
% 
% 
% % MASK FOOD AND ARENA:
% [raw_flyCount, mid_flyCount] = deal([]);
% for vid = 1:nvids
%     % number of flies labeled on each frame:
%     raw_flyCount = [raw_flyCount; sum(data(vid).occupancy_matrix)'];
%     
%     % all data for head tracked location
%     headData = squeeze(data(vid).tracks(:,2,:,:));
% %     headData = squeeze(data(vid).tracks(:,1,:,:)); % TODO REVERT TO HEADPOINT
%     %save the unadjusted data in the same format
%     data(vid).x_loc_raw = squeeze(headData(:,1,:)).*scaleFactor;
%     data(vid).y_loc_raw = squeeze(headData(:,2,:)).*scaleFactor;
%     nframes = size(headData,1); 
%     
%     % x-y coordinates of flies for each frame
%     x_loc = squeeze(headData(:,1,:)).*scaleFactor;
%     y_loc = squeeze(headData(:,2,:)).*scaleFactor;
% 
%     % REMOVE FOOD TRACKED POINTS HERE 
%     Xdim = size(x_loc);
%     Ydim = size(y_loc);
%     %resize the data:
%     X = reshape(x_loc,[numel(x_loc),1]);
%     Y = reshape(y_loc,[numel(y_loc),1]);
%     % Find points within the masked region and turn to NaN
%     if exist('mask','var')
%         for ii = 1:length(mask)
%             [in,on] = inpolygon(X,Y, mask(ii).n(:,1),mask(ii).n(:,2));   % Logical Matrix
%             inon = in | on;                                    % Combine ‘in’ And ‘on’
%             X(inon) = NaN;
%             Y(inon) = NaN;
%         end
%     end
%     % Remove any labeled points outside the arena:
%     loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)>=r;
%     X(loc) = NaN; Y(loc) = NaN;
% 
%     % Resize the data to OG structure:
%     X = reshape(X,Xdim);
%     Y = reshape(Y,Ydim);
% 
%     data(vid).x_loc_mid = X; % save for later convenience
%     data(vid).y_loc_mid = Y;
%     
%     % New fly count measure:
%     mid_flyCount = [mid_flyCount; sum(~isnan(X),2)];
% end
% 
% %% --------------------------- Mask the overtracked points ---------------------------
% removalRadius = 10;
% % Find the video with the greatest avg over-count of fly points:
% for vid = 1:nvids
%     ROI = occupancy.frameROI(vid,:);
%     numberCount(vid) = median(mid_flyCount(ROI(1):ROI(2)));
% end
% % find the top 4 frames and display them
% ptsPath = [analysisDir expName arenaSel ' FalsePoints.mat'];
% if ~isfile(ptsPath)
%     pullFrames = 5;
%     frameList = [];
%     vid_idx = find(numberCount == max(numberCount)); vid_idx = vid_idx(1);
%     ROI = occupancy.frameROI(vid_idx,:);
%     numlist = mid_flyCount(ROI(1):ROI(2));
%     [~,Idx] = (sort(numlist));
%     frameList = Idx(end-pullFrames+1:end); %pull the four highest overcounts
% 
%     % pull up the images in order and click on points that ARE NOT FLIES:
%     movieInfo = VideoReader([baseFolder vidFolder '\' expName '_' num2str(vid_idx) '.avi']); %read in video
%     for ii = 1:length(frameList)
%         frame = frameList(ii);
%         img = read(movieInfo,frame);
%         fig = figure;
%         imshow(img);
%             hold on
%             x = data(vid_idx).x_loc_mid(frame,:); x(isnan(x)) = [];
%             y = data(vid_idx).y_loc_mid(frame,:); y(isnan(y)) = [];
%         scatter(x,y, 10, 'y')
%         title(['select all points that are NOT flies ' num2str(ii) '/' num2str(pullFrames)])
%         pointLabels(ii).coord = labelWrongPoints(fig);
%     end
% 
%     % Now determine how they are related and if there is consistent overlap
%     % that can be targeted for deletion...
%     fig = figure; set(fig, 'color', 'k')
%     imshow(img)
%     hold on
%     for ii = 1:length(frameList)
%         scatter(pointLabels(ii).coord(1,:), pointLabels(ii).coord(2,:),20, 'filled')
%     end
%     title('Select points with overlap for masking')
%     % save_figure(fig, [analysisDir expName arenaSel ' overtracking point IDs'], '-png');
%     pointLabels(1).finalRound = labelWrongPoints(fig);
%     nmaskpoints = length(frameList);
% 
%     % visualize the selected ROIs 
%     fig = figure; set(fig, 'color', 'k');
%     imshow(img)
%     hold on
%     for ii = 1:nmaskpoints
%         scatter(pointLabels(ii).coord(1,:), pointLabels(ii).coord(2,:),20, 'filled')
%     end
%     R = removalRadius*ones(size(pointLabels(1).finalRound,2),1);
%     viscircles(pointLabels(1).finalRound', R ,'Color','r');
%     save_figure(fig,[analysisDir,expName,arenaSel,' overtracking points selected for deletion'], '-png');
%     
%     % save the food masks:
%     save(ptsPath, 'frameList','vid_idx', 'pointLabels','nmaskpoints')
% else
%     load(ptsPath)
% end
% 
% % how many flies are in those circles?? TODO integrate the last OG fly count...
% flyCount = [];
% for vid = 1:nvids
% %resize the data:
%     x_loc = data(vid).x_loc_mid; 
%     y_loc = data(vid).y_loc_mid;
%     
%     % how many flies OG
%     Xdim = size(x_loc); 
%     Ydim = size(y_loc);
%     X = reshape(x_loc,[numel(x_loc),1]); 
%     Y = reshape(y_loc,[numel(y_loc),1]);
% 
%     % Find points within the masked region and turn to NaN
%     if ~isempty(pointLabels(1).finalRound)
%         for ii = 1:size(pointLabels(1).finalRound,2)
%             loc = (((X-pointLabels(1).finalRound(1,ii)).^2 + (Y-pointLabels(1).finalRound(2,ii)).^2).^0.5)<=removalRadius;
%             X(loc) = NaN; Y(loc) = NaN;
%         end
%     end
%     % Resize the data to OG structure:
%     X = reshape(X,Xdim);
%     Y = reshape(Y,Ydim);
%     flyCount = [flyCount; sum(~isnan(X),2)];
%     
%     data(vid).x_loc = X; % save for later convience
%     data(vid).y_loc = Y;
% end
% 
% % save long-term fly count data:
% occupancy.time = linspace(1,(length(mid_flyCount)/3)/60,length(mid_flyCount));
% occupancy.raw_flyCount = raw_flyCount;
% occupancy.mid_flyCount = mid_flyCount;
% occupancy.flycount = flyCount;
% 
% %%  ------------------------------ Visualization --------------------------------------
% %read in video
% movieInfo = VideoReader([baseFolder,vidFolder,'\',expName,'_',num2str(vid_idx),'.avi']); 
% frame = frameList(1);
% img = read(movieInfo,frame);
% vid = vid_idx;
% nrows = 4;
% ncols = 3;
% sb(1).idx = [1,2,4,5]; %imshow image
% sb(2).idx = [7,8,10,11]; %time course frame count
% sb(3).idx = 3:3:12; %histogram
%     
% fig = figure; set(fig, 'pos', [300 100 1171 826]); % X-off, Y-off, width, height
% % ARENA IMAGE
% subplot(nrows,ncols,sb(1).idx)
%     imshow(img)
%     axis tight square
%     hold on
%     for well = 1:4
%         kolor = pullFoodColor(wellLabels{well});
%         scatter(wellcenters(1,well),wellcenters(2,well), 75,...
%             'MarkerFaceColor', kolor, 'MarkerEdgeColor', 'w') 
%     end 
%     % raw points:
%     x = data(vid).x_loc_raw(frame,:);
%     y = data(vid).y_loc_raw(frame,:);
%     x(isnan(x)) = []; % remove empty tracks
%     y(isnan(y)) = [];
%     scatter(x,y, 20, 'k')
%     scatter(x,y, 15, Color('grey'),'filled')
%     
%     % intermediate points:
%     x = data(vid).x_loc_mid(frame,:);
%     y = data(vid).y_loc_mid(frame,:);
%     x(isnan(x)) = []; % remove empty tracks
%     y(isnan(y)) = [];
%     scatter(x,y, 15, Color('orangered'),'filled')
%     
%     % intermediate points:
%     x = data(vid).x_loc(frame,:);
%     y = data(vid).y_loc(frame,:);
%     x(isnan(x)) = []; % remove empty tracks
%     y(isnan(y)) = [];
%     scatter(x,y, 15, Color('teal'),'filled')
% 
% % OVERTRACKING OVER TIME
% subplot(nrows,ncols,sb(2).idx)
%     hold on
%     plot(occupancy.time, raw_flyCount, 'color', Color('grey'))
%     plot(occupancy.time, mid_flyCount,'color', Color('orangered'))
%     plot(occupancy.time, flyCount, 'color', Color('teal'))
%     hline(nflies, 'w')
%     ylabel('Fly count'); xlabel('time (min)')
%     axis tight
%     yyaxis right
%     plot(occupancy.time,full_temp,'color', 'y')
%     ylabel('temp (\circC)')
%     xlim([0,occupancy.time(end)+10])
% 
% % FLY COUNT HISTOGRAM
% subplot(nrows,ncols,sb(3).idx)
%     yyaxis right
%     hold on
%     h = histogram(raw_flyCount);
%     h.FaceColor = Color('grey');
%     h = histogram(mid_flyCount);
%     h.FaceColor = Color('orangered'); 
%     h = histogram(flyCount);
%     h.FaceColor = Color('teal');
%     vline(nflies, 'w')
%     xlabel('Number of flies')
%     ylabel('Frame count')
% 
% % LABELS AND FORMATTING
% fig = formatFig(fig, true, [nrows, ncols], sb);
% l = legend({['SLEAP ' num2str(mean(raw_flyCount)) ' avg'],...
%             ['wells & arena masks ' num2str(mean(mid_flyCount)) ' avg'],...
%             ['wells, arena & manual ' num2str(mean(flyCount)) ' avg']});
% set(l,'color','k','textcolor','w','edgecolor','k','Position', [0.0169 0.8981 0.3122 0.0877]);
% subplot(nrows,ncols,sb(2).idx)
% yyaxis left
% set(gca, 'ycolor', 'w')
% subplot(nrows,ncols,sb(3).idx)
% yyaxis left
% set(gca, 'ycolor', 'k')
% 
% % save and export figure:
% if strcmpi(questdlg('Append figure to summary pdf?'),'Yes')
%     export_fig(fig, expPDF, '-pdf', '-nocrop', '-r300' , '-painters', '-rgb','-append');
% end  
% save_figure(fig, [analysisDir expName arenaSel ' quality control'], '-png');
% 
% 
% % ----------------------------- zoom-in-take on the tracking correction -------------
% % % visual confirmation that the selected points are near the well:
% % AllPoints = [];
% % for vid = 1:nvids
% %    X = reshape(data(vid).x_loc,numel(data(vid).x_loc),1); 
% %    Y = reshape(data(vid).y_loc,numel(data(vid).y_loc),1); 
% %    X(isnan(X)) = [];
% %    Y(isnan(Y)) = []; 
% %    AllPoints = [AllPoints; X , Y];
% % end
% % fig = getfig; set(fig, 'color', 'k');
% % hist2d(AllPoints(:,1),AllPoints(:,2), 'probability', 'tile')
% % axis tight; axis square
% % set(gca, 'visible', 'off')
% % c = colorbar;
% % c.Color = [1,1,1];
% % hold on 
% % % c = drawcircle('Center', wellcenters(:,5)', 'Radius', r);
% % % hold on
% % % viscircles(wellcenters(:,1:4)',[radii,radii,radii,radii])
% % % b = wellcenters(:,well); % well 1 points
% % save_figure(fig, [analysisDir expName arenaSel ' cropped out food locations quality control'], '-png');
% 
% %% Temp vs fly count
% 
% fig = figure;
% idx = discretize(full_temp,1:30);
% boxplot(flyCount, idx,'Plotstyle', 'traditional','colors', 'w')
% h_line(nflies, 'Teal', '-',3)
% xlabel('Temperature (\circC)')
% ylabel('Tracking Fly Count')
% fig = formatFig(fig, true);
% 
% save_figure(fig, [analysisDir expName arenaSel ' fly count v temperature'], '-png',true);
% 
% 
% %% ------------------- Save preformatted data for QuadStep2 ------------------------
% disp('Saving data file...')
% clearvars('-except',initial_vars{:})
% save([analysisDir expName ' preformed data'])
% disp('Formatted data saved')
% disp('Done')

%% Plot work vs temp for full experiment 
% 
% x = tempLog(:,3); % temperature
% y = tempLog(:,4); % work
% 
% 
% fig = figure;
%     scatter(x,y,45, 'w')
%     xlabel('temperature (\circC)')
%     ylabel('work (%)')
%     formatFig(fig,true);
% save_figure(fig, [analysisDir expName ' temperature vs work'], '-png',true);
% 
% 
% fig = figure;
%     yyaxis left
%     plot(x,'color','b','linewidth',1.5)
%     ylabel('temperature (\circC)')
%     yyaxis right
%     plot(y,'color','w','linewidth',1.5)
%     ylabel('work (%)')
%     xlabel('time (au)')
%     formatFig(fig,true);
% save_figure(fig, [analysisDir expName ' temperature overlay with work'], '-png',true);
% 
% 
% 
% fig = figure;
%     yyaxis left
%     plot(smooth(x,180,'moving'),'color','b','linewidth',1.5)
%     ylabel('temperature (\circC)')
%     yyaxis right
%     plot(smooth(y,180,'moving'),'color','w','linewidth',1.5)
%     ylabel('work (%)')
%     xlabel('time (au)')
%     formatFig(fig,true);
% save_figure(fig, [analysisDir expName ' smoothed temperature overlay with work'], '-png',true);
% 
% 





















