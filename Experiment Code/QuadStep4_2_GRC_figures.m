
fig_type = '-pdf';
blkbgd = false;


%% Occupancy 'null' distribution for no food trials
clearvars('-except',initial_vars{:})

% link each data set to it's 'null' data set: (make this more involved and
% automated later -- gui driven, maybe?)

null_pair = [1,2]; % first idx is the test trial and the second is the null for each row

% Plot out the quadrant data:

np = 2; %null-pair idx
exp = 1; % active trial
% find the quads with the highest and lowest occupancy over the course of
% the experiment:

sSpan = 180;

dummy = [];
plotData = [];
for i = 1:4
    a = grouped(np).fullquad.(quadOrder{i}).all;
    dummy(i,:) = sum(a,1,'omitnan');
    plotData(:,:,i) = a;
end
% find min and max occupancy quadrants: 
[~, lowerIDX] = min(dummy);
[~, upperIDX] = max(dummy);

[minOcc,maxOcc] = deal([]);
for i = 1:num.trial(np)
    minOcc(:,i) = squeeze(plotData(:,i,lowerIDX(i)));
    maxOcc(:,i) = squeeze(plotData(:,i,upperIDX(i)));
end
    y_err = smooth(mean((minOcc-maxOcc)./2,2,'omitnan'),sSpan,'moving');
    y_avg = smooth(mean([minOcc,maxOcc],2,'omitnan'),sSpan, 'moving');



fig = getfig('',1); 
    hold on
    % plot the null distribution data
    y_err = smooth(mean((minOcc-maxOcc)./2,2,'omitnan'),sSpan,'moving');
    y_avg = smooth(mean([minOcc,maxOcc],2,'omitnan'),sSpan, 'moving');
    kolor = grouped(np).color;
    time = grouped(np).time;
    y1 = smooth(mean(minOcc,2),sSpan, 'moving');
    y2 = smooth(mean(maxOcc,2),sSpan, 'moving');
    plot_error_fills(true, time, y_avg,y_err,kolor,fig_type,0.5);
    % plot(time,y1,'color',kolor)
    % plot(time,y2,'color',kolor)
    % plot the paired food trial on top:
    x = grouped(exp).time;
    kolor = grouped(exp).color;
    y = grouped(exp).fullquad.food.avg;
    y_err = grouped(exp).fullquad.food.std./sqrt(num.trial(exp));
    plot_error_fills(true, x, y, y_err, kolor, fig_type);
    plot(x,y,'color',kolor,'linewidth', 1)
    % formatting
    xlabel('time (min)')
    ylabel('food quadrant occupancy (%)')
    xlim([0 700])
    formatFig(fig,false);
    
save_figure(fig,[figDir, 'full quad occ over time'],fig_type);


%% 




%% FIGURE: STATIC caviar trials -- plot heatmap of fly position within arena
clearvars('-except',initial_vars{:})
save_path = createFolder([saveDir 'COM/']);
autoSave = true;

[foreColor,backColor] = formattingColors(blkbgd);


% Find the occupancy for each bin:

n = 26; % number of spatial bins
autoLim = false;
axis_limits = [0, 0.01];

% Parameters: 
expList = 1:num.exp;
temp_list = [15 17 20 25 27 33 35]; % temps that we have temp hold data for...
temp_match_proto = 'Large_temp_sweep_15_35'; % this is where we want to match the time points for the experiment

% Find the time points that correspond to the desired temp period
tp = getTempTurnPoints(temp_match_proto);
temp_path = getCloudPath;
temp_path = temp_path(1:end-5);
temp = load([temp_path 'LTS 15-35 temp data.mat']); % pull in the fictive temp region information (and berlin caviar LTS data)
LTS = temp.LTS_temp; clear temp

% heating and cooling separated for each desired temp bin
% find the time index locations for the periods where the behavior should
% be compiled to match the fictive temp protocol

% grouped(exp).position.loc(rate, tempbin).data(trial).pos  

% find the time point index: 

for exp = 1:length(temp_list) % could also do this as auto find of the avg temp for the trial...
    temp = temp_list(exp); % this is the temp for the temp hold
    [~,idx] = min(abs(LTS.temp_list-temp)); %  bin index for this temp
    nRates = length(LTS.temp_rates);
    % find the frame numbers for the selected temp
    for rr = 1:nRates % for heating and cooling, respectively
        frames = LTS.loc(rr,idx).frames;
        % find the center of the arena for this exp type
        Cx = mean(grouped(exp).position.well_pos.x(5,:)); %center X
        Cy = mean(grouped(exp).position.well_pos.y(5,:)); %center Y
        for trial = 1:num.trial(exp)
            con_type = data(exp).con_type(trial);
            if any(con_type==[1 2]) % check if it matches plate 1 
            x = grouped(exp).position.trial(trial).x;
            y = grouped(exp).position.trial(trial).y;
            r = conversion(con_type).R*conversion(con_type).pix2mm; % pixel radius of this arena
            % find the edges of the spatial bins 
            x_edge = linspace(Cx-r,Cx+r,n);
            y_edge = linspace(Cy-r,Cy+r,n);
            % find which fly locations go to which spatial bin
            nanLoc = isnan(x)| isnan(y);
            x(nanLoc) = [];
            y(nanLoc) = []; % this also reorganizes as a linear bin
            xInd = discretize(x,x_edge);
            yInd = discretize(y,y_edge);

            % WORKING HERE 6.16 ESD
            % find the number of flies within each spatial bin:
            for row = 1:n
                for col = 1:n
                    nflies(row,col) = sum(yInd==row & xInd==col);
                end
            end
            % turn to prob and not direct occupancy
            plotData(i,rr).data = nflies./sum(sum(nflies));
            
            max_occ = max([max_occ,max(max(plotData(i,rr).data))]);
    
            % Find the wells within the binned space
            % wellX = (grouped(i).position.well_pos.x(1:4,:)); 
            % wellY = (grouped(i).position.well_pos.y(1:4,:));
            % wellX = wellX(:);
            % wellY = wellY(:);
            xInd = discretize(0,x_edge);
            yInd = discretize(0,y_edge);
    
            plotData(i,rr).wells = [xInd,yInd];




% Set Temperature
for temp = temp_list 

plotData = [];
max_occ = [];
% GROUP DATA
for i = expList
    for trial = 1:num.trial
        % r = 
        % get the 'square' units for partitioning space
        Cx = mean(grouped(i).position.well_pos.x(5,:)); %center X
        Cy = mean(grouped(i).position.well_pos.y(5,:)); %center Y
        x_edge = linspace(Cx-r,Cx+r,n);
        y_edge = linspace(Cy-r,Cy+r,n);
    
        % determine what a circle around the arena would look like:
        % r needs to be transformed into unit square space...
        square_unit = mean(diff(x_edge)); % pixel size for one bin
        circ_r = r/square_unit; % arena radius in bin size
        circ_X = discretize(Cx, x_edge);
        circ_Y = discretize(Cy, y_edge);
    
        % find the temp bin that corresponds to the selected temp:
        [~,idx] = min(abs(grouped(i).position.temp_list-temp));
        nRates = length(grouped(i).position.temp_rates);
        nflies = [];
        for rr = 1:nRates
            % find x and y that are within each 'box'
            x = grouped(i).position.loc(rr,idx).x;
            y = grouped(i).position.loc(rr,idx).y;
            nanLoc = isnan(x)| isnan(y);
            x(nanLoc) = [];
            y(nanLoc) = [];
    
            xInd = discretize(x,x_edge);
            yInd = discretize(y,y_edge);
        
            % find the number of flies within each spatial bin:
            for row = 1:n
                for col = 1:n
                    nflies(row,col) = sum(yInd==row & xInd==col);
                end
            end
            % turn to prob and not direct occupancy
            plotData(i,rr).data = nflies./sum(sum(nflies));
            
            max_occ = max([max_occ,max(max(plotData(i,rr).data))]);
    
            % Find the wells within the binned space
            % wellX = (grouped(i).position.well_pos.x(1:4,:)); 
            % wellY = (grouped(i).position.well_pos.y(1:4,:));
            % wellX = wellX(:);
            % wellY = wellY(:);
            xInd = discretize(0,x_edge);
            yInd = discretize(0,y_edge);
    
            plotData(i,rr).wells = [xInd,yInd];
        end
end

disp(['Max occupancy: ' num2str(max_occ)])

% PLOT 
fig_W = 20 + (400*nRates);

for i = expList
    fig = getfig('',false,[fig_W, 340]); 
    for rr = 1:nRates
        subplot(1,nRates,rr)
        hold on
        imagesc(plotData(i,rr).data); hold on
        scatter(plotData(i,rr).wells(:,1),plotData(i,rr).wells(:,2),10,'r','filled')
        axis tight;
        axis square;
        % h = drawcircle('Center',[circ_X,circ_Y],'Radius',circ_r,'StripeColor',foreColor);
        v = viscircles([circ_X,circ_Y],circ_r, 'color', foreColor);
    end
    formatFig(fig, blkbgd,[1,nRates]);
    for rr = 1:nRates
        subplot(1,nRates,rr)
        set(gca,'XColor',backColor,'Ycolor',backColor,'XTick', [],'YTick', [])
        t_str = [num2str(grouped(i).position.temp_rates(rr)) '\circC/min | ' num2str(temp) '\circC'];
        title({grouped(i).name; t_str},'color',foreColor,'fontsize', 12)
        
        % set(gca,'ColorScale','log')

        c = colorbar;
        c.Label.String = 'Occupancy Probability';
        c.Label.Color = foreColor;
        c.Color = foreColor;
        if autoLim
            clim([0,max_occ]) 
        else
            clim(axis_limits)
        end
    end
    colormap(flipud(gray))
    % save the figure to a folder specific to that cohort?
    save_figure(fig,[save_path grouped(i).name ' ' num2str(temp) ' deg'], fig_type,autoSave,true);
end

end
