 


function [newdata, conditions] = insert_averages_in_holes(testdata, conditions, parameters, missing_values)
% newdata = insert_averages_in_holes(testdata)
% Insert an average line of data into the fictrac data when there was a
% dropped frame during the control or stimulus periods
% Inputs:
% 'testdata' [raw data from fictrac (23 or 24 columns)]
% 'conditions' [structure with the start|end fictrac frame number info]
% 'missing_values' [structure with condition location for missing numbers]
% Outputs:
% 'newdata' [fictrac data with averages stuck into the dropped frames loc.]
% prints any further missing locations
%  Load variables
Type = {'Stim', 'Control'};

switch nargin
    case 3 %manual
        param = str2double(cell2mat(inputdlg('Stim (1) or Control (2)?')));
        rep =   str2double(cell2mat(inputdlg('Rep number?')));    % missing rep number
        cond = str2double(cell2mat(inputdlg('Condition number?')));    % missing cond number    
    case 4 %automated
        param = missing_values(1);
        rep = missing_values(2);
        cond = missing_values(3);    
end
% Check if the missing data is that there are more or less data points than
% desired:

% Position of start and end points within the testdata ...
%(not garanteed to be the exact location due to dropped frames at start of recording)
frame_numbers = testdata(:,1); 
search_buffer = 10;
search_begin = find((conditions.fictracdata.(Type{param}).start_num(cond,rep)-search_buffer) == frame_numbers);
search_end = find((conditions.fictracdata.(Type{param}).end_num(cond,rep)+search_buffer) == frame_numbers);
%             search_begin = find((Start.(Type{param})(cond,jj)-search_buffer) == frame_numbers);
%             search_end = find((End.(Type{param})(cond,jj)+search_buffer) == frame_numbers);

% pull up the frame number info:
frame_starts = frame_numbers(search_begin:search_end, 1);
% find the differences between frame numbers, look for >1
missing_location = find(diff(frame_starts) > 1);
%if there are no frame number gaps but a shortage of a frame
if sum(missing_location) == 0
    a = length(frame_starts) - (2*search_buffer);
    switch param
        case 1 %stim missing data
            if a == (parameters.OL_time*parameters.fps)
              conditions.fictracdata.(Type{param}).end_num(cond,rep) = ...
                  conditions.fictracdata.(Type{param}).end_num(cond,rep)+1;
              fprintf(['\n Stim | Rep ' num2str(rep) ' | Cond ' num2str(cond) ...
                       ' adjusted from ' num2str(a-1) ' to ' num2str(a) ' frames\n'])
            end 
        case 2 %control missing data
             if a == (parameters.OL_time*parameters.fps)
              conditions.fictracdata.(Type{param}).end_num(cond,rep) = ...
                  conditions.fictracdata.(Type{param}).end_num(cond,rep)+1;
              fprintf(['\n Control | Rep ' num2str(rep) ' | Cond ' num2str(cond) ...
                       ' adjusted from ' num2str(a-1) ' to ' num2str(a) ' frames\n'])
            end 
    end
    newdata = testdata;
%missing point is in the middle of the stimulus
else
    points_previous = missing_location(1);
    points_next = missing_location(1)+1;
    frame_one = frame_starts(points_previous, 1);
    frame_two = frame_starts(points_next, 1);
    gap_location = find(frame_numbers == frame_one);
    num_missing_points = frame_two-frame_one-1;
    average_start = find(frame_numbers == frame_one);
    average_end = find(frame_numbers == frame_two);
    testdata_average = mean(testdata(average_start:average_end,:));

    % multiply the averages for the correct numbers of missing rows
    for ii = 1:num_missing_points
        td_average(ii,:) = testdata_average;
        td_average(ii,1) = frame_one + ii;
    end
    %insert the averaged gap fillers into a new variable along with original data
    newdata = [testdata(1:gap_location,:); td_average; testdata(gap_location+1:end,:)];
    %check for an empty missing location:
    frame_starts = newdata(search_begin:search_end,1);
    missing_location = find(diff(frame_starts) > 1);
    fprintf('\n Missing Location: \n -----------------\n')
    disp(missing_location)    

end 

% figure; hold all
% plot(frame_starts, 'b--o')
% hline(frame_numbers(search_begin))
% hline(frame_numbers(search_end))

end
