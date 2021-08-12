
function fig = DLC_plot_ROIForFit(STEP)
% fig = DLC_plot_ROIForFit(input, stance, ROI)
% plots the selected frames in black
% 
% ES Dickinson,
% University of Washington, 2020

LW = 1;
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
for leg = 1:length(leg_colors)
    C(leg,:) = Color(leg_colors{leg});
end

fig = getfig('',1);
for n = 1:length(STEP)
    if STEP(n).qual == false; continue; end
    
    min_loc = STEP(n).stance(:,2); %stance end
    max_loc = STEP(n).stance(:,1); %stance start
    ROI =  STEP(n).ROI;
    % Pull up the data:
    input = STEP(n).input;
    plot3(input(:,1),input(:,2),input(:,3), 'color', C(n,:), 'linewidth', LW)
    hold on
    scatter3(input(min_loc,1), input(min_loc,2), input(min_loc,3),100, 'g', 'filled') % start of stance
    scatter3(input(max_loc,1), input(max_loc,2), input(max_loc,3),100, 'm', 'filled') % start of swing
    
    % plot the selected region that will be used to fit the sphere:
    for ii = 1:size(ROI,1)
        roi = ROI(ii,1):ROI(ii,2);
        scatter3(input(roi,1), input(roi,2), input(roi,3), 10, 'k', 'filled')
    end
    grid on

end
legend({'tarsus position', 'stance end', 'stance start', 'ROI'})
title('Selected ROI for sphere fitting')
xlabel('X')
ylabel('Y')
zlabel('Z')


%     % OLD
%     min_loc = stance(:,2); %stance end
%     max_loc = stance(:,1); %stance start
% 
% fig = getfig('',1);
%     plot3(input(:,1),input(:,2),input(:,3))
%     hold on
%     scatter3(input(min_loc,1), input(min_loc,2), input(min_loc,3),100, 'g') % start of stance
%     scatter3(input(max_loc,1), input(max_loc,2), input(max_loc,3),100, 'm') % start of swing
%     
%     % plot the selected region that will be used to fit the sphere:
%     for n = 1:size(ROI,1)
%         roi = ROI(n,1):ROI(n,2);
%         scatter3(input(roi,1), input(roi,2), input(roi,3), 10, 'k', 'filled')
%     end
%     grid on
%     legend({'tarsus position', 'stance end', 'stance start', 'ROI'})
end