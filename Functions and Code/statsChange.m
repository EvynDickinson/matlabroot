% % intact walking fly analysis
function [diff, err] = statsChange(inputdata,cROI,sROI)

    a = inputdata;
    pre = nanmean(a(cROI,:));
    post = nanmean(a(sROI,:));
    strt = median(pre);
    stp = median(post);
    diff = stp-strt;
    
%     error calc
    b = pre-post;
    err = std(b);
%     
%     a = inputdata;
%     pre = nanmean(a(cROI,:),2);
%     post = nanmean(a(sROI,:),2);
%     strt = mean(pre);
%     stp = mean(post);
%     diff = stp-strt;
%     
%     % error calc
%     b = pre-post;
%     err = std(b);
end

% % interneuron joint angle analysis
% function diff = statsChange(inputdata,cROI,sROI)
% 
%     a = inputdata;
%     pre = mean(a(cROI,:));
%     post = mean(a(sROI,:));
%     strt = mean(pre);
%     stp = mean(post);
%     diff = stp-strt;
% 
% end


