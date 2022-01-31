
function [threshHigh, threshLow] = getTempThresholds
% [threshHigh, threshLow] = getTempThresholds;
%
% Returns temperature thresholds from a list
% as selected by the user
%
%
% ES Dickinson, Jan 2022


tempList = {'Low (8-20)', 'Mid (6-25)','High (5-30)','Other'};
UserChoice = listdlg('ListString', tempList,'ListSize', [100, 150],'SelectionMode','single');

switch UserChoice
    case 1 %'Low (8-20)'
        threshHigh = 19.88;
        threshLow = 8.02;
    case 2 %'Mid (6-25)'
        threshHigh = 24.5;
        threshLow = 6.5;
    case 3 %'High (5-30)'
        threshHigh = 29.88;
        threshLow = 5.02;
    case 4 %'Other'
        prompt = {'Low temp threshold:','High temp threshold:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'5','13'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        threshLow = str2double(answer{1});
        threshHigh = str2double(answer{2});
        
    case 5 %'' or cancel
        return
end
