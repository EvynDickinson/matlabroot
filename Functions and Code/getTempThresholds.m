
function [threshHigh, threshLow] = getTempThresholds
% [threshHigh, threshLow] = getTempThresholds;
%
% Returns temperature thresholds from a list
% as selected by the user
%
%
% ES Dickinson, Jan 2022


tempList = {'(8-20)', '(6-25)', '(6-26)', '(7-23)', '(5-30)','Other'};
UserChoice = tempList{listdlg('ListString', tempList,'ListSize', [100, 150],'SelectionMode','single')};

switch UserChoice
    case '(8-20)'
        threshHigh = 19.88;
        threshLow = 8.02;
    case '(6-25)'
        threshHigh = 24.5;
        threshLow = 6.5;
    case '(5-30)'
        threshHigh = 29.88;
        threshLow = 5.02;
    case '(6-26)'
        threshHigh = 26;
        threshLow = 6;
    case '(7-23)'
        threshHigh = 23;
        threshLow = 7;
    case 'Other'
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
