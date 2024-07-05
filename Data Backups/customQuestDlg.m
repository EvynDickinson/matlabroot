

function actionChoice = customQuestDlg()
    options = {'Raw Datal', 'Single Trial', 'Pooled Trial', 'Grouped Data', 'Cancel'};
    fig = uifigure('Position', [600 600 300 250], 'Name', 'Select an option');

    lbl = uilabel(fig, 'Position', [20 150 260 30], 'Text', 'Select the data type to access:');

    for i = 1:numel(options)
        uibutton(fig, 'Position', [20 150-30*i 260 30], 'Text', options{i}, ...
                 'ButtonPushedFcn', @(btn, event) buttonCallback(fig, options{i}));
    end
end

function choice = buttonCallback(fig, option)
    disp(['You selected: ' option]);
    choice = option;
    delete(fig);  % Close the dialog
end
