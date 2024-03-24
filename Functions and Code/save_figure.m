
function results = save_figure(fig_handle, figure_name, type, autoSave, closeFig,fig_quality)
% 
% results = save_figure(fig_handle, figure_name, type, autoSave, closeFig)
% 
% Export the input figure to the given location and name
% with the following settings: 
% '-png', '-nocrop', '-r300' , '-painters', '-rgb'
% Default type is '-png' but can be specified
%
% Inputs:
% 'fig_handle' [handle for figure being saved]
% 'figure_name' [path and name for saving location of fig]
% 'type' ['-pdf' or '-png' output type]
% 'autoSave' [true|false : save fig without user input]
%
% Outputs: 
% 'results' [logical true|false if figure is saved]
%     
% % ES Dickinson, University of Washington, Jan 2019    

% Defaults:
if nargin == 2
    type = '-pdf'; % pdf filetype
end
if nargin < 4
    autoSave = false; % manually choose to save/not save
end
if nargin<5
    closeFig = true;
end
if nargin<6
    fig_quality = '-r300';
end

% Save figure:
if autoSave==true
    export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', fig_quality , '-painters', '-rgb');
    if closeFig
        close(fig_handle)
    end
    fprintf('\nSaved:')
    disp(figure_name) 
    results = true;
elseif autoSave == false
    switch questdlg('Save Image?', 'Figure', 'Save Figure', 'Close Figure', 'Cancel', 'Save Figure')
        case 'Save Figure'
            export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', fig_quality , '-painters', '-rgb');
            if closeFig
                close(fig_handle)
            end
            fprintf('\nSaved:')
            disp(figure_name) 
            results = true;
        case 'Close Figure'
            if closeFig
                close(fig_handle)
            end
            results = false;
        case 'Cancel'
            results = 'Cancel';
    end
end


% export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', '-r300' , '-painters', '-rgb');
% close(fig_handle)
% fprintf('\nSaved:')
% disp(figure_name) 
% results = true;


%  switch questdlg('Save Image?', 'Figure', 'Save Figure', 'Close Figure', 'Cancel', 'Save Figure')
%     case 'Save Figure'
%         savefig(fig_handle, figure_name)
%     case 'Cancel'
%         results = false;
% end
   
% switch questdlg('Save Image?', 'Figure', 'Save Figure', 'Close Figure', 'Cancel', 'Save Figure')
%     case 'Save Figure'
%         export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', '-r300' , '-painters', '-rgb');
%         close(fig_handle)
%         fprintf('\nSaved:')
%         disp(figure_name) 
%         results = true;
%     case 'Close Figure'
%         close(fig_handle)
%         results = false;
%     case 'Cancel'
%         results = false;
% end

end