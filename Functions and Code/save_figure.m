
function results = save_figure(fig_handle, figure_name, type)
% 
% results = save_figure(fig_handle, figure_name)
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
%
% Outputs: 
% 'results' [logical true|false if figure is saved]
%     
% % ES Dickinson, University of Washington, Jan 2019    

%default 
if nargin == 2
    type = '-pdf';
end

% export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', '-r300' , '-painters', '-rgb');
% close(fig_handle)
% fprintf('\nSaved:')
% disp(figure_name) 
% results = true;



switch questdlg('Save Image?', 'Figure', 'Save Figure', 'Close Figure', 'Cancel', 'Save Figure')
    case 'Save Figure'
        export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', '-r300' , '-painters', '-rgb');
        close(fig_handle)
        fprintf('\nSaved:')
        disp(figure_name) 
        results = true;
    case 'Close Figure'
        close(fig_handle)
        results = false;
    case 'Cancel'
        results = 'Cancel';
end


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