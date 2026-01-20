
function results = save_figure(fig_handle, figure_name, type, autoSave, closeFig, fig_quality)
% 
% results = save_figure(fig_handle, figure_name, type, autoSave, closeFig, figure_quality)
%
% PURPOSE
% Export the input figure to the given location and name
% with the following settings: 
% '-png', '-nocrop', '-r300' , '-painters', '-rgb'
%
% INPUTS:
%   'fig_handle' : handle for figure being saved
%   'figure_name' : path and name for saving location of fig
%           e.g. 'S:\Evyn\DATA\this_figure_name'
%   'type' : saved format for the figure
%           common options: '-pdf' or '-png' output type (**current version saves both**)
%   'autoSave' : true|false : save fig without user input
%   'closeFig' : true|false : close the figure after saving
%   'figure_quality' : pixels per inch figure saving quality (format:  '-rNUM')
%           '-r80' is a quick save but '-r300' is good quality 
%           (default: '-r300' on PC  &  '-r100' on MAC)
%
% OUTPUTS: 
%   'results' : logical true|false if figure is saved
%
% UPDATE 12.4.25: AUTO SAVE BOTH A PDF AND A PNG FILE 
%     
% ES DICKINSON, 2019    

%%

% Defaults:
if nargin < 4
    autoSave = false; % manually choose to save/not save
end
if nargin<5
    closeFig = true;
end
if nargin<6 % DETERMINE FIGURE QUALITY
    if ismac
        fig_quality = '-r100'; %this massively increases the image saving time
    else
        fig_quality = '-r300'; 
end

warning off

% Save figure:
if autoSave==true
    % save PDF first: 
    export_fig(fig_handle, [figure_name '.pdf'], '-pdf', '-nocrop', fig_quality , '-painters', '-rgb');
    export_fig(fig_handle, [figure_name '.png'], '-png', '-nocrop', fig_quality , '-painters', '-rgb');
    if closeFig
        close(fig_handle)
    end
    fprintf('\nSaved:')
    disp(figure_name) 
    results = true;
elseif autoSave == false
    switch questdlg('Save Image?', 'Figure', 'Save Figure', 'Close Figure', 'Cancel', 'Save Figure')
        case 'Save Figure'
            export_fig(fig_handle, [figure_name '.pdf'], '-pdf', '-nocrop', fig_quality , '-painters', '-rgb');
            export_fig(fig_handle, [figure_name '.png'], '-png', '-nocrop', fig_quality , '-painters', '-rgb');
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



% LEGACY (PRE 12/4/25) CODE: 
% ------------------------------------------------------------------------------------------------------------
% % Defaults:
% if nargin == 2
%     type = '-pdf'; % pdf filetype
% end
% if nargin < 4
%     autoSave = false; % manually choose to save/not save
% end
% if nargin<5
%     closeFig = true;
% end
% if nargin<6
%     fig_quality = '-r300';
% end
% 
% if ismac && nargin<6
%     fig_quality = '-r100'; %this massively increases the image saving time
% end
% 
% 
% % Save figure:
% if autoSave==true
%     export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', fig_quality , '-painters', '-rgb');
%     if closeFig
%         close(fig_handle)
%     end
%     fprintf('\nSaved:')
%     disp(figure_name) 
%     results = true;
% elseif autoSave == false
%     switch questdlg('Save Image?', 'Figure', 'Save Figure', 'Close Figure', 'Cancel', 'Save Figure')
%         case 'Save Figure'
%             export_fig(fig_handle, [figure_name '.' type(2:end)], type, '-nocrop', fig_quality , '-painters', '-rgb');
%             if closeFig
%                 close(fig_handle)
%             end
%             fprintf('\nSaved:')
%             disp(figure_name) 
%             results = true;
%         case 'Close Figure'
%             if closeFig
%                 close(fig_handle)
%             end
%             results = false;
%         case 'Cancel'
%             results = 'Cancel';
%     end
% end

% ------------------------------------------------------------------------------------------------------------

