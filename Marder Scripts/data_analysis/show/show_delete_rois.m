function f(roi_indices)

% get some handles we'll need
figure_h=gcbf;
image_axes_h=findobj(figure_h,'Tag','image_axes_h');  
image_h=findobj(figure_h,'Tag','image_h');
colorbar_axes_h=findobj(figure_h,'Tag','colorbar_axes_h');

% do it
old_n_rois=get_userdata(figure_h,'n_rois');
old_roi_struct=get(colorbar_axes_h,'UserData');
keep=repmat(1,[old_n_rois 1]);
keep(roi_indices)=0;
keep=logical(keep);
toss=~keep;

% kill widgets
delete(old_roi_struct.border_h(toss));
delete(old_roi_struct.label_h(toss));

% modify roi_struct
roi_struct=struct('roi_ids',old_roi_struct.roi_ids(keep),...
                  'border_h',old_roi_struct.border_h(keep),...
                  'label_h',old_roi_struct.label_h(keep));

% update n_rois
n_rois=old_n_rois-length(roi_indices);

% modify other things
selected_roi_index=get_userdata(figure_h,'selected_roi_index');
if any(roi_indices==selected_roi_index)
  % gray the rename, delete roi menus
  delete_roi_menu_h=findobj(figure_h,'Tag','delete_roi_menu_h');
  set(delete_roi_menu_h,'Enable','off');
  rename_roi_menu_h=findobj(figure_h,'Tag','rename_roi_menu_h');
  set(rename_roi_menu_h,'Enable','off');
  % set the state
  selected_roi_index=NaN;
else
  selected_roi_index=selected_roi_index-...
                     sum(roi_indices<selected_roi_index);
end

% commit changes
set(colorbar_axes_h,'UserData',roi_struct);
set_userdata(figure_h,'n_rois',n_rois);
set_userdata(figure_h,'selected_roi_index',selected_roi_index);
  
% modify ancillary crap  
if n_rois==0
  % need to set image erase mode to none, since now there are no 
  % more lines in front of the image
  set(image_h,'EraseMode','none');
  delete_all_rois_menu_h=...
    findobj(figure_h,'Tag','delete_all_rois_menu_h');
  set(delete_all_rois_menu_h,'Enable','off');
  save_rois_to_file_menu_h=...
    findobj(figure_h,'Tag','save_rois_to_file_menu_h');
  set(save_rois_to_file_menu_h,'Enable','off');  
  hide_rois_menu_h=findobj(figure_h,'Tag','hide_rois_menu_h');
  set(hide_rois_menu_h,'Enable','off');
  show_set_hide_rois(0);
  % disable the select ROI button
  select_button_h=...
    findobj(figure_h,'Tag','select_button_h');
  set(select_button_h,'Enable','off');
  % enable the move all ROIs button
  move_all_button_h=...
    findobj(figure_h,'Tag','move_all_button_h');
  set(move_all_button_h,'Enable','off');  
  % change the mode to elliptic_roi
  show_set_mode('elliptic_roi');
end

