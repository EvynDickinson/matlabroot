function f()

% get handles we'll need
ifig_h=gcbf;
colorbar_axes_h=findobj(ifig_h,'Tag','colorbar_axes_h');
roi_ids=get_userdata(colorbar_axes_h,'roi_ids');
border_h=get_userdata(colorbar_axes_h,'border_h');
label_h=get_userdata(colorbar_axes_h,'label_h');
image_h=findobj(ifig_h,'Tag','image_h');
rename_roi_menu_h=findobj(ifig_h,'Tag','rename_roi_menu_h');
delete_roi_menu_h=findobj(ifig_h,'Tag','delete_roi_menu_h');

% get the old roi
old_selected_roi_index=get_userdata(ifig_h,'selected_roi_index');

if isnan(old_selected_roi_index)
  % nothing to do
else
  % make the old ROI label blue, make it's height 1
  set(label_h(old_selected_roi_index),'Color',[0 0 1]);
  pos=get(label_h(old_selected_roi_index),'Position');
  pos(3)=1;
  set(label_h(old_selected_roi_index),'Position',pos);
  % make the old ROI border blue, make it's height 1
  set(border_h(old_selected_roi_index),'Color',[0 0 1]);
  zdata=get(border_h(old_selected_roi_index),'ZData');
  zdata=repmat(1,size(zdata));
  set(border_h(old_selected_roi_index),'ZData',zdata);
  % these two lines are necessary to prevent the unselected boxes from 
  % disappearing.  I don't know why they disappear, but there you go.
  set(image_h,'Selected','on');
  set(image_h,'Selected','off');  
  % gray the rename, delete roi menus
  set(rename_roi_menu_h,'Enable','off');
  set(delete_roi_menu_h,'Enable','off');
  % no ROI selected anymore
  set_userdata(ifig_h,'selected_roi_index',NaN);
  % disable the edit mode
  edit_roi_button_h=findobj(ifig_h,'Tag','edit_roi_button_h');
  set(edit_roi_button_h,'Enable','off');
end



