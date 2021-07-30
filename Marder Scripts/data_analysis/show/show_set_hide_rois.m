function f(new_hide_rois)

if new_hide_rois
  % hide them
  colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
  label_h=get_userdata(colorbar_axes_h,'label_h');
  border_h=get_userdata(colorbar_axes_h,'border_h');
  for j=1:length(label_h)
    set(label_h(j),'Visible','off');
    set(border_h(j),'Visible','off');
  end
  set_userdata(gcbf,'hide_rois',1);
  hide_rois_menu_h=findobj(gcbf,'Tag','hide_rois_menu_h');
  set(hide_rois_menu_h,'Label','Show ROIs');
  % set the image erase mode
  image_h=findobj(gcbf,'Tag','image_h');
  set(image_h,'EraseMode','none');
else
  % show them
  colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
  label_h=get_userdata(colorbar_axes_h,'label_h');
  border_h=get_userdata(colorbar_axes_h,'border_h');
  for j=1:length(label_h)
    set(label_h(j),'Visible','on');
    set(border_h(j),'Visible','on');
  end
  set_userdata(gcbf,'hide_rois',0);
  hide_rois_menu_h=findobj(gcbf,'Tag','hide_rois_menu_h');
  set(hide_rois_menu_h,'Label','Hide ROIs');
  % set the image erase mode
  image_h=findobj(gcbf,'Tag','image_h');
  set(image_h,'EraseMode','normal');
end
