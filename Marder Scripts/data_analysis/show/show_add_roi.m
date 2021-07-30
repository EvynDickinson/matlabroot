function f(this_line_h)

% get the figure handle
figure_h=gcbf;

% get some handles we'll need
image_axes_h=findobj(figure_h,'Tag','image_axes_h');  
image_h=findobj(figure_h,'Tag','image_h');
colorbar_axes_h=findobj(figure_h,'Tag','colorbar_axes_h');

% what is the ID for the new ROI?
n_rois=get_userdata(figure_h,'n_rois');
roi_struct=get(colorbar_axes_h,'UserData');
if n_rois==0
  this_roi_id=1;
else
  this_roi_id=roi_struct.roi_ids(n_rois)+1;
end

% calc the COM for this ROI
x=get(this_line_h,'XData');
y=get(this_line_h,'YData');
this_border=zeros(2,length(x));
this_border(1,:)=x;
this_border(2,:)=y;
com=border_com(this_border);

% figure out what the label for this ROI will be
counter=this_roi_id;
while 1
  tentative_label=show_id_to_alpha(counter);
  if ~show_label_in_use(tentative_label)
    break;
  end
  counter=counter+1;
end

% make the label for this ROI
this_label_h=...
  text('Parent',image_axes_h,...
       'Position',[com(1) com(2) 2],...
       'String',tentative_label,...
       'HorizontalAlignment','center',...
       'VerticalAlignment','middle',...
       'Color',[0 0 1],...
       'Tag','label_h',...
       'Clipping','on',...         
       'ButtonDownFcn','show_callback');

% change roi_struct
roi_struct.roi_ids(n_rois+1,1)=this_roi_id;
roi_struct.border_h(n_rois+1,1)=this_line_h;
roi_struct.label_h(n_rois+1,1)=this_label_h;

% commit the changes to the figure variables
set_userdata(figure_h,'n_rois',n_rois+1);
set(colorbar_axes_h,'UserData',roi_struct);

% select the just-created ROI
show_select_roi(n_rois+1);

% if this is the first ROI, need to do some stuff
if n_rois==0
  % need to set image erase mode to normal, since now there's something
  % other than the image in that image axes
  set(image_h,'EraseMode','normal');
  % enable the delete ROIS menu item
  delete_all_rois_menu_h=...
    findobj(figure_h,'Tag','delete_all_rois_menu_h');
  set(delete_all_rois_menu_h,'Enable','on');
  % enable the save ROIS to file menu item
  save_rois_to_file_menu_h=...
    findobj(figure_h,'Tag','save_rois_to_file_menu_h');
  set(save_rois_to_file_menu_h,'Enable','on');
  % enable the hide ROIs menu item
  hide_rois_menu_h=...
    findobj(figure_h,'Tag','hide_rois_menu_h');
  set(hide_rois_menu_h,'Enable','on');
  show_set_hide_rois(0);
  % enable the select ROI button
  select_button_h=...
    findobj(figure_h,'Tag','select_button_h');
  set(select_button_h,'Enable','on');
  % enable the move all ROIs button
  move_all_button_h=...
    findobj(figure_h,'Tag','move_all_button_h');
  set(move_all_button_h,'Enable','on');
end  


