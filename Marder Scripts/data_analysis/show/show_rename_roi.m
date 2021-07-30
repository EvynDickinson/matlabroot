function f(roi_index)

% get some handles we'll need
figure_h=gcbf;
colorbar_axes_h=findobj(figure_h,'Tag','colorbar_axes_h');
roi_label_hs=get_userdata(colorbar_axes_h,'label_h');
this_label_h=roi_label_hs(roi_index);

% get the current label string
this_label_string=get(this_label_h,'String');

% throw up the dialog box
new_label_string=...
  inputdlg({ 'New ROI name:' },...
           'Rename ROI...',...
           1,...
           { this_label_string },...
           'off');
if isempty(new_label_string) 
  break; 
end

% break out the returned cell array                
new_label_string=new_label_string{1};

% if new value is not taken, change label
if ~show_label_in_use(new_label_string)
  set(this_label_h,'String',new_label_string);
  % these lines are necessary to prevent the unselected boxes from 
  % disappearing.  I don't know why they disappear, but there you go.
  image_h=findobj(gcbf,'Tag','image_h');
  set(image_h,'Selected','on');
  set(image_h,'Selected','off');  
end

