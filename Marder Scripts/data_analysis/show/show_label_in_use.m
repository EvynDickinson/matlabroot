function retval=f(str)

% get some handles we'll need
colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
roi_label_hs=get_userdata(colorbar_axes_h,'label_h');

% see if str is already a ROI label
retval=0;
for i=1:length(roi_label_hs)
  if strcmp(get(roi_label_hs(i),'String'),str)
    retval=1;
    break;
  end
end
