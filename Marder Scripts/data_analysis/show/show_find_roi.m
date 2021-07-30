function roi_index=f(h)

% if h doesn't represent an roi decoration, 0 is returned

% get handles we'll need
ifig_h=gcbf;
colorbar_axes_h=findobj(ifig_h,'Tag','colorbar_axes_h');
border_h=get_userdata(colorbar_axes_h,'border_h');
label_h=get_userdata(colorbar_axes_h,'label_h');
n_rois=get_userdata(ifig_h,'n_rois');

% return the ROI index
if n_rois>0
  new_selected_roi=find(border_h==h|label_h==h);
  if isempty(new_selected_roi)
    roi_index=0;
  else
    roi_index=new_selected_roi(1);
  end
else
  roi_index=0;
end