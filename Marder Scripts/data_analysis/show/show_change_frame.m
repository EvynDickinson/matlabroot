function f(new_frame_index)

% get state variables
figure_h=gcbf;
image_axes_h=findobj(figure_h,'Tag','image_axes_h');
indexed_data=get(image_axes_h,'UserData');
n_frames=size(indexed_data,3);
image_h = findobj(figure_h,'Tag','image_h');
frame_index_edit_h=findobj(figure_h,'Tag','frame_index_edit_h');
% change things to accord with the new frame index
if (new_frame_index>=1) & (new_frame_index<=n_frames)
  set_userdata(figure_h,'frame_index',new_frame_index);
  set(frame_index_edit_h,'String',sprintf('%d',new_frame_index));
  this_frame = indexed_data(:,:,new_frame_index);
  set(image_h,'CData',this_frame);
end
