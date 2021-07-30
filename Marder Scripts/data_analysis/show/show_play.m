function f(direction)

% get handles of figs
ifig_h=gcbf;

% play the movie
image_axes_h=findobj(ifig_h,'Tag','image_axes_h');
indexed_data=get(image_axes_h,'UserData');
start_frame_index=get_userdata(ifig_h,'frame_index');
n_frames=size(indexed_data,3);
frame_index_edit_h=findobj(ifig_h,'Tag','frame_index_edit_h');
image_h = findobj(ifig_h,'Tag','image_h');
if (get_userdata(ifig_h,'n_rois')>0)
  set(image_h,'EraseMode','none');
end
fps=get_userdata(ifig_h,'fps');
spf=1/fps;
if (direction>0)
  frame_sequence=start_frame_index:n_frames;
else
  frame_sequence=start_frame_index:-1:1;
end
set_userdata(ifig_h,'stop_button_hit',0);
for frame_index=frame_sequence
  tic;
  this_frame = indexed_data(:,:,frame_index);
  set(image_h,'CData',this_frame);
  set(frame_index_edit_h,'String',num2str(frame_index));
  drawnow;
  while (toc < spf)
  end
  if get_userdata(ifig_h,'stop_button_hit')
    break;
  end
end
set_userdata(ifig_h,'stop_button_hit',0);
set_userdata(ifig_h,'frame_index',frame_index);
if (get_userdata(ifig_h,'n_rois')>0)
  set(image_h,'EraseMode','normal');
end