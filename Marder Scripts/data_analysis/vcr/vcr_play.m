function f(direction)
   
figure_h = gcbf;
data_var_name = get_userdata(figure_h,'data_var_name');
axes_h=findobj(figure_h,'Tag','axes_h');
indexed_data=get(axes_h,'UserData');
start_frame_index=get_userdata(figure_h,'frame_index');
n_frames=size(indexed_data,4);
frame_index_edit_h=findobj(figure_h,'Tag','frame_index_edit_h');
image_h = findobj(figure_h,'Tag','image_h');
fps=get_userdata(figure_h,'fps');
spf=1/fps;
if (direction>0)
  frame_sequence=start_frame_index:n_frames;
else
  frame_sequence=start_frame_index:-1:1;
end
%set(frame_index_edit_h,'String',' ');
%drawnow;
set_userdata(figure_h,'stop_button_hit',0);
for frame_index=frame_sequence
  tic;
  this_frame = indexed_data(:,:,:,frame_index);
  set(image_h,'CData',this_frame);
  set(frame_index_edit_h,'String',num2str(frame_index));
  drawnow;
  while (toc < spf)
  end
  if get_userdata(figure_h,'stop_button_hit')
    break;
  end
end
set_userdata(figure_h,'stop_button_hit',0);
set_userdata(figure_h,'frame_index',frame_index);


