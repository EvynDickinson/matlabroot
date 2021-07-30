function f()

figure_h = gcbf;
key=get(figure_h,'CurrentCharacter');
switch key
  case { ',', '<' }
    frame_index=get_userdata(figure_h,'frame_index',1);
    show_change_frame(frame_index-1);
  case { '.', '>' }
    frame_index=get_userdata(figure_h,'frame_index',1);
    show_change_frame(frame_index+1);
  case 'p'
    show_play(+1);
end
