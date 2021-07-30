function f(varargin)

switch(get(gcbo,'Tag'))
  case 'FPS_edit_h',
    new_fps_string=get(gcbo,'String');
    new_fps=str2num(new_fps_string);
    if (new_fps>0)
      set_userdata(gcbf,'fps',new_fps);
    end
    set(gcbo,'String',sprintf('%6.2f',get_userdata(gcbf,'fps')));
  case 'frame_index_edit_h',
    new_frame_index_string=get(gcbo,'String');
    new_frame_index=str2num(new_frame_index_string);
    n_frames=size(get(findobj(gcbf,'Tag','axes_h'),'UserData'),4);
    if (new_frame_index>=1) & (new_frame_index<=n_frames)
      vcr_change_frame(new_frame_index);
    else
      set(gcbo,'String',sprintf('%d',get_userdata(gcbf,'frame_index')));
    end
  case 'to_start_button_h',
    vcr_change_frame(1);
  case 'play_backward_button_h',
    vcr_play(-1);    
  case 'frame_backward_button_h',
    vcr_change_frame(get_userdata(gcbf,'frame_index')-1);    
  case 'stop_button_h',
    set_userdata(gcbf,'stop_button_hit',1);
  case 'frame_forward_button_h',
    vcr_change_frame(get_userdata(gcbf,'frame_index')+1);        
  case 'play_forward_button_h',
    vcr_play(+1);    
  case 'to_end_button_h',
    axes_h=findobj(gcbf,'Tag','axes_h');
    vcr_change_frame(size(get(axes_h,'UserData'),4));  
  case 'copy_as_bitmap_menu_h',
    print -dbitmap;
end  

