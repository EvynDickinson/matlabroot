function f()

%fprintf(1,'key\n');

% handle the key
key=get(gcbf,'CurrentCharacter');
switch key
  case { ',' , '<' }
    frame_index=get_userdata(gcbf,'frame_index');
    show_change_frame(frame_index-1);
  case { '.' , '>' }
    frame_index=get_userdata(gcbf,'frame_index');
    show_change_frame(frame_index+1);
  case 'p'
    show_play(+1);
  case { char(8) , char(127) }  
    selected_roi_index=get_userdata(gcbf,'selected_roi_index');
    if ~isnan(selected_roi_index)
      show_delete_rois(selected_roi_index);
    end
end
