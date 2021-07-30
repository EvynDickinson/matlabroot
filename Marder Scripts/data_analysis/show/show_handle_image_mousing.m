function f()

% get handles of figs
ifig_h=gcbf;

% get the current mode
mode=get_userdata(ifig_h,'mode');
% switch on the current roi mode appropriately
switch(mode)
  case 'rect_roi'
    show_draw_rect_roi('start');
  case 'elliptic_roi'
    show_draw_elliptic_roi('start');
  case 'edit_roi'
    % do nothing
    %show_edit_roi;
  case 'select'
    roi_index=show_find_roi(gcbo);
    if roi_index==0
      show_unselect_selected_roi;
    else
      show_select_roi(roi_index);
      show_move_roi('start',roi_index);
    end
  case 'zoom'
    sel_type=get(gcbf,'SelectionType');
    switch sel_type
      case 'extend'
        show_zoom_out;
      case 'alternate'
        show_zoom_out;
      case 'open'
        show_zoom_out;
      otherwise
        show_draw_zoom_rect('start');
    end
  case 'move_all'
    show_move_all('start');
end

