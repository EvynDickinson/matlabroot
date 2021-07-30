function f(varargin)

% the big switch
tag=get(gcbo,'Tag');
switch(tag)
  % the text-entry boxes
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
    image_axes_h=findobj(gcbf,'Tag','image_axes_h');
    n_frames=size(get(image_axes_h,'UserData'),3);
    if (new_frame_index>=1) & (new_frame_index<=n_frames)
      show_change_frame(new_frame_index);
    else
      set(gcbo,'String',sprintf('%d',get_userdata(gcbf,'frame_index')));
    end
  % VCR buttons
  case 'to_start_button_h',
    show_change_frame(1);
  case 'play_backward_button_h',
    show_play(-1);    
  case 'frame_backward_button_h',
    show_change_frame(get_userdata(gcbf,'frame_index')-1);    
  case 'stop_button_h',
    set_userdata(gcbf,'stop_button_hit',1);
  case 'frame_forward_button_h',
    show_change_frame(get_userdata(gcbf,'frame_index')+1);        
  case 'play_forward_button_h',
    show_play(+1);    
  case 'to_end_button_h',
    image_axes_h=findobj(gcbf,'Tag','image_axes_h');
    n_frames=size(get(image_axes_h,'UserData'),3);
    show_change_frame(n_frames);
  % image and decorations
  case 'image_h'
    show_handle_image_mousing;
  case 'border_h'  
    show_handle_image_mousing;
  case 'label_h'  
    show_handle_image_mousing;
  % mode buttons  
  case 'rect_roi_button_h'    
    show_set_mode('rect_roi');
  case 'elliptic_roi_button_h'
    show_set_mode('elliptic_roi');
  case 'edit_roi_button_h'
    show_set_mode('edit_roi');
  case 'select_button_h'    
    show_set_mode('select');
  case 'zoom_button_h'
    show_set_mode('zoom');
  case 'move_all_button_h'
    show_set_mode('move_all');
  % action buttons
  case 'chart_button_h'
    serial_number=get_userdata(gcbf,'serial_number');
    [roi_sums,roi_n_pels]=show_calc_roi_sums;
    colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
    roi_ids=get_userdata(colorbar_axes_h,'roi_ids');
    colorbar_h=findobj(gcbf,'Tag','colorbar_h');
    glom=get(colorbar_h,'UserData');
    roi_label_hs=get_userdata(colorbar_axes_h,'label_h');
    n_rois=length(roi_label_hs);
    roi_labels=cell(n_rois,1);
    for i=1:n_rois
      roi_labels{i}=get(roi_label_hs(i),'String');
    end
    chart_update(serial_number,...
                 glom.t_o,roi_sums,roi_n_pels,...
                 roi_ids,roi_labels,...
                 glom.t_e,glom.e_phys);
  case 'polka_button_h'
    serial_number=get_userdata(gcbf,'serial_number');
    [roi_sums,roi_n_pels]=show_calc_roi_sums;
    colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
    roi_ids=get_userdata(colorbar_axes_h,'roi_ids');
    colorbar_h=findobj(gcbf,'Tag','colorbar_h');
    glom=get(colorbar_h,'UserData');
    roi_label_hs=get_userdata(colorbar_axes_h,'label_h');
    n_rois=length(roi_label_hs);
    roi_labels=cell(n_rois,1);
    for i=1:n_rois
      roi_labels{i}=get(roi_label_hs(i),'String');
    end
    polka_update(serial_number,...
                 glom.t_o,roi_sums,roi_n_pels,...
                 roi_ids,roi_labels,...
                 glom.t_e,glom.e_phys);
        
  %
  % Color menu items
  %
  case 'min_max_menu_h'
    show_handle_colorbar_menus(tag); 
  case 'five_95_menu_h'
    show_handle_colorbar_menus(tag); 
  case 'abs_max_menu_h'
    show_handle_colorbar_menus(tag); 
  case 'ninety_symmetric_menu_h'
    show_handle_colorbar_menus(tag); 
  case 'colorbar_edit_bounds_menu_h'
    show_handle_colorbar_menus(tag); 
  case 'gray_menu_h',
    show_set_cmap_name('gray');
  case 'bone_menu_h',
    show_set_cmap_name('bone');
  case 'hot_menu_h',
    show_set_cmap_name('hot');
  case 'jet_menu_h',
    show_set_cmap_name('jet');
  case 'red_green_menu_h'
    show_set_cmap_name('red_green');
  case 'red_blue_menu_h'
   show_set_cmap_name('red_blue');
  case 'brighten_menu_h',
    cmap=get(gcbf,'Colormap');
    new_cmap=brighten(cmap,0.1);
    set(gcbf,'Colormap',new_cmap);
  case 'darken_menu_h',
    cmap=get(gcbf,'Colormap');
    new_cmap=brighten(cmap,-0.1);
    set(gcbf,'Colormap',new_cmap);
  case 'revert_menu_h',
    cmap_name=get_userdata(gcbf,'cmap_name');
    eval(sprintf('set(gcbf,''Colormap'',%s(256));',cmap_name));

  %
  % ROIs menu items
  %
  case 'save_rois_to_file_menu_h'
    show_save_rois_to_file;
  case 'open_rois_menu_h'
    show_load_rois_from_file;    
  case 'rename_roi_menu_h'
    selected_roi_index=get_userdata(gcbf,'selected_roi_index');
    if ~isnan(selected_roi_index)
      show_rename_roi(selected_roi_index);
    end
  case 'delete_roi_menu_h'
    selected_roi_index=get_userdata(gcbf,'selected_roi_index');
    if ~isnan(selected_roi_index)
      show_delete_rois(selected_roi_index);
    end
  case 'delete_all_rois_menu_h'
    n_rois=get_userdata(gcbf,'n_rois');
    roi_list=(1:n_rois)';
    show_delete_rois(roi_list);
  case 'hide_rois_menu_h'
    hide_rois=get_userdata(gcbf,'hide_rois');
    show_set_hide_rois(~hide_rois);

  %
  % Print menu items
  %
  case 'page_setup_menu_h'
    pagesetupdlg(gcbf);
  case 'print_setup_menu_h'
    % this suffers from a MATLAB bug
    eval(sprintf('print -dsetup -f%.16f',gcbf));
  case 'print_preview_menu_h'
    printpreview(gcbf);
  case 'print_menu_h'
    eval(sprintf('print -v -f%.16f',gcbf));
end


