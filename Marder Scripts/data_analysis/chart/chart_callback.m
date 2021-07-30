function f(varargin)

% get handle of fig
pfig_h=gcbf;

% the big switch
tag=get(gcbo,'Tag');
%fprintf(1,sprintf('%s\n',tag));
switch(tag)
  % x-axis menu items
  case 'frame_number_menu_h'
    chart_change_plot_x_axis('frame_number');
  case 'time_sec_menu_h'
    chart_change_plot_x_axis('time_sec');        
  case 'time_msec_menu_h'
    chart_change_plot_x_axis('time_msec');        
  case 'frequency_hz_menu_h'
    chart_change_plot_x_axis('frequency_hz');        
  % y-axis menu items
  case 'autoscale_menu_h'
    old_autoscale=get_userdata(pfig_h,'autoscale');
    new_autoscale=~old_autoscale;
    set_userdata(pfig_h,'autoscale',new_autoscale);
    if (new_autoscale)
      set(gcbo,'Checked','on');
      chart_update;
    else
      set(gcbo,'Checked','off');
    end
  case 'edit_bounds_menu_h'
    % get the current y min and y max strings
    y_min_string=get_userdata(pfig_h,'y_min_string');
    y_max_string=get_userdata(pfig_h,'y_max_string');    
    % throw up the dialog box
    bounds=inputdlg({ 'Y Axis Maximum:' , 'Y Axis Minimum:' },...
                    'Edit Bounds...',...
                    1,...
                    { y_max_string , y_min_string },...
                    'off');
    if isempty(bounds) 
      break; 
    end
    % break out the returned cell array                
    new_y_max_string=bounds{1};
    new_y_min_string=bounds{2};
    % convert all these strings to real numbers
    y_min=str2num(y_min_string);
    y_max=str2num(y_max_string);
    new_y_min=str2num(new_y_min_string);
    new_y_max=str2num(new_y_max_string);
    % if new values are kosher, change plot limits
    if ~isempty(new_y_min) & ~isempty(new_y_max) & ...
       isfinite(new_y_min) & isfinite(new_y_max) & ...
       (new_y_max>new_y_min)
      set_userdata(pfig_h,'y_min_string',new_y_min_string);
      set_userdata(pfig_h,'y_max_string',new_y_max_string);
      optical_axes_h=findobj(pfig_h,'Tag','optical_axes_h');
      set(optical_axes_h,'YLim',[new_y_min new_y_max]);
    end
  case 'mean_menu_h'  
    chart_change_plot_y_axis('mean');
  case 'mean_ac_menu_h'  
    chart_change_plot_y_axis('mean_ac');
  case 'sum_menu_h'  
    chart_change_plot_y_axis('sum');
  case 'sum_ac_menu_h'  
    chart_change_plot_y_axis('sum_ac');
  case 'dff_menu_h'
    chart_change_plot_y_axis('dff');    
  % debleaching menu items
  case 'no_debleach_menu_h'
    chart_change_debleach('none');
  case 'linear_debleach_menu_h'
    chart_change_debleach('1');
  case 'quadratic_debleach_menu_h'
    chart_change_debleach('2');
  case 'cubic_debleach_menu_h'
    chart_change_debleach('3');
  case 'quartic_debleach_menu_h'
    chart_change_debleach('4');
  case 'quintic_debleach_menu_h'
    chart_change_debleach('5');
  case 'lpf_debleach_menu_h'
    chart_change_debleach('lpf');
  % traces menu items
  case 'select_all_menu_h'
    chart_handle_select_all_none(1);
  case 'select_none_menu_h'
    chart_handle_select_all_none(0);
  case 'next_trace_menu_h'  
    chart_handle_step_trace(+1);
  case 'previous_trace_menu_h'  
    chart_handle_step_trace(-1);    
  % print menu items
  case 'page_setup_menu_h'
    pagesetupdlg(gcbf);
  case 'print_setup_menu_h'
    % this suffers from a MATLAB bug
    eval(sprintf('print -dsetup -f%.16f',gcbf));
  case 'print_preview_menu_h'
    printpreview(gcbf);
  case 'print_menu_h'
    eval(sprintf('print -v -f%.16f',gcbf));
  % ROI listbox
  case 'roi_listbox_h'
    chart_handle_roi_listbox_action;
end  



