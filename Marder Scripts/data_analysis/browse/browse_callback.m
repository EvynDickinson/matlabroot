function f(varargin)

% get handle of fig
fig_h=gcbf;

% need something special to catch axes events
axes_hs=get_userdata(fig_h,'axes_hs');
if any(gcbo==axes_hs)
  browse_draw_zoom_limits('start',gcbo);
end

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
    old_autoscale=get_userdata(fig_h,'autoscale');
    new_autoscale=~old_autoscale;
    set_userdata(fig_h,'autoscale',new_autoscale);
    if (new_autoscale)
      set(gcbo,'Checked','on');
      chart_update;
    else
      set(gcbo,'Checked','off');
    end
  case 'edit_bounds_menu_h'
    % get the current y min and y max strings
    y_min_string=get_userdata(fig_h,'y_min_string');
    y_max_string=get_userdata(fig_h,'y_max_string');    
    % throw up the dialog box
    bounds=inputdlg({ 'Y Axis Maximum:' , 'Y Axis Minimum:' },...
                    'Edit Bounds...',...
                    1,...
                    { y_max_string , y_min_string },...
                    'off');
    if ~isempty(bounds) 
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
        set_userdata(fig_h,'y_min_string',new_y_min_string);
        set_userdata(fig_h,'y_max_string',new_y_max_string);
        optical_axes_h=findobj(fig_h,'Tag','optical_axes_h');
        set(optical_axes_h,'YLim',[new_y_min new_y_max]);
      end
    end

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
    
  % bottom buttons
  case 'to_start_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    t0=t(1);
    set(axes_hs,'xlim',[t0 t0+tw]);
  case 'page_left_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    t0=t(1);
    tl_new=tl-tw;
    if tl_new(1)<t0
      tl_new=[t0 t0+tw];
    end
    set(axes_hs,'xlim',tl_new);    
  case 'step_left_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    t0=t(1);
    tl_new=tl-0.1*tw;
    if tl_new(1)<t0
      tl_new=[t0 t0+tw];
    end
    set(axes_hs,'xlim',tl_new);    
  case 'step_right_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    tf=t(end);
    tl_new=tl+0.1*tw;
    if tl_new(2)>tf
      tl_new=[tf-tw tf];
    end
    set(axes_hs,'xlim',tl_new);    
  case 'page_right_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    tf=t(end);
    tl_new=tl+tw;
    if tl_new(2)>tf
      tl_new=[tf-tw tf];
    end
    set(axes_hs,'xlim',tl_new);    
  case 'to_end_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    tf=t(end);
    set(axes_hs,'xlim',[tf-tw tf]);
  case 'zoom_way_out_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=[t(1) t(end)];
    set(axes_hs,'xlim',tl);
  case 'zoom_out_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    t0=t(1);  tf=t(end);
    if tl(2)~=tf
      tl_new=[tl(1) tl(1)+2*tw];
      if tl_new(2)>tf
        tl_new(2)=tf;
      end
      set(axes_hs,'xlim',tl_new);    
    elseif tl(1)~=t0
      tl_new=[tf-2*tw tf];
      if tl_new(1)<t0
        tl_new(1)=t0;
      end
      set(axes_hs,'xlim',tl_new);
    end
  case 'zoom_in_button_h'
    axes_hs=get_userdata(fig_h,'axes_hs');
    x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
    t=get_userdata(x_axis_menu_h,'t');
    tl=xlim(axes_hs(1));
    tw=tl(2)-tl(1);
    tl_new=[tl(1) tl(1)+tw/2];
    set(axes_hs,'xlim',tl_new);        
end  



