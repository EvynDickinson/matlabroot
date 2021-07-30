function f(new_debleach_type)

% get handles of figs
pfig_h=gcbf;

% get the current debleaching type
old_debleach_type=get_userdata(pfig_h,'debleach_type');

% uncheck the old menu item
switch (old_debleach_type)
  case 'none'
    menu_h=findobj(pfig_h,'Tag','no_debleach_menu_h');
  case '1'
    menu_h=findobj(pfig_h,'Tag','linear_debleach_menu_h');
  case '2'
    menu_h=findobj(pfig_h,'Tag','quadratic_debleach_menu_h');
  case '3'
    menu_h=findobj(pfig_h,'Tag','cubic_debleach_menu_h');
  case '4'
    menu_h=findobj(pfig_h,'Tag','quartic_debleach_menu_h');
  case '5'
    menu_h=findobj(pfig_h,'Tag','quintic_debleach_menu_h');
  case 'lpf'
    menu_h=findobj(pfig_h,'Tag','lpf_debleach_menu_h');
end        
set(menu_h,'Checked','off');

% check the new menu item
switch (new_debleach_type)
  case 'none'
    menu_h=findobj(pfig_h,'Tag','no_debleach_menu_h');
  case '1'
    menu_h=findobj(pfig_h,'Tag','linear_debleach_menu_h');
  case '2'
    menu_h=findobj(pfig_h,'Tag','quadratic_debleach_menu_h');
  case '3'
    menu_h=findobj(pfig_h,'Tag','cubic_debleach_menu_h');
  case '4'
    menu_h=findobj(pfig_h,'Tag','quartic_debleach_menu_h');
  case '5'
    menu_h=findobj(pfig_h,'Tag','quintic_debleach_menu_h');
  case 'lpf'
    menu_h=findobj(pfig_h,'Tag','lpf_debleach_menu_h');
end        
set(menu_h,'Checked','on');

% set the figure var
set_userdata(pfig_h,'debleach_type',new_debleach_type);

% update the plot
chart_update;








