function f(new_plot_x_axis)

% get handles of figs
pfig_h=gcbf;

% get the current units
old_plot_x_axis=get_userdata(pfig_h,'plot_x_axis');

% set the figure var
set_userdata(pfig_h,'plot_x_axis',new_plot_x_axis);

% uncheck the old menu item
menu_h=findobj(pfig_h,'Tag',sprintf('%s_menu_h',old_plot_x_axis));
set(menu_h,'Checked','off');

% check the new menu item
menu_h=findobj(pfig_h,'Tag',sprintf('%s_menu_h',new_plot_x_axis));
set(menu_h,'Checked','on');

% update the plot
chart_update;

