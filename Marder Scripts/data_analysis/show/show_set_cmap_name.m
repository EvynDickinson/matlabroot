function f(new_cmap_name)

% get handles of figs
ifig_h=gcbf;

% need to remember the old cmap_name so that we can
% uncheck that menu item
old_cmap_name=get_userdata(ifig_h,'cmap_name');

% set the chosen cmap_name
set_userdata(ifig_h,'cmap_name',new_cmap_name);

% uncheck the old menu item
menu_h=...
  findobj(ifig_h,'Tag',sprintf('%s_menu_h',old_cmap_name));
set(menu_h,'Checked','off');

% check the new menu item
menu_h=...
  findobj(ifig_h,'Tag',sprintf('%s_menu_h',new_cmap_name));
set(menu_h,'Checked','on');

% set the colormap
eval(sprintf('set(ifig_h,''Colormap'',%s(256));',new_cmap_name));
