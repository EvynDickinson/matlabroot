function f(new_mode)

% get handles of figs
ifig_h=gcbf;

% need to remember the old mode so that we can
% uncheck that menu item
old_mode=get_userdata(ifig_h,'mode');

% untoggle the old mode button
old_button_h=...
  findobj(ifig_h,'Tag',sprintf('%s_button_h',old_mode));
set(old_button_h,'Value',0);

% check the new menu item
new_button_h=...
  findobj(ifig_h,'Tag',sprintf('%s_button_h',new_mode));
set(new_button_h,'value',1);

% set the chosen mode
set_userdata(ifig_h,'mode',new_mode);
