function f(tl_new)

% get handle of fig
fig_h=gcbf;
axes_hs=get_userdata(fig_h,'axes_hs');
tl=xlim(axes_hs(1));
x_axis_menu_h=findobj(fig_h,'tag','x_axis_menu_h');
t=get_userdata(x_axis_menu_h,'t');
t0=t(1);  tf=t(end);
if tl_new(1)>tl_new(2)
  tl_new=fliplr(tl_new);
end
if tl_new(1)<t0
  tl_new(1)=t0;
end
if tl_new(2)>tf
  tl_new(2)=tf;
end
if tl_new(1)~=tl_new(2)
  set(axes_hs,'xlim',tl_new);
end
