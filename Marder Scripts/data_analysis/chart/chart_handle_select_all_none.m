function f(all_or_none)

% all_or_none==0 means none
% all_or_none==1 means all

% get stuff we need
pfig_h=gcbf;
roi_listbox_h=findobj(pfig_h,'Tag','roi_listbox_h');
n_rois=get_userdata(pfig_h,'n_rois');

% do it
trace_menu_h=get_userdata(pfig_h,'trace_menu_h');
trace_on=get_userdata(pfig_h,'trace_on');
if all_or_none
  % all
  % leave the scrolling alone  
  listboxtop=get(roi_listbox_h,'ListboxTop');
  % now turn all traces on
  trace_on(:)=1;
  set(roi_listbox_h,'Value',1:n_rois);
  % leave the scrolling alone
  set(roi_listbox_h,'ListboxTop',listboxtop);
else
  % none
  trace_on(:)=0;
  set(roi_listbox_h,'Value',[]);
end

% commit the change
set_userdata(pfig_h,'trace_on',trace_on);

% update the plot
chart_update;
