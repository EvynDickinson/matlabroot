function f(di)

% di is +1 for next trace
% di is -1 for prev trace

% get state vars we need
pfig_h=gcbf;
n_rois=get_userdata(pfig_h,'n_rois');
roi_listbox_h=findobj(pfig_h,'Tag','roi_listbox_h');
trace_on=get_userdata(pfig_h,'trace_on');

% figure out which trace to turn off, which to turn on, do it
old=find(trace_on);
if length(old)==1
  % inc/dec to get new trace
  new=old+di;
  % wrap
  if new<1
    new=n_rois;
  end
  if new>n_rois
    new=1;
  end
  % set the state
  trace_on(old)=0;
  trace_on(new)=1;
  % change the listbox
  set(roi_listbox_h,'Value',new);
  % commit the change to state
  set_userdata(pfig_h,'trace_on',trace_on);
  % update the plot
  chart_update;
end