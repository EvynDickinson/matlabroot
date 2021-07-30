function f()

% get whos selected
pfig_h=gcbf;
old_trace_on=get_userdata(pfig_h,'trace_on');
new_trace_on_indices=get(gcbo,'Value');
new_trace_on=zeros(size(old_trace_on));
new_trace_on(new_trace_on_indices)=1;

% commit
set_userdata(pfig_h,'trace_on',new_trace_on);

% update the plot
chart_update;

