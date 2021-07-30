function f(serial_number,varargin)

% create a color sequence if this is the first run
persistent color_sequence;
if isempty(color_sequence)
  % this creates a sequence of colors to be used for traces
  % the sequence is specified in the HSV colorspace
  % all traces have the same satuaration and brightness
  % the hue sequence is determined by taking 256 evenly spaced samples
  %   from 0 to 1, and then shuffling them in such a way that each
  %   hue gets mapped to the 'bit-reversed' hue.
  % This means that subsequent colors tend to be far apart in hue, 
  %   which is the desired effect.
  n_colors=256;
  saturation=1;
  brightness=0.9;
  indices=uint8((0:(n_colors-1))');
  for j=1:n_colors
    indices(j)=bit_reverse(indices(j));
  end
  color_sequence_hsv=...
    [ double(indices)/n_colors ...
      repmat(saturation,[n_colors 1]) ...
      repmat(brightness,[n_colors 1]) ];
  color_sequence=hsv2rgb(color_sequence_hsv);
  %figure;
  %colormap(color_sequence);
  %colorbar;
end
  
% handle args
if nargin>1
  t_o=varargin{1};
end
if nargin>2
  roi_sums=varargin{2};
end
if nargin>3
  roi_n_pels=varargin{3};
end
if nargin>4
  roi_ids=varargin{4};
end
if nargin>5
  roi_labels=varargin{5};
end
if nargin>6
  t_e=varargin{6};
end
if nargin>7
  e_phys=varargin{7};
end

% is gcbf a chart_figure_h, or does a chart_figure_h exist with the 
% given serial number?
chart_figure_h=NaN;
if strcmp(get(gcbf,'Tag'),'chart_figure_h')
  % gcbf _is_ a chart
  chart_figure_h=gcbf;
else
  % gcbf is _not_ a chart, ergo need to find a chart w/ same serial 
  %   number
  hs=findobj('Tag','chart_figure_h');
  for j=1:length(hs)
    this_guys_serial_number=get_userdata(hs(j),'serial_number');
    if this_guys_serial_number==serial_number
      chart_figure_h=hs(j);
      break;
    end
  end
end

% if not, create one and proceed
if isnan(chart_figure_h)
  chart_figure_h=...
    chart_create(t_o,roi_sums,roi_n_pels,...
                 roi_ids,roi_labels,...
                 t_e,e_phys,serial_number);
end

% if we have non-SN args, replace the current chart data w/ them and 
% proceed
if nargin>1
  % filter out ROIs w/ zero pels
  roi_sums=roi_sums(:,roi_n_pels>0);
  roi_ids=roi_ids(roi_n_pels>0);
  roi_labels={roi_labels{roi_n_pels>0}};
  roi_n_pels=roi_n_pels(roi_n_pels>0);
  % figure out how old rois map to new ones, as best we can
  old_n_rois=get_userdata(chart_figure_h,'n_rois');
  old_roi_ids=get_userdata(chart_figure_h,'roi_ids');
  old_trace_on=get_userdata(chart_figure_h,'trace_on');
  n_rois=length(roi_ids);
  trace_on=repmat(0,[n_rois 1]);  % default is off
  if n_rois>0 & old_n_rois>0
    matches=...
      repmat(roi_ids,[1 old_n_rois])==repmat(old_roi_ids',[n_rois 1]);
    [r,c]=find(matches);
    trace_on(r)=old_trace_on(c);
  end
  % if none are on, turn first one on
  if n_rois>0 & sum(trace_on)==0
    trace_on(1)=1;
  end
  % commit the new roi info to the figure state 
  set_userdata(chart_figure_h,'n_rois',n_rois);
  set_userdata(chart_figure_h,'roi_ids',roi_ids);
  set_userdata(chart_figure_h,'roi_labels',roi_labels);
  set_userdata(chart_figure_h,'trace_on',trace_on);
  % commit the actual roi data to the figure state
  data_glom=struct('t_o',t_o,...
                   'roi_sums',roi_sums,...
                   'roi_n_pels',roi_n_pels);
  optical_axes_h=findobj(chart_figure_h,'Tag','optical_axes_h');
  set(optical_axes_h,'UserData',data_glom);
  % commit the e-phys data
  data_glom=struct('t_e',t_e,'e_phys',e_phys);
  xaxis_menu_h=findobj(chart_figure_h,'Tag','xaxis_menu_h');
  set(xaxis_menu_h,'UserData',data_glom);
  clear data_glom;
  % put the roi labels into the ROI listbox
  roi_listbox_h=findobj(chart_figure_h,'Tag','roi_listbox_h');
  set(roi_listbox_h,'String',roi_labels);
  set(roi_listbox_h,'Value',find(trace_on));
  one_max_menu_h=findobj(chart_figure_h,'Tag','one_max_menu_h');
end

% rescale the roi_sums
[x_o,x_e,y]=chart_scale_traces(chart_figure_h);
n_frames=length(x_o);
n_samples=length(x_e);
n_rois=size(y,2);

% get the e-phys data
xaxis_menu_h=findobj(chart_figure_h,'Tag','xaxis_menu_h');
glom=get(xaxis_menu_h,'UserData');
t_e=glom.t_e;
e_phys=glom.e_phys;
clear glom;
n_signals=size(e_phys,2);
% the data in e_phys break down like this:
% v1 
% i1 
% v2 
% i2 
% e1 
% e2
% notscan
% notready
% axoclamp_command
% shutter_command



%
% plot the optical traces
%

% write optical data into the traces
trace_h=get_userdata(chart_figure_h,'trace_h');
delete(trace_h);
trace_on=get_userdata(chart_figure_h,'trace_on');
roi_ids=get_userdata(chart_figure_h,'roi_ids');
optical_axes_h=findobj(chart_figure_h,'Tag','optical_axes_h');
trace_h=zeros(n_rois,1);
for j=1:n_rois
  if trace_on(j)
    trace_h(j)=line('Parent',optical_axes_h,...
                    'XData',x_o,...
                    'YData',y(:,j),...
                    'Color',color_sequence(roi_ids(j),:),...
                    'Visible','on');
  else
    trace_h(j)=line('Parent',optical_axes_h,...
                    'Visible','off');
  end
end
set_userdata(chart_figure_h,'trace_h',trace_h);

% deal with autoscaling
autoscale=get_userdata(chart_figure_h,'autoscale');
n_traces_on=sum(trace_on);
if autoscale & n_traces_on>0
  y_on=y(:,logical(trace_on));  
  y_lim=tight_y_limits(y_on(:));
  y_max_string=sprintf('%.4e',y_lim(2));
  y_min_string=sprintf('%.4e',y_lim(1));
  y_max=str2num(y_max_string);
  y_min=str2num(y_min_string);
  optical_axes_h=findobj(chart_figure_h,'Tag','optical_axes_h');
  set(optical_axes_h,'YLim',[y_min y_max]);
  set_userdata(chart_figure_h,'y_min_string',y_min_string);  
  set_userdata(chart_figure_h,'y_max_string',y_max_string);
end

% change the optical axes y axis label
axes(optical_axes_h);
plot_y_axis=get_userdata(chart_figure_h,'plot_y_axis');
switch(plot_y_axis)
  case 'mean'
    ylabel('Mean');
  case 'mean_ac'
    ylabel('Mean AC');
  case 'sum'
    ylabel('Sum');
  case 'sum_ac'
    ylabel('Sum AC');
  case 'dff'
    ylabel('dF/F (%)');
end

% redraw the optical traces legend
%legend(optical_axes_h,'off');
%if n_traces_on>0
%  roi_id_strings=cell(n_rois,1);
%  for j=1:n_rois
%    roi_id_strings{j}=num2str(roi_ids(j));
%  end
%  legend(optical_axes_h,...
%         trace_h(logical(trace_on)),...
%         roi_id_strings(logical(trace_on)));
%end



%
% plot the e-phys traces
%

% get the axes handles
voltage_axes_h=findobj(chart_figure_h,'Tag','voltage_axes_h');
current_axes_h=findobj(chart_figure_h,'Tag','current_axes_h');
extra_axes_h=findobj(chart_figure_h,'Tag','extra_axes_h');
ttl_axes_h=findobj(chart_figure_h,'Tag','ttl_axes_h');

% set the axes and color for each e-phys trace
axes_h=[ voltage_axes_h ; ...
         current_axes_h ; ...
         voltage_axes_h ; ...
         current_axes_h ; ...
         extra_axes_h ; ...
         extra_axes_h ; ...
         ttl_axes_h ; ...
         ttl_axes_h ; ...
         ttl_axes_h ; ...
         ttl_axes_h ];
e_phys_colors=[ 0 0 1 ;
                0 0 1 ;
                0 .5 0 ;
                0 .5 0 ;
                0 0 1 ;
                0 .5 0 ;
                0 0 1 ;
                0 .5 0 ;
                0 0 1 ;
                0 .75 .75 ];
e_phys_trace_on=[ 1 1 1 1 1 1 0 0 1 1]';

% write the data into the traces
e_phys_trace_h=get_userdata(chart_figure_h,'e_phys_trace_h');
delete(e_phys_trace_h);
e_phys_trace_h=zeros(n_signals,1);
for j=n_signals:-1:1
  if e_phys_trace_on(j)
    e_phys_trace_h(j)=line('Parent',axes_h(j),...
                           'XData',x_e,...
                           'YData',e_phys(:,j),...
                           'Color',e_phys_colors(j,:));
  else
    e_phys_trace_h(j)=line('Parent',axes_h(j),...
                           'Visible','off');
  end                       
end

% commit the e-phys traces
set_userdata(chart_figure_h,'e_phys_trace_h',e_phys_trace_h);

% set the y limits for each axes
if ~isempty(e_phys)
  set(voltage_axes_h,'YLim',tight_y_limits([e_phys(:,1) e_phys(:,3)]));
  set(current_axes_h,'YLim',tight_y_limits([e_phys(:,2) e_phys(:,4)]));
  set(extra_axes_h,'YLim',tight_y_limits(e_phys(:,5:6)));
end
set(ttl_axes_h,'YLim',[-1 6],...
               'YTick',[0 5]);

% set the y axis label for each e-phys plot
axes(voltage_axes_h); ylabel('Voltage (mV)');
axes(current_axes_h); ylabel('Current (nA)');
axes(extra_axes_h); ylabel('Extra (uV)');
axes(ttl_axes_h); ylabel('TTL (V)');


           
%
% X axis stuff for all plots
%

% determine the x axis limits
if isempty(x_e)
  x_min=x_o(1);
  x_max=x_o(n_frames);
else
  x_min=min(x_o(1),x_e(1));
  x_max=max(x_o(n_frames),x_e(n_samples));
end

% change the x axis limits
set(optical_axes_h,'XLim',[x_min x_max]);
set(voltage_axes_h,'XLim',[x_min x_max]);
set(current_axes_h,'XLim',[x_min x_max]);
set(extra_axes_h,'XLim',[x_min x_max]);
set(ttl_axes_h,'XLim',[x_min x_max]);

% if plot_x_axis==frame number, make sure first frame and last frame 
%   have tick mark
plot_x_axis=get_userdata(chart_figure_h,'plot_x_axis');
%if strcmp(plot_x_axis,'frame_number')
%  x_tick=get(optical_axes_h,'XTick');
%  x_tick=x_tick(x_tick>=1 & x_tick<=n_frames);
%  if x_tick(1)~=1
%    x_tick=[1 x_tick];
%  end
%  if x_tick(length(x_tick))~=n_frames)
%    x_tick=[x_tick n_frames];
%  end
%  set(optical_axes_h,'XTick',x_tick);        
%end

% change the x axis label
axes(ttl_axes_h);
switch(plot_x_axis)
  case 'frame_number'
    xlabel('Frame Number');
  case 'time_sec'
    xlabel('Time (s)');
  case 'time_msec'
    xlabel('Time (ms)');
end

% other crap
next_trace_menu_h=...
  findobj(chart_figure_h,'Tag','next_trace_menu_h');
previous_trace_menu_h=...
  findobj(chart_figure_h,'Tag','previous_trace_menu_h');
if n_traces_on==1
  set(next_trace_menu_h,'Enable','on');
  set(previous_trace_menu_h,'Enable','on');
else
  set(next_trace_menu_h,'Enable','off');
  set(previous_trace_menu_h,'Enable','off');
end
