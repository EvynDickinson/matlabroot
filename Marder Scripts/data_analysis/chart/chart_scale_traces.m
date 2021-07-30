function [x_o,x_e,y_optical] = f(chart_figure_h)

% get the data
optical_axes_h=findobj(chart_figure_h,'Tag','optical_axes_h');
data_glom=get(optical_axes_h,'UserData');
t_o=data_glom.t_o;
y_optical=data_glom.roi_sums;
n_pels_mask=data_glom.roi_n_pels;
clear data_glom;
xaxis_menu_h=findobj(chart_figure_h,'Tag','xaxis_menu_h');
glom=get(xaxis_menu_h,'UserData');
t_e=glom.t_e;
clear glom;

% get dims
n_frames=size(y_optical,1);
n_rois=size(y_optical,2);

% calculate the time-average of the y_optical
y_optical_time_avg=mean(y_optical,1);

% debleach the data, if necessary
debleach_type=get_userdata(chart_figure_h,'debleach_type');
switch debleach_type
  case 'none'
    % do nothing
  case 'lpf'  
    % lpf debleach w/ 500 ms filter
    lpf_time=500;  % ms
    ts=(t_o(n_frames)-t_o(1))/(n_frames-1);
    for j=1:n_rois
      roi_fast=divide_gauss_trunc(y_optical(:,j),lpf_time/ts);
      y_optical(:,j)=roi_fast*y_optical_time_avg(j);
    end
  otherwise
    % use poly filter
    debleach_order=str2num(debleach_type);
    frame_index=(1:n_frames)';
    for j=1:n_rois
      p=polyfit(frame_index,y_optical(:,j),debleach_order);
      bleaching=polyval(p,frame_index);
      y_optical(:,j)=(y_optical(:,j)./bleaching)*...
                       y_optical_time_avg(j);
    end
end

% change the units of the y_optical data, if necessary
switch get_userdata(chart_figure_h,'plot_y_axis')
  case 'sum'
    % do nothing
  case 'sum_ac'
    y_optical=y_optical-repmat(y_optical_time_avg,[n_frames 1]);
  case 'mean'
    y_optical=y_optical./repmat(n_pels_mask,[n_frames 1]);
  case 'mean_ac'
    y_optical=(y_optical-repmat(y_optical_time_avg,[n_frames 1]))./...
              repmat(n_pels_mask,[n_frames 1]);
  case 'dff'
    y_optical=100*(y_optical./repmat(y_optical_time_avg,[n_frames 1])-1);
end

% figure out what the x axis should be
plot_x_axis=get_userdata(chart_figure_h,'plot_x_axis');
switch(plot_x_axis)
  case 'frame_number'
    x_o=(1:n_frames)';
    x_e=(n_frames-1)/(t_o(n_frames)-t_o(1))*(t_e-t_o(1))+1;
  case 'time_sec'
    x_o=t_o/1000;
    x_e=t_e/1000;
  case 'time_msec'
    x_o=t_o;
    x_e=t_e;
end


