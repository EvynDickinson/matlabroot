function f(cb_min_string,cb_max_string)

% get handles of figs
ifig_h=gcbf;

% find the relevant objects
colorbar_axes_h=findobj(ifig_h,'Tag','colorbar_axes_h');
colorbar_h=findobj(ifig_h,'Tag','colorbar_h');
image_axes_h=findobj(ifig_h,'Tag','image_axes_h');
image_h=findobj(ifig_h,'Tag','image_h');

% change the figure strings
set_userdata(ifig_h,'colorbar_min_string',cb_min_string);
set_userdata(ifig_h,'colorbar_max_string',cb_max_string);     

% change the axes and colorbar
cb_min=str2num(cb_min_string);
cb_max=str2num(cb_max_string);
cb_increment=(cb_max-cb_min)/256;
set(colorbar_axes_h,'YLim',[cb_min cb_max]);
set(colorbar_h,'YData',[cb_min+0.5*cb_increment...
                        cb_max-0.5*cb_increment]);

% recalculate indexed_data, set in figure
data=get(image_h,'UserData');
epsilon=1e-6;
new_indexed_data=uint8(floor((256-epsilon)*(double(data)-cb_min)/...
                                           (cb_max-cb_min)));
set(image_axes_h,'UserData',new_indexed_data);

% change the displayed image
frame_index=get_userdata(ifig_h,'frame_index');
this_frame=new_indexed_data(:,:,frame_index);
set(image_h,'CData',this_frame);

