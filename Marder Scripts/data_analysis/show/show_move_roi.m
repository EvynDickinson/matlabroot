function f(action,roi_index)

persistent ifig_h;
persistent image_axes_h;
persistent anchor;
persistent roi_border_h;
persistent roi_label_h;
persistent border_x_anchor;
persistent border_y_anchor;
persistent label_pos_anchor;

switch(action)
  case 'start'
    ifig_h=gcbf;
    image_axes_h=findobj(ifig_h,'Tag','image_axes_h');
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    anchor=point;
    colorbar_axes_h=findobj(ifig_h,'Tag','colorbar_axes_h');
    roi_border_hs=get_userdata(colorbar_axes_h,'border_h');
    roi_border_h=roi_border_hs(roi_index);
    border_x_anchor=get(roi_border_h,'XData');
    border_y_anchor=get(roi_border_h,'YData');
    roi_label_hs=get_userdata(colorbar_axes_h,'label_h');
    roi_label_h=roi_label_hs(roi_index);
    label_pos_anchor=get(roi_label_h,'Position');
    % set the callbacks for the drag
    set(ifig_h,'WindowButtonMotionFcn',...
               'show_move_roi(''move'')');
    set(ifig_h,'WindowButtonUpFcn',...
               'show_move_roi(''stop'')');
  case 'move'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    set(roi_border_h,'XData',border_x_anchor+(point(1)-anchor(1)));
    set(roi_border_h,'YData',border_y_anchor+(point(2)-anchor(2)));
    set(roi_label_h,'Position',label_pos_anchor+[point-anchor 0]);
  case 'stop'
    % change the move and buttonup calbacks
    set(ifig_h,'WindowButtonMotionFcn','show_update_pointer');
    set(ifig_h,'WindowButtonUpFcn','');
    % now do the stuff we'd do for a move also
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    set(roi_border_h,'XData',border_x_anchor+(point(1)-anchor(1)));
    set(roi_border_h,'YData',border_y_anchor+(point(2)-anchor(2)));
    set(roi_label_h,'Position',label_pos_anchor+[point-anchor 0]);
    % clear the persistents
    ifig_h=[];
    image_axes_h=[];
    anchor=[];
    roi_border_h=[];
    roi_label_h=[];
    border_x_anchor=[];
    border_y_anchor=[];
    label_pos_anchor=[];
end  % switch






