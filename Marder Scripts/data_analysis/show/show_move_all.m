function f(action)

persistent ifig_h;
persistent image_axes_h;
persistent n_rois;
persistent roi_borders;
persistent roi_border_line_hs;
persistent roi_label_hs;
persistent anchor;

switch(action)
  case 'start'
    % set up the persistents
    ifig_h=gcbf;
    image_axes_h=findobj(ifig_h,'Tag','image_axes_h');
    colorbar_axes_h=findobj(ifig_h,'Tag','colorbar_axes_h');
    roi_border_line_hs=get_userdata(colorbar_axes_h,'border_h');
    n_rois=length(roi_border_line_hs);
    roi_borders=cell(n_rois,1);
    for i=1:n_rois
      x=get(roi_border_line_hs(i),'XData');
      y=get(roi_border_line_hs(i),'YData');
      this_roi_border=zeros(2,length(x));
      this_roi_border(1,:)=x;
      this_roi_border(2,:)=y;
      roi_borders{i}=this_roi_border;
    end
    % hide the ROI labels
    roi_label_hs=get_userdata(colorbar_axes_h,'label_h');
    for i=1:n_rois
      set(roi_label_hs(i),'Visible','off');
    end    
    % now do real stuff    
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    anchor=point;
    % set the callbacks for the drag
    set(ifig_h,'WindowButtonMotionFcn',...
               'show_move_all(''move'')');
    set(ifig_h,'WindowButtonUpFcn',...
               'show_move_all(''stop'')');
  case 'move'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    vect=point-anchor;
    % move all the borders & labels
    for i=1:n_rois
      this_roi_border=roi_borders{i};
      set(roi_border_line_hs(i),'XData',vect(1)+this_roi_border(1,:));
      set(roi_border_line_hs(i),'YData',vect(2)+this_roi_border(2,:));
    end
  case 'stop'
    % change the move and buttonup calbacks
    set(ifig_h,'WindowButtonMotionFcn','show_update_pointer');
    set(ifig_h,'WindowButtonUpFcn','');
    % now do the stuff we'd do for a move also
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    vect=point-anchor;
    % move all the borders & labels
    for i=1:n_rois
      this_roi_border=roi_borders{i};
      set(roi_border_line_hs(i),'XData',vect(1)+this_roi_border(1,:));
      set(roi_border_line_hs(i),'YData',vect(2)+this_roi_border(2,:));
    end
    % move the labels, make visible again
    for i=1:n_rois
      this_position=get(roi_label_hs(i),'Position');
      set(roi_label_hs(i),'Position',[vect 0]+this_position);
      set(roi_label_hs(i),'Visible','on');
    end        
    % clear the persistents
    ifig_h=[];
    image_axes_h=[];
    n_rois=[];
    roi_borders=[];
    roi_border_line_hs=[];
    roi_label_hs=[];
    anchor=[];
    % since the borders are stored in the lines, we're done
end  % switch






