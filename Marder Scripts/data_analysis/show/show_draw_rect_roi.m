function f(action)

persistent ifig_h;
persistent image_axes_h;
persistent anchor;
persistent rect_h;

switch(action)
  case 'start'
    ifig_h=gcbf;
    image_axes_h=findobj(ifig_h,'Tag','image_axes_h');
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    anchor=point;
    % create a new rectangle
    rect_h=...
      line('Parent',image_axes_h,...
           'Color',[1 0 0],...
           'Tag','border_h',...
           'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)],...
           'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)],...
           'ZData',[2 2 2 2 2],...
           'ButtonDownFcn','show_callback');
    % set the callbacks for the drag
    set(ifig_h,'WindowButtonMotionFcn',...
               'show_draw_rect_roi(''move'')');
    set(ifig_h,'WindowButtonUpFcn',...
               'show_draw_rect_roi(''stop'')');
  case 'move'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    set(rect_h,...
        'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)]);
    set(rect_h,...
        'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)]);
  case 'stop'
    % change the move and buttonup calbacks
    set(ifig_h,'WindowButtonMotionFcn','show_update_pointer');
    set(ifig_h,'WindowButtonUpFcn','');
    % now do the stuff we'd do for a move also
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    set(rect_h,...
        'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)]);
    set(rect_h,...
        'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)]);
    % clear the persistents
    ifig_h=[];
    image_axes_h=[];
    % now add the roi to the list
    show_add_roi(rect_h);
end  % switch






