function f(action)

persistent ifig_h;
persistent image_axes_h;
persistent anchor;
persistent rect_h;

% init the persistents that point to the image figure and the image 
% axes, if this is the first time through this function
if isempty(ifig_h)
  ifig_h=gcbf;
end  
if isempty(image_axes_h)
  image_axes_h=findobj(ifig_h,'Tag','image_axes_h');
end  

switch(action)
  case 'start'
    point=show_nearest_visible_corner(get(image_axes_h,'CurrentPoint'));
    anchor=point;
    % create a new rectangle
    rect_h=...
      line('Parent',image_axes_h,...
           'Color',[0 0 1],...
           'Tag','border_h',...
           'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)],...
           'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)],...
           'ZData',[1 1 1 1 1],...
           'ButtonDownFcn','show_callback');
    % set the callbacks for the drag
    set(ifig_h,'WindowButtonMotionFcn',...
               'show_draw_zoom_rect(''move'')');
    set(ifig_h,'WindowButtonUpFcn',...
               'show_draw_zoom_rect(''stop'')');
  case 'move'
    point=show_nearest_visible_corner(get(image_axes_h,'CurrentPoint'));
    set(rect_h,...
        'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)]);
    set(rect_h,...
        'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)]);
  case 'stop'
    % change the move and buttonup calbacks
    set(ifig_h,'WindowButtonMotionFcn','show_update_pointer');
    set(ifig_h,'WindowButtonUpFcn','');
    % now do the stuff we'd do for a move also
    point=show_nearest_visible_corner(get(image_axes_h,'CurrentPoint'));
    set(rect_h,...
        'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)]);
    set(rect_h,...
        'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)]);
    % clear the persistents
    ifig_h=[];
    image_axes_h=[];
    % do the zoom
    delete(rect_h)
    show_zoom_in(point,anchor);
end  % switch






