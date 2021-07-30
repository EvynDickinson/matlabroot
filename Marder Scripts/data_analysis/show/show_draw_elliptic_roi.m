function f(action)

persistent ifig_h;
persistent image_axes_h;
persistent anchor;
persistent ellipse_h;
persistent n_segs;
persistent sintheta;
persistent costheta;

% init the persistents that point to the image figure and the image axes, if
% this is the first time through this function
if isempty(ifig_h)
  ifig_h=gcbf;
end  
if isempty(image_axes_h)
  image_axes_h=findobj(ifig_h,'Tag','image_axes_h');
end  
if isempty(n_segs)
  n_segs=40;
  theta=linspace(0,2*pi,n_segs+1);
  sintheta=sin(theta);
  costheta=cos(theta);
end

% the big switcheroo  
switch(action)
  case 'start'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    anchor=point;
    % create a new 'ellipse'
    ellipse_h=...
      line('Parent',image_axes_h,...
           'Color',[1 0 0],...
           'Tag','border_h',...
           'XData',repmat(anchor(1),[1 n_segs+1]),...
           'YData',repmat(anchor(2),[1 n_segs+1]),...
           'ZData',repmat(2,[1 n_segs+1]),...
           'ButtonDownFcn','show_callback');
    % set the callbacks for the drag
    set(ifig_h,'WindowButtonMotionFcn',...
               'show_draw_elliptic_roi(''move'')');
    set(ifig_h,'WindowButtonUpFcn',...
               'show_draw_elliptic_roi(''stop'')');
  case 'move'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    center=(point+anchor)/2;
    a=abs(point(1)-anchor(1))/2; b=abs(point(2)-anchor(2))/2;
    dx=a*costheta;
    dy=b*sintheta;
    set(ellipse_h,'XData',center(1)+dx);
    set(ellipse_h,'YData',center(2)+dy);
  case 'stop'
    % change the move and buttonup calbacks
    set(ifig_h,'WindowButtonMotionFcn','show_update_pointer');
    set(ifig_h,'WindowButtonUpFcn','');
    % now do the stuff we'd do for a move also
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    center=(point+anchor)/2;
    a=abs(point(1)-anchor(1))/2; b=abs(point(2)-anchor(2))/2;
    dx=a*costheta;
    dy=b*sintheta;
    set(ellipse_h,'XData',center(1)+dx);
    set(ellipse_h,'YData',center(2)+dy);
    % clear the persistents
    ifig_h=[];
    image_axes_h=[];
    % now add the roi to the list
    show_add_roi(ellipse_h);
end  % switch






