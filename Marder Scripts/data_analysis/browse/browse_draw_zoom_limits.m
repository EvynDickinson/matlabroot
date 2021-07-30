function f(action)

persistent fig_h;
persistent axes_h;
persistent anchor;
persistent anchor_line_h;
persistent point_line_h;

switch(action)
  case 'start'
    fig_h=gcbf;
    axes_hs=get_userdata(fig_h,'axes_hs');
    axes_h=axes_hs(find(axes_hs==gcbo));
    cp=get(axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    anchor=point;
    % create new limit lines
    yl=ylim(axes_h);
    anchor_line_h=...
      line('Parent',axes_h,...
           'Color',[0.25 0.25 0.25],...
           'Tag','anchor_line_h',...
           'XData',[anchor(1) anchor(1)],...
           'YData',yl,...
           'ZData',[2 2]);
    point_line_h=...
      line('Parent',axes_h,...
           'Color',[0.25 0.25 0.25],...
           'Tag','point_line_h',...
           'XData',[point(1) point(1)],...
           'YData',yl,...
           'ZData',[2 2]);
    % set the callbacks for the drag
    set(fig_h,'WindowButtonMotionFcn',...
               'browse_draw_zoom_limits(''move'')');
    set(fig_h,'WindowButtonUpFcn',...
               'browse_draw_zoom_limits(''stop'')');
  case 'move'
    cp=get(axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    set(point_line_h,...
        'XData',[point(1) point(1)]);
  case 'stop'
    % change the move and buttonup calbacks
    set(fig_h,'WindowButtonMotionFcn','');
    set(fig_h,'WindowButtonUpFcn','');
    % now do the stuff we'd do for a move also
    cp=get(axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    set(point_line_h,...
        'XData',[point(1) point(1)]);
    % clear the persistents
    fig_h=[];
    axes_h=[];
    % now do the zoom
    browse_tlim([anchor(1) point(1)]);
    % now delete the lines
    delete(anchor_line_h);
    delete(point_line_h);    
end  % switch






