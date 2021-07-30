function f(point,anchor)

% get the figure handle
figure_h=gcbf;

% replace the point and anchor with min_corner and max_corner
min_corner=[min(point(1),anchor(1)) min(point(2),anchor(2))];
max_corner=[max(point(1),anchor(1)) max(point(2),anchor(2))];

% make sure there aren't zero pels in the rectangle
if ((min_corner(1)<max_corner(1))&(min_corner(2)<max_corner(2)))
  image_axes_h=findobj(figure_h,'Tag','image_axes_h');
  image_h=findobj(figure_h,'Tag','image_h');
  n_rois=get_userdata(figure_h,'n_rois');
  if n_rois==0
    set(image_h,'EraseMode','normal');
  end
  set(image_axes_h,'XLim',[min_corner(1) max_corner(1)]);
  set(image_axes_h,'YLim',[min_corner(2) max_corner(2)]);
  if n_rois==0
    set(image_h,'EraseMode','none');
  end
  fprintf(1,'Current view:  xlim:[%f %f]  ylim:[%f %f]\n',...
          min_corner(1),max_corner(1),min_corner(2),max_corner(2));
end    
