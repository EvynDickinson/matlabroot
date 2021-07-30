function f()

persistent above_image;

% init above_image
if isempty(above_image)
  mode=get_userdata(gcbf,'mode');
  if strcmp(mode,'select') | strcmp(mode,'move_all')
    above_image=0;
    set(gcbf,'Pointer','arrow');
  else  
    image_axes_h=findobj(gcbf,'Tag','image_axes_h');
    cp=get(image_axes_h,'CurrentPoint');
    cp=cp(1,1:2);
    xlim=get(image_axes_h,'XLim');
    ylim=get(image_axes_h,'YLim');
    above_image=xlim(1)<=cp(1) & ...
                cp(1)<=xlim(2) & ...
                ylim(1)<=cp(2) & ...
                cp(2)<=ylim(2);
    if above_image
      set(gcbf,'Pointer','crosshair');
    else
      set(gcbf,'Pointer','arrow');
    end
  end
end

% if we are above the image now and weren't before, or vice-versa, change
% the pointer appropriately
mode=get_userdata(gcbf,'mode');
if ~( strcmp(mode,'select') | strcmp(mode,'move_all') )
  image_axes_h=findobj(gcbf,'Tag','image_axes_h');
  cp=get(image_axes_h,'CurrentPoint');
  cp=cp(1,1:2);
  xlim=get(image_axes_h,'XLim');
  ylim=get(image_axes_h,'YLim');
  above_image_now=xlim(1)<=cp(1) & ...
                  cp(1)<=xlim(2) & ...
                  ylim(1)<=cp(2) & ...
                  cp(2)<=ylim(2);
  if above_image_now~=above_image
    if above_image_now
      set(gcbf,'Pointer','crosshair');
    else
      set(gcbf,'Pointer','arrow');
    end
    above_image=above_image_now;
  end
end

