function f(action)

persistent ifig_h;
persistent image_axes_h;
persistent anchor;
persistent radius;
persistent A_phi;
persistent ellipse_h;
persistent n_segs;
%persistent theta_hat;
%persistent ones_mat;
persistent circle;

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
  theta_hat=[sin(theta); cos(theta)];
  ones_mat=repmat([1;1],1,n_segs+1);
  circle=ones_mat+theta_hat;
end

% the big switcheroo  
switch(action)
  case 'start'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2)'; 
    anchor=point;
    radius=[0 ; 0];
    phi=0;
    cos_phi=cos(phi); sin_phi=sin(phi);
    A_phi=[cos_phi -sin_phi ; sin_phi cos_phi];
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
    point=cp(1,1:2)';
    ap=point-anchor;
    if strcmp(get(ifig_h,'SelectionType'),'extend')
      % shift drag
      phi=atan2(ap(2),ap(1))-atan2(radius(2),radius(1));
      cos_phi=cos(phi); sin_phi=sin(phi);
      A_phi=[cos_phi -sin_phi ; sin_phi cos_phi];
    else
      % normal drag
      radius=A_phi\(ap/2);
    end
    ellipse=repmat(anchor,[1 n_segs+1])+A_phi*diag(radius)*circle;
    set(ellipse_h,'XData',ellipse(1,:));
    set(ellipse_h,'YData',ellipse(2,:));
  case 'stop'
    % change the move and buttonup calbacks
    set(ifig_h,'WindowButtonMotionFcn','show_update_pointer');
    set(ifig_h,'WindowButtonUpFcn','');
    % now do the stuff we'd do for a move also
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2)';
    ap=point-anchor;
    if strcmp(get(ifig_h,'SelectionType'),'extend')
      % shift drag
      phi=atan2(ap(2),ap(1))-atan2(radius(2),radius(1));
      cos_phi=cos(phi); sin_phi=sin(phi);
      A_phi=[cos_phi -sin_phi ; sin_phi cos_phi];
    else
      % normal drag
      radius=A_phi\(ap/2);
    end
    ellipse=repmat(anchor,[1 n_segs+1])+A_phi*diag(radius)*circle;
    set(ellipse_h,'XData',ellipse(1,:));
    set(ellipse_h,'YData',ellipse(2,:));
    % clear the persistents
    ifig_h=[];
    image_axes_h=[];
    % now add the roi to the list
    show_add_roi(ellipse_h);
end  % switch

