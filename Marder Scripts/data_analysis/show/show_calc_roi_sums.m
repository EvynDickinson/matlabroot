function [y,n_pels_mask] = f()

% get handle of fig
ifig_h=gcbf;

% get the data
image_h=findobj(ifig_h,'Tag','image_h');
data=get(image_h,'UserData');

% get dims
n_rows=size(data,1);
n_cols=size(data,2);
n_frames=size(data,3);

% reshape the data for fast indexing with masks
n_ppf=n_rows*n_cols;
data=reshape(data,[n_ppf n_frames]);

% get the roi info
n_rois=get_userdata(ifig_h,'n_rois');
colorbar_axes_h=findobj(ifig_h,'Tag','colorbar_axes_h');
border_h=get_userdata(colorbar_axes_h,'border_h');

% translate the borders into masks
mask=zeros(n_rows,n_cols,n_rois);
template=zeros(n_rows,n_cols);
for k=1:n_rois
  x=get(border_h(k),'XData');
  y=get(border_h(k),'YData');
  mask(:,:,k)=roipoly(template,x,y);
end

% reshape the mask for fast mask indexing
mask=reshape(mask,[n_ppf n_rois]);

% calc the raw y
y=zeros(n_frames,n_rois);
n_pels_mask=zeros(1,n_rois);
for j=1:n_rois
  this_mask=mask(:,j);
  n_pels_mask(j)=sum(this_mask);
  y(:,j)=sum(data(logical(this_mask),:),1)';
end

