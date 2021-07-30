function f()

% get the ROI info from the figure state
colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
border_h=get_userdata(colorbar_axes_h,'border_h');
label_h=get_userdata(colorbar_axes_h,'label_h');
n_rois=length(border_h);

% throw up the dialog box
[filename,pathname]=uiputfile('*.rpb','Save ROIs to File...');
if isnumeric(filename) | isnumeric(pathname)
  % this happens if user hits Cancel
  return;
end
full_filename=strcat(pathname,filename);

%
% Write the ROI borders to the file
%

% open the file for writing
fid=fopen(filename,'w','ieee-be');
if (fid == -1)
  errordlg(sprintf('Unable to open file %s',filename),...
           'File Error');
  return;
end

% write the number of ROIs
count=fwrite(fid,n_rois,'uint32');
if (count ~= 1)
  errordlg(sprintf('Error writing ROIs to file %s',filename),...
           'File Error');
  fclose(fid);
  delete(filename);
  return;
end

% for each ROI, write a label and a vertex list
for j=1:n_rois
  % first the label
  label_string=get(label_h(j),'String');
  n_chars=length(label_string);
  count=fwrite(fid,n_chars,'uint32');
  count=fwrite(fid,label_string,'uchar');
  if (count ~= n_chars)
    errordlg(sprintf('Error writing ROIs to file %s',filename),...
             'File Error');
    fclose(fid);
    delete(filename);
    return;
  end
  % then the vertex list
  x=get(border_h(j),'XData');
  y=get(border_h(j),'YData');
  n_vertices=length(x);  % NB: last one is a repeat of first
  % we want each vertex to be a col
  vl=zeros(2,n_vertices);
  vl(1,:)=x;
  vl(2,:)=y;
  % write it
  count=fwrite(fid,n_vertices,'uint32');
  count=fwrite(fid,vl,'float32');
  if (count ~= 2*n_vertices)
    errordlg(sprintf('Error writing ROIs to file %s',filename),...
             'File Error');
    fclose(fid);
    delete(filename);
    return;
  end
end  

% close the file
fclose(fid);

