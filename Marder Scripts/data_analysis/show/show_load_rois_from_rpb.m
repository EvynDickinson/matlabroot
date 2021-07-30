function f(filename,pathname)

%
% load in the ROI data from the file, w/ error checking
%

% open the file
full_filename=strcat(pathname,filename);
fid=fopen(full_filename,'r','ieee-be');
if (fid == -1)
  errordlg(sprintf('Unable to open file %s',filename),...
           'File Error');
  return;
end

% read the number of rois
[n_rois,count]=fread(fid,1,'uint32');
%n_rois
if (count ~= 1)
  errordlg(sprintf('Error loading ROIs from file %s',filename),...
           'File Error');
  fclose(fid);
  return;
end

% dimension cell arrays to hold the ROI labels and vertex lists
labels=cell(n_rois,1);
borders=cell(n_rois,1);

% for each ROI, read the label and the vertex list
for j=1:n_rois
  % the label
  [n_chars,count]=fread(fid,1,'uint32');
  if (count ~= 1)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  [temp,count]=fread(fid,[1 n_chars],'uchar');
  if (count ~= n_chars)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  labels{j}=char(temp);
  % the vertex list
  [n_vertices,count]=fread(fid,1,'uint32');
  if (count ~= 1)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  this_border=zeros(2,n_vertices);
  [this_border,count]=fread(fid,[2 n_vertices],'float32');
  %this_border
  if (count ~= 2*n_vertices)
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'Show File Error');
    fclose(fid);
    return;
  end
  borders{j}=this_border;
end

% close the file
fclose(fid);

% clear the old roi borders & labels
colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
glom=get(colorbar_axes_h,'UserData');
delete(glom.border_h);
delete(glom.label_h);

% generate new ROIs state
image_axes_h=findobj(gcbf,'Tag','image_axes_h');
roi_ids=(1:n_rois)';
label_h=zeros(n_rois,1);
border_h=zeros(n_rois,1);
for j=1:n_rois
  this_border=borders{j};
  com=border_com(this_border);
  label_h(j)=...
    text('Parent',image_axes_h,...
         'Position',[com(1) com(2) 1],...
         'String',labels{j},...
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Color',[0 0 1],...
         'Tag','label_h',...
         'Clipping','on',...
         'ButtonDownFcn','show_callback');
  border_h(j)=...
    line('Parent',image_axes_h,...
         'Color',[0 0 1],...
         'Tag','border_h',...
         'XData',this_border(1,:),...
         'YData',this_border(2,:),...
         'ZData',repmat(1,[1 size(this_border,2)]),...
         'ButtonDownFcn','show_callback');
end
   
% write the new ROI into the figure state
set_userdata(gcbf,'n_rois',n_rois);
set_userdata(gcbf,'selected_roi_index',NaN);
set(colorbar_axes_h,'UserData',struct('roi_ids',roi_ids,....
                                      'border_h',border_h,...
                                      'label_h',label_h));
% modify ancillary crap
if n_rois>0
  % need to set image erase mode to normal, since now there's something
  % other than the image in that image axes
  image_h=findobj(gcbf,'Tag','image_h');
  set(image_h,'EraseMode','normal');
  % disable the delete (selected) roi menu, since no ROI
  % is currently sellected
  delete_roi_menu_h=findobj(gcbf,'Tag','delete_roi_menu_h');
  set(delete_roi_menu_h,'Enable','off');
  % enable the delete ROIS menu item
  delete_all_rois_menu_h=...
    findobj(gcbf,'Tag','delete_all_rois_menu_h');
  set(delete_all_rois_menu_h,'Enable','on');
  % enable the save ROIS to file menu item
  save_rois_to_file_menu_h=...
    findobj(gcbf,'Tag','save_rois_to_file_menu_h');
  set(save_rois_to_file_menu_h,'Enable','on');
  % enable the hide ROIs menu item
  hide_rois_menu_h=...
    findobj(gcbf,'Tag','hide_rois_menu_h');
  set(hide_rois_menu_h,'Enable','on');
  show_set_hide_rois(0);
  % enable the select ROI button
  select_button_h=...
    findobj(gcbf,'Tag','select_button_h');
  set(select_button_h,'Enable','on');
  % enable the move all ROIs button
  move_all_button_h=...
    findobj(gcbf,'Tag','move_all_button_h');
  set(move_all_button_h,'Enable','on');
end  
