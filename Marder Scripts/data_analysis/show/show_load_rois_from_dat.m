function f(filename,pathname)

% what follows is basically load_binary_data, but with errors in dialog
% boxes instead of to stderr
full_filename=strcat(pathname,filename);
fid=fopen(full_filename,'r','ieee-be');
if (fid == -1)
  errordlg(sprintf('Unable to open file %s',filename),...
           'File Error');
  return;
end
[rank,count]=fread(fid,1,'float64');
if (count ~= 1)
  errordlg(sprintf('Error loading ROIs from file %s',filename),...
           'File Error');
  fclose(fid);
  return;
end
[dims,count]=fread(fid,rank,'float64');
if (count ~= rank)
  errordlg(sprintf('Error loading ROIs from file %s',filename),...
           'Show File Error');
  fclose(fid);
  return;
end
if (rank<=2)
  [mask,count]=fread(fid,dims','float64');
  if (count ~= prod(dims))
    errordlg(sprintf('Error loading ROIs from file %s',filename),...
             'File Error');
    fclose(fid);
    return;
  end
elseif (rank==3)
  for i=1:dims(3)
    [mask(:,:,i),count]=fread(fid,dims(1:2)','float64');
    if (count ~= prod(dims(1:2)))
      errordlg(sprintf('Error loading ROI %d from from file %s',...
                       i,filename),...
               'File Error');
      fclose(fid);
      return;
    end 
  end
else
  errordlg(sprintf('File %s is not a valid stack of ROIs',...
                   filename),...
           'File Error');
  fclose(fid);
  return;
end
fclose(fid);

% make sure the mask is the right shape
image_axes_h=findobj(gcbf,'Tag','image_axes_h');
indexed_data=get(image_axes_h,'UserData');
n_rows=size(indexed_data,1);
n_cols=size(indexed_data,2);
clear indexed_data;
if size(mask,1)~=n_rows | size(mask,2)~=n_cols
  errordlg(...
    sprintf('ROIs in %s are not the same size as the optical data',...
            filename),...
    'File Error');
  fclose(fid);
  return;
end

% clear the old roi borders & labels
colorbar_axes_h=findobj(gcbf,'Tag','colorbar_axes_h');
glom=get(colorbar_axes_h,'UserData');
delete(glom.border_h);
delete(glom.label_h);

% generate new ROIs state
n_rois=size(mask,3);
roi_ids=(1:n_rois)';
label_h=zeros(n_rois,1);
border_h=zeros(n_rois,1);
for j=1:n_rois
  this_mask=mask(:,:,j);
  this_border=mask_border(this_mask);
  com=border_com(this_border);
  label_h(j)=...
    text('Parent',image_axes_h,...
         'Position',[com(1) com(2) 1],...
         'String',show_id_to_alpha(roi_ids(j)),...
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
