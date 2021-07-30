function f(tag)

% get handles of figs
ifig_h=gcbf;

% switch on the tag
switch(tag)
  case 'min_max_menu_h'
    image_h=findobj(ifig_h,'Tag','image_h');
    data=get(image_h,'UserData');
    d_min=min(data(:));
    d_max=max(data(:));
    cb_min_string=sprintf('%.4e',d_min);
    cb_max_string=sprintf('%.4e',d_max);
    show_set_colorbar_bounds(cb_min_string,cb_max_string);     
  case 'five_95_menu_h'
    image_h=findobj(ifig_h,'Tag','image_h');
    data=get(image_h,'UserData');
    data=data(:);
    n_bins=1000;
    [h,t]=hist(data,n_bins);
    ch=cumsum(h);
    % figure; plot(t,ch);
    cb_min=crossing_times(t,ch,0.05*ch(n_bins));
    cb_max=crossing_times(t,ch,0.95*ch(n_bins));
    cb_min_string=sprintf('%.4e',cb_min);
    cb_max_string=sprintf('%.4e',cb_max);
    show_set_colorbar_bounds(cb_min_string,cb_max_string);
  case 'abs_max_menu_h'
    image_h=findobj(ifig_h,'Tag','image_h');
    data=get(image_h,'UserData');
    cb_max=max(abs(data(:)));
    cb_min_string=sprintf('%.4e',-cb_max);
    cb_max_string=sprintf('%.4e',+cb_max);
    show_set_colorbar_bounds(cb_min_string,cb_max_string);     
  case 'ninety_symmetric_menu_h'
    % need to fix this, since what it does now is useless
    image_h=findobj(ifig_h,'Tag','image_h');
    data=get(image_h,'UserData');
    data=abs(data(:));
    n_bins=1000;
    [h,t]=hist(data,n_bins);
    ch=cumsum(h);
    % figure; plot(t,ch);
    cb_max=crossing_times(t,ch,0.9*ch(n_bins));
    cb_min=-cb_max;
    cb_min_string=sprintf('%.4e',cb_min);
    cb_max_string=sprintf('%.4e',cb_max);
    show_set_colorbar_bounds(cb_min_string,cb_max_string);
  case 'colorbar_edit_bounds_menu_h'
    % get the current y min and y max strings
    cb_min_string=get_userdata(ifig_h,'colorbar_min_string');
    cb_max_string=get_userdata(ifig_h,'colorbar_max_string');    
    % throw up the dialog box
    bounds=inputdlg({ 'Colorbar Maximum:' , 'Colorbar Minimum:' },...
                    'Edit Colorbar Bounds...',...
                    1,...
                    { cb_max_string , cb_min_string },...
                    'off');
    if isempty(bounds) 
      break; 
    end
    % break out the returned cell array                
    new_cb_max_string=bounds{1};
    new_cb_min_string=bounds{2};
    % convert all these strings to real numbers
    cb_min=str2num(cb_min_string);
    cb_max=str2num(cb_max_string);
    new_cb_min=str2num(new_cb_min_string);
    new_cb_max=str2num(new_cb_max_string);
    % if new values are kosher, change colorbar bounds
    if ~isempty(new_cb_min) & ~isempty(new_cb_max) & ...
       isfinite(new_cb_min) & isfinite(new_cb_max) & ...
       (new_cb_max>new_cb_min)
      show_set_colorbar_bounds(new_cb_min_string,new_cb_max_string);
    end
end  



