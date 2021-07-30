function f()

% throw up the dialog box
[filename,pathname]=uigetfile('*.rpb','Load ROIs from File...');
if isnumeric(filename) | isnumeric(pathname)
  % this happens if user hits Cancel
  return;
end

% depending on the extension, call the appropriate loader
len=length(filename);
if strcmp(filename(len-3:len),'.rpb')
  show_load_rois_from_rpb(filename,pathname);
else
  show_load_rois_from_dat(filename,pathname);
end